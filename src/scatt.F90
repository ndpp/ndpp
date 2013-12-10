module scatt_class

  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, header, print_ascii_array
  use sab
  use scattdata_class,  only: scattdata
  use search,           only: binary_search
  use string,           only: to_str

#ifdef HDF5
  use hdf5_interface
#endif

  implicit none

  contains

!===============================================================================
! CALC_SCATT Calculates the group-to-group transfer matrices and scattering
! moments for the current nuclide
!===============================================================================

    subroutine calc_scatt(nuc, energy_bins, scatt_type, order, &
                          scatt_mat, mu_bins, thin_tol, E_grid)
      type(Nuclide), pointer, intent(in)  :: nuc            ! Nuclide
      real(8), intent(in)                 :: energy_bins(:) ! Energy groups
      integer, intent(in)                 :: scatt_type     ! Scattering output type
      integer, intent(in)                 :: order          ! Scattering data order
      real(8), allocatable, intent(inout) :: scatt_mat(:,:,:) ! Unionized Scattering Matrices
      integer, intent(in)                 :: mu_bins        ! Number of angular points
                                                            ! to use during f_{n,MT} conversion
      real(8), intent(in)                 :: thin_tol       ! Thinning tolerance
      real(8), allocatable, intent(in)    :: E_grid(:)      ! Incoming Energy Grid

      type(DistEnergy), pointer :: edist
      type(Reaction),   pointer :: rxn
      integer :: num_tot_rxn
      integer :: i_rxn, i_nested_rxn, imu
      real(8) :: max_err, dmu
      type(ScattData), allocatable, target  :: rxn_data(:)
      type(ScattData), pointer :: mySD
      real(8), allocatable     :: mu_out(:) ! The tabular output mu grid

      ! This routine will parse through each nuc % reaction entry.
      ! For each, it will determine if the rxn is a scattering reaction, and if
      ! so it will set up the required memory spaces (the scatt_data object)
      ! and then perform the distribution calculations, storing the results in
      ! scatt_data.

      ! First we have to find the total number of reactions, including the
      ! nested edists, so we can allocate rxn_data correctly.

      num_tot_rxn = 0
      do i_rxn = 1, nuc % n_reaction
        rxn => nuc % reactions(i_rxn)
        num_tot_rxn = num_tot_rxn + 1
        if (rxn % has_energy_dist) then
          edist => rxn % edist
          do while (associated(edist % next))
            edist => edist % next
            num_tot_rxn = num_tot_rxn + 1
          end do
        end if
      end do

      ! Size rxn_data according to num_tot_rxn
      allocate(rxn_data(num_tot_rxn))

      ! Allocate all the rxn_data objects with the correct energy distro, if
      ! needed.
      i_nested_rxn = 0
      do i_rxn = 1, nuc % n_reaction
        i_nested_rxn = i_nested_rxn + 1

        rxn => nuc % reactions(i_rxn)
        mySD => rxn_data(i_nested_rxn)
        edist => rxn % edist
        call mySD % init(nuc, rxn, edist, energy_bins, scatt_type, order, &
          mu_bins)
        if (associated(rxn % edist)) then
          do while (associated(edist % next))
            edist => edist % next
            i_nested_rxn = i_nested_rxn + 1
            mySD => rxn_data(i_nested_rxn)
            call mySD % init(nuc, rxn, edist, energy_bins, scatt_type, order, &
              mu_bins)
          end do
        end if
      end do

      do i_rxn = 1, num_tot_rxn
        mySD => rxn_data(i_rxn)
        ! Convert the angular distributions from the ACE data to Tabular format
        call mySD % convert_distro()
        nullify(mySD)
      end do

      ! Calculate our equi-width mu points
      if (scatt_type == SCATT_TYPE_TABULAR) then
        allocate(mu_out(order))
        dmu = TWO / real(order - 1, 8)
        do imu = 1, order
          mu_out(imu) = -ONE + real(imu - 1, 8) * dmu
        end do
      end if

      ! Now combine the results on to E_grid
      call calc_scatt_grid(nuc, mu_out, rxn_data, E_grid, scatt_mat)

      ! Now clear rxn_datas members
      do i_rxn = 1, num_tot_rxn
        call rxn_data(i_rxn) % clear()
      end do

    end subroutine calc_scatt

!===============================================================================
! CALC_SCATTSAB Calculates the group-to-group transfer matrices and scattering
! moments for the current thermal scattering table
!===============================================================================

    subroutine calc_scattsab(sab, energy_bins, scatt_type, order, &
                             scatt_mat, mu_bins, thin_tol, E_grid)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: energy_bins(:) ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), allocatable, intent(inout)   :: scatt_mat(:,:,:) ! Unionized Scattering Matrices
      integer, intent(in)                   :: mu_bins        ! Number of angular points
                                                              ! to use during f_{n,MT} conversion
      real(8), intent(in)                   :: thin_tol       ! Thinning tolerance
      real(8), allocatable, intent(in)      :: E_grid(:)      ! Incoming Energy Grid

      real(8), allocatable :: mu_out(:)      ! Tabular output mu grid
      real(8), allocatable :: sab_int_el(:,:,:) ! Integrated SAB elastic data [L, G, NEin]
      real(8), allocatable :: sab_int_inel(:,:,:) ! Integrated SAB inelastic data [L, G, NEin]
      real(8), allocatable :: sig_el(:)      ! Elastic microscopic x/s on E_grid
      real(8), allocatable :: sig_inel(:)    ! Inelastic microscopic x/s on E_grid
      integer :: groups    ! Number of energy groups
      integer :: imu       ! angle bin counter
      real(8) :: dmu       ! angle spacing

      ! Initialize memory spaces
      groups = size(energy_bins) - 1
      allocate(sab_int_el(order + 1, groups, size(E_grid)))
      allocate(sab_int_inel(order + 1, groups, size(E_grid)))

      ! This routine will go through the sab table and for the
      ! inelastic and elastic reactions, it will call the proper routines
      ! to populate the data in scatt_mat.  The routines to be called here
      ! will live in the sab module

      ! Perform integral for elastic and inelastic depending on output type
      ! Calculate our equi-width mu points
      if (scatt_type == SCATT_TYPE_LEGENDRE) then
        call integrate_sab_el(sab, E_grid, energy_bins, scatt_type, order, &
                              sab_int_el, sig_el)
        call integrate_sab_inel(sab, E_grid, energy_bins, scatt_type, order, &
                                sab_int_inel, sig_inel)

      else if (scatt_type == SCATT_TYPE_TABULAR) then
        allocate(mu_out(order))
        dmu = TWO / real(order - 1, 8)
        do imu = 1, order
          mu_out(imu) = -ONE + real(imu - 1, 8) * dmu
        end do
        !!!TD - create tabular versions of the above  routines
        ! call integrate_sab_el_tab(sab, energy_bins, scatt_type, order, sab_int(1))
        ! call integrate_sab_inel_tab(sab, energy_bins, scatt_type, order, sab_int(2))
      end if


      ! Now combine the results on to E_grid
      call combine_sab_grid(sab_int_el, sab_int_inel, sig_el, sig_inel, &
                            scatt_mat)

      ! Clear the space so its ready next time
      deallocate(sab_int_el)
      deallocate(sab_int_inel)
    end subroutine calc_scattsab

!===============================================================================
! CALC_SCATT_GRID Combines all the scattering data points on to one single
! energy grid.
!===============================================================================

    subroutine calc_scatt_grid(nuc, mu_out, rxn_data, E_grid, scatt_mat)
      type(Nuclide), pointer, intent(in)   :: nuc   ! The nuclide of interest
      real(8), intent(inout)               :: mu_out(:) ! The tabular output mu grid
      type(ScattData), intent(inout), target :: rxn_data(:) ! The converted distros
      real(8), allocatable, intent(in)     :: E_grid(:) ! Ein grid
      real(8), allocatable, intent(out)    :: scatt_mat(:,:,:) ! Output scattering matrix

      integer :: iE, NE                          ! Ein counter, # Ein
      integer :: irxn, Nrxn                      ! reaction counter, # Reactions
      integer :: groups, order                   ! # Groups, # Orders
      type(ScattData), pointer, SAVE :: mySD => NULL() ! Current working ScattData object
      real(8) :: norm_tot                        ! Sum of all normalization consts
      integer :: iE_print                 ! iE range to print status of
      integer :: iE_pct, last_iE_pct      ! Current and previous pct complete


      groups = rxn_data(1) % groups
      order = rxn_data(1) % order
      NE = size(E_grid)
      Nrxn = size(rxn_data)
      iE_print = NE / 20
      last_iE_pct = -1

      ! Allocate the scatt_mat according to the groups, order and number of E pts
      allocate(scatt_mat(order + 1, groups, NE))

      ! Step through each Ein and reactions and sum the scattering distros @ Ein
      !$omp parallel do schedule(dynamic,50) num_threads(omp_threads) default(shared),private(iE,mySD,norm_tot,irxn)
      do iE = 1, NE
#ifndef OPENMP
        if (iE_print > 0) then
          if ((mod(iE, iE_print) == 1) .or. (iE == NE)) then
            iE_pct = 100 * iE / NE
            if (iE_pct /= last_iE_pct) then
              message = "    Evaluation " // &
                trim(to_str(100 * iE / NE)) // "% Complete"
              call write_message(7)
            end if
            last_iE_pct = iE_pct
          end if
        end if
#endif
        if (E_grid(iE) <= rxn_data(1) % E_bins(size(rxn_data(1) % E_bins))) then
          scatt_mat(:, :, iE) = ZERO
          norm_tot = ZERO
          do irxn = 1, Nrxn
            mySD => rxn_data(irxn)
            ! If we do not have a scatter reaction, don't score it.
            if (.not. mySD % is_init) cycle
            ! Add the scattering distribution to the union scattering grid
            scatt_mat(:, :, iE) = scatt_mat(:, :, iE) + &
              mySD % interp_distro(mu_out, nuc, E_grid(iE), norm_tot)
          end do

          ! Normalize for later multiplication in the MC code
          if (norm_tot == ZERO) norm_tot = ONE
          scatt_mat(:, :, iE) = scatt_mat(:, :, iE) / norm_tot
        else
          ! This step is taken so that interpolation works OK if the MC code
          ! has a particle with an energy == the top energy group value.
          ! With this step, it has something to interpolate to, and that
          ! value is the same as the Ein value, which will be more accurate.
          scatt_mat(:, :, iE) = scatt_mat(:, :, iE - 1)
        end if
      end do
      !$omp end parallel do

    end subroutine calc_scatt_grid

!===============================================================================
! PRINT_SCATT prints the scattering data to the specified output file
! in the specified format.
!===============================================================================

  subroutine print_scatt(name, lib_format, data, E_grid, tol)
    character(len=*),     intent(in) :: name        ! (hdf5 specific) name of group
    integer,              intent(in) :: lib_format  ! Library output type
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8),              intent(in) :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep

    if (lib_format == ASCII) then
      call print_scatt_ascii(data, E_grid, tol)
    else if (lib_format == BINARY) then
      call print_scatt_bin(data, E_grid, tol)
    else if (lib_format == HUMAN) then
      call print_scatt_human(data, E_grid, tol)
    else if (lib_format == H5) then
      call print_scatt_hdf5(name, data, E_grid, tol)
    end if

  end subroutine print_scatt

!===============================================================================
! PRINT_SCATT_ASCII prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_scatt_ascii(data, E_grid, tol)
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    integer :: gmin, gmax, iE

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line:
    ! <size of incoming energy array, # E pts>
    ! <incoming energy array>
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>

    ! Begin writing:

    ! # energy points
    write(UNIT_NUC,'(I20)') size(E_grid)

    ! <incoming energy array>
    call print_ascii_array(E_grid, UNIT_NUC)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(E_grid)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(data, dim = 2)
        if (data(1, gmin, iE) > tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(data, dim = 2), 1, -1
        if (data(1, gmax, iE) > tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC, '(I20,I20)') 0,0
      else
        write(UNIT_NUC, '(I20,I20)') gmin,gmax
        call print_ascii_array(reshape(data(:, gmin : gmax, iE), (/ &
          size(data, dim=1) * (gmax - gmin + 1)/)), UNIT_NUC)
      end if
    end do

  end subroutine print_scatt_ascii

!===============================================================================
! PRINT_SCATT_HUMAN prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_scatt_human(data, E_grid, tol)
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    integer :: g, gmin, gmax, iE

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line:
    ! <size of incoming energy array, # E pts>
    ! <incoming energy array>
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>

    ! Begin writing:

    ! # energy points
    write(UNIT_NUC,'(I20)') size(E_grid)

    ! <incoming energy array>
    call print_ascii_array(E_grid, UNIT_NUC)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(E_grid)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(data, dim = 2)
        if (data(1, gmin, iE) > tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(data, dim = 2), 1, -1
        if (data(1, gmax, iE) > tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',E_grid(iE), &
          '   gmin = ', 0, '   gmax = ', 0
      else
        write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',E_grid(iE), &
          '   gmin = ', gmin, '   gmax = ', gmax
        do g = gmin, gmax
          write(UNIT_NUC,'(A,I5)') 'outgoing group = ', g
          call print_ascii_array(data(:, g, iE), UNIT_NUC)
        end do
      end if
    end do

  end subroutine print_scatt_human

!===============================================================================
! PRINT_SCATT_BIN prints the scattering data to the specified output file
! in a native Fortran stream format.
!===============================================================================

  subroutine print_scatt_bin(data, E_grid, tol)
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    integer :: g, gmin, gmax, iE

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format:
    ! <size of incoming energy array, # E pts>
    ! <incoming energy array>
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>

    ! Begin writing:

    ! # energy points
    write(UNIT_NUC)  size(E_grid)

    ! <incoming energy array>
    write(UNIT_NUC) E_grid

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(E_grid)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(data, dim = 2)
        if (data(1, gmin, iE) > tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(data, dim = 2), 1, -1
        if (data(1, gmax, iE) > tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC) 0, 0
      else
        write(UNIT_NUC) gmin, gmax
        do g = gmin, gmax
          write(UNIT_NUC) data(:, g, iE)
        end do
      end if

    end do

  end subroutine print_scatt_bin

!===============================================================================
! PRINT_SCATT_HDF5 prints the scattering data to the specified output group
! with the HDF5 library.
!===============================================================================

  subroutine print_scatt_hdf5(name, data, E_grid, tol)
    character(len=*),     intent(in) :: name        ! name of nuclide for group
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    integer :: g, gmin, gmax, iE

#ifdef HDF5
    integer(HID_T) :: orig_group, scatt_group
    character(MAX_FILE_LEN) :: group_name, iE_name, nuc_name
    integer :: period_loc

   ! Create a new hdf5 group for the scatter.
    nuc_name = trim(adjustl(name))
    period_loc = scan(nuc_name, '.')
    nuc_name(period_loc : period_loc) = '/'
    group_name = "/" // trim(adjustl(nuc_name)) // "/scatt"
    orig_group = temp_group
    call hdf5_open_group(group_name)
    scatt_group = temp_group

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format:
    ! <size of incoming energy array, # E pts>
    ! <incoming energy array>
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>

    ! Begin writing:

    ! # energy points !!! Maybe not necessary with HDF5, can I get the size
    ! easily???
    call hdf5_write_integer(temp_group, 'NEin', size(E_grid))

    ! <incoming energy array>
    call hdf5_write_double_1Darray(temp_group, 'Ein', E_grid, size(E_grid))

    call hdf5_close_group()

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(E_grid)
      iE_name = trim(adjustl(group_name)) // "/iE" // trim(adjustl(to_str(iE)))
      call hdf5_open_group(iE_name)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(data, dim = 2)
        if (data(1, gmin, iE) > tol) then
          call hdf5_close_group()
          temp_group = scatt_group
          exit
        end if
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(data, dim = 2), 1, -1
        if (data(1, gmax, iE) > tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        call hdf5_write_integer(temp_group, 'gmin', 0)
        call hdf5_write_integer(temp_group, 'gmax', 0)
      else
        call hdf5_write_integer(temp_group, 'gmin', gmin)
        call hdf5_write_integer(temp_group, 'gmax', gmax)
        call hdf5_write_double_2Darray(temp_group, 'data', data(:, gmin : gmax, iE), &
          (/size(data(:, gmin : gmax, iE), dim = 1), &
          size(data(:, gmin : gmax, iE), dim = 2)/))
      end if
      call hdf5_close_group()
      temp_group = scatt_group
    end do
    call hdf5_close_group()
    temp_group = orig_group
#endif
  end subroutine print_scatt_hdf5

end module scatt_class
