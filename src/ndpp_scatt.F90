module ndpp_scatt
  
  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global,           only: nuclides, message
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, header, print_ascii_array
  use scattdata_class,  only: scattdata
  use search,           only: binary_search
  use string,           only: to_str

  implicit none
  
  contains
 
!===============================================================================
! CALC_SCATT Calculates the group-to-group transfer matrices and scattering
! moments for the current nuclide
!===============================================================================
    
    subroutine calc_scatt(nuc, energy_bins, scatt_type, order, &
                          scatt_mat, mu_bins, thin_tol)
      type(Nuclide), pointer, intent(in)  :: nuc            ! Nuclide
      real(8), intent(in)                 :: energy_bins(:) ! Energy groups
      integer, intent(in)                 :: scatt_type     ! Scattering output type
      integer, intent(in)                 :: order          ! Scattering data order
      real(8), allocatable, intent(inout) :: scatt_mat(:,:,:) ! Unionized Scattering Matrices
      integer, intent(in)                 :: mu_bins        ! Number of angular points
                                                            ! to use during f_{n,MT} conversion
      real(8), intent(in)                 :: thin_tol       ! Thinning tolerance
      
      real(8), pointer     :: E_grid(:)      ! Common Energy Grid
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
        if (rxn % has_energy_dist) then
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
      
      ! Create the energy grid in which the distributions will be calculated on
      E_grid => nuc % energy
      
      ! Combine the reactions to a union grid
      if (scatt_type == SCATT_TYPE_TABULAR) then
        allocate(mu_out(order))
        dmu = TWO / real(order - 1, 8)
        do imu = 1, order
          mu_out(imu) = -ONE + real(imu - 1, 8) * dmu
        end do
      end if
      
      ! Now combine the results on to the nuc % energy grid.
      call calc_scatt_grid(nuc, mu_out, rxn_data, E_grid, scatt_mat)
      
      ! Now clear rxn_datas members
      do i_rxn = 1, num_tot_rxn
        call rxn_data(i_rxn) % clear()
      end do
      
    end subroutine calc_scatt
 
!===============================================================================
! CALC_SCATT_GRID Combines all the scattering data points on to one single 
! energy grid.
!===============================================================================
    
    subroutine calc_scatt_grid(nuc, mu_out, rxn_data, E_grid, scatt_mat)
      type(Nuclide), pointer, intent(in)   :: nuc   ! The nuclide of interest
      real(8), intent(inout)               :: mu_out(:) ! The tabular output mu grid
      type(ScattData), intent(inout), target  :: rxn_data(:) ! The converted distros
      real(8), pointer, intent(out)        :: E_grid(:) ! Ein grid
      real(8), allocatable, intent(out)    :: scatt_mat(:,:,:) ! Output scattering matrix
      
      integer :: iE, NE
      integer :: irxn, Nrxn
      integer :: groups, order
      type(ScattData), pointer :: mySD => NULL()
      
      groups = rxn_data(1) % groups
      order = rxn_data(1) % order
      NE = size(E_grid)
      Nrxn = size(rxn_data)
      
      ! Allocate the scatt_mat according to the groups, order and number of E pts
      
      allocate(scatt_mat(order, groups, NE))
      
      ! Step through each Ein and reactions and sum the scattering distros @ Ein
      do iE = 1, NE
        scatt_mat(:, :, iE) = ZERO
        do irxn = 1, Nrxn
          mySD => rxn_data(irxn)
          if (.not. mySD % is_init) cycle
          ! If the incoming energy point is outside the range of this reaction,
          ! dont score it
          if (E_grid(iE) < mySD % E_grid(1)) cycle
          ! Add the scattering distribution to the union scattering grid
          scatt_mat(:, :, iE) = scatt_mat(:, :, iE) + &
            mySD % interp_distro(mu_out, nuc, E_grid(iE))
        end do
      end do
    end subroutine calc_scatt_grid
    
!===============================================================================
! MIN_VALUE_LOCS Finds the reaction in rxn_data with the next highest Ein, and
! stores the result in hits
!===============================================================================
  
    subroutine min_value_locs(rxn_data, loc, done, hits)
      type(ScattData), intent(in)  :: rxn_data(:) ! Rxn information
      integer, intent(in)    :: loc(:)  ! Location of last-used energy
      logical, intent(in)    :: done(:) ! Flag to state whether a rxn's grid is complete
      logical, intent(inout) :: hits(:) ! The reactions which are the lowest E
      
      integer :: i
      real(8) :: min_val
      
      hits = .false.
      min_val = INFINITY
      
      ! Find the next smallest energy value
      do i = 1, size(rxn_data)
        if (.not. done(i)) then
          if (rxn_data(i) % E_grid(loc(i)) < min_val) &
            min_val = rxn_data(i) % E_grid(loc(i))
        end if
      end do
      
      ! Find all locs  that are within a certain tolerance of the min value
      do i = 1, size(rxn_data)
        if (.not. done(i)) then
          if (abs(rxn_data(i) % E_grid(loc(i)) - min_val) < FP_PRECISION) then
            hits(i) = .true.
          end if
        end if
      end do
    
    end subroutine min_value_locs 
  
!===============================================================================
! PRINT_SCATT prints the scattering data to the specified output file
! in the specified format.
!===============================================================================  
  
  subroutine print_scatt(lib_format, data, E_grid, tol)
    integer,              intent(in) :: lib_format  ! Library output type
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print 
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized Total Energy in
    real(8),              intent(in) :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    
    if (lib_format == ASCII) then
      call print_scatt_ascii(data, E_grid, tol)
    else if (lib_format == BINARY) then
      call print_scatt_bin(data, E_grid, tol)
    else if (lib_format == HDF5) then
      ! TBI
    end if
    
  end subroutine print_scatt
    
!===============================================================================
! PRINT_SCATT_ASCII prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================
     
  subroutine print_scatt_ascii(data, E_grid, tol)
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print 
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized Total Energy in
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    integer :: gmin, gmax, iE
    
    character(MAX_LINE_LEN) :: line
  
    ! Assumes that the file and header information is already printed 
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line: 
    ! <size of incoming energy array, size(E_grid)>
    ! <incoming energy array>
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    
    ! Begin writing:
    
    ! <size(E_grid)>
    line = ''
    write(line,'(I20)') size(E_grid)
    write(UNIT_NUC,'(A)') trim(line)
    
    ! <incoming energy array>
    call print_ascii_array(E_grid, UNIT_NUC)
    
!~     ! < \Sigma_{s,g',l}(Ein) ordered by: order, groups, Ein>
!~     call print_ascii_array(pack(data, .true.), UNIT_NUC)
    
    
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
        call print_ascii_array((/ZERO, ZERO/), UNIT_NUC)
      else
        call print_ascii_array((/real(gmin,8), real(gmax,8), &
          reshape(data(:, gmin : gmax, iE), (/ &
          size(data, dim=1) * (gmax - gmin + 1)/))/), UNIT_NUC)
      end if
    end do
  
    
  end subroutine print_scatt_ascii
  
  !===============================================================================
! PRINT_SCATT_BIN prints the scattering data to the specified output file
! in in native Fortran stream format.
!===============================================================================
     
  subroutine print_scatt_bin(data, E_grid, tol)
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print 
                                                    ! (order x g x Ein)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized Total Energy in
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    integer :: gmin, gmax, iE
    
    ! Assumes that the file and header information is already printed 
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format: 
    ! <size of incoming energy array, size(E_grid)>
    ! <incoming energy array>
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    
    ! Begin writing:
    
    ! <size(E_grid)>
    write(UNIT_NUC)  size(E_grid)
    
    ! <incoming energy array>
    write(UNIT_NUC) E_grid
    
!~     ! < \Sigma_{s,g',l}(Ein) ordered by: order, groups, Ein>
!~     write(UNIT_NUC) data

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
        write(UNIT_NUC) (/ZERO, ZERO/)
      else
        write(UNIT_NUC) (/real(gmin,8), real(gmax,8), &
          reshape(data(:, gmin : gmax, iE), (/ &
          size(data, dim=1) * (gmax - gmin + 1)/))/)
      end if
      
    end do
    
  end subroutine print_scatt_bin
    
end module ndpp_scatt
