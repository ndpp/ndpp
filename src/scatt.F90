module scatt

  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, print_ascii_array
  use sab
  use scattdata_header, only: scattdata
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

    subroutine calc_scatt(nuc, energy_bins, scatt_type, order, mu_bins, &
                          thin_tol, nuscatt, E_grid, scatt_mat, nuscatt_mat)
      type(Nuclide), pointer, intent(in)  :: nuc            ! Nuclide
      real(8), intent(in)                 :: energy_bins(:) ! Energy groups
      integer, intent(in)                 :: scatt_type     ! Scattering output type
      integer, intent(inout)              :: order          ! Scattering data order
      integer, intent(in)                 :: mu_bins        ! Number of angular points
                                                            ! to use during f_{n,MT} conversion
      real(8), intent(in)                 :: thin_tol       ! Thinning tolerance
      logical, intent(in)                 :: nuscatt        ! Whether or not to include nuscatt
      real(8), allocatable, intent(inout) :: E_grid(:)      ! Incoming Energy Grid
      real(8), allocatable, intent(inout) :: scatt_mat(:,:,:) ! Unionized Scattering Matrices
      real(8), allocatable, intent(inout) :: nuscatt_mat(:,:,:) ! Unionized Nu-Scattering Matrices

      type(DistEnergy), pointer :: edist
      type(Reaction),   pointer :: rxn
      integer :: num_tot_rxn
      integer :: i_rxn, i_nested_rxn, imu
      real(8) :: max_err, dmu
      type(ScattData), allocatable, target  :: rxn_data(:)
      type(ScattData), pointer :: mySD
      real(8), allocatable     :: mu_out(:) ! The tabular output mu grid
      type(ScattData), pointer :: inittedSD => null()
      real(8) :: inel_thresh

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
      inel_thresh = E_grid(size(E_grid))
      do i_rxn = 1, nuc % n_reaction
        i_nested_rxn = i_nested_rxn + 1

        rxn => nuc % reactions(i_rxn)
        mySD => rxn_data(i_nested_rxn)
        edist => rxn % edist
        call mySD % init(nuc, rxn, edist, energy_bins, scatt_type, order, mu_bins)
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
        ! Get the threshold for inelastic collisions
        if (mySD % is_init) then
          inittedSD => rxn_data(i_rxn)
          rxn => mySD % rxn
          edist => rxn % edist
          if ((rxn % MT >= N_N1) .and. (rxn % MT < N_NC)) then
            if (edist % data(1) < inel_thresh) then
              inel_thresh = edist % data(1)
            end if
          else if ((rxn % MT == N_2N) .or. (rxn % MT == N_3N) .or.&
                   (rxn % MT == N_4N) .or. (rxn % MT == N_NC)) then
            if (mySD % E_grid(1) < inel_thresh) then
              inel_thresh = mySD % E_grid(1)
            end if
          end if
        end if
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

      ! Add in incoming energies from all reaction channel distributions to
      ! the Ein grid
      call combine_Eins(rxn_data, E_grid)

      ! Add points to aid in interpolation if elastic scattering points
      call add_elastic_Eins(nuc % awr, energy_bins, E_grid)

      ! Expand the Incoming energy grid (E_grid) to include points which will
      ! significantly improve linear interpolation of inelastic level scatter
      ! results.  These results are kind of like stair functions but with
      ! near-linear (depends on f(mu)) ramps inbetween each `step'.  For now
      ! we will combat this by putting EXTEND_PTS per current point above the
      ! threshold for inelastic level scatter to begin
      call extend_grid(inel_thresh, E_grid)

      ! Finally, we need to make sure there are at least EXTEND_PTS
      ! incoming energy points per group.  If not, then add them in.
      call add_pts_per_group(energy_bins, E_grid)

      ! Now combine the results on to E_grid
      call calc_scatt_grid(nuc, mu_out, rxn_data, E_grid, inittedSD % order, &
        energy_bins, nuscatt, scatt_mat, nuscatt_mat)

      ! Now clear rxn_datas members
      do i_rxn = 1, num_tot_rxn
        call rxn_data(i_rxn) % clear()
      end do

    end subroutine calc_scatt

!===============================================================================
! COMBINE_EINS ensures that the incoming energies from all of the different
! reaction channels' incoming energy grids are included in our final grid so
! as not to miss any highly varying information.
!===============================================================================

    subroutine combine_Eins(rxn_data, Ein)
      type(ScattData), target, intent(in) :: rxn_data(:) ! Reaction data
      real(8), allocatable, intent(inout) :: Ein(:)      ! Incoming Energy Grid

      real(8), allocatable      :: new_grid(:)
      real(8), allocatable      :: tmp_grid(:)
      type(ScattData), pointer  :: mySD
      integer :: i_rxn

      allocate(new_grid(1))
      new_grid = Ein(1)

      write(*,*) 'Before', size(Ein)

      do i_rxn = 1, size(rxn_data)
        mySD => rxn_data(i_rxn)
        if (mySD % is_init) then
          call merge(mySD % E_grid, new_grid, tmp_grid)
          deallocate(new_grid)
          allocate(new_grid(size(tmp_grid)))
          new_grid = tmp_grid
        end if
        nullify(mySD)
      end do

      deallocate(tmp_grid)
      allocate(tmp_grid(size(Ein)))
      tmp_grid = Ein
      call merge(new_grid, tmp_grid, Ein)
      deallocate(tmp_grid)
      deallocate(new_grid)

      write(*,*) 'After', size(Ein)

    end subroutine combine_Eins

    subroutine add_elastic_Eins(awr, E_bins, Ein)
      real(8), intent(in)                 :: awr       ! Atomic Weight Ratio
      real(8), intent(in)                 :: E_bins(:) ! Energy groups
      real(8), allocatable, intent(inout) :: Ein(:)    ! Incoming Energy Grid

      real(8), allocatable  :: new_pts(:)
      real(8), allocatable  :: old_grid(:)
      integer               :: g
      real(8)               :: alpha
      integer               :: num_pts

      allocate(old_grid(size(Ein)))
      old_grid = Ein

      alpha = ((awr - ONE) / (awr + ONE))**2

      allocate(new_pts(size(E_bins) - 1))

      num_pts = 0
      do g = 1, size(E_bins) - 1
        new_pts(g) = E_bins(g) / alpha
        ! Now we increment new_pts if the new point is inside our energy range
        ! of interest, and if not, quit this loop.  That way we dont add points
        ! which are above the peak energy we need.
        if (new_pts(g) < E_bins(size(E_bins))) then
          num_pts = num_pts + 1
        else
          exit
        end if
      end do

      call merge(new_pts, old_grid, Ein)
      deallocate(new_pts)
      deallocate(old_grid)
    end subroutine add_elastic_Eins


!===============================================================================
! EXPAND_GRID Expands the Incoming energy grid to include points which will
! significantly improve linear interpolation of inelastic level scatter
! results.  These results are kind of like stair functions but with
! near-linear (depends on f(mu)) ramps inbetween each `step'.  So
! to combat this, we will put an Ein point at the start and end of this
! `ramp'.
!===============================================================================

    subroutine expand_grid(rxn_data, E_bins, Ein)
      type(ScattData), target, intent(in) :: rxn_data(:) ! Reaction data
      real(8), intent(in)                 :: E_bins(:)   ! Energy groups
      real(8), allocatable, intent(inout) :: Ein(:)      ! Incoming Energy Grid

      real(8), allocatable :: new_grid(:,:)
      integer :: groups, i_rxn, g, iE, iMT
      type(DistEnergy), pointer :: edist => null()
      type(Reaction),   pointer :: rxn   => null()
      type(ScattData), pointer  :: mySD  => null()
      real(8) :: Eo, Ei, Ap12, D1, D2, thresh


      groups = size(E_bins) - 1

      allocate(new_grid(2 * groups , N_NC - N_N1))
      ! Set them to this so when we merge them all, duplicates are removed
      new_grid = E_bins(1)

      iMT = 0
      do i_rxn = 1, size(rxn_data)
        mySD => rxn_data(i_rxn)
        if (.not. mySD % is_init) then
          cycle
        end if
        rxn => mySD % rxn
        Ap12 = (mySD % awr + ONE) * (mySD % awr + ONE)
        ! Lets make sure that we have an inelastic level, and that law 3
        ! is used to describe it (which it totally better, but why not check)
        if (((rxn % MT >= N_N1) .and. (rxn % MT < N_NC)) .and. &
            (mySD % law == 3)) then
          edist => mySD % edist
          ! Find which Ein puts the lowest energy point at each group boundary
          D1 = edist % data(1)
          D2 = edist % data(2)
          iMT = iMT + 1
          iE = 0
          do g = 1, groups
            Eo = E_bins(g)
            ! Eo_low and Eo_high yield two possibilities on the group boundary:
            !!! INCORRECT
            Ei = sqrt((Ap12*(-(D1*D2)+Eo)+Ap12*Ap12*D2*(D1*D2+Eo) - &
                      TWO*sqrt(Ap12*Ap12*D2*Eo*(-D1+Ap12*D1*D2+Ap12*Eo))) / &
                      ((Ap12 * D2 - ONE) * (Ap12 * D2 - ONE)))
            ! Lets check for validity (Ei above threshold)
            if (Ei > D1) then
              iE = iE + 1
              new_grid(iE, iMT) = Ei
            end if
            !!! INCORRECT
            Ei = sqrt((Ap12*(-(D1*D2)+Eo)+Ap12*Ap12*D2*(D1*D2+Eo) + &
                      TWO*sqrt(Ap12*Ap12*D2*Eo*(-D1+Ap12*D1*D2+Ap12*Eo))) / &
                      ((Ap12 * D2 - ONE) * (Ap12 * D2 - ONE)))
            ! Lets check for validity (Ei above threshold)
            if (Ei > D1) then
              iE = iE + 1
              new_grid(iE, iMT) = Ei
            end if
          end do

        end if

      end do

    end subroutine expand_grid

!===============================================================================
! ADD_PTS_PER_GROUP Checks to make sure there are EXTEND_PTS per group available.
! If not, then EXTEND_PTS are added to the group.
!===============================================================================

    subroutine add_pts_per_group(E_bins, Ein)
      real(8), intent(in)                 :: E_bins(:)   ! Energy groups
      real(8), allocatable, intent(inout) :: Ein(:)      ! Incoming Energy Grid

      real(8), allocatable :: new_grid(:)
      real(8), allocatable :: temp_grid(:)
      integer :: groups, g, lo, hi, i, j
      real(8) :: dE


      groups = size(E_bins) - 1

      allocate(new_grid(EXTEND_PTS * groups))
      ! Set them to this so when we merge them all, duplicates are removed
      new_grid = E_bins(1)

      hi = 1
      i = 1
      do g = 1, groups
        lo = hi
        hi = binary_search(Ein, size(Ein), E_bins(g + 1))
        if (hi - lo < EXTEND_PTS) then
          dE = (Ein(hi) - Ein(lo)) / real(EXTEND_PTS,8)
          do j = 0, EXTEND_PTS - 1
            new_grid(i + j) = Ein(lo) + real(j,8) * dE
          end do
          i = i + EXTEND_PTS
        end if
      end do

      ! Now we can merge in new_grid(:i-1) with Ein to get our new grid
      if (i > 1) then
        ! Otherwise, nothing needs to happen, E_grid does not change
        call merge(new_grid(:i - 1), Ein, temp_grid)
        deallocate(Ein)
        allocate(Ein(size(temp_grid)))
        Ein = temp_grid
        deallocate(temp_grid)
      end if
      deallocate(new_grid)

    end subroutine add_pts_per_group



!===============================================================================
! EXTEND_GRID increases the number of points present above a threshold (min)
! The number of points to increase is EXTEND_PTS for every point on Ein.
!===============================================================================

    subroutine extend_grid(min, a)
      real(8), intent(in)                 :: min
      real(8), allocatable, intent(inout) :: a(:)
      real(8), allocatable :: temp(:)
      integer :: i, j, k, l
      real(8) :: dE

      ! First lets find where min is
      k = binary_search(a, size(a), min)
      allocate(temp(k + EXTEND_PTS * (size(a) - k) - 1))
      temp(1:k-1) = a(1:k-1)

      j = k
      do i = k, size(a) - 1
        dE = (a(i + 1) - a(i)) / real(EXTEND_PTS,8)
        do l = 0, EXTEND_PTS - 1
          temp(j + l) = a(i) + real(l,8) * dE
        end do
        j = j + EXTEND_PTS
      end do
      temp(size(temp)) = a(size(a))

      deallocate(a)
      allocate(a(size(temp)))
      a = temp
      deallocate(temp)

    end subroutine extend_grid

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

    subroutine calc_scatt_grid(nuc, mu_out, rxn_data, E_grid, order, E_bins, &
                               nuscatt, scatt_mat, nuscatt_mat)
      type(Nuclide), pointer, intent(in)   :: nuc   ! The nuclide of interest
      real(8), intent(inout)               :: mu_out(:) ! The tabular output mu grid
      type(ScattData), intent(inout), target :: rxn_data(:) ! The converted distros
      real(8), allocatable, intent(in)     :: E_grid(:) ! Ein grid
      integer, intent(in)                  :: order     ! Angular order
      real(8), intent(in)                  :: E_bins(:) ! Energy groups
      logical, intent(in)                  :: nuscatt   ! Include nuscatter?
      real(8), allocatable, intent(out)    :: scatt_mat(:,:,:) ! Output scattering matrix
      real(8), allocatable, intent(out)    :: nuscatt_mat(:,:,:) ! Output nu-scattering matrix

      integer :: iE, NE              ! Ein counter, # Ein
      integer :: irxn, Nrxn          ! reaction counter, # Reactions
      integer :: groups              ! # Groups
      type(ScattData), pointer, SAVE :: mySD => NULL() ! Current working ScattData object
      real(8) :: norm_tot            ! Sum of all normalization consts
      integer :: iE_print            ! iE range to print status of
      integer :: iE_pct, last_iE_pct ! Current and previous pct complete
      real(8), allocatable :: temp_scatt(:,:) ! calculated scattering matrix


      groups = size(E_bins) - 1
      NE = size(E_grid)
      Nrxn = size(rxn_data)
      iE_print = NE / 20
      last_iE_pct = -1

      ! Allocate the scatt_mat according to the groups, order and number of E pts
      allocate(scatt_mat(order, groups, NE))

      if (nuscatt) then
        allocate(nuscatt_mat(order, groups, NE))
        nuscatt_mat = ZERO
      end if

      allocate(temp_scatt(order, groups))

      ! Step through each Ein and reactions and sum the scattering distros @ Ein
      !$omp parallel do schedule(dynamic,50) num_threads(omp_threads) &
      !$omp default(shared), private(iE,mySD,norm_tot,irxn,temp_scatt)
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
        if (E_grid(iE) <= E_bins(size(E_bins))) then
          scatt_mat(:, :, iE) = ZERO
          norm_tot = ZERO
          do irxn = 1, Nrxn
            mySD => rxn_data(irxn)
            ! If we do not have a scatter reaction, don't score it.
            if (.not. mySD % is_init) cycle
            ! Some reactions, ENDF-VII.0's Ca-40 e.g., have two angdist Ein
            ! points in a row being the same value. Check for this and just skip
            ! the first point
            if (iE < mySD % NE) then
              if (mySD % E_grid(iE) >= mySD % E_grid(iE + 1)) then
                cycle
              end if
            end if
            ! Add the scattering distribution to the union scattering grid
            temp_scatt = mySD % interp_distro(mu_out, nuc, E_grid(iE), norm_tot)
            scatt_mat(:, :, iE) = scatt_mat(:, :, iE) + temp_scatt
            if (nuscatt) then
              nuscatt_mat(:, :, iE) = nuscatt_mat(:, :, iE) + &
                real(mySD % rxn % multiplicity, 8) * temp_scatt
            end if
          end do

          ! Normalize for later multiplication in the MC code
          if (norm_tot == ZERO) norm_tot = ONE
          scatt_mat(:, :, iE) = scatt_mat(:, :, iE) / norm_tot

          if (nuscatt) then
            nuscatt_mat(:, :, iE) = nuscatt_mat(:, :, iE) / norm_tot
          end if
        else
          ! This step is taken so that interpolation works OK if the MC code
          ! has a particle with an energy == the top energy group value.
          ! With this step, it has something to interpolate to, and that
          ! value is the same as the Ein value, which will be more accurate.
          scatt_mat(:, :, iE) = scatt_mat(:, :, iE - 1)
          if (nuscatt) then
            nuscatt_mat(:, :, iE) = nuscatt_mat(:, :, iE - 1)
          end if
        end if
      end do
      !$omp end parallel do

    end subroutine calc_scatt_grid

!===============================================================================
! PRINT_SCATT prints the scattering data to the specified output file
! in the specified format.
!===============================================================================

  subroutine print_scatt(name, lib_format, E_grid, tol, data, nudata)
    character(len=*),     intent(in) :: name        ! (hdf5 specific) name of group
    integer,              intent(in) :: lib_format  ! Library output type
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8),              intent(in) :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, optional, intent(in) :: nudata(:,:,:) ! Nu-Scatt data to print
                                                    !            (order x g x Ein)

    if (present(nudata)) then
      if (lib_format == ASCII) then
        call print_scatt_ascii(E_grid, tol, data, nudata)
      else if (lib_format == BINARY) then
        call print_scatt_bin(E_grid, tol, data, nudata)
      else if (lib_format == HUMAN) then
        call print_scatt_human(E_grid, tol, data, nudata)
      else if (lib_format == H5) then
        call print_scatt_hdf5(E_grid, tol, data, nudata)
      end if
    else
      if (lib_format == ASCII) then
        call print_scatt_ascii(E_grid, tol, data)
      else if (lib_format == BINARY) then
        call print_scatt_bin(E_grid, tol, data)
      else if (lib_format == HUMAN) then
        call print_scatt_human(E_grid, tol, data)
      else if (lib_format == H5) then
        call print_scatt_hdf5(E_grid, tol, data)
      end if
    end if

  end subroutine print_scatt

!===============================================================================
! PRINT_SCATT_ASCII prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_scatt_ascii(E_grid, tol, data, nudata)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, optional, intent(in) :: nudata(:,:,:) ! Nu-Scatt data to print
                                                                ! (order x g x Ein)

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

    if (present(nudata)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(E_grid)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nudata, dim = 2)
          if (nudata(1, gmin, iE) > tol) exit
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nudata, dim = 2), 1, -1
          if (nudata(1, gmax, iE) > tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          write(UNIT_NUC, '(I20,I20)') 0,0
        else
          write(UNIT_NUC, '(I20,I20)') gmin,gmax
          call print_ascii_array(reshape(nudata(:, gmin : gmax, iE), (/ &
            size(nudata, dim=1) * (gmax - gmin + 1)/)), UNIT_NUC)
        end if
      end do
    end if

  end subroutine print_scatt_ascii

!===============================================================================
! PRINT_SCATT_HUMAN prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_scatt_human(E_grid, tol, data, nudata)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, optional, intent(in) :: nudata(:,:,:) ! Nu-Scatt data to print
                                                                ! (order x g x Ein)

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

    if (present(nudata)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(E_grid)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nudata, dim = 2)
          if (nudata(1, gmin, iE) > tol) exit
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nudata, dim = 2), 1, -1
          if (nudata(1, gmax, iE) > tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',E_grid(iE), &
            '   gmin = ', 0, '   gmax = ', 0
        else
          write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',E_grid(iE), &
            '   gmin = ', gmin, '   gmax = ', gmax
          do g = gmin, gmax
            write(UNIT_NUC,'(A,I5)') 'outgoing group = ', g
            call print_ascii_array(nudata(:, g, iE), UNIT_NUC)
          end do
        end if
      end do
    end if

  end subroutine print_scatt_human

!===============================================================================
! PRINT_SCATT_BIN prints the scattering data to the specified output file
! in a native Fortran stream format.
!===============================================================================

  subroutine print_scatt_bin(E_grid, tol, data, nudata)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, optional, intent(in) :: nudata(:,:,:) ! Nu-Scatt data to print
                                                                ! (order x g x Ein)

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

    if (present(nudata)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(E_grid)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nudata, dim = 2)
          if (nudata(1, gmin, iE) > tol) exit
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nudata, dim = 2), 1, -1
          if (nudata(1, gmax, iE) > tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          write(UNIT_NUC) 0, 0
        else
          write(UNIT_NUC) gmin, gmax
          do g = gmin, gmax
            write(UNIT_NUC) nudata(:, g, iE)
          end do
        end if
      end do
    end if

  end subroutine print_scatt_bin

!===============================================================================
! PRINT_SCATT_HDF5 prints the scattering data to the specified output group
! with the HDF5 library.
!===============================================================================

  subroutine print_scatt_hdf5(E_grid, tol, data, nudata)
    real(8), allocatable, intent(in) :: E_grid(:)   ! Unionized E_{in} grid
    real(8), intent(in)              :: tol         ! Minimum grp-to-grp prob'y
                                                    ! to keep
    real(8), allocatable, intent(in) :: data(:,:,:) ! Scatt data to print
                                                    ! (order x g x Ein)
    real(8), allocatable, optional, intent(in) :: nudata(:,:,:) ! Nu-Scatt data to print
                                                                ! (order x g x Ein)

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

    if (present(nudata)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(E_grid)
        iE_name = trim(adjustl(group_name)) // "/iE" // trim(adjustl(to_str(iE)))
        call hdf5_open_group(iE_name)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nudata, dim = 2)
          if (nudata(1, gmin, iE) > tol) then
            call hdf5_close_group()
            temp_group = scatt_group
            exit
          end if
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nudata, dim = 2), 1, -1
          if (nudata(1, gmax, iE) > tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          call hdf5_write_integer(temp_group, 'nu_gmin', 0)
          call hdf5_write_integer(temp_group, 'nu_gmax', 0)
        else
          call hdf5_write_integer(temp_group, 'nu_gmin', gmin)
          call hdf5_write_integer(temp_group, 'nu_gmax', gmax)
          call hdf5_write_double_2Darray(temp_group, 'nu_data', nudata(:, gmin : gmax, iE), &
            (/size(nudata(:, gmin : gmax, iE), dim = 1), &
            size(nudata(:, gmin : gmax, iE), dim = 2)/))
        end if
        call hdf5_close_group()
        temp_group = scatt_group
      end do
    end if

    call hdf5_close_group()
    temp_group = orig_group
#endif
  end subroutine print_scatt_hdf5

end module scatt
