module scatt

  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, print_ascii_array, &
                              print_ascii_integer_array
  use sab
  use scattdata_header, only: scattdata
  use search,           only: binary_search
  use string,           only: to_str

#ifdef HDF5
  use hdf5_interface
#endif

#ifdef OPENMP
  use omp_lib
#endif

  implicit none

  contains

!===============================================================================
! CALC_SCATT Calculates the group-to-group transfer matrices and scattering
! moments for the current nuclide
!===============================================================================

    subroutine calc_scatt(nuc, energy_bins, scatt_type, order, mu_bins, &
                          nuscatt, Ein_el, Ein_inel, el_mat, inel_mat, nuinel_mat, &
                          normalization)
      type(Nuclide), pointer, intent(in)  :: nuc            ! Nuclide
      real(8), intent(in)                 :: energy_bins(:) ! Energy groups
      integer, intent(in)                 :: scatt_type     ! Scattering output type
      integer, intent(inout)              :: order          ! Scattering data order
      integer, intent(in)                 :: mu_bins        ! Number of angular points
                                                            ! to use during f_{n,MT} conversion
      logical, intent(in)                 :: nuscatt        ! Whether or not to include nuscatt
      real(8), allocatable, intent(inout) :: Ein_el(:)      ! Incoming Energy Grid for elastic
      real(8), allocatable, intent(inout) :: Ein_inel(:)    ! Incoming Energy Grid for inelastic
      real(8), allocatable, intent(inout) :: el_mat(:,:,:)   ! Unionized Elastic Matrices
      real(8), allocatable, intent(inout) :: inel_mat(:,:,:) ! Unionized Inelastic Matrices
      real(8), allocatable, intent(inout) :: nuinel_mat(:,:,:) ! Unionized Nu-Inelastic Matrices
      real(8), allocatable, intent(inout) :: normalization(:)  ! Norm. Constant for inelastic Data

      type(DistEnergy), pointer :: edist
      type(Reaction),   pointer :: rxn
      integer :: num_tot_rxn
      integer :: i_rxn, i_nested_rxn, imu
      real(8) :: dmu
      type(ScattData), allocatable, target  :: rxn_data(:)
      type(ScattData), pointer :: mySD
      real(8), allocatable     :: mu_out(:) ! The tabular output mu grid
      type(ScattData), pointer :: inittedSD => null()
      real(8) :: inel_thresh
      real(8) :: cutoff

      cutoff = ZERO

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
      inel_thresh = energy_bins(size(energy_bins))
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
        ! Also get the threshold for inelastic collisions
        if (mySD % is_init) then
          inittedSD => rxn_data(i_rxn)
          rxn => mySD % rxn
          edist => rxn % edist
          ! Find the threshold of the inelastic reactions so we can tell
          ! where to begin placing additional points for interpolation.
          if (rxn % MT == ELASTIC) then
            cutoff = mySD % freegas_cutoff
          else
            inel_thresh = nuc % energy(rxn % threshold)
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

      ! Build our incoming energy grid to use. There will be one for elastic
      ! and another for inelastic
      call create_Ein_grid(rxn_data, energy_bins, nuc % energy, nuc % awr, &
                           nuc % kT, cutoff, inel_thresh, Ein_el, Ein_inel)

      ! Now combine the results on to our Energy grids for each case
      call calc_elastic_grid(nuc, mu_out, rxn_data, Ein_el, inittedSD % order, &
                             energy_bins, el_mat)

      call calc_inelastic_grid(nuc, mu_out, rxn_data, Ein_inel, &
                               inittedSD % order, energy_bins, nuscatt, &
                               inel_mat, nuinel_mat, normalization)

      ! Now clear rxn_datas members
      do i_rxn = 1, num_tot_rxn
        call rxn_data(i_rxn) % clear()
      end do

    end subroutine calc_scatt

!===============================================================================
! CREATE_EIN_GRID creates the incoming energy grid to be provided to in the
! output library.  It effectively takes the nuclide's energy grid then calls
! combine_eins, add_elastic_Eins, add_inelastic_Eins, add_pts_per_group,
! and add_one_more_point.
!===============================================================================

    subroutine create_Ein_grid(rxn_data, E_bins, nuc_grid, awr, kT, cutoff, &
                               thresh, Ein_el, Ein_inel)
      type(ScattData), target, intent(in) :: rxn_data(:) ! Reaction data
      real(8), intent(in)                 :: E_bins(:)   ! Group structure
      real(8), allocatable, intent(in)    :: nuc_grid(:) ! Nuclidic Ein points
      real(8), intent(in)                 :: kT          ! Library Temperature
      real(8), intent(in)                 :: awr         ! Atomic weight ratio
      real(8), intent(in)                 :: cutoff      ! Free Gas cutoff E
      real(8), intent(in)                 :: thresh      ! Inelastic threshold
      real(8), allocatable, intent(inout) :: Ein_el(:)   ! Incoming Elastic Energy Grid
      real(8), allocatable, intent(inout) :: Ein_inel(:) ! Incoming Inelastic Energy Grid

      integer :: iEmax, iEthresh

      ! Create energy grid to use after limiting nuc_grid to
      ! the maximum energy in E_bins
      if (E_bins(size(E_bins)) >= nuc_grid(size(nuc_grid))) then
        iEmax = size(nuc_grid)
      else
        iEmax = binary_search(nuc_grid, size(nuc_grid), E_bins(size(E_bins)))
      end if

      call merge(nuc_grid(1:iEMax), E_bins, Ein_el)

      ! Add in incoming energies from all reaction channel distributions to
      ! the Ein grid
      call combine_Eins(rxn_data, Ein_el)

      ! Now set inelastic and elastic to be the same (at least for Ein >  thresh)
      iEthresh = binary_search(Ein_el, size(Ein_el), thresh)
      allocate(Ein_inel(1 + size(Ein_el(iEthresh:))))
      Ein_inel = Ein_el(1)
      Ein_inel(2:) = Ein_el(iEthresh:)

      ! Add points to aid in interpolation if elastic scattering points
      call add_elastic_Eins(awr, kT, cutoff, E_bins, Ein_el)

      ! Expand the Incoming energy grid (E_grid) to include points which will
      ! improve interpolation of inelastic level scatter results.
      ! These results are kind of like stair functions but with
      ! near-linear (depends on f(mu)) ramps inbetween each `step'.  For now
      ! we will combat this by putting EXTEND_PTS per current point above the
      ! threshold for inelastic level scatter to begin
      call add_inelastic_Eins(thresh, Ein_inel)

      ! Finally add in one point above energy_bins to give MC code something to
      ! interpolate to if Ein==E_bins(size(E_bins))
      call add_one_more_point(Ein_el)
      call add_one_more_point(Ein_inel)

    end subroutine create_Ein_grid

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
      integer :: i_rxn, iEmax
      real(8) :: max_grp, min_grp

      allocate(new_grid(1))
      new_grid = Ein(1)
      do i_rxn = 1, size(rxn_data)
        mySD => rxn_data(i_rxn)
        if (mySD % is_init) then
          min_grp = mySD % E_bins(1)
          max_grp = mySD % E_bins(size(mySD % E_bins))
          ! Find maximum point that is within our group boundaries before merging
          if (min_grp >= mySD % E_grid(size(mySD % E_grid))) then
            cycle
          else if (max_grp <= mySD % E_grid(1)) then
            cycle
          else
            iEmax = binary_search(mySD % E_grid, size(mySD % E_grid), max_grp)
          end if
          ! Now combine with the rest
          call merge(mySD % E_grid(1: iEmax), new_grid, tmp_grid)
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

    end subroutine combine_Eins

!===============================================================================
! ADD_ELASTIC_EINS adds in incoming energy points to the Ein grid so that the
! continuous scattering distributions can be more accurately reconstructed by
! linear interpolation.  This specific routine adds in points which better
! characterize the behavior as the elastic outgoing energies pass from one
! outgoing group to another.  A total of EXTEND_PTS are added for each of these
! transitions.
!===============================================================================

    subroutine add_elastic_Eins(awr, kT, cutoff, E_bins, Ein)
      real(8), intent(in)                 :: awr       ! Atomic Weight Ratio
      real(8), intent(in)                 :: kT        ! Library Temperature
      real(8), intent(in)                 :: cutoff    ! Free-Gas Cutoff Energy
      real(8), intent(in)                 :: E_bins(:) ! Energy groups
      real(8), allocatable, intent(inout) :: Ein(:)    ! Incoming Energy Grid

      real(8), allocatable :: new_pts(:)
      real(8), allocatable :: old_grid(:)
      integer              :: i, g        ! EXTEND_PTS and Group loop indices
      real(8)              :: alpha
      integer              :: num_pts
      real(8)              :: Ehi, Elo
      real(8)              :: newE       ! Current new incoing energy point
      real(8)              :: dElo, dEhi ! interval between each Ein point
      real(8)              :: lo_shift

      allocate(old_grid(size(Ein)))
      old_grid = Ein

      alpha = ((awr - ONE) / (awr + ONE))**2

      allocate(new_pts(EXTEND_PTS * size(E_bins)))
      new_pts = ZERO

      ! lo_shift is just the delta-Ein below the group boundary to add points
      ! for.
      lo_shift = TWO * kT * (awr + ONE) / awr

      num_pts = 0
      ! Add points for upscattering up until we get to freegas cutoff
      ! (at which point there will no longer be upscatter)
      if (cutoff /= ZERO) then
        do g = 1, size(E_bins) - 1
          Ehi = E_bins(g + 1)
          Elo = E_bins(g)
          if (Ehi <= cutoff) then
            dElo = log(Ehi / (Ehi - lo_shift)  ) / real(EXTEND_PTS, 8)
            do i = -EXTEND_PTS, -1
              newE = Ehi * exp(real(i, 8) * dElo)
              if (newE >= Elo) then
                num_pts = num_pts + 1
                new_pts(num_pts) = newE
              else
                cycle
              end if
            end do
          else if (Elo < cutoff) then
            ! Do something similar to theabove, but adding points before the
            ! freegas cutoff instead of at the group boundary
            Ehi = cutoff
            dElo = log(Ehi / (Ehi - lo_shift)) / real(EXTEND_PTS, 8)
            do i = -EXTEND_PTS, -1
              newE = Ehi * exp(real(i, 8) * dElo)
              if (newE > Elo) then
                num_pts = num_pts + 1
                new_pts(num_pts) = newE
              else
                cycle
              end if
            end do
          end if
        end do

        call merge(new_pts(1: num_pts), old_grid, Ein)
        deallocate(new_pts)
        deallocate(old_grid)
        allocate(old_grid(size(Ein)))
        old_grid = Ein
        allocate(new_pts(EXTEND_PTS * size(E_bins) - 1))
        new_pts = ZERO
        num_pts = 0
      end if

      ! Now add points for the downscatter that occurs as Ein approaches
      ! a group boundary from higher energies

      ! 1/alpha is what you get by saying maximum point is Eg/alpha on log scale
      ! multiplying dEhi by two essentially doubles the range we want to apply over.
      dEhi = 7.0_8 * log(ONE / alpha) / real(EXTEND_PTS, 8)

      do g = 1, size(E_bins) - 1
        if (E_bins(g) == ZERO) then
          cycle
        end if

        Ehi = E_bins(g + 1)

        do i = 1, EXTEND_PTS - 1
          newE = E_bins(g) * exp(real(i, 8) * dEhi)
          if (newE < Ehi) then
            num_pts = num_pts + 1
            new_pts(num_pts) = newE
          else
            exit
          end if
        end do
      end do

      call merge(new_pts(1: num_pts), old_grid, Ein)
      deallocate(new_pts)
      deallocate(old_grid)
    end subroutine add_elastic_Eins

!===============================================================================
! ADD_ONE_MORE_POINT adds ... one more point. To the top, to provide MC code with
! something to interpolate to if particle is exactly the top energy boundary.
!===============================================================================

    subroutine add_one_more_point(Ein)
      real(8), allocatable, intent(inout) :: Ein(:)      ! Incoming Energy Grid

      real(8), allocatable :: temp_grid(:)
      real(8)              :: Ehi
      integer              :: NEin

      Ehi = Ein(size(Ein))
      NEin = size(Ein)

      allocate(temp_grid(NEin + 1))
      temp_grid(1: NEin) = Ein(1:NEin)
      temp_grid(NEin + 1) = Ehi * (ONE + 1.0E-3)

      deallocate(Ein)
      allocate(Ein(size(temp_grid)))
      Ein = temp_grid
      deallocate(temp_grid)

    end subroutine add_one_more_point

!===============================================================================
! ADD_INELASTIC_EINS increases the number of points present above a threshold (min)
! The number of points to increase is EXTEND_PTS for every point on Ein.
!===============================================================================

    subroutine add_inelastic_Eins(min, a)
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
        dE = log(a(i + 1) / a(i)) / real(EXTEND_PTS,8)
        do l = 0, EXTEND_PTS - 1
          temp(j + l) = a(i) * exp(real(l,8) * dE)
        end do
        j = j + EXTEND_PTS
      end do
      temp(size(temp)) = a(size(a))

      deallocate(a)
      allocate(a(size(temp)))
      a = temp
      deallocate(temp)

    end subroutine add_inelastic_eins

!===============================================================================
! CALC_SCATTSAB Calculates the group-to-group transfer matrices and scattering
! moments for the current thermal scattering table
!===============================================================================

    subroutine calc_scattsab(sab, energy_bins, scatt_type, order, &
                             scatt_mat, mu_bins, E_grid)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: energy_bins(:) ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), allocatable, intent(inout)   :: scatt_mat(:,:,:) ! Unionized Scattering Matrices
      integer, intent(in)                   :: mu_bins        ! Number of angular points
                                                              ! to use during f_{n,MT} conversion
      real(8), allocatable, intent(in)      :: E_grid(:)      ! Incoming Energy Grid

      real(8), allocatable :: mu_out(:)      ! Tabular output mu grid
      real(8), allocatable :: sab_int_el(:,:,:) ! Integrated SAB elastic data [L, G, NEin]
      real(8), allocatable :: sab_int_inel(:,:,:) ! Integrated SAB inelastic data [L, G, NEin]
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
                              sab_int_el)
        call integrate_sab_inel(sab, E_grid, energy_bins, scatt_type, order, &
                                sab_int_inel)

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
      call combine_sab_grid(sab_int_el, sab_int_inel, scatt_mat)

      ! Clear the space so its ready next time
      deallocate(sab_int_el)
      deallocate(sab_int_inel)
    end subroutine calc_scattsab

!===============================================================================
! CALC_ELASTIC_GRID Combines all the scattering data points on to one single
! energy grid for elastic scattering
!===============================================================================

    subroutine calc_elastic_grid(nuc, mu_out, rxn_data, Ein, order, E_bins, &
                                 scatt_mat)
      type(Nuclide), pointer, intent(in)   :: nuc   ! The nuclide of interest
      real(8), intent(inout)               :: mu_out(:) ! The tabular output mu grid
      type(ScattData), intent(inout), target :: rxn_data(:) ! The converted distros
      real(8), allocatable, intent(in)     :: Ein(:) ! Elastic Ein grid
      integer, intent(in)                  :: order     ! Angular order
      real(8), intent(in)                  :: E_bins(:) ! Energy groups
      real(8), allocatable, intent(out)    :: scatt_mat(:,:,:) ! Output scatt matrix

      integer :: iE, NE              ! Ein counter, # Ein
      integer :: irxn, Nrxn          ! reaction counter, # Reactions
      integer :: groups              ! # Groups
      type(ScattData), pointer, SAVE :: mySD => NULL() ! Current working ScattData object
      real(8) :: norm_tot            ! Sum of all normalization consts
      integer :: iE_print            ! iE range to print status of
      integer :: iE_pct, last_iE_pct ! Current and previous pct complete
      integer :: tid                 ! Thread id

      groups = size(E_bins) - 1
      NE = size(Ein)
      Nrxn = size(rxn_data)
      iE_print = NE / 20
      last_iE_pct = -1

      ! Allocate the scatt_mat according to the groups, order and number of E pts
      allocate(scatt_mat(order, groups, NE))

      ! Step through each Ein and reactions and sum the scattering distros @ Ein
!$omp parallel do schedule(dynamic,100) num_threads(omp_threads) &
!$omp default(shared), private(iE,mySD,norm_tot,irxn)
      do iE = 1, NE
#ifdef OPENMP
        tid = OMP_GET_THREAD_NUM()
#else
        tid = 0
#endif

        if (.not. mpi_enabled .and. tid == 0) then
          if (iE_print > 0) then
            if ((mod(iE, iE_print) == 1) .or. (iE == NE)) then
              iE_pct = 100 * iE / NE
              if (iE_pct /= last_iE_pct) then
                message = "    Elastic Evaluation " // &
                  trim(to_str(100 * iE / NE)) // "% Complete"
                call write_message(7)
              end if
              last_iE_pct = iE_pct
            end if
          end if
        end if
        if (Ein(iE) <= E_bins(size(E_bins))) then
          scatt_mat(:, :, iE) = ZERO
          norm_tot = ZERO
          do irxn = 1, Nrxn
            mySD => rxn_data(irxn)
            ! If we do not have a scatter reaction, don't score it.
            if (.not. mySD % is_init) cycle
            if (mySD % rxn % MT /= ELASTIC) cycle

            ! Add the scattering distribution to the union scattering grid
            scatt_mat(:, :, iE) = mySD % interp_distro(mu_out, nuc, Ein(iE), norm_tot)
          end do

        else
          ! This step is taken so that interpolation works OK if the MC code
          ! has a particle with an energy == the top energy group value.
          ! With this step, it has something to interpolate to, and that
          ! value is the same as the Ein value, which will be more accurate.
          scatt_mat(:, :, iE) = scatt_mat(:, :, iE - 1)

        end if
      end do
!$omp end parallel do

    end subroutine calc_elastic_grid

!===============================================================================
! CALC_INELASTIC_GRID Combines all the scattering data points on to one single
! energy grid for inelastic scattering
!===============================================================================

    subroutine calc_inelastic_grid(nuc, mu_out, rxn_data, Ein, order, &
                                   E_bins, nuscatt, scatt_mat, nuscatt_mat, norm_tot)
      type(Nuclide), pointer, intent(in)   :: nuc   ! The nuclide of interest
      real(8), intent(inout)               :: mu_out(:) ! The tabular output mu grid
      type(ScattData), intent(inout), target :: rxn_data(:) ! The converted distros
      real(8), allocatable, intent(in)     :: Ein(:)    ! Inelastic Ein grid
      integer, intent(in)                  :: order     ! Angular order
      real(8), intent(in)                  :: E_bins(:) ! Energy groups
      logical, intent(in)                  :: nuscatt   ! Include nuscatter?
      real(8), allocatable, intent(out)    :: scatt_mat(:,:,:) ! Output inelastic matrix
      real(8), allocatable, intent(out)    :: nuscatt_mat(:,:,:) ! Output nu-inelastic matrix
      real(8), allocatable, intent(out)    :: norm_tot(:)

      integer :: iE, NE              ! Ein counter, # Ein
      integer :: irxn, Nrxn          ! reaction counter, # Reactions
      integer :: groups              ! # Groups
      type(ScattData), pointer, SAVE :: mySD => NULL() ! Current working ScattData object
      integer :: iE_print            ! iE range to print status of
      integer :: iE_pct, last_iE_pct ! Current and previous pct complete
      real(8), allocatable :: temp_scatt(:,:) ! calculated scattering matrix
      integer :: tid                 ! Thread id

      groups = size(E_bins) - 1
      NE = size(Ein)
      Nrxn = size(rxn_data)
      iE_print = NE / 20
      last_iE_pct = -1

      ! Allocate the scatt_mat according to the groups, order and number of E pts
      allocate(scatt_mat(order, groups, NE))
      allocate(norm_tot(NE))

      if (nuscatt) then
        allocate(nuscatt_mat(order, groups, NE))
        nuscatt_mat = ZERO
      end if

      allocate(temp_scatt(order, groups))

      ! Step through each Ein and reactions and sum the scattering distros @ Ein
!$omp parallel do schedule(dynamic,100) num_threads(omp_threads) &
!$omp default(shared), private(iE,mySD,irxn,temp_scatt)
      do iE = 1, NE
#ifdef OPENMP
        tid = OMP_GET_THREAD_NUM()
#else
        tid = 0
#endif

        if (.not. mpi_enabled .and. tid == 0) then
          if (iE_print > 0) then
            if ((mod(iE, iE_print) == 1) .or. (iE == NE)) then
              iE_pct = 100 * iE / NE
              if (iE_pct /= last_iE_pct) then
                message = "    Inelastic Evaluation " // &
                  trim(to_str(100 * iE / NE)) // "% Complete"
                call write_message(7)
              end if
              last_iE_pct = iE_pct
            end if
          end if
        end if
        if (Ein(iE) <= E_bins(size(E_bins))) then
          scatt_mat(:, :, iE) = ZERO
          norm_tot(iE) = ZERO
          do irxn = 1, Nrxn
            mySD => rxn_data(irxn)
            ! If we do not have a scatter reaction, don't score it.
            if (.not. mySD % is_init) cycle
            if (mySD % rxn % MT == ELASTIC) cycle

            ! Add the scattering distribution to the union scattering grid
            temp_scatt = mySD % interp_distro(mu_out, nuc, Ein(iE), norm_tot(iE))
            scatt_mat(:, :, iE) = scatt_mat(:, :, iE) + temp_scatt
            if (nuscatt) then
              nuscatt_mat(:, :, iE) = nuscatt_mat(:, :, iE) + &
                real(mySD % rxn % multiplicity, 8) * temp_scatt
            end if
          end do

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

    end subroutine calc_inelastic_grid

!===============================================================================
! PRINT_SCATT prints the scattering data to the specified output file
! in the specified format.
!===============================================================================

  subroutine print_scatt(lib_format, grp_index_el, grp_index_inel, Ein_el, &
                         Ein_inel, print_tol, thin_tol, el_mat, inel_mat, &
                         nuinel_mat)
    integer,              intent(in) :: lib_format   ! Library output type
    integer,              intent(in) :: grp_index_el(:) ! energy_group locations in Ein_el
    integer,              intent(in) :: grp_index_inel(:) ! energy_group locations in Ein_inel
    real(8), allocatable, intent(in) :: Ein_el(:)    ! Elastic Ein grid
    real(8), allocatable, intent(in) :: Ein_inel(:)  ! Inelastic Ein grid
    real(8),              intent(in) :: print_tol    ! Minimum grp-to-grp prob'y
                                                     ! to print
    real(8),              intent(in) :: thin_tol     ! Minimum grp-to-grp prob'y
                                                     ! to print
    real(8), allocatable, intent(in) :: el_mat(:,:,:)    ! Elastic data to print
    real(8), allocatable, intent(in) :: inel_mat(:,:,:)  ! Inelastic data to print
    real(8), allocatable, optional, intent(in) :: nuinel_mat(:,:,:) ! Nu-Inel data to print

    if (present(nuinel_mat)) then
      if (lib_format == ASCII) then
        call print_scatt_ascii(grp_index_el, grp_index_inel, Ein_el, &
                               Ein_inel, print_tol, el_mat, inel_mat, &
                               nuinel_mat)
      else if (lib_format == BINARY) then
        call print_scatt_bin(grp_index_el, grp_index_inel, Ein_el, &
                             Ein_inel, print_tol, el_mat, inel_mat, &
                             nuinel_mat)
      else if (lib_format == HUMAN) then
        call print_scatt_human(grp_index_el, grp_index_inel, Ein_el, &
                               Ein_inel, print_tol, el_mat, inel_mat, &
                               nuinel_mat)
#ifdef HDF5
      else if (lib_format == H5) then
        call print_scatt_hdf5(grp_index_el, grp_index_inel, Ein_el, &
                              Ein_inel, print_tol, el_mat, inel_mat, &
                              nuinel_mat)
#endif
      end if
    else
      if (lib_format == ASCII) then
        call print_scatt_ascii(grp_index_el, grp_index_inel, Ein_el, &
                               Ein_inel, print_tol, el_mat, inel_mat)
      else if (lib_format == BINARY) then
        call print_scatt_bin(grp_index_el, grp_index_inel, Ein_el, &
                             Ein_inel, print_tol, el_mat, inel_mat)
      else if (lib_format == HUMAN) then
        call print_scatt_human(grp_index_el, grp_index_inel, Ein_el, &
                               Ein_inel, print_tol, el_mat, inel_mat)
#ifdef HDF5
      else if (lib_format == H5) then
        call print_scatt_hdf5(grp_index_el, grp_index_inel, Ein_el, &
                              Ein_inel, print_tol, el_mat, inel_mat)
#endif
      end if
    end if

  end subroutine print_scatt

!===============================================================================
! PRINT_SCATT_ASCII prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_scatt_ascii(grp_index_el, grp_index_inel, Ein_el, &
                               Ein_inel, print_tol, el_mat, inel_mat, &
                               nuinel_mat)
    integer,              intent(in) :: grp_index_el(:) ! energy_group locations in Ein_el
    integer,              intent(in) :: grp_index_inel(:) ! energy_group locations in Ein_inel
    real(8), allocatable, intent(in) :: Ein_el(:)    ! Elastic Ein grid
    real(8), allocatable, intent(in) :: Ein_inel(:)  ! Inelastic Ein grid
    real(8),              intent(in) :: print_tol    ! Minimum grp-to-grp prob'y
                                                     ! to print
    real(8), allocatable, intent(in) :: el_mat(:,:,:)    ! Elastic data to print
    real(8), allocatable, intent(in) :: inel_mat(:,:,:)  ! Inelastic data to print
    real(8), allocatable, optional, intent(in) :: nuinel_mat(:,:,:) ! Nu-Inel data to print

    integer :: gmin, gmax, iE

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line:
    ! FOR ELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! FOR INELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! Repeat data for nuscatt

    ! Begin writing:
    ! ELASTIC
    ! # energy points
    write(UNIT_NUC,'(I20)') size(Ein_el)

    ! <incoming energy array>
    call print_ascii_array(Ein_el, UNIT_NUC)

    ! # Group Indices
    call print_ascii_integer_array(grp_index_el, UNIT_NUC)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(Ein_el)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(el_mat, dim = 2)
        if (el_mat(1, gmin, iE) > print_tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(el_mat, dim = 2), 1, -1
        if (el_mat(1, gmax, iE) > print_tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC, '(I20,I20)') 0,0
      else
        write(UNIT_NUC, '(I20,I20)') gmin,gmax
        call print_ascii_array(reshape(el_mat(:, gmin : gmax, iE), (/ &
          size(el_mat, dim=1) * (gmax - gmin + 1)/)), UNIT_NUC)
      end if
    end do

    ! INELASTIC
    ! # energy points
    write(UNIT_NUC,'(I20)') size(Ein_inel)

    ! <incoming energy array>
    call print_ascii_array(Ein_inel, UNIT_NUC)

    ! # Group Indices
    call print_ascii_integer_array(grp_index_inel, UNIT_NUC)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(Ein_inel)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(inel_mat, dim = 2)
        if (inel_mat(1, gmin, iE) > print_tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(inel_mat, dim = 2), 1, -1
        if (inel_mat(1, gmax, iE) > print_tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC, '(I20,I20)') 0,0
      else
        write(UNIT_NUC, '(I20,I20)') gmin,gmax
        call print_ascii_array(reshape(inel_mat(:, gmin : gmax, iE), (/ &
          size(inel_mat, dim=1) * (gmax - gmin + 1)/)), UNIT_NUC)
      end if
    end do

    if (present(nuinel_mat)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(Ein_inel)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nuinel_mat, dim = 2)
          if (nuinel_mat(1, gmin, iE) > print_tol) exit
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nuinel_mat, dim = 2), 1, -1
          if (nuinel_mat(1, gmax, iE) > print_tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          write(UNIT_NUC, '(I20,I20)') 0,0
        else
          write(UNIT_NUC, '(I20,I20)') gmin,gmax
          call print_ascii_array(reshape(nuinel_mat(:, gmin : gmax, iE), (/ &
            size(nuinel_mat, dim=1) * (gmax - gmin + 1)/)), UNIT_NUC)
        end if
      end do
    end if

  end subroutine print_scatt_ascii

!===============================================================================
! PRINT_SCATT_HUMAN prints the scattering data to the specified output file
! in an ASCII format.
!===============================================================================

  subroutine print_scatt_human(grp_index_el, grp_index_inel, Ein_el, &
                               Ein_inel, print_tol, el_mat, inel_mat, &
                               nuinel_mat)
    integer,              intent(in) :: grp_index_el(:) ! energy_group locations in Ein_el
    integer,              intent(in) :: grp_index_inel(:) ! energy_group locations in Ein_inel
    real(8), allocatable, intent(in) :: Ein_el(:)    ! Elastic Ein grid
    real(8), allocatable, intent(in) :: Ein_inel(:)  ! Inelastic Ein grid
    real(8),              intent(in) :: print_tol    ! Minimum grp-to-grp prob'y
                                                     ! to print
    real(8), allocatable, intent(in) :: el_mat(:,:,:)    ! Elastic data to print
    real(8), allocatable, intent(in) :: inel_mat(:,:,:)  ! Inelastic data to print
    real(8), allocatable, optional, intent(in) :: nuinel_mat(:,:,:) ! Nu-Inel data to print

    integer :: g, gmin, gmax, iE

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format with at max 4 entries per line:
    ! FOR ELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! FOR INELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! Repeat data for nuscatt

    ! Begin writing:
    ! ELASTIC
    ! # energy points
    write(UNIT_NUC,'(I20)') size(Ein_el)

    ! <incoming energy array>
    call print_ascii_array(Ein_el, UNIT_NUC)

    ! # Group Indices
    call print_ascii_integer_array(grp_index_el, UNIT_NUC)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(Ein_el)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(el_mat, dim = 2)
        if (el_mat(1, gmin, iE) > print_tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(el_mat, dim = 2), 1, -1
        if (el_mat(1, gmax, iE) > print_tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',Ein_el(iE), &
          '   gmin = ', 0, '   gmax = ', 0
      else
        write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',Ein_el(iE), &
          '   gmin = ', gmin, '   gmax = ', gmax
        do g = gmin, gmax
          write(UNIT_NUC,'(A,I5)') 'outgoing group = ', g
          call print_ascii_array(el_mat(:, g, iE), UNIT_NUC)
        end do
      end if
    end do

    ! INELASTIC
    ! # energy points
    write(UNIT_NUC,'(I20)') size(Ein_inel)

    ! <incoming energy array>
    call print_ascii_array(Ein_inel, UNIT_NUC)

    ! # Group Indices
    call print_ascii_integer_array(grp_index_inel, UNIT_NUC)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(Ein_inel)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(inel_mat, dim = 2)
        if (inel_mat(1, gmin, iE) > print_tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(inel_mat, dim = 2), 1, -1
        if (inel_mat(1, gmax, iE) > print_tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',Ein_inel(iE), &
          '   gmin = ', 0, '   gmax = ', 0
      else
        write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',Ein_inel(iE), &
          '   gmin = ', gmin, '   gmax = ', gmax
        do g = gmin, gmax
          write(UNIT_NUC,'(A,I5)') 'outgoing group = ', g
          call print_ascii_array(inel_mat(:, g, iE), UNIT_NUC)
        end do
      end if
    end do

    if (present(nuinel_mat)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(Ein_inel)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nuinel_mat, dim = 2)
          if (nuinel_mat(1, gmin, iE) > print_tol) exit
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nuinel_mat, dim = 2), 1, -1
          if (nuinel_mat(1, gmax, iE) > print_tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',Ein_inel(iE), &
            '   gmin = ', 0, '   gmax = ', 0
        else
          write(UNIT_NUC, '(A,1PE20.12,A,I5,A,I5)') 'Ein = ',Ein_inel(iE), &
            '   gmin = ', gmin, '   gmax = ', gmax
          do g = gmin, gmax
            write(UNIT_NUC,'(A,I5)') 'outgoing group = ', g
            call print_ascii_array(nuinel_mat(:, g, iE), UNIT_NUC)
          end do
        end if
      end do
    end if

  end subroutine print_scatt_human

!===============================================================================
! PRINT_SCATT_BIN prints the scattering data to the specified output file
! in a native Fortran stream format.
!===============================================================================

  subroutine print_scatt_bin(grp_index_el, grp_index_inel, Ein_el, &
                             Ein_inel, print_tol, el_mat, inel_mat, &
                             nuinel_mat)
    integer,              intent(in) :: grp_index_el(:) ! energy_group locations in Ein_el
    integer,              intent(in) :: grp_index_inel(:) ! energy_group locations in Ein_inel
    real(8), allocatable, intent(in) :: Ein_el(:)    ! Elastic Ein grid
    real(8), allocatable, intent(in) :: Ein_inel(:)  ! Inelastic Ein grid
    real(8),              intent(in) :: print_tol    ! Minimum grp-to-grp prob'y
                                                     ! to print
    real(8), allocatable, intent(in) :: el_mat(:,:,:)    ! Elastic data to print
    real(8), allocatable, intent(in) :: inel_mat(:,:,:)  ! Inelastic data to print
    real(8), allocatable, optional, intent(in) :: nuinel_mat(:,:,:) ! Nu-Inel data to print

    integer :: g, gmin, gmax, iE

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format:
    ! FOR ELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! FOR INELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! Repeat data for nuscatt

    ! Begin writing:
    ! ELASTIC
    ! # energy points
    write(UNIT_NUC) size(Ein_el)

    ! <incoming energy array>
    write(UNIT_NUC) Ein_el

    ! Group Indices
    write(UNIT_NUC) grp_index_el

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(Ein_el)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(el_mat, dim = 2)
        if (el_mat(1, gmin, iE) > print_tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(el_mat, dim = 2), 1, -1
        if (el_mat(1, gmax, iE) > print_tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC) 0, 0
      else
        write(UNIT_NUC) gmin, gmax
        do g = gmin, gmax
          write(UNIT_NUC) el_mat(:, g, iE)
        end do
      end if
    end do

    ! INELASTIC
    ! # energy points
    write(UNIT_NUC) size(Ein_inel)

    ! <incoming energy array>
    write(UNIT_NUC) Ein_inel

    ! Group Indices
    write(UNIT_NUC) grp_index_inel

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    do iE = 1, size(Ein_inel)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(inel_mat, dim = 2)
        if (inel_mat(1, gmin, iE) > print_tol) exit
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(inel_mat, dim = 2), 1, -1
        if (inel_mat(1, gmax, iE) > print_tol) exit
      end do
      if (gmin > gmax) then ! we have effectively all zeros
        write(UNIT_NUC) 0, 0
      else
        write(UNIT_NUC) gmin, gmax
        do g = gmin, gmax
          write(UNIT_NUC) inel_mat(:, g, iE)
        end do
      end if
    end do

    if (present(nuinel_mat)) then
      ! < \nu-\Sigma_{s,g',l}(Ein) array as follows for each Ein:
      ! g'_min, g'_max, for g' in g'_min to g'_max: \nu-\Sigma_{s,g',1:L}(Ein)>
      do iE = 1, size(Ein_inel)
        ! find gmin by checking the P0 moment
        do gmin = 1, size(nuinel_mat, dim = 2)
          if (nuinel_mat(1, gmin, iE) > print_tol) exit
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nuinel_mat, dim = 2), 1, -1
          if (nuinel_mat(1, gmax, iE) > print_tol) exit
        end do
        if (gmin > gmax) then ! we have effectively all zeros
          write(UNIT_NUC) 0, 0
        else
          write(UNIT_NUC) gmin, gmax
          do g = gmin, gmax
            write(UNIT_NUC) nuinel_mat(:, g, iE)
          end do
        end if
      end do
    end if

  end subroutine print_scatt_bin

!===============================================================================
! PRINT_SCATT_HDF5 prints the scattering data to the specified output group
! with the HDF5 library.
!===============================================================================
#ifdef HDF5
  subroutine print_scatt_hdf5(grp_index_el, grp_index_inel, Ein_el, &
                              Ein_inel, print_tol, el_mat, inel_mat, &
                              nuinel_mat)
    integer,              intent(in) :: grp_index_el(:) ! energy_group locations in Ein_el
    integer,              intent(in) :: grp_index_inel(:) ! energy_group locations in Ein_inel
    real(8), allocatable, intent(in) :: Ein_el(:)    ! Elastic Ein grid
    real(8), allocatable, intent(in) :: Ein_inel(:)  ! Inelastic Ein grid
    real(8),              intent(in) :: print_tol    ! Minimum grp-to-grp prob'y
                                                     ! to print
    real(8), allocatable, intent(in) :: el_mat(:,:,:)    ! Elastic data to print
    real(8), allocatable, intent(in) :: inel_mat(:,:,:)  ! Inelastic data to print
    real(8), allocatable, optional, intent(in) :: nuinel_mat(:,:,:) ! Nu-Inel data to print

    integer :: g, gmin, gmax, iE

    integer(HID_T) :: orig_group, temp_group, scatt_group
    character(MAX_FILE_LEN) :: group_name, iE_name, nuc_name
    integer :: period_loc

   ! Create a new hdf5 group for the scatter.
    group_name = "/scatt"
    call hdf5_open_group(h5_file, group_name, temp_group)

    ! Assumes that the file and header information is already printed
    ! (including # of groups and bins, and thinning tolerance)
    ! Will follow this format:
    ! FOR ELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! FOR INELASTIC:
    ! NEin
    ! Ein[:]
    ! Group Indices
    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! Repeat data for nuscatt

    ! Begin writing:

    ! Group Indices
    call hdf5_write_integer_1Darray(temp_group, 'Group_Index', grp_index, &
                                    size(grp_index))

    ! # energy points !!! Maybe not necessary with HDF5, can I get the size
    ! easily???
    call hdf5_write_integer(temp_group, 'NEin', size(E_grid))

    ! <incoming energy array>
    call hdf5_write_double_1Darray(temp_group, 'Ein', E_grid, size(E_grid))

    call hdf5_close_group(temp_group)

    ! < \Sigma_{s,g',l}(Ein) array as follows for each Ein:
    ! g'_min, g'_max, for g' in g'_min to g'_max: \Sigma_{s,g',1:L}(Ein)>
    ! THis would be better served with a data set.
    do iE = 1, size(E_grid)
      iE_name = trim(adjustl(group_name)) // "/iE" // trim(adjustl(to_str(iE)))
      call hdf5_open_group(iE_name)
      ! find gmin by checking the P0 moment
      do gmin = 1, size(data, dim = 2)
        if (data(1, gmin, iE) > print_tol) then
          call hdf5_close_group()
          temp_group = scatt_group
          exit
        end if
      end do
      ! find gmax by checking the P0 moment
      do gmax = size(data, dim = 2), 1, -1
        if (data(1, gmax, iE) > print_tol) exit
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
          if (nudata(1, gmin, iE) > print_tol) then
            call hdf5_close_group()
            temp_group = scatt_group
            exit
          end if
        end do
        ! find gmax by checking the P0 moment
        do gmax = size(nudata, dim = 2), 1, -1
          if (nudata(1, gmax, iE) > print_tol) exit
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
  end subroutine print_scatt_hdf5
#endif

end module scatt
