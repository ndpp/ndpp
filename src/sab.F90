module sab

  use ace_header,  only: SAlphaBeta
  use array_merge, only: merge
  use constants
  use error,       only: fatal_error, warning
  use global
  use legendre
  use search,      only: binary_search
  use string,      only: to_str

  implicit none

  contains

!===============================================================================
! INTEGRATE_SAB_EL performs integration of the thermal scattering elastic
! collisions and places the result in sab_int
!===============================================================================

    subroutine integrate_sab_el(sab, ein_grid, e_bins, scatt_type, order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      real(8) :: sig       ! Micros. x/s
      integer :: iEin      ! incoming energy counter
      integer :: groups    ! shorthand for number of energy groups
      integer :: g, imu, l ! indices for: group, angle, legendre order
      real(8) :: mu        ! cosine of angle of scatter
      real(8) :: Ein       ! Incoming energy
      integer :: isab      ! index on sab Ein grid
      real(8) :: f         ! fraction for interpolation
      real(8) :: wgt       ! weighting of each angle

      ! Initialize data
      sab_int = ZERO
      groups = size(e_bins) - 1

      ! Exit if there are no elastic reactions
      if (sab % threshold_elastic == ZERO) then
        return
      end if

      if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
        wgt = ONE / real(sab % n_elastic_mu, 8)
      end if

      !$omp parallel do schedule(dynamic,20) num_threads(omp_threads) &
      !$omp default(shared),private(iEin, g, imu, l, mu, Ein, isab, f)
      do iEin = 1, size(ein_grid)
        Ein = ein_grid(iEin)
        ! Find the index and interpolation factor
        if (Ein < sab % elastic_e_in(1)) then
          cycle
        else if (Ein >= sab % threshold_elastic) then
          cycle
        else
          isab = binary_search(sab % elastic_e_in, sab % n_elastic_e_in, Ein)
          f = (Ein - sab % elastic_e_in(isab)) / &
            (sab % elastic_e_in(isab + 1) - sab % elastic_e_in(isab))
        end if

        ! set outgoing group, g.  No outgoing group integration necessary
        ! because Eout does not change for elastic; Ein == Eout
        if (Ein < e_bins(1)) then
          cycle
        else if (Ein > e_bins(size(e_bins))) then
          cycle
        else
          g = binary_search(e_bins, size(e_bins), Ein)
        end if

        ! Get x/s for normalizing
        if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          sig = sab % elastic_P(isab) / Ein
        else if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
          sig = (ONE - f) * sab % elastic_P(isab) + &
            f * sab % elastic_P(isab + 1)
        end if

        if (sab % n_elastic_mu == 0) then
          ! Perform calculations for coherent elastic (u = 1-E_{bragg}/E)
          mu = ONE - sab % elastic_e_in(isab) / Ein
          do l = 1, order + 1
            sab_int(l, g, iEin) = sab_int(l, g, iEin) + calc_pn(l - 1, mu)
          end do
        else if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
          ! Loop over outgoing angles and add legendre of each
          do imu = 1, sab % n_elastic_mu
            ! Find our interpolated mu
            mu = (ONE - f) * sab % elastic_mu(imu, isab) + &
              f * sab % elastic_mu(imu, isab + 1)
            ! And add the discrete point to our integration
            do l = 1, order + 1
              sab_int(l, g, iEin) = sab_int(l, g, iEin) + wgt * &
                calc_pn(l - 1, mu)
            end do
          end do
        else
          ! pass
        end if
        sab_int(:, :, iEin) = sig * sab_int(:, :, iEin)
      end do

    end subroutine integrate_sab_el

!===============================================================================
! INTEGRATE_SAB_INEL performs integration of the thermal scattering inelastic
! collisions and places the result in sab_int by calling INTEGRATE_SAB_INEL_DISC
! or INTEGRATE_SAB_INEL_CONT when necessary
!===============================================================================

    subroutine integrate_sab_inel(sab, ein_grid, e_bins, scatt_type, order, &
                                  sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      if ((sab % secondary_mode == SAB_SECONDARY_EQUAL) .or. &
          (sab % secondary_mode == SAB_SECONDARY_SKEWED)) then
          call integrate_sab_inel_disc(sab, ein_grid, e_bins, scatt_type, &
                                       order, sab_int)
      else if (sab % secondary_mode == SAB_SECONDARY_CONT) then
          call integrate_sab_inel_cont(sab, ein_grid, e_bins, scatt_type, &
                                       order, sab_int)
      end if
    end subroutine integrate_sab_inel

!===============================================================================
! INTEGRATE_SAB_INEL_DISC performs integration of the thermal scattering
! inelastic collisions and places the result in sab_int for discrete inelastic
! (skewed or equi-probable) secondary energy distributions
!===============================================================================

    subroutine integrate_sab_inel_disc(sab, ein_grid, e_bins, scatt_type, &
                                       order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      real(8) :: sig       ! Micros. x/s
      integer :: iEin      ! incoming energy counter
      integer :: iEout     ! outgoing energy counter
      integer :: groups    ! shorthand for number of energy groups
      integer :: g, imu, l ! indices for: group, angle, legendre order
      real(8) :: mu        ! cosine of angle of scatter
      real(8) :: Ein, Eout ! Incoming & outgoing energies
      integer :: isab      ! index on sab Ein grid
      real(8) :: f         ! fraction for interpolation
      real(8), allocatable :: wgt(:) ! weighting (based on skewed or equal mode)

      ! Initialize data
      sab_int = ZERO
      groups = size(e_bins) - 1

      ! First lets set up our weighting
      if (sab % secondary_mode == SAB_SECONDARY_EQUAL) then
        allocate(wgt(sab % n_inelastic_e_out))
        wgt = ONE / (real(sab % n_inelastic_e_out, 8) * &
                     real(sab % n_inelastic_mu, 8))
      else if (sab % secondary_mode == SAB_SECONDARY_SKEWED) then
        if (sab % n_inelastic_e_out > 4) then
          allocate(wgt(sab % n_inelastic_e_out))
          ! 0.1, 0.4, equally-likely, 0.4, 0.1 (all normalized to 1)
          wgt(1) = 0.1_8
          wgt(2) = 0.4_8
          wgt(3: sab % n_inelastic_e_out - 2) = ONE
          wgt(sab % n_inelastic_e_out - 1) = 0.4_8
          wgt(sab % n_inelastic_e_out) = 0.1_8
          wgt = wgt / (sum(wgt) * real(sab % n_inelastic_mu, 8))
        else
          ! Leave this in until I know what to do (this should be the continuous)
          message = "Number of Inelastic Outgoing Energies Less Than 4,&
                    & but Skewed Weighting Requested by Data!"
          call fatal_error()
        end if
      end if

      !$omp parallel do schedule(dynamic,20) num_threads(omp_threads) &
      !$omp default(shared),private(iEin, iEout, g, imu, l, mu, Ein, Eout, isab, f)
      do iEin = 1, size(ein_grid)
        Ein = ein_grid(iEin)
        ! Find the index and interpolation factor
        if (Ein < sab % inelastic_e_in(1)) then
          isab = 1
          f = ZERO
        else if (Ein > sab % threshold_inelastic) then
          cycle
        else if (Ein == sab % threshold_inelastic) then
          isab = size(sab % inelastic_e_in) - 1
          f = ONE
        else
          isab = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, Ein)
          f = (Ein - sab % inelastic_e_in(isab)) / &
            (sab % inelastic_e_in(isab + 1) - sab % inelastic_e_in(isab))
        end if

        ! Get x/s for normalizing
        sig = (ONE - f) * sab % inelastic_sigma(isab) + &
          f * sab % inelastic_sigma(isab + 1)

        ! Integrate over outgoing energy (outer) and outgoing mu (inner)

        ! Loop through each equally likely energy, find its group
        ! and sum legendre moments of discrete mu to it.
        ! Eout and mu need to be interpolated to
        do iEout = 1, sab % n_inelastic_e_out
          ! Interpolate to our Eout
          Eout = (ONE - f) * sab % inelastic_e_out(iEout, isab) + &
            f * sab % inelastic_e_out(iEout, isab + 1)

          ! Get the outgoing group, E_bins(1) should be zero, but thats OK
          if (Eout < e_bins(1)) then
            cycle
          else if (Eout >= e_bins(size(e_bins))) then
            cycle
          else
            g = binary_search(e_bins, size(e_bins), Eout)
          end if

          ! Now integrate all the discrete mu values
          do imu = 1, sab % n_inelastic_mu
            ! Find our interpolated mu
            mu = (ONE - f) * sab % inelastic_mu(imu, iEout, isab) + &
              f * sab % inelastic_mu(imu, iEout, isab + 1)
            ! Probably could write a function to do this using Legendre recursion
            ! avoiding the loop over order. Oh well, this should be pretty cheap
            do l = 1, order + 1
              sab_int(l, g, iEin) = sab_int(l, g, iEin) + &
                calc_pn(l - 1, mu) * wgt(iEout)
            end do
          end do
        end do
        sab_int(:, :, iEin) = sig * sab_int(:, :, iEin)
      end do
    end subroutine integrate_sab_inel_disc

!===============================================================================
! INTEGRATE_SAB_INEL_DISC performs integration of the thermal scattering
! inelastic collisions and places the result in sab_int for discrete inelastic
! (skewed or equi-probable) secondary energy distributions
!===============================================================================

    subroutine integrate_sab_inel_cont(sab, ein_grid, e_bins, scatt_type, &
                                       order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      real(8) :: sig       ! Micros. x/s
      integer :: iEin      ! incoming energy counter
      integer :: groups    ! shorthand for number of energy groups
      integer :: g, imu, l ! indices for: group, angle, legendre order
      real(8) :: mu        ! cosine of angle of scatter
      real(8) :: Ein       ! Incoming Energy
      integer :: isab      ! index on sab Ein grid
      real(8) :: f         ! fraction for interpolation
      real(8), allocatable :: distro(:,:,:) ! E'-mu double integral on S(a,b) Ein grid
      integer :: iE, iE_lo, iE_hi
      real(8) :: f_lo, f_hi   ! Interpolation for lo and hi
      integer :: NEout
      real(8), allocatable :: pdf(:)
      real(8), pointer :: Eout_arr(:), mu_arr(:,:)
      real(8) :: mult

      ! Initialize data
      sab_int = ZERO
      groups = size(e_bins) - 1
      allocate(distro(order + 1, groups, sab % n_inelastic_e_in))
      distro = ZERO

      ! Before we do the integration at the ein_grid points, we will
      ! calculate the energy-angle double integral over the S(a,b) data at
      ! the S(a,b) data's own Ein,Eout grid.  Then when that is done we can
      ! interpolate these results to the ein_grid provided to this routine.

      ! So, start with the energy-angle double integral of the S(a,b) data
      !$omp parallel do schedule(dynamic,20) num_threads(omp_threads) &
      !$omp default(private),shared(sab,groups,order,e_bins,distro,sab_int)
      do iEin = 1, sab % n_inelastic_e_in
        NEout = sab % inelastic_data(iEin) % n_e_out
        allocate(pdf(NEout))
        Eout_arr => sab % inelastic_data(iEin) % e_out
        mu_arr => sab % inelastic_data(iEin) % mu

        ! First we normalize the pdf
        do iE_lo = 1, NEout - 1
          pdf(iE_lo) = sab % inelastic_data(iEin) % e_out_pdf(iE_lo) * &
            (Eout_arr(iE_lo + 1) - Eout_arr(iE_lo))
        end do
        pdf(NEout) = ZERO

        ! Now go through each group and integrate within the energy bounds of
        ! that group.
        do g = 1, groups
          ! Find iE_lo, and add first term to distro
          if (e_bins(g) < Eout_arr(1)) then
            ! We need to skip this lower bound
            iE_lo = 1
          else if (e_bins(g) >= Eout_arr(NEout)) then
            ! In this case, the lower group boundary is above all energies;
            ! this means the group is outside the Eout range, and thus this group
            ! has a zero distribution
            distro(:, g, iEin) = ZERO
            cycle
          else
            ! The group boundary is inbetween 2 Eout pts. Interpolate.
            iE_lo = binary_search(Eout_arr, NEout, e_bins(g))
            f_lo = (e_bins(g) - Eout_arr(iE_lo)) / &
              (Eout_arr(iE_lo + 1) - Eout_arr(iE_lo))
            mult = f_lo * pdf(iE_lo)
            do imu = 1, sab % n_inelastic_mu
              ! Find our interpolated mu
              mu = (ONE - f_lo) * mu_arr(imu, iE_lo) + &
                f_lo * mu_arr(imu, iE_lo + 1)
              do l = 1, order + 1
                distro(l, g, iEin) = distro(l, g, iEin) + &
                  calc_pn(l - 1, mu) * mult
              end do
            end do
            iE_lo = iE_lo + 1
          end if

          ! Find iE_hi and add the last term to distro
          if (e_bins(g + 1) < Eout_arr(1)) then
            ! We can skip this group completely then
            distro(:, g, iEin) = ZERO
            cycle
          else if (e_bins(g + 1) >= Eout_arr(NEout)) then
            ! The upper grp boundary is above Eout, and so its value is zero
            ! We therefore dont need to add its contribution to the integral
            iE_hi = NEout - 1
          else
            ! The group boundary is inbetween 2 Eout pts. Interpolate.
            iE_hi = binary_search(Eout_arr, NEout, e_bins(g + 1))
            f_hi = (e_bins(g + 1) - Eout_arr(iE_hi)) / &
              (Eout_arr(iE_hi + 1) - Eout_arr(iE_hi))
            mult = f_hi * pdf(iE_hi)
            do imu = 1, sab % n_inelastic_mu
              ! Find our interpolated mu
              mu = (ONE - f_hi) * mu_arr(imu, iE_hi) + &
                f_hi * mu_arr(imu, iE_hi + 1)
              do l = 1, order + 1
                distro(l, g, iEin) = distro(l, g, iEin) + &
                  calc_pn(l - 1, mu) * mult
              end do
            end do
            iE_hi = iE_hi - 1
          end if

          ! Now we can do the intermediate points
          do iE = iE_lo, iE_hi
            mult = pdf(iE)
            do imu = 1, sab % n_inelastic_mu
              do l = 1, order + 1
                distro(l, g, iEin) = distro(l, g, iEin) + &
                  calc_pn(l - 1, mu_arr(imu, iE)) * pdf(iE)
              end do
            end do
          end do

          ! Finish the division from trapezoidal integration
          distro(:, g, iEin) = distro(:, g, iEin)  / real(sab % n_inelastic_mu, 8)
        end do
        deallocate(pdf)
      end do

      ! Now interpolate to the ein_grid, while also finding sig
      !$omp parallel do schedule(dynamic,20) num_threads(omp_threads) &
      !$omp default(shared),private(iEin,Ein,isab,f)
      do iEin = 1, size(ein_grid)
        Ein = ein_grid(iEin)
        ! Find the index and interpolation factor
        if (Ein <= sab % inelastic_e_in(1)) then
          isab = 1
          f = ZERO
          ! Get x/s for normalizing
          sig = sab % inelastic_sigma(isab)
          ! Find the sab integral
          sab_int(:,:,iEin) = distro(:,:,isab) * sig
        else if (Ein >= sab % threshold_inelastic) then
          sab_int(:,:,iEin) = ZERO
          cycle
        else
          isab = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, Ein)
          f = (Ein - sab % inelastic_e_in(isab)) / &
            (sab % inelastic_e_in(isab + 1) - sab % inelastic_e_in(isab))
          ! Get x/s for normalizing
          sig = (ONE - f) * sab % inelastic_sigma(isab) + &
            f * sab % inelastic_sigma(isab + 1)
          ! Find the sab integral
          sab_int(:,:,iEin) = ((ONE - f) * distro(:,:,isab) + &
            f * distro(:,:,isab + 1)) * sig
        end if
      end do
    end subroutine integrate_sab_inel_cont

!===============================================================================
! COMBINE_SAB_GRID takes elastic and inelastic integrals and puts them on the
! same grid by weighting by the x/s
!===============================================================================

  subroutine combine_sab_grid(sab_int_el, sab_int_inel, scatt_mat)
    real(8), intent(in) :: sab_int_el(:,:,:)   ! Integrated SAB elastic data [L, G, NEin]
    real(8), intent(in) :: sab_int_inel(:,:,:) ! Integrated SAB inelastic data [L, G, NEin]
    real(8), allocatable, intent(inout) :: scatt_mat(:,:,:) ! Unionized Scattering Matrices

    integer :: iE, NE
    real(8) :: norm_sum

    ! set up our scatt_mat space
    allocate(scatt_mat(size(sab_int_el, 1), size(sab_int_el, 2), &
             size(sab_int_el, 3)))

    NE = size(sab_int_el, dim=3)

    !$omp parallel do schedule(dynamic,20) num_threads(omp_threads) &
    !$omp default(shared),private(iE, norm_sum)
    do iE = 1, NE

      ! Combine the elastic and inelastic data to the same grid and
      ! normalize the results since the PDF weighting doesn't seem
      ! to be guaranteed to be normalized, at least in lwtr.26t
      ! Never-the-less, since Emin is forced to 0 by NDPP, and Emax is assumed
      ! sufficiently high, we know the PDF should normalize to 1,
      scatt_mat(:, :, iE) = (sab_int_el(:, :, iE) + sab_int_inel(:, :, iE))
      norm_sum = sum(scatt_mat(1, :, iE))
      if (norm_sum > ZERO) then
        norm_sum = ONE / norm_sum
        scatt_mat(:, :, iE) = scatt_mat(:, :, iE) * norm_sum
      else
        scatt_mat(:, :, iE) = ZERO
      end if
    end do

    ! Finally, since we want the interpolation between the very last two
    ! points to work out right (the threshold scatt_mat will be zero, which is
    ! not physical), we will set the threshold scatt_mat equal to the scatt_mat
    ! for the iE just before this
    scatt_mat(:, :, NE) = scatt_mat(:, :, NE - 1)

  end subroutine combine_sab_grid

!===============================================================================
! SAB_EGRID Creates the thermal scattering energy grid for NDPP integration
!===============================================================================

    subroutine sab_egrid(sab, energy_bins, Ein)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! SAB Library
      real(8), intent(in)                   :: energy_bins(:) ! Grp bounds
      real(8), allocatable, intent(inout)   :: Ein(:)         ! Resultant Egrid

      real(8), allocatable   :: e_grid_tmp(:) ! Temporary array for bldg E-grid
      integer                :: i_max_ein     ! Index of maximum value we want from Ein
      real(8)                :: max_ein       ! maximum value we want from Ein
      integer                :: num_pts       ! Final number of pts for Ein
      integer                :: i, j, iE      ! Counters for final indexing
      real(8)                :: dE
      integer                :: g1, g2, g
      real(8)                :: Eo1, Eo2, Ei1, Ei2
      real(8), allocatable   :: e_grid_tmp2(:)

      ! First take the incoming energy grid of inelastic and elastic, and combine
      ! with the energy group structure, energy_bins.
      if (allocated(sab % elastic_e_in)) then
        call merge(sab % inelastic_e_in, sab % elastic_e_in, e_grid_tmp)
        call merge(e_grid_tmp, energy_bins, Ein)
        ! Now I want to cap this at the maximum value
        max_ein = max(sab % inelastic_e_in(sab % n_inelastic_e_in), &
                      sab % elastic_e_in(sab % n_elastic_e_in))
      else
        call merge(sab % inelastic_e_in, energy_bins, Ein)
        ! Now I want to cap this at the maximum value
        max_ein = sab % inelastic_e_in(sab % n_inelastic_e_in)
      end if

      if (allocated(e_grid_tmp)) deallocate(e_grid_tmp)

      ! We want to have an incoming point whenever an inelastic Eout
      ! point crosses a group boundary.
      if (sab % secondary_mode /= SAB_SECONDARY_CONT) then
        do i = 1, sab % n_inelastic_e_in - 1
          Ei1 = sab % inelastic_e_in(i)
          Ei2 = sab % inelastic_e_in(i + 1)
          do j = 1, sab % n_inelastic_e_out
            ! This is slow, but will due until we can sort the resultant array
            allocate(e_grid_tmp(size(energy_bins)))
            num_pts = 0
            Eo1 = sab % inelastic_e_out(j, i)
            ! Ifs to check bounds!
            g1 = binary_search(energy_bins, size(energy_bins), Eo1)
            Eo2 = sab % inelastic_e_out(j, i + 1)
            ! Ifs to check bounds!
            g2 = binary_search(energy_bins, size(energy_bins), Eo2)
            ! Switch data if *2 point is not the higher energy point
            if (Eo2 < Eo1) then
              g = g1
              g2 = g1
              g1 = g
            end if
            ! If g1 /= g2 then we have a critical point in the middle.
            ! Find the Ein which makes this Eout equal to Eg using log interp
            do g = g1 + 1, g2
              num_pts = num_pts + 1
              ! For linear interpolation:
              e_grid_tmp(num_pts) = (energy_bins(g) - Eo1) / (Eo2 - Eo1) * &
                (Ei2 - Ei1) + Ei1
              ! Logarithmic Interpolation:
              !e_grid_tmp(num_pts) = Ei11 * ((energy_bins(g) - Eo1) / (Eo2 - Eo1) * &
              !  log (Ei2 / Ei1))
            end do
            ! Now combine the results on to the Ein grid
            if (num_pts > 0) then
              ! Put relevant pts of e_grid_tmp on to Ein grid, stored in e_grid_tmp2
              call merge(e_grid_tmp(1: num_pts), Ein, e_grid_tmp2)
              deallocate(Ein)
              ! Now put them back in to Ein
              allocate(Ein(size(e_grid_tmp2)))
              Ein = e_grid_tmp2
              deallocate(e_grid_tmp2)
            end if
            deallocate(e_grid_tmp)
          end do
        end do
      end if

      ! Now find where the threshold energy is for s(a,b) on this new grid
      i_max_ein = binary_search(Ein, size(Ein), max_ein)

      if (allocated(e_grid_tmp)) deallocate(e_grid_tmp)
      allocate(e_grid_tmp(size(Ein)))
      e_grid_tmp = Ein
      deallocate(Ein)

      ! Now expand the number of points per currently existing Ein point
      if (SAB_EPTS_PER_BIN == 0) then
        allocate(Ein(i_max_ein))
        Ein = e_grid_tmp(1:i_max_ein)
      else
        num_pts = (i_max_ein - 1) * (EXTEND_PTS) + i_max_ein

        allocate(Ein(num_pts))
        j = 0
        do iE = 1, i_max_ein - 1
          dE = (log(e_grid_tmp(iE + 1) / e_grid_tmp(iE))) / real(EXTEND_PTS + 1, 8)
          j = j + 1
          Ein(j) = e_grid_tmp(iE)
          do i = 1, EXTEND_PTS
            j = j + 1
            Ein(j) = Ein(j - 1) * exp(dE)
          end do
        end do
        Ein(num_pts) = e_grid_tmp(i_max_ein)
      end if

    end subroutine sab_egrid

end module sab
