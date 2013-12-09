module sab

  use ace_header, only: SAlphaBeta
  use constants
  use error,      only: fatal_error, warning
  use global
  use legendre
  use search,     only: binary_search
  use string,     only: to_str

  implicit none

  contains

!===============================================================================
! INTEGRATE_SAB_EL performs integration of the thermal scattering elastic
! collisions and places the result in sab_int
!===============================================================================

    subroutine integrate_sab_el(sab, ein_grid, e_bins, scatt_type, order, &
                                sab_int, sig)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]
      real(8), allocatable, intent(inout)   :: sig(:)         ! Micros. x/s

      integer :: iEin      ! incoming energy counter
      integer :: groups    ! shorthand for number of energy groups
      integer :: g, imu, l ! indices for: group, angle, legendre order
      real(8) :: mu        ! cosine of angle of scatter
      real(8) :: Ein       ! Incoming energy
      integer :: isab      ! index on sab Ein grid
      real(8) :: f         ! fraction for interpolation

      ! Initialize data
      sab_int = ZERO
      groups = size(e_bins) - 1
      allocate(sig(size(ein_grid)))

      ! We have to exit if there are no elastic reactions
      if (sab % threshold_elastic == ZERO) then
        return
      end if
      !!!TD: Parallelize this over Ein

      do iEin = 1, size(ein_grid)
        Ein = ein_grid(iEin)
        ! Find the index and interpolation factor
        if (Ein < sab % elastic_e_in(1)) then
          isab = 1
          f = ZERO
        else if (Ein > sab % threshold_elastic) then
          return
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
          return
        else
          g = binary_search(e_bins, size(e_bins), Ein)
        end if

        ! Get x/s for normalizing
        if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          sig(iEin) = sab % elastic_P(isab) / Ein
        else if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
          sig(iEin) = (ONE - f) * sab % elastic_P(isab) + &
            f * sab % elastic_P(isab + 1)
        end if

        if (sab % n_elastic_mu == 0) then
          ! Perform calculations for coherent elastic (u = 1-E_{bragg}/E)
          mu = ONE - sab % elastic_e_in(isab) / Ein
          do l = 1, order + 1
            sab_int(l, g, iEin) = sab_int(l, g, iEin) + calc_pn(l - 1, mu)
          end do
          sab_int(:, g, iEin) = sig(iEin) * sab_int(:, g, iEin)
        else if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
          ! Loop over outgoing angles and add legendre of each
          do imu = 1, sab % n_elastic_mu
            ! Find our interpolated mu
            mu = (ONE - f) * sab % elastic_mu(imu, isab) + &
              f * sab % elastic_mu(imu, isab + 1)
            ! Probably could write a function to do this using Legendre recursion
            ! avoiding the loop over order. Oh well, this should be pretty cheap
            do l = 1, order + 1
              sab_int(l, g, iEin) = sab_int(l, g, iEin) + calc_pn(l - 1, mu)
            end do
          end do
          sab_int(:, g, iEin) = sig(iEin) * sab_int(:, g, iEin)
        else
          ! pass
        end if
      end do

    end subroutine integrate_sab_el

!===============================================================================
! INTEGRATE_SAB_INEL performs integration of the thermal scattering inelastic
! collisions and places the result in sab_int
!===============================================================================

    subroutine integrate_sab_inel(sab, ein_grid, e_bins, scatt_type, order, &
                                  sab_int, sig)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]
      real(8), allocatable, intent(inout)   :: sig(:)         ! Micros. x/s

      integer :: iEin      ! incoming energy counter
      integer :: iEout     ! outgoing energy counter
      integer :: groups    ! shorthand for number of energy groups
      integer :: g, imu, l ! indices for: group, angle, legendre order
      real(8) :: mu        ! cosine of angle of scatter
      real(8) :: Ein, Eout ! Incoming & outgoing energies
      integer :: isab      ! index on sab Ein grid
      real(8) :: f         ! fraction for interpolation
      real(8), allocatable :: Eout_wgt(:) ! weighting (based on skewed or equal mode)

      ! Initialize data
      sab_int = ZERO
      groups = size(e_bins) - 1
      allocate(sig(size(ein_grid)))

      ! First lets set up our weighting
      allocate(Eout_wgt(sab % n_inelastic_e_out))
      if (sab % secondary_mode == SAB_SECONDARY_EQUAL) then
        Eout_wgt = ONE / real(sab % n_inelastic_e_out, 8)
      else if (sab % secondary_mode == SAB_SECONDARY_SKEWED) then
        if (sab % n_inelastic_e_out > 4) then
          ! 0.1, 0.4, equally-likely, 0.4, 0.1 (all normalized to 1)
          Eout_wgt(1) = 0.1_8 / real(sab % n_inelastic_e_out, 8)
          Eout_wgt(2) = 0.4_8 / real(sab % n_inelastic_e_out, 8)
          Eout_wgt(3: sab % n_inelastic_e_out - 2) = ONE / &
            real(sab % n_inelastic_e_out, 8)
          Eout_wgt(sab % n_inelastic_e_out - 1) = &
            0.4_8 / real(sab % n_inelastic_e_out, 8)
          Eout_wgt(sab % n_inelastic_e_out) = &
            0.1_8 / real(sab % n_inelastic_e_out, 8)
        else
          ! Leave this in until I know what to do (if any data has this)
          message = "Number of Inelastic Outgoing Energies Less Than 4,&
                    & but Skewed Weighting Requested by Data!"
          call fatal_error()
        end if
      end if

      !!!TD: Parallelize this over Ein
      do iEin = 1, size(ein_grid)
        Ein = ein_grid(iEin)
        ! Find the index and interpolation factor
        if (Ein < sab % inelastic_e_in(1)) then
          isab = 1
          f = ZERO
        else if (Ein > sab % threshold_inelastic) then
          return
        else
          isab = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, Ein)
          f = (Ein - sab % inelastic_e_in(isab)) / &
            (sab % inelastic_e_in(isab + 1) - sab % inelastic_e_in(isab))
        end if

        ! Get x/s for normalizing
        sig(iEin) = (ONE - f) * sab % inelastic_sigma(isab) + &
          f * sab % inelastic_sigma(isab + 1)

        ! Integrate over outgoing energy (outer) and outgoing mu (inner)
        ! How to do it depends upon if the second energy mode is
        ! equally-likely bins, or skewed (or in the future, continuous)
        if ((sab % secondary_mode == SAB_SECONDARY_EQUAL) .or. &
          (sab % secondary_mode == SAB_SECONDARY_SKEWED)) then

          ! Loop through each equally likely energy, find its group
          ! and sum legendre moments of discrete mu to it.
          ! Eout and mu need to be interpolated to
          do iEout = 1, sab % n_inelastic_e_out
            ! Interpolate to our Eout
            Eout = (ONE - f) * sab % inelastic_e_out(iEout, isab) + &
              f * sab % inelastic_e_out(iEout, isab + 1)

            ! Get the outgoing group
            if (Eout < e_bins(1)) then
              cycle
            else if (Eout >= e_bins(size(e_bins))) then
              cycle
            else
              g = binary_search(e_bins, size(e_bins), Eout)
            end if

            ! Now integrate all the discrete mu values (equal weight)
            do imu = 1, sab % n_inelastic_mu
              ! Find our interpolated mu
              mu = (ONE - f) * sab % inelastic_mu(imu, iEout, isab) + &
                f * sab % inelastic_mu(imu, iEout, isab + 1)
              ! Probably could write a function to do this using Legendre recursion
              ! avoiding the loop over order. Oh well, this should be pretty cheap
              do l = 1, order + 1
                sab_int(l, g, iEin) = sab_int(l, g, iEin) + &
                  calc_pn(l - 1, mu) * Eout_wgt(iEout)
              end do
            end do
          end do
        else
          ! pass (continuous goes here)
        end if
        sab_int(:, :, iEin) = sig(iEin) * sab_int(:, :, iEin)
      end do

    end subroutine integrate_sab_inel

!===============================================================================
! COMBINE_SAB_GRID takes elastic and inelastic integrals and puts them on the
! same grid by weighting by the x/s
!===============================================================================

  subroutine combine_sab_grid(sab_int_el, sab_int_inel, sig_el, sig_inel,  &
                              scatt_mat)
    real(8), intent(in) :: sab_int_el(:,:,:)   ! Integrated SAB elastic data [L, G, NEin]
    real(8), intent(in) :: sab_int_inel(:,:,:) ! Integrated SAB inelastic data [L, G, NEin]
    real(8), intent(in) :: sig_el(:)           ! Elastic microscopic x/s on E_grid
    real(8), intent(in):: sig_inel(:)          ! Inelastic microscopic x/s on E_grid
    real(8), allocatable, intent(inout) :: scatt_mat(:,:,:) ! Unionized Scattering Matrices

    integer :: iE

    ! set up our scatt_mat space
    allocate(scatt_mat(size(sab_int_el, 1), size(sab_int_el, 2), &
             size(sab_int_el, 3)))

    do iE = 1, size(sab_int_el, 3)
      scatt_mat(:, :, iE) = (sab_int_el(:, :, iE) + sab_int_inel(:, :, iE)) / &
        (sig_el(iE) + sig_inel(iE))
    end do
  end subroutine combine_sab_grid


end module sab
