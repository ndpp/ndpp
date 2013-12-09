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

      do iEin = 1, size(ein_grid)
        Ein = ein_grid(iEin)
        if (Ein < sab % elastic_e_in(1)) then
          cycle
        else if (Ein > sab % threshold_elastic) then
          return
        else
          isab = binary_search(sab % elastic_e_in, sab % n_elastic_e_in, Ein)
        end if

        ! set outgoing group, g.  No outgoing group integration necessary
        ! because Eout does not change for elastic; Ein == Eout
        if (Ein <= e_bins(1)) then
          g = 1
        else if (Ein > e_bins(size(e_bins))) then
          return
        else
          g = binary_search(e_bins, size(e_bins), Ein)
        end if

        ! Get x/s for normalizing
        f = (Ein - sab % elastic_e_in(isab)) / &
          (sab % elastic_e_in(isab + 1) - sab % elastic_e_in(isab))
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
              sab_int(l, g, iEin) = sab_int(l, g, iEin) + &
                calc_pn(l - 1, sab % elastic_mu(imu, isab))
            end do
          end do
          sab_int(:, g, iEin) = sig(iEin) * sab_int(:, g, iEin)
        else
          ! pass
        end if
write(*,*) Ein, sab_int(:, g, iEin)
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

      integer :: iEin  ! incoming energy counter

    end subroutine integrate_sab_inel

end module sab
