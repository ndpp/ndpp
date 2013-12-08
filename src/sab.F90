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

    subroutine integrate_sab_el(sab, e_bins, scatt_type, order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      integer :: iEin      ! incoming energy counter
      real(8) :: sigma     ! cross-section
      integer :: groups    ! shorthand for number of energy groups
      integer :: g, imu, l ! indices for: group, angle, legendre order
      real(8) :: mu        ! angle for coherent elastic
      real(8) :: Ein       ! Incoming energy

      ! Initialize data
      SAB_INT = ZERO
      groups = size(e_bins) - 1

      ! We have to exit if there is no elastic reactions
      if (sab % threshold_elastic == ZERO) then
        return
      end if

      do iEin = 1, sab % n_elastic_e_in
        ! For now, Ein = grid point, need to add finer grid spacing though
        Ein = sab % elastic_e_in(iEin)
        ! Get x/s for normalizing
        sigma = sab % elastic_P(iEin)
        if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          sigma = sigma / Ein
        end if

        ! set outgoing group, g.  No outgoing group integration necessary
        ! because Eout does not change for elastic; Ein == Eout
        g = binary_search(e_bins, size(e_bins), Ein)

        if (sab % n_elastic_mu == 0) then
          ! Perform calculations for coherent elastic (u = 1-E_{bragg}/E)
!!!TD this one requires a funny interpolation, actually.
! This makes me want to implement a custom energy grid before we get here
          mu = ONE - sab % elastic_e_in(iEin) / Ein
          do l = 1, order + 1
            sab_int(l, g, iEin) = sab_int(l, g, iEin) + calc_pn(l - 1, mu)
          end do
        else if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
          ! Loop over outgoing angles and add legendre of each
          ! For u/zr, there are 69 Ein, but 0 imu.. basically ITCA does not
          ! exist. Is this a bug in ace.F90, or does it simply mean
          ! we should treat angles as isotropic?
          do imu = 1, sab % n_elastic_mu
            ! Probably could write a function to do this using Legendre recursion
            ! avoiding the loop over order. Oh well, this should be pretty cheap
            do l = 1, order + 1
              sab_int(l, g, iEin) = sab_int(l, g, iEin) + &
                calc_pn(l - 1, sab % elastic_mu(imu, iEin))
            end do
          end do
          sab_int(:, g, iEin) = sigma * sab_int(:, g, iEin)
        else
          ! pass
        end if
write(*,*) Ein, g, sab_int(:, g, iEin)
      end do
write(17, *) " "
write(17, *) sab_int

    end subroutine integrate_sab_el

!===============================================================================
! INTEGRATE_SAB_INEL performs integration of the thermal scattering inelastic
! collisions and places the result in sab_int
!===============================================================================

    subroutine integrate_sab_inel(sab, e_bins, scatt_type, order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      integer :: iEin  ! incoming energy counter

    end subroutine integrate_sab_inel

end module sab
