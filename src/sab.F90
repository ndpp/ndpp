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

    subroutine integrate_sab_el(sab, ein_grid, e_bins, scatt_type, order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
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
      integer :: isab      ! index on sab Ein grid

      ! Initialize data
      sab_int = ZERO
      groups = size(e_bins) - 1

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
        sigma = sab % elastic_P(isab)
        if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          sigma = sigma / Ein
        end if

        if (sab % n_elastic_mu == 0) then
          ! Perform calculations for coherent elastic (u = 1-E_{bragg}/E)
          mu = ONE - sab % elastic_e_in(isab) / Ein
          do l = 1, order + 1
            sab_int(l, g, iEin) = sab_int(l, g, iEin) + calc_pn(l - 1, mu)
          end do
          sab_int(:, g, iEin) = sigma * sab_int(:, g, iEin)
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
write(*,*) Ein, sab_int(:, g, iEin)
      end do

    end subroutine integrate_sab_el

!===============================================================================
! INTEGRATE_SAB_INEL performs integration of the thermal scattering inelastic
! collisions and places the result in sab_int
!===============================================================================

    subroutine integrate_sab_inel(sab, ein_grid, e_bins, scatt_type, order, sab_int)
      type(SAlphaBeta), pointer, intent(in) :: sab            ! Nuclide
      real(8), intent(in)                   :: ein_grid(:)    ! Pre-set incoming E grid
      real(8), intent(in)                   :: e_bins(:)      ! Energy groups
      integer, intent(in)                   :: scatt_type     ! Scattering output type
      integer, intent(in)                   :: order          ! Scattering data order
      real(8), intent(inout)                :: sab_int(:,:,:) ! Integrated SAB data [L, G, NEin]

      integer :: iEin  ! incoming energy counter

    end subroutine integrate_sab_inel

end module sab
