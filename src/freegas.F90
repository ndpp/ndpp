module freegas

  use constants
  use error,            only: fatal_error, warning
  use global,           only: nuclides, message
  use legendre
  use search,           only: binary_search
  use string,           only: to_str

  implicit none

  contains

!===============================================================================
! INTEGRATE_FREEGAS_LEG Finds Legendre moments of the energy-angle
! distribution of an elastic collision using the free-gas scattering kernel.
!===============================================================================

  subroutine integrate_freegas_leg(Ein, A, kT, fEmu, mu, E_bins, order, distro)
    real(8), intent(in)  :: Ein         ! Incoming energy of neutron
    real(8), intent(in)  :: A           ! Atomic-weight-ratio of target
    real(8), intent(in)  :: kT          ! Target Temperature (MeV)
    real(8), intent(in)  :: fEmu(:)     ! Energy-angle distro to act on
    real(8), intent(in)  :: mu(:)       ! fEmu angular grid
    real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
    integer, intent(in)  :: order       ! Number of moments to find
    real(8), intent(out) :: distro(:,:) ! Resultant integrated distribution

    integer :: g             ! outgoing energy group index
    real(8) :: p0_1g_norm    ! normalization constant so that P0 = 1.0
    integer :: l             ! Scattering order
    real(8) :: Eout_lo, Eout_hi ! Low and High bounds of Eout integration
    real(8) :: Elo, Ehi      ! Low and High bounds of Eout integration
                             ! for each grp
    real(8) :: alphaEin      ! Eout lower limit for target-at-rest elastic

    ! This routine does the double integration of the free-gas kernel
    ! using an adaptive simpsons integration scheme for both Eout and mu.
    ! mu is the inner integration.

    ! Find alphaEin; we will use this for setting a breakpoint in our adaptive
    ! Eout integration
    alphaEin = (A - ONE) / (A + ONE)
    alphaEin = alphaEin * alphaEin * Ein

    ! Set the normalization constant to zero so we can tally it w/ each group
    p0_1g_norm = ZERO

    ! Calculate the lower and upper bounds of integration
    call calc_FG_Eout_bounds(A, kT, Ein, Eout_lo, Eout_hi)

    do g = 1, size(E_bins) - 1
      if ((E_bins(g) < Eout_hi) .and. (E_bins(g + 1) > Eout_lo)) then
        ! Now lets set the lower and upper bounds of integration
        ! for this group which
        ! will progress through the rest of these steps
        if (Eout_lo > E_bins(g)) then
          Elo = Eout_lo
        else
          Elo = E_bins(g)
        end if
        if (Eout_hi < E_bins(g + 1)) then
          Ehi = Eout_hi
        else
          Ehi = E_bins(g + 1)
        end if

        ! Integrate the tails of the distribution (low grp boundary to Elo,
        ! Ehi to high group boundary)
        ! We do this because it is essentially free anyways (the tails
        ! should be smooth and thus very few points are needed), and has
        ! led to small but sensible improvements in accuracy.
        do l = 1, order
          distro(l, g) = &
            adaptiveSimpsons_Eout(A, kT, Ein, l - 1, fEmu, mu, &
            E_bins(g), Elo) + &
            adaptiveSimpsons_Eout(A, kT, Ein, l - 1, fEmu, mu, &
            Ehi, E_bins(g + 1))
        end do

        ! If alphaEin exists in the group of interest use it as a
        ! adaptive simpsons break point
        if ((Elo < alphaEin) .and. (alphaEin < Ehi)) then
          do l = 1, order
            distro(l, g) = distro(l, g) + &
              adaptiveSimpsons_Eout(A, kT, Ein, l - 1, fEmu, mu, &
              Elo, alphaEin)
          end do
          Elo = alphaEin
        end if

        ! If Ein exists in the group of interest use it as a
        ! adaptive simpsons break point
        if ((Elo < Ein) .and. (Ein < Ehi)) then
          do l = 1, order
            distro(l, g) = distro(l, g) + &
              adaptiveSimpsons_Eout(A, kT, Ein, l - 1, fEmu, mu, &
              Elo, Ein)
          end do
          Elo = Ein
        end if

        ! Do the remainder (Elo will have been set to whatever the start of
        ! the remainder is by the time we get here)
        do l = 1, order
          distro(l, g) = distro(l, g) + &
            adaptiveSimpsons_Eout(A, kT, Ein, l - 1, fEmu, mu, &
            Elo, Ehi)
        end do

      else ! What the heck, do the integral anyways, should be pretty cheap.
        do l = 1, order
          distro(l, g) = &
            adaptiveSimpsons_Eout(A, kT, Ein, l - 1, fEmu, mu, &
            E_bins(g), E_bins(g + 1))
        end do
      end if

      ! Tally the normalization constant.
      p0_1g_norm = p0_1g_norm + distro(1, g)
    end do

    ! And normalize
    distro = distro / p0_1g_norm
  end subroutine integrate_freegas_leg

!===============================================================================
! CALC_FG_EOUT_BOUNDS determines the outgoing energy (Eout) integration bounds
! to use such that the negligible 'tails' of the S(a,b) distribution can
! be ignored, saving computational time.
!===============================================================================

  pure subroutine calc_FG_Eout_bounds(A, kT, Ein, Eout_lo, Eout_hi)
    real(8), intent(in) :: A    ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(out) :: Eout_lo ! Low bound of Eout
    real(8), intent(out) :: Eout_hi ! High bound of Eout

    real(8) :: alpha            ! Multiplicative factor which represents
                                ! maximum energy loss in a pure TAR collision.

    alpha = ((A - ONE) / (A + ONE))**2

    ! alpha*Ein is the minimum Eout of a TAR elastic collision.
    ! Let's take 5% of that energy, to give room for the fact that
    ! the target is not at rest and this isn't a hard boundary anymore.
    ! This... is willy nilly.
    Eout_lo = 0.005_8 * alpha * Ein

    ! The only way for upscatter to occur is with the target nuclide having
    ! a large maxwellian energy. So, add on larger than the mean to Ein.
    ! Again, this is willy nilly.
    if (Ein > 300.0_8 * kT / A) then
      Eout_hi = 12.0_8 * kT * (A + ONE) / A + 1.5_8 * Ein
    else
      Eout_hi = 12.0_8 * kT * (A + ONE) / A + TWO * Ein
    end if

  end subroutine calc_FG_Eout_bounds

!===============================================================================
! CALC_SAB calculates the value of S(a,b) for the free-gas scattering kernel.
! Note that beta is provided, but alpha is calculated within the function.
!===============================================================================

  pure function calc_sab(A, kT, Ein, Eout, beta, mu) result(sab)
    real(8), intent(in) :: A    ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(in) :: Eout ! Outgoing energy of neutron
    real(8), intent(in) :: beta ! Energy Transfer
    real(8), intent(in) :: mu   ! Angle in question

    real(8) :: sab   ! The result of this function
    real(8) :: alpha ! Momentum Transfer
    real(8) :: lterm ! The leading term from the rest of the integral with S(a,b)
                     ! calculated in this routine to save FLOPs elsewhere.

    ! These are limits used by NJOY99 to help with numerical stability.
    real(8), parameter :: alpha_min = 1.0E-6_8
    real(8), parameter :: sab_min   = -225.0_8
    real(8), parameter :: lterm_min = 2.0E-10_8

    ! Find the leading term (the part not inside S(a,b), but still in eqn)
    lterm = sqrt(Eout / Ein) / kT * ((A + ONE) / A) ** 2

    ! Calculate momentum transfer, alpha
    alpha = (Ein + Eout - TWO * mu * sqrt(Ein * Eout)) / (A * kT)
    if (alpha < alpha_min) then
      alpha = alpha_min
    end if

    ! Find the argument to the exponent in S(a,b), for testing against
    ! sab_min
    sab = -(alpha + beta)**2 / (4.0_8 * alpha) ! The sab exp argument

    if (sab < sab_min) then
      sab = ZERO
    else
      ! We have an acceptable value, plug in and move on.
      sab = lterm * exp(sab) / (sqrt(4.0_8 * PI * alpha))
      if (sab < lterm_min) then
        sab = ZERO
      end if
    end if
  end function calc_sab

!===============================================================================
! BRENT_MU uses Brents method of root finding to determine where the S(a,b)
! function crosses the provided threshold as a function of mu.
!===============================================================================

  pure function brent_mu(awr, kT, Ein, Eout, beta, thresh, lo, hi) result(mu_val)
    real(8), intent(in) :: awr    ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(in) :: Eout ! Outgoing energy of neutron
    real(8), intent(in) :: beta ! Energy Transfer
    real(8), intent(in) :: thresh ! SAB Threshold to search for
    real(8), intent(in) :: lo   ! Low mu value boundary
    real(8), intent(in) :: hi   ! High mu value boundary

    real(8) :: mu_val
    real(8) :: a, b, c, d, fa, fb, fc, s, fs
    real(8) :: tmpval
    integer :: i     ! Iteration counter
    logical :: mflag ! Method flag

    a = lo
    b = hi
    c = ZERO
    d = INFINITY

    fa = calc_sab(awr, kT, Ein, Eout, beta, a) - thresh
    fb = calc_sab(awr, kT, Ein, Eout, beta, b) - thresh

    fc = ZERO
    s  = ZERO
    fs = ZERO

    ! Check the bounds, exit if we are not in bounds
    if (fa * fb >= ZERO) then
      if (fa < fb) then
        mu_val = a
      else
        mu_val = b
      end if
      return
    end if

    ! Now, we will switch a and b if abs(fa) < abs(fb)
    if (abs(fa) < abs(fb)) then
      tmpval = a
      a = b
      b = tmpval
      tmpval = fa
      fa = fb
      fb = tmpval
    end if

    ! Set up our initial run through

    c = a
    fc = fa
    mflag = .true.
    i = 0

    do while ((fb /= ZERO) .and. (abs(a - b) > BRENT_MU_THRESH))
      if ((fa /= fc) .and. (fb /= fc)) then
        ! Inverse quadratic interpolation
        s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / &
          (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb)
      else
        ! Secant Rule
        s = b - fb * (b - a) / (fb - fa)
      end if

      tmpval = (3.0_8 * a + b) * 0.25_8
      if ((.not. (((s > tmpval) .and. (s < b)) .or. &
        ((s < tmpval) .and. (s > b)))) .or. &
        (mflag .and. (abs(s - b) >= (0.5_8 * abs(b - c)))) .or. &
        (.not. mflag .and. (abs(s - b) >= (abs(c - d) * 0.5_8)))) then
          s = 0.5_8 * (a + b)
          mflag = .true.
      else
        if ((mflag .and. (abs(b - c) < BRENT_MU_THRESH)) .or. &
          (.not. mflag .and. (abs(c - d) < BRENT_MU_THRESH))) then
          s = (a + b) * 0.5_8
          mflag = .true.
        else
            mflag = .false.
        end if
      end if

      fs = calc_sab(awr, kT, Ein, Eout, beta, s) - thresh
      d = c
      c = b
      fc = fb

      if (fa * fs < ZERO) then
        b = s
        fb = fs
      else
        a = s
        fa = fs
      end if

      ! Now swap a and b if we need to again
      if (abs(fa) < abs(fb)) then
        tmpval = a
        a = b
        b = tmpval
        tmpval = fa
        fa = fb
        fb = tmpval
      end if

      i = i + 1
    end do

    mu_val =  b

  end function brent_mu

!===============================================================================
! FIND_FG_MU calculates the angular boundaries (for a given Ein, Eout pair) to
! use for the integration over the change in angle ($\mu$) variable.  This is
! needed because S(a,b) approaches a delta function as the incoming energy
! increases (to the order of eV - keV).  Delta functions, obviously, are very
! difficult to numerically integrate, so this routine puts the angular points
! where they are needed, saving FLOPs.
!===============================================================================

  pure subroutine find_FG_mu(A, kT, Ein, Eout, mu)
    real(8), intent(in) :: A    ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(in) :: Eout ! Outgoing energy of neutron
    real(8), intent(inout) :: mu(:) ! Angle points to analyze

    real(8) :: mu_max    ! mu which gives the peak of sab(mu)
    real(8) :: beta      ! Energy Transfer
    real(8) :: alpha_max ! Value of alpha corresponding to mu_max
    real(8) :: sab_max   ! Maximum sab value in [-1,1]
    real(8) :: sab_minthresh ! Minimum value of sab to consider
    real(8) :: mu_lo, mu_hi  ! Low, hi mu pts

    ! Calculate beta
    beta = (Eout - Ein) / kT

    ! Calculate the alpha corresponding to the max of s(a,b)
    ! This was derived by ds(a,b)/da = 0
    alpha_max = sqrt(beta * beta + ONE) - ONE
    ! Find the mu values corresponding to alpha_max
    mu_max = (Ein + Eout - alpha_max * A * kT) / (TWO * sqrt(Ein * Eout))

    ! Now that I have mu_max, lets chack the values to the left and right
    ! side to try and find when we can start inspecting sab

    ! First, if mu_max is outside of [-1,1], then this is a relatively
    ! flat profile and we should be integrating all of the range
    if (abs(mu_max) > ONE) then
      mu_lo = -ONE
      mu_hi = ONE
    else
      ! Lets first start checking the low side.
      ! Our threshold is based on the value of sab
      ! at the maximum, so lets find the maximum value
      sab_max = calc_sab(A, kT, Ein, Eout, beta, mu_max)
      sab_minthresh = sab_max * SAB_THRESHOLD
      if (calc_sab(A, kT, Ein, Eout, beta, -ONE) > sab_minthresh) then
        mu_lo = -ONE
      else
        mu_lo = brent_mu(A, kT, Ein, Eout, beta, sab_minthresh, -ONE, mu_max)
      end if
      if (calc_sab(A, kT, Ein, Eout, beta, ONE) > sab_minthresh) then
        mu_hi = ONE
      else
        mu_hi = brent_mu(A, kT, Ein, Eout, beta, sab_minthresh, mu_max, ONE)
      end if
    end if

    ! Now set the return values
    mu(1) = mu_lo
    mu(2) = mu_hi

  end subroutine find_FG_mu

!===============================================================================
! CALC_FGK calculates the value of the free-gas scattering kernel.
!===============================================================================

  function calc_fgk(awr, kT, Ein, Eout, l, mu, fEmu, global_mu) &
    result(fgk)

    real(8), intent(in) :: awr    ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(in) :: Eout ! Outgoing energy of neutron
    integer, intent(in) :: l    ! Legendre order
    real(8), intent(in) :: mu   ! Angle in question
    real(8), intent(in) :: fEmu(:) ! Energy-angle distro to act on
    real(8), intent(in) :: global_mu(:)   ! fEmu angular grid

    real(8) :: fgk   ! The result of this function
    real(8) :: alpha ! Momentum Transfer
    real(8) :: beta  ! Energy Transfer
    real(8) :: lterm ! The leading term from the rest of the integral with S(a,b)
                     ! calculated in this routine to save FLOPs elsewhere.
    integer :: i                    ! Global mu index
    real(8) :: interp, fEmu_val     ! interpolation fraction and interpolated value
                                    ! of the angular distribution

    ! Find fEmu val to use
    if (mu <= global_mu(1)) then
        i = 1
      else if (mu >= global_mu(size(global_mu))) then
        i = size(global_mu) - 1
      else
        i = binary_search(global_mu, size(global_mu), mu)
      end if
    interp = (mu - global_mu(i)) / (global_mu(i + 1) - global_mu(i))
    fEmu_val = (ONE - interp) * fEmu(i) + interp * fEmu(i + 1)

    ! Find the leading term (the part not inside S(a,b), but still in eqn)
    lterm = fEmu_val * sqrt(Eout / Ein) / kT * ((awr + ONE) / awr) ** 2

    ! Calculate momentum transfer, alpha
    alpha = (Ein + Eout - TWO * mu * sqrt(Ein * Eout)) / (awr * kT)

    ! Calculate energy transfer, beta
    beta = (Eout - Ein) / kT

    if (alpha < 1.0E-6_8) alpha = 1.0E-6_8

    ! Find the argument to the exponent in S(a,b), for testing against
    ! sab_min
    fgk = -(alpha + beta)**2 / (4.0_8 * alpha) ! The sab exp argument
    if (fgk <= -2300.0_8) then
      ! This is to avoid a floating point exception b/c exp(fgk) was too small
      ! Really would be better treated by putting a maximum Eout value
      ! on our iterations.
      fgk = ZERO
    else
      fgk = lterm * exp(fgk) / (sqrt(4.0_8 * PI * alpha)) * calc_pn(l, mu)
    end if

  end function calc_fgk

!===============================================================================
! ADAPTIVESIMPSONS_MU and ADAPTIVESIMPSONSAUX_MU in conjunction perform adaptive
! simpsons integration over the angular variable (mu) between bounds a and b.
! The code for these routines was adapted from that found at:
! http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method#C
!===============================================================================

  function adaptiveSimpsons_mu(awr, kT, Ein, Eout, l, fEmu, global_mu, &
    a, b) result(integral)

    real(8), intent(in) :: awr  ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(in) :: Eout ! Outgoing energy of neutron
    integer, intent(in) :: l    ! Legendre order
    real(8), intent(in) :: fEmu(:) ! Energy-angle distro to act on
    real(8), intent(in) :: global_mu(:)   ! fEmu angular grid
    real(8), intent(in) :: a
    real(8), intent(in) :: b

    real(8) :: c, h, fa, fb, fc, S
    real(8) :: integral

    c = (a + b)* 0.5_8
    h = (b - a)

    fa = calc_fgk(awr, kT, Ein, Eout, l, a, fEmu, global_mu)
    fb = calc_fgk(awr, kT, Ein, Eout, l, b, fEmu, global_mu)
    fc = calc_fgk(awr, kT, Ein, Eout, l, c, fEmu, global_mu)

    S = (h / 6.0_8) * (fa + 4.0_8 * fc + fb)
    integral =  adaptiveSimpsonsAux_mu(awr, kT, Ein, Eout, l, fEmu, global_mu, &
      a, b, ADAPTIVE_MU_TOL, S, fa, fb, fc, ADAPTIVE_MU_ITS)

  end function adaptiveSimpsons_mu

  recursive function adaptiveSimpsonsAux_mu(awr, kT, Ein, Eout, l, &
    fEmu, global_mu, a, b, eps, S, fa, fb, fc, bottom) result(val)

    real(8), intent(in) :: awr  ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    real(8), intent(in) :: Eout ! Outgoing energy of neutron
    integer, intent(in) :: l    ! Legendre order
    real(8), intent(in) :: fEmu(:) ! Energy-angle distro to act on
    real(8), intent(in) :: global_mu(:)   ! fEmu angular grid
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: eps
    real(8), intent(in) :: S
    real(8), intent(in) :: fa
    real(8), intent(in) :: fb
    real(8), intent(in) :: fc
    integer, intent(in) :: bottom

    real(8) :: c, d, h, e, fd, fe, Sleft, Sright, S2
    real(8) :: val

    c = 0.5_8 * (a + b)
    h = b - a
    d = 0.5_8 * (a + c)
    e = 0.5_8 * (c + b)
    fd = calc_fgk(awr, kT, Ein, Eout, l, d, fEmu, global_mu)
    fe = calc_fgk(awr, kT, Ein, Eout, l, e, fEmu, global_mu)

    Sleft  = (h / 12.0_8) * (fa + 4.0_8 * fd + fc)
    Sright = (h / 12.0_8) * (fc + 4.0_8 * fe + fb)
    S2 = Sleft + Sright

    if ((bottom <= 0) .or. (abs(S2 - S) <= 15.0_8 * eps)) then
      val = S2 + (S2 - S) / 15.0_8
    else
      val = adaptiveSimpsonsAux_mu(awr, kT, Ein, Eout, l, fEmu, global_mu, &
              a, c, 0.5_8 * eps, Sleft,  fa, fc, fd, bottom - 1) + &
            adaptiveSimpsonsAux_mu(awr, kT, Ein, Eout, l, fEmu, global_mu, &
              c, b, 0.5_8 * eps, Sright, fc, fb, fe, bottom - 1)
    end if

  end function adaptiveSimpsonsAux_mu


!===============================================================================
! ADAPTIVESIMPSONS_EOUT and ADAPTIVESIMPSONSAUX_EOUT in conjunction perform
! adaptive simpsons integration over the outgoing Energy (Eout) between bounds
! a and b.  The code for these routines was adapted from that found at:
! http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method#C
!===============================================================================

  function adaptiveSimpsons_Eout(awr, kT, Ein, l, fEmu, global_mu, &
    a, b) result(integral)

    real(8), intent(in) :: awr  ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    integer, intent(in) :: l    ! Legendre order
    real(8), intent(in) :: fEmu(:) ! Energy-angle distro to act on
    real(8), intent(in) :: global_mu(:)   ! fEmu angular grid
    real(8), intent(in) :: a
    real(8), intent(in) :: b

    real(8) :: c, h, fa, fb, fc, S
    real(8) :: mu_a(2), mu_b(2), mu_c(2)
    real(8) :: integral

    c = 0.5_8 * (a + b)
    h = b - a

    call find_FG_mu(awr, kT, Ein, a, mu_a)
    call find_FG_mu(awr, kT, Ein, b, mu_b)
    call find_FG_mu(awr, kT, Ein, c, mu_c)

    fa = adaptiveSimpsons_mu(awr, kT, Ein, a, l, &
            fEmu, global_mu, mu_a(1), mu_a(2))
    fb = adaptiveSimpsons_mu(awr, kT, Ein, b, l, &
            fEmu, global_mu, mu_b(1), mu_b(2))
    fc = adaptiveSimpsons_mu(awr, kT, Ein, c, l, &
            fEmu, global_mu, mu_c(1), mu_c(2))

    S = (h / 6.0_8) * (fa + 4.0_8 * fc + fb)
    integral =  adaptiveSimpsonsAux_Eout(awr, kT, Ein, l, fEmu, global_mu, &
      a, b, ADAPTIVE_EOUT_TOL, S, fa, fb, fc, ADAPTIVE_EOUT_ITS)
  end function adaptiveSimpsons_Eout

  recursive function adaptiveSimpsonsAux_Eout(awr, kT, Ein, l, &
    fEmu, global_mu, a, b, eps, S, fa, fb, fc, bottom) result(val)

    real(8), intent(in) :: awr  ! Atomic-weight-ratio of target
    real(8), intent(in) :: kT   ! Target Temperature (MeV)
    real(8), intent(in) :: Ein  ! Incoming energy of neutron
    integer, intent(in) :: l    ! Legendre order
    real(8), intent(in) :: fEmu(:) ! Energy-angle distro to act on
    real(8), intent(in) :: global_mu(:)   ! fEmu angular grid
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: eps
    real(8), intent(in) :: S
    real(8), intent(in) :: fa
    real(8), intent(in) :: fb
    real(8), intent(in) :: fc
    integer, intent(in) :: bottom

    real(8) :: c, d, h, e, fd, fe, Sleft, Sright, S2
    real(8) :: mu_d(2), mu_e(2)
    real(8) :: val

    c = 0.5_8 * (a + b)
    d = 0.5_8 * (a + c)
    e = 0.5_8 * (c + b)
    h = b - a

    call find_FG_mu(awr, kT, Ein, d, mu_d)
    call find_FG_mu(awr, kT, Ein, e, mu_e)

    fd = adaptiveSimpsons_mu(awr, kT, Ein, d, l, &
            fEmu, global_mu, mu_d(1), mu_d(2))
    fe = adaptiveSimpsons_mu(awr, kT, Ein, e, l, &
            fEmu, global_mu, mu_e(1), mu_e(2))

    Sleft  = (h / 12.0_8) * (fa + 4.0_8 * fd + fc)
    Sright = (h / 12.0_8) * (fc + 4.0_8 * fe + fb)
    S2 = Sleft + Sright
    if ((bottom <= 0) .or. (abs(S2 - S) <= 15.0_8*eps)) then
      val = S2 + (S2 - S) / 15.0_8
    else
      val = adaptiveSimpsonsAux_Eout(awr, kT, Ein, l, fEmu, global_mu, &
              a, c, 0.5_8 * eps, Sleft, fa, fc, fd, bottom - 1) + &
            adaptiveSimpsonsAux_Eout(awr, kT, Ein, l, fEmu, global_mu, &
              c, b, 0.5_8 * eps, Sright, fc, fb, fe, bottom - 1)
    end if
  end function adaptiveSimpsonsAux_Eout

end module freegas