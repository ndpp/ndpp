module chidata_header

  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use fission,          only: nu_total, nu_delayed
  use global,           only: nuclides
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, header, print_ascii_array
  use search,           only: binary_search
  use string,           only: to_str

#ifdef HDF5
  use hdf5_interface
#endif

  implicit none

!===============================================================================
! CHIDATA Stores the data for each fission energy spectrum of a given nuclide
!===============================================================================

  type :: ChiData
    logical              :: is_init = .false. ! Initialization Status
    integer              :: NE = 0            ! Number of Ein values
    integer              :: NR = 0            ! Number of interpolation regions
    real(8), allocatable :: E_grid(:)         ! Ein values
    real(8), pointer     :: E_bins(:)         ! Energy grp boundaries from input
    integer              :: groups     =  0   ! Number of outgoing energy groups
    type(Nuclide),  pointer   :: nuc   => NULL() ! Working nuclide
    type(Reaction), pointer   :: rxn   => NULL() ! Working reaction
    type(DistEnergy), pointer :: edist => NULL() ! Reaction's energy distro
    real(8), pointer     :: sigma(:) => NULL() ! Reaection cross-section set
    integer              :: law               ! Energy distribution law
    logical              :: is_delayed        ! Prompt or delayed reaction
    integer              :: precursor_grp = 0 ! If delayed, the precursor group #

    ! Type-Bound procedures
    contains
      procedure :: init  => chi_init  ! Sets NE, allocates spaces, etc
      procedure :: clear => chi_clear ! Deallocates this object
      procedure :: beta =>  chi_beta  ! Calculates Beta value
      procedure :: prob =>  chi_prob  ! Calculates probability of rxn value
      procedure :: integrate => chi_integrate ! Performs integration of edist
  end type ChiData

contains

!===============================================================================
! CHI_INIT Initializes the scatt_data object.  This entails finding the energy
! grid, allocating accordingly and storing this information.
!===============================================================================

    subroutine chi_init(self, nuc, rxn, sigma, edist, E_bins, is_delayed, &
                        precursor_grp)
      class(ChiData), intent(inout)         :: self  ! Working object
      type(Nuclide), pointer, intent(in)    :: nuc   ! Nuclide we are working on
      type(Reaction), pointer, intent(in)   :: rxn   ! Reaction of interest
      real(8), pointer, intent(in)          :: sigma(:) ! X/S of reaction of interest
      type(Distenergy), pointer, intent(in) :: edist ! The energy distribution to use
      real(8), target, intent(in)           :: E_bins(:)  ! Energy group bounds
      logical, intent(in)                   :: is_delayed ! Is prompt (F) or delayed (T)
      integer, optional, intent(in)         :: precursor_grp ! Precursor grp

      integer :: lc ! Edist % data location counter

      self % is_init = .false.

      ! Store nuclide
      self % nuc => nuc

      ! Store if we are delayed or not
      self % is_delayed = is_delayed
      ! Store the precursor group
      if (is_delayed) then
        if (present(precursor_grp)) then
          self % precursor_grp = precursor_grp
        else
          call fatal_error("Precursor Group Must Be Provided For Delayed Chi Data!")
        end if
        ! Store the reaction type
        self % sigma => null()
        self % rxn => null()
      else
        self % precursor_grp = 0
        ! Store the reaction type
        self % sigma => sigma
        self % rxn => rxn
      end if

      ! Allocate the group-dependent attributes
      self % E_bins => E_bins
      self % groups = size(E_bins) - 1

      ! Store energy dist information
      self % edist => edist
      self % law = edist % law
      self % NR = int(edist % data(1))
      ! Error checking:
      !if (self % NR /= 0) then
      !  call fatal_error("NR /= 0 for Law " // trim(to_str(self % law)))
      !end if
      self % NE = int(edist % data(2 + 2 * self % NR))
      allocate(self % E_grid(self % NE))
      lc = 2 + 2 * self % NR
      ! Get Ein values from the energy distribution
      self % E_grid = edist % data(lc + 1 : lc + self % NE)

      ! The final initialization
      self % is_init = .true.
    end subroutine chi_init

!===============================================================================
! CHI_CLEAR Clears (deallocates) the Chi object.
!===============================================================================

    subroutine chi_clear(self)
      class(ChiData), intent(inout) :: self  ! Working object

      self % NE = 0
      self % NR = 0
      if (allocated(self % E_grid)) deallocate(self % E_grid)
      nullify(self % E_bins)
      self % groups = 0
      nullify(self % nuc)
      nullify(self % rxn)
      nullify(self % sigma)
      nullify(self % edist)
      self % law = 0
      self % precursor_grp = 0
      self % is_init = .false.
    end subroutine chi_clear

!===============================================================================
! CHI_BETA Calculates Beta (or 1-beta for prompt) at a given incoming energy
!===============================================================================

  function chi_beta(self, Ein) result(beta)
    class(ChiData), intent(in) :: self  ! Working object
    real(8), intent(in)        :: Ein   ! Incoming energy
    real(8)                    :: beta  ! beta value to return

    beta = nu_delayed(self % nuc, Ein) / nu_total(self % nuc, Ein)

  end function chi_beta

!===============================================================================
! CHI_PROB Calculates the probability of the working reaction compared to all
! the other fission reactions of this nuclide at a given incoming energy.
! For delayed neutron channels, this becomes the yield of that precursor group
!===============================================================================

  function chi_prob(self, Ein) result(prob)
    class(ChiData), intent(in) :: self  ! Working object
    real(8), intent(in)        :: Ein   ! Incoming energy
    real(8)                    :: prob  ! probability value to return

    integer :: NR, NE, lc, j ! Various counters
    real(8) :: f             ! Interpolant for prompt

    if (self % is_delayed) then
      ! Get the yield data
      lc = 1
      ! Our lc of interest depends on the data of previous groups, so we have to
      ! 'fast-forward' through the previous groups to get to our data.
      do j = 1, self % nuc % n_precursor
        ! determine number of interpolation regions and energies
        NR = int(self % nuc % nu_d_precursor_data(lc + 1))
        NE = int(self % nuc % nu_d_precursor_data(lc + 2 + 2 * NR))
        if (j == self % precursor_grp) then
          exit
        else
          ! Update the lc pointer
          lc = lc + 2 + 2*NR + 2*NE + 1
        end if
      end do

      ! Determine delayed neutron precursor yield for our precursor group
      prob = &
        interpolate_tab1(self % nuc % nu_d_precursor_data(lc+1:lc+2+2*NR+2*NE),Ein)
    else
      ! Do the prompt data
      ! calculate probabilty of self nested distribution over the others
      if (Ein < self % nuc % energy(1)) then
        j = 1
        f = ZERO
      elseif (Ein >= self % nuc % energy(self % nuc % n_grid)) then
        j = self % nuc % n_grid - 1
        f = ONE
      else
        j = binary_search(self % nuc % energy, self % nuc % n_grid, Ein)
        f = (Ein - self % nuc % energy(j)) / &
          (self % nuc % energy(j + 1) - self % nuc % energy(j))
      end if
      ! check for rare case where two energy points are the same
      if (self % nuc % energy(j) == self % nuc % energy(j + 1)) then
        j = j + 1
      end if

      ! Get the probability of law validity. If block needed because of Pu-240 issue.
      if (j < self % rxn % threshold) then
        prob = ZERO
      else
        ! Find the probability using linear interpolation
        prob = ((ONE - f) * self % sigma(j - self % rxn % threshold + 1) + &
                f * self % sigma(j - self % rxn % threshold + 2)) / &
               ((ONE - f) * self % nuc % fission(j) + f * self % nuc % fission(j + 1))
      end if
      if (associated(self % edist % next) .and. self % edist % p_valid % n_regions > 0) then
        prob = prob * interpolate_tab1(self % edist % p_valid, Ein)
      end if
    end if

  end function chi_prob

!===============================================================================
! CHI_INTEGRATE Integrates the fission energy spectrum over the outgoing groups
!===============================================================================

  function chi_integrate(self, Ein) result(chis)
    class(ChiData), intent(in) :: self    ! Working object
    real(8), intent(in)        :: Ein     ! Incoming energy
    real(8), allocatable       :: chis(:) ! Resultant chi distribution

    type(DistEnergy), pointer :: edist
    integer :: NR, NE, NP, iE, INTTp, INTT, ND, INTT_in
    integer :: lEout_min
    real(8) :: T              ! Fission Spectra Theta Value
    real(8) :: U              ! Restriction Energy
    real(8) :: I              ! Non-dimensional Fiss Spectra param.
    real(8) :: x, x0          ! Spectra constants
    integer :: lc             ! location in the data array
    integer :: g              ! E group indices
    real(8) :: Egp1           ! upper bound of integral
    real(8) :: Eg             ! lower bound of integral
    real(8) :: Watt_a, Watt_b ! Watt spectrum values
    real(8) :: interp         ! interpolation value of the energy point
    real(8) :: runsum         ! Running sum of chi integration for law 4

    allocate(chis(self % groups))
    chis = ZERO
    edist => self % edist

    ! Determine which secondary energy distribution law to use
    select case (self % law)
    case (1)
      ! =======================================================================
      ! TABULAR EQUIPROBABLE ENERGY BINS
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // &
                   " Not Yet Supported.")
    case (3)
      ! =======================================================================
      ! INELASTIC LEVEL SCATTERING
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // &
                   " Not Yet Supported.")
    case (4, 61)
      ! =======================================================================
      ! CONTINUOUS TABULAR DISTRIBUTION AND
      ! CORRELATED ENERGY AND ANGLE DISTRIBUTION

      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      if (NR == 1) then
        call warning("Assuming linear-linear interpolation when sampling &
                     &continuous tabular distribution")
      else if (NR > 1) then
        call fatal_error("Multiple interpolation regions not supported while &
                     &attempting to sample continuous tabular distribution.")
      end if
      ! Set interpolation type (assuming lin-lin)
      INTT_in = LINEAR_LINEAR

      ! find energy bin and calculate interpolation factor -- if the energy is
      ! outside the range of the tabulated energies, choose the first or last
      ! bins
      lc = 2 + 2*NR
      if (Ein < edist % data(lc+1)) then
        iE = 1
        x = ZERO
      elseif (Ein >= edist % data(lc+NE)) then
        iE = NE-1
        x = ONE
      else
        iE = binary_search(edist % data(lc+1:lc+NE), NE, Ein)
        x = (Ein - edist%data(lc+iE)) / (edist%data(lc+iE+1) - edist%data(lc+iE))
      end if

      if (x > 0.5_8) then
        iE = iE + 1
      end if

      ! determine location of outgoing energies, pdf, cdf for E(l)
      lc = int(edist % data(2 + 2*NR + NE + iE))

      ! determine type of interpolation and number of discrete lines
      INTTp = int(edist % data(lc + 1))
      NP    = int(edist % data(lc + 2))
      if (INTTp > 10) then
        INTT = mod(INTTp,10)
        ND = (INTTp - INTT)/10
      else
        INTT = INTTp
        ND = 0
      end if
      if (ND > 0) then
        ! discrete lines present
        call fatal_error("Discrete lines in continuous tabular distributed not &
                         &yet supported")
      end if

      ! Loop through energy groups
      lc = lc + 3
      lEout_min = lc
      runsum = ZERO
      do g = 1, self % groups
        ! Find the location in edist%data corresponding to the upper
        ! bound of the energy group.
        ! The lower bound is already known: lEout_min
        do iE = lEout_min, NP + lc - 2
          if (edist % data(iE + 1) > self % E_bins(g + 1)) exit
        end do
        if (iE == NP + lc - 1) iE = iE - 1
        ! Since we have the CDF, the integral of chis(g) is simply:
        ! cdf_chis(E_high) - cdf_chis(E_low)
        ! Calculate cdf_chis(E_high),set equal to chis(g)
        ! Will do this with LINEAR_LINEAR, regardless of what INTT says
        interp = (self % E_bins(g + 1) - edist % data(iE)) / &
          (edist % data(iE + 1) - edist % data(iE))

        chis(g) = (edist % data(iE + 2 * NP) + interp * &
          (edist % data(iE + 1 + 2 * NP) - edist % data(iE + 2 * NP)))

        ! The above put the upper bound of the integral in to chis(g);
        ! still need to take off the lower bound.
        ! We have been keeping track of this in runsum - the lower bound
        ! is simply the sum of integrals of all lower groups
        chis(g) = chis(g) - runsum
        ! Update runsum
        runsum = runsum + chis(g)
        ! Move locator up in the world
        lEout_min = iE
      end do

    case (5)
      ! =======================================================================
      ! GENERAL EVAPORATION SPECTRUM
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported.")
    case (7)
      ! =======================================================================
      ! MAXWELL FISSION SPECTRUM
      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))

      ! determine nuclear temperature from tabulated function
      T = interpolate_tab1(edist % data, Ein)

      ! determine restriction energy
      lc = 2 + 2*NR + 2*NE
      U = edist % data(lc + 1)
      if (Ein - U <= ZERO) return
      x = (Ein - U) / T
      I = sqrt(T*T*T) * (sqrt(0.25_8*PI) * erf(x) - x * exp(-x))
      do g = 1, self % groups
        ! The upper energy boundary
        Egp1 = self % E_bins(g+1)
        if (Egp1 > Ein - U) Egp1 = U
        chis(g) = 0.5_8 * (sqrt(PI * T) * erf(sqrt(Egp1/T)) * exp(Egp1/T) - &
          TWO * sqrt(Egp1)) * T * exp(-Egp1/T)
        ! The lower energy boundary
        Eg = self % E_bins(g)
        if (Eg > Ein - U) Eg = U
        chis(g) = chis(g) - &
          (0.5_8 * (sqrt(PI * T) * erf(sqrt(Eg/T))*exp(Eg/T) - &
          TWO * sqrt(Eg)) * T * exp(-Eg/T))
        chis(g) = chis(g) / I
      end do

    case (9)
      ! =======================================================================
      ! EVAPORATION SPECTRUM
      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))

      ! determine nuclear temperature from tabulated function
      T = interpolate_tab1(edist % data, Ein)

      ! determine restriction energy
      ! These were derived using a sage worksheet; note that in the equation for I,
      ! the MCNP5 manual has exp(x) instead of exp(-x) which is found on the T2
      ! website (http://t2.lanl.gov/nis/endf/intro25.html); the T2 formulation
      ! makes sense and is correct.
      lc = 2 + 2*NR + 2*NE
      U = edist % data(lc + 1)
      x = (Ein - U) / T
      if (Ein - U <= ZERO) return
      do g = 1, self % groups
        Egp1 = self % E_bins(g+1)
        Eg = self % E_bins(g)
        if (Egp1 > (Ein - U)) Egp1 = Ein - U
        if (Eg > (Ein - U)) Eg = Ein - U
        chis(g) = (Egp1*exp(x) + T*exp(x))*exp(-Egp1 / T)
        chis(g) = chis(g) - (Eg*exp(x) + T*exp(x))*exp(-Eg / T)
        chis(g) = chis(g) / (T*(x - exp(x) + ONE))
      end do
    case (11)
      ! =======================================================================
      ! ENERGY-DEPENDENT WATT SPECTRUM
      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))

      ! Get interpolation type
      if (NR == 0) then
        INTT_in = LINEAR_LINEAR
      else
        write(*,*) 'Error, INTT /= LINEAR_LINEAR', edist % law
      end if

      ! determine Watt parameter 'a' from tabulated function
      Watt_a = interpolate_tab1(edist % data, Ein)
      ! determine Watt parameter 'b' from tabulated function
      lc = 2 + 2*(NR + NE)
      Watt_b = interpolate_tab1(edist % data, Ein, lc + 1)
      ! read number of interpolation regions and incoming energies for
      ! parameter 'a'
      NR = int(edist % data(lc + 1))
      NE = int(edist % data(lc + 2 + 2*NR))
      ! determine restriction energy
      lc = lc + 2 + 2*(NR + NE)
      U = edist % data(lc + 1)
      x = (Ein - U) / Watt_a
      if (Ein - U <= ZERO) return
      x0 = Watt_a * Watt_b * 0.25_8
      I = 0.25_8 * sqrt(PI * Watt_a ** 3 * Watt_b) * exp(x0) * &
        (erf(sqrt(x) - sqrt(x0)) + erf(sqrt(x) + sqrt(x0))) - &
        Watt_a * exp(-x * sinh(Watt_a * Watt_b * x))
      !Reuse Watt_b, and x
      Watt_b = sqrt(Watt_b)
      x = sqrt(PI * Watt_a) * Watt_b * exp(0.25_8 * Watt_a * Watt_b**2)
      do g = 1, self % groups
        Egp1 = self % E_bins(g+1)
        if (Egp1 > U) Egp1 = U
        chis(g) = &
          (-x * erf((Watt_a * Watt_b - TWO * sqrt(Egp1)/(TWO * Watt_a))) + &
          x * erf((Watt_a * Watt_b + TWO * sqrt(Egp1)/(TWO * Watt_a))) - &
          TWO * (exp(TWO * Watt_b * sqrt(Egp1)) * exp(-(Watt_a * Watt_b * sqrt(Egp1))/Watt_a)))
        Eg = self % E_bins(g)
        if (Eg > U) Eg = U
        chis(g) = chis(g) - &
          (-x * erf((Watt_a * Watt_b - TWO * sqrt(Eg)/(TWO * Watt_a))) + &
          x * erf((Watt_a * Watt_b + TWO * sqrt(Eg)/(TWO * Watt_a))) - &
          TWO * (exp(TWO * Watt_b * sqrt(Eg)) * exp(-(Watt_a * Watt_b * sqrt(Eg))/Watt_a)))
        chis(g) = 0.25_8 * Watt_a *chis(g) / I
      end do
    case (12)
      ! =======================================================================
      ! Madland-Nix Fission Spectrum
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported.")
    case (44)
      ! =======================================================================
      ! KALBACH-MANN CORRELATED SCATTERING
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported.")
    case (66)
      ! =======================================================================
      ! N-BODY PHASE SPACE DISTRIBUTION
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported.")
    case (67)
      ! =======================================================================
      ! LABORATORY ENERGY-ANGLE LAW
      call warning("Energy Distribution Type " // trim(to_str(edist % law)) // " Not Yet Supported.")
    end select

    ! Normalize chis(g) in case the interpolation rule causes it to be > 1.0
    I = ZERO
    do g = 1, self % groups
      I = I + chis(g)
    end do
    if (I /= ONE) then
      I = ONE / I
      do g = 1, self % groups
        chis(g) = chis(g) * I
      end do
    end if

  end function chi_integrate

end module chidata_header