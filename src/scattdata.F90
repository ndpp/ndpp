module scattdata_class

  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use freegas,          only: integrate_freegas_leg
  use global
  use interpolation,    only: interpolate_tab1
  use legendre
  use output,           only: write_message, header, print_ascii_array
  use search,           only: binary_search
  use string,           only: to_str

  implicit none

!===============================================================================
! JAGGED1D and JAGGED2D is a type which allows for jagged 1-D or 2-D array.
! This is needed for ScattData % distro (2D) and ScattData % Eouts (1D)
!===============================================================================

  type :: jagged2d
    real(8), allocatable :: data(:,:)
  end type jagged2d

  type :: jagged1d
    real(8), allocatable :: data(:)
  end type jagged1d


!===============================================================================
! SCATTDATA Stores the data for each reaction of the nuclide in question.
!===============================================================================

  type :: ScattData
    logical              :: is_init = .false. ! Initialization status
    integer              :: NE = 0            ! Number of Ein values
    real(8), allocatable :: E_grid(:)         ! Ein values

    type(jagged2d), allocatable :: distro(:)  ! Output distribution
                                              ! # distro pts x # E_out x # E_in
    type(jagged1d), allocatable :: Eouts(:)   ! Output Energy values
                                              ! # E_out x # E_in
    type(jagged1d), allocatable :: pdfs(:)    ! Output Energy PDF
                                              ! # E_out x # E_in
    type(jagged1d), allocatable :: cdfs(:)    ! Output Energy CDF
                                              ! # E_out x # E_in
    integer, allocatable :: INTT(:)           ! Interpolation type for each Ein
    real(8), allocatable :: mu(:)             ! mu pts to find f(mu) values at
    real(8), pointer     :: E_bins(:)         ! Energy grp boundaries from input
    integer              :: scatt_type = -1   ! Type of format to store the data
    integer              :: order      =  0   ! Order of the data storage format
    integer              :: groups     =  0   ! Number of outgoing energy groups
    real(8)              :: awr        =  ZERO   ! atomic weight ratio
    type(Reaction),   pointer :: rxn   => NULL() ! My reaction
    type(DistEnergy), pointer :: edist => NULL() ! My reaction's combined dist
    type(DistAngle),  pointer :: adist => NULL() ! My reaction's angle dist
    real(8)              :: freegas_cutoff = ZERO ! Free gas cutoff energy
    real(8)              :: kT         = ZERO ! kT of library
    integer              :: law               ! Scattering law

    ! Type-Bound procedures
    contains
      procedure :: init  => scatt_init  ! Sets NE, allocates spaces, etc
      procedure :: clear => scatt_clear ! Deallocates this object
      procedure :: convert_distro => scatt_convert_distro  ! Converts from ACE to Tabular
      procedure :: interp_distro => scatt_interp_distro ! Interpolate the distro
                                                        ! to the given energy pt
  end type ScattData

  contains

!===============================================================================
! SCATT_INIT Initializes the scatt_data object.  This entails finding the energy
! grid, allocating accordingly and storing this information.
!===============================================================================

    subroutine scatt_init(this, nuc, rxn, edist, E_bins, scatt_type, order, &
      mu_bins)

      class(ScattData), intent(inout)       :: this  ! The object to initialize
      type(Nuclide), intent(in)             :: nuc   ! Nuclide we are working on
      type(Reaction), target, intent(inout) :: rxn   ! The reaction of interest
      type(Distenergy), pointer, intent(in) :: edist ! The energy distribution to use
      ! Edist is intended to specify which of the nested distros we are
      ! actually using. t can be null.
      real(8), target, intent(in)           :: E_bins(:)  ! Energy group bounds
      integer, intent(in) :: scatt_type ! Type of format to store the data
      integer, intent(in) :: order      ! Order of the data storage format
      integer, intent(in) :: mu_bins    ! Number of angular pts in tabular rep.

      integer :: i          ! loop counter
      real(8) :: dmu        ! mu spacing
      integer :: NP, NR     ! Number of outgoing energy values and number of interp regions
      integer :: lc         ! Location of outgoing energy data in edist % data

      ! This test will leave this as uninitialized (and is_init == .false.),
      ! which will be used as a flag when the reactions are combined.
      ! Test reactions to ensure we have a scattering reaction.
      this % is_init = .false.
      if (.not. is_valid_scatter(rxn % MT, rxn % scatter_in_cm)) return

      ! Now, check edist, if passed, and ensure it is of the right law type
      ! before proceeding
      if (associated(edist)) then
        if ((edist % law  /= 3) .and. (edist % law  /= 44) .and. &
          (edist % law  /= 61) .and. (edist % law /= 9)) return
      end if
      ! We survived the above check and thus have a scattering reaction.

      ! Save the order and scattering type
      this % scatt_type = scatt_type
      if (scatt_type == SCATT_TYPE_LEGENDRE) then
        this % order = order + 1
      else
        this % order = order
      end if
      ! Store the reaction type
      this % rxn => rxn
      ! Store the atomic weight ratio
      this % awr = nuc % awr
      ! Store library kT
      this % kT = nuc % kT

      ! Set freegas cutoff;
      ! only do if we are dealing with elastic reactions for now.
      ! Non-elastic reactions will use the default values of zero  .
      if (this % rxn % MT == ELASTIC) then
        ! Store the freegas cutoff energy (MeV)
        this % freegas_cutoff = nuc % freegas_cutoff
      end if

      ! Set distributions
      if (rxn % has_angle_dist) then
        this % adist => rxn % adist
        if (associated(edist)) then ! Like for inelastic level scattering
          this % edist => edist
        else
          this % edist => null()
        end if
      else if (associated(edist)) then
        if ((edist % law == 3) .or. (edist % law == 9)) then
          ! We have either an isotropic inelastic level scatter or an evaporation
          ! spectrum
          this % edist => edist
          ! set up the isotropic adist
          rxn % adist % n_energy = 2
          if (allocated(rxn % adist % energy)) deallocate(rxn % adist % energy)
          allocate(rxn % adist % energy(2))
          ! Set the upper value to that in E_bins(max)
          rxn % adist % energy(2) = E_bins(size(E_bins))
          ! Set the lower value; but, if the threshold of this reaction is above
          ! the lower bound, use that instead
          if (nuc % energy(rxn % threshold) > E_bins(1)) then
            rxn % adist % energy(1) = nuc % energy(rxn % threshold)
          else
            rxn % adist % energy(1) = E_bins(1)
          end if
          ! Set the type to Isotropic
          if (allocated(rxn % adist % type)) deallocate(rxn % adist % type)
          allocate(rxn % adist % type(2))
          rxn % adist % type = ANGLE_ISOTROPIC
          if (allocated(rxn % adist % location)) deallocate(rxn % adist % location)
          allocate(rxn % adist % location(2))
          rxn % adist % location = 0
          if (allocated(rxn % adist % data)) deallocate(rxn % adist % data)
          allocate(rxn % adist % data(2))
          rxn % adist % data = ZERO
          this % adist => rxn % adist
        else
          this % adist => null()
          this % edist => edist ! We dont use rxn % edist because this allows
                                ! ScattData type to be used for nested distros.
        end if
      else
        ! This is an isotropic distribution, to reduce the amount of code later,
        ! lets set up rxn % adist to reflect that and point to it.
        ! This is isotropic. Force rxn % adist to be one that we can read later
        rxn % adist % n_energy = 2
        allocate(rxn % adist % energy(2))
        ! Set the upper value to that in E_bins(max)
        rxn % adist % energy(2) = E_bins(size(E_bins))
        ! Set the lower value; but, if the threshold of this reaction is above
        ! the lower bound, use that instead
        if (nuc % energy(rxn % threshold) > E_bins(1)) then
          rxn % adist % energy(1) = nuc % energy(rxn % threshold)
        else
          rxn % adist % energy(1) = E_bins(1)
        end if
        ! Set the type to Isotropic
        allocate(rxn % adist % type(2))
        rxn % adist % type = ANGLE_ISOTROPIC
        allocate(rxn % adist % location(2))
        rxn % adist % location = 0
        allocate(rxn % adist % data(2))
        rxn % adist % data = ZERO
        this % adist => rxn % adist
        rxn % scatter_in_cm = .true.
        this % edist => null()
      end if

      ! Get the number of energy points; doing so depends on where the info
      ! exists.
      if (associated(this % adist)) then
        this % NE = this % adist % n_energy
        allocate(this % E_grid(this % NE))
        this % E_grid = this % adist % energy
        allocate(this % distro(this % NE))
        NP = 1
        do i = 1, this % NE
          ! There is no Energy-out dependence to the distribution, so we only
          ! need a 1 in the Eout dimension.
          allocate(this % distro(i) % data(mu_bins, NP))
          this % distro(i) % data = ZERO
        end do
      else if ((associated(this % edist)) .and. (this % edist % law /= 3)) then
        NR = int(edist % data(1))
        this % NE = int(edist % data(2 + 2*NR))
        allocate(this % E_grid(this % NE))
        lc = 2 + 2*NR
        this % E_grid(1 : this % NE) = edist % data(lc + 1 : lc + this % NE)
        allocate(this % distro(this % NE))
        do i = 1, this % NE
          ! Find the number of outgoing energy points
          lc = int(edist % data(2+2*NR+this%NE+i))
          NP = int(edist % data(lc + 2))
          allocate(this % distro(i) % data(mu_bins, NP))
          this % distro(i) % data = ZERO
        end do
      end if

      ! Assign the mu points
      allocate(this % mu(mu_bins))
      dmu = TWO / real(mu_bins - 1, 8)
      do i = 1, mu_bins - 1
        this % mu(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      this % mu(size(this % mu)) = ONE

      ! Allocate the Eouts, PDF, CDF jagged container and the interpolation type
      allocate(this % Eouts(this % NE))
      allocate(this % pdfs(this % NE))
      allocate(this % cdfs(this % NE))
      allocate(this % INTT(this % NE))

      ! Allocate the group-dependent attributes
      this % E_bins => E_bins
      this % groups = size(E_bins) - 1

      ! Store the law
      if (associated(this % edist)) then
        this % law = this % edist % law
      else
        this % law = 0
      end if

      ! The final initialization. If we made it here then this bad boy was
      ! successful
      this % is_init = .true.
    end subroutine scatt_init

!===============================================================================
! SCATT_CLEAR Clears (deallocates) the scatt_data object.
!===============================================================================

    subroutine scatt_clear(this)
      class(ScattData), intent(inout) :: this ! The object to clear

      integer :: i ! loop counter

      ! Set attributes to initial values
      this % NE         =  0
      this % scatt_type = -1
      this % order      =  0
      this % groups     =  0
      this % awr        =  ZERO
      this % freegas_cutoff = ZERO
      this % kT         =  ZERO
      this % law        = -1

      ! Reset pointers
      nullify(this % rxn)
      nullify(this % edist)
      nullify(this % adist)
      nullify(this % E_bins)
      ! Deallocate the attribute arrays
      if (allocated(this % E_grid)) then
        deallocate(this % E_grid)
        do i = 1, size(this % distro)
          deallocate(this % distro(i) % data)
          if (allocated(this % Eouts(i) % data)) &
            deallocate(this % Eouts(i) % data)
          if (allocated(this % pdfs(i) % data)) &
            deallocate(this % pdfs(i) % data)
          if (allocated(this % cdfs(i) % data)) &
            deallocate(this % cdfs(i) % data)
        end do
        deallocate(this % distro)
        deallocate(this % Eouts)
        deallocate(this % pdfs)
        deallocate(this % cdfs)
        deallocate(this % INTT)
        deallocate(this % mu)
      end if
      ! Finalize the clear
      this % is_init = .false.
    end subroutine scatt_clear

!===============================================================================
! SCATT_CONVERT_DISTRO Converts the ACE data energy/angle distribution to a
! tabular representation stored in this % distro(:) % data
!===============================================================================

    subroutine scatt_convert_distro(this)
      class(ScattData), intent(inout) :: this ! The object to act on

      integer :: iE ! incoming energy grid index

      ! Check to see if this SD is initialized (if it is not, then it is not
      ! a scattering reaction, or is an invalid law type)
      if (.not. this % is_init) return

      if ((.not. associated(this % edist)) .and. &
        (.not. associated(this % adist))) then
        message = "No distribution associated with this ScattData &
          &object."
        call fatal_error()
      end if

      ! Step through each incoming energy value to do these calculations
      !$omp parallel do schedule(dynamic,50) num_threads(omp_threads) default(shared),private(iE)
      do iE = 1, this % NE
        this % distro(iE) % data = ZERO
        if ((this % law == 0) .or. (this % law == 3)) then
          ! angle distribution only. N_Nx will enter through here.
          call convert_file4(iE, this % mu, this % adist, &
            this % Eouts(iE) % data, this % INTT(iE), &
            this % distro(iE) % data(:, 1))
        else
          ! combined energy/angle distribution.
          call convert_file6(iE, this % mu, this % edist, &
            this % Eouts(iE) % data, this % INTT(iE), this % pdfs(iE) % data, &
            this % cdfs(iE) % data, this % distro(iE) % data)
        end if
      end do
      !$omp end parallel do

    end subroutine scatt_convert_distro

!===============================================================================
! SCATT_INTERP_DISTRO calculates the value of the distribution (times its
! probability) at the optionally given incoming energy.  If no incoming energy
! is given, then no interpolation is performed, and instead we can just use the
! data directly from the given index (iE)
!===============================================================================

    function scatt_interp_distro(this, mu_out, nuc, Ein, norm_tot) &
      result(distro)

      class(ScattData), target, intent(in) :: this ! Working ScattData object
      real(8), intent(in)                  :: mu_out(:) ! The tabular output mu grid
      type(Nuclide), intent(in), pointer   :: nuc  ! Working nuclide
      real(8), intent(inout)               :: norm_tot ! Running total of the micro scattering x/s

      real(8), intent(in)       :: Ein     ! Incoming energy to interpolate on

      type(Reaction), pointer   :: rxn     ! The reaction of interest
      real(8), allocatable :: distro(:,:)  ! the output distribution
      real(8), allocatable :: distro_int(:,:) ! the distribution at Ein before cm2lab
      real(8) :: f, p_valid, sigS          ! interpolation, probability of this
                                           ! (nested) reaction, and \sigma_s(E)
      integer :: iE                        ! incoming energy index (searched)
                                           ! (lower bound if Ein is provided)
      integer :: nuc_iE                    ! Energy index on the nuclide's x/s
      real(8), pointer :: sigS_array(:) => null() ! sigS pointer
      real(8) :: Enorm                    ! Range of Energy space represented by
                                          ! this reaction
      integer :: g                        ! Group index

      ! Set up the results memory
      allocate(distro(this % order, this % groups))
      distro = ZERO
      Enorm  = ZERO

      ! Set rxn, so we can save some characters throughout this function.
      rxn => this % rxn

      ! Point the cross-section pointer to the right place
      if (rxn % MT == ELASTIC) then
        sigS_array => nuc % elastic
      else
        sigS_array => rxn % sigma
      end if

      ! Get sigS and the integrated distro
      if (((Ein <= nuc % energy(rxn % threshold)) .and. &
          (rxn % threshold > 1)) .or. &
        (Ein > this % E_bins(size(this % E_bins)))) then
        ! This is a catch-all, our energy was below the threshold or above the
        ! max group value,
        ! distro should be set to zero and we shall just exit this function
        distro = ZERO
        return
      else if (Ein >= nuc % energy(nuc % n_grid)) then
        ! We are above the global energy grid, so take the highest value
        sigS = sigS_array(size(sigS_array))

        iE = this % NE

        allocate(distro_int(size(this % distro(iE) % data, dim=1), &
          size(this % distro(iE) % data, dim=2)))
        distro_int = this % distro(iE) % data

        distro = integrate_distro(this, Ein, iE, ONE, distro_int, Enorm)
      else
        if (Ein <= nuc % energy(1)) then
          nuc_iE = 1
        else
          nuc_iE = binary_search(nuc % energy, nuc % n_grid, Ein)
        end if
        ! check for rare case where two energy points are the same
        if (nuc % energy(nuc_iE) == nuc % energy(nuc_iE + 1)) &
          nuc_iE = nuc_iE + 1
        ! calculate interpolation factor
        f = (Ein - nuc % energy(nuc_iE)) / &
          (nuc % energy(nuc_iE + 1) - nuc % energy(nuc_iE))
        ! Adjust nuc_iE to point to the sigS_array array
        nuc_iE = nuc_iE - rxn % threshold + 1
        sigS = (ONE - f) * sigS_array(nuc_iE) + f * sigS_array(nuc_iE + 1)

        ! Some reactions have an upper limit on reactions (set by NJOY thresholds)
        ! So, if sigS is zero, (and we already checked for below threshold)
        ! we know that we are above the threshold and can just stop processing
        ! this particular Ein.
        ! One would think that we can stop processing this MT altogether, but
        ! its possible that a x/s was just set to 0 for this particular pt because
        ! its value was below the threshold (perhaps in a resonance dip?)
        if (sigs == ZERO) then
          distro = ZERO
          return
        end if
        ! Search on the angular distribution's energy grid to find what energy
        ! index Ein is at.
        if (Ein < this % E_grid(1)) then
          iE = 1
        else
          iE = binary_search(this % E_grid, this % NE, Ein)
        end if

        ! Interpolate the distribution
        f = (Ein - this % E_grid(iE)) / &
          (this % E_grid(iE + 1) - this % E_grid(iE))
        ! Do on nearest neighbor, or with linear interpolation?
        if (INTERP_NEAREST) then
          if (f >= 0.5_8) iE = iE + 1
          allocate(distro_int(size(this % distro(iE) % data, dim=1), &
            size(this % distro(iE) % data, dim=2)))
          distro_int = this % distro(iE) % data

          distro = integrate_distro(this, Ein, iE, ONE, distro_int, Enorm)

        else if (this % law == 3) then
          allocate(distro_int(size(this % distro(iE) % data, dim=1), &
            size(this % distro(iE) % data, dim=2)))
          distro_int = this % distro(iE) % data

          distro = integrate_distro(this, Ein, iE, ONE, distro_int, Enorm)

        else ! Linear interpolation it is
          ! Do the lower distribution
          allocate(distro_int(size(this % distro(iE) % data, dim=1), &
            size(this % distro(iE) % data, dim=2)))
          distro_int = this % distro(iE) % data
          distro = integrate_distro(this, Ein, iE, (ONE - f), distro_int, Enorm)
          deallocate(distro_int)
          ! Do the upper distribution
          allocate(distro_int(size(this % distro(iE + 1) % data, dim=1), &
            size(this % distro(iE + 1) % data, dim=2)))
          distro_int = this % distro(iE + 1) % data
          distro = distro + integrate_distro(this, Ein, iE + 1, f, distro_int, &
            Enorm)
        end if

      end if

      ! Get the probability value
      if (associated(this % edist)) then
        if (this % edist % p_valid % n_regions > 0) then
          p_valid = interpolate_tab1(this % edist % p_valid, Ein)
        else
          p_valid = ONE
        end if
      else
        p_valid = ONE
      end if

      ! Combine the results, normalizing by the total probability of transfer
      ! from all energies to the energy range represented in the outgoing groups

      ! Combine the terms in to one before multiplying
      sigS = sigS * p_valid * real(rxn % multiplicity, 8)
      do g = 1, this % groups
        distro(:, g) = distro(:, g) * sigS
      end do
      ! Add this contribution to the normalization constant
      ! divide by multiplicity since we dont want to take that in to account
      ! when normalizing
      norm_tot = norm_tot + sigS / real(rxn % multiplicity, 8)

    end function scatt_interp_distro

!===============================================================================
!===============================================================================
! FUNCTIONS TO SUPPORT SCATTDATA
!===============================================================================
!===============================================================================

!===============================================================================
! INTEGRATE_DISTRO finds the energy/angle boundaries of a distribution, then
! integrates over mu/Eout to produce the tabular or legendre distributions
! requested at iE.
!===============================================================================

    function integrate_distro(this, Ein, iE, f, distro_int, Enorm) &
      result(result_distro)

      class(ScattData), target, intent(in) :: this ! Working ScattData object
      real(8), intent(in) :: Ein      ! Incoming energy to interpolate on
      real(8), intent(in) :: f        ! Interpolation factor for this distro_int
      integer, intent(in) :: iE       ! incoming energy index (searched)
      real(8), intent(in) :: distro_int(:,:) ! the distribution at Ein before cm2lab
      real(8), intent(out) :: Enorm   ! Range of Energy space represented by
                                      ! this reaction
      real(8), allocatable :: result_distro(:,:)  ! the output distribution

      ! the distribution in the lab frame
      real(8) :: distro_lab(size(distro_int, dim=1), size(distro_int, dim=2))
      integer :: bins(2, this % groups)   ! Start/end energy/angle bin indices
      real(8) :: vals(2, this % groups)   ! values on start/end energy/angle
      real(8) :: interp(2, this % groups) ! interpolation on start/end energy/angle
      real(8) :: temp_Enorm               ! Storage for Enorm (for summing to Enorm)

      ! Set up the results memory
      allocate(result_distro(this % order, this % groups))
      result_distro = ZERO

      select case (this % scatt_type)
      case (SCATT_TYPE_LEGENDRE)
        if ((associated(this % adist)) .and. &
          (.not. associated(this % edist))) then
          if (Ein < this % freegas_cutoff) then
            distro_lab = distro_int
            call integrate_freegas_leg(Ein, this % awr, this % kT, &
              distro_lab(:, 1), this % mu, this % E_bins, this % order, &
              result_distro)
            temp_Enorm = ONE
          else
            if (this % rxn % scatter_in_cm) then
              call cm2lab(Ein, this % rxn % Q_value, this % awr, this % mu, &
                distro_int, distro_lab)
            else
              distro_lab = distro_int
            end if
            call calc_mu_bounds(this % awr, this % rxn % Q_value, Ein, &
              this % E_bins, this % mu, interp, vals, bins, temp_Enorm)
            ! integrate according to scatt_type
            call integrate_file4_lab_leg(distro_lab(:, 1), this % mu, &
              interp, vals, bins, this % order, result_distro)
          end if
        else if (associated(this % edist)) then
          if (this % rxn % scatter_in_cm) then
            if (this % law == 3) then
                distro_lab = distro_int
                call law3_scatter_cm_leg(distro_lab(:, 1), this % edist, Ein, &
                  this % awr, this % E_bins, this % mu, this % order, &
                  result_distro, temp_Enorm)
            else
              distro_lab = distro_int
              !call integrate_file6_cm_leg(distro_lab, this % mu, Ein, this % awr, &
              !  this % Eouts(iE) % data, this % INTT(iE), this % pdfs(iE) % data, &
              !  this % E_bins, this % order, result_distro, temp_Enorm)
              call integrate_file6_lab_leg(distro_lab, this % mu, &
                this % Eouts(iE) % data, this % INTT(iE), this % pdfs(iE) % data, &
                this % E_bins, this % order, result_distro, temp_Enorm)
            end if
          else
            distro_lab = distro_int
            if (associated(this % adist)) then
              ! Here, we have angular info in adist, but it goes to a single energy
              ! specified in edist. This is for level inelastic scattering and some
              ! (n,xn) reactions.
              if (this % edist % law == 9) then
                ! For these cases, we need to find out which outgoing groups
                ! get the isotropic angular distribution data, which will be
                ! the same for each group.
                ! The group transfer should also be normalized by the prob. of
                ! transfer to that group.
                call law9_scatter_lab_leg(distro_lab(:, 1), this % edist, Ein, &
                  this % E_bins, this % mu, this % order, result_distro, temp_Enorm)
              end if
            else
              call integrate_file6_lab_leg(distro_lab, this % mu, &
                this % Eouts(iE) % data, this % INTT(iE), this % pdfs(iE) % data, &
                this % E_bins, this % order, result_distro, temp_Enorm)
            end if
          end if
        end if
      case (SCATT_TYPE_TABULAR)

      end select

      ! Multiply by f for interpolation.
      result_distro = result_distro * f
      Enorm = Enorm + temp_Enorm * f

    end function integrate_distro

!===============================================================================
! CONVERT_FILE4 performs the overall control for converting file 4 ACE data to
! the required tabular format.
!===============================================================================

    subroutine convert_file4(iE, mu, adist, Eouts, INTT, distro)
      integer, intent(in)  :: iE            ! Energy index to act on
      real(8), intent(in)  :: mu(:)         ! tabular mu points
      type(DistAngle), pointer, intent(in) :: adist    ! My angle dist
      real(8), allocatable, intent(inout)  :: Eouts(:) ! Energy out grid @ Ein
      integer,  intent(inout)              :: INTT     ! Energy out INTT grid
      real(8), intent(inout) :: distro(:) ! resultant distro (# pts)

      real(8), pointer :: data(:) => null() ! Shorthand for adist % data
      integer :: lc           ! Location inside adist % data
      integer :: idata, idata_prev ! Loop counters
      integer :: imu          ! mu loop indices
      integer :: interp, NP   ! Tabular format data (interp type and # pts)
      real(8) :: r            ! Interpolation parameter

      data => adist % data

      lc   = adist % location(iE)

      ! Check what type of distribution we have

      select case(adist % type(iE))
        case (ANGLE_ISOTROPIC)
          distro = 0.5_8
        case (ANGLE_32_EQUI)
          idata_prev = lc + 1
          do imu = 1, size(mu)
            do idata = idata_prev, lc + 1 + NUM_EP
              if (data(idata) >= mu(imu)) then
                ! Find the area of the bin, and use that.
                if (imu == 1) then
                  ! If we are at the start, then use the next data point.
                  distro(imu) = R_NUM_EP /  (data(idata + 1) - data(idata))
                else
                  ! Look backwards to figure it out
                  distro(imu) = R_NUM_EP / (data(idata) - data(idata - 1))
                end if
                idata_prev = idata
                exit
              end if
            end do
          end do
        case (ANGLE_TABULAR)
          interp = int(data(lc + 1))
          NP = int(data(lc + 2))
          lc = lc + 3
          idata_prev = lc
          if (interp == HISTOGRAM) then
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NP - 1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match - set to previous value of PDF (since hist.)
                  distro(imu) = data(idata - 1 + NP)
                  idata_prev = idata
                  exit
                else if (abs(data(idata) - mu(imu)) <= FP_PRECISION) then
                  ! Found a match - set to current value of PDF (since hist.)
                  distro(imu) = data(idata + NP)
                  idata_prev = idata
                  exit
                end if
              end do
            end do
          else if (interp == LINEAR_LINEAR) then
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NP - 1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match - interpolate
                  r = (mu(imu) - data(idata - 1)) / &
                    (data(idata) - data(idata - 1))
                  distro(imu) = data(idata + NP -1) + r * &
                    (data(idata + NP) - data(idata + NP - 1))
                  idata_prev = idata
                  exit
                else if (abs(data(idata) - mu(imu)) <= FP_PRECISION) then
                  ! Found a match - just set to current value of PDF
                  distro(imu) = data(idata + NP)
                  idata_prev = idata
                  exit
                end if
              end do
            end do
          end if
      end select

      ! Finally, set Eouts and INTT
      allocate(Eouts(2))
      Eouts(1) = ZERO
      Eouts(2) = INFINITY
      INTT = HISTOGRAM

    end subroutine convert_file4

!===============================================================================
! CONVERT_FILE6 performs the overall control for converting file 6 ACE data to
! the required tabular format.
!===============================================================================

    subroutine convert_file6(iE, mu, edist, Eouts, INTT, pdf, cdf, distro)
      integer, intent(in)  :: iE            ! Energy index to act on
      real(8), intent(in)  :: mu(:)         ! tabular mu points
      type(DistEnergy), pointer, intent(in) :: edist    ! My energy dist
      real(8), allocatable, intent(inout)   :: Eouts(:) ! Energy out grid @ Ein
      integer, intent(inout)                :: INTT     ! Energy out INTT grid
      real(8), allocatable, intent(inout)   :: pdf(:)   ! Eout PDF
      real(8), allocatable, intent(inout)   :: cdf(:)   ! Eout CDF
      real(8), intent(inout) :: distro(:,:) ! resultant distro (pts x NEout)

      real(8), pointer :: data(:) => null() ! Shorthand for adist % data
      integer :: lcin, lc     ! Locations inside edist % data
      integer :: idata, idata_prev ! Loop counters
      integer :: imu          ! mu loop indices
      integer :: iEout        ! outgoing loop indices and number of points
      integer :: interp, NP   ! Tabular format data (interp type and # pts)
      integer :: NPang        ! Number of angular points
      real(8) :: r            ! Interpolation parameter
      integer :: NR, NE       ! edist % data navigation information
      real(8) :: KMR, KMA, KMconst ! Various data for Kalbach-Mann

      ! Exit if using an unsupported law
      if ((edist % law /= 44) .and. (edist % law /= 61)) return

      data => edist % data

      NR = int(data(1))
      if (NR > 0) then
        message = "Multiple interpolation regions not supported while &
             &attempting to sample Kalbach-Mann distribution."
        call fatal_error()
      end if
      NE = int(data(2 + 2*NR))

      lc = int(data(2 + 2*NR + NE + iE)) ! start of LDAT for iE

      ! determine type of interpolation
      INTT = int(edist % data(lc + 1))
      if (INTT > 10) INTT = mod(INTT,10)

      ! Get number of outgoing energy points at Ein
      NP = int(data(lc + 2))
      allocate(Eouts(NP))
      Eouts = data(lc + 2 + 1 : lc + 2 + NP)
      allocate(pdf(NP))
      allocate(cdf(NP))
      ! Get the Eout PDF and CDF
      pdf = data(lc + 2 + NP + 1 : lc + 2 + 2 * NP)
      cdf = data(lc + 2 + 2 * NP + 1 : lc + 2 + 3 * NP)

      if (edist % law  == 44) then
        lc = lc + 2
        do iEout = 1, NP
          ! Get the KM parameters
          KMR = data(lc + 3 * NP + iEout)
          KMA = data(lc + 4 * NP + iEout)
          ! Calculate the leading term
          KMconst = 0.5_8 * KMA / sinh(KMA)
          distro(:, iEout) = KMconst * (cosh(KMA * mu(:)) + KMR * sinh(KMA * mu(:)))
        end do
      else if (edist % law  == 61) then
        lcin = lc + 2
        do iEout = 1, NP
          lc = int(data(lcin + 3*NP + iEout))

          ! Check if isotropic
          if (lc == 0) then
            distro(:, iEout) = 0.5_8
            cycle
          end if

          interp = int(data(lc + 1))
          NPang = int(data(lc + 2))
          lc = lc + 3

          ! Calculate a PDF value for each mu pt.
          if (interp == HISTOGRAM) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev , lc + NPang -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, take value at mu(imu - 1)
                  distro(imu, iEout) = data(idata + NPang - 1)
                  idata_prev = idata
                  exit
                else if (abs(data(idata) - mu(imu)) <= FP_PRECISION) then
                  ! Found a match, take value at mu(imu)
                  distro(imu, iEout) = data(idata + NPang)
                  idata_prev = idata
                  exit
                end if
              end do
            end do
          else if (interp == LINEAR_LINEAR) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NPang -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, interpolate value
                  r = (mu(imu) - data(idata -1)) / (data(idata) - data(idata-1))
                  distro(imu, iEout) = data(idata + NPang - 1) + r * &
                    (data(idata + NPang) - data(idata - 1 + NPang))
                  idata_prev = idata
                  exit
                else if (abs(data(idata) - mu(imu)) <= FP_PRECISION) then
                  ! Found a match, take value at mu(imu)
                  distro(imu, iEout) = data(idata + NPang)
                  idata_prev = idata
                  exit
                end if
              end do
            end do
          else
            message = "Unknown interpolation type: " // trim(to_str(interp))
            call fatal_error()
          end if

          ! Multiply by the PDF from ENDF
          distro(:, iEout) = distro(:, iEout)

        end do
      end if

    end subroutine convert_file6

!===============================================================================
! CM2LAB Converts a tabular center-of-mass distribution to a laboratory frame
! of reference tabular distribution.
!===============================================================================

    subroutine cm2lab(Ein, Q, awr, mu, data, distro_out)
      real(8), intent(in) :: Ein     ! Incoming energy
      real(8), intent(in) :: Q       ! Reaction Q-value
      real(8), intent(in) :: awr     ! Atomic weight ratio
      real(8), intent(in) :: mu(:)   ! Angular grid
      real(8), intent(in) :: data(:,:) ! The distribution to convert
      real(8), intent(out) :: distro_out(:,:) ! The distribution to convert

      real(8) :: R, Rinv, R2       ! Reduced Effective Mass, 1/R, and R^2
      integer :: imu, iEout        ! mu indices, outgoing energy index
      real(8) :: interp            ! Interpolation parameter b.t wi pts
      integer :: mu_bins           ! Number of mu bins in dim 1 of data
      real(8) :: wi                ! CM cosine corresponding to mu(imu)
      real(8) :: fwi               ! CM distro value at wi
      real(8) :: dwdui             ! CM to Lab Jacobian at wi
      real(8) :: S                 ! Sqrt(mu^2+R^2-1) temp variable
      real(8) :: mu2               ! mu(i)^2 temp variable
      integer :: iwi               ! location in mu(:) of wi
      real(8) :: dmu               ! mu(:) spacing
      real(8) :: Rinv2

      mu_bins = size(data, dim = 1)
      dmu = mu(2) - mu(1)

      ! First check to make sure this is energetically possible
      R2 = -Q * (awr + ONE) / awr
      if (Ein <= R2) then
        distro_out(:,:) = ZERO ! To make sure its energetically impossible
        return
      end if

      ! Now carry on
      ! From equation 234 in Methods for Processing ENDF/B-VII (pg 2798)
      R = awr * sqrt((ONE - R2 / Ein))
      R2 = R * R
      Rinv = ONE / R
      Rinv2 = Rinv * Rinv

      if ((awr > ONE) .and. (R < ONE)) then
        distro_out(:,:) = ZERO ! To make sure its energetically impossible
        return
      end if

      if (R2 < ONE) then
        do imu = 1, mu_bins
          ! Calculat the w to mu conversion terms
          mu2   = mu(imu) * mu(imu)
          if (mu2 < 0.01_8) then
            ! Switch R and Rinv (same for R2 and Rinv2)
            S     = sqrt(mu2 + Rinv2 - ONE)
            wi    = R * (mu2 - ONE + mu(imu) * S)
            dwdui = R * (TWO * mu(imu) + S + mu2 / S)
          else
            S     = sqrt(mu2 + R2 - ONE)
            wi    = Rinv * (mu2 - ONE + mu(imu) * S)
            dwdui = Rinv * (TWO * mu(imu) + S + mu2 / S)
          end if
          ! Now find fwi at wi (linear interpolation)
          iwi    = int((wi + ONE) / dmu) + 1
          if (imu < mu_bins) then
            interp = (wi - mu(iwi)) / (mu(iwi + 1) - mu(iwi))
            ! Now we can get fwi easily, lets do Eout loop now to not have
            ! to recalculate all the above paramters at every imu
            do iEout = 1, size(data, dim = 2)
              fwi = (ONE - interp) * data(iwi, iEout) + &
                    interp * data(iwi + 1, iEout)
              distro_out(imu, iEout) = fwi * dwdui
            end do
          else
            do iEout = 1, size(data, dim = 2)
              fwi = data(iwi, iEout)
              distro_out(imu, iEout) = fwi * dwdui
            end do
          end if
        end do
      else
        do imu = 1, mu_bins
          ! Calculat the w to mu conversion terms
          mu2   = mu(imu) * mu(imu)
          S     = sqrt(mu2 + R2 - ONE)
          wi    = Rinv * (mu2 - ONE + mu(imu) * S)
          dwdui = Rinv * (TWO * mu(imu) + S + mu2 / S)
          ! Now find fwi at wi (linear interpolation)
          iwi    = int((wi + ONE) / dmu) + 1
          if (imu < mu_bins) then
            interp = (wi - mu(iwi)) / (mu(iwi + 1) - mu(iwi))
            ! Now we can get fwi easily, lets do Eout loop now to not have
            ! to recalculate all the above paramters at every imu
            do iEout = 1, size(data, dim = 2)
              fwi = (ONE - interp) * data(iwi, iEout) + &
                    interp * data(iwi + 1, iEout)
              distro_out(imu, iEout) = fwi * dwdui
            end do
          else
            do iEout = 1, size(data, dim = 2)
              fwi = data(iwi, iEout)
              distro_out(imu, iEout) = fwi * dwdui
            end do
          end if
        end do
      end if

    end subroutine cm2lab

!===============================================================================
! CALC_MU_BOUNDS Calculates the mu-values corresponding to the energy-out
! bins for a given input energy.
!===============================================================================

    subroutine calc_mu_bounds(awr, Q, Ein, E_bins, mu, interp, vals, bins, Enorm)

      real(8), intent(in)  :: awr         ! Atomic-weight ratio
      real(8), intent(in)  :: Q           ! Q-Value of this reaction
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: mu(:)       ! tabular mu values
      real(8), intent(out) :: interp(:,:) ! Outgoing mu interpolants
                                          ! corresponding to the energy groups
      real(8), intent(out) :: vals(:,:)   ! Outgoing mu values corresponding
                                          ! to the energy groups
      integer, intent(out) :: bins(:,:)   ! Outgoing mu indices corresponding
                                          ! to the energy groups
      real(8), intent(out) :: Enorm       ! Fraction of possible energy space
                                          ! of this Ein reaction represented by
                                          ! the energy group structure of the
                                          ! problem

      real(8) :: Emin, Emax       ! Max/Min E transfer of this reaction
      real(8) :: R                ! The Reduced Mass (takes in to account Qval)
      real(8) :: mu_low, mu_high  ! Low and high angular points
      integer :: g                ! Group index variable
      integer :: imu              ! angle bin counter

      R = awr * sqrt((ONE + Q * (awr + ONE) / (awr * Ein)))

      do g = 1, size(E_bins) - 1
        ! Calculate the values of mu corresponding to this energy group
        ! These come from eqs. 232-233 in Methods for Processing ENDF/B-VII,
        ! also on pg 2798

        if (R < ONE) then
          ! Set R = 1/R, to conserve the min being alpha*Ein, but
          ! avoid the sqrt(-1) issues.
          if (E_bins(g) == ZERO) then
            mu_low = -ONE
          else
            mu_low  = 0.5_8 * ((ONE + ONE / awr) * sqrt(E_bins(g) / Ein) + &
              (ONE - ONE / (R * R)) / (ONE + ONE / awr) * sqrt(Ein / E_bins(g)))
          end if
          mu_high = 0.5_8 * ((ONE + ONE / awr) * sqrt(E_bins(g + 1) / Ein) + &
            (ONE - ONE / (R * R)) / (ONE + ONE / awr) * sqrt(Ein / E_bins(g + 1)))
        else if (E_bins(g) > ZERO) then
          mu_low  = 0.5_8 * ((ONE + awr) * sqrt(E_bins(g) / Ein) + &
            (ONE - R * R) / (ONE + awr) * sqrt(Ein / E_bins(g)))
          mu_high = 0.5_8 * ((ONE + awr) * sqrt(E_bins(g + 1) / Ein) + &
            (ONE - R * R) / (ONE + awr) * sqrt(Ein / E_bins(g + 1)))
        else
          mu_low = -ONE
          mu_high = 0.5_8 * ((ONE + awr) * sqrt(E_bins(g + 1) / Ein) + &
            (ONE - R * R) / (ONE + awr) * sqrt(Ein / E_bins(g + 1)))
        end if

        ! Find the index in mu corresponding to mu_low
        if (mu_low <= -ONE) then
          bins(MU_LO, g) = 1
          vals(MU_LO, g) = -ONE
          interp(MU_LO, g) = ZERO
        else if (mu_low >= ONE) then
          bins(MU_LO, g) = size(mu) - 1
          vals(MU_LO, g) = ONE
          interp(MU_LO, g) = ONE
        else
          imu = binary_search(mu, size(mu), mu_low)
          bins(MU_LO, g) = imu
          vals(MU_LO, g) = mu_low
          interp(MU_LO, g) = (mu_low - mu(imu)) / (mu(imu + 1) - mu(imu))
        end if

        ! Find the index in mu corresponding to mu_high
        if (mu_high <= -ONE) then
          bins(MU_HI, g) = 1
          vals(MU_HI, g) = -ONE
          interp(MU_HI, g) = ZERO
        else if (mu_high >= ONE) then
          bins(MU_HI, g) = size(mu) - 1
          vals(MU_HI, g) = ONE
          interp(MU_HI, g) = ONE
        else
          imu = binary_search(mu, size(mu), mu_high)
          bins(MU_HI, g) = imu
          vals(MU_HI, g) = mu_high
          interp(MU_HI, g) = (mu_high - mu(imu)) / (mu(imu + 1) - mu(imu))
        end if

      end do

      ! Now set the normalization term (this is the fraction of all energy
      ! space that this Ein -> g transfer represents out of all g.
      Emin = Ein * (ONE + R * R - TWO * R) / ((ONE + awr) * (ONE + awr))
      Emax = Ein * (ONE + R * R + TWO * R) / ((ONE + awr) * (ONE + awr))

      ! Three cases to consider for assigning the normalization:
      ! 1) Emin and Emax within this group (transfer probability = 1)
      ! 2) Emin and Emax outside this group (transfer probability = 0)
      ! 3) Emin, Emax span multiple groups
      ! These correspond to the if branches below.

      if (((Emin >= E_bins(1)) .and. (Emin <= E_bins(size(E_bins)))) .and. &
        ((Emax >= E_bins(1)) .and. (Emax <= E_bins(size(E_bins))))) then

        Enorm = ONE
      else if (((Emin < E_bins(1)) .or. (Emin > E_bins(size(E_bins)))) .and. &
        ((Emax < E_bins(1)) .or. (Emax > E_bins(size(E_bins))))) then

        Enorm = ZERO
      else
        Enorm = (min(Emax, E_bins(size(E_bins))) - max(Emin, E_bins(1))) / &
          (Emax - Emin)
      end if

    end subroutine calc_mu_bounds

!===============================================================================
! LAW3_SCATTER_CM_LEG Finds which group contains the inelastic level outgoing energy
! and sets interp, vals, bins, and Enorm accordingly for later integration.
!===============================================================================

    subroutine law3_scatter_cm_leg(fmu, edist, Ein, awr, E_bins, mu, order, &
                                   distro, Enorm)
      real(8), intent(in)  :: fmu(:)      ! Angle distro to act on
      type(DistEnergy), pointer, intent(in) :: edist ! My energy dist
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: awr         ! atomic weight ratio
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: mu(:)       ! fmu angular grid
      integer, intent(in)  :: order       ! Number of moments to find
      real(8), intent(out) :: distro(:,:) ! Resultant integrated distribution
      real(8), intent(out) :: Enorm       ! Fraction of possible energy space
                                          ! of this Ein reaction represented by
                                          ! the energy group structure of the
                                          ! problem

      integer :: g, imu_c   ! Group index variable, mu point index
      real(8) :: Eo_cm, Eo  ! CM Outgoing Energy, Lab outgoing energy
      real(8) :: dEo        ! change in outgoing energy points
      integer :: iE         ! Outgoing energy counter
      real(8) :: mu_l, mu_c ! lab angle, cm angle
      real(8) :: c          ! c from eq 359
      integer :: NUM_E = 10 ! Number of Eout points per group
      integer :: l          ! Scattering order counter
      real(8) :: proby, f   ! fmu at mu_c, interpolant between fmu points
      real(8) :: integ      ! the integrand
      real(8) :: J          ! the jacobian
      real(8) :: Eo_lo, Eo_hi ! Laboratory Eout bounds
      integer :: g_lo, g_hi ! Lower and upper groups
      real(8), allocatable :: E_bnds(:)
      real(8) :: ap1inv     ! 1/(AWR+1)

      ap1inv = ONE / (awr + ONE)

      Eo_cm = edist % data(2) * (Ein - edist % data(1))
      Eo_lo = Ein * (sqrt(Eo_cm / Ein) - ap1inv)**2
      Eo_hi = Ein * (sqrt(Eo_cm / Ein) + ap1inv)**2

      ! To find Enorm (the total prob'y encapsulated by lower and upper group bnds)
      ! we will solve for mu_l at E_bins(1) and E_bins(size(E_bins)),
      ! and find the fraction of mu_l=[-1,1] that is spanned by this domain.
      ! If this is larger than one, then the energy groups more than cover it.
      c = ap1inv * sqrt(Ein / E_bins(size(E_bins)))
      Enorm = 0.5_8 * (ONE + c * c - Eo_cm / E_bins(size(E_bins))) / c
      if (E_bins(1) == ZERO) then
        Enorm = Enorm + ONE
      else
        c = ap1inv * sqrt(Ein / E_bins(1))
        Enorm = Enorm - 0.5_8 * (ONE + c * c - Eo_cm / E_bins(1)) / c
      end if
      if (abs(Enorm) > TWO) then
        Enorm = ONE
      else
        Enorm = 0.5_8 * Enorm ! Divide by total span of mu_l (1-(-1))
      end if

      g_lo = binary_search(E_bins, size(E_bins), Eo_lo)
      g_hi = binary_search(E_bins, size(E_bins), Eo_hi)
      allocate(E_bnds(g_lo : g_hi + 1))
      E_bnds(g_lo) = Eo_lo
      E_bnds(g_lo + 1: g_hi) = E_bins(g_lo + 1: g_hi)
      E_bnds(g_hi + 1) = Eo_hi

      do g = g_lo, g_hi
        Eo = E_bnds(g)
        dEo = (E_bnds(g + 1) - E_bnds(g)) / real(NUM_E - 1, 8)
        do iE = 1, NUM_E
          ! Adjust Eo to be the top value so that the if Eo==Eo_hi check
          ! will work despite FP precision issues
          if (iE == NUM_E) Eo = E_bnds(g + 1)
          ! Now that we have Ein, Eo, and Eo_cm, find c, mu_l, mu_c, J
          c = ap1inv * sqrt(Ein / Eo)
          J = sqrt(Eo / Eo_cm)
          if (Eo == Eo_lo) then
            mu_l = 0.5_8 * (ONE + c * c - Eo_cm / Eo) / c
            mu_c = -ONE
          else if (Eo == Eo_hi) then
            mu_l = 0.5_8 * (ONE + c * c - Eo_cm / Eo) / c
            mu_c = ONE
          else
            mu_l = 0.5_8 * (ONE + c * c - Eo_cm / Eo) / c
            mu_c = J * (mu_l - c)
          end if
          ! Now we can find fmu(mu=mu_c), using lin-lin interpolation
          if (mu_c == -ONE) then
            proby = fmu(1)
          else if (mu_c == ONE) then
            proby = fmu(size(mu))
          else
            imu_c = binary_search(mu, size(mu), mu_c)
            f = (mu_c - mu(imu_c)) / (mu(imu_c + 1) - mu(imu_c))
            proby = (ONE - f) * fmu(imu_c) + f * fmu(imu_c + 1)
          end if
          integ = proby * J
          if ((iE > 1) .and. (iE < NUM_E)) then
            integ = TWO * integ
          end if
          do l = 1, order
            distro(:, g) = distro(:, g) + integ * calc_pn(l - 1, mu_l)
          end do
          Eo = Eo + dEo
        end do
        distro(:, g) = distro(:, g) * dEo * 0.5_8
      end do

    end subroutine law3_scatter_cm_leg

!===============================================================================
! INTEGRATE_FILE6_CM_LEG Finds which group contains the inelastic level outgoing energy
! and sets interp, vals, bins, and Enorm accordingly for later integration.
!===============================================================================

    subroutine integrate_file6_cm_leg(fEmu, mu, Ein, awr, Eout, INTT, thispdf, &
                                      E_bins, order, distro, Enorm)
      real(8), intent(in)  :: fEmu(:,:)   ! Energy-angle distro to act on
      real(8), intent(in)  :: mu(:)       ! fEmu angular grid
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: awr         ! atomic weight ratio
      real(8), intent(in)  :: Eout(:)     ! Outgoing energies
      integer, intent(in)  :: INTT        ! Interpolation type (Hist || Lin-Lin)
      real(8), intent(in)  :: thispdf(:)  ! Outgoing E-dist PDF from mySD
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      integer, intent(in)  :: order       ! Number of moments to find
      real(8), intent(out) :: distro(:,:) ! Resultant integrated distribution
      real(8), intent(out) :: Enorm       ! Fraction of possible energy space
                                          ! of this Ein reaction represented by
                                          ! the energy group structure of the
                                          ! problem

      integer :: g, imu_c   ! Group index variable, mu point index
      integer :: imu        ! Lab mu index
      real(8) :: Eo_cm, Eo  ! CM Outgoing Energy, Lab outgoing energy
      real(8) :: dEo        ! change in outgoing energy points
      integer :: iE         ! Outgoing energy counter
      real(8) :: mu_c ! cm angle
      real(8) :: mu_l_min   ! Minimum lab angle
      real(8) :: dmu        ! change in lab angle
      real(8) :: c          ! c from eq 359
      integer :: NUM_E = 10 ! Number of Eout points per group
      integer :: l          ! Scattering order counter
      real(8) :: proby, f   ! fmu at mu_c, interpolant between fmu points
      real(8) :: integ      ! the integrand
      real(8) :: J          ! the jacobian
      real(8) :: Eo_lo, Eo_hi ! Laboratory Eout bounds
      real(8), allocatable :: fEl(:) ! Integrated (over mu) fEmu
      real(8), allocatable :: pdf(:) ! Eout probability distrib.
      real(8)              :: pEo    ! probability of this Eo value
      integer              :: iEo    ! Index of Eout array of current Eo
      real(8)              :: fEo    ! interpolant for Eo on Eout grid
      integer :: g_lo, g_hi ! Lower and upper groups
      real(8), allocatable :: E_bnds(:)
      real(8), allocatable :: fmu(:), mu_l(:)

      Enorm = ONE
      allocate(fEl(order))
      allocate(fmu(size(mu)))
      allocate(mu_l(size(mu)))

      ! First lets normalize the PDF
      allocate(pdf(size(Eout)))
      pdf = thispdf
      do iE = 1, size(Eout) - 1
        pdf(iE) = pdf(iE) * (Eout(iE + 1) - Eout(iE))
      end do

      ! Set up lower and upper lab energy boundaries
      Eo_lo = Eout(1) + (Ein - TWO * (awr + ONE) * sqrt(Ein * Eout(1))) &
        / ((awr + ONE)*(awr + ONE))
      Eo_hi = Eout(size(Eout)) + (Ein + TWO * (awr + ONE) * &
        sqrt(Ein * Eout(size(Eout)))) / ((awr + ONE)*(awr + ONE))

      g_lo = binary_search(E_bins, size(E_bins), Eo_lo)
      g_hi = binary_search(E_bins, size(E_bins), Eo_hi)
      allocate(E_bnds(g_lo : g_hi + 1))
      E_bnds(g_lo) = Eo_lo
      E_bnds(g_lo + 1: g_hi) = E_bins(g_lo + 1: g_hi)
      E_bnds(g_hi + 1) = Eo_hi

      do g = g_lo, g_hi
        Eo = E_bnds(g)
        dEo = (E_bnds(g + 1) - E_bnds(g)) / real(NUM_E - 1, 8)
        Eo = Eo - dEo
        do iE = 1, NUM_E
          Eo = Eo + dEo
          fEl = ZERO
          fmu = ZERO
          c = sqrt(Ein / Eo) / (awr + ONE)
          mu_l_min = (ONE + c * c - Eout(size(Eout)) / Eo) / (TWO * c)
          if (mu_l_min < -ONE) then
            mu_l_min = -ONE
          else if (mu_l_min >= ONE) then
            cycle
          end if
          dmu = (ONE - mu_l_min) / real(size(mu) - 1, 8)
          do imu = 1, size(mu)
            mu_l(imu) = mu_l_min + dmu * real(imu - 1, 8)
            Eo_cm = Eo * (ONE + c * c - TWO * c * mu_l(imu))
            if (Eo_cm <= ZERO) cycle
            if (Eo_cm <= Eout(1)) then
              iEo = 1
            else if (Eo_cm >= Eout(size(Eout))) then
              iEo = size(Eout) - 1
            else
              iEo = binary_search(Eout, size(Eout), Eo_cm)
            end if
            if (INTT == HISTOGRAM) then
              fEo = ZERO
              pEo = pdf(iEo)
            else
              fEo = (Eo_cm - Eout(iEo)) / (Eout(iEo + 1) - Eout(iEo))
              pEo = (ONE - fEo) * pdf(iEo) + fEo * pdf(iEo + 1)
            end if
            J = sqrt(Eo / Eo_cm)
            if (mu_l(imu) == -ONE) then
              mu_c = -ONE
            else if (mu_l(imu) == ONE) then
              mu_c = ONE
            else
              mu_c = (mu_l(imu) - c) * J
              if (abs(mu_c) > ONE) cycle
            end if
            if ((mu_c < -ONE) .or. (mu_c > ONE)) cycle
            imu_c = binary_search(mu, size(mu), mu_c)
            f = (mu_c - mu(imu_c)) / (mu(imu_c + 1) - mu(imu_c))
            ! Do lower Eout point
            proby = (ONE - fEo) * &
              ((ONE - f) * fEmu(imu_c, iEo) + f * fEmu(imu_c + 1, iEo))
            ! And now the upper
            proby = proby + fEo * &
              ((ONE - f) * fEmu(imu_c, iEo + 1) + f * fEmu(imu_c + 1, iEo + 1))
            integ = proby * J
            fmu(imu) = integ
          end do
          do imu = 1, size(mu) - 1
            fEl = fEl + &
              calc_int_pn_tablelin(order, mu_l(imu), mu_l(imu + 1), &
              fmu(imu), fmu(imu + 1))
          end do
          if ((iE /= 1) .and. (iE /= NUM_E)) then
            distro(:, g) = distro(:, g) + TWO * pEo * fEl(:)
          else
            distro(:, g) = distro(:, g) + pEo * fEl(:)
          end if
        end do
        distro(:, g) = distro(:, g) * dEo *0.5_8
      end do

    end subroutine integrate_file6_cm_leg

!===============================================================================
! LAW9_SCATTER_LAB_LEG Finds the legendre moments of the incoming adistro and
! assigns the results to each of the outgoing energy groups (according to a
! law 9 evaporation spectrum) according to the probabilities.
!===============================================================================

    subroutine law9_scatter_lab_leg(fmu, edist, Ein, E_bins, mu, order, distro, &
      Enorm)

      real(8), intent(in)  :: fmu(:)      ! Angle distro to act on
      type(DistEnergy), pointer, intent(in) :: edist    ! My energy dist
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: mu(:)       ! fEmu angular grid
      integer, intent(in)  :: order       ! Number of moments to find
      real(8), intent(out) :: distro(:,:) ! Resultant integrated distribution
      real(8), intent(out) :: Enorm       ! Fraction of possible energy space
                                          ! of this Ein reaction represented by
                                          ! the energy group structure of the
                                          ! problem

      integer :: g, NR, NE, lc, imu
      real(8) :: T, U, x, I, Egp1, Eg, pE_xfer

      pE_xfer = ZERO
      Enorm = ZERO

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
      I = T * T * (ONE - exp(-x) * (ONE + x))
      if (U == Ein) return
      if (U < ZERO) U = Ein - U
      do g = 1, size(E_bins) - 1
        Egp1 = E_bins(g + 1)
        if (Egp1 > (Ein - U)) Egp1 = Ein - U
        Eg = E_bins(g)
        if (Eg > (Ein - U)) Eg = Ein - U
        pE_xfer = (T * exp(Egp1 / T) - Egp1 - T) * exp(-Egp1 / T)
        pE_xfer = pE_xfer - (T * exp(Eg / T) - Eg - T) * exp(-Eg / T)
        pE_xfer = T * pE_xfer / I

        ! Now set the angular distributions according to pE_xfer
        do imu = 1, size(mu) - 1
          distro(:, g) = distro(:, g) + &
            calc_int_pn_tablelin(order, mu(imu), mu(imu + 1), &
              fmu(imu), fmu(imu + 1)) * pE_xfer
        end do
        ! Add pE_xfer to running tally in Enorm
        Enorm = Enorm + pE_xfer
      end do

    end subroutine law9_scatter_lab_leg

!===============================================================================
! INTEGRATE_FILE*_LAB_LEG Finds Legendre moments of the energy-angle
! distribution in fEmu over each of the outgoing energy groups and stores the
! result in distro. The FILE4 version does this for only an angular distribution,
! while the FILE6 version does the same for a combined energy-angle distribution
!===============================================================================

    subroutine integrate_file4_lab_leg(fEmu, mu, interp, vals, &
      bins, order, distro)

      real(8), intent(in)  :: fEmu(:)        ! Energy-angle distro to act on
      real(8), intent(in)  :: mu(:)          ! fEmu angular grid
      real(8), intent(in)  :: interp(:,:)    ! interpolants of mu values
      real(8), intent(in)  :: vals(:,:)      ! mu values
      integer, intent(in)  :: bins(:,:)      ! indices of fEmu corresponding to
                                             ! the group boundaries

      integer, intent(in)  :: order          ! Number of moments to find
      real(8), intent(out) :: distro(:,:)    ! Resultant integrated distribution

      integer :: g                           ! outgoing energy group index
      integer :: imu                         ! angle bin index
      real(8) :: flo, fhi                    ! pdf low and high values

      ! Integrating over each of the bins.
      do g = 1, size(distro, dim = 2)
        if (bins(MU_LO, g) /= bins(MU_HI, g)) then
          ! Integrate part of the moment from the low point to the index above
          ! the low point
          flo = fEmu(bins(MU_LO, g)) + interp(MU_LO, g) * &
            (fEmu(bins(MU_LO, g) + 1) - fEmu(bins(MU_LO, g)))
          distro(:, g) = &
            calc_int_pn_tablelin(order, vals(MU_LO, g), &
            mu(bins(MU_LO, g) + 1), flo, fEmu(bins(MU_LO, g) + 1))
          ! Integrate the inbetween pts
          do imu = bins(MU_LO, g) + 1, bins(MU_HI, g) -1
            distro(:, g) = distro(:, g) + &
              calc_int_pn_tablelin(order, mu(imu), &
              mu(imu + 1), fEmu(imu), fEmu(imu + 1))
          end do
          ! Integrate part of the moment from the index below the high point to
          ! the high point
          fhi = fEmu(bins(MU_HI, g)) + interp(MU_HI, g) * &
            (fEmu(bins(MU_HI, g) + 1) - fEmu(bins(MU_HI, g)))
          distro(:, g) = distro(:, g) + &
            calc_int_pn_tablelin(order, mu(bins(MU_HI, g)), vals(MU_HI, g), &
            fEmu(bins(MU_HI, g)), fhi)

        else
          ! The points are all within the same bin, can get flo and fhi directly
          flo = fEmu(bins(MU_LO, g)) + interp(MU_LO, g) * &
            (fEmu(bins(MU_LO, g) + 1) - fEmu(bins(MU_LO, g)))
          fhi = fEmu(bins(MU_HI, g)) + interp(MU_HI, g) * &
            (fEmu(bins(MU_HI, g) + 1) - fEmu(bins(MU_HI, g)))
          distro(:, g) = distro(:, g) + &
            calc_int_pn_tablelin(order, vals(MU_LO, g), vals(MU_HI, g), &
            flo, fhi)

        end if
      end do

    end subroutine integrate_file4_lab_leg

    subroutine integrate_file6_lab_leg(fEmu, mu, Eout, INTT, thispdf, &
      E_bins, order, distro, Enorm)

      real(8), intent(in)    :: fEmu(:,:)     ! Energy-angle distro to act on
      real(8), intent(in)    :: mu(:)         ! fEmu angular grid
      real(8), intent(in)    :: Eout(:)       ! Outgoing energies
      integer, intent(in)    :: INTT          ! Interpolation type (Hist || Lin-Lin)
      real(8), intent(in)    :: thispdf(:)    ! Outgoing E-dist PDF from mySD
      real(8), intent(in)    :: E_bins(:)     ! Energy group boundaries
      integer, intent(in)    :: order         ! Number of moments to find
      real(8), intent(out)   :: distro(:,:)   ! Resultant integrated distro
      real(8), intent(inout) :: Enorm         ! Energy normalization, will be one.

      real(8), allocatable   :: fEmu_int(:,:) ! Integrated (over E) fEmu
      real(8), allocatable   :: pdf(:)        ! local version of pdf to mess with
      integer :: g            ! Energy group index
      integer :: imu          ! angular grid index
      integer :: iE           ! outgoing energy grid index
      integer :: NEout        ! Number of outgoing energies
      integer :: iE_lo, iE_hi ! tmp val of bins(MU_LO, g)
      real(8) :: f_lo, f_hi   ! Interpolation for lo and hi

      allocate(fEmu_int(size(mu), size(E_bins) - 1))
      fEmu_int = ZERO
      NEout = size(Eout)

      Enorm = ONE

      ! First lets normalize the PDF
      allocate(pdf(NEout))
      do iE = 1, NEout - 1
        pdf(iE) = thispdf(iE) * (Eout(iE + 1) - Eout(iE))
      end do

      if (NEout > 1) then
        ! This branch will perform integration of the outgoing energy of fEmu
        ! over each energy group in E_bins.
        ! Since sum(p(E)*deltaE) = 1.0, we know how to do our integration.
        ! sum(f(E)*p(E)*deltaE)
        do g = 1, size(E_bins) - 1
          ! Find iE_lo, and add the first term to fEmu_int
          if (E_bins(g) < Eout(1)) then
            ! We need to skip this lower bound
            iE_lo = 1

          else if (E_bins(g) >= Eout(NEout)) then
            ! In this case, the lower group boundary is above all energies;
            ! this means the group is outside the Eout range, and thus this group
            ! has a zero distribution
            distro(:, g) = ZERO
            cycle

          else
            ! The group boundary is inbetween 2 Eout pts. Interpolate.
            iE_lo = binary_search(Eout, NEout, E_bins(g))
            f_lo = (E_bins(g) - Eout(iE_lo)) / (Eout(iE_lo + 1) - Eout(iE_lo))

            fEmu_int(:, g) = fEmu_int(:, g) + f_lo * pdf(iE_lo) * fEmu(:, iE_lo)
            iE_lo = iE_lo + 1

          end if
          ! Find iE_hi and add the last term to fEmu_int
          if (E_bins(g + 1) < Eout(1)) then
            ! We can skip this group completely then
            distro(:, g) = ZERO
            cycle
          else if (E_bins(g + 1) >= Eout(NEout)) then
            ! The upper grp boundary is above Eout, and so its value is zero
            ! We therefore dont need to add its contribution to the integral
            iE_hi = NEout - 1
          else
            ! The group boundary is inbetween 2 Eout pts. Interpolate.
            iE_hi = binary_search(Eout, NEout, E_bins(g + 1))
            f_hi = (E_bins(g + 1) - Eout(iE_hi)) / (Eout(iE_hi + 1) - Eout(iE_hi))

            fEmu_int(:, g) = fEmu_int(:, g) + f_hi * pdf(iE_hi) * &
              fEmu(:, iE_hi)
            iE_hi = iE_hi - 1
          end if

          ! Now we can do the intermediate points
          do iE = iE_lo, iE_hi
            fEmu_int(:, g) = fEmu_int(:, g) + pdf(iE) * fEmu(:, iE)
          end do

          ! Perform the legendre expansion, and divide by two from trapezoid int
          do imu = 1, size(mu) - 1
            distro(:, g) = distro(:, g) + &
              calc_int_pn_tablelin(order, mu(imu), mu(imu + 1), &
                fEmu_int(imu, g), fEmu_int(imu + 1, g))
          end do

        end do

      else
        ! We do not need to integrate at all, just set distro = fEmu if its'
        ! Eout is in that group (where each group is (Emin, Emax])
        do g = 1, size(E_bins) - 1
          if ((Eout(1) > E_bins(g)) .and. (Eout(1) <= E_bins(g + 1))) then
            ! Perform the legendre expansion
            do imu = 1, size(mu) - 1
              distro(:, g) = distro(:, g) + &
                calc_int_pn_tablelin(order, mu(imu), mu(imu + 1), &
                  fEmu(imu, 1), fEmu(imu + 1, 1))
            end do
          else
            distro(:, g) = ZERO
          end if
        end do
      end if

    end subroutine integrate_file6_lab_leg

!===============================================================================
! INTEGRATE_FILE*_LAB_TAB Finds the tabular distribution representation of
! the energy-angle distrib in fEmu over each of the outgoing energy groups and
! stores the result in distro. The FILE4 version does this for only an angular
! distribution while the FILE6 version does the same for a combined energy-angle
! distribution
!===============================================================================

!!! NOT YET IMPLEMENTED

!===============================================================================
! IS_VALID_SCATTER determines if a given MT number is that of a scattering event
!===============================================================================

  function is_valid_scatter(MT, cm) result(scatter_event)

    integer, intent(in) :: MT ! Reaction channel
    logical, intent(in) :: cm ! Is scatter in CM?
    logical             :: scatter_event

    if ((MT == ELASTIC) .or. &
        ((MT == N_2N) .or. (MT == N_3N) .or. (MT == N_4N)) .or. &
        ((MT >= N_N1) .and. (MT <= N_NC))) then
      scatter_event = .true.
    else
      scatter_event = .false.
    end if

if ((MT == N_2N) .or. (MT == N_3N) .or. (MT == N_4N)) scatter_event = .false.
if (MT == N_NC) scatter_event = .false.
if (MT == ELASTIC) scatter_event = .false.

  end function is_valid_scatter

end module scattdata_class