module scattdata_header

  use ace_header
  use array_merge
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
      ! actually using. It can be null.
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
      if (.not. is_valid_scatter(rxn % MT)) return

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
          if (edist % law == 3) then
            ! Set the edist to null (rxn already has Q_value so we dont need
            ! any data from edist)
            this % edist => null()
          else
            this % edist => edist
          end if
          this % law = edist % law
        else
          this % edist => null()
          this % law = 0
        end if
      else if (associated(edist)) then
        if ((edist % law == 3) .or. (edist % law == 9)) then
          ! We have either an isotropic inelastic level scatter or an evaporation
          ! spectrum
          ! If the reaction is isotropic inelastic, then we should set up an
          ! isotropic reaction, but not point this % edist to edist.
          ! Doin so allows integrate_distro to correctly enter the
          ! integrate_file4_cm_* path.
          if (edist % law == 9) then
            this % edist => edist
          end if
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
        this % law = edist % law
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
        this % law = 0
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

      ! The final initialization
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
        if ((this % law == 0) .or. (this % law == 3) .or. (this % law == 9)) then
          ! angle distribution only. elastic, inelastic level, &
          ! N_Nx will enter through here.
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

    function scatt_interp_distro(this, mu_out, nuc, Ein) result(distro)
      class(ScattData), target, intent(in) :: this ! Working ScattData object
      real(8), intent(in)                  :: mu_out(:) ! The tabular output mu grid
      type(Nuclide), intent(in), pointer   :: nuc  ! Working nuclide

      real(8), intent(in)       :: Ein     ! Incoming energy to interpolate on

      type(Reaction), pointer   :: rxn     ! The reaction of interest
      real(8), allocatable :: distro(:,:)  ! the output distribution
      real(8), allocatable :: distro_int(:,:) ! the distribution at Ein before cm2lab
      real(8) :: f, p_valid, sigS          ! interpolation, probability of this
                                           ! (nested) reaction, and \sigma_s(E)
      integer :: iE                        ! incoming energy index (searched)
                                           ! (lower bound if Ein is provided)
      integer :: nuc_iE                    ! Energy index on the nuclide's x/s
      real(8), pointer :: sigS_array(:)    ! sigS pointer

      ! Set up the results memory
      allocate(distro(this % order, this % groups))
      distro = ZERO

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
        ! distro should be left as zero
        ! and we shall just exit this function
        return
      else if (Ein >= nuc % energy(nuc % n_grid)) then
        ! We are above the global energy grid, so take the highest value
        sigS = sigS_array(size(sigS_array))

        iE = this % NE

        allocate(distro_int(size(this % distro(iE) % data, dim=1), &
          size(this % distro(iE) % data, dim=2)))
        distro_int = this % distro(iE) % data

        distro = integrate_distro(this, Ein, iE - 1)
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
        if (sigS == ZERO) then
          return
        end if
        ! Search on the angular distribution's energy grid to find what energy
        ! index Ein is at.
        if (Ein < this % E_grid(1)) then
          iE = 1
        else
          iE = binary_search(this % E_grid, this % NE, Ein)
        end if

        ! Some reactions, ENDF-VII.0's O-16, and Ca-40 e.g., have two angdist Ein
        ! points in a row being the same value. Check for this and just skip
        ! the first point
        if (this % E_grid(iE) >= this % E_grid(iE + 1)) then
          return
        end if

        distro = integrate_distro(this, Ein, iE)
      end if

      ! Get the probability value
      if (associated(this % edist)) then
        p_valid = interpolate_tab1(this % edist % p_valid, Ein)
      else
        p_valid = ONE
      end if

      if (rxn % MT /= ELASTIC) then
        ! Scale the results
        distro = distro * sigS * p_valid
      end if

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

    function integrate_distro(this, Ein, iE) result(result_distro)

      class(ScattData), target, intent(in) :: this ! Working ScattData object
      real(8), intent(in) :: Ein      ! Incoming energy to interpolate on
      integer, intent(in) :: iE       ! incoming energy index (searched)

      real(8), allocatable :: Eout(:)
      real(8), allocatable :: pdf(:)
      integer              :: INTT
      real(8), allocatable :: fEmu(:,:)
      real(8), allocatable :: result_distro(:,:)  ! the output distribution

      ! the distribution in the lab frame
      real(8), allocatable :: distro(:,:)

      ! Set up the results memory
      allocate(result_distro(this % order, this % groups))
      result_distro = ZERO

      select case (this % scatt_type)
      case (SCATT_TYPE_LEGENDRE)
        if ((associated(this % adist)) .and. &
          (.not. associated(this % edist))) then
          if ((Ein < this % freegas_cutoff) .and. &
              (this % rxn % MT == ELASTIC)) then
            call integrate_freegas_leg(Ein, this % awr, this % kT, &
              this % distro(iE) % data(:, 1), this % mu, this % E_bins, &
              this % order, result_distro)
          else
            if (this % rxn % scatter_in_cm) then
              call integrate_file4_cm_leg(this % distro(iE) % data(:,1), Ein, &
                                          this % awr, this % rxn % Q_value, &
                                          this % E_bins, this % mu, this % order, &
                                          result_distro)
            else
              ! As a notification of future issues:
              message = "File 4 Reaction Found With Lab Angle Distribution and &
                        &No Energy Distribution!"
              call fatal_error()
            end if
          end if

        else if (associated(this % edist)) then
          if (this % rxn % scatter_in_cm) then
            ! Need to come up with the input distro via unit-base interpolation
            call unitbase(this, Ein, iE, Eout, pdf, INTT, fEmu)
            call integrate_file6_cm_leg(fEmu, this % mu, Ein, this % awr, &
              Eout, INTT, pdf, this % E_bins, this % order, result_distro)
          else
            if (associated(this % adist)) then
              ! Here, we have angular info in adist, but it goes to a single energy
              ! specified in edist. This is for level inelastic scattering and some
              ! (n,xn) reactions.
              if (this % law == 9) then
                ! For these cases, we need to find out which outgoing groups
                ! get the isotropic angular distribution data, which will be
                ! the same for each group.
                ! The group transfer should also be normalized by the prob. of
                ! transfer to that group.
                call law9_scatter_lab_leg(this % distro(iE) % data(:, 1), &
                                          this % edist, Ein, this % E_bins, &
                                          this % mu, this % order, result_distro)

              else
                ! As a notification of future issues:
                message = " Associated Edist and Adist, but not law 9: " // &
                  to_str(this % edist % law) //", "//to_str(this % rxn % MT)
              end if
            else
              ! Need to come up with the input distro via unit-base interpolation
              call unitbase(this, Ein, iE, Eout, pdf, INTT, fEmu)
              call integrate_file6_lab_leg(fEmu, this % mu, Eout, INTT, pdf, &
                this % E_bins, this % order, result_distro)
            end if
          end if
        end if
      case (SCATT_TYPE_TABULAR)

      end select

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
          else if (interp == LINEAR_LOG) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NPang -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, interpolate value
                  r = (log(mu(imu)) - log(data(idata - 1)))/ &
                      (log(data(idata)) - log(data(idata - 1)))
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
          else if (interp == LOG_LINEAR) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NPang -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, interpolate value
                  r = (mu(imu) - data(idata -1)) / (data(idata) - data(idata-1))
                  distro(imu,iEout) = exp((ONE-r) * log(data(idata + NPang)) + &
                                      r * log(data(idata + NPang - 1)))
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
          else if (interp == LOG_LOG) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NPang -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, interpolate value
                  r = (log(mu(imu)) - log(data(idata - 1)))/ &
                      (log(data(idata)) - log(data(idata - 1)))
                  distro(imu,iEout) = exp((ONE-r) * log(data(idata + NPang)) + &
                                      r * log(data(idata + NPang - 1)))
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

        end do
      end if

    end subroutine convert_file6

!===============================================================================
! INTEGRATE_FILE4_CM_LEG Calculates the legendre Moments of a file 4 distrib.
!===============================================================================

    subroutine integrate_file4_cm_leg(fw, Ein, awr, Q, E_bins, w, order, distro)
      real(8), intent(in)  :: fw(:)       ! CM Angle distro to act on
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: awr         ! atomic weight ratio
      real(8), intent(in)  :: Q           ! Energy level of excited nucleus
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: w(:)        ! fw angular grid
      integer, intent(in)  :: order       ! Number of moments to find
      real(8), intent(out) :: distro(:,:) ! Resultant integrated distribution

      integer :: g          ! Group index variable, mu point index
      real(8) :: R          ! Reduced Mass
      real(8) :: wlo, whi   ! CM mu-bounds of integration
      integer :: ilo, ihi   ! Grid locations of lo and hi pts
      integer :: iw         ! loop counter for w
      real(8) :: ulo, uhi   ! Lab angles for input in legendre
      real(8) :: flo, fhi   ! f(w) at low and high points
      real(8) :: interp     ! Interpolation factor between w points in fw
      real(8) :: onepawr2   ! (1+AWR)^2
      real(8) :: onepR2     ! ( 1 + R^2)
      real(8) :: inv2REin   ! 1/(2*R*Ein)
      real(8) :: dw         ! Bin spacing of w array
      integer :: l          ! Legendre moment index

      dw = w(2) - w(1)

      R = awr * sqrt((ONE + Q * (awr + ONE) / (awr * Ein)))
      onepawr2 = (ONE + awr)**2
      onepR2 = ONE + R * R
      inv2REin = 0.5_8 / (R * Ein)

      do g = 1, size(E_bins) - 1
        ! Get bounds of integration, starting with low
        wlo = (E_bins(g) * onepawr2 - Ein * onepR2) * inv2REin
        ! Check to make sure wlo is in bounds
        if (wlo < -ONE) then
          wlo = -ONE
        else if (wlo > ONE) then
          wlo = ONE
        end if
        ilo = int((wlo + ONE) / dw) + 1
        ! Repeat for high end
        whi = (E_bins(g + 1) * onepawr2 - Ein * onepR2) * inv2REin
        ! Check to make sure whi is in bounds
        if (whi < -ONE) then
          whi = -ONE
        else if (whi > ONE) then
          whi = ONE
        end if
        ihi = int((whi + ONE) / dw) + 1

        ! Now we can skip groups we do not need to consider.
        ! We will cycle if whi == wlo = - 1 since that means we have not yet reached
        ! the 'active' groups.  We will return if whi == wlo = 1 since that means
        ! we have already passed the `active` groups.
        if (wlo == whi) then
          if (wlo == -ONE) then
            cycle
          else if (wlo == ONE) then
            return
          end if
        end if

        ! Get corresponding values of fw at wlo and whi
        ! Do low
        if (ilo == size(w)) then
          flo = fw(size(w))
        else
          interp = (wlo - w(ilo)) / (w(ilo + 1) - w(ilo))
          flo = (ONE - interp) * fw(ilo) + interp * fw(ilo + 1)
        end if

        ! Do high
        if (ihi == size(w)) then
          fhi = fw(size(w))
        else
          interp = (whi - w(ihi)) / (w(ihi + 1) - w(ihi))
          fhi = (ONE - interp) * fw(ihi) + interp * fw(ihi + 1)
        end if

        ! Now we simply integrate f(w) * Pl(mu(w))dw between wlo and whi
        if (ilo /= ihi) then
          ! Integrate part of the moment from the low point to the index above
          ! the low point
          ulo = tolab(R, wlo)
          uhi = tolab(R, w(ilo + 1))
          do l = 1, order
            distro(l, g) = (w(ilo + 1) - wlo) * (flo * calc_pn(l - 1, ulo) + &
              fw(ilo + 1) * calc_pn(l - 1, uhi))
          end do
          ! Integrate the inbetween pts
          do iw = ilo + 1, ihi -1
            ulo = uhi
            uhi = tolab(R, w(iw + 1))
            do l = 1, order
              distro(l, g) = distro(l, g) + (w(iw + 1) - w(iw)) * &
                (fw(iw) * calc_pn(l - 1, ulo) + fw(iw + 1) * calc_pn(l - 1, uhi))
            end do
          end do
          ! Integrate part of the moment from the index below the high point to
          ! the high point
          ulo = uhi
          uhi = tolab(R, whi)
          do l = 1, order
            distro(l, g) = distro(l, g) + (whi - w(ihi)) * &
              (fw(ihi) * calc_pn(l - 1, ulo) + fhi * calc_pn(l - 1, uhi))
          end do
        else
          ! The points are all within the same bin, can get flo and fhi directly
          ulo = tolab(R, wlo)
          uhi = tolab(R, whi)
          do l = 1, order
            distro(l, g) = (whi - wlo) * &
              (flo * calc_pn(l - 1, ulo) + fhi * calc_pn(l - 1, uhi))
          end do
        end if

        ! Now perform multiplication by 1/2, the last step of trapezoidal rule
        distro(:, g) = 0.5_8 * distro(:, g)

      end do

    end subroutine integrate_file4_cm_leg

!===============================================================================
! INTEGRATE_FILE6_CM_LEG Integrated the Center-of-Mass combined energy/angle
! distribution over all outgoing groups.
!===============================================================================

    subroutine integrate_file6_cm_leg(fEmu, mu, Ein, awr, Eout, INTT, thispdf, &
                                      E_bins, order, distro)
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

      integer :: g, imu_c   ! Group index variable, mu point index
      integer :: imu        ! Lab mu index
      real(8) :: Eo_cm, Eo  ! CM Outgoing Energy, Lab outgoing energy
      real(8) :: dEo        ! change in outgoing energy points
      integer :: iE         ! Outgoing energy counter
      real(8) :: mu_c       ! cm angle
      real(8) :: mu_l_min   ! Minimum lab angle
      real(8) :: dmu        ! change in lab angle
      real(8) :: c          ! c from eq 359
      real(8) :: proby, f   ! fmu at mu_c, interpolant between fmu points
      real(8) :: integ      ! the integrand
      real(8) :: J          ! the jacobian
      real(8) :: Eo_lo, Eo_hi ! Laboratory Eout bounds
      real(8), allocatable :: fEl(:) ! Integrated (over mu) fEmu
      real(8)              :: pEo    ! probability of this Eo value
      integer              :: iEo    ! Index of Eout array of current Eo
      real(8)              :: fEo    ! interpolant for Eo on Eout grid
      integer :: g_lo, g_hi ! Lower and upper groups
      real(8), allocatable :: E_bnds(:)
      real(8) :: ap1inv     ! 1/(AWR+1)
      real(8), allocatable :: fmu(:), mu_l(:)
      real(8), allocatable :: pdf(:) ! Eout PDF
      real(8) :: deltamu    ! mu bin spacing

      deltamu = mu(2) - mu(1)

      ! First lets normalize the PDF
      allocate(pdf(size(Eout)))
      pdf = thispdf
      if (Eout(size(Eout)) == Eout(size(Eout) - 1)) then
      ! Deal with Zr-90 (others?) who have the same Eout points at end
        pdf(size(Eout) - 1) = ZERO
      end if

      allocate(fEl(order))
      allocate(fmu(size(mu)))
      allocate(mu_l(size(mu)))

      ap1inv = ONE / (awr + ONE)

      ! Set up lower and upper lab energy boundaries
      Eo_lo = Eout(1) + (Ein - TWO * (awr + ONE) * sqrt(Ein * Eout(1))) &
        * ap1inv * ap1inv
      Eo_lo = 1E-12_8
      Eo_hi = Eout(size(Eout)) + (Ein + TWO * (awr + ONE) * &
        sqrt(Ein * Eout(size(Eout)))) * ap1inv * ap1inv

      if (Eo_lo <= E_bins(1)) then
        g_lo = 1
      else if (Eo_lo >= E_bins(size(E_bins))) then
        return
      else
        g_lo = binary_search(E_bins, size(E_bins), Eo_lo)
      end if
      if (Eo_hi <= E_bins(1)) then
        return
      else if (Eo_hi >= E_bins(size(E_bins))) then
        g_hi = size(E_bins) - 1
        allocate(E_bnds(g_lo : g_hi + 1))
        E_bnds(g_lo) = Eo_lo
        E_bnds(g_lo + 1: g_hi) = E_bins(g_lo + 1: g_hi)
        E_bnds(g_hi + 1) = E_bins(g_hi)
      else
        g_hi = binary_search(E_bins, size(E_bins), Eo_hi)
        allocate(E_bnds(g_lo : g_hi + 1))
        E_bnds(g_lo) = Eo_lo
        E_bnds(g_lo + 1: g_hi) = E_bins(g_lo + 1: g_hi)
        E_bnds(g_hi + 1) = Eo_hi
      end if

      do g = g_lo, g_hi
        Eo = E_bnds(g)
        dEo = (E_bnds(g + 1) - E_bnds(g)) / real(NE_PER_GRP - 1, 8)
        Eo = Eo - dEo
        do iE = 1, NE_PER_GRP
          Eo = Eo + dEo
          fEl = ZERO
          fmu = ZERO
          c = ap1inv * sqrt(Ein / Eo)
          mu_l_min = (ONE + c * c - Eout(size(Eout)) / Eo) / (TWO * c)
          if (mu_l_min < -ONE) then
            mu_l_min = -ONE
          else if (abs(mu_l_min - ONE) < 1E-10_8) then
            mu_l_min = ONE
          else if (mu_l_min > ONE) then
            cycle
          end if
          dmu = (ONE - mu_l_min) / real(size(mu) - 1, 8)
          do imu = 1, size(mu)
            mu_l(imu) = mu_l_min + dmu * real(imu - 1, 8)
            Eo_cm = Eo * (ONE + c * c - TWO * c * mu_l(imu))
            if (Eo_cm <= ZERO) then
              cycle
            else if (Eo_cm <= Eout(1)) then
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
              if (Eout(iEo + 1) == Eout(iEo)) then
                ! Take care of Zr-90 (and others?) where Law 61 has two
                ! Eout datapoints which are equal; since these seem to always
                ! be at the end, we can just say that fEo = 0, pE0 = 0
                fEo = ZERO
                pEo = pdf(iEo)
              else
                fEo = (Eo_cm - Eout(iEo)) / (Eout(iEo + 1) - Eout(iEo))
                pEo = (ONE - fEo) * pdf(iEo) + fEo * pdf(iEo + 1)
              end if
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

            if (abs(mu_c - ONE) < 1E-10_8) then
              imu_c = size(mu) - 1
              f = ONE
            else
              imu_c = int((mu_c + ONE) / deltamu) + 1
              f = (mu_c - mu(imu_c)) / (mu(imu_c + 1) - mu(imu_c))
            end if
            ! Do lower Eout point
            proby = (ONE - fEo) * &
              ((ONE - f) * fEmu(imu_c, iEo) + f * fEmu(imu_c + 1, iEo))
            ! And now the upper
            proby = proby + fEo * &
              ((ONE - f) * fEmu(imu_c, iEo + 1) + f * fEmu(imu_c + 1, iEo + 1))
            integ = proby * J * pEo
            fmu(imu) = integ
          end do

          do imu = 1, size(mu) - 1
            fEl = fEl + &
              calc_int_pn_tablelin(order, mu_l(imu), mu_l(imu + 1), &
              fmu(imu), fmu(imu + 1))
          end do

          if ((iE /= 1) .and. (iE /= NE_PER_GRP)) then
            distro(:, g) = distro(:, g) + TWO * fEl(:)
          else
            distro(:, g) = distro(:, g) + fEl(:)
          end if
        end do
        distro(:, g) = distro(:, g) * dEo * 0.5_8
      end do

      ! Lets do a normalization:
      ! 40090.7*c is up to 35% from being normalized, oddly
      fEo = ZERO
      do g = g_lo, g_hi
        fEo = fEo + distro(1, g)
      end do
      if (fEo > ZERO) fEo = ONE / fEo
      do g = g_lo, g_hi
        distro(:, g) =  distro(:, g) * fEo
      end do

    end subroutine integrate_file6_cm_leg

!===============================================================================
! LAW9_SCATTER_LAB_LEG Finds the legendre moments of the incoming adistro and
! assigns the results to each of the outgoing energy groups (according to a
! law 9 evaporation spectrum) according to the probabilities.
!===============================================================================

    subroutine law9_scatter_lab_leg(fmu, edist, Ein, E_bins, mu, order, distro)

      real(8), intent(in)  :: fmu(:)      ! Angle distro to act on
      type(DistEnergy), pointer, intent(in) :: edist    ! My energy dist
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: mu(:)       ! fEmu angular grid
      integer, intent(in)  :: order       ! Number of moments to find
      real(8), intent(out) :: distro(:,:) ! Resultant integrated distribution

      integer :: g, NR, NE, lc, imu
      real(8) :: T, U, x, I, Egp1, Eg, pE_xfer

      pE_xfer = ZERO

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
      if (Ein - U <= ZERO) return
      do g = 1, size(E_bins) - 1
        Egp1 = E_bins(g + 1)
        Eg = E_bins(g)
        if (Egp1 > (Ein - U)) Egp1 = Ein - U
        if (Eg > (Ein - U)) Eg = Ein - U
        pE_xfer = (exp(-Egp1 / T) * (T + Egp1)) - (exp(-Eg / T) * (T + Eg))
        pE_xfer= -T * pE_xfer / I

        ! Now set the angular distributions according to pE_xfer
        !!! This can be done separately and not for every group....
        !!! Change this when you have some more time.
        !!! (No big rush - this routine is not used frequently so a speed-up of
        !!! any kind isn't worth much)
        do imu = 1, size(mu) - 1
          distro(:, g) = distro(:, g) + &
            calc_int_pn_tablelin(order, mu(imu), mu(imu + 1), &
              fmu(imu), fmu(imu + 1)) * pE_xfer
        end do
      end do

    end subroutine law9_scatter_lab_leg

!===============================================================================
! INTEGRATE_FILE6_LAB_LEG Finds Legendre moments of the energy-angle
! distribution in fEmu over each of the outgoing energy groups and stores the
! result in distro.
!===============================================================================

    subroutine integrate_file6_lab_leg(fEmu, mu, Eout, INTT, thispdf, &
      E_bins, order, distro)
      real(8), intent(in)    :: fEmu(:,:)     ! Energy-angle distro to act on
      real(8), intent(in)    :: mu(:)         ! fEmu angular grid
      real(8), intent(in)    :: Eout(:)       ! Outgoing energies
      integer, intent(in)    :: INTT          ! Interpolation type (Hist || Lin-Lin)
      real(8), intent(in)    :: thispdf(:)    ! Outgoing E-dist PDF from mySD
      real(8), intent(in)    :: E_bins(:)     ! Energy group boundaries
      integer, intent(in)    :: order         ! Number of moments to find
      real(8), intent(out)   :: distro(:,:)   ! Resultant integrated distro

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

      ! First lets normalize the PDF
      allocate(pdf(NEout))
      pdf = thispdf
      do iE = 1, NEout - 1
        pdf(iE) = thispdf(iE) * (Eout(iE + 1) - Eout(iE))
      end do
      if (Eout(size(Eout)) == Eout(size(Eout) - 1)) then
      ! Deal with Zr-90 (others?) who have the same Eout points at end
        pdf(size(Eout) - 1) = ZERO
      end if

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
! TOLAB Calculates the Lab angle u given the CM angle w and reduced mass R
!===============================================================================

  pure function tolab(R, w) result(u)
    real(8), intent(in) :: R ! Reduced Mass
    real(8), intent(in) :: w ! Center-of-Mass Angle
    real(8)             :: u ! Resultant Lab angle

    real(8) :: f ! Interpolant

    if (R > ONE) then
      u = (ONE + R * w) / sqrt(ONE + R * R + TWO * R * w)
    else if (R == ONE) then
      ! divide by zero error at w=-1, avoid this
      if (w == -ONE) then
        u = -ONE
      else
        u = (ONE + R * w) / sqrt(ONE + R * R + TWO * R * w)
      end if
    else ! R < ONE
      if (w < -R) then
        ! Avoid unphysical results for w=[-1,-R]:
        ! Assume a linear shape to u(w) from w=-1 to w=-R
        ! First find u(-R)
        u = sqrt(ONE - R * R)
        ! Now do the linear interpolation
        f = (w - (-ONE)) / (-R - ONE)
        u = (ONE - f) * (-ONE) + f * u
      else
        u = (ONE + R * w) / sqrt(ONE + R * R + TWO * R * w)
      end if
    end if

  end function tolab

!===============================================================================
! IS_VALID_SCATTER determines if a given MT number is that of a scattering event
!===============================================================================

  function is_valid_scatter(MT) result(scatter_event)

    integer, intent(in) :: MT ! Reaction channel
    logical             :: scatter_event

    scatter_event = .false.
    if ((MT == ELASTIC) .or. ((MT >= N_2ND) .and. (MT <= N_NC))) then
      if (MT /= N_FISSION .and. MT /= N_F .and. MT /= N_NF .and. &
          MT /= N_2NF .and. MT /= N_3NF) then
        scatter_event = .true.
      end if
    end if

  end function is_valid_scatter

!===============================================================================
! UNITBASE performs the whole unit base conversion.
!===============================================================================

  subroutine unitbase(this, Ein, iE, Eout, pdf, INTT, fEmu)
    class(ScattData), target, intent(in) :: this ! Working ScattData object
    real(8), intent(in) :: Ein      ! Incoming energy to interpolate on
    integer, intent(in) :: iE       ! incoming energy index (searched)
    real(8), allocatable, intent(out) :: Eout(:)
    real(8), allocatable, intent(out) :: pdf(:)
    integer,              intent(out) :: INTT
    real(8), allocatable, intent(out) :: fEmu(:,:)

    real(8), allocatable :: ub1(:), ub2(:)

    call cast_to_unitbase(this % Eouts(iE) % data, this % pdfs(iE) % data, &
                          this % cdfs(iE) % data, this % INTT(iE), ub1)
    call cast_to_unitbase(this % Eouts(iE+1) % data, this % pdfs(iE+1) % data, &
                          this % cdfs(iE+1) % data, this % INTT(iE+1), ub2)

    call interp_unitbase(Ein, ub1, this % Eouts(iE) % data, &
                         this % pdfs(iE) % data, this % INTT(iE), &
                         this % distro(iE) % data, this % E_grid(iE), &
                         ub2, this % Eouts(iE+1) % data, &
                         this % pdfs(iE+1) % data, this % INTT(iE+1), &
                         this % distro(iE+1) % data, this % E_grid(iE+1), &
                         Eout, pdf, INTT, fEmu)

  end subroutine unitbase


!===============================================================================
! CAST_TO_UNITBASE converts an outgoing energy distribution to a PDF and CDF
! which ranges from 0 to 1 instead of from Eo_lo to Eo_hi
!===============================================================================

  subroutine cast_to_unitbase(Eout, pdf, cdf, INTT, ub_grid)
    real(8), allocatable, intent(in)  :: Eout(:)
    real(8), allocatable, intent(in)  :: pdf(:)
    real(8), allocatable, intent(in)  :: cdf(:)
    integer,              intent(in)  :: INTT
    real(8), allocatable, intent(out) :: ub_grid(:)

    integer :: ilo  ! Low index to start from
    integer :: i, j ! Loop counters
    real(8) :: inv_dE ! Difference between high and low energy points for scaling
    real(8), allocatable :: ub_temp(:)

    ! Some of these Eout distributions have the first two points with 0 probability
    ! (and the first will be an energy of zero)
    ! Find out of this is the case so we can discard the first one
    ilo = 1
    !!! Skipping for now, interp_unitbase is consistent with skipping this.
    !if (INTT == HISTOGRAM) then
    !  if (pdf(1) == ZERO) then
    !    ilo = 2
    !  end if
    !else
    !  if (cdf(2) == ZERO) then
    !    ilo = 2
    !  end if
    !end if

    ! Set up our storage location
    allocate(ub_temp(size(Eout) - ilo + 1))
    ub_temp = ZERO

    inv_dE = ONE / (Eout(size(Eout)) - Eout(ilo))

    j = 1
    do i = ilo, size(Eout) - 1
      ub_temp(j) = (Eout(i) - Eout(ilo)) * inv_dE
      j = j + 1
    end do
    ub_temp(j) = ONE

    if (ub_temp(j - 1) == ONE) then
      allocate(ub_grid(size(ub_temp) - 1))
      ub_grid = ub_temp(1:size(ub_temp) - 1)
    else
      allocate(ub_grid(size(ub_temp)))
      ub_grid = ub_temp
    end if
    deallocate(ub_temp)


  end subroutine cast_to_unitbase

!===============================================================================
! INTERP_UNITBASE interpolates between two unit-base grids to produce one grid
! in energy space (vice unit-base [0,1] space)
!===============================================================================

  subroutine interp_unitbase(Ein, ub1, Eout1, pdf1, INTT1, fEmu1, Ei1, ub2, &
                             Eout2, pdf2, INTT2, fEmu2, Ei2, Eout, pdf, INTT, fEmu)
    real(8),              intent(in)  :: Ein
    real(8), allocatable, intent(in)  :: ub1(:)
    real(8), allocatable, intent(in)  :: Eout1(:)
    real(8), allocatable, intent(in)  :: pdf1(:)
    integer,              intent(in)  :: INTT1
    real(8), allocatable, intent(in)  :: fEmu1(:,:)
    real(8),              intent(in)  :: Ei1
    real(8), allocatable, intent(in)  :: ub2(:)
    real(8), allocatable, intent(in)  :: Eout2(:)
    real(8), allocatable, intent(in)  :: pdf2(:)
    integer,              intent(in)  :: INTT2
    real(8), allocatable, intent(in)  :: fEmu2(:,:)
    real(8),              intent(in)  :: Ei2
    real(8), allocatable, intent(out) :: Eout(:)
    real(8), allocatable, intent(out) :: pdf(:)
    integer,              intent(out) :: INTT
    real(8), allocatable, intent(out) :: fEmu(:,:)

    real(8), allocatable :: ub(:)
    real(8) :: p1, p2 ! Value of pdf for 1 and 2 at ub point of interest
    real(8) :: dE1, dE2
    real(8) :: f, r
    integer :: i, j, k

    ! This routine has to create a new pdf grid with an Eout grid that has been adjusted

    ! First merge the two ub grids
    call merge(ub1, ub2, ub)
    allocate(Eout(size(ub)))

    Eout = ZERO
    allocate(pdf(size(ub)))
    pdf = ZERO
    allocate(fEmu(size(fEmu1(:,1)), size(ub)))
    fEmu = ZERO

    ! Find our interpolant
    f = (Ein - Ei1) / (Ei2 - Ei1)

    ! Find the energy widths of the two sets
    dE1 = (Eout1(size(Eout1)) - Eout1(1))
    dE2 = (Eout2(size(Eout2)) - Eout2(1))
    ! Now step through ub, and interpolate to find new pdf values
    do i = 1, size(ub)
      ! Find p1 (via interpolation)
      j = binary_search(ub1, size(ub1), ub(i))
      if (INTT1 == HISTOGRAM) then
        r = ZERO
      else if (INTT1 == LINEAR_LINEAR .or. INTT1 == LOG_LINEAR) then
        r = (ub(i) - ub1(j)) / (ub1(j+1) - ub1(j))
      else if (INTT1 == LINEAR_LOG .or. INTT1 == LOG_LOG) then
        r = log(ub(i) / ub1(j)) / log(ub1(j+1) / ub1(j))
      end if
      ! Now calculate the pdf at ub(1) for data-set1
      if (INTT1 == HISTOGRAM .or. INTT1 == LINEAR_LINEAR .or. &
          INTT1 == LINEAR_LOG) then
        p1 = (ONE - r) * pdf1(j) + r * pdf1(j+1)
      else if (INTT1 == LOG_LINEAR .or. INTT1 == LOG_LOG) then
        p1 = exp((ONE - r) * log(pdf1(j)) + r * log(pdf1(j + 1)))
      end if
      ! And add the portion of our interpolant for this distrib while we have r and j
      do k = 1, size(fEmu1(:,1))
        fEmu(k,i) = (ONE - f) * ((ONE - r) * fEmu1(k,j) + r * fEmu1(k,j+1))
      end do

      ! Repeat above for data set 2
      j = binary_search(ub2, size(ub2), ub(i))
      if (INTT1 == HISTOGRAM) then
        r = ZERO
      else if (INTT1 == LINEAR_LINEAR .or. INTT1 == LOG_LINEAR) then
        r = (ub(i) - ub2(j)) / (ub2(j+1) - ub2(j))
      else if (INTT1 == LINEAR_LOG .or. INTT1 == LOG_LOG) then
        r = log(ub(i) / ub2(j)) / log(ub2(j+1) / ub2(j))
      end if
      ! Now calculate the pdf at ub(1) for data-set1
      if (INTT1 == HISTOGRAM .or. INTT1 == LINEAR_LINEAR .or. &
          INTT1 == LINEAR_LOG) then
        p2 = (ONE - r) * pdf2(j) + r * pdf2(j+1)
      else if (INTT1 == LOG_LINEAR .or. INTT1 == LOG_LOG) then
        p2 = exp((ONE - r) * log(pdf2(j)) + r * log(pdf2(j + 1)))
      end if
      ! And add the portion of our interpolant for this distrib while we have r and j
      do k = 1, size(fEmu1(:,1))
        fEmu(k,i) = fEmu(k,i) + f * ((ONE - r) * fEmu2(k,j) + r * fEmu2(k,j+1))
      end do

      ! Now we have our p1 and p2, lets interpolate (with f now) to combine.
      ! I'm not sure how else to do this besides LIN-LIN
      pdf(i) = (ONE - f) * p1 + f * p2

      ! Now lets set what the Eout value actually is for this pdf
      Eout(i) = (ONE - f) * (Eout1(1) + dE1 * ub(i)) + &
        f * (Eout2(1) + dE2 * ub(i))

    end do

    ! Since we used LIN-LIN for interpolating, lets set that here too for
    ! consistency with future routines.
    INTT = LINEAR_LINEAR
  end subroutine interp_unitbase

end module scattdata_header
