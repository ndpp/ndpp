module scattdata_header
  
  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global,           only: nuclides, message
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
      type(Nuclide), intent(in)    :: nuc   ! Nuclide we are working on
      type(Reaction), target, intent(in)   :: rxn   ! The reaction of interest
      type(Distenergy), pointer, intent(in) :: edist ! The energy distribution to use
      ! Edist is intended to specify which of the nested distros we are 
      ! actually using. t can be null.
      real(8), target, intent(in) :: E_bins(:)  ! Energy group boundaries
      integer, intent(in) :: scatt_type ! Type of format to store the data
      integer, intent(in) :: order      ! Order of the data storage format
      integer, intent(in) :: mu_bins    ! Number of angular pts in tabular rep.
      
      integer :: i          ! loop counter
      real(8) :: dmu        ! mu spacing
      integer :: NP, NR     ! Number of outgoing energy values and number of interp regions
      integer :: lc         ! Location of outgoing energy data in edist % data
      
      ! Test reactions to ensure we have a scattering reaction.
      ! These tests will leave this as uninitialized (and is_init == .false.),
      ! which will be used as a flag when the reactions are combined.
      if (rxn % MT == N_FISSION .or. rxn % MT == N_F .or. rxn % MT == N_NF &
        .or. rxn % MT == N_2NF .or. rxn % MT == N_3NF) return

      ! some materials have gas production cross sections with MT > 200 that
      ! are duplicates. Also MT=4 is total level inelastic scattering which
      ! should be skipped.
      if (rxn % MT >= 200 .or. rxn % MT == N_LEVEL) return
      
      ! We survived the above check and thus have a scattering reaction.
      
      ! Save the order and scattering type
      this % scatt_type = scatt_type
      this % order = order
      ! Store the reaction type
      this % rxn => rxn
      ! Store the atomic weight ratio
      this % awr = nuc % awr
      
      ! Set distributions
      if (rxn % has_angle_dist) then
        this % adist => rxn % adist
        this % edist => null()
      else if (associated(edist)) then
        this % adist => null()
        this % edist => edist ! We dont use rxn % edist because this allows 
                              ! ScattData type to be used for nested distros.
      else
        this % adist => null()
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
      else if (associated(this % edist)) then
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
      else
        ! Must be isotropic. Just use the top and bottom of the energy bins
        this % NE = 2
        allocate(this % E_grid(this % NE))
        ! Set the upper value
        this % E_grid(this % NE) = E_bins(size(E_bins))
        ! Set the lower value; but, if the threshold of this reaction is above 
        ! the lower bound, use that instead
        if (rxn % threshold > E_bins(1)) then
          this % E_grid(1) = nuc % energy(rxn % threshold)
        else
          this % E_grid(1) = E_bins(1)
        end if
        allocate(this % distro(this % NE))
        NP = 1
        do i = 1, this % NE
          ! There is no Energy-out dependence to the distribution, so we only
          ! need a 1 in the Eout dimension.
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
      
      ! Allocate the Eouts jagged container and the interpolation type array
      allocate(this % Eouts(this % NE))
      allocate(this % INTT(this % NE))
      
      ! Allocate the group-dependent attributes
      this % E_bins => E_bins
      this % groups = size(E_bins) - 1
      
      ! The final initialization. If we made it here then this bad boy was 
      ! successful.
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
      this % awr        = ZERO
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
        end do
        deallocate(this % distro)
        deallocate(this % Eouts)
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
      
      ! Step through each incoming energy value to do these calculations
      do iE = 1, this % NE
        if (associated(this % edist)) then
          ! combined energy/angle distribution - have to integrate over
          ! the energy part of the data, AND the angle part.
          call convert_file6(iE, this % mu, this % edist, &
            this % Eouts(iE) % data, this % INTT(iE), this % distro(iE) % data)
        else if (associated(this % adist)) then
          ! angle distribution only. Only have to integrate the angle.
          call convert_file4(iE, this % mu, this % adist, &
            this % Eouts(iE) % data, this % INTT(iE), &
            this % distro(iE) % data(:, 1))
        end if 
      end do
      
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
      real(8) :: f, p_valid, sigS          ! interpolation, probability of this
                                           ! (nested) reaction, and \sigma_s(E)
      integer :: iE                        ! incoming energy index (searched)
                                           ! (lower bound if Ein is provided)
      integer :: nuc_iE                    ! Energy index on the nuclide's x/s
      real(8), pointer :: sigS_array(:) => null()   ! sigS pointer
      real(8), pointer :: my_distro(:, :) => null() ! the distribution chosen
                                                    ! by the binary search on Ein
                                                    
      integer :: bins(2, this % groups)   ! Start/end energy/angle bin indices
      real(8) :: vals(2, this % groups)   ! values on start/end energy/angle
      real(8) :: interp(2, this % groups) ! interpolation on start/end energy/angle
      
      ! Set up the result memory 
      if (this % scatt_type == SCATT_TYPE_LEGENDRE) then
        allocate(distro(this % order + 1, this % groups))
      else
        allocate(distro(this % order, this % groups))
      end if
      
      ! Set rxn, so we can save some characters throughout this function.
      rxn => this % rxn
      
      ! Point the cross-section pointer to the right place
      if (rxn % MT == ELASTIC) then
        sigS_array => nuc % elastic
      else
        sigS_array => rxn % sigma
      end if
      
      ! Get sigS - this code is pulled out of the big if(present) loop since both
      ! branches need it.
      if (Ein < nuc % energy(rxn % threshold)) then
        ! This is a catch-all, our energy was below the threshold, distro
        ! should be set to zero and we shall just exit this function
        distro = ZERO
        return
      else if (Ein <= nuc % energy(1)) then
        ! We are below the lowest global energy value so take the first entry
        sigS = sigS_array(rxn % threshold)
      else if (Ein >= nuc % energy(nuc % n_grid)) then
        ! We are above the global energy grid, so take the highest value
        sigS = sigS_array(nuc % n_grid)
      else
        nuc_iE = binary_search(nuc % energy, nuc % n_grid, &
          Ein)
        ! check for rare case where two energy points are the same
        if (nuc % energy(nuc_iE) == nuc % energy(nuc_iE + 1)) &
          nuc_iE = nuc_iE + 1
        ! calculate interpolation factor
        f = (Ein - nuc % energy(nuc_iE))/ &
          (nuc % energy(nuc_iE + 1) - nuc % energy(nuc_iE))
        ! Adjust nuc_iE to point to the sigS_array array
        nuc_iE = nuc_iE - rxn % threshold + 1
        sigS = (ONE - f) * sigS_array(nuc_iE) + f * sigS_array(nuc_iE + 1)
      end if
      
      ! Get the probability value
      if (associated(this % edist)) then
        p_valid = interpolate_tab1(this % edist % p_valid, Ein)
      else
        p_valid = ONE
      end if
      
      ! Search on the angular distribution's energy grid to find what energy
      ! index Ein is at.
      iE = binary_search(this % E_grid, this % NE, Ein)
      
      ! Interpolate the distribution
      !!! For now we will just `interpolate' based on whichever
      !!! distribution is the closest to the requested energy
      f = (Ein - this % E_grid(iE))/ &
        (this % E_grid(iE + 1) - this % E_grid(iE))
      if (f >= 0.5_8) iE = iE + 1
      my_distro => this % distro(iE) % data
      
      ! We know which distribution to work with, now it is time to:
      ! 1) convert from CM to Lab, if necessary
      ! 2) calculate the angular boundaries for integration
      ! 3) calculate the energy boundaries for integration
      ! 4) integrate according to scatt_type
      
      ! 1) convert from CM to Lab, if necessary
      if (rxn % scatter_in_cm) &
        call cm2lab(this % awr, this % rxn % Q_value, Ein, this % mu, my_distro)
      
      select case (this % scatt_type)
        case (SCATT_TYPE_LEGENDRE)
          if (associated(this % adist)) then
            ! 2) calculate the angular boundaries for integration
            call calc_mu_bounds(this % awr, this % rxn % Q_value, Ein, this % E_bins, &
              this % mu, interp, vals, bins)
            ! 4) integrate according to scatt_type
            call integrate_energyangle_file4_leg(my_distro(:, 1), this % mu, &
              interp, vals, bins, this % order, distro)
          else if (associated(this % edist)) then
            ! 3) calculate the energy boundaries for integration
            call calc_E_bounds(this % E_bins, this % Eouts(iE) % data, &
              this % INTT(iE), interp, bins)
            ! 4) integrate according to scatt_type
            call integrate_energyangle_file6_leg(my_distro, this % mu, &
              this % Eouts(iE) % data, this % E_bins, interp, bins, &
              this % order, distro)
          end if
        case (SCATT_TYPE_TABULAR)
          if (associated(this % adist)) then
            ! 2) calculate the angular boundaries for integration
            call calc_mu_bounds(this % awr, this % rxn % Q_value, Ein, this % E_bins, &
              this % mu, interp, vals, bins)
            ! 4) integrate according to scatt_type
!~             call integrate_energyangle_file4_tab(my_distro(:, 1), this % mu, &
!~               interp, vals, bins, this % order, distro)
          else if (associated(this % edist)) then
            ! 3) calculate the energy boundaries for integration
            call calc_E_bounds(this % E_bins, this % Eouts(iE) % data, &
              this % INTT(iE), interp, bins)
            ! 4) integrate according to scatt_type
            
          end if
      end select
      
      ! Combine the results
      distro = sigS * p_valid * distro * rxn % multiplicity
      
    end function scatt_interp_distro

!===============================================================================
!===============================================================================
! FUNCTIONS TO SUPPORT SCATTDATA
!===============================================================================
!===============================================================================

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
      
      ! Check what type of distribution we have
      lc   = adist % location(iE)
      
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

    subroutine convert_file6(iE, mu, edist, Eouts, INTT, distro)
      integer, intent(in)  :: iE            ! Energy index to act on
      real(8), intent(in)  :: mu(:)         ! tabular mu points
      type(DistEnergy), pointer, intent(in) :: edist    ! My energy dist
      real(8), allocatable, intent(inout)   :: Eouts(:) ! Energy out grid @ Ein
      integer, intent(inout)                :: INTT     ! Energy out INTT grid
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
      
      if (edist % law  == 44) then
        lc = lc + 2
        do iEout = 1, NP
          KMR = data(lc + 3 * NP + iEout)
          KMA = data(lc + 4 * NP + iEout)
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
        end do
      end if
      
    end subroutine convert_file6

!===============================================================================
! CM2LAB Converts a tabular center-of-mass distribution to a laboratory frame 
! of reference tabular distribution.
!===============================================================================

    subroutine cm2lab(awr, Q, Ein, mu, data)
      real(8), intent(in) :: awr   ! Atomic Weight Ratio for this nuclide
      real(8), intent(in) :: Q     ! Binding Energy of reaction, for finding R
      real(8), intent(in) :: Ein   ! Incoming energy
      real(8), intent(in) :: mu(:) ! Angular grid
      real(8), intent(inout) :: data(:,:) ! The distribution to convert
      
      real(8) :: R, Rinv, R2       ! Reduced Effective Mass, 1/R, and R^2
      integer :: imu_min           ! Location in mu of the minimum feasible
                                   ! cm to lab mu point.
      
      real(8) :: mu_l(size(data, dim = 1))  ! CM angular points corresponding to 
                                            ! mu(:), if mu was in Lab.
      real(8) :: tempdistro(size(data, dim = 1)) ! Temporary storage of the 
                                                 ! converted distro
      real(8) :: tempsqrt
      integer :: imu, imu_l, iEout      ! mu indices, outgoing energy index
      real(8) :: interp                 ! Interpolation parameter
      integer :: mu_bins                ! Number of mu bins in dim 1 of data
      real(8) :: mu_min                 ! Minimum value of mu for R < 1
      
      mu_bins = size(data, dim = 1)
      
      ! From equation 234 in Methods for Processing ENDF/B-VII (pg 2798)
      R = sqrt(awr * awr  * (ONE + Q * (awr + ONE) / (awr * Ein)))

      Rinv = ONE / R
      R2   = R * R
        
      ! Calculate the lab (mu) and CM (mu_l) mu grid points.
      do imu = 1, mu_bins
        if ((R < 1) .and. (mu(imu) == -ONE)) then
          mu_l(imu) = sqrt(ONE - R2)
        else
          tempsqrt = sqrt(ONE + R2 + TWO * R * mu(imu))
          mu_l(imu) = (ONE + R * mu(imu)) / tempsqrt
        end if
      end do
      
      do iEout = 1, size(data, dim = 2)
        ! Convert the CM distro to the laboratory system
        do imu = 1, mu_bins
          if ((R < 1) .and. (mu(imu) == -ONE)) then
            ! Is this really zero??
            tempdistro(imu) = ZERO
          else
            tempsqrt = sqrt(mu_l(imu) * mu_l(imu) + R2 - ONE)
            tempdistro(imu) = data(imu, iEout) * (TWO * mu_l(imu) + tempsqrt + &
              mu_l(imu) * mu_l(imu) / tempsqrt) * Rinv
          end if
        end do
        
        ! Now put this back on to the original mu grid.
        do imu = 1 , mu_bins - 1
          if (mu(imu) < mu_l(1)) then
            imu_l = 1
          else
            imu_l = binary_search(mu_l, mu_bins, mu(imu))
          end if
          ! Get the interpolation parameter
          interp = (mu(imu) -  mu_l(imu_l)) / (mu_l(imu_l + 1) - mu_l(imu_l))
          data(imu, iEout) = tempdistro(imu_l) + interp * &
            (tempdistro(imu_l + 1) - tempdistro(imu_l))
        end do
        data(mu_bins, iEout) = tempdistro(mu_bins)
        
        do imu = 1, mu_bins
          if (data(imu, iEout) < ZERO) data(imu, iEout) = ZERO
        end do
      end do
            
    end subroutine cm2lab

!===============================================================================
! CALC_MU_BOUNDS Calculates the mu-values corresponding to the energy-out
! bins for a given input energy.
!===============================================================================

    subroutine calc_mu_bounds(awr, Q, Ein, E_bins, mu, interp, vals, bins)
      
      real(8), intent(in)  :: awr         ! Atomic-weight ratio
      real(8), intent(in)  :: Q           ! Reaction Q-Value
      real(8), intent(in)  :: Ein         ! Incoming energy
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: mu(:)       ! tabular mu values
      real(8), intent(out) :: interp(:,:) ! Outgoing mu interpolants
                                          ! corresponding to the energy groups
      real(8), intent(out) :: vals(:,:)   ! Outgoing mu values corresponding
                                          ! to the energy groups
      integer, intent(out) :: bins(:,:)   ! Outgoing mu indices corresponding
                                          ! to the energy groups
      
      real(8) :: mu_low, mu_high  ! Low and high angular points
      real(8) :: R            ! The Reduced Mass (takes in to account Qval)
      real(8) :: alpha        ! Energy-transfer constant
      integer :: g            ! Group index variable
      integer :: imu
      
      ! From equation 234 in Methods for Processing ENDF/B-VII (pg 2798)
      R = sqrt(awr * awr  * (ONE + Q * (awr + ONE) / (awr * Ein)))
      
      do g = 1, size(E_bins) - 1
        ! Calculate the values of mu corresponding to this energy group
        ! These come from eqs. 232-233 in Methods for Processing ENDF/B-VII, 
        ! also on pg 2798
        
        ! First check to see if the energy transfer is possible
        ! (the pg 2798 eqns dont work if R<1)
        ! If it isnt, then we set the boundary to -1
        alpha = (R - 1) / (R + 1)
        alpha = alpha * alpha
        if (E_bins(g) < alpha * Ein) then
          mu_low = -ONE
        else
          mu_low  = 0.5_8 * ((ONE + awr) * sqrt(E_bins(g) / Ein) + &
            (ONE - R * R) / (ONE + awr) * sqrt(Ein / E_bins(g)))
        end if
        if (E_bins(g + 1) < alpha * Ein) then
          mu_high = -ONE
        else
          mu_high = 0.5_8 * ((ONE + awr) * sqrt(E_bins(g + 1) / Ein) + &
            (ONE - R * R) / (ONE + awr) * sqrt(Ein / E_bins(g + 1)))
        end if
        
        ! Find the index in mu corresponding to mu_low
        if (mu_low < -ONE) then
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
    end subroutine calc_mu_bounds

!===============================================================================
! CALC_E_BOUNDS Calculates the energy points corresponding to the energy-out
! bins for a given input energy.
!===============================================================================

    subroutine calc_E_bounds(E_bins, Eout, INTT, interp, bins)
      real(8), intent(in)  :: E_bins(:)   ! Energy group boundaries
      real(8), intent(in)  :: Eout(:)     ! Output energies
      integer, intent(in)  :: INTT        ! Output energies interpolation type
      real(8), intent(out) :: interp(:,:) ! Outgoing E interpolants
                                          ! corresponding to the energy groups
      integer, intent(out) :: bins(:,:)   ! Outgoing E indices corresponding
                                          ! to the energy groups
      
      integer :: g                        ! Group index variable
      
      do g = 1, size(E_bins) - 1
        if (E_bins(g) < Eout(1)) then
          bins(MU_LO, g)   = 0
          interp(MU_LO, g) = ZERO
        else if (E_bins(g) > Eout(size(Eout))) then
          bins(MU_LO, g)   = 0
          interp(MU_LO, g) = ZERO
        else
          bins(MU_LO, g)   = binary_search(Eout, size(Eout), E_bins(g))
          interp(MU_LO, g) = (E_bins(g) - Eout(bins(MU_LO, g))) / & 
            (Eout(bins(MU_LO, g) + 1) - Eout(bins(MU_LO, g)))
        end if
        if (E_bins(g + 1) < Eout(1)) then
          bins(MU_HI, g)   = 0
          interp(MU_HI, g) = ZERO
        else if (E_bins(g + 1) > Eout(size(Eout))) then
          bins(MU_HI, g)   = 0
          interp(MU_HI, g) = ZERO
        else
          bins(MU_HI, g)   = binary_search(Eout, size(Eout), E_bins(g + 1))
          interp(MU_HI, g) = (E_bins(g + 1) - Eout(bins(MU_HI, g))) / & 
            (Eout(bins(MU_HI, g) + 1) - Eout(bins(MU_HI, g)))
        end if
        ! adjust the interpolant so that if using histogram it is correctly 0
        if (INTT == HISTOGRAM) then
          interp(MU_LO, g) = ZERO
          interp(MU_HI, g) = ZERO
        end if
      end do
      
    end subroutine calc_E_bounds

!===============================================================================
! INTEGRATE_ENERGYANGLE_*_LEG Finds Legendre moments of the energy-angle 
! distribution in fEmu over each of the outgoing energy groups and stores the
! result in distro. The FILE4 version does this for only an angular distribution,
! while the FILE6 version does the same for a combined energy-angle distribution
!===============================================================================

    subroutine integrate_energyangle_file4_leg(fEmu, mu, interp, vals, &
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
          distro(:, g) = ZERO
        end if
      end do
      
    end subroutine integrate_energyangle_file4_leg
    
    subroutine integrate_energyangle_file6_leg(fEmu, mu, Eout, E_bins, interp, &
      bins, order, distro)
      
      real(8), intent(in)  :: fEmu(:,:)      ! Energy-angle distro to act on
      real(8), intent(in)  :: mu(:)          ! fEmu angular grid
      real(8), intent(in)  :: Eout(:)        ! Outgoing energies
      real(8), intent(in)  :: E_bins(:)      ! Energy group boundaries
      real(8), intent(in)  :: interp(:,:)    ! interpolants of E values
      integer, intent(in)  :: bins(:,:)      ! indices of fEmu corresponding to
                                             ! the group boundaries
      integer, intent(in)  :: order          ! Number of moments to find
      real(8), intent(out) :: distro(:,:)    ! Resultant integrated distribution
      
      real(8), allocatable :: fEmu_int(:,:)  ! Integrated (over E) fEmu
      integer :: g         ! Energy group index
      integer :: imu       ! angular grid index
      integer :: iE        ! outgoing energy grid index
      real(8) :: Elo, Ehi  ! interpolated value of the energy
      real(8) :: flo, fhi  ! interpolated value of fEmu(imu,Eout)
      
      allocate(fEmu_int(size(mu), size(bins, dim = 2)))
      
      ! This routine will perform integration of the outgoing energy of fEmu
      ! over each energy group in E_bins. This will be done with trapezoidal
      ! integration !!!But can later be improved to something else, like what
      ! NJOY does that I dont quite understand yet!!!.  
      ! Trapezoidal integration = 1/2 * (b-a)*(f(b)+f(a))
      do g = 1, size(E_bins ) - 1
        ! Do the lower energyout distribution
        if (bins(MU_LO, g) /= 0) then
          Elo = E_bins(g)
          Ehi = Eout(bins(MU_LO, g))
          do imu = 1, size(mu)
            flo = fEmu(imu, bins(MU_LO, g)) + interp(MU_LO, g) * &
              (fEmu(imu, bins(MU_LO, g) + 1) - fEmu(imu, bins(MU_LO, g)))
            fhi = fEmu(imu, bins(MU_LO, g) + 1)
            fEmu_int(imu, g) = (Ehi - Elo) * (fhi + flo)
          end do
        end if
        ! Do the intermediate energyout distributions
        do iE = bins(MU_LO, g) + 1, bins(MU_HI, g) - 1
          Elo = Eout(iE)
          Ehi = Eout(iE + 1)
          do imu = 1, size(mu)
            flo = fEmu(imu, iE)
            fhi = fEmu(imu, iE + 1)
            fEmu_int(imu, g) = (Ehi - Elo) * (fhi + flo)
          end do
        end do
        ! Do the upper energyout distribution
        if (bins(MU_HI, g) /= 0) then
          Elo = Eout(bins(MU_HI, g))
          Ehi = E_bins(g + 1)
          do imu = 1, size(mu)
            flo = fEmu(imu, bins(MU_HI, g))
            fhi = fEmu(imu, bins(MU_HI, g)) + interp(MU_HI, g) * &
              (fEmu(imu, bins(MU_HI, g) + 1) - fEmu(imu, bins(MU_HI, g)))
            fEmu_int(imu, g) = (Ehi - Elo) * (fhi + flo)
          end do
        end if
        ! Apply the 1/2 term from trapezoidal integration
        fEmu_int(:, g) = 0.5_8 * fEmu_int(:, g)
        
        ! Perform the legendre expansion
        do imu = 1, size(mu)
          distro(:, g) = distro(:, g) + &
            calc_int_pn_tablelin(order, mu(imu), &
            mu(imu + 1), fEmu_int(imu, g), fEmu_int(imu + 1, g))
        end do
      end do
      
    end subroutine integrate_energyangle_file6_leg
    
!===============================================================================
! INTEGRATE_ENERGYANGLE_*_TAB Finds the tabular distribution representation of
! the energy-angle distrib in fEmu over each of the outgoing energy groups and 
! stores the result in distro. The FILE4 version does this for only an angular 
! distribution while the FILE6 version does the same for a combined energy-angle
! distribution
!===============================================================================
 

end module scattdata_header
