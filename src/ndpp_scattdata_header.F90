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

! Module constants
  integer, parameter :: DISTRO_PTS = 2001

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
    real(8), allocatable :: mu(:)    ! mu pts to find f(mu) values at
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
    
    subroutine scatt_init(this, nuc, rxn, edist, E_bins, scatt_type, order)
      class(ScattData), intent(inout)       :: this  ! The object to initialize
      type(Nuclide), intent(in), pointer    :: nuc   ! Nuclide we are working on
      type(Reaction), intent(in), pointer   :: rxn   ! The reaction of interest
      type(Distenergy), pointer, intent(in) :: edist ! The energy distribution to use
      ! Edist is intended to specify which of the nested distros we are 
      ! actually using. t can be null.
      real(8), target, intent(in) :: E_bins(:)  ! Energy group boundaries
      integer, intent(in) :: scatt_type ! Type of format to store the data
      integer, intent(in) :: order      ! Order of the data storage format
      
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
      if (rxn % has_angle_dist) this % adist => rxn % adist
      this % edist => edist ! We dont use rxn % edist because this allows 
                            ! the ScattData type to be used for nested distros.
      
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
          allocate(this % distro(i) % data(DISTRO_PTS, NP))
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
          lc = int(edist % data(2 + 2*NR + int(edist % data(2 + 2*NR)) + i))
          NP   = int(edist % data(lc + 2))
          allocate(this % distro(i) % data(DISTRO_PTS, NP))
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
          this % E_grid(1) = rxn % threshold
        else
          this % E_grid(1) = E_bins(1)
        end if
        allocate(this % distro(this % NE))
        NP = 1
        do i = 1, this % NE
          ! There is no Energy-out dependence to the distribution, so we only
          ! need a 1 in the Eout dimension.
          allocate(this % distro(i) % data(DISTRO_PTS, NP))
          this % distro(i) % data = ZERO
        end do
      end if
      
      ! Assign the mu points
      allocate(this % mu(DISTRO_PTS))
      dmu = TWO / real(DISTRO_PTS - 1, 8)
      do i = 1, DISTRO_PTS - 1
        this % mu(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point is exactly ONE
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

    function scatt_interp_distro(this, mu_out, nuc, iE, Ein) result(distro)
      class(ScattData), target, intent(in)  :: this ! Working ScattData object
      real(8), intent(in), pointer          :: mu_out(:) ! The tabular output mu grid
      type(Nuclide), intent(in), pointer    :: nuc  ! Working nuclide
      integer, intent(in)           :: iE   ! incoming energy index
                                            ! (lower bound if Ein is provided)
      real(8), intent(in), optional :: Ein  ! Incoming energy to interpolate on
      
      type(Reaction), pointer   :: rxn   ! The reaction of interest
      real(8), allocatable :: distro(:,:)
      real(8) :: f, p_valid, sigS, Ein2
      integer :: global_iE
      real(8), pointer :: sigma(:) => null()
      real(8), pointer :: my_distro(:, :) => null()
      integer, allocatable :: mu_bins(:, :)
      real(8), allocatable :: mu_interp(:, :)
      
      rxn => this % rxn
      allocate(distro(this % order, this % groups))
      if (rxn % MT == ELASTIC) then
        sigma => nuc % elastic
      else
        sigma => rxn % sigma
      end if
      
      ! Get sigS - this code is pulled out of the big if(present) loop since both
      ! branches need it.
      if (present(Ein)) then
        Ein2 = Ein
      else
        Ein2 = this % E_grid(iE)
      end if
      if (Ein2 < nuc % energy(rxn % threshold)) then
        ! This is a catch-all, our energy was below the threshold, distro
        ! should be set to zero and we shall just exit this function
        distro = ZERO
        return
      else if (Ein2 <= nuc % energy(1)) then
        ! We are below the lowest global energy value so take the first entry
        sigS = sigma(rxn % threshold)
      else if (Ein2 >= nuc % energy(nuc % n_grid)) then
        ! We are above the global energy grid, so take the highest value
        sigS = sigma(nuc % n_grid)
      else
        global_iE = binary_search(nuc % energy, nuc % n_grid, &
          Ein2)
        ! check for rare case where two energy points are the same
        if (nuc % energy(global_iE) == nuc % energy(global_iE + 1)) &
          global_iE = global_iE + 1
        ! calculate interpolation factor
        f = (Ein2 - nuc % energy(global_iE))/ &
          (nuc % energy(global_iE + 1) - nuc % energy(global_iE))
        ! Adjust global_iE to point to the sigma array
        global_iE = global_iE - rxn % threshold + 1
        sigS = (ONE - f) * sigma(global_iE) + f * sigma(global_iE + 1)
      end if
      
      if (present(Ein)) then
        ! Get the probability value
        if (associated(this % edist)) then
          p_valid = interpolate_tab1(this % edist % p_valid, Ein)
        else
          p_valid = ONE
        end if
        
        ! Interpolate the distribution
        !!! For now we will just `interpolate' based on whichever
        !!! distribution is the closest to the requested energy
        global_iE = binary_search(this % E_grid(iE : this % NE), &
          this % NE - iE + 1, Ein)
        f = (Ein - this % E_grid(global_iE))/ &
          (this % E_grid(global_iE + 1) - this % E_grid(global_iE))
        if (f < 0.5_8) then
          my_distro => this % distro(global_iE) % data
        else
          my_distro => this % distro(global_iE + 1) % data
        end if
        
      else
        ! No interpolation is needed just multiply and return
        if (associated(this % edist)) then
          p_valid = interpolate_tab1(this % edist % p_valid, this % E_grid(iE))
        else
          p_valid = ONE
        end if
        
        my_distro => this % distro(global_iE) % data
      end if
      
      ! We know which distribution to work with, now it is time to:
      ! 1) convert from CM to Lab, if necessary
      ! 2) calculate the angular boundaries for integration
      ! 3) integrate according to scatt_type
      
      ! 1) convert from CM to Lab, if necessary
      call cm2lab(this % awr, this % rxn % Q_value, Ein2, this % mu, my_distro)
      
      ! 2) calculate the angular boundaries for integration
      allocate(mu_bins(2, this % groups))
      allocate(mu_interp(2, this % groups))
      call calc_mu_bins(this % awr, this % rxn % Q_value, Ein2, this % E_bins, &
        this % mu, mu_interp, mu_bins)
        
      ! 3) integrate according to scatt_type
      if (this % scatt_type == SCATT_TYPE_LEGENDRE) then
        call integrate_energyangle_leg(my_distro, this % mu, mu_interp, &
          mu_bins, this % order, distro)
      else if (this % scatt_type == SCATT_TYPE_TABULAR) then
        call integrate_energyangle_tab(my_distro, this % mu, mu_interp, &
          mu_bins, mu_out, distro)
      end if
      
      ! Combine the results
      !!! Where does multiplicity come in to play??
      distro = sigS * p_valid * distro
      
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
                  !!! If we are at the start, then use the next data point.
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
      integer,  intent(inout)               :: INTT     ! Energy out INTT grid
      real(8), intent(inout) :: distro(:,:) ! resultant distro (pts x NEout)
      
      real(8), pointer :: data(:) => null() ! Shorthand for adist % data
      integer :: lcin, lc     ! Locations inside edist % data
      integer :: idata, idata_prev ! Loop counters
      integer :: imu          ! mu loop indices
      integer :: iEout, NEout ! outgoing loop indices and number of points
      integer :: interp, NP   ! Tabular format data (interp type and # pts)
      real(8) :: r            ! Interpolation parameter
      integer :: NR, NE       ! edist % data navigation information
      real(8) :: KMR, KMA, KMconst ! Various data for Kalbach-Mann
      
      data => edist % data  
      
      NR = int(data(1))
      NE = int(data(2 + 2*NR))
      
      lc = int(data(2 + 2 * NR + NE + iE)) ! start of LDAT for iE
      NEout = int(data(lc + 2))
      allocate(Eouts(NEout))
      Eouts = data(lc + 2 + 1 : lc + 2 + NP)
      
      ! determine type of interpolation
      INTT = int(edist % data(lc + 1))
      if (INTT > 10) INTT = mod(INTT,10)
      
      if (edist % law  == 44) then
        lc = lc + 2
        do iEout = 1, NEout
          KMR = data(lc + 3 * NEout + iEout)
          KMA = data(lc + 4 * NEout + iEout)
          KMconst = 0.5_8 * KMA / sinh(KMA)
          distro(:, iEout) = KMconst * (cosh(KMA * mu(:)) + KMR * sinh(KMA * mu(:)))
        end do
      else if (edist % law  == 61) then
        lcin = lc + 2
        do iEout = 1, NEout
          lc = int(data(lcin + 3*NP + iEout))
          ! Check if isotropic
          if (lc == 0) then
            distro(:, iEout) = 0.5_8
            cycle
          end if
          
          interp = int(data(lc + 1))
          NP = int(data(lc + 2))
          lc = lc + 3
          
          ! Calculate a PDF value for each mu pt.
          if (interp == HISTOGRAM) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev , lc + NP -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, take value at mu(imu - 1)
                  distro(imu, iEout) = data(idata + NP - 1)
                  idata_prev = idata
                  exit
                else if (abs(data(idata) - mu(imu)) <= FP_PRECISION) then
                  ! Found a match, take value at mu(imu)
                  distro(imu, iEout) = data(idata + NP)
                  idata_prev = idata
                  exit
                end if
              end do
            end do
          else if (interp == LINEAR_LINEAR) then
            idata_prev = lc
            do imu = 1, size(mu)
              do idata = idata_prev, lc + NP -1
                if ((data(idata) - mu(imu)) > FP_PRECISION) then
                  ! Found a match, interpolate value
                  r = (mu(imu) - data(idata -1)) / (data(idata) - data(idata-1))
                  distro(imu, iEout) = data(idata + NP - 1) + r * &
                    (data(idata + NP) - data(idata - 1 + NP))
                  idata_prev = idata
                  exit
                else if (abs(data(idata) - mu(imu)) <= FP_PRECISION) then
                  ! Found a match, take value at mu(imu)
                  distro(imu, iEout) = data(idata + NP)
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

    pure subroutine cm2lab(awr, Q, Ein, mu, data)
      real(8), intent(in) :: awr   ! Atomic Weight Ratio for this nuclide
      real(8), intent(in) :: Q     ! Binding Energy of reaction, for finding R
      real(8), intent(in) :: Ein   ! Incoming energy
      real(8), intent(in) :: mu(:) ! Angular grid
      real(8), intent(inout) :: data(:,:) ! The distribution to convert
      
      real(8) :: R, mu_min ! Reduced Mass and minimum mu value for integration
      
      real(8) :: mu_cm(DISTRO_PTS) ! CM angular points corresponding to 
                                           ! mu(:), if mu was in Lab.
      real(8) :: tempsqrt
      integer :: imu, iEout  ! mu index, outgoing energy index
      
      ! Calculate the effective mass ratio
      if (Q /= ZERO) then
        R = awr * sqrt(ONE - (awr + ONE) * Q / (awr * Ein))
      else 
        R = awr
      end if
      
      ! Calculate the lower bounds of integration (-1 for everything besides
      ! H-1)
      if (R < ONE) then 
        mu_min = sqrt(ONE - R * R)
      else
        mu_min = -ONE
      end if
      
      ! Calculate the lab (mu) and CM (mu_cm) mu grid points.
      do imu = 1, DISTRO_PTS
        mu_cm(imu) = (mu(imu) * mu(imu) - ONE + mu(imu) * sqrt(R * R + &
          mu(imu) * mu(imu) - ONE)) / R
      end do
      
      do iEout = 1, size(data, dim = 2)
        ! Convert the CM distro to the laboratory system
        do imu = 1, DISTRO_PTS
          if (mu(imu) >= mu_min) then
            tempsqrt = sqrt(mu(imu) * mu(imu) + R * R - ONE)
            data(imu, iEout) = data(imu, iEout) * (TWO * mu(imu) + tempsqrt + &
              mu(imu) * mu(imu) / tempsqrt) / R
          end if
        end do
      end do
      
    end subroutine cm2lab

!===============================================================================
! CALC_MU_BINS Calculates the mu-values corresponding to the energy-out
! bins for a given input energy.
!===============================================================================

    subroutine calc_mu_bins(awr, Q, Ein, E_bins, mu, mu_interp, mu_bins)
      real(8), intent(in)  :: awr            ! Atomic-weight ratio
      real(8), intent(in)  :: Q              ! Reaction Q-Value
      real(8), intent(in)  :: Ein            ! Incoming energy
      real(8), intent(in)  :: E_bins(:)      ! Energy group boundaries
      real(8), intent(in)  :: mu(:)          ! tabular mu values
      real(8), intent(out) :: mu_interp(:,:) ! Outgoing mu values corresponding
                                             ! to the energy groups
      integer, intent(out) :: mu_bins(:,:)   ! Outgoing mu values corresponding
                                             ! to the energy groups
      
      real(8) :: mu_low, mu_high  ! Low and high angular points
      real(8) :: R            ! The Reduced Mass (takes in to account Qval)
      integer :: g            ! Group index variable
      integer :: imu
      
      ! From equation 234 in Methods for Processing ENDF/B-VII (pg 2798)
      R = sqrt(awr * awr  * (ONE + Q * (awr + ONE) / (awr * Ein)))
      
      do g = 1, size(E_bins) - 1
        ! Calculate the values of mu corresponding to this energy group
        mu_low = 0.5_8 * ((R + ONE) * sqrt(E_bins(g)/ Ein) - &
          (R - ONE) * sqrt(Ein / E_bins(g)))
        mu_high = 0.5_8 * ((R + ONE) * sqrt(E_bins(g + 1)/ Ein) - &
          (R - ONE) * sqrt(Ein / E_bins(g + 1)))
        
        ! Find the index in mu corresponding to mu_low
        if (mu_low < -ONE) then
          mu_bins(MU_LO, g) = 1
          mu_interp(MU_LO, g) = ZERO
        else if (mu_low > ONE) then
          mu_bins(MU_LO, g) = size(mu) - 1
          mu_interp(MU_LO, g) = ONE
        else
          imu = binary_search(mu, size(mu), mu_low)
          mu_bins(MU_LO, g) = imu
          mu_interp(MU_LO, g) = (mu_low - mu(imu - 1)) / (mu(imu) - mu(imu - 1))
        end if
        
        ! Find the index in mu corresponding to mu_high
        if (mu_high < -ONE) then
          mu_bins(MU_HI, g) = 1
          mu_interp(MU_HI, g) = ZERO
        else if (mu_high > ONE) then
          mu_bins(MU_HI, g) = size(mu) - 1
          mu_interp(MU_HI, g) = ONE
        else
          imu = binary_search(mu, size(mu), mu_high)
          mu_bins(MU_HI, g) = imu
          mu_interp(MU_HI, g) = (mu_high - mu(imu - 1)) / (mu(imu) - mu(imu - 1))
        end if
      end do
    end subroutine calc_mu_bins

!===============================================================================
! INTEGRATE_ENERGYANGLE_LEG Finds Legendre moments of the energy-angle 
! distribution in fEmu over each of the outgoing energy groups and stores the
! result in distro.
!===============================================================================

    subroutine integrate_energyangle_leg(fEmu, mu, mu_interp, mu_bins, order, &
      distro)
      
      real(8), intent(in)  :: fEmu(:,:)      ! Energy-angle distro to act on
      real(8), intent(in)  :: mu(:)          ! fEmu angular grid
      real(8), intent(in)  :: mu_interp(:,:) ! Interpolants on the group bounds
      integer, intent(in)  :: mu_bins(:,:)   ! indices of fEmu corresponding to
                                             ! the group boundaries
      integer, intent(in)  :: order          ! Number of moments to find
      real(8), intent(out) :: distro(:,:)    ! Resultant integrated distribution
      
      integer :: g                           ! outgoing energy group index
      
      ! Loop through energy groups
      do g = 1, size(distro, dim = 2)
        
      end do
      
    end subroutine integrate_energyangle_leg

!===============================================================================
! INTEGRATE_ENERGYANGLE_TAB Finds the tabular distribution representation of the
! energy-angle distribution in fEmu over each of the outgoing energy groups and 
! stores the result in distro.
!===============================================================================

    subroutine integrate_energyangle_tab(fEmu, mu, mu_interp, mu_bins, mu_out, &
      distro)
      
      real(8), intent(in)  :: fEmu(:,:)      ! Energy-angle distro to act on      
      real(8), intent(in)  :: mu(:)          ! fEmu angular grid
      real(8), intent(in)  :: mu_interp(:,:) ! Interpolants on the group bounds
      integer, intent(in)  :: mu_bins(:,:)   ! indices of fEmu corresponding to
                                             ! the group boundaries
      real(8), intent(in)  :: mu_out(:)      ! Number of moments to find
      real(8), intent(out) :: distro(:,:)    ! Resultant integrated distribution
      
      
      
    end subroutine integrate_energyangle_tab

end module scattdata_header
