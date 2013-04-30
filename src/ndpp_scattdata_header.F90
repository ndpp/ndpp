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
  implicit none

! Module constants
  integer, parameter :: CM_ORDER = 20

!===============================================================================
! SCATT_DATA Stores the data for each reaction of the nuclide in question.
!===============================================================================
  
  type :: ScattData
    logical :: is_init = .false.            ! Initialization status
    integer              :: NE = 0          ! Number of Ein values
    real(8), allocatable :: E_grid(:)       ! Ein values
    real(8), allocatable :: distro(:,:,:)   ! Output distribution
                                            ! # orders x g_out x # E points
    real(8), allocatable :: mu_bins(:,:)    ! mu_min and mu_max for each group
    real(8), pointer     :: E_bins(:)       ! Energy group boundaries from input
    integer              :: scatt_type = -1 ! Type of format to store the data
    integer              :: order      =  0 ! Order of the data storage format
    integer              :: groups     =  0 ! Number of outgoing energy groups      
    real(8)              :: awr        =  ZERO   ! atomic weight ratio
    type(Reaction),   pointer :: rxn   => NULL() ! My reaction
    type(DistEnergy), pointer :: edist => NULL() ! My reaction's combined dist
    type(DistAngle),  pointer :: adist => NULL() ! My reaction's angle dist
    
    ! Type-Bound procedures
    contains
      procedure :: init  => scatt_init  ! Sets NE, allocates spaces, etc
      procedure :: clear => scatt_clear ! Deallocates this object
      procedure :: calc_distro => scatt_calc_distro ! Calculate the distribution
      procedure :: interp_distro => scatt_interp_distro ! Interpolates the distro
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
      
      integer :: NR ! Number of interpolation regions (from ACE)
      integer :: lc ! ACE data locator variable
      integer :: iE ! energy grid point index of the scattering data
      integer :: global_iE ! energy grid index on the nnuclide's global data
      real(8) :: interp_val ! Interpolated cross-section value
      real(8) :: f          ! Interpolation factor
      
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
      else if (associated(this % edist)) then
        NR = int(edist % data(1))
        this % NE = int(edist % data(2 + 2*NR))
        allocate(this % E_grid(this % NE))
        lc = 2 + 2*NR
        this % E_grid(1 : this % NE) = edist % data(lc + 1 : lc + this % NE)
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
      end if
      
      ! Allocate the group-dependent attributes
      this % E_bins => E_bins
      this % groups = size(E_bins) - 1
      allocate(this % distro(CM_ORDER, this % groups, this % NE))
      allocate(this % mu_bins(2,this % groups))
      
      ! The final initialization. If we made it here then this bad boy was 
      ! successful.
      this % is_init = .true.
    end subroutine scatt_init
    
!===============================================================================
! SCATT_CLEAR Clears (deallocates) the scatt_data object.
!===============================================================================

    subroutine scatt_clear(this)
      class(ScattData), intent(inout) :: this ! The object to clear
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
        deallocate(this % distro)
        deallocate(this % mu_bins)
      end if
      
      this % is_init = .false.
    end subroutine scatt_clear
    
!===============================================================================
! SCATT_CALC_DISTRO Controls the flow of execution necessary to calculate the 
! angular distribution as requested, and store it in this % distro.
!===============================================================================

    subroutine scatt_calc_distro(this)
      class(ScattData), intent(inout) :: this ! The object to act on
      
      integer :: iE ! Ein grid index
      
      ! Step through each incoming energy value to do these calculations
      do iE = 1, this % NE
        ! Calculate the minimum mu and maximum mu for each energy group
        ! BUT only do this if we do not have a combined energy/angle
        ! distribution!  IF we do have that case, then mu_bins are -1 and 1, 
        ! per the Methods for Processing ENDF guide.
        if (associated(this % edist)) then
          if (this % edist % law /= 44 .and. this % edist % law /= 61) then
            call calc_mu_bins(this % E_grid(iE), this % awr, &
              this % rxn % Q_value, this % E_bins, this % mu_bins)
          else
            this % mu_bins = ONE
          end if
        else
          call calc_mu_bins(this % E_grid(iE), this % awr, &
            this % rxn % Q_value, this % E_bins, this % mu_bins)
        end if
        
        ! Calculate the energy-angle integrals (where necessary)
        ! of the outgoing distro
        if (this % rxn % scatter_in_cm) then
          ! If we are in center-of-mass, store and then convert
          ! Set order to 20 so that the conversion from CM to LAB has a minimal
          ! loss in accuracy. (The lower orders of Lab depend on the higher
          ! orders of CM)
          this % distro(:, :, iE) = integrate_distro(this % scatt_type, &
            20, this % mu_bins, iE, this % adist, this % edist)
        else
          ! If we are in lab, just store the result
          ! Since we arent converting the moments later to Lab, we only need to
          ! get this % order number of orders.
          this % distro(:, :, iE) = integrate_distro(this % scatt_type, &
            this % order, this % mu_bins, iE, this % adist, this % edist)
        end if
        
      end do
      
      if (this % rxn % scatter_in_cm) then
        ! If we are in center-of-mass then convert to lab
        if (this % scatt_type == SCATT_TYPE_LEGENDRE) then
          call cm2lab_legendre_numeric(this % awr, this % rxn % Q_value, &
            this % E_grid, this % order, CM_ORDER, this % distro)
        else
          write(*,*) 'Not yet implemented.'
        end if
      end if
      
    end subroutine scatt_calc_distro

!===============================================================================
! SCATT_INTERP_DISTRO calculates the value of the distribution (times its 
! probability) at the optionally given incoming energy.  If no incoming energy
! is given, then no interpolation is performed, and instead we can just use the
! data directly from the given index (iE)
!===============================================================================

    function scatt_interp_distro(this, nuc, iE, Ein) result(distro)
      class(ScattData), intent(in)  :: this ! The ScattData object to work with
      type(Nuclide), intent(in), pointer :: nuc   ! Nuclide we are working on
      integer, intent(in)           :: iE   ! incoming energy index
                                            ! (lower bound if Ein is provided)
      real(8), intent(in), optional :: Ein  ! Incoming energy to interpolate on
      
      type(Reaction), pointer   :: rxn   ! The reaction of interest
      real(8), allocatable :: distro(:,:)
      real(8) :: f, p_valid, sigS, Ein2
      integer :: global_iE
      real(8), pointer :: sigma(:) => NULL()
      
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
        ! This is a catch-all, our energy was below the threshold, so set the 
        ! probability to 0
        sigS = ZERO
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
          distro = this % distro(:,:, global_iE)
        else
          distro = this % distro(:,:, global_iE + 1)
        end if
        
      else
        ! No interpolation is needed just multiply and return
        if (associated(this % edist)) then
          p_valid = interpolate_tab1(this % edist % p_valid, this % E_grid(iE))
        else
          p_valid = ONE
        end if
        
        distro = this % distro(:, :, iE)
      end if
      
      ! Combine the results
      distro = sigS * p_valid * distro
      
    end function scatt_interp_distro


!===============================================================================
!===============================================================================
! FUNCTIONS TO SUPPORT SCATTDATA
!===============================================================================
!===============================================================================

!===============================================================================
! CALC_mu_bins Calculates the mu-values corresponding to the energy-out
! bins for a given input energy.
!===============================================================================

    pure subroutine calc_mu_bins(Ein, awr, Q, E_bins, mu_bins)
      real(8), intent(in)  :: Ein            ! Incoming energy
      real(8), intent(in)  :: awr            ! Atomic-weight ratio
      real(8), intent(in)  :: Q              ! Reaction Q-Value
      real(8), intent(in)  :: E_bins(:)      ! Energy group boundaries
      real(8), intent(out) :: mu_bins(:,:) ! Outgoing mu values corresponding
                                             ! to the energy groups
      
      real(8) :: R            ! The Reduced Mass (takes in to account Qval)
      integer :: g            ! Group index variable
      
      ! From equation 234 in Methods for Processing ENDF/B-VII (pg 2798)
      R = sqrt(awr * awr  * (ONE + Q * (awr + ONE) / (awr * Ein)))
      
      do g = 1, size(E_bins) - 1
        mu_bins(MU_LO,g) = 0.5_8 * ((R + ONE) * sqrt(E_bins(g)/ Ein) - &
          (R - ONE) * sqrt(Ein / E_bins(g)))
        mu_bins(MU_HI,g) = 0.5_8 * ((R + ONE) * sqrt(E_bins(g + 1)/ Ein) - &
          (R - ONE) * sqrt(Ein / E_bins(g + 1)))
        if (mu_bins(MU_LO,g) < -ONE) mu_bins(MU_LO,g) = -ONE
        if (mu_bins(MU_HI,g) < -ONE) mu_bins(MU_HI,g) = -ONE
        if (mu_bins(MU_LO,g) >  ONE) mu_bins(MU_LO,g) =  ONE
        if (mu_bins(MU_HI,g) >  ONE) mu_bins(MU_HI,g) =  ONE
      end do     
mu_bins(MU_LO,:) = -ONE      
mu_bins(MU_HI,:) =  ONE
    end subroutine calc_mu_bins

!===============================================================================
! INTEGRATE_DISTRO Performs the integration of the scattering distribution for
! the current energy point for all outgoing energy groups.
!===============================================================================

    function integrate_distro(scatt_type, order, mu_bins, iE, &
      adist, edist) result(distro)
      integer, intent(in)  :: scatt_type   ! Type of format to store the data
      integer, intent(in)  :: order        ! Order of the data storage format
      real(8), intent(in)  :: mu_bins(:,:) ! mu_min and mu_max for each group
      integer, intent(in)  :: iE           ! Energy index to act on
      type(DistAngle), pointer, intent(in)   :: adist ! My angle dist
      type(DistEnergy), pointer , intent(in) :: edist ! My combined dist
      
      real(8), allocatable :: distro(:,:) ! resultant distro (order x groups)
      
      allocate(distro(order, size(mu_bins, dim = 2)))
      distro = ZERO

      ! This routine needs to:
      ! 1) Determine which type of distro we have to deal with
      ! 2) Depending on the type of law and output type, 
      ! integrate over just angle or energy and angle
      
      select case(scatt_type)
        case (SCATT_TYPE_LEGENDRE)
          if (associated(edist)) then
            if (edist % law == 44 .or. edist % law == 61) then
              ! combined energy/angle distribution - have to integrate over
              ! the energy part of the data, AND the angle part.
!~               call integrate_file6_legendre()
            else
              ! angle distribution only. Only have to integrate the angle.
              call integrate_file4_legendre(order, mu_bins, iE, adist, distro)
            end if
          else if (associated(adist)) then
            ! angle distribution only. Only have to integrate the angle.
            call integrate_file4_legendre(order, mu_bins, iE, adist, distro)
          end if 
        case (SCATT_TYPE_TABULAR)
          if (associated( edist)) then
            if (edist % law == 44 .or. edist % law == 61) then
              ! combined energy/angle distribution - have to integrate over
              ! the energy part of the data, AND the angle part.
            else
              ! angle distribution only. Only have to integrate the angle.
              
            end if
          else if (associated(adist)) then
            ! angle distribution only. Only have to integrate the angle.
          end if 
      end select
      
      
      
    end function integrate_distro

end module scattdata_header
