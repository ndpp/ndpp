module ndpp_scatt
  
  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global,           only: nuclides, message
  use interpolation,    only: interpolate_tab1
  use output,           only: write_message, header, print_ascii_array
  use scattdata_header, only: scattdata
  use search,           only: binary_search
  use string,           only: to_str

  implicit none
  
  contains
 
!===============================================================================
! CALC_SCATT Calculates the group-to-group transfer matrices and scattering
! moments for the current nuclide
!===============================================================================
    
    subroutine calc_scatt(nuc, energy_bins, scatt_type, order, &
                          scatt_union, E_union, thin_tol)
      type(Nuclide), pointer, intent(in) :: nuc            ! Nuclide
      real(8), intent(in)                :: energy_bins(:) ! Energy groups
      integer, intent(in)                :: scatt_type     ! Scattering output type
      integer, intent(in)                :: order          ! Scattering data order
      real(8), allocatable,intent(inout) :: scatt_union(:,:,:) ! Unionized Scattering Matrices
      real(8), allocatable,intent(inout) :: E_union(:)     ! Unionized Energy Grid
      real(8), intent(in)                :: thin_tol       ! Thinning tolerance
      
      type(DistEnergy), pointer :: edist
      type(Reaction),   pointer :: rxn
      integer :: num_tot_rxn
      integer :: i_rxn, i_nested_rxn, imu
      integer :: iE
      real(8) :: max_err, dmu
      type(ScattData), allocatable, target  :: rxn_data(:)
      type(ScattData), pointer :: mySD
      real(8), pointer  :: mu_out(:) ! The tabular output mu grid
      
      ! This routine will parse through each nuc % reaction entry.
      ! For each, it will determine if the rxn is a scattering reaction, and if
      ! so it will set up the required memory spaces (the scatt_data object)
      ! and then perform the distribution calculations, storing the results in
      ! scatt_data.
      ! We will also keep a tally of which reactions were scatters to simplify 
      ! the collection of results
      
      ! First we have to find the total number of reactions, including the
      ! nested edists, so we can allocate rxn_data correctly.
      
      num_tot_rxn = 0
      do i_rxn = 1, nuc % n_reaction
        rxn => nuc % reactions(i_rxn)
        num_tot_rxn = num_tot_rxn + 1
        if (rxn % has_energy_dist) then
          edist => rxn % edist
          do while (associated(edist % next))
            edist => edist % next
            num_tot_rxn = num_tot_rxn + 1
          end do
        end if
      end do
      
      ! Size rxn_data according to num_tot_rxn
      allocate(rxn_data(num_tot_rxn))
      
      ! Allocate all the rxn_data objects with the correct energy distro, if
      ! needed.
      i_nested_rxn = 0
      do i_rxn = 1, nuc % n_reaction
        i_nested_rxn = i_nested_rxn + 1
        
        rxn => nuc % reactions(i_rxn)
        mySD => rxn_data(i_nested_rxn)
        edist => rxn % edist
        call mySD % init(nuc, rxn, edist, energy_bins, scatt_type, order)
        if (rxn % has_energy_dist) then
          do while (associated(edist % next))
            edist => edist % next
            i_nested_rxn = i_nested_rxn + 1
            mySD => rxn_data(i_nested_rxn)
            call mySD % init(nuc, rxn, edist, energy_bins, scatt_type, order)
          end do
        end if
      end do  
      
      do i_rxn = 1, num_tot_rxn
        mySD => rxn_data(i_rxn)
        ! Convert the angular distributions from the ACE data to Tabular format
        call mySD % convert_distro()
      end do
      
      ! Combine the reactions to a union grid
      if (scatt_type == SCATT_TYPE_TABULAR) then
        allocate(mu_out(order))
        dmu = TWO / real(order - 1, 8)
        do imu = 1, order
          mu_out(imu) = -ONE + real(imu - 1, 8) * dmu
        end do
      end if
      call union_grid(nuc, mu_out, rxn_data, E_union, scatt_union)

    end subroutine calc_scatt
 
!===============================================================================
! UNION_GRID Combines all the scattering data points on to one single energy
! grid, (the union grid).
!===============================================================================
    
    subroutine union_grid(nuc, mu_out, rxn_data, E_union, scatt_union)
      type(Nuclide), pointer, intent(in)   :: nuc   ! The nuclide of interest
      real(8), intent(inout), pointer      :: mu_out(:) ! The tabular output mu grid
      type(ScattData), intent(in), target  :: rxn_data(:)
      real(8), allocatable, intent(out)    :: E_union(:)
      real(8), allocatable, intent(out)    :: scatt_union(:,:,:)
      
      real(8), allocatable :: E_temp(:)
      real(8), allocatable :: scatt_temp(:,:,:)
      integer, allocatable :: scatt_loc(:)      ! Location of last-used energy
      logical, allocatable :: scatt_done(:)     ! Flag to state whether a rxn's grid is complete
      logical, allocatable :: scatt_hits(:)     ! Flag for when value is @ the lowest E
      integer :: iE, NE, maxNE, rxnE
      integer :: irxn
      integer :: groups, order
      type(ScattData), pointer :: mySD => NULL()
      
      groups = rxn_data(1) % groups
      order = rxn_data(1) % order
      
      ! Find the maximum number of energy points to allocate scatt_temp and E_temp
      maxNE = 0
      do irxn = 1, size(rxn_data)
        maxNE = maxNE + rxn_data(irxn) % NE
      end do
      allocate(E_temp(maxNE))
      E_temp = ZERO
      allocate(scatt_temp(order, groups, maxNE))
      scatt_temp = ZERO
      
      ! Allocate the searching variables, scatt_loc, _done, _hits
      allocate(scatt_loc(size(rxn_data)))
      scatt_loc = 1
      allocate(scatt_done(size(rxn_data)))
      scatt_done = .false.
      allocate(scatt_hits(size(rxn_data)))
      scatt_hits = .false.
      
      ! Grab the next energy data points, put them in the grid, interpolate the rest
      do iE = 1, maxNE
        ! Find the reactions at the smallest Energy points, store in scatt_hits
        call min_value_locs(rxn_data, scatt_loc, scatt_done, scatt_hits)
        ! Get the hits first since they involve no interpolation
        do irxn = 1, size(rxn_data)
          if (.not. scatt_done(irxn)) then
            if (scatt_hits(irxn)) then
              rxnE = scatt_loc(irxn)
              mySD => rxn_data(irxn)
              ! Add the energy point to the union grid
              E_temp(iE) = mySD % E_grid(rxnE)
              ! Add the scattering distribution to the union scattering grid
!~               scatt_temp(:,:,iE) = scatt_temp(:,:,iE) + &
!~                 mySD % interp_distro(mu_out, nuc, rxnE)
              ! Increment the location counter and set to done if the counter is
              ! out of bounds
              scatt_loc(irxn) = scatt_loc(irxn) + 1
              if (scatt_loc(irxn) > mySD % NE) scatt_done(irxn) = .true. 
            end if
          end if
        end do
        
        ! Now interpolate those which are not the next data point
        do irxn = 1, size(rxn_data)
          if (.not. scatt_done(irxn)) then
            if(.not. scatt_hits(irxn)) then
              rxnE = scatt_loc(irxn)
              mySD => rxn_data(irxn)
              ! If the iE energy point is outside the range of this reaction,
              ! don't score it
              if (E_temp(iE) < mySD % E_grid(1)) cycle
              ! Add the scattering distribution to the union scattering grid
!~               scatt_temp(:,:,iE) = scatt_temp(:,:,iE) + &
!~                 mySD % interp_distro(mu_out, nuc, rxnE, E_temp(iE))
            end if
          end if
        end do
        
        ! Check to see if all reactions have reached their last data point
        ! If so, then we can complete the unionizing
        if (all(scatt_done)) exit
      end do
      
      ! Allocate and assign the final union grid spaces
      allocate(scatt_union(order, groups, iE))
      allocate(E_union(iE))
      scatt_union(:,:,:) = scatt_temp(:,:, 1 : iE)
      E_union(:) = E_temp(1 : iE)
      
      deallocate(scatt_temp)
      deallocate(E_temp)
      deallocate(scatt_loc)
      deallocate(scatt_done)
      deallocate(scatt_hits)
      
    end subroutine union_grid

!===============================================================================
! MIN_VALUE_LOCS Finds the reaction in rxn_data with the next highest Ein, and
! stores the result in hits
!===============================================================================
  
    subroutine min_value_locs(rxn_data, loc, done, hits)
      type(ScattData), intent(in)  :: rxn_data(:) ! Rxn information
      integer, intent(in)    :: loc(:)  ! Location of last-used energy
      logical, intent(in)    :: done(:) ! Flag to state whether a rxn's grid is complete
      logical, intent(inout) :: hits(:) ! The reactions which are the lowest E
      
      integer :: i
      real(8) :: min_val
      
      hits = .false.
      min_val = INFINITY
      
      ! Find the next smallest energy value
      do i = 1, size(rxn_data)
        if (.not. done(i)) then
          if (rxn_data(i) % E_grid(loc(i)) < min_val) &
            min_val = rxn_data(i) % E_grid(loc(i))
        end if
      end do
      
      ! Find all locs  that are within a certain tolerance of the min value
      do i = 1, size(rxn_data)
        if (.not. done(i)) then
          if (abs(rxn_data(i) % E_grid(loc(i)) - min_val) < FP_PRECISION) then
            hits(i) = .true.
          end if
        end if
      end do
    
    end subroutine min_value_locs  
    
end module ndpp_scatt
