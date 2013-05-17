program test_scattdata
  use ace_header
  use constants
  use scattdata_header
  
  implicit none

!===============================================================================
  
  write(*,*) '***************************************************'
  write(*,*) 'Testing ScattData'
  write(*,*) '***************************************************'
  
  ! Test clear routine
  call test_clear()
  
  ! Test initialization routine
  call test_init()
  
  write(*,*)
  write(*,*) '***************************************************'
  write(*,*) 'Testing ScattData Passed!'
  write(*,*) '***************************************************'

!===============================================================================  
  contains
  
!===============================================================================
! TEST_CLEAR Tests the clearing of the scattdata object.
!===============================================================================

    subroutine test_clear()
      type(ScattData)  :: mySD           ! Testing object
!~       type(Nuclide), target     :: nuc   ! Nuclide we are working on
      type(Reaction), target    :: rxn   ! The reaction of interest
      type(DistAngle),  pointer :: adist ! The angular distribution to use
      type(DistEnergy), pointer :: edist ! The energy distribution to use
      real(8), allocatable, target  :: E_bins(:)  ! Energy group boundaries
      integer :: scatt_type ! Type of format to store the data
      integer :: order      ! Order of the data storage format
      integer :: mu_bins    ! Number of angular pts in tabular rep.
      
      integer :: i          ! Loop counter
            
      write(*,*) 
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing ScattData % clear'
      
      ! The members of scattData are all public. Manually set these to some
      ! value and ensure clear sets them back the way they should be.
      mySD % is_init = .true.
      mySD % NE = 2
      allocate(mySD % E_grid(2))
      allocate(mySD % distro(2))
      allocate(mySD % Eouts(2))
      do i = 1, 2
        allocate(mySD % distro(i) % data(2,2))
        allocate(mySD % Eouts(i) % data(2))
      end do
      allocate(mySD % INTT(2))
      allocate(mySD % mu(2))
      allocate(E_bins(2))
      mySD % E_bins => E_bins
      mySD % scatt_type = 2
      mySD % order = 2
      mySD % groups = 2
      mySD % awr = TWO
      mySD % rxn => rxn
      mySD % edist => edist
      mySD % adist => adist
      
      ! Now test that everything is cleared after calling clear
      call mySD % clear()
      if (mySD % is_init .or. mySD % NE /= 0 .or. allocated(mySD % E_grid) .or. &
        allocated(mySD % distro) .or. allocated(mySD % Eouts) .or. &
        allocated(mySD % INTT) .or. allocated(mySD % mu) .or. &
        associated(mySD % E_bins) .or. mySD % scatt_type /= -1 .or. &
        mySD % order /= 0 .or. mySD % groups /= 0 .or. mySD % awr /= ZERO .or. &
        associated(mySD % rxn) .or. associated(mySD % edist) .or. &
        associated(mySD % adist)) then
        
        write(*,*) 'ScattData % clear FAILED!'
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'ScattData % clear Passed!'
      write(*,*) '---------------------------------------------------'
      
      deallocate(E_bins)
    end subroutine test_clear
  
!===============================================================================
! TEST_INIT Tests the initialization of the scattdata object. This test is 
! written assuming scattdata % clear() has already been tested and passes.
!===============================================================================    
    subroutine test_init()
      type(ScattData)  :: mySD           ! Testing object
      type(Nuclide), target     :: nuc   ! Nuclide we are working on
      type(Reaction), target    :: rxn   ! The reaction of interest
      type(Distenergy), pointer :: edist ! The energy distribution to use
      real(8), allocatable, target  :: E_bins(:)  ! Energy group boundaries
      integer :: scatt_type ! Type of format to store the data
      integer :: order      ! Order of the data storage format
      integer :: mu_bins    ! Number of angular pts in tabular rep.
      
      integer :: i          ! loop counters
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing ScattData % init'
      
      ! Set the init parameters to simple values
      rxn % MT = ELASTIC
      rxn % has_angle_dist = .true.
      rxn % adist % n_energy = 2
      allocate(rxn % adist % energy(2))
      rxn % adist % energy = (/1.0E-11_8, 20.0_8/)
      edist => null()
      nuc % awr = ONE
      allocate(E_bins(2))
      E_bins = rxn % adist % energy
      scatt_type = SCATT_TYPE_LEGENDRE
      order = 5
      mu_bins = 3
      
      ! The parameters currently set up will pass, call init and ensure 
      ! is_init is true and that the parameters are initialized as expected
      
      
      ! Next, test a situation with a rxn % MT which will cause an automatic
      ! return of init
      ! This includes all of:
      ! (n,fiss), (n,f), (n,nf), (n,2nf), (n,3nf), MT>=200, (n,level)
      do i = 1, 7
        select case (i)
          case (1) 
            rxn % MT = N_FISSION
          case (2) 
            rxn % MT = N_F
          case (3) 
            rxn % MT = N_NF
          case (4) 
            rxn % MT = N_2NF
          case (5) 
            rxn % MT = N_3NF
          case (6) 
            rxn % MT = 200
          case (7) 
            rxn % MT = N_LEVEL
        end select
        call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)
        if (mySD % is_init) then
          write(*,*) 'ScattData % init FAILED! (CASE = ',i,')'
          stop 10
        end if
        call mySD % clear()
      end do
      ! Restore initial conditions
      rxn % MT = ELASTIC
      
      
      
      write(*,*)
      write(*,*) 'ScattData % init Passed!'
      write(*,*) '---------------------------------------------------'
    end subroutine test_init
    
end program test_scattdata

