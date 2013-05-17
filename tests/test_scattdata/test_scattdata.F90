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
      type(Reaction), target    :: rxn   ! The reaction of interest
      type(DistAngle),  target :: adist ! The angular distribution to use
      type(DistEnergy), target :: edist ! The energy distribution to use
      real(8), allocatable, target  :: E_bins(:)  ! Energy group boundaries
      
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
      type(ScattData)           :: mySD  ! Testing object
      type(Nuclide)             :: nuc   ! Nuclide we are working on
      type(Reaction), target    :: rxn   ! The reaction of interest
      type(DistEnergy), pointer :: edist ! The energy distribution to use
      type(DistEnergy), target  :: myedist 
      real(8), allocatable, target  :: E_bins(:)  ! Energy group boundaries
      integer :: scatt_type ! Type of format to store the data
      integer :: order      ! Order of the data storage format
      integer :: mu_bins    ! Number of angular pts in tabular rep.
      
      integer :: i          ! loop counters
      integer :: NR, NE, NP ! edist data
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing ScattData % init'
      
      ! Set the init parameters to simple values
      rxn % MT = ELASTIC
      rxn % has_angle_dist = .true.
      rxn % adist % n_energy = 2
      allocate(rxn % adist % energy(2))
      rxn % adist % energy = (/2.0E-11_8, 20.0_8/)
      edist => null()
      nuc % awr = ONE
      allocate(nuc % energy(3))
      nuc % energy = (/1.0E-11_8, ONE, 20.0_8/)
      allocate(E_bins(2))
      E_bins = (/1.0E-11_8, 20.0_8/)
      scatt_type = SCATT_TYPE_LEGENDRE
      order = 5
      mu_bins = 3
      
      ! First test a situation with a rxn % MT which will cause an automatic
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
      write(*,*) 'Testing File 4 Distribution'
      ! Next test a file 4 distribution initialization (adist)
      call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)
      ! To reduce the scope of each if-block, we will check the scalars,
      ! then the allocatables, then the pointers
      ! First check the scalars
      if ((.not. mySD % is_init) .or. &
        mySD % scatt_type /= SCATT_TYPE_LEGENDRE .or. &
        mySD % order /= 5 .or. mySD % awr /= ONE .or. mySD % NE /= 2) then
        write(*,*) 'ScattData % init FAILED! (Scalars)'
        stop 10
      end if
      ! Allocatables, one at a time
      if (.not. allocated(mySD % E_grid)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Not Allocated)'
        stop 10
      else if (size(mySD % E_grid) /= 2) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Size)'
        stop 10
      else if ((mySD % E_grid(1) /= 2E-11_8) .or. &
        (mySD % E_grid(2) /= 20.0_8)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % INTT)) then
        write(*,*) 'ScattData % init FAILED! (INTT Not Allocated)'
        stop 10
      else if (size(mySD % INTT) /= 2) then
        write(*,*) 'ScattData % init FAILED! (INTT Incorrect Size)'
        stop 10
      end if
      if (.not. allocated(mySD % mu)) then
        write(*,*) 'ScattData % init FAILED! (mu Not Allocated)'
        stop 10
      else if (size(mySD % mu) /= 3) then
        write(*,*) 'ScattData % init FAILED! (mu Incorrect Size)'
        stop 10
      else if ((mySD % mu(1) /= -ONE) .or. (mySD % mu(2) /= ZERO) .or. &
        (mySD % mu(3) /= ONE)) then
        write(*,*) 'ScattData % init FAILED! (mu Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % distro)) then
        write(*,*) 'ScattData % init FAILED! (distro Not Allocated)'
        stop 10
      else if (size(mySD % distro) /= 2) then
        write(*,*) 'ScattData % init FAILED! (distro Incorrect Size)'
        stop 10
      else if ((size(mySD % distro(1) % data, dim = 1) /= 3) .or. &
        (size(mySD % distro(1) % data, dim = 2) /= 1) .or. &
        (size(mySD % distro(2) % data, dim = 1) /= 3) .or. &
        (size(mySD % distro(2) % data, dim = 2) /= 1)) then
        write(*,*) 'ScattData % init FAILED! (distro % data Incorrect Size)'
        stop 10
      else if ((any(mySD % distro(1) % data /= ZERO)) .or. &
        (any(mySD % distro(2) % data /= ZERO))) then
        write(*,*) 'ScattData % init FAILED! (distro % data Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % Eouts)) then
        write(*,*) 'ScattData % init FAILED! (Eouts Not Allocated)'
        stop 10
      else if (size(mySD % Eouts) /= 2) then
        write(*,*) 'ScattData % init FAILED! (Eouts Incorrect Size)'
        stop 10
      end if
      ! Pointers
      if ((.not. associated(mySD % rxn, rxn)) .or. &
        (.not. associated(mySD % adist, rxn % adist)) .or. &
        associated(mySD % edist) .or. &
        (.not. associated(mySD % E_bins, E_bins))) then
        write(*,*) 'ScattData % init FAILED! (Associated)'
        stop 10
      end if
      call mySD % clear()
      
      ! Test a file 6 distribution case (edist passed)
      write(*,*)
      write(*,*) 'Testing File 6 Distribution'
      rxn % has_angle_dist = .false.
      rxn % adist % n_energy = 0
      deallocate(rxn % adist % energy)
      ! Now set up myedist with two Eins, each with 2 Eouts
      NR = 0
      NE = 2
      NP = 2
      allocate(myedist % data(5+2*NR+NE + 1+5*NP + 1 + 2+5*NP - 1))
      myedist % data(1:2+2*NR+NE) = (/ZERO, TWO, 2.0E-11_8, 20.0_8/)  ! NR, NE, Ein(:) 
      myedist % data(3+2*NR+NE : 4+2*NR+NE) = (/4+2*NR+NE, 5+2*NR+NE + 1+5*NP/) ! lc(Ein)
      ! Set Law 44 data for first Ein
      ! INTT', NP, Eout(:), PDF(:), CDF(:), R(:), A(:)
      myedist % data(5+2*NR+NE : ) = &
        (/ONE, real(NP,8), 1.0E-11_8, 20.0_8, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
        ZERO, ZERO/) 
      ! Set Law 44 data for second Ein
      ! INTT', NP, Eout(:), PDF(:), CDF(:), R(:), A(:)
      myedist % data(5+2*NR+NE + 2+5*NP : ) = &
        (/ONE, real(NP,8), 1.0E-11_8, 20.0_8, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
        ZERO, ZERO/) 
      edist => myedist
      
      call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)
      ! To reduce the scope of each if-block, we will check the scalars,
      ! then the allocatables, then the pointers
      ! First check the scalars
      if ((.not. mySD % is_init) .or. &
        mySD % scatt_type /= SCATT_TYPE_LEGENDRE .or. &
        mySD % order /= 5 .or. mySD % awr /= ONE .or. mySD % NE /= 2) then
        write(*,*) 'ScattData % init FAILED! (Scalars)'
        stop 10
      end if
      ! Allocatables, one at a time
      if (.not. allocated(mySD % E_grid)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Not Allocated)'
        stop 10
      else if (size(mySD % E_grid) /= 2) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Size)'
        stop 10
      else if ((mySD % E_grid(1) /= 2E-11_8) .or. &
        (mySD % E_grid(2) /= 20.0_8)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % INTT)) then
        write(*,*) 'ScattData % init FAILED! (INTT Not Allocated)'
        stop 10
      else if (size(mySD % INTT) /= 2) then
        write(*,*) 'ScattData % init FAILED! (INTT Incorrect Size)'
        stop 10
      end if
      if (.not. allocated(mySD % mu)) then
        write(*,*) 'ScattData % init FAILED! (mu Not Allocated)'
        stop 10
      else if (size(mySD % mu) /= 3) then
        write(*,*) 'ScattData % init FAILED! (mu Incorrect Size)'
        stop 10
      else if ((mySD % mu(1) /= -ONE) .or. (mySD % mu(2) /= ZERO) .or. &
        (mySD % mu(3) /= ONE)) then
        write(*,*) 'ScattData % init FAILED! (mu Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % distro)) then
        write(*,*) 'ScattData % init FAILED! (distro Not Allocated)'
        stop 10
      else if (size(mySD % distro) /= 2) then
        write(*,*) 'ScattData % init FAILED! (distro Incorrect Size)'
        stop 10
      else if ((size(mySD % distro(1) % data, dim = 1) /= 3) .or. &
        (size(mySD % distro(1) % data, dim = 2) /= 2) .or. &
        (size(mySD % distro(2) % data, dim = 1) /= 3) .or. &
        (size(mySD % distro(2) % data, dim = 2) /= 2)) then
        write(*,*) 'ScattData % init FAILED! (distro % data Incorrect Size)'
        stop 10
      else if ((any(mySD % distro(1) % data /= ZERO)) .or. &
        (any(mySD % distro(2) % data /= ZERO))) then
        write(*,*) 'ScattData % init FAILED! (distro % data Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % Eouts)) then
        write(*,*) 'ScattData % init FAILED! (Eouts Not Allocated)'
        stop 10
      else if (size(mySD % Eouts) /= 2) then
        write(*,*) 'ScattData % init FAILED! (Eouts Incorrect Size)'
        stop 10
      end if
      ! Pointers
      if ((.not. associated(mySD % rxn, rxn)) .or. &
        (.not. associated(mySD % edist, myedist)) .or. &
        associated(mySD % adist) .or. &
        (.not. associated(mySD % E_bins, E_bins))) then
        write(*,*) 'ScattData % init FAILED! (Associated)'
        stop 10
      end if
      call mySD % clear()
      
      ! Test an isotropic distribution case (adist and edist dont exist)
      write(*,*)
      write(*,*) 'Testing Isotropic Distribution'
      edist => null()
      call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)
      ! To reduce the scope of each if-block, we will check the scalars,
      ! then the allocatables, then the pointers
      ! First check the scalars
      if ((.not. mySD % is_init) .or. &
        mySD % scatt_type /= SCATT_TYPE_LEGENDRE .or. &
        mySD % order /= 5 .or. mySD % awr /= ONE .or. mySD % NE /= 2) then
        write(*,*) 'ScattData % init FAILED! (Scalars)'
        stop 10
      end if
      ! Allocatables, one at a time
      if (.not. allocated(mySD % E_grid)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Not Allocated)'
        stop 10
      else if (size(mySD % E_grid) /= 2) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Size)'
        stop 10
      else if ((mySD % E_grid(1) /= 1E-11_8) .or. &
        (mySD % E_grid(2) /= 20.0_8)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % INTT)) then
        write(*,*) 'ScattData % init FAILED! (INTT Not Allocated)'
        stop 10
      else if (size(mySD % INTT) /= 2) then
        write(*,*) 'ScattData % init FAILED! (INTT Incorrect Size)'
        stop 10
      end if
      if (.not. allocated(mySD % mu)) then
        write(*,*) 'ScattData % init FAILED! (mu Not Allocated)'
        stop 10
      else if (size(mySD % mu) /= 3) then
        write(*,*) 'ScattData % init FAILED! (mu Incorrect Size)'
        stop 10
      else if ((mySD % mu(1) /= -ONE) .or. (mySD % mu(2) /= ZERO) .or. &
        (mySD % mu(3) /= ONE)) then
        write(*,*) 'ScattData % init FAILED! (mu Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % distro)) then
        write(*,*) 'ScattData % init FAILED! (distro Not Allocated)'
        stop 10
      else if (size(mySD % distro) /= 2) then
        write(*,*) 'ScattData % init FAILED! (distro Incorrect Size)'
        stop 10
      else if ((size(mySD % distro(1) % data, dim = 1) /= 3) .or. &
        (size(mySD % distro(1) % data, dim = 2) /= 1) .or. &
        (size(mySD % distro(2) % data, dim = 1) /= 3) .or. &
        (size(mySD % distro(2) % data, dim = 2) /= 1)) then
        write(*,*) 'ScattData % init FAILED! (distro % data Incorrect Size)'
        stop 10
      else if ((any(mySD % distro(1) % data /= ZERO)) .or. &
        (any(mySD % distro(2) % data /= ZERO))) then
        write(*,*) 'ScattData % init FAILED! (distro % data Incorrect Values)'
        stop 10
      end if
      if (.not. allocated(mySD % Eouts)) then
        write(*,*) 'ScattData % init FAILED! (Eouts Not Allocated)'
        stop 10
      else if (size(mySD % Eouts) /= 2) then
        write(*,*) 'ScattData % init FAILED! (Eouts Incorrect Size)'
        stop 10
      end if
      ! Pointers
      if ((.not. associated(mySD % rxn, rxn)) .or. &
        associated(mySD % edist) .or. associated(mySD % adist) .or. &
        (.not. associated(mySD % E_bins, E_bins))) then
        write(*,*) 'ScattData % init FAILED! (Associated)'
        stop 10
      end if
      call mySD % clear()
      
      ! Finally, repeat the isotropic case, but set a threshold energy
      ! above E_bins(1) so that the E_grid(1) == rxn % threshold
      write(*,*)
      write(*,*) 'Testing Isotropic Distribution w/ Threshold'
      rxn % threshold = 2
      call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)      
      ! To reduce the scope of each if-block, we will check the scalars,
      ! then the allocatables, then the pointers
      ! First check the scalars
      if (.not. mySD % is_init) then
        write(*,*) 'ScattData % init FAILED! (Scalars)'
        stop 10
      end if
      ! Allocatables, one at a time
      if (.not. allocated(mySD % E_grid)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Not Allocated)'
        stop 10
      else if (size(mySD % E_grid) /= 2) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Size)'
        stop 10
      else if ((mySD % E_grid(1) /= ONE) .or. &
        (mySD % E_grid(2) /= 20.0_8)) then
        write(*,*) 'ScattData % init FAILED! (E_Grid Incorrect Values)'
        stop 10
      end if   
      call mySD % clear()
      
      write(*,*)
      write(*,*) 'ScattData % init Passed!'
      write(*,*) '---------------------------------------------------'
      
      ! Clear all the data set up for this test.
    end subroutine test_init
    
end program test_scattdata


