program test_scatt
  use ace_header
  use constants
  use scattdata_class
  use ndpp_scatt
  
  implicit none
  
  REAL(8), PARAMETER :: TEST_TOL = 1E-10_8

!===============================================================================
  
  write(*,*) '***************************************************'
  write(*,*) 'Testing Scattering Capabilities'
  write(*,*) '***************************************************'
  
  ! Test clear routine
  call test_clear()
  
  ! Test initialization routine
  call test_init()
  
  ! Test convert_file4 routine
  call test_convert_file4()
  
  ! Test convert_file6 routine
  call test_convert_file6()
  
  ! Test convert_distro
  call test_convert_distro()
  
  ! Test cm2lab
  call test_cm2lab()
  
  ! Test calc_mu_bounds
  call test_calc_mu_bounds()
  
  ! Test calc_E_bounds
  call test_calc_E_bounds()
  
  ! Test integrate_energyangle_file4_leg
  call test_integrate_file4_leg()
  
  ! Test integrate_energyangle_file6_leg
  call test_integrate_file6_leg()
  
  ! Test interp_distro
  call test_interp_distro()
  
  ! Perform integral test of scattering capability
  call test_calc_scatt()
  
  write(*,*)
  write(*,*) '***************************************************'
  write(*,*) 'Test of Scattering Capabilities Passed!'
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
      rxn % threshold = 1
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
      ! (n,fiss), (n,f), (n,nf), (n,2nf), (n,3nf), MT>=200
      do i = 1, 6
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
      
      write(*,*) 'Testing File 4 Distribution'
      ! Next test a file 4 distribution initialization (adist)
      call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)
      ! To reduce the scope of each if-block, we will check the scalars,
      ! then the allocatables, then the pointers
      ! First check the scalars
      if ((.not. mySD % is_init) .or. &
        mySD % scatt_type /= SCATT_TYPE_LEGENDRE .or. &
        mySD % order /= 6 .or. mySD % awr /= ONE .or. mySD % NE /= 2) then
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
      write(*,*) 'Testing File 6 Distribution'
      rxn % has_angle_dist = .false.
      rxn % adist % n_energy = 0
      deallocate(rxn % adist % energy)
      ! Now set up myedist with two Eins, each with 2 Eouts
      NR = 0
      NE = 2
      NP = 2
      myedist%law = 44
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
        mySD % order /= 6 .or. mySD % awr /= ONE .or. mySD % NE /= 2) then
        write(*,*) 'ScattData % init FAILED! (Scalars)'
        write(*,*) mySD % is_init, mySD % scatt_type, mySD % order
        write(*,*) mySD % awr, mySD % NE
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
      write(*,*) 'Testing Isotropic Distribution'
      edist => null()
      call mySD % init(nuc, rxn, edist, E_bins, scatt_type, order, mu_bins)
      ! To reduce the scope of each if-block, we will check the scalars,
      ! then the allocatables, then the pointers
      ! First check the scalars
      if ((.not. mySD % is_init) .or. &
        mySD % scatt_type /= SCATT_TYPE_LEGENDRE .or. &
        mySD % order /= 6 .or. mySD % awr /= ONE .or. mySD % NE /= 2) then
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
        associated(mySD % edist) .or. &
        (.not. associated(mySD % E_bins, E_bins))) then
        write(*,*) 'ScattData % init FAILED! (Associated)'
        stop 10
      end if
      call mySD % clear()
      
      ! Finally, repeat the isotropic case, but set a threshold energy
      ! above E_bins(1) so that the E_grid(1) == rxn % threshold
      write(*,*) 'Testing Isotropic Distribution w/ Threshold'
      call rxn % adist % clear()
      rxn % threshold = 2
      ! Clear the adist that was written in the previous init
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
   
!===============================================================================
! TEST_CONVERT_FILE4 Tests the conversion of File 4 ACE data to a tabular 
! distribution.
!===============================================================================    
    
    subroutine test_convert_file4()
      integer :: num_pts       ! Number of angular points to store
      integer :: iE            ! Energy index to act on
      real(8), allocatable :: mu(:)     ! tabular mu points
      real(8), allocatable :: distro(:) ! resultant distro (# pts)
      type(DistAngle), target  :: myadist   ! adist to be pointed to (intel issue)
      type(DistAngle), pointer :: adist ! My angle dist
      real(8), allocatable :: Eouts(:)  ! Energy out grid @ Ein
      integer :: INTT          ! Energy out INTT grid
      real(8) :: dmu           ! delta-mu for mu grid
      integer :: i             ! Generic loop counter  
      real(8), allocatable :: distro_ref(:)
      integer :: NP            ! Number of pts for ANGLE_TABULAR distribution
      
      adist => myadist
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing convert_file4'
      write(*,*)
      
      ! Allocate/Set mu and distro
      num_pts = 5
      allocate(mu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do i = 1, num_pts - 1
        mu(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      mu(num_pts) = ONE
      allocate(distro(num_pts))
      distro = ZERO
      allocate(distro_ref(num_pts))
      
      ! Set adist (one energy point, at 1 MeV, with data that corresponds to the
      ! type of data we will need to test.
      adist % n_energy = 1
      allocate(adist % energy(1))
      adist % energy(1) = ONE
      allocate(adist % type(1))
      allocate(adist % location(1))
      iE = 1
      INTT = -1
      
      ! Test an isotropic distribution
      write(*,*) 'Testing Isotropic Distribution'
      adist % type(iE) = ANGLE_ISOTROPIC
      adist % location(1) = 0
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check Eouts and INTT
      if (allocated(Eouts)) then
        if (Eouts(1) /= ZERO .or. Eouts(2) /= INFINITY) then
          write(*,*) 'convert_file4 FAILED! (Isotropic Eouts Values)'
          stop 10
        end if
      else
        write(*,*) 'convert_file4 FAILED! (Isotropic Eouts Allocation)'
        stop 10
      end if
      if (INTT /= HISTOGRAM) then
        write(*,*) 'convert_file4 FAILED! (Isotropic INTT Value)'
        stop 10
      end if
      ! Check distro value
      distro_ref = 0.5_8
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Isotropic Distro Values)'
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
        
      ! Test Equiprobable distribution. Do not need to check Eouts and INTT
      ! here since these were already tested above.
      ! Will test two cases:
      ! 1) an isotropic distribution
      ! 2) linear distribution
      ! Set quantities which are the same for both:
      adist % type(iE) = ANGLE_32_EQUI
      adist % location(iE) = 1
      allocate(adist % data(NUM_EP + 2))
      ! Move on to isotropic distribution
      write(*,*) 'Testing Equiprobable Distribution (Isotropic)'
      ! Set an isotropic distribution
      dmu = TWO * R_NUM_EP
      do i = 1, NUM_EP
        adist % data(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      adist % data(NUM_EP + 1) = ONE
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      distro_ref = 0.5_8
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Equiprob Isotropic Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      
      ! Repeat for linear distribution
      ! Set the distribution to be an equiprobable representation of a linear
      ! pdf.  This is calculated in supporting_calcs.xlsx
      write(*,*) 'Testing Equiprobable Distribution (Linear)'
      adist % data = &
        (/ZERO, -ONE, -0.6464466094_8, -0.5_8, -0.3876275643_8, &
        -0.2928932188_8, -0.209430585_8, -0.1339745962_8, -0.0645856533_8, &
        ZERO, 0.0606601718_8, 0.1180339887_8, 0.17260394_8, 0.2247448714_8, &
        0.2747548784_8, 0.3228756555_8, 0.3693063938_8, 0.4142135624_8, &
        0.4577379737_8, 0.5_8, 0.5411035007_8, 0.5811388301_8, 0.6201851746_8, &
        0.6583123952_8, 0.6955824958_8, 0.7320508076_8, 0.767766953_8, &
        0.8027756377_8, 0.8371173071_8, 0.8708286934_8, 0.9039432765_8, &
        0.9364916731_8, 0.9685019685_8, ONE/)
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      distro_ref = (/8.8388347646636875E-002_8, 0.21338834765811932_8, &
        0.48385358672217688_8, 0.73943449322968258_8, 0.99212549203273326_8/)
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Equiprob Linear Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      deallocate(adist % data)
      
      ! Test linear distribution.
      ! Will test three cases (each with histogram and lin-lin interpolation):
      ! 1) Isotropic
      ! 2) Linear, with grid width larger than mu
      ! 3) Linear, with grid width smaller than mu
      
      ! Set isotropic tabular representation. Only include the two end points
      adist % type(iE) = ANGLE_TABULAR
      adist % location(iE) = 1
      NP = 2
      allocate(adist % data(3 + 2 * NP))
      adist % data(1)   = ZERO
      adist % data(3)   = NP
      adist % data(4:5) = (/-ONE, ONE/)
      adist % data(6:7) = 0.5_8
      ! Calc histogram case
      write(*,*) 'Testing Tabular Distribution (Isotropic, Histogram)'
      adist % data(2) = HISTOGRAM
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      distro_ref = 0.5_8
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Isotropic, Histogram Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      ! Calc lin-lin case
      write(*,*) 'Testing Tabular Distribution (Isotropic, Lin-Lin)'
      adist % data(2) = LINEAR_LINEAR
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      distro_ref = 0.5_8
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Isotropic, Lin-Lin Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      
      ! Set linear distribution case with only end-points
      adist % data(6:7) = (/ZERO, ONE/)
      ! Calc histogram case
      write(*,*) 'Testing Tabular Distribution (Linear-1, Histogram)'
      adist % data(2) = HISTOGRAM
      distro_ref = (/ZERO, ZERO, ZERO, ZERO, ONE/)
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Linear-1, Histogram Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      ! Calc linear-linear case
      write(*,*) 'Testing Tabular Distribution (Linear-1, Lin-Lin)'
      adist % data(2) = LINEAR_LINEAR
      distro_ref = (/ZERO, 0.25_8, 0.5_8, 0.75_8, ONE/)
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Linear-1, Lin-Lin Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      deallocate(adist % data)
      
      ! Set linear distribution case with 11 points within
      NP = 11
      allocate(adist % data(3 + 2 * NP))
      adist % data(1)   = ZERO
      adist % data(3)   = NP
      adist % data(4:) = (/-ONE, -0.8_8, -0.6_8, -0.4_8, -0.2_8, ZERO, &
                            0.2_8, 0.4_8, 0.6_8, 0.8_8, ONE/)
      adist % data(4+NP:) = (/ZERO, 0.1_8, 0.2_8, 0.3_8, 0.4_8, 0.5_8, &
                                  0.6_8, 0.7_8, 0.8_8, 0.9_8, ONE/)
      ! Calc histogram case
      write(*,*) 'Testing Tabular Distribution (Linear-2, Histogram)'
      adist % data(2) = HISTOGRAM
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      distro_ref = (/ZERO, 0.2_8, 0.5_8, 0.7_8, ONE/)
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Linear-2, Histogram Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      ! Calc lin-lin case
      write(*,*) 'Testing Tabular Distribution (Linear-2, Lin-Lin)'
      adist % data(2) = LINEAR_LINEAR
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check distro value
      distro_ref = (/ZERO, 0.25_8, 0.5_8, 0.75_8, ONE/)
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Linear-2, Lin-Lin Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      
      ! Test invalid interpolation type. Just expect graceful failure (distro
      ! still zero)
      write(*,*) 'Testing Tabular Distribution (Linear-2, Invalid)'
      adist % data(2) = 17 ! 17 is arbitrarily chosen.
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check Eouts and INTT
      if (allocated(Eouts)) then
        if (Eouts(1) /= ZERO .or. Eouts(2) /= INFINITY) then
          write(*,*) 'convert_file4 FAILED! (Linear-2, Invalid Eouts Values)'
          stop 10
        end if
      else
        write(*,*) 'convert_file4 FAILED! (Linear-2, Invalid Eouts Allocation)'
        stop 10
      end if
      if (INTT /= HISTOGRAM) then
        write(*,*) 'convert_file4 FAILED! (Linear-2, Invalid INTT Value)'
        stop 10
      end if
      ! Check distro value
      distro_ref = ZERO
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Tabular Linear-2, Invalid Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      
      ! Final test, invalid angular distribution type. 
      ! Just expect graceful failure (distro still zero)
      write(*,*) 'Testing Invalid Distribution'
      adist % type(iE) = 17
      call convert_file4(iE, mu, adist, Eouts, INTT, distro)
      ! Check Eouts and INTT
      if (allocated(Eouts)) then
        if (Eouts(1) /= ZERO .or. Eouts(2) /= INFINITY) then
          write(*,*) 'convert_file4 FAILED! (Invalid Eouts Values)'
          stop 10
        end if
      else
        write(*,*) 'convert_file4 FAILED! (Invalid Eouts Allocation)'
        stop 10
      end if
      if (INTT /= HISTOGRAM) then
        write(*,*) 'convert_file4 FAILED! (Invalid INTT Value)'
        stop 10
      end if
      ! Check distro value
      distro_ref = ZERO
      if (any(distro /= distro_ref)) then
        write(*,*) 'convert_file4 FAILED! (Invalid Distro Values)'
        write(*,*) mu
        write(*,*) distro
        stop 10
      end if
      ! reset INTT, Eouts, and distro
      INTT = -1
      deallocate(Eouts)
      distro = ZERO
      
      write(*,*)
      write(*,*) 'convert_file4 Passed!'
      write(*,*) '---------------------------------------------------'
      
      ! Clean-up
      call adist % clear()
      deallocate(mu)
      deallocate(distro)
      deallocate(distro_ref)
      
    end subroutine test_convert_file4

!===============================================================================
! TEST_CONVERT_FILE6 Tests the conversion of File 6 ACE data to a tabular 
! distribution.
!===============================================================================    
    
    subroutine test_convert_file6()
      integer                   :: iE          ! Energy index to act on
      real(8), allocatable      :: mu(:)       ! tabular mu points
      type(DistEnergy), pointer :: edist       ! My energy dist
      type(DistEnergy), target  :: myedist     ! My energy dist
      real(8), allocatable      :: Eouts(:)    ! Energy out grid @ Ein
      integer                   :: INTT        ! Energy out INTT grid
      real(8), allocatable      :: distro(:,:) ! resultant distro (pts x NEout)
      real(8)                   :: dmu         ! delta-mu
      integer                   :: num_pts     ! number of mu pts
      real(8), allocatable      :: distro_ref(:,:) ! resultant distro (pts x NEout)
      integer                   :: i           ! loop counter
      
      ! Law 44/61 data holders
      integer :: NR, NEin, INTTp, NPEout
      real(8), allocatable :: Ein(:), P(:), L(:), Eout(:), PDF(:), CDF(:), &
                              LC(:,:), R(:), A(:), CSOUT(:,:), PDFang(:,:), &
                              CDFang(:,:)
      integer, allocatable :: NPang(:), JJ(:)
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing convert_file6'
      write(*,*)
      
      ! Set the quantities which dont change with all the varying tests
      edist => myedist
      
      ! Set the law-independent information
      NR = 0
      NEin = 2
      allocate(Ein(NEin))
      Ein = (/ONE, TWO/)
      allocate(P(NEin))
      P = ONE
      NPEout = 2
      allocate(Eout(NPEout))
      Eout = (/0.5_8, ONE/)
      allocate(PDF(NPEout))
      PDF = (/0.5_8, 0.5_8/)
      allocate(CDF(NPEout))
      CDF = (/ZERO, ONE/)
      
      ! Set the storage grid
      iE = 2
      num_pts = 5
      allocate(mu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do i = 1, num_pts - 1
        mu(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      mu(num_pts) = ONE
      allocate(distro(num_pts, NPEout))
      distro = ZERO
      allocate(distro_ref(num_pts, NPEout))
      
      ! This test will test laws 44 (KM) and 61 (tabular).  
      ! Included that is testing law 61 with isotropic, histogram, lin-lin,
      ! and an unknown interpolation type (although a fatal error will be 
      ! called, so I'll have to figure out how to deal with that.
      ! We will also need to try INTT < 10 and > 10.
      
      ! Test Kalbach-Mann (will do with INTT <10 and >10).
      ! Doing INTT < 10 first.
      write(*,*) 'Testing Law 44 (INTTp < 10)'
      ! First set up edist
      edist % law = 44
      ! Set L, INTTp, R, A
      allocate(L(NEin))
      L = (/6.0_8, 18.0_8/)
      INTTp = 1
      allocate(R(NPEout))
      R = (/ONE, ZERO/)
      allocate(A(NPEout))
      A = (/ONE, 0.5_8/)
      allocate(edist % data((3+2+NR+NEin + 2*(2+5*NPEout))))
      ! Set edist. The first Ein's R and A are twice the normal values so it is
      ! clear if we accidentally hit them.
      edist % data = (/real(NR,8), real(NEin,8), Ein, L, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, TWO * R, TWO * A, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, R, A/)
      ! Set the reference solution
      distro_ref(:, 1) = 0.5_8 * (/0.1565176427_8, 0.2580539668_8, &
        0.4254590641_8, 0.7014634088_8, 1.1565176427_8/)
      distro_ref(:, 2) = 0.5_8 * (/0.5409883534_8, 0.4948293954_8, &
        0.4797586878_8, 0.4948293954_8, 0.5409883534_8/)
      call convert_file6(iE, mu, edist, Eouts, INTT, distro)
      ! Check results
      if ((any((distro(:,1) - distro_ref(:,1)) > TEST_TOL)) .or. &
        (any((distro(:,2) - distro_ref(:,2)) > TEST_TOL))) then
        write(*,*) 'convert_file6 FAILED! (Invalid Distro Values - Law 44 INTT=1)'
        write(*,*) distro
        write(*,*) distro_ref
        stop 10
      end if
      if (any(Eouts /= Eout)) then
        write(*,*) 'convert_file6 FAILED! (Invalid Eouts Values - Law 44 INTT=1)'
        write(*,*) Eouts
        write(*,*) Eout
        stop 10
      end if
      if (INTT /= INTTp) then
        write(*,*) 'convert_file6 FAILED! (Invalid INTT Value - Law 44 INTT=1)'
        write(*,*) INTT, INTTp
        stop 10
      end if
      ! Reset values
      deallocate(Eouts)
      distro = ZERO
      INTT = -1
      
      ! Now set INTTp > 10, and make sure we get correct results
      write(*,*) 'Testing Law 44 (INTTp > 10)'
      INTTp = 12
      edist % data = (/real(NR,8), real(NEin,8), Ein, L, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, TWO * R, TWO * A, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, R, A/)
      call convert_file6(iE, mu, edist, Eouts, INTT, distro)
      ! Check results
      if ((any((distro(:,1) - distro_ref(:,1)) > TEST_TOL)) .or. &
        (any((distro(:,2) - distro_ref(:,2)) > TEST_TOL))) then
        write(*,*) 'convert_file6 FAILED! (Invalid Distro Values - Law 44 INTT=2)'
        write(*,*) distro
        write(*,*) distro_ref
        stop 10
      end if
      if (any(Eouts /= Eout)) then
        write(*,*) 'convert_file6 FAILED! (Invalid Eouts Values - Law 44 INTT=2)'
        write(*,*) Eouts
        write(*,*) Eout
        stop 10
      end if
      if (INTT /= 2) then
        write(*,*) 'convert_file6 FAILED! (Invalid INTT Value - Law 44 INTT=2)'
        write(*,*) INTT, 2
        stop 10
      end if
      ! Reset values
      deallocate(Eouts)
      distro = ZERO
      INTT = -1
      deallocate(edist % data)
      
      ! Test Law 61
      ! Set the things which dont depend on type of LDAT
      edist % law = 61
      INTTp = 1
      ! Start with isotropic
      write(*,*) 'Testing Law 61, Isotropic'
      L = (/6.0_8, 16.0_8/)
      allocate(LC(NPEout, NEin))
      LC(:, 1) = (/0, 0/) ! LC == 0 signifies isotropic
      LC(:, 2) = (/0, 0/) ! LC == 0 signifies isotropic
      allocate(edist % data(26))
      edist % data = (/real(NR,8), real(NEin,8), Ein, L, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, LC(:, 1), &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, LC(:, 2)/)
      ! Set the reference solution
      distro_ref(:, 1) = 0.25_8
      distro_ref(:, 2) = 0.25_8
      call convert_file6(iE, mu, edist, Eouts, INTT, distro)
      ! Check results
      if ((any((distro(:,1) - distro_ref(:,1)) > TEST_TOL)) .or. &
        (any((distro(:,2) - distro_ref(:,2)) > TEST_TOL))) then
        write(*,*) 'convert_file6 FAILED! (Invalid Distro Values - Law 61)'
        write(*,*) distro
        write(*,*) distro_ref
        stop 10
      end if
      if (any(Eouts /= Eout)) then
        write(*,*) 'convert_file6 FAILED! (Invalid Eouts Values - Law 61)'
        write(*,*) Eouts
        write(*,*) Eout
        stop 10
      end if
      if (INTT /= INTTp) then
        write(*,*) 'convert_file6 FAILED! (Invalid INTT Value - Law 61)'
        write(*,*) INTT, INTTp
        stop 10
      end if
      ! Reset values
      deallocate(Eouts)
      distro = ZERO
      INTT = -1
      deallocate(edist % data)
      
      ! Law 61, tabular 
      ! for iE = 1, will use isotropic hist and lin-lin distros
      ! for iE = 2, will use linear hist and lin-lin distros
      write(*,*) 'Testing Law 61, Tabular'
      allocate(JJ(NPEout))
      JJ = (/HISTOGRAM, LINEAR_LINEAR/)
      allocate(NPang(NPEout))
      NPang = (/2, 11/)
      allocate(CSOUT(NPang(2),NPEout))
      CSOUT = ZERO
      CSOUT(:, 1) = (/-ONE, ONE/)
      CSOUT(:, 2) = (/-ONE, -0.8_8, -0.6_8, -0.4_8, -0.2_8, ZERO, 0.2_8, &
                      0.4_8, 0.6_8, 0.8_8, ONE/)
      allocate(PDFang(NPang(2),NPEout))
      PDFang = ZERO
      PDFang(:, 1) = (/0.5_8, 0.5_8/)
      PDFang(:, 2) = (/ZERO, 0.1_8, 0.2_8, 0.3_8, 0.4_8, 0.5_8, &
                       0.6_8, 0.7_8, 0.8_8, 0.9_8, ONE/)
      allocate(CDFang(NPang(2),NPEout))
      CDFang = ZERO ! Values are unnecessary, just need it to take up space.
      ! Set L
      L = (/6.0_8, 32.0_8/)
      ! Set LC
      LC(:, 1) = (/16.0_8, 24.0_8/)
      LC(:, 2) = (/42.0_8, 77.0_8/)
      ! Build edist % data
      allocate(edist % data(166))
      edist % data = (/real(NR,8), real(NEin,8), Ein, L, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, LC(:, 1), &
        real(JJ(1),8), real(NPang(1),8), CSOUT(1:NPang(1),1), &
        PDFang(1:NPang(1),1), CDFang(1:NPang(1),1), &
        real(JJ(2),8), real(NPang(1),8), CSOUT(1:NPang(1),1), &
        PDFang(1:NPang(1),1), CDFang(1:NPang(1),1), &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, LC(:, 2), &
        real(JJ(1),8), real(NPang(2),8), CSOUT(:,2), PDFang(:,2), CDFang(:,2), &
        real(JJ(2),8), real(NPang(2),8), CSOUT(:,2), PDFang(:,2), CDFang(:,2)/)
      
      ! Test the isotropic values
      iE = 1
      call convert_file6(iE, mu, edist, Eouts, INTT, distro)
      ! Check results
      if ((any(distro(:,1) /= 0.25_8)) .or. (any(distro(:,2) /= 0.25_8))) then
        write(*,*) 'convert_file6 FAILED! (Invalid Distro Values - Law 61 Iso)'
        write(*,*) distro
        stop 10
      end if
      if (any(Eouts /= Eout)) then
        write(*,*) 'convert_file6 FAILED! (Invalid Eouts Values - Law 61 Iso)'
        write(*,*) Eouts
        write(*,*) Eout
        stop 10
      end if
      if (INTT /= INTTp) then
        write(*,*) 'convert_file6 FAILED! (Invalid INTT Value - Law 61 Iso)'
        write(*,*) INTT, INTTp
        stop 10
      end if
      ! Reset values
      deallocate(Eouts)
      distro = ZERO
      INTT = -1
      
      ! Test the linear values
      iE = 2
      distro_ref(:, 1) = 0.5_8 * (/ZERO, 0.2_8, 0.5_8, 0.7_8, ONE/)
      distro_ref(:, 2) = 0.5_8 * (/ZERO, 0.25_8, 0.5_8, 0.75_8, ONE/)
      call convert_file6(iE, mu, edist, Eouts, INTT, distro)
      ! Check results
      if ((any((distro(:,1) - distro_ref(:,1)) > TEST_TOL)) .or. &
        (any((distro(:,2) - distro_ref(:,2)) > TEST_TOL))) then
        write(*,*) 'convert_file6 FAILED! (Invalid Distro Values - Law 61 Lin)'
        write(*,*) distro(:,1)
        write(*,*) distro_ref(:,1)
        write(*,*) distro(:,2)
        write(*,*) distro_ref(:,2)
        stop 10
      end if
      if (any(Eouts /= Eout)) then
        write(*,*) 'convert_file6 FAILED! (Invalid Eouts Values - Law 61 Lin)'
        write(*,*) Eouts
        write(*,*) Eout
        stop 10
      end if
      if (INTT /= INTTp) then
        write(*,*) 'convert_file6 FAILED! (Invalid INTT Value - Law 61 Lin)'
        write(*,*) INTT, INTTp
        stop 10
      end if
      ! Reset values
      deallocate(Eouts)
      distro = -ONE
      INTT = -1
      
      ! Invalid Law
      write(*,*) 'Testing Invalid Law'
      edist % law = 7
      call convert_file6(iE, mu, edist, Eouts, INTT, distro)
      ! Check results
      if ((any(distro(:,1) /= -ONE)) .or. (any(distro(:,2) /= -ONE))) then
        write(*,*) 'convert_file6 FAILED! (Invalid Distro Values - Invalid Law)'
        write(*,*) distro(:,1)
        write(*,*) distro_ref(:,1)
        write(*,*) distro(:,2)
        write(*,*) distro_ref(:,2)
        stop 10
      end if
      if (allocated(Eouts)) then
        write(*,*) 'convert_file6 FAILED! (Allocated Eouts Values - Invalid Law)'
        stop 10
      end if
      if (INTT /= -1) then
        write(*,*) 'convert_file6 FAILED! (Invalid INTT Value - Invalid Law)'
        write(*,*) INTT
        stop 10
      end if      
      
      write(*,*)
      write(*,*) 'convert_file6 Passed!'
      write(*,*) '---------------------------------------------------'
      
      call edist % clear()
      nullify(edist)
      
    end subroutine test_convert_file6

!===============================================================================
! TEST_CONVERT_DISTRO Tests the functionality of scatt_convert_distro
!===============================================================================

    subroutine test_convert_distro()
      type(ScattData)           :: mySD    ! Testing object
      type(DistAngle), target   :: myadist ! adist to be pointed to (intel issue)
      type(DistAngle), pointer  :: adist   ! My angle dist
      type(DistEnergy), pointer :: edist   ! My energy dist
      type(DistEnergy), target  :: myedist ! My energy dist
      integer                   :: num_pts ! Number of angular pts in mu
      real(8)                   :: dmu     ! delta-mu
      integer                   :: NE      ! # incoming energy pts
      integer                   :: i       ! Loop counter
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing convert_distro'
      write(*,*)
      
      ! convert_distro simply takes a ScattData object and enters either 
      ! convert_file6 or convert_file4 depending on the type of distribution
      ! that is associated with it.
      ! There are four code paths for this routine, which are envoked by the 
      ! following situations: 
      ! 1) an energy-angle distribution is passed
      ! 2) an angle distribution is passed
      ! 3) both energy-angle and angle distributions are passed 
      ! 4) no distributions are passed.
      ! 3) and 4) yield in fatal errors and therefore cant be tested here.
      ! So we will test 1) and 2) by setting up edists and adists and ensuring
      ! the correct routine was called. The validity of the convert_file* 
      ! results does not need to be checked here since other unit tests fill
      ! that role.
      
      ! Begin by setting up mySD
      mySD % is_init = .true.
      NE = 2
      mySD % NE = NE
      ! Allocate/Set mu and distro
      num_pts = 5
      allocate(mySD % mu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do i = 1, num_pts - 1
        mySD % mu(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      mySD % mu(num_pts) = ONE
      allocate(mySD % distro(NE))
      allocate(mySD % distro(1) % data(num_pts, 2))
      allocate(mySD % distro(2) % data(num_pts, 2))
      mySD % distro (1) % data = ZERO
      mySD % distro (2) % data = ZERO
      ! Set other info needed in mySD
      allocate(mySD % Eouts(NE))
      allocate(mySD % INTT(NE))
      mySD % INTT = 0
      ! Set the angular distribution info 
      adist => myadist
      adist % n_energy = NE
      allocate(adist % energy(NE))
      adist % energy(1) = ONE
      adist % energy(2) = TWO
      allocate(adist % type(NE))
      allocate(adist % location(NE))
      adist % type = ANGLE_ISOTROPIC
      adist % location = 0
      allocate(adist % data(1))
      ! Set energy-angle distribution info
      edist => myedist
      edist % law = 44
      allocate(edist % data(31))
      edist % data = (/ZERO, TWO, ONE, TWO, 6.0_8, 18.0_8, ONE, TWO, 0.5_8, &    
        ONE, 0.5_8, 0.5_8, ZERO, ONE, TWO, ZERO, TWO, ONE, ONE, TWO, 0.5_8, &
        ONE, 0.5_8, 0.5_8, ZERO, ONE, ONE, ZERO, ONE, 0.5_8/)
      
      ! Test passing energy-angle distribution
      write(*,*) 'Testing File 6 Distribution'
      mySD % edist => edist
      call mySD % convert_distro()
      ! Check we got to convert_file6 by ensuring Eouts is not ZERO and INFINITY
      ! as it would be for file 4
      if (all(mySD % Eouts(1) % data(:) == (/ZERO, INFINITY/)) .and. &
        all(mySD % Eouts(2) % data(:) == (/ZERO, INFINITY/))) then
        write(*,*) 'convert_distro FAILED! (File 6)'
        stop 10
      end if
      ! Reset the info
      mySD % distro (1) % data = ZERO
      mySD % distro (2) % data = ZERO
      mySD % INTT = 0
      deallocate(mySD % Eouts(1) % data, mySD % Eouts(2) % data)
      nullify(mySD % edist)
      
      ! Test passing angle distribution
      write(*,*) 'Testing File 4 Distribution'
      mySD % adist => adist
      call mySD % convert_distro()
      ! Check we got to convert_file4 by ensuring Eouts are ZERO and INFINITY
      if (all(mySD % Eouts(1) % data(:) /= (/ZERO, INFINITY/)) .and. &
        all(mySD % Eouts(2) % data(:) /= (/ZERO, INFINITY/))) then
        write(*,*) 'convert_distro FAILED! (File 4)'
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'convert_distro Passed!'
      write(*,*) '---------------------------------------------------'
      ! clean up
      call edist % clear()
      call adist % clear()
      call mySD  % clear()
      
    end subroutine test_convert_distro

!===============================================================================
! TEST_CM2LAB Tests the conversion of a CM distro to the lab frame
!===============================================================================

    subroutine test_cm2lab()
      real(8)              :: awr             ! Atomic Weight Ratio
      real(8)              :: Q               ! Binding Energy of reaction
      real(8)              :: Ein             ! Incoming energy
      real(8), allocatable :: mu(:)           ! Angular grid
      real(8), allocatable :: distro(:,:)     ! The distribution to convert
      real(8), allocatable :: distro_in(:,:)  ! The input distro for each test
      real(8), allocatable :: distro_ref(:,:,:) ! The reference sol'n
      real(8)              :: dmu             ! delta-mu
      integer              :: num_pts         ! Number of angular pts
      integer              :: NP              ! Number of energy grid
      integer              :: i               ! loop counter
      integer              :: nz_loc          ! loc of 1st non-zero ref value
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing cm2lab'
      write(*,*)
      
      ! In this test we will take an angular distribution and convert it to CM
      ! for the cases of R > 1 and R < 1.  To properly test the presence of
      ! multiple Eout distributions, we will include an isotropic and linear
      ! distribution.
      
      ! Set the storage grid and distro values
      num_pts = 201
      NP = 2
      allocate(mu(num_pts))
      allocate(distro(num_pts, NP))
      allocate(distro_in(num_pts, NP))
      allocate(distro_ref(num_pts, NP, 2))
      dmu = TWO / real(num_pts - 1, 8)
      do i = 1, num_pts
        mu(i) = -ONE + real(i - 1, 8) * dmu
        if (i == num_pts) mu(i) = ONE
        ! distro_in(:, 1) is isotropic
        distro_in(i, 1) = 0.5_8
        ! Set distro_in(:, 1) ref soln for R > 1 (from sage workbook)
        distro_ref(i, 1, 1) = 0.00212765957446809_8*(mu(i)*mu(i)) / &
          sqrt((mu(i)*mu(i)) + 55224.0_8) + &
          0.00425531914893617_8*mu(i) + 0.00212765957446809_8 * &
          sqrt((mu(i)*mu(i)) + 55224.0_8)

        ! Set distro_in(:, 1) reference solution for R < 1(comes from sage workbook)
        if ((mu(i)*mu(i)) > 0.00166530611099991_8) then
          distro_ref(i, 1, 2) = 0.500416847233746_8 * (mu(i)*mu(i)) / &
            sqrt((mu(i)*mu(i)) - 0.00166530611099991_8) + &
            1.00083369446749_8 * mu(i) + &
            0.500416847233746_8 * sqrt((mu(i)*mu(i)) - 0.00166530611099991_8)
        else
          nz_loc = i
        end if
        
        ! distro_in(:, 2) is linear with mu
        distro_in(i, 2) = 0.5_8 * (mu(i) + ONE)
        ! Set distro(:, 2) ref soln for R > 1 (from sage workbook)
        distro_ref(i, 2, 1) = (0.00212765957446809_8*(mu(i)*mu(i)) + &
          0.00212765957446809_8*sqrt((mu(i)*mu(i)) + 55224.0_8)*mu(i) + &
          0.497872340425532_8)*(0.00425531914893617_8*(mu(i)*mu(i)) / &
          sqrt((mu(i)*mu(i)) + 55224.0_8) + &
          0.00851063829787234_8*mu(i) + 0.00425531914893617_8 * &
          sqrt((mu(i)*mu(i)) + 55224.0_8))
        ! Set distro(:, 2) reference solution (comes from sage workbook)
        if ((mu(i)*mu(i)) > 0.00166530611099991_8) then
          distro_ref(i, 2, 2) = (0.500416847233746_8*(mu(i)*mu(i)) + &
            0.500416847233746_8*sqrt((mu(i)*mu(i)) - 0.00166530611099991_8) * &
            mu(i) - 0.000416847233745687_8)*(1.00083369446749_8*(mu(i)*mu(i)) / &
            sqrt((mu(i)*mu(i)) - 0.00166530611099991_8) + &
            2.00166738893498_8*mu(i) + 1.00083369446749_8*sqrt((mu(i)*mu(i)) - &
            0.0016653061109999_8))
        end if
      end do
      
      distro_ref(1:nz_loc, 1, 2) = ZERO
      distro_ref(1:nz_loc, 2, 2) = ZERO
      
      ! Set Ein, will keep it a constant throughout this test.
      Ein = ONE
      
      ! Test R > 1
      write(*,*) 'Testing R > 1'
      awr = 235.0_8
      Q = ZERO
      distro = distro_in
      call cm2lab(awr, Q, Ein, mu, distro_in, distro)
      ! Check results
      if ((any(abs(distro(:,1) - distro_ref(:,1,1)) > 1.0E-6_8)) .or. &
        (any(abs(distro(:,2) - distro_ref(:,2,1)) > 1.0E-3_8))) then
        write(*,*) 'cm2lab FAILED! (Invalid Distro Values - R > 1)'
        write(*,*) maxval(abs(distro(:,1)-distro_ref(:,1,1)))
        write(*,*) maxval(abs(distro(:,2)-distro_ref(:,2,1)))
        stop 10
      end if
      
      ! Test R < 1
      write(*,*) 'Testing R < 1'
      awr = 0.999167_8
      Q = ZERO
      distro = distro_in
      call cm2lab(awr, Q, Ein, mu, distro_in, distro)
      ! Check results
      ! To not subject all domains of the solution to the same error bounds
      ! check the zero, first non-zero point, and remaining points separate.
      ! Check the first zeros
      if ((any(distro(:nz_loc,1) /= ZERO)) .or. &
        (any(distro(:nz_loc,2) /= ZERO))) then
        write(*,*) 'cm2lab FAILED! (Invalid Distro Values - R < 1, Zeros)'
        stop 10
      end if
      ! Check the first non-zero value, since this will be the most likely to
      ! have the highest error
      if ((abs(distro(nz_loc+1,1) - distro_ref(nz_loc+1,1,2)) > 2E-3_8) .or. &
        (abs(distro(nz_loc+1,2) - distro_ref(nz_loc+1,2,2)) > 3E-4_8)) then
        write(*,*) 'cm2lab FAILED! (Invalid Distro Values - R < 1, First Non-zero)'
        write(*,*) distro(nz_loc+1,1), distro_ref(nz_loc+1,1,2)
        write(*,*) distro(nz_loc+1,2), distro_ref(nz_loc+1,2,2)
        stop 10
      end if
      ! Finally, check the rest   
      if ((any(abs(distro(nz_loc+2:,1) - &
        distro_ref(nz_loc+2:,1,2)) > 4.0E-3_8)) .or. &
        (any(abs(distro(nz_loc+2:,2) - &
        distro_ref(nz_loc+2:,2,2)) > 3.0E-4_8))) then
        write(*,*) 'cm2lab FAILED! (Invalid Distro Values - R < 1, Last Non-Zeros)'
        write(*,*) mu(nz_loc+2:)
        write(*,*) 
        write(*,*) distro(nz_loc+2:,2)
        write(*,*) 
        write(*,*) distro_ref(nz_loc+2:,2,2)
        write(*,*) 
        write(*,*) (abs(distro(nz_loc+2:,2) - distro_ref(nz_loc+2:,2,2)))
        write(*,*) maxval(abs(distro(nz_loc+2:,2) - distro_ref(nz_loc+2:,2,2)))
        write(*,*) maxloc(abs(distro(nz_loc+2:,2) - distro_ref(nz_loc+2:,2,2)))
        stop 10
      end if
      
      write(*,*) 'Max Error in R < 1 Isotropic Distribution: '
      write(*,*) maxval(abs(distro(nz_loc+2:,1) - distro_ref(nz_loc+2:,1,2)))
      write(*,*) 'Max Error in R < 1 Linear Distribution: '
      write(*,*) maxval(abs(distro(nz_loc+2:,2) - distro_ref(nz_loc+2:,2,2)))
    
      write(*,*)
      write(*,*) 'cm2lab Passed!'
      write(*,*) '---------------------------------------------------'
      
    end subroutine test_cm2lab
    
!===============================================================================
! TEST_CALC_MU_BOUNDS Tests the calculation of the angular boundary 
!===============================================================================    
    
    subroutine test_calc_mu_bounds()
      real(8), allocatable :: mu(:)       ! Angular bin points
      real(8), allocatable :: E_bins(:)   ! Outgoing Energy Bins
      real(8), allocatable :: interp(:,:) ! Outgoing mu interpolants
      real(8), allocatable :: vals(:,:)   ! Outgoing mu values
      integer, allocatable :: bins(:,:)   ! Outgoing mu indices
      real(8)              :: Q           ! Reaction Q-Value
      real(8)              :: awr         ! Atomic-weight ratio
      real(8)              :: Ein         ! Incoming energy
      real(8)              :: dmu         ! change in mu
      integer              :: num_pts     ! Number of angular pts to use
      integer              :: i           ! temp loop counter
      integer              :: num_Eout    ! Number of Eout bins
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing calc_mu_bounds'
      write(*,*)
      
      ! Allocate/Set mu
      num_pts = 21
      allocate(mu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do i = 1, num_pts - 1
        mu(i) = -ONE + real(i - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      mu(num_pts) = ONE
      
      ! Allocate E_bins, interp, vals, bins
      num_Eout = 3
      allocate(E_bins(num_Eout))
      allocate(interp(2, num_Eout - 1))
      allocate(  vals(2, num_Eout - 1))
      allocate(  bins(2, num_Eout - 1))
      interp = ZERO
      vals   = ZERO
      bins   = 0
      
      ! Set the awr and Q-value to that for hydrogen (R<1 is the toughest case)
      awr = 0.999167_8 ! H-1 value
      Q   = ZERO       
      
      ! The possible situations one could encounter for the energy vars are:
      ! 1) Ein > all of E_bins
      ! 2) Ein < all of E_bins
      ! 3) Ein within one bin
      ! The result of calc_mu_bounds should always return a value b/t [-1,1].  
      
      ! 1) Ein > all of E_bins
      write(*,*) 'Testing Ein > all E_bins'
      Ein = 20.0_8
      E_bins = (/1E-11_8, ONE, TWO/)
      call calc_mu_bounds(awr, Q, Ein, E_bins, mu, interp, vals, bins)
      ! For transfer to the lower group, hand calcs show that:
      ! mu_low = -1 and mu_high = 0.22537631014397342822
      if ((any((interp(:,1) - &
        (/ZERO, 0.2537631014397342822_8/)) > TEST_TOL)) .or. &
        (any((vals(:,1) - &
        (/-ONE, 0.22537631014397342822_8/)) > TEST_TOL)) .or. &
        (any(bins(:,1) /= (/1, 13/)))) then
        write(*,*) 'calc_mu_bounds FAILED! (Ein > all of E_bins, MU_LO)'
        stop 10
      end if
      ! For transfer to the higher group, hand calcs show that:
      ! mu_low = 0.22537631014397342822 and mu_high = 0.31741314579775205019
      if ((any((interp(:,2) - (/0.2537631014397342822_8, &
        0.1741314579775205019_8/)) > TEST_TOL)) .or. &
        (any((vals(:,2) - (/0.22537631014397342822_8, &
        0.31741314579775205019_8/)) > TEST_TOL)) .or. &
        (any(bins(:,2) /= (/13, 14/)))) then
        write(*,*) 'calc_mu_bounds FAILED! (Ein > all of E_bins, MU_HI)'
        stop 10
      end if
      ! Reset values
      interp = ZERO
      vals   = ZERO
      bins   = 0
      
      ! 2) Ein < all of E_bins
      write(*,*) 'Testing Ein < all E_bins'
      Ein = 1E-11_8
      E_bins = (/ONE, TWO, 20.0_8/)
      call calc_mu_bounds(awr, Q, Ein, E_bins, mu, interp, vals, bins)
      ! Check results
      ! Since Ein < min(E_bins), we know that all mu transfers are to mu=1
      ! (and thus improbable). This means interp == 1, bins == size(mu) - 1,
      ! and vals == 1
      if ((any(interp /= ONE)) .or. (any(vals /= ONE)) .or. &
        (any(bins /= num_pts-1))) then
        write(*,*) 'calc_mu_bounds FAILED! (Ein < min(E_bins))'
        stop 10
      end if
      ! Reset values
      interp = ZERO
      vals   = ZERO
      bins   = 0
      
      ! 3) Ein within one bin
      write(*,*) 'Testing Ein within one E_bin'
      Ein = 1.5_8
      E_bins = (/ONE, TWO, 20.0_8/)
      call calc_mu_bounds(awr, Q, Ein, E_bins, mu, interp, vals, bins)
      ! Check results
      ! Since Ein < E_bins(2), we know that in the 1st group, the upper mu
      ! transfer is to mu=1 (like in the previous case).
      ! Next, a hand calculation shows that mu=0.81666661634070423168 for Eout=1
      if ((any((interp(:,1) - &
        (/0.1666661634070423168_8, ONE/)) > TEST_TOL)) .or. &
        (any((vals(:,1) - (/0.81666661634070423168_8, ONE/)) > TEST_TOL)) .or. &
        (any(bins(:,1) /= (/19, num_pts - 1/)))) then
        write(*,*) 'calc_mu_bounds FAILED! (Ein within one bin, MU_LO)'
        stop 10
      end if
      ! Since Ein < E_bins(2:3), we know that the 2nd group mu transfers are to 
      ! mu=1 (and thus improbable). This means interp == 1, bins == num_pts-1,
      ! and vals == 1
      if ((any(interp(:,2) /= ONE)) .or. (any(vals(:,2) /= ONE)) .or. &
        (any(bins(:,2) /= num_pts-1))) then
        write(*,*) 'calc_mu_bounds FAILED! (Ein within one bin, MU_HI)'
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'calc_mu_bounds Passed!'
      write(*,*) '---------------------------------------------------'
      
    end subroutine test_calc_mu_bounds
    
!===============================================================================
! TEST_CALC_E_BOUNDS Tests the calculation of the energy integration points 
!===============================================================================  

    subroutine test_calc_E_bounds()
      real(8), allocatable :: E_bins(:)   ! Energy group boundaries
      real(8), allocatable :: Eout(:)     ! Output energies
      integer              :: INTT        ! Output energies interpolation type
      real(8), allocatable :: interp(:,:) ! Outgoing E interpolants
                                          ! corresponding to the energy groups
      integer, allocatable :: bins(:,:)   ! Outgoing E indices corresponding
                                          ! to the energy groups
      integer              :: NEout       ! Number of Eout pts
      integer              :: NEbins      ! Number of E_bins
      ! Reference solution holders
      real(8), allocatable :: interp_ref(:,:) 
      integer, allocatable :: bins_ref(:,:)  
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing calc_E_bounds'
      write(*,*)
      
      ! This test will examine E_bins situated at various positions relative to
      ! Eout (above, below, between, etc).  
      ! We will have to do this with linear and histogram interpolations.
      
      ! Allocate and set up Eout
      NEout  = 3
      allocate(Eout(NEout))
      Eout = (/ONE, TWO, 3.0_8/)
      
      ! Allocate and set up E_bins
      NEbins = 5
      allocate(E_bins(NEbins))
      E_bins = (/0.25_8, 0.75_8, TWO, 2.5_8, 4.0_8/)
      
      ! Allocate interp, vals, bins, and their reference solution spaces
      allocate(interp    (2, NEbins - 1))
      allocate(bins      (2, NEbins - 1))
      allocate(interp_ref(2, NEbins - 1))
      allocate(bins_ref  (2, NEbins - 1))
      
      ! First test linear interpolation
      write(*,*) 'Testing Linear Interpolation'
      INTT = LINEAR_LINEAR
      ! Set reference solution
      bins_ref(MU_LO, :) = (/-1, -1, 2, 2/)
      bins_ref(MU_HI, :) = (/-1, 2, 2, -2/)
      interp_ref(MU_LO, :) = (/ZERO, ZERO, ZERO, 0.5_8/)
      interp_ref(MU_HI, :) = (/ZERO, ZERO, 0.5_8, ZERO/)
      call calc_E_bounds(E_bins, Eout, INTT, interp, bins)
      ! Check results
      if (any(bins /= bins_ref)) then
        write(*,*) 'calc_E_bounds FAILED! (Incorrect Bins - INTT = Linear)'
        stop 10
      end if
      if (any(interp /= interp_ref)) then
        write(*,*) 'calc_E_bounds FAILED! (Incorrect Interps - INTT = Linear)'
        stop 10
      end if
      
      ! Test Histogram interpolation.
      write(*,*) 'Testing Histogram Interpolation'
      INTT = HISTOGRAM
      ! Set reference solution
      bins_ref(MU_LO, :) = (/-1, -1, 2, 2/)
      bins_ref(MU_HI, :) = (/-1, 2, 2, -2/)
      interp_ref = -interp_ref
      call calc_E_bounds(E_bins, Eout, INTT, interp, bins)
      ! Check results
      if (any(bins /= bins_ref)) then
        write(*,*) 'calc_E_bounds FAILED! (Incorrect Bins - INTT = Histogram)'
        stop 10
      end if
      if (any(interp /= interp_ref)) then
        write(*,*) 'calc_E_bounds FAILED! (Incorrect Interps - INTT = Histogram)'
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'calc_E_bounds Passed!'
      write(*,*) '---------------------------------------------------'

    end subroutine test_calc_E_bounds
     
!===============================================================================
! TEST_INTEGRATE_FILE4_LEG Tests the ability to integrate a set of file 4 
! distributions to produce Legendre moments.
!===============================================================================

    subroutine test_integrate_file4_leg()
      real(8), allocatable :: fEmu(:)           ! Energy-angle distro to act on
      real(8), allocatable :: mu(:)             ! fEmu angular grid
      real(8), allocatable :: interp(:,:)       ! interpolants of mu values
      real(8), allocatable :: vals(:,:)         ! mu values
      integer, allocatable :: bins(:,:)         ! indices of fEmu corresponding to
                                                ! the group boundaries
      integer              :: order             ! Number of moments to find
      real(8), allocatable :: distro(:,:)       ! Resultant integrated distribution
      
      real(8), allocatable :: distro_ref(:,:)   ! Reference solution
      integer              :: imu               ! mu_val indices
      integer              :: num_pts, num_G    ! Number of mu pts and energy grps
      real(8)              :: dmu               ! mu spacing
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing integrate_energyangle_file4_leg'
      write(*,*)
      
      ! This test will do find the legendre moments of two cases:
      ! 1) Linear anisotropic distribution with only one energy group and 
      !    no truncation
      ! 2) Linear anisotropic distribution with 3 energy groups across the
      !    angular domain.
      
      ! For both cases, we will use a mu with 5 grid pts. Set mu and fEmu now
      ! Allocate/Set mu
      num_pts = 5
      allocate(mu(num_pts))
      allocate(fEmu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do imu = 1, num_pts - 1
        mu(imu) = -ONE + real(imu - 1, 8) * dmu
        fEmu(imu) = 0.5_8 * (mu(imu) + ONE)
      end do
      ! Set the end point to exactly ONE
      mu(num_pts) = ONE
      fEmu(num_pts) = ONE
      ! Set the number of orders
      order = 6
      
      ! Start with case 1
      write(*,*) 'Testing One Energy Group'
      ! Set variables and allocate as needed
      num_G = 1
      allocate(distro(order, num_G))
      distro = ZERO
      allocate(interp(2, num_G))
      interp(:,1) = (/ZERO, ONE/)
      allocate(vals(2, num_G))
      vals(:,1) = (/-ONE, ONE/)
      allocate(bins(2, num_G))
      bins(:, 1) = (/1, num_pts-1/)
      ! Set reference solution (analytical linear aniso moments)
      allocate(distro_ref(order, num_G))
      distro_ref(:,1) = (/ONE, ONE / 3.0_8, ZERO, ZERO, ZERO, ZERO/)
      ! Run the function
      call integrate_energyangle_file4_leg(fEmu, mu, interp, vals, bins, &
        order - 1, distro)
      ! Test the results
      if (any(abs(distro - distro_ref) > TEST_TOL)) then
        write(*,*) 'integrate_energy_angle_file4_leg FAILED! (Case 1)'
        stop 10
      end if
      ! Clear results
      deallocate(distro, distro_ref, interp, vals, bins)
      
      ! Case 2, three groups
      write(*,*) 'Testing Three Energy Groups'
      ! Set variables and allocate as needed
      num_G = 3
      allocate(distro(order, num_G))
      distro = ZERO
      allocate(interp(2, num_G))
      interp(:, 1) = (/ZERO, 0.5_8/)
      interp(:, 2) = (/0.5_8, 0.5_8/)
      interp(:, 3) = (/0.5_8, ONE/)
      allocate(vals(2, num_G))
      vals(:, 1) = (/-ONE, -0.75_8/)
      vals(:, 2) = (/-0.75_8, 0.25_8/)
      vals(:, 3) = (/0.25_8, ONE/)
      allocate(bins(2, num_G))
      bins(:, 1) = (/1, 1/) ! 
      bins(:, 2) = (/1, 3/)
      bins(:, 3) = (/3, num_pts-1/)
      ! Set reference solution (analytical linear aniso moments w/ 3 groups
      ! as defined above)
      ! This reference solution is found in the Sage notebook:
      ! integrate_file4_leg_reference.sws
      allocate(distro_ref(order, num_G))
      distro_ref(:, 1) = (/0.015625_8, -0.0130208333333_8, 0.008544921875_8, &
        -0.00341796875_8, -0.0010503133138_8, 0.00387191772461_8/)
      distro_ref(:, 2) = (/0.375_8, -0.0520833333333_8, -0.13671875_8, &
        0.0400390625_8, 0.0531209309896_8, -0.00587463378906_8/)
      distro_ref(:, 3) = (/0.609375_8, 0.3984375_8, 0.128173828125_8, &
        -0.03662109375_8, -0.0520706176758_8, 0.00200271606445_8/)
      ! Run the function
      call integrate_energyangle_file4_leg(fEmu, mu, interp, vals, bins, &
        order - 1, distro)
      ! Test the results
      if (any(abs(distro - distro_ref) > TEST_TOL)) then
        write(*,*) 'integrate_energy_angle_file4_leg FAILED! (Case 2)'
        write(*,*) abs(distro(:,1)-distro_ref(:,1))
        write(*,*)
        write(*,*) abs(distro(:,2)-distro_ref(:,2))
        write(*,*)
        write(*,*) abs(distro(:,3)-distro_ref(:,3))
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'integrate_energyangle_file4_leg Passed!'
      write(*,*) '---------------------------------------------------'
      
    end subroutine test_integrate_file4_leg

!===============================================================================
! TEST_INTEGRATE_FILE6_LEG Tests the ability to integrate a set of file 6 
! energy/angle distributions to produce Legendre moments.
!===============================================================================

    subroutine test_integrate_file6_leg()
      real(8), allocatable :: fEmu(:,:)       ! Energy-angle distro to act on
      real(8), allocatable :: mu(:)           ! fEmu angular grid
      real(8), allocatable :: Eout(:)         ! Outgoing energies
      real(8), allocatable :: E_bins(:)       ! Energy grp boundaries
      real(8), allocatable :: interp(:,:)     ! interpolants of E values
      integer, allocatable :: bins(:,:)       ! indices of fEmu corresponding to
                                              ! the group boundaries
      integer              :: order           ! Number of moments to find
      real(8), allocatable :: distro(:,:)     ! Resultant integrated distribution
      
      real(8), allocatable :: distro_ref(:,:) ! Reference solution
      integer              :: imu             ! mu_val indices
      integer              :: num_pts, num_G  ! Number of mu pts and energy grps
      integer              :: num_Eout, iE    ! Number of outgoing E pts, and index
      real(8)              :: dmu             ! mu spacing
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing integrate_energyangle_file6_leg'
      write(*,*)
      
      ! In this test, we will set up angular distributions which span the entire
      ! angular space, and set an array of these up in energy. We then
      ! will have energy groups to integrate the distribution over.
      ! These groups will be such that we test all the portions of the 
      ! file6 integration routine. This yields the following cases:
      ! 1) only one Eout distribution
      ! 2) all groups within two adjacent energy out points
      ! 3) groups spanning multiple Eout points
      
      ! Set values constant for all cases
      num_pts = 5
      allocate(mu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do imu = 1, num_pts - 1
        mu(imu) = -ONE + real(imu - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      mu(num_pts) = ONE
      ! Set the number of orders
      order = 6
      
      ! Case 1 - only one Eout
      write(*,*) 'Testing One Eout'
      ! Set energy group structure
      num_G = 1
      num_Eout = 1
      allocate(E_bins(num_G+1))
      E_bins = (/ONE, TWO/)
      ! Set angular distribution to linearly anisotropic
      allocate(fEmu(num_pts, num_Eout))
      do iE = 1, num_Eout
        do imu = 1, num_pts
          fEmu(imu, iE) = 0.5_8 * (mu(imu) + ONE)
        end do
      end do
      ! Set Eout
      allocate(Eout(num_Eout))
      Eout = 1.5_8
      ! Set the bins and interp vals
      allocate(bins(2, num_G))
      bins(:, 1) = (/1, 1/)
      allocate(interp(2, num_G))
      interp(:, 1) = (/0.5_8, ONE/)
      ! Allocate and ready distro
      allocate(distro(order, num_G))
      distro = ZERO
      ! Set the reference solution
      allocate(distro_ref(order, num_G))
      distro_ref(:,1) = (/ONE, ONE / 3.0_8, ZERO, ZERO, ZERO, ZERO/)
      ! Run the calcs!
      call integrate_energyangle_file6_leg(fEmu, mu, Eout, E_bins, interp, &
        bins, order - 1, distro)
      ! Test the results
      if (any(abs(distro - distro_ref) > TEST_TOL)) then
        write(*,*) 'integrate_energy_angle_file6_leg FAILED! (Case 1)'
        write(*,*) abs(distro(:,1)-distro_ref(:,1))
        stop 10
      end if
      ! Deallocs to get ready for case 2
      deallocate(E_bins, fEmu, Eout, bins, interp, distro, distro_ref)
      
      ! Case 2, all groups within one set of Eouts
      write(*,*) 'Testing Multiple Groups within two Eout points'
      ! Set energy group structure
      num_G = 2
      num_Eout = 2
      allocate(E_bins(num_G+1))
      E_bins = (/ONE, TWO, 3.0_8/)
      ! Set angular distribution to isotropic at the lower Eout, 
      ! linearly anisotropic at the higher Eout
      allocate(fEmu(num_pts, num_Eout))
      fEmu(:, 1) = 0.5_8
      do imu = 1, num_pts
        fEmu(imu, 2) = 0.5_8 * (mu(imu) + ONE)
      end do
      ! Set Eout
      allocate(Eout(num_Eout))
      Eout = (/ONE, 3.0_8/)
      ! Set the bins and interp vals
      allocate(bins(2, num_G))
      bins(:, 1) = (/1, 1/)
      bins(:, 2) = (/1, 1/)
      allocate(interp(2, num_G))
      interp(:, 1) = (/ ZERO, 0.5_8/)
      interp(:, 2) = (/0.5_8,   ONE/)
      ! Allocate and ready distro
      allocate(distro(order, num_G))
      distro = ZERO
      ! Set the reference solution
      allocate(distro_ref(order, num_G))
      ! These are calculated in the Sage worksheet associated with this test.
      distro_ref(:, 1) = (/ONE, ONE / 12.0_8, ZERO, ZERO, ZERO, ZERO/)
      distro_ref(:, 2) = (/ONE, ONE / 4.0_8, ZERO, ZERO, ZERO, ZERO/)
      ! Run the calcs!
      call integrate_energyangle_file6_leg(fEmu, mu, Eout, E_bins, interp, &
        bins, order - 1, distro)
      ! Test the results
      if (any(abs(distro - distro_ref) > TEST_TOL)) then
        write(*,*) 'integrate_energy_angle_file6_leg FAILED! (Case 2)'
        write(*,*) abs(distro(:,1)-distro_ref(:,1))
        write(*,*) abs(distro(:,2)-distro_ref(:,2))
        stop 10
      end if
      ! Deallocs to get ready for case 3
      deallocate(E_bins, fEmu, Eout, bins, interp, distro, distro_ref)
      
      ! Case 3, groups spanning multiple Eout points
      ! Will do with 1 group that surrounds 3 Eouts
      write(*,*) 'Testing Groups spanning multiple Eout points'
      ! Set energy structures
      num_G = 1
      allocate(E_bins(num_G+1))
      ! Set Eout
      num_Eout = 3
      allocate(Eout(num_Eout))
      Eout = (/ONE, TWO, 3.0_8/)
      E_bins = (/ZERO, 4.0_8/)
      ! Set angular distribution to isotropic at the lower Eout, 
      ! increasing linearly anisotropic at the middle Eout
      ! decreasing linearly anisotropic at the higher Eout
      allocate(fEmu(num_pts, num_Eout))
      fEmu(:, 1) = 0.5_8
      do imu = 1, num_pts
        fEmu(imu, 2) = 0.5_8 * (mu(imu) + ONE)
        fEmu(imu, 3) = -0.5_8 * mu(imu) + 0.5_8
      end do
      ! Set the bins and interp vals
      allocate(bins(2, num_G))
      bins(:, 1) = (/-1, -(num_Eout-1)/)
      allocate(interp(2, num_G))
      interp(:, 1) = (/ZERO, ZERO/)
      ! Allocate and ready distro
      allocate(distro(order, num_G))
      distro = ZERO
      ! Set the reference solution
      allocate(distro_ref(order, num_G))
      ! These are calculated in the Sage worksheet associated with this test.
      distro_ref(:, 1) = (/ONE, ONE / 12.0_8, ZERO, ZERO, ZERO, ZERO/)
      ! Run the calcs!
      call integrate_energyangle_file6_leg(fEmu, mu, Eout, E_bins, interp, &
        bins, order - 1, distro)
      ! Test the results
      if (any(abs(distro - distro_ref) > TEST_TOL)) then
        write(*,*) 'integrate_energy_angle_file6_leg FAILED! (Case 3)'
        write(*,*) abs(distro(:,1)-distro_ref(:,1))
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'integrate_energyangle_file6_leg Passed!'
      write(*,*) '---------------------------------------------------'
      
    end subroutine test_integrate_file6_leg

!===============================================================================
! TEST_INTERP_DISTRO Tests the functionality of scatt_interp_distro.
!===============================================================================

    subroutine test_interp_distro()
      type(ScattData)           :: mySD    ! Testing object
      type(Nuclide), target     :: mynuc   ! Testing nuclide
      type(Nuclide), pointer    :: nuc     ! Testing nuclide
      type(Reaction), allocatable, target :: rxn(:) ! Elastic and inelastic rxn
      type(DistEnergy), target  :: myedist ! My energy dist
      real(8), allocatable      :: mu_out(:)
      real(8), allocatable      :: distro_out(:,:)
      real(8), allocatable      :: distro_ref(:,:)
      integer                   :: num_pts ! Number of angular pts in mu
      real(8)                   :: dmu     ! delta-mu
      real(8)                   :: Ein
      integer                   :: imu     ! Loop counter
      
      ! interp_distro takes an array of distributions in this % distro,
      ! and an incoming energy, stores the result in the result(distro).
      ! The major parts of my testing will be, with legendre results requested,
      ! to pass in elastic and inelastic reactions, an incoming energy below
      ! the threshold energy, above the minimum energy in the 
      ! nuc % energy grid and one w/in the energy range.  
      ! We then have to test something with nested distributions 
      ! (and thus edist % p_valid exists),
      ! CM and Lab, and a file 4 and file 6 distribution.
      ! Our distributions for this will be an isotropic and linear distro.
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing interp_distro'
      write(*,*)
      
      ! Set up the reactions I need (elastic and inelastic)
      allocate(rxn(2))
      ! The first: elastic
      rxn(1) % MT = ELASTIC
      rxn(1) % Q_value = ZERO
      rxn(1) % multiplicity = 1
      rxn(1) % has_angle_dist = .true.
      rxn(1) % has_energy_dist = .false.
      rxn(1) % scatter_in_cm = .true.
      rxn(1) % threshold = 1
      ! Second, inelastic
      rxn(2) % MT = N_LEVEL
      rxn(2) % Q_value = ZERO
      rxn(2) % multiplicity = 2
      rxn(2) % has_angle_dist = .false.
      rxn(2) % has_energy_dist = .true.
      rxn(2) % scatter_in_cm = .true.
      rxn(2) % threshold = 2
      allocate(rxn(2) % sigma(2))
      rxn(2) % sigma = (/0.5_8, ONE/)
      rxn(2) % edist => myedist
      myedist % law = 44
      myedist % p_valid % n_pairs = 2
      allocate(myedist % p_valid % x(2))
      myedist % p_valid % x = (/ONE, TWO/) 
      allocate(myedist % p_valid % y(2))
      myedist % p_valid % y = (/0.5_8, ONE/)
      
      ! Set nuc
      nuc => mynuc
      nuc % awr = TWO
      nuc % n_grid = 2
      allocate(nuc % energy(nuc % n_grid))
      allocate(nuc % elastic(nuc % n_grid))
      nuc % energy =  (/ONE, TWO/)
      nuc % elastic = (/0.5_8, ONE/)
      
      ! Build up mySD
      mySD % scatt_type = SCATT_TYPE_LEGENDRE
      mySD % order = 5
      mySD % groups = 1
      mySD % awr = nuc % awr
      num_pts = 5001
!~       num_pts = 11
      allocate(mySD % mu(num_pts))
      dmu = TWO / real(num_pts - 1, 8)
      do imu = 1, num_pts - 1
        mySD % mu(imu) = -ONE + real(imu - 1, 8) * dmu
      end do
      ! Set the end point to exactly ONE
      mySD % mu(num_pts) = ONE
      allocate(mySD % E_bins(mySD % groups + 1))
      mySD % E_bins = (/1E-11_8, 20.0_8/)
      
      ! Now lets do an elastic distribution
      mySD % rxn => rxn(1)
      mySD % adist => rxn(1) % adist
      mySD % NE = 2
      allocate(mySD % E_grid(mySD % NE))
      mySD % E_grid = (/ONE, TWO/)
      allocate(mySD % distro(mySD % NE))
      allocate(mySD % distro(1) % data(num_pts, mySD % groups))
      mySD % distro(1) % data(:,1) = 0.5_8
      allocate(mySD % distro(2) % data(num_pts, mySD % groups))
      ! Set as linear
      do imu = 1, num_pts
        mySD % distro(2) % data(imu,1) = 0.5_8 * (mySD % mu(imu) + ONE)
      end do
      Ein = 1.5_8
      ! Calculate the resultant distro
      allocate(distro_out(mySD % order + 1, mySD % groups))
      distro_out = mySD % interp_distro(mu_out, nuc, Ein)
      ! Set the reference solution
      ! The reference is simply the linear distribution converted to lab,
      ! multiplied by sigS (which is 0.75 due to interpolation). 
      ! These results come from the sage notebook for this test.
      allocate(distro_ref(mySD % order + 1, mySD % groups))
      distro_ref(:,1) = (/0.75_8, 0.4625_8, 0.177758602769303_8, &
        0.0428571428571428_8, 0.00554988817179086_8, ZERO/)
!~       if (any(abs(distro_out - distro_ref) > 1.0E-7_8)) then
!~         write(*,*) 'interp_distro FAILED! (Elastic)'
!~         write(*,*) distro_out
!~         write(*,*) distro_ref
!~         write(*,*) maxval(abs(distro_out(:,1)-distro_ref(:,1)))
!~         write(*,*) maxloc(abs(distro_out(:,1)-distro_ref(:,1)))
!~         stop 10
!~       end if
      
      ! Lets do the inelastic distribution, with Ein < rxn%threshold
      mySD % rxn => rxn(2)
      nullify(mySD % adist)
      mySD % edist => myedist
      distro_out = ZERO
      distro_out = mySD % interp_distro(mu_out, nuc, Ein)
      ! Set the reference solution
      distro_ref = ZERO
      if (any(abs(distro_out - distro_ref) > TEST_TOL)) then
        write(*,*) 'interp_distro FAILED! (Inelastic, < Threshold)'
        write(*,*) distro_out
        write(*,*) distro_ref
        write(*,*) maxval(abs(distro_out(:,1)-distro_ref(:,1)))
        write(*,*) maxloc(abs(distro_out(:,1)-distro_ref(:,1)))
        stop 10
      end if
      
      ! Inelastic, but now with Ein > max incoming energy
      rxn % threshold = 1 ! Set it back to the beginning so we dont need to care
      Ein = 3.0_8
      allocate(mySD % Eouts(mySD % NE))
      allocate(mySD % Eouts(1) % data(1))
      mySD % Eouts(1) % data(1) = ONE
      allocate(mySD % Eouts(2) % data(1))
      mySD % Eouts(2) % data(1) = TWO
      allocate(mySD % INTT(mySD % NE))
      mySD % INTT = LINEAR_LINEAR
      distro_out = ZERO
      distro_out = mySD % interp_distro(mu_out, nuc, Ein)
      ! Set the reference solution (2/.75 is to use the same reference
      ! as the first (elastic) case, but modifying it for the multiplicative
      ! constants out front.
      distro_ref(:,1) = TWO / 0.75_8 * (/0.75_8, 0.4625_8, 0.177758602769303_8, &
        0.0428571428571428_8, 0.00554988817179086_8, ZERO/)
      if (any(abs(distro_out - distro_ref) > 1.0E-7_8)) then
        write(*,*) 'interp_distro FAILED! (Inelastic, >  Max Ein)'
        write(*,*) distro_out
        write(*,*) distro_ref
        write(*,*) maxval(abs(distro_out(:,1)-distro_ref(:,1)))
        write(*,*) maxloc(abs(distro_out(:,1)-distro_ref(:,1)))
        stop 10
      end if
      
      ! Inelastic, but now with Ein within the range
      Ein = 1.1_8
      ! What the heck, lets turn back on CM2Lab
      rxn % scatter_in_cm = .true.
      distro_out = ZERO
      distro_out = mySD % interp_distro(mu_out, nuc, Ein)
      ! Set the reference solution
      distro_ref(:,1) = (/0.605_8, 0.201666666666667_8, 0.0314322417310492_8, &
        ZERO, -0.000696052807637037_8, ZERO/)
      if (any(abs(distro_out - distro_ref) > 1.0E-7_8)) then
        write(*,*) 'interp_distro FAILED! (Inelastic, w/in Ein range)'
        write(*,*) distro_out
        write(*,*) distro_ref
        write(*,*) maxval(abs(distro_out(:,1)-distro_ref(:,1)))
        write(*,*) maxloc(abs(distro_out(:,1)-distro_ref(:,1)))
        stop 10
      end if
      
      ! Test tabular response type
      distro_out = ZERO
      mySD % scatt_type = SCATT_TYPE_TABULAR
      distro_out = mySD % interp_distro(mu_out, nuc, Ein)
      ! Set the reference solution
      distro_ref = ZERO
      if (any(abs(distro_out - distro_ref(1:5,:)) > TEST_TOL)) then
        write(*,*) 'interp_distro FAILED! (Tabular)'
        write(*,*) distro_out
        write(*,*) distro_ref
        write(*,*) maxval(abs(distro_out(:,1)-distro_ref(:,1)))
        write(*,*) maxloc(abs(distro_out(:,1)-distro_ref(:,1)))
        stop 10
      end if
      
      write(*,*)
      write(*,*) 'interp_distro Test Passed!'
      write(*,*) '---------------------------------------------------'
      
    end subroutine test_interp_distro

!===============================================================================
! TEST_SCATT Performs an integral test of the scattdata module and 
! ndpp_scatt%calc_scatt()
!===============================================================================

    subroutine test_calc_scatt()
      type(Nuclide), target     :: mynuc            ! Nuclide to work with 
      type(Nuclide), pointer    :: nuc => null()    ! Nuclide passed to calc_scatt
      type(Reaction), pointer   :: rxn => null()    ! Temporary pointer to use
      type(DistAngle), pointer  :: adist => null()  ! temporary ptr
      type(DistEnergy), pointer :: edist => null()  ! temporary ptr
      real(8), allocatable      :: energy_bins(:)   ! Energy group boundaries
      integer                   :: scatt_type       ! Type of data to obtain
      integer                   :: order            ! Order of data to obtain
      integer                   :: mu_bins          ! # of angular bins to use
      real(8), allocatable      :: results(:,:,:)   ! Output of calc_scatt
      real(8), allocatable      :: reference(:,:,:) ! Reference solution
      real(8)                   :: thin_tol         ! Thinning tolerance (NYI)
      integer                   :: NE_Ein, NEsig    ! # of Eins and sigS
      type(DistEnergy), allocatable, target :: myedist(:)   ! Edists to point to
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Testing calc_scatt'
      write(*,*)
      
      ! In this test we will set up a nuclide with four reactions:
      ! 1) no distribution (ISO), 2) an adist, and 
      ! 3) an edist with 2 nested edists (first being valid, second not being
      ! a valid law). 4) one that is not scattering.
      ! We will use 2 groups, 2 Eins per reaction, and a nuc%energy array with 
      ! 3 energy entries. We will do this once for legendre and tabular
      ! (But for now, only legendre since tabular is NYI)
      
      ! This test is intended as an integral test of calc_scatt and scattdata to
      ! increase confidence at NDPPs ability to handle the many complicated 
      ! nuclides in ENDF data.
      
      ! Set some problem constants
      NE_Ein = 2
      NEsig = 3
      
      ! Lets set up our nuclide (Only set values we will need)
      nuc => mynuc
      nuc % awr = 100.0_8 ! The larger it is the fewer mu pts we need to
                          ! adequately match the reference
      nuc % n_grid = NEsig
      allocate(nuc % energy(NEsig))
      nuc % energy = (/ONE, TWO, 3.0_8/)
      allocate(nuc % elastic(NEsig))
      nuc % elastic = (/0.25_8, 0.5_8, ONE/)
      nuc % n_reaction = 4
      allocate(nuc % reactions(nuc % n_reaction))
      
      ! Set up reactions
      ! 1) No distribution (isotropic)
      rxn => nuc % reactions(1)
      rxn % MT = ELASTIC
      rxn % Q_value = ZERO
      rxn % multiplicity = 1
      rxn % threshold = 1
      rxn % scatter_in_cm = .true.
      allocate(rxn % sigma(NEsig))
      rxn % has_angle_dist  = .false.
      rxn % has_energy_dist = .false.
      
      ! 2) an adist
      rxn => nuc % reactions(2)
      rxn % MT = N_LEVEL
      rxn % Q_value = ZERO
      rxn % multiplicity = 1
      rxn % threshold = 1
      rxn % scatter_in_cm = .true.
      allocate(rxn % sigma(NEsig))
      rxn % sigma = (/ONE, 0.5_8, 0.25_8/)
      rxn % has_angle_dist  = .true.
      rxn % has_energy_dist = .false.
      adist => rxn % adist
      adist % n_energy = NE_Ein
      allocate(adist % energy(NE_Ein))
      adist % energy = (/1.5_8, 2.5_8/)
      allocate(adist % type(NE_Ein))
      adist % type = ANGLE_TABULAR
      allocate(adist % location(NE_Ein))
      adist % location(1) = 1
      adist % location(2) = 7
      allocate(adist % data(14))
      adist % data(1:7)  = (/ZERO, real(LINEAR_LINEAR,8), TWO, -ONE, ONE, &
        0.5_8, 0.5_8/)
      adist % data(8:13) = (/real(LINEAR_LINEAR,8), TWO, -ONE, ONE, &
        0.5_8, 0.5_8/)
      nullify(adist)
        
      ! 3) an edist with 2 nested edists (first being valid, second not being
      ! a valid law)
      rxn => nuc % reactions(3)
      rxn % MT = N_LEVEL
      rxn % Q_value = -0.02_8
      rxn % multiplicity = 2
      rxn % threshold = 2
      rxn % scatter_in_cm = .false.
      allocate(rxn % sigma(NEsig - rxn % threshold + 1))
      rxn % sigma = (/ONE, TWO/)
      rxn % has_angle_dist  = .false.
      rxn % has_energy_dist = .true.
      allocate(myedist(2))
      rxn % edist => myedist(1)
      edist => rxn % edist
      edist % law = 44
      edist % p_valid % n_pairs = 2
      allocate(edist % p_valid % x(2))
      edist % p_valid % x = (/ONE, TWO/) 
      allocate(edist % p_valid % y(2))
      edist % p_valid % y = (/0.5_8, ONE/)
      allocate(edist % data(31))
      edist % data = (/ZERO, TWO, TWO, 3.0_8, 6.0_8, 18.0_8, TWO, TWO, ONE, &    
        TWO, 0.5_8, 0.5_8, ZERO, ONE, 0.5_8, ZERO, 0.5_8, 0.5_8, TWO, TWO, &
        TWO, 3.0_8, 0.5_8, 0.5_8, ZERO, ONE, ONE, ZERO, ONE, 0.5_8/)
      edist % next => myedist(2)
      edist % next % law = 66 ! Invalid so we shouldn't need any other info
      edist % next % next => null()
      nullify(edist)
      
      ! 4) Not scattering
      rxn => nuc % reactions(4)
      rxn % MT = N_FISSION
      rxn % has_energy_dist = .false.
      
      ! Set problem parameters for Legendre case
      allocate(energy_bins(3))
      energy_bins = (/ONE, TWO, 3.0_8/)
      scatt_type = SCATT_TYPE_LEGENDRE
      order = 5
      mu_bins = 3001
      thin_tol = 1.0E-8_8
      
      ! Create the reference solution
      allocate(reference(order + 1, size(energy_bins) - 1, NEsig))
      reference(:, :, 1) = ZERO ! Ein == lower bound, therefore no distribution
      reference(:, 1, 2) = (/3.0_8, 0.08864337353812592746_8, &
        0.03257903545734936728_8, 0.00057911919634860452_8, &
        0.00012837354419521826_8, 1.39180690734675999921E-6_8/)
      reference(:, 2, 2) = ZERO ! No upscattering (Ein==lower bound of grp 2)
      reference(:, 1, 3) = ZERO ! No scattering from Ein=3 to group 1
      reference(:, 2, 3) = (/5.25_8, 0.63440390433224552687_8, &
        0.1543723225381074322_8, 0.01712913596389337176_8, &
        0.00201270147011421999_8, 0.00017011895592263393_8/)
      
      ! Lets run the test
      write(*,*) 'Testing Legendre Case'
      call calc_scatt(nuc, energy_bins, scatt_type, order, results, mu_bins, &
        thin_tol)
      write(*,*) 'results, grp 1, Ein1', results(:,1,1)
      write(*,*) 'results, grp 2, Ein1', results(:,2,1)
      write(*,*) 'results, grp 1, Ein2', results(:,1,2)
      write(*,*) 'results, grp 2, Ein2', results(:,2,2)
      write(*,*) 'results, grp 1, Ein3', results(:,1,3)
      write(*,*) 'results, grp 2, Ein3', results(:,2,3)
      write(*,*) maxval(results(:, 1, 2) - reference(:, 1, 2))
      write(*,*) maxval(results(:, 2, 3) - reference(:, 2, 3))
      write(*,*)
      write(*,*) 'calc_scatt Test Passed!'
      write(*,*) '---------------------------------------------------'
      
    end subroutine test_calc_scatt

end program test_scatt
