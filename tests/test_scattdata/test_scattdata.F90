program test_scattdata
  use ace_header
  use constants
  use scattdata_header
  
  implicit none
  
  REAL(8), PARAMETER :: TEST_TOL = 1E-10_8

!===============================================================================
  
  write(*,*) '***************************************************'
  write(*,*) 'Testing ScattData'
  write(*,*) '***************************************************'
  
  ! Test clear routine
  call test_clear()
  
  ! Test initialization routine
  call test_init()
  
  ! Test convert_file4 routine
  call test_convert_file4()
  
  ! Test convert_file6 routine
  call test_convert_file6()
  
  
  ! Test calc_mu_bounds
  call test_calc_mu_bounds()
  
  ! Test calc_E_bounds
  call test_calc_E_bounds()
  
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
      NR = ZERO
      NEin = TWO
      allocate(Ein(NEin))
      Ein = (/ONE, TWO/)
      allocate(P(NEin))
      P = ONE
      NPEout = TWO
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
      R = (/ONE, ZERO/)
      A = (/ONE, 0.5_8/)
      ! allocate edist % data to ((3+2+NR+NEin + 2*(2+5*NPEout)))
      allocate(edist % data((3+2+NR+NEin + 2*(2+5*NPEout))))
      ! Set edist. The first Ein's R and A are twice the normal values so it is
      ! clear if we accidentally hit them.
      edist % data = (/real(NR,8), real(NEin,8), Ein, L, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, TWO * R, TWO * A, &
        real(INTTp,8), real(NPEout,8), Eout, PDF, CDF, R, A/)
      ! Set the reference solution
      distro_ref(:, 1) = (/0.1565176427_8, 0.2580539668_8, 0.4254590641_8, &
        0.7014634088_8, 1.1565176427_8/)
      distro_ref(:, 2) = (/0.5409883534_8, 0.4948293954_8, 0.4797586878_8, &
        0.4948293954_8, 0.5409883534_8/)
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
      distro_ref(:, 1) = 0.5_8
      distro_ref(:, 2) = 0.5_8
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
      if ((any(distro(:,1) /= 0.5_8)) .or. (any(distro(:,2) /= 0.5_8))) then
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
      distro_ref(:, 1) = (/ZERO, 0.2_8, 0.5_8, 0.7_8, ONE/)
      distro_ref(:, 2) = (/ZERO, 0.25_8, 0.5_8, 0.75_8, ONE/)
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
      bins_ref(MU_LO, :) = (/0, 0, 2, 2/)
      bins_ref(MU_HI, :) = (/0, 2, 2, 0/)
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
      bins_ref(MU_LO, :) = (/0, 0, 2, 2/)
      bins_ref(MU_HI, :) = (/0, 2, 2, 0/)
      interp_ref = ZERO
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
     
end program test_scattdata
