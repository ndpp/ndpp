module legendre
  
  use ace_header
  use constants
  use jacobian!,   only: calc_jacobian_legendre
  use error,      only: fatal_error
  use global,     only: message
  
  implicit none

  contains

!===============================================================================
! CM2LAB_LEGENDRE_NUMERICAL Converts the center-of-mass distribution in 
! ScattData to a laboratory frame of reference distribution.  This is performed 
! with an numerical expansion then reintegration.
! This works for the Legendre output option, a similar function exists in 
! ndpp_tabular for tabular output option.
!===============================================================================

    pure subroutine cm2lab_legendre_numeric(awr, Q, E_grid, order_L, order_CM, &
      moments)
      real(8), intent(in) :: awr ! Atomic Weight Ratio for this nuclide
      real(8), intent(in) :: Q   ! Binding Energy of reaction, for finding R
      real(8), intent(in) :: E_grid(:) ! Incoming energy grid
      integer, intent(in) :: order_L   ! Legendre order of the lab moments (user input)
      integer, intent(in) :: order_CM  ! Legendre order of the CM moments (set at 20)
      real(8), intent(inout) :: moments(:,:,:) ! The moments to convert
      
      real(8) :: R, mumin ! Reduced Mass and minimum mu value for integration
      integer :: groups, NE
      integer :: g, iE, l, imu
      real(8) :: mu(CM2LAB_NUMERIC_PTS)
      real(8) :: mu_cm(CM2LAB_NUMERIC_PTS)
      real(8) :: dmu, tempsqrt
      real(8) :: distro(CM2LAB_NUMERIC_PTS)
      
      groups = size(moments, dim = 2)
      NE = size(moments, dim = 3)
      
      do iE = 1, NE
        ! Calculate the effective mass ratio
        if (Q /= ZERO) then
          R = awr * sqrt(ONE - (awr + ONE) * Q / (awr * E_grid(iE)))
        else 
          R = awr
        end if
        
        ! Calculate the lower bounds of integration (-1 for everything besides
        ! H-1)
        if (R < ONE) then 
          mumin = sqrt(ONE - R * R)
        else
          mumin = -ONE
        end if
        
        ! Calculate the lab (mu) and CM (mu_cm) \mu grid points.
        dmu = (ONE - mumin) / CM2LAB_NUMERIC_PTS
        do imu = 1, CM2LAB_NUMERIC_PTS
          mu(imu) = mumin + real(imu - 1, 8) * dmu
          mu_cm(imu) = (mu(imu) * mu(imu) - ONE + mu(imu) * sqrt(R * R + &
            mu(imu) * mu(imu) - ONE)) / R
        end do
        
        ! Expand the Center-of-mass distribution on to the mu_cm grid
        do g = 1, groups
          do imu = 1, CM2LAB_NUMERIC_PTS
            distro(imu) = ZERO
            do l = 1, order_CM
              distro(imu) = distro(imu) + (real(l - 1, 8) + 0.5_8) * &
                moments(l, g, iE) * calc_pn(l - 1, mu_cm(imu))
            end do
          end do
          
          ! Convert the CM distro to the laboratory system
          do imu = 1, CM2LAB_NUMERIC_PTS
            tempsqrt = sqrt(mu(imu) * mu(imu) + R * R - ONE)
            distro(imu) = distro(imu) * (TWO * mu(imu) + tempsqrt + &
              mu(imu) * mu(imu) / tempsqrt) / R
          end do
          
          ! Finally, re-integrate this tabular distribution and store the results
          ! as the laboratory moments
          moments(:, g, iE) = ZERO
          do imu = 1, CM2LAB_NUMERIC_PTS - 1
            moments(:, g, iE) = moments(:, g, iE) + &
              calc_int_pn_tablelin(order_L, mu(imu), mu(imu + 1), &
              distro(imu), distro(imu + 1))
          end do
        end do
        
      end do
      
    end subroutine cm2lab_legendre_numeric

!===============================================================================
! CM2LAB_LEGENDRE_JACOB Converts the center-of-mass distribution in ScattData to a 
! laboratory frame of reference distribution.  This is performed with an
! analytical expression making use of a Jacobian.
! This works for the Legendre output option, a similar function exists in 
! ndpp_tabular for tabular output option.
!===============================================================================

    pure subroutine cm2lab_legendre_jacob(awr, Q, E_grid, order_L, order_CM, moments)
      real(8), intent(in) :: awr ! Atomic Weight Ratio for this nuclide
      real(8), intent(in) :: Q   ! Binding Energy of reaction, for finding R
      real(8), intent(in) :: E_grid(:) ! Incoming energy grid
      integer, intent(in) :: order_L   ! Legendre order of the lab moments (user input)
      integer, intent(in) :: order_CM  ! Legendre order of the CM moments (set at 20)
      real(8), intent(inout) :: moments(:,:,:) ! The moments to convert
      
      real(8), allocatable :: lab_moments(:)
      real(8), allocatable :: jacob(:,:)
      real(8) :: R ! Reduced Mass
      integer :: groups, NE
      integer :: g, iE, l, lcm
      
      ! This problem is solved by setting up the following matrix calculation:
      ! Lab_Moments(vector) = Jacobian(Matrix) * CM_Moments(vector)
      ! The Jacobian was derived via the SAGE computer algebra system.
      
      groups = size(moments, dim = 2)
      NE = size(moments, dim = 3)
      
      allocate(lab_moments(order_L))
      allocate(jacob(order_L, order_CM))
      
      do iE = 1, NE
        ! From equation 234 in Methods for Processing ENDF/B-VII (pg 2798)
        R = sqrt(awr * awr  * (ONE + Q * (awr + ONE) / (awr * E_grid(iE))))
        call calc_jacobian_legendre(R, order_L, order_CM, jacob)
        do g = 1, groups
          do l = 1, order_L + 1
            do lcm = 1, order_CM + 1
              lab_moments(l) = lab_moments(l) + jacob(lcm, l) * &
                moments(lcm, g, iE)
            end do
          end do
          ! Now put the new lab moments back in to the original array
          moments(1:order_L,g,iE) = lab_moments
        end do
      end do
      
      deallocate(lab_moments)
      deallocate(jacob)
      
    end subroutine cm2lab_legendre_jacob

!===============================================================================
! INTEGRATE_FILE4_LEGENDRE performs the overall control for calculating the
! Legendre moments of a file 4 (angle-only) distribution.
!===============================================================================

    subroutine integrate_file4_legendre(order, mu_bins, iE, adist, distro)
      integer, intent(in)  :: order         ! Order of the expansion
      real(8), intent(in)  :: mu_bins(:,:)  ! mu_min and mu_max for each group
      integer, intent(in)  :: iE            ! Energy index to act on
      type(DistAngle), pointer, intent(in)   :: adist ! My angle dist
      real(8), intent(inout) :: distro(:,:) ! resultant distro (order x groups)
      
      integer :: lc           ! Location inside adist % data
      integer :: istart, ilc  ! adist % data indices
      integer :: g            ! Outgoing group index
      integer :: groups       ! Number of outgoing groups
      integer :: interp, NP   ! Tabular format data (interp type and # pts)
      real(8) :: f, interpval ! Interpolation parameter and interpolated value.
      
      groups = size(distro, dim = 2)

      ! Check what type of distribution we have
      lc   = adist % location(iE)
      
      select case(adist % type(iE))
        case (ANGLE_ISOTROPIC)
          do g = 1, groups
            distro(:, g) = calc_leg_iso(order, mu_bins(MU_LO, g), mu_bins(MU_HI, g))
          end do
        case (ANGLE_32_EQUI)
          istart = lc
          do g = 1, groups
            if (mu_bins(MU_LO, g) /= mu_bins(MU_HI, g)) then
              do ilc = istart, istart + 32
                if ((mu_bins(MU_LO, g) >= adist % data(ilc)) .and. &
                  (mu_bins(MU_LO, g) <= adist % data(ilc + 1))) then
                  ! Intermediate point
                  distro(:, g) = distro(:, g) + calc_int_pn_equi(order, &
                    mu_bins(MU_LO, g), adist%data(ilc + 1), ONE/32.0_8)
                else if ((adist % data(ilc) >= mu_bins(MU_LO, g)) .and. &
                  (adist % data(ilc) <= mu_bins(MU_HI, g))) then
                  ! Entire adist bin is within mu_bin
                  distro(:, g) = distro(:, g) + &
                    calc_int_pn_equi(order, adist%data(ilc), &
                    adist%data(ilc + 1), ONE/32.0_8)
                else if ((adist % data(ilc + 1) >= mu_bins(MU_HI, g))) then
                  ! right side of mu_bin is greater than an adist bin.
                  distro(:, g) = distro(:, g) + &
                    calc_int_pn_equi(order, adist%data(ilc), mu_bins(MU_HI, g), &
                    ONE/32.0_8)
                  istart = ilc
                  exit
                end if
              end do
            end if
          end do
        case (ANGLE_TABULAR)
          istart = lc + 3
          interp = int(adist % data(lc + 1))
          NP = int(adist % data(lc + 2))
          if (interp == HISTOGRAM) then
            do g = 1, groups
              if (mu_bins(MU_LO, g) /= mu_bins(MU_HI, g)) then
                do ilc = istart, istart + NP - 2
                  if ((mu_bins(MU_LO, g) >= adist % data(ilc)) .and. &
                    (mu_bins(MU_LO, g) <= adist % data(ilc + 1))) then
                    ! Intermediate point
                    distro(:, g) = distro(:, g) + &
                      calc_int_pn_tablehist(order, mu_bins(MU_LO, g), &
                      adist%data(ilc + 1), adist%data(ilc + NP))
                  else if ((adist % data(ilc) >= mu_bins(MU_LO, g)) .and. &
                    (adist % data(ilc) <= mu_bins(MU_HI, g))) then
                    ! Entire adist bin is within mu_bin
                    distro(:, g) = distro(:, g) + &
                      calc_int_pn_tablehist(order, adist%data(ilc), &
                      adist%data(ilc + 1), adist%data(ilc + NP))
                  else if ((adist % data(ilc + 1) >= mu_bins(MU_HI, g))) then
                    ! right side of mu_bin is greater than an adist bin.
                    distro(:, g) = distro(:, g) + &
                      calc_int_pn_tablehist(order, adist%data(ilc), &
                      mu_bins(MU_HI, g), adist%data(ilc + NP))
                    istart = ilc
                    exit
                  end if
                end do
              end if
            end do
          else if (interp == LINEAR_LINEAR) then
            do g = 1, groups
              if (mu_bins(MU_LO, g) /= mu_bins(MU_HI, g)) then
                do ilc = istart, istart + NP - 2
                  if ((mu_bins(MU_LO, g) >= adist % data(ilc)) .and. &
                    (mu_bins(MU_LO, g) <= adist % data(ilc + 1))) then
                    ! Intermediate point
                    ! Calculate interpolated value
                    f = (mu_bins(MU_LO, g) - adist % data(ilc)) / &
                      (adist % data(ilc + 1) - adist % data(ilc))
                    interpval = (ONE - f) * adist % data(ilc + NP) + f * &
                      adist % data(ilc + NP + 1)
                    distro(:, g) = distro(:, g) + &
                      calc_int_pn_tablelin(order, mu_bins(MU_LO, g), &
                      adist%data(ilc + 1), interpval, adist%data(ilc + 1 + NP))
                  else if ((adist % data(ilc) >= mu_bins(MU_LO, g)) .and. &
                    (adist % data(ilc) <= mu_bins(MU_HI, g))) then
                    ! Entire adist bin is within mu_bin
                    distro(:, g) = distro(:, g) +  calc_int_pn_tablelin(order, &
                      adist%data(ilc), adist%data(ilc + 1), adist%data(ilc + NP), &
                      adist%data(ilc + 1 + NP))
                  else if ((adist % data(ilc + 1) >= mu_bins(MU_HI, g))) then
                    ! right side of mu_bin is greater than start of this adist bin.
                    ! Calculate interpolated value
                    f = (mu_bins(MU_HI, g) - adist % data(ilc)) / &
                      (adist % data(ilc + 1) - adist % data(ilc))
                    interpval = (ONE - f) * adist % data(ilc + NP) + f * &
                      adist % data(ilc + NP + 1)
                    distro(:, g) = distro(:, g) + calc_int_pn_tablelin(order, &
                      adist%data(ilc), mu_bins(MU_HI, g), adist%data(ilc + NP), &
                      interpval)
                    istart = ilc
                    exit
                  end if
                end do
              end if
            end do
          end if
      end select
    end subroutine integrate_file4_legendre

!===============================================================================
! INTEGRATE_FILE6_LEGENDRE performs the overall control for calculating the
! Legendre moments of a file 6 (combined energy-angle) distribution.
!===============================================================================

    subroutine integrate_file6_legendre(order, iE, E_bins, edist, distro)
      integer, intent(in)  :: order      ! Order of the expansion
      integer, intent(in)  :: iE         ! Energy index to act on
      real(8), intent(in)  :: E_bins(:)  ! Energy group boundaries
      type(DistEnergy), pointer, intent(in)   :: edist ! My energy-angle dist
      real(8), intent(inout) :: distro(:,:) ! resultant distro (order x groups)
      
      integer :: iEout, iEout_min, iEout_max  ! Indices of outgoing energy
      integer :: g            ! Outgoing group index
      integer :: groups       ! Number of outgoing groups
      integer :: NP           ! Number of Eout points
      real(8) :: flo, fhi     ! Interpolation parameter for the high and low pts
      real(8), allocatable :: Eout(:) ! Eout values from edist % data at Ein
      real(8), allocatable :: PDF(:)  ! f(Eout) values from edist % data at Ein (PDF)
      real(8), allocatable :: R(:)    ! Kalbach-Mann R parameter
      real(8), allocatable :: A(:)    ! Kalbach-Mann A parameter
      integer, allocatable :: adata_lc(:) ! Location of Law 61 angular distribution
      integer, allocatable :: angtype(:)  ! Type of Law 61 angular distribution
      integer :: angtype_lo, angtype_hi ! Same, but for the lo and hi points
      real(8) :: Rlo, Alo, Rhi, Ahi   ! Interpolated KM A and R
      real(8) :: PDFlo, PDFhi ! Interpolated PDF of Eout values
      integer :: INTT         ! Interpolation type between Eout data
      integer :: lc, NR, ND   ! edist%data location, num of interp regions, and
                              ! number of discrete lines
      integer :: lc1, lc2, NP1, NP2
      
      ! This routine will step through every Eout point in the edist (@ E_in),
      ! calculate the moments, and integrate them.  This will effectively follow
      ! this equation:
      ! \Sigma_{s,n,MT,l,g}(E_{in}) \= 
      !   \int_{E_g}^{E_{g-1}}\int_{-1}^{1}f(E,E',\mu) P_l{\mu} d\mu dE'
      
      ! The inner angular integral can be performed using either:
      ! calc_int_pn_tablehist, calc_int_pn_tablelin, or calc_leg_km, depending
      ! on the law and angular interpolation type.
      
      ! The calculation will proceed by:
      ! 1) Locating the Ein e-mu distro in edist
      ! 2) For each group:
      !  2.a) Find the Eout locatins in edistro which correspond to the
      !       boundaries of this group.
      !  2.b) Interpolate to get distribution parameters on the bin boundaries.
      !  2.c) Trapezoidal integration of the points within each group multiplied
      !       by the PDF of that energy point
      
      groups = size(distro, dim = 2)
      
      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      
      ! find energy bin and calculate interpolation factor -- if the energy is
      ! outside the range of the tabulated energies, choose the first or last
      ! bins
      lc = 2 + 2*NR
      
      ! determine location of outgoing energy data
      lc = int(edist % data(2 + 2*NR + int(edist % data(2 + 2*NR)) + iE))

      ! determine type of interpolation and number of discrete lines
      INTT = int(edist % data(lc + 1))
      NP   = int(edist % data(lc + 2))
      if (INTT > 10) then
        INTT = mod(INTT,10)
        ND = (INTT - INTT)/10
      else
        ND = 0
      end if
      if (ND > 0) then
        ! discrete lines present
        message = "Discrete lines in continuous tabular distributed not &
             &yet supported"
        call fatal_error()
      end if
      
      ! Now that I have NP, allocate and set Eout and the PDF for Eout
      allocate(Eout(NP))
      allocate(PDF(NP))
      lc = lc + 3
      Eout = edist % data(lc : lc + NP - 1)
      PDF  = edist % data(lc + NP : lc + 2 * NP - 1)
      ! allocate our A, R if Kalbach-Mann.
      if (edist % law == 44) then ! Kalbach-Mann
        allocate(R(NP))
        allocate(A(NP))
        R = edist % data(lc + 3 * NP : lc + 4 * NP - 1)
        A = edist % data(lc + 4 * NP : lc + 5 * NP - 1)
      else if (edist % law == 61) then ! Correlated energy-angle in tabular form
        allocate(adata_lc(NP))
        allocate(angtype(NP))
        adata_lc = int(abs(edist % data(lc + 3 * NP : lc + 4 * NP - 1)))
        do g = 1, NP
          if (edist % data(lc + 3 * NP - 1 + g) == 0) then
            angtype(g) = ANGLE_TABULAR
          else
            angtype(g) = ANGLE_ISOTROPIC
          end if
        end do
      end if
      
      ! Loop through energy groups
      iEout_min = lc
      do g = 1, groups
        ! Find the location in edist%data corresponding to the upper
        ! bound of this energy group.
        ! The lower bound is already known: iEout_min
        do iEout_max = iEout_min, NP - 1
          if (Eout(iEout_max + 1) > E_bins(g + 1)) exit
        end do
        if (iE == NP - 1) iEout_max = iEout_max - 1
        ! Find the interpolation factor for the low and high Eout points
        if (INTT == LINEAR_LINEAR) then
          flo = (E_bins(g) - Eout(iEout_min)) / &
            (Eout(iEout_min + 1) - Eout(iEout_min))
          fhi = (E_bins(g + 1) - Eout(iEout_max)) / &
            (Eout(iEout_max + 1) - Eout(iEout_max))
        else if (INTT == HISTOGRAM) then
          flo = ZERO
          fhi = ZERO
        end if
        ! Perform the integration over each Eout.
        select case(edist % law)
          case (44) ! Kalbach-Mann
            ! Find the end point R and A values
            Rlo = R(iEout_min) + flo * (R(iEout_min + 1) - R(iEout_min))
            Alo = A(iEout_min) + flo * (A(iEout_min + 1) - A(iEout_min))
            PDFlo = PDF(iEout_min) + flo * (PDF(iEout_min + 1) - PDF(iEout_min))
            Rhi = R(iEout_max) + fhi * (R(iEout_max + 1) - R(iEout_max))
            Ahi = A(iEout_max) + fhi * (A(iEout_max + 1) - A(iEout_max))
            PDFlo = PDF(iEout_max) + flo * (PDF(iEout_max + 1) - PDF(iEout_max))
            ! Trapezoidal integration: 0.5*(b-a)*(f(a)+f(b))
            ! Add the integral for the end points to our distribution
            distro(:, g) = (Eout(iEout_min + 1) - Eout(iEout_min)) * &
              (calc_leg_km(order, Rlo, Alo) * PDFlo + &
              calc_leg_km(order, R(iEout_min + 1), A(iEout_min + 1)) * &
              PDF(iEout_min + 1))
            distro(:, g) = distro(:, g) + &
              (Eout(iEout_max + 1) - Eout(iEout_max)) * &
              (calc_leg_km(order, Rhi, Ahi) * PDFhi + &
              calc_leg_km(order, R(iEout_max), A(iEout_max)) * PDF(iEout_max))
            ! Do all the intermediate points in a similar fashion
            do iEout = iEout_min + 1, iEout_max - 1
              distro(:, g) = distro(:, g) + (Eout(iEout + 1) - Eout(iEout)) * &
                (calc_leg_km(order, R(iEout), A(iEout)) * PDF(iEout) + &
                calc_leg_km(order, R(iEout + 1), A(iEout + 1)) * PDF(iEout + 1))
            end do
            distro(:, g) = 0.5_8 * distro(:, g)
          case (61) ! Correlated Energy-Angle
            ! Find the end point mu_bins and their PDF
            ! There is no interpolation, per se, in energy.  For consistency with
            ! OpenMC we will take the angular distribution closest to our group
            ! boundary Eout values.
            
            ! Perform trapezoidal integration of these Eout points within this group.
            ! Find the distribution to use for flo, and apply it
            if (flo < 0.5_8) then ! use the lower energy point's angular data
              lc1 = adata_lc(iEout_min) + 1
              NP1 = int(edist % data(lc1 + 1))
              lc2 = adata_lc(iEout_min + 1) + 1
              NP2 = int(edist % data(lc2 + 1))
            else ! Use the higher energy points angular data
              lc1 = adata_lc(iEout_min + 1) + 1
              NP1 = int(edist % data(lc1 + 1))
              lc2 = adata_lc(iEout_min + 1) + 1
              NP2 = int(edist % data(lc2 + 1))
            end if
            distro(:, g) = (Eout(iEout_min + 1) - Eout(iEout_min)) * &
              (integrate_law61_legendre(order, angtype_lo, &
              edist % data(lc1 : lc1 + 3 + 2 * NP1)) * PDFlo + &
              integrate_law61_legendre(order, angtype(iEout_min + 1), &
              edist % data(lc2 : lc2 + 3 + 2 * NP2)) * PDF(iEout_min + 1))
            ! Do the upper point in a similar fashion
            if (fhi < 0.5_8) then ! use the lower energy point's angular data
              lc1 = adata_lc(iEout_max) + 1
              NP1 = int(edist % data(lc1 + 1))
              lc2 = adata_lc(iEout_max + 1) + 1
              NP2 = int(edist % data(lc2 + 1))
            else ! Use the higher energy points angular data
              lc1 = adata_lc(iEout_max + 1) + 1
              NP1 = int(edist % data(lc1 + 1))
              lc2 = adata_lc(iEout_max + 1) + 1
              NP2 = int(edist % data(lc2 + 1))
            end if
            distro(:, g) = distro(:, g) + &
              (Eout(iEout_max + 1) - Eout(iEout_max)) * &
              (integrate_law61_legendre(order, angtype_hi, &
              edist % data(lc1 : lc1 + 3 + 2 * NP1)) * PDFhi + &
              integrate_law61_legendre(order, angtype(iEout_max + 1), &
              edist % data(lc2 : lc2 + 3 + 2 * NP2)) * PDF(iEout_max + 1))
            
            ! Do all the intermediate points in a similar fashion, 
            ! but there is no need to 'interpolate' now.
            do iEout = iEout_min + 1, iEout_max - 1
              lc1 = adata_lc(iEout) + 1
              NP1 = int(edist % data(lc1 + 1))
              lc2 = adata_lc(iEout + 1) + 1
              NP2 = int(edist % data(lc2 + 1))
              distro(:, g) = distro(:, g) + (Eout(iEout + 1) - Eout(iEout)) * &
                (integrate_law61_legendre(order, angtype(iEout), &
                edist % data(lc1 : lc1 + 3 + 2 * NP1)) * PDF(iEout) + &
                integrate_law61_legendre(order, angtype(iEout + 1), &
                edist % data(lc2 : lc2 + 3 + 2 * NP2)) * PDF(iEout + 1))
            end do
            distro(:, g) = 0.5_8 * distro(:, g)
        end select
        
        iEout_min = iEout
      end do
      
    end subroutine integrate_file6_legendre

!===============================================================================
! INTEGRATE_LAW61_LEGENDRE Integrates a single Eout data set of law 61 data.  
! This function very closely resembles integrate_file4_legendre because
! the data to operate on is exactly the same, it just exists in a different
! location for law 61.
!===============================================================================

    pure function integrate_law61_legendre(order, ang_type, adata) &
      result(distro)
      
      integer, intent(in)  :: order         ! Order of the expansion
      integer, intent(in)  :: ang_type      ! Angular distribution type
      real(8), intent(in)  :: adata(:)      ! The angular distribution at Eout
      
      real(8), allocatable :: distro(:) ! resultant distro (order x groups)
      
      integer :: lc           ! Location inside adist % adata
      integer :: istart, ilc  ! adist % adata indices
      integer :: interp, NP   ! Tabular format adata (interp type and # pts)
      
      allocate(distro(order))
      lc = 1
      ! Check what type of distribution we have
      select case(ang_type)
        case (ANGLE_ISOTROPIC)
          !!! This doesnt make sense. How can we have an isotropic distro
          ! across all angle ranges??
          distro(:) = calc_leg_iso(order,-ONE,ONE)
        case (ANGLE_TABULAR)
          istart = lc + 3
          interp = int(adata(lc + 1))
          NP = int(adata(lc + 2))
          if (interp == HISTOGRAM) then
            do ilc = istart, istart + NP - 2
              distro(:) = distro(:) + &
                calc_int_pn_tablehist(order, adata(ilc), &
                adata(ilc + 1), adata(ilc + NP))
            end do

          else if (interp == LINEAR_LINEAR) then
            do ilc = istart, istart + NP - 2
              distro(:) = distro(:) +  calc_int_pn_tablelin(order, &
                adata(ilc), adata(ilc + 1), adata(ilc + NP), &
                adata(ilc + 1 + NP))
            end do
          end if
      end select
    end function integrate_law61_legendre

!===============================================================================
! The CALC_INT_PN_* functions calculate the analytic integral of one portion of
! the group (i.e., from one set of points in the ACE table).
! The CALC_LEG_* functions calculate the legendre moments of an entire group.
!===============================================================================

!===============================================================================
! CALC_LEG_ISO calculates the Legendre expansion of an isotropic
! distribution. This function returns all n orders. It operates on an entire
! energy group. Since this function is called repeatedly,
! neither n or x is checked to see if they are in the applicable range. 
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range, 
! the return value will be 1.0_8.
!===============================================================================

    pure function calc_leg_iso(n, xlow, xhigh) result(integrals)
      integer, intent(in)  :: n     ! Legendre orders requested (inclusive)
      real(8), intent(in)  :: xlow  ! Low boundary of range of integration
      real(8), intent(in)  :: xhigh ! High boundary of range of integration
      
      real(8) :: integrals(n + 1) ! The Legendres evaluated at x (0:nth order)
      
      real(8)              :: flow, fhigh  ! the height of the function, 0.5 for iso
      
      flow = 0.5_8
      fhigh = flow
      
      ! Calculate the moments using the line integral method, but with
      ! flow == fhigh, which forces it to be flat.
      integrals =  calc_int_pn_tablelin(n, xlow, xhigh, flow, fhigh)
      
    end function calc_leg_iso
    
!===============================================================================
! CALC_LEG_KM calculates the integral of the Legendre polynomial times a 
! Kalbach correlation (LAW 44).  This function returns 
! all n orders. It is intended to operate on one outgoing energy group.
! Since this function is called repeatedly,
! neither n or x is checked to see if they are in the applicable range. 
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range, 
! the return value will be ONE
!===============================================================================
    
    pure function calc_leg_km(n, R, A) result (integrals)
      integer, intent(in)  :: n     ! Legendre orders requested (inclusive)
      
      real(8), intent(in)  :: R      ! Kalbach-Mann parameter R
      real(8), intent(in)  :: A      ! Kalbach-Mann parameter A
      
      real(8) :: integrals(n  + 1) ! The Legendres evaluated at x (0:nth order)
      
      real(8) :: cothA, RcothA, values
      
      integer :: l     ! Index of the Legendre order
      
      integrals = ZERO
      
      cothA = 1.0_8 / tanh(A)
      RcothA = R * cothA
      ! These were derived from the following Mathematica commands:
      ! f := A / (2*Sinh[A]) * (Cosh[A*x] + R * Sinh[A*x])
      ! FortranForm[HornerForm[Integrate[f*LegendreP[l,x],{x,-1,1}]]]
      ! where l is the order of interest.
      do l = 0, n
        select case (l)
          case (0)
            values = ONE
          case(1)
            values = (-R + A*RcothA)/A
          case(2)
            values = (3.0_8 + A*(A - 3.0_8*cothA))/A**2
          case(3)
            values = (-15.0_8*R + A*(15.0_8*RcothA + A*(-6.0_8*R + A*RcothA)))/A**3
          case(4)
            values = (105.0_8 + A*(A*(45 + A*(A - 10.0_8*cothA)) - 105.0_8*cothA))/A**4
          case(5)
            values = (-945.0_8*R + A*(945.0_8*RcothA + A*(-420.0_8*R + A*(105.0_8*RcothA + &
              A*(-15.0_8*R + A*RcothA)))))/A**5
          case(6)
            values = (10395.0_8 + A*(A*(4725.0_8 + A*(A*(210.0_8 + A*(A - 21.0_8*cothA)) - &
            1260.0_8*cothA)) - 10395.0_8*cothA))/A**6
          case(7)
            values = (-135135.0_8*R + A*(135135.0_8*RcothA + A*(-62370.0_8*R + A*(17325.0_8*RcothA + &
              A*(-3150.0_8*R + A*(378.0_8*RcothA + A*(-28.0_8*R + A*RcothA)))))))/A**7
          case(8)
            values = (2027025.0_8 + A*(A*(945945.0_8 + A*(A*(51975.0_8 + A*(A*(630.0_8 + A* &
              (A - 36.0_8*cothA)) - 6930.0_8*cothA)) - 270270.0_8*cothA)) - 2027025.0_8*cothA))/A**8
          case(9)
            values = (-34459425.0_8*R + A*(34459425.0_8*RcothA + A*(-16216200.0_8*R + &
              A*(4729725.0_8*RcothA + A*(-945945.0_8*R + A*(135135.0_8*RcothA + &
              A*(-13860.0_8*R + A*(990.0_8*RcothA + A*(-45.0_8*R + A*RcothA)))))))))/A**9
          case(10)
            values = (654729075.0_8 + A*(A*(310134825.0_8 + A*(A*(18918900.0_8 + A*(A* &
              (315315.0_8 + A*(A*(1485.0_8 + A*(A - 55.0_8*cothA)) - 25740.0_8*cothA)) - &
              2837835.0_8*cothA)) - 91891800.0_8*cothA)) - 654729075.0_8*cothA))/A**10
          case(11)
            values = (-13749310575.0_8*R + A*(13749310575.0_8*RcothA + A*(-6547290750.0_8*R + &
              A*(1964187225.0_8*RcothA + &
              A*(-413513100.0_8*R + A*(64324260.0_8*RcothA + A*(-7567560.0_8*R + A*(675675.0_8*RcothA + &
              A*(-45045.0_8*R + A*(2145.0_8*RcothA + A*(-66.0_8*R + A*RcothA)))))))))))/A**11
          case(12)
            values = (316234143225.0_8 + A*(A*(151242416325.0_8 + A*(A*(9820936125.0_8 + A*(A* &
              (192972780.0_8 + A*(A*(1351350.0_8 + A*(A*(3003.0_8 + A*(A - 78.0_8*cothA)) - 75075.0_8*cothA)) - &
              18378360.0_8*cothA)) - 1571349780.0_8*cothA)) - 45831035250.0_8*cothA)) - &
              316234143225.0_8*cothA))/A**12
          case(13)
            values = (-7905853580625.0_8*R + A*(7905853580625.0_8*RcothA + &
              A*(-3794809718700.0_8*R + A*(1159525191825.0_8*RcothA + A*(-252070693875.0_8*R + &
              A*(41247931725.0_8*RcothA + A*(-5237832600.0_8*R + A*(523783260.0_8*RcothA + &
              A*(-41351310.0_8*R + A*(2552550.0_8*RcothA + A*(-120120.0_8*R + A*(4095.0_8*RcothA + &
              A*(-91.0_8*R + A*RcothA)))))))))))))/A**13
          case(14)
            values = (213458046676875.0_8 + A*(A*(102776096548125.0_8 + A*(A*(6957151150950.0_8 + A* &
              (A*(151242416325.0_8 + A*(A*(1309458150.0_8 + A*(A*(4594590.0_8 + A*(A*(5460.0_8 + A* &
              (A - 105.0_8*cothA)) - 185640.0_8*cothA)) - 87297210.0_8*cothA)) - 15713497800.0_8*cothA)) - &
              1159525191825.0_8*cothA)) - 31623414322500.0_8*cothA)) - 213458046676875.0_8*cothA))/A**14
          case(15)
            values = (-6190283353629375.0_8*R + A*(6190283353629375.0_8*RcothA + &
              A*(-2988412653476250.0_8*R + A*(924984868933125.0_8*RcothA + A*(-205552193096250.0_8*R + &
              A*(34785755754750.0_8*RcothA + A*(-4638100767300.0_8*R + A*(496939367925.0_8*RcothA + &
              A*(-43212118950.0_8*R + A*(3055402350.0_8*RcothA + A*(-174594420.0_8*R + A*(7936110.0_8*RcothA + &
              A*(-278460.0_8*R + A*(7140.0_8*RcothA + A*(-120.0_8*R + A*RcothA)))))))))))))))/A**15
          case(16)
            values = (191898783962510625.0_8 + A*(A*(92854250304440625.0_8 + A*(A*(6474894082531875.0_8 + &
              A*(A*(150738274937250.0_8 + A*(A*(1490818103775.0_8 + A*(A*(6721885170.0_8 + A*(A*(13226850.0_8 + &
              A*(A*(9180.0_8 + A*(A - 136.0_8*cothA)) - 406980.0_8*cothA)) - 333316620.0_8*cothA)) - &
              110430970650.0_8*cothA)) - 16564645597500.0_8*cothA)) - 1109981842719750.0_8*cothA)) - &
              28887988983603750.0_8*cothA)) - 191898783962510625.0_8*cothA))/A**16
          case(17)
            values = (-6332659870762850625.0_8*R + A*(6332659870762850625.0_8*RcothA + &
              A*(-3070380543400170000.0_8*R + A*(959493919812553125.0_8*RcothA + A*(-216659917377028125.0_8*R + &
              A*(37554385678684875.0_8*RcothA + A*(-5179915266025500.0_8*R + A*(581419060472250.0_8*RcothA + &
              A*(-53835098191875.0_8*R + A*(4141161399375.0_8*RcothA + A*(-265034329560.0_8*R + A* &
              (14054850810.0_8*RcothA + A*(-611080470.0_8*R + A*(21366450.0_8*RcothA + A* &
              (-581400.0_8*R + A*(11628.0_8*RcothA + A*(-153.0_8*R + A*RcothA)))))))))))))))))/A**17
          case(18)
            values = (221643095476699771875.0_8 + A*(A*(107655217802968460625.0_8 + &
              A*(A*(7675951358500425000.0_8 + A*(A*(187771928393424375.0_8 + A*(A*(2034966711652875.0_8 + &
              A*(A*(10767019638375.0_8 + A*(A*(28109701620.0_8 + A*(A*(33575850.0_8 + &
              A*(A*(14535.0_8 + A*(A - 171.0_8*cothA)) - 813960.0_8*cothA)) - 1081142370.0_8*cothA)) - &
              602350749000.0_8*cothA)) - 161505294575625.0_8*cothA)) - 21459648959248500.0_8*cothA)) - &
              1343291487737574375.0_8*cothA)) - 33774185977401870000.0_8*cothA)) - &
              221643095476699771875.0_8*cothA))/A**18
          case(19)
            values = (-8200794532637891559375.0_8*R + A*(8200794532637891559375.0_8*RcothA + &
              A*(-3989575718580595893750.0_8*R + A*(1255977541034632040625.0_8*RcothA + &
              A*(-287080580807915895000.0_8*R + A*(50661278966102805000.0_8*RcothA + A*(-7164221267933730000.0_8*R + &
              A*(831561397170879375.0_8*RcothA + A*(-80473683597181875.0_8*R + A*(6557114959770375.0_8*RcothA + &
              A*(-452214824811750.0_8*R + A*(26428139112375.0_8*RcothA + A*(-1305093289500.0_8*R + A* &
              (54057118500.0_8*RcothA + A*(-1853386920.0_8*R + A*(51482970.0_8*RcothA + A* &
              (-1119195.0_8*R + A*(17955.0_8*RcothA + A*(-190.0_8*R + A*RcothA)))))))))))))))))))/A**19
          case(20)
            values = (319830986772877770815625.0_8 + A*(A*(155815096120119939628125.0_8 + &
              A*(A*(11303797869311688365625.0_8 + A*(A*(287080580807915895000.0_8 + &
              A*(A*(3326245588683517500.0_8 + A*(A*(19671344879311125.0_8 + A*(A*(61665657928875.0_8 + &
              A*(A*(100391791500.0_8 + A*(A*(77224455.0_8 + A*(A*(21945.0_8 + A*(A - 210.0_8*cothA)) - &
              1514205.0_8*cothA)) - 3088978200.0_8*cothA)) - 2710578370500.0_8*cothA)) - &
              1192202719958250.0_8*cothA)) - 277187132390293125.0_8*cothA)) - &
              33774185977401870000.0_8*cothA)) - 2009564065655411265000.0_8*cothA)) - &
              49204767195827349356250.0_8*cothA)) - 319830986772877770815625.0_8*cothA))/A**20   
          case default
            values = ONE
        end select
        integrals(l + 1) = integrals(l + 1) + values
      end do
    end function calc_leg_km
       
!===============================================================================
! CALC_INT_PN_EQUI calculates the Legendre expansion of an equiprobable
! angular representation (for one bin). This function returns 
! all n orders. It is intended to operate on one set of points from an equiprob
! bin. Since this function is called repeatedly,
! neither n or x is checked to see if they are in the applicable range. 
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range, 
! the return value will be 1.0_8.
!===============================================================================

    pure function calc_int_pn_equi(n, xlow, xhigh, proby) result(integrals)
      integer, intent(in)  :: n     ! Legendre orders requested (inclusive)
      real(8), intent(in)  :: xlow  ! Low boundary of range of integration
      real(8), intent(in)  :: xhigh ! High boundary of range of integration
      real(8), intent(in)  :: proby ! Probabiltiy of this bin (expect 1/32)
      
      real(8) :: integrals(n + 1) ! The Legendres evaluated at x (0:nth order)
      
      real(8)              :: flow, fhigh  ! the height of the function
      
      flow = proby / (xhigh - xlow)
      fhigh = flow
      
      ! Calculate the moments using the line integral method, but with
      ! flow == fhigh, which forces it to be flat.
      integrals = calc_int_pn_tablelin(n, xlow, xhigh, flow, fhigh)
      
    end function calc_int_pn_equi    
    
!===============================================================================
! CALC_INT_PN_TABLEHIST calculates the Legendre expansion of a tabular
! angular distribution with histogram interpolation. This function returns 
! all n orders. It is intended to operate on one set of points from the table.
! Since this function is called repeatedly,
! neither n or x is checked to see if they are in the applicable range. 
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range, 
! the return value will be 1.0_8.
!===============================================================================

    pure function calc_int_pn_tablehist(n, xlow, xhigh, flow) result(integrals)
      integer, intent(in)  :: n     ! Legendre orders requested (inclusive)
      real(8), intent(in)  :: xlow  ! Low boundary of range of integration
      real(8), intent(in)  :: xhigh ! High boundary of range of integration
      real(8), intent(in)  :: flow  ! Line y-value at low boundary
      
      real(8) :: integrals(n + 1) ! The Legendres evaluated at x (0:nth order)
      
      real(8) :: fhigh ! the height of the function
      
      fhigh = flow
      
      ! Calculate the moments using the line integral method, but with
      ! flow == fhigh, which forces it to be flat.
      integrals = calc_int_pn_tablelin(n, xlow, xhigh, flow, fhigh)
      
    end function calc_int_pn_tablehist

!===============================================================================
! CALC_INT_PN_TABLELIN calculates the Legendre expansion of a tabular angular 
! distribution with linear-linear interpolation. This function returns 
! all n orders. It is intended to operate on one set of points from the table.
! Since this function is called repeatedly,
! neither n or x is checked to see if they are in the applicable range. 
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range, 
! the return value will be 1.0_8.
!===============================================================================
  
    pure function calc_int_pn_tablelin(n, xlow, xhigh, flow, fhigh) &
      result(integrals)
      integer, intent(in)  :: n     ! Legendre orders requested (inclusive)
      real(8), intent(in)  :: xlow  ! Low boundary of range of integration
      real(8), intent(in)  :: xhigh ! High boundary of range of integration
      real(8), intent(in)  :: flow  ! Line y-value at low boundary
      real(8), intent(in)  :: fhigh ! Line y-value at high boundary
      
      real(8) :: integrals(n + 1)       ! The Legendres evaluated at x (0:n)
      
      integer :: l     ! Index of the Legendre order
      real(8) :: values
      
      ! The values below for integrals(l+1) are the integral of (f(x)*P_l(x)) where
      ! f(x) is a straight line connecting flow and fhigh at xlow and xhigh.
      integrals = ZERO
      do l = 0, n
        select case (l)
          case (0)
            values = 0.5_8*((fhigh+flow)*xlow**2-TWO*flow*xhigh*xlow) &
              /(xhigh-xlow)+0.5_8*((fhigh+flow)*xhigh**2-TWO*fhigh*xhigh*xlow) &
              /(xhigh-xlow)
          case (1)
            values = (ONE/6.0_8*((TWO*fhigh + flow)*xhigh**3 - &
              3.0_8*fhigh*xhigh**2*xlow)/(xhigh - xlow) + &
              ONE/6.0_8*((fhigh + TWO*flow)*xlow**3 - 3.0_8*flow*xhigh*xlow**2) &
              /(xhigh - xlow))
          case (2)
            values = ONE/8.0_8*((3.0_8*fhigh+flow)*xhigh**4- &
              2.0_8*(fhigh+flow)*xhigh**2-4.0_8*(fhigh*xhigh**3-fhigh*xhigh)*xlow) &
              /(xhigh-xlow)+ONE/8.0_8*((fhigh+3.0_8*flow)*xlow**4-4.0_8*flow* &
              xhigh*xlow**3-2.0_8*(fhigh+flow)*xlow**2+4.0_8*flow*xhigh*xlow)/ &
              (xhigh-xlow)
          case (3)
            values = (ONE/8.0_8*((4.0_8*fhigh+flow)*xhigh**5-2.0_8* &
              (2.0_8*fhigh+flow)*xhigh**3-(5.0_8*fhigh*xhigh**4-6.0_8*fhigh*xhigh**2)* &
              xlow)/(xhigh-xlow)+ONE/8.0_8*((fhigh+4.0_8*flow)*xlow**5-5.0_8* &
              flow*xhigh*xlow**4-2.0_8*(fhigh+2.0_8*flow)*xlow**3+6.0_8*flow* &
              xhigh*xlow**2)/(xhigh-xlow))
          case (4)
            values = ONE/48.0_8*(7.0_8*(5.0_8*fhigh+flow)*xhigh**6- &
              15.0_8*(3.0_8*fhigh+flow)*xhigh**4+9.0_8*(fhigh+flow)*xhigh**2- &
              6.0_8*(7.0_8*fhigh*xhigh**5-10.0_8*fhigh*xhigh**3+3.0_8*fhigh*xhigh)* &
              xlow)/(xhigh-xlow)+ONE/48.0_8*(7.0_8*(fhigh+5.0_8*flow)*xlow**6- &
              42.0_8*flow*xhigh*xlow**5-15.0_8*(fhigh+3.0_8*flow)*xlow**4+60.0_8* &
              flow*xhigh*xlow**3+9.0_8*(fhigh+flow)*xlow**2-18.0_8*flow*xhigh*xlow) &
              /(xhigh-xlow)
          case (5)
            values = ONE/16.0_8*(3.0_8*(6.0_8*fhigh+flow)*xhigh**7- &
              7.0_8*(4.0_8*fhigh+flow)*xhigh**5+5.0_8*(2.0_8*fhigh+flow)*xhigh**3- &
              (21.0_8*fhigh*xhigh**6-35.0_8*fhigh*xhigh**4+15.0_8*fhigh*xhigh**2)*xlow) &
              /(xhigh-xlow)+ONE/16.0_8*(3.0_8*(fhigh+6.0_8*flow)*xlow**7-21.0_8* &
              flow*xhigh*xlow**6-7.0_8*(fhigh+4.0_8*flow)*xlow**5+35.0_8*flow* &
              xhigh*xlow**4+5.0_8*(fhigh+2.0_8*flow)*xlow**3-15.0_8*flow* &
              xhigh*xlow**2)/(xhigh-xlow)
          case (6)
            values = (ONE/128.0_8*(33.0_8*(7.0_8*fhigh + flow)*xhigh**8 - &
              84.0_8*(5.0_8*fhigh + flow)*xhigh**6 + 70.0_8*(3.0_8*fhigh + flow)* &
              xhigh**4 - 20.0_8*(fhigh + flow)*xhigh**2 - 8.0_8*(33.0_8*fhigh* &
              xhigh**7 - 63.0_8*fhigh*xhigh**5 + 35.0_8*fhigh*xhigh**3 - &
              5.0_8*fhigh*xhigh)*xlow)/(xhigh - xlow) + ONE/128.0_8*(33.0_8* &
              (fhigh + 7.0_8*flow)*xlow**8 - 264.0_8*flow*xhigh*xlow**7 - &
              84.0_8*(fhigh + 5.0_8*flow)*xlow**6 +504.0_8*flow*xhigh*xlow**5 + &
              70.0_8*(fhigh + 3.0_8*flow)*xlow**4 - 280.0_8*flow*xhigh*xlow**3 - &
              20.0_8*(fhigh + flow)*xlow**2 + 40.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (7)
            values = (ONE/384.0_8*(143.0_8*(8.0_8*fhigh + flow)*xhigh**9 - &
              396.0_8*(6.0_8*fhigh + flow)*xhigh**7 + 378.0_8*(4.0_8*fhigh + flow)* &
              xhigh**5 - 140.0_8*(2.0_8*fhigh + flow)*xhigh**3 -3.0_8*(429.0_8* &
              fhigh*xhigh**8 - 924.0_8*fhigh*xhigh**6 + 630.0_8*fhigh*xhigh**4 - &
              140.0_8*fhigh*xhigh**2)*xlow)/(xhigh - xlow) + ONE/384.0_8*(143.0_8* &
              (fhigh +8.0_8*flow)*xlow**9 - 1287.0_8*flow*xhigh*xlow**8 - 396.0_8* &
              (fhigh + 6.0_8*flow)*xlow**7 +2772.0_8*flow*xhigh*xlow**6 + 378.0_8* &
              (fhigh + 4.0_8*flow)*xlow**5-1890.0_8*flow*xhigh*xlow**4 - 140.0_8* &
              (fhigh+2.0_8*flow)*xlow**3 +420.0_8*flow*xhigh*xlow**2)/(xhigh-xlow))
          case (8)
            values = (ONE/256.0_8*(143.0_8*(9.0_8*fhigh + flow)*xhigh**10- &
              429.0_8*(7.0_8*fhigh + flow)*xhigh**8 +462.0_8*(5.0_8*fhigh + flow)* &
              xhigh**6 - 210.0_8*(3.0_8*fhigh + flow)*xhigh**4 + 35.0_8*(fhigh+ &
              flow)*xhigh**2 - 2.0_8*(715.0_8*fhigh*xhigh**9 - 1716.0_8*fhigh* &
              xhigh**7 +1386.0_8*fhigh*xhigh**5 - 420.0_8*fhigh*xhigh**3 + 35.0_8* &
              fhigh*xhigh)*xlow)/(xhigh-xlow) + ONE/256.0_8*(143.0_8*(fhigh + &
              9.0_8*flow)*xlow**10 - 1430.0_8*flow*xhigh*xlow**9 -429.0_8*(fhigh+ &
              7.0_8*flow)*xlow**8 + 3432.0_8*flow*xhigh*xlow**7 + 462.0_8*(fhigh+ &
              5.0_8*flow)*xlow**6 - 2772.0_8*flow*xhigh*xlow**5 - 210.0_8*(fhigh+ &
              3.0_8*flow)*xlow**4 +840.0_8*flow*xhigh*xlow**3 + 35.0_8*(fhigh + &
              flow)*xlow**2 -70.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (9)
            values = (ONE/384.0_8*(143.0_8*(8.0_8*fhigh + flow)*xhigh**9- &
              396.0_8*(6.0_8*fhigh + flow)*xhigh**7 +378.0_8*(4.0_8*fhigh + flow)* &
              xhigh**5 - 140.0_8*(2.0_8*fhigh + flow)*xhigh**3 -3.0_8*(429.0_8* &
              fhigh*xhigh**8 - 924.0_8*fhigh*xhigh**6 + 630.0_8*fhigh*xhigh**4 - &
              140.0_8*fhigh*xhigh**2)*xlow)/(xhigh - xlow) + ONE/384.0_8*(143.0_8* &
              (fhigh +8.0_8*flow)*xlow**9 - 1287.0_8*flow*xhigh*xlow**8 - 396.0_8* &
              (fhigh + 6.0_8*flow)*xlow**7 +2772.0_8*flow*xhigh*xlow**6 + 378.0_8* &
              (fhigh + 4.0_8*flow)*xlow**5 -1890.0_8*flow*xhigh*xlow**4 - 140.0_8* &
              (fhigh+2.0_8*flow)*xlow**3 +420.0_8*flow*xhigh*xlow**2)/(xhigh-xlow))
          case (10)
            values = (ONE/3072.0_8*(4199.0_8*(11.0_8*fhigh + flow)*xhigh**12- &
              14586.0_8*(9.0_8*fhigh +flow)*xhigh**10 + 19305.0_8*(7.0_8*fhigh+ &
              flow)*xhigh**8 - 12012.0_8*(5.0_8*fhigh +flow)*xhigh**6 + 3465.0_8* &
              (3.0_8*fhigh + flow)*xhigh**4 - 378.0_8*(fhigh +flow)*xhigh**2 - &
              12.0_8*(4199.0_8*fhigh*xhigh**11 - 12155.0_8*fhigh*xhigh**9 + &
              12870.0_8*fhigh*xhigh**7 - 6006.0_8*fhigh*xhigh**5 + 1155.0_8*fhigh*xhigh**3 - &
              63.0_8*fhigh*xhigh)*xlow)/(xhigh - xlow) + ONE/3072.0_8*(4199.0_8*(fhigh + &
              11.0_8*flow)*xlow**12 - 50388.0_8*flow*xhigh*xlow**11 - 14586.0_8*(fhigh + &
              9.0_8*flow)*xlow**10 + 145860.0_8*flow*xhigh*xlow**9 + 19305.0_8*(fhigh + &
              7.0_8*flow)*xlow**8 - 154440.0_8*flow*xhigh*xlow**7 - 12012.0_8*(fhigh + &
              5.0_8*flow)*xlow**6 + 72072.0_8*flow*xhigh*xlow**5 + 3465.0_8*(fhigh + &
              3.0_8*flow)*xlow**4-13860.0_8*flow*xhigh*xlow**3 - 378.0_8*(fhigh + &
              flow)*xlow**2 +756.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (11)
            values = (ONE/1024.0_8*(2261.0_8*(12.0_8*fhigh + flow)*xhigh**13 - &
              8398.0_8*(10.0_8*fhigh +flow)*xhigh**11 + 12155.0_8*(8.0_8*fhigh + &
              flow)*xhigh**9 - 8580.0_8*(6.0_8*fhigh +flow)*xhigh**7 + 3003.0_8* &
              (4.0_8*fhigh + flow)*xhigh**5 - 462.0_8*(2.0_8*fhigh +flow)*xhigh**3 - &
              (29393.0_8*fhigh*xhigh**12 - 92378.0_8*fhigh*xhigh**10 +109395.0_8* &
              fhigh*xhigh**8 - 60060.0_8*fhigh*xhigh**6 + 15015.0_8*fhigh*xhigh**4 - &
              1386.0_8*fhigh*xhigh**2)*xlow)/(xhigh - xlow) + ONE/1024.0_8*(2261.0_8*(fhigh + &
              12.0_8*flow)*xlow**13 - 29393.0_8*flow*xhigh*xlow**12 - 8398.0_8*(fhigh + & 
              10.0_8*flow)*xlow**11 + 92378.0_8*flow*xhigh*xlow**10 + 12155.0_8*(fhigh + &
              8.0_8*flow)*xlow**9 - 109395.0_8*flow*xhigh*xlow**8 - 8580.0_8* &
              (fhigh + 6.0_8*flow)*xlow**7+ 60060.0_8*flow*xhigh*xlow**6 + 3003.0_8* &
              (fhigh + 4.0_8*flow)*xlow**5 -15015.0_8*flow*xhigh*xlow**4 - 462.0_8* &
              (fhigh + 2.0_8*flow)*xlow**3 +1386.0_8*flow*xhigh*xlow**2)/(xhigh - xlow))
          case (12)
            values = (ONE/2048.0_8*(7429.0_8*(13.0_8*fhigh + flow)*xhigh**14 - &
              29393.0_8*(11.0_8*fhigh +flow)*xhigh**12 + 46189.0_8*(9.0_8*fhigh + flow)* &
              xhigh**10 - 36465.0_8*(7.0_8*fhigh +flow)*xhigh**8 + 15015.0_8*(5.0_8* &
              fhigh + flow)*xhigh**6 - 3003.0_8*(3.0_8*fhigh +flow)*xhigh**4 + 231.0_8* &
              (fhigh + flow)*xhigh**2 - 2.0_8*(52003.0_8*fhigh*xhigh**13-176358.0_8* &
              fhigh*xhigh**11 + 230945.0_8*fhigh*xhigh**9 - 145860.0_8*fhigh*xhigh**7 +&
              45045.0_8*fhigh*xhigh**5 - 6006.0_8*fhigh*xhigh**3 + 231.0_8*fhigh*xhigh)* &
              xlow)/(xhigh-xlow) + ONE/2048.0_8*(7429.0_8*(fhigh + 13.0_8*flow)*xlow**14 - &
              104006.0_8*flow*xhigh*xlow**13 - 29393.0_8*(fhigh + 11.0_8*flow)*xlow**12 + &
              352716.0_8*flow*xhigh*xlow**11 + 46189.0_8*(fhigh + 9.0_8*flow)*xlow**10 - &
              461890.0_8*flow*xhigh*xlow**9 - 36465.0_8*(fhigh + 7.0_8*flow)*xlow**8 + &
              291720.0_8*flow*xhigh*xlow**7 + 15015.0_8*(fhigh + 5.0_8*flow)*xlow**6 - &
              90090.0_8*flow*xhigh*xlow**5 - 3003.0_8*(fhigh + 3.0_8*flow)*xlow**4 + &
              12012.0_8*flow*xhigh*xlow**3 + 231.0_8*(fhigh + flow)*xlow**2 - &
              462.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (13)
            values = (ONE/6144.0_8*(37145.0_8*(14.0_8*fhigh + flow)*xhigh**15 - 156009.0_8*(12.0_8*fhigh + &
              flow)*xhigh**13 + 264537.0_8*(10.0_8*fhigh + flow)*xhigh**11 - 230945.0_8*(8.0_8*fhigh + &
              flow)*xhigh**9 + 109395.0_8*(6.0_8*fhigh + flow)*xhigh**7 - 27027.0_8*(4.0_8*fhigh + &
              flow)*xhigh**5 + 3003.0_8*(2.0_8*fhigh + flow)*xhigh**3 - 3.0_8*(185725.0_8*fhigh*xhigh**14 &
              - 676039.0_8*fhigh*xhigh**12 + 969969.0_8*fhigh*xhigh**10 - 692835.0_8*fhigh*xhigh**8 + &
              255255.0_8*fhigh*xhigh**6 - 45045.0_8*fhigh*xhigh**4 + &
              3003.0_8*fhigh*xhigh**2)*xlow)/(xhigh - xlow) + ONE/6144.0_8*(37145.0_8*(fhigh + &
              14.0_8*flow)*xlow**15 - 557175.0_8*flow*xhigh*xlow**14 - 156009.0_8*(fhigh + &
              12.0_8*flow)*xlow**13 + 2028117.0_8*flow*xhigh*xlow**12 + 264537.0_8*(fhigh + &
              10.0_8*flow)*xlow**11 - 2909907.0_8*flow*xhigh*xlow**10 - 230945.0_8*(fhigh + &
              8.0_8*flow)*xlow**9 + 2078505.0_8*flow*xhigh*xlow**8 + 109395.0_8*(fhigh + &
              6.0_8*flow)*xlow**7 - 765765.0_8*flow*xhigh*xlow**6 - 27027.0_8*(fhigh + &
              4.0_8*flow)*xlow**5 + 135135.0_8*flow*xhigh*xlow**4 + 3003.0_8*(fhigh + 2.0_8*flow)*xlow**3 &
              - 9009.0_8*flow*xhigh*xlow**2)/(xhigh - xlow))
          case (14)
            values = (ONE/32768.0_8*(334305.0_8*(15.0_8*fhigh + flow)*xhigh**16 - 1485800.0_8*(13.0_8*fhigh + &
              flow)*xhigh**14 + 2704156.0_8*(11.0_8*fhigh + flow)*xhigh**12 - 2586584.0_8*(9.0_8*fhigh + &
              flow)*xhigh**10 + 1385670.0_8*(7.0_8*fhigh + flow)*xhigh**8 - 408408.0_8*(5.0_8*fhigh + &
              flow)*xhigh**6 + 60060.0_8*(3.0_8*fhigh + flow)*xhigh**4 - 3432.0_8*(fhigh + &
              flow)*xhigh**2 - 16.0_8*(334305.0_8*fhigh*xhigh**15 - 1300075.0_8*fhigh*xhigh**13 + &
              2028117.0_8*fhigh*xhigh**11 - 1616615.0_8*fhigh*xhigh**9 + 692835.0_8*fhigh*xhigh**7 - &
              153153.0_8*fhigh*xhigh**5 + 15015.0_8*fhigh*xhigh**3 - &
              429.0_8*fhigh*xhigh)*xlow)/(xhigh - xlow) + ONE/32768.0_8*(334305.0_8*(fhigh + &
              15.0_8*flow)*xlow**16 - 5348880.0_8*flow*xhigh*xlow**15 - 1485800.0_8*(fhigh + &
              13.0_8*flow)*xlow**14 + 20801200.0_8*flow*xhigh*xlow**13 + 2704156.0_8*(fhigh + &
              11.0_8*flow)*xlow**12 - 32449872.0_8*flow*xhigh*xlow**11 - 2586584.0_8*(fhigh + &
              9.0_8*flow)*xlow**10 + 25865840.0_8*flow*xhigh*xlow**9 + 1385670.0_8*(fhigh + &
              7.0_8*flow)*xlow**8 - 11085360.0_8*flow*xhigh*xlow**7 - 408408.0_8*(fhigh + &
              5.0_8*flow)*xlow**6 + 2450448.0_8*flow*xhigh*xlow**5 + 60060.0_8*(fhigh + &
              3.0_8*flow)*xlow**4 - 240240.0_8*flow*xhigh*xlow**3 - 3432.0_8*(fhigh + flow)*xlow**2 + &
              6864.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (15)
            values = (ONE/32768.0_8*(570285.0_8*(16.0_8*fhigh + flow)*xhigh**17 - &
              2674440.0_8*(14.0_8*fhigh + flow)*xhigh**15 + 5200300.0_8*(12.0_8*fhigh + flow)* &
              xhigh**13 - 5408312.0_8*(10.0_8*fhigh+ flow)*xhigh**11 + 3233230.0_8* &
              (8.0_8*fhigh + flow)*xhigh**9 - 1108536.0_8*(6.0_8*fhigh +flow)*xhigh**7 + &
              204204.0_8*(4.0_8*fhigh + flow)*xhigh**5 - 17160.0_8*(2.0_8*fhigh +flow)* &
              xhigh**3 - (9694845.0_8*fhigh*xhigh**16 - 40116600.0_8*fhigh*xhigh**14 + &
              67603900.0_8*fhigh*xhigh**12 - 59491432.0_8*fhigh*xhigh**10+29099070.0_8* &
              fhigh*xhigh**8 - 7759752.0_8*fhigh*xhigh**6 + 1021020.0_8*fhigh*xhigh**4 - &
              51480.0_8*fhigh*xhigh**2)*xlow)/(xhigh - xlow) + ONE/32768.0_8*(570285.0_8*(fhigh + &
              16.0_8*flow)*xlow**17 - 9694845.0_8*flow*xhigh*xlow**16 - 2674440.0_8*(fhigh + &
              14.0_8*flow)*xlow**15 + 40116600.0_8*flow*xhigh*xlow**14 + 5200300.0_8*(fhigh + &
              12.0_8*flow)*xlow**13 - 67603900.0_8*flow*xhigh*xlow**12 - 5408312.0_8*(fhigh + &
              10.0_8*flow)*xlow**11 + 59491432.0_8*flow*xhigh*xlow**10 + 3233230.0_8*(fhigh + &
              8.0_8*flow)*xlow**9 - 29099070.0_8*flow*xhigh*xlow**8 - 1108536.0_8*(fhigh + &
              6.0_8*flow)*xlow**7 + 7759752.0_8*flow*xhigh*xlow**6 + 204204.0_8*(fhigh + &
              4.0_8*flow)*xlow**5 - 1021020.0_8*flow*xhigh*xlow**4 - 17160.0_8*(fhigh + &
              2.0_8*flow)*xlow**3 + 51480.0_8*flow*xhigh*xlow**2)/(xhigh - xlow))
          case (16)
            values = (ONE/196608.0_8*(5892945.0_8*(17.0_8*fhigh + flow)*xhigh**18 - &
              29084535.0_8*(15.0_8*fhigh +flow)*xhigh**16 + 60174900.0_8*(13.0_8*fhigh + flow)*xhigh**14 - &
              67603900.0_8*(11.0_8*fhigh + flow)*xhigh**12 + 44618574.0_8*(9.0_8*fhigh + flow)*xhigh**10 - &
              17459442.0_8*(7.0_8*fhigh + flow)*xhigh**8 + 3879876.0_8*(5.0_8*fhigh + flow)*xhigh**6 - &
              437580.0_8*(3.0_8*fhigh + flow)*xhigh**4 + 19305.0_8*(fhigh + flow)*xhigh**2 - &
              6.0_8*(17678835.0_8*fhigh*xhigh**17 - 77558760.0_8*fhigh*xhigh**15 + &
              140408100.0_8*fhigh*xhigh**13 - 135207800.0_8*fhigh*xhigh**11 + &
              74364290.0_8*fhigh*xhigh**9 - 23279256.0_8*fhigh*xhigh**7 + 3879876.0_8*fhigh*xhigh**5 - &
              291720.0_8*fhigh*xhigh**3 + 6435.0_8*fhigh*xhigh)*xlow)/(xhigh - xlow) + &
              ONE/196608.0_8*(5892945.0_8*(fhigh + 17.0_8*flow)*xlow**18 - &
              106073010.0_8*flow*xhigh*xlow**17 - 29084535.0_8*(fhigh + 15.0_8*flow)*xlow**16 + &
              465352560.0_8*flow*xhigh*xlow**15 + 60174900.0_8*(fhigh + 13.0_8*flow)*xlow**14 - &
              842448600.0_8*flow*xhigh*xlow**13 - 67603900.0_8*(fhigh + 11.0_8*flow)*xlow**12 + &
              811246800.0_8*flow*xhigh*xlow**11 + 44618574.0_8*(fhigh + 9.0_8*flow)*xlow**10 - &
              446185740.0_8*flow*xhigh*xlow**9 - 17459442.0_8*(fhigh + 7.0_8*flow)*xlow**8 + &
              139675536.0_8*flow*xhigh*xlow**7 + 3879876.0_8*(fhigh + 5.0_8*flow)*xlow**6 - &
              23279256.0_8*flow*xhigh*xlow**5 - 437580.0_8*(fhigh + 3.0_8*flow)*xlow**4 + &
              1750320.0_8*flow*xhigh*xlow**3 + 19305.0_8*(fhigh + flow)*xlow**2 - &
              38610.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (17)
            values = (ONE/65536.0_8*(3411705.0_8*(18.0_8*fhigh + flow)*xhigh**19 - &
              17678835.0_8*(16.0_8*fhigh +flow)*xhigh**17 + 38779380.0_8*(14.0_8*fhigh+flow)* &
              xhigh**15 -46802700.0_8*(12.0_8*fhigh + flow)*xhigh**13 + 33801950.0_8*(10.0_8*fhigh + &
              flow)*xhigh**11 - 14872858.0_8*(8.0_8*fhigh + flow)*xhigh**9 + 3879876.0_8*(6.0_8*fhigh + &
              flow)*xhigh**7 - 554268.0_8*(4.0_8*fhigh + flow)*xhigh**5 + 36465.0_8*(2.0_8*fhigh + &
              flow)*xhigh**3 - (64822395.0_8*fhigh*xhigh**18 - 300540195.0_8*fhigh*xhigh**16 + &
              581690700.0_8*fhigh*xhigh**14 - 608435100.0_8*fhigh*xhigh**12 + &
              371821450.0_8*fhigh*xhigh**10 - 133855722.0_8*fhigh*xhigh**8 + &
              27159132.0_8*fhigh*xhigh**6 - 2771340.0_8*fhigh*xhigh**4 + &
              109395.0_8*fhigh*xhigh**2)*xlow)/(xhigh - xlow) + ONE/65536.0_8*(3411705.0_8*(fhigh + &
              18.0_8*flow)*xlow**19 - 64822395.0_8*flow*xhigh*xlow**18 - 17678835.0_8*(fhigh + &
              16.0_8*flow)*xlow**17 + 300540195.0_8*flow*xhigh*xlow**16 + 38779380.0_8*(fhigh + &
              14.0_8*flow)*xlow**15 - 581690700.0_8*flow*xhigh*xlow**14 - 46802700.0_8*(fhigh + &
              12.0_8*flow)*xlow**13 + 608435100.0_8*flow*xhigh*xlow**12 + 33801950.0_8*(fhigh + &
              10.0_8*flow)*xlow**11 - 371821450.0_8*flow*xhigh*xlow**10 - 14872858.0_8*(fhigh + &
              8.0_8*flow)*xlow**9 + 133855722.0_8*flow*xhigh*xlow**8 + 3879876.0_8*(fhigh + &
              6.0_8*flow)*xlow**7 - 27159132.0_8*flow*xhigh*xlow**6 - 554268.0_8*(fhigh + &
              4.0_8*flow)*xlow**5 + 2771340.0_8*flow*xhigh*xlow**4 + 36465.0_8*(fhigh + &
              2.0_8*flow)*xlow**3 - 109395.0_8*flow*xhigh*xlow**2)/(xhigh - xlow))
          case (18)
            values = (ONE/262144.0_8*(23881935.0_8*(19.0_8*fhigh + flow)*xhigh**20 - &
              129644790.0_8*(17.0_8*fhigh +flow)*xhigh**18 + 300540195.0_8*(15.0_8*fhigh + flow)*xhigh**16 - &
              387793800.0_8*(13.0_8*fhigh + flow)*xhigh**14 + 304217550.0_8*(11.0_8*fhigh + &
              flow)*xhigh**12 - 148728580.0_8*(9.0_8*fhigh + flow)*xhigh**10 + 44618574.0_8*(7.0_8*fhigh &
              + flow)*xhigh**8 - 7759752.0_8*(5.0_8*fhigh + flow)*xhigh**6 + 692835.0_8*(3.0_8*fhigh + &
              flow)*xhigh**4 - 24310.0_8*(fhigh + flow)*xhigh**2 - &
              4.0_8*(119409675.0_8*fhigh*xhigh**19 - 583401555.0_8*fhigh*xhigh**17 + &
              1202160780.0_8*fhigh*xhigh**15 - 1357278300.0_8*fhigh*xhigh**13 + &
              912652650.0_8*fhigh*xhigh**11 - 371821450.0_8*fhigh*xhigh**9 + &
              89237148.0_8*fhigh*xhigh**7 - 11639628.0_8*fhigh*xhigh**5 + 692835.0_8*fhigh*xhigh**3 - &
              12155.0_8*fhigh*xhigh)*xlow)/(xhigh - xlow) + ONE/262144.0_8*(23881935.0_8*(fhigh + &
              19.0_8*flow)*xlow**20 - 477638700.0_8*flow*xhigh*xlow**19 - 129644790.0_8*(fhigh + &
              17.0_8*flow)*xlow**18 + 2333606220.0_8*flow*xhigh*xlow**17 + 300540195.0_8*(fhigh + &
              15.0_8*flow)*xlow**16 - 4808643120.0_8*flow*xhigh*xlow**15 - 387793800.0_8*(fhigh + &
              13.0_8*flow)*xlow**14 + 5429113200.0_8*flow*xhigh*xlow**13 + 304217550.0_8*(fhigh + &
              11.0_8*flow)*xlow**12 - 3650610600.0_8*flow*xhigh*xlow**11 - 148728580.0_8*(fhigh + &
              9.0_8*flow)*xlow**10 + 1487285800.0_8*flow*xhigh*xlow**9 + 44618574.0_8*(fhigh + &
              7.0_8*flow)*xlow**8 - 356948592.0_8*flow*xhigh*xlow**7 - 7759752.0_8*(fhigh + &
              5.0_8*flow)*xlow**6 + 46558512.0_8*flow*xhigh*xlow**5 + 692835.0_8*(fhigh + &
              3.0_8*flow)*xlow**4 - 2771340.0_8*flow*xhigh*xlow**3 - 24310.0_8*(fhigh + flow)*xlow**2 &
              + 48620.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case (19)
            values = (ONE/786432.0_8*(126233085.0_8*(20.0_8*fhigh + flow)*xhigh**21- &
              716458050.0_8*(18.0_8*fhigh +flow)*xhigh**19 + 1750204665.0_8*(16.0_8* &
              fhigh + flow)*xhigh**17 -2404321560.0_8*(14.0_8*fhigh + flow)*xhigh**15 + &
              2035917450.0_8*(12.0_8*fhigh +flow)*xhigh**13 - 1095183180.0_8*(10.0_8*fhigh+ &
              flow)*xhigh**11 +371821450.0_8*(8.0_8*fhigh + flow)*xhigh**9 - &
              76488984.0_8*(6.0_8*fhigh + flow)*xhigh**7 +8729721.0_8*(4.0_8*fhigh + flow)* &
              xhigh**5 - 461890.0_8*(2.0_8*fhigh + flow)*xhigh**3 -3.0_8*(883631595.0_8* &
              fhigh*xhigh**20 - 4537567650.0_8*fhigh*xhigh**18 +9917826435.0_8*fhigh* &
              xhigh**16 - 12021607800.0_8*fhigh*xhigh**14 +8822308950.0_8*fhigh*xhigh**12 - &
              4015671660.0_8*fhigh*xhigh**10 +1115464350.0_8*fhigh*xhigh**8 - 178474296.0_8* &
              fhigh*xhigh**6 +14549535.0_8*fhigh*xhigh**4 - 461890.0_8*fhigh*xhigh**2)*xlow)/ &
              (xhigh - xlow) +ONE/786432.0_8*(126233085.0_8*(fhigh + 20.0_8*flow)*xlow**21 - &
              2650894785.0_8*flow*xhigh*xlow**20 - 716458050.0_8*(fhigh + 18.0_8*flow)*xlow**19 + &
              13612702950.0_8*flow*xhigh*xlow**18 + 1750204665.0_8*(fhigh + 16.0_8*flow)*xlow**17 - &
              29753479305.0_8*flow*xhigh*xlow**16 - 2404321560.0_8*(fhigh + 14.0_8*flow)*xlow**15 + &
              36064823400.0_8*flow*xhigh*xlow**14 + 2035917450.0_8*(fhigh + 12.0_8*flow)*xlow**13 - & 
              26466926850.0_8*flow*xhigh*xlow**12 - 1095183180.0_8*(fhigh + 10.0_8*flow)*xlow**11 + &
              12047014980.0_8*flow*xhigh*xlow**10 + 371821450.0_8*(fhigh + 8.0_8*flow)*xlow**9 - &
              3346393050.0_8*flow*xhigh*xlow**8 - 76488984.0_8*(fhigh + 6.0_8*flow)*xlow**7 + &
              535422888.0_8*flow*xhigh*xlow**6 + 8729721.0_8*(fhigh + 4.0_8*flow)*xlow**5 - &
              43648605.0_8*flow*xhigh*xlow**4 - 461890.0_8*(fhigh + 2.0_8*flow)*xlow**3 + &
              1385670.0_8*flow*xhigh*xlow**2)/(xhigh - xlow))
          case (20)
            values = (ONE/524288.0_8*(149184555.0_8*(21.0_8*fhigh + flow)*xhigh**22- &
              883631595.0_8*(19.0_8*fhigh +flow)*xhigh**20 + 2268783825.0_8*(17.0_8* &
              fhigh + flow)*xhigh**18 -3305942145.0_8*(15.0_8*fhigh + flow)*xhigh**16 + &
              3005401950.0_8*(13.0_8*fhigh +flow)*xhigh**14 - 1764461790.0_8*(11.0_8* &
              fhigh + flow)*xhigh**12 +669278610.0_8*(9.0_8*fhigh + flow)*xhigh**10 - &
              159352050.0_8*(7.0_8*fhigh + flow)*xhigh**8+ 22309287.0_8*(5.0_8*fhigh + flow)* &
              xhigh**6 - 1616615.0_8*(3.0_8*fhigh + flow)*xhigh**4 +46189.0_8*(fhigh + flow)* &
              xhigh**2 - 2.0_8*(1641030105.0_8*fhigh*xhigh**21 -8836315950.0_8*fhigh* &
              xhigh**19 + 20419054425.0_8*fhigh*xhigh**17 -26447537160.0_8*fhigh* &
              xhigh**15 + 21037813650.0_8*fhigh*xhigh**13 -10586770740.0_8*fhigh* &
              xhigh**11 + 3346393050.0_8*fhigh*xhigh**9 -637408200.0_8*fhigh*xhigh**7 + &
              66927861.0_8*fhigh*xhigh**5 - 3233230.0_8*fhigh*xhigh**3+ 46189.0_8*fhigh*xhigh)* &
              xlow)/(xhigh - xlow) + ONE/524288.0_8*(149184555.0_8*(fhigh +21.0_8*flow)*xlow**22- &
              3282060210.0_8*flow*xhigh*xlow**21 - 883631595.0_8*(fhigh +19.0_8*flow)*xlow**20 + &
              17672631900.0_8*flow*xhigh*xlow**19 + 2268783825.0_8*(fhigh +17.0_8*flow)*xlow**18 - &
              40838108850.0_8*flow*xhigh*xlow**17 - 3305942145.0_8*(fhigh +15.0_8*flow)*xlow**16 + &
              52895074320.0_8*flow*xhigh*xlow**15 + 3005401950.0_8*(fhigh +13.0_8*flow)*xlow**14 - &
              42075627300.0_8*flow*xhigh*xlow**13 - 1764461790.0_8*(fhigh +11.0_8*flow)*xlow**12 + &
              21173541480.0_8*flow*xhigh*xlow**11 + 669278610.0_8*(fhigh+9.0_8*flow)*xlow**10 - &
              6692786100.0_8*flow*xhigh*xlow**9 - 159352050.0_8*(fhigh +7.0_8*flow)*xlow**8 + &
              1274816400.0_8*flow*xhigh*xlow**7 + 22309287.0_8*(fhigh +5.0_8*flow)*xlow**6 - &
              133855722.0_8*flow*xhigh*xlow**5 - 1616615.0_8*(fhigh +3.0_8*flow)*xlow**4 + &
              6466460.0_8*flow*xhigh*xlow**3 + 46189.0_8*(fhigh + flow)*xlow**2- &
              92378.0_8*flow*xhigh*xlow)/(xhigh - xlow))
          case default
            values = ONE
        end select
        integrals(l + 1) = integrals(l + 1) + values
      end do    
    end function calc_int_pn_tablelin

!===============================================================================
! CALC_PN calculates the n-th order Legendre polynomial at the value of x.
! Since this function is called repeatedly during the neutron transport process,
! neither n or x is checked to see if they are in the applicable range. 
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range, 
! the return value will be 1.0.  This is to replace the one in math.F90 when
! done. But is kept here to avoid merge issues in the mean time.
!===============================================================================
  
  pure function calc_pn(n,x) result(pnx)

    integer, intent(in) :: n   ! Legendre order requested
    real(8), intent(in) :: x   ! Independent variable the Legendre is to be 
                               ! evaluated at; x must be in the domain [-1,1]
    real(8)             :: pnx ! The Legendre poly of order n evaluated at x
    
    select case(n)
    case(1)
      pnx = x
    case(2)
      pnx = 1.5_8 * x * x - 0.5_8
    case(3)
      pnx = 2.5_8 * x * x * x - 1.5_8 * x
    case(4)
      pnx = 4.375_8 * (x ** 4) - 3.75_8 * x * x + 0.375_8
    case(5)
      pnx = 7.875_8 * (x ** 5) - 8.75_8 * x * x * x + 1.875 * x
    case(6)
      pnx = 14.4375_8 * (x ** 6) - 19.6875_8 * (x ** 4) + &
        6.5625_8 * x * x - 0.3125_8
    case(7) 
      pnx = 26.8125_8 * (x ** 7) - 43.3125_8 * (x ** 5) + &
        19.6875_8 * x * x * x - 2.1875_8 * x
    case(8)
      pnx = 50.2734375_8 * (x ** 8) - 93.84375_8 * (x ** 6) + &
        54.140625 * (x ** 4) - 9.84375_8 * x * x + 0.2734375_8
    case(9)
      pnx = 94.9609375_8 * (x ** 9) - 201.09375_8 * (x ** 7) + &
        140.765625_8 * (x ** 5) - 36.09375_8 * x * x * x + 2.4609375_8 * x
    case(10)
      pnx = 180.42578125_8 * (x ** 10) - 427.32421875_8 * (x ** 8) + &
        351.9140625_8 * (x ** 6) - 117.3046875_8 * (x ** 4) + &
        13.53515625_8 * x * x - 0.24609375_8
    ! Cases 11-20 come from Mathematica
    case(11) 
      pnx = (-2.70703125_8 * x + 58.65234375_8 * x**3 - 351.9140625_8 * x**5 + &
        854.6484375_8 * x**7 - 902.12890625_8 * x**9 + 344.44921875_8 * x**11)
    case(12)
      pnx = (0.225585937_8 - 17.595703125_8 * x**2 + 219.946289062_8 * x**4 - &
        997.08984375_8 * x**6 + 2029.790039062_8 * x**8 - 1894.470703125_8 * x**10 + &
        660.194335937_8 * x**12)
    case(13)
      pnx = (2.932617187_8 * x - 87.978515625_8 * x**3 + 747.817382812_8 * x**5 - &
        2706.38671875_8 * x**7 + 4736.176757812_8 * x**9 - 3961.166015625_8 * x**11 + &
        1269.604492187_8 * x**13)
    case(14)
      pnx = (-0.209472656_8 + 21.994628906_8 * x**2 - 373.908691406_8 * x**4 + &
        2368.088378906_8 * x**6 - 7104.265136719_8 * x**8 + 10893.206542969_8 * x**10 &
        - 8252.429199219_8 * x**12 + 2448.522949219_8 * x**14)
    case(15)
      pnx = (-3.142089844_8 * x + 124.636230469_8 * x**3 - 1420.853027344_8 * x**5 + &
        7104.265136719_8 * x**7 - 18155.344238281_8 * x**9 + 24757.287597656_8 * x**11 - &
        17139.660644531_8 * x**13 + 4733.811035156_8 * x**15)
    case(16)
      pnx = (0.196380615_8 - 26.707763672_8 * x**2 + 592.022094727_8 * x**4 - &
        4972.985595703_8 * x**6 + 20424.762268066_8 * x**8 - 45388.360595703_8 * x**10 + &
        55703.897094727_8 * x**12 - 35503.582763672_8 * x**14 + 9171.758880615_8 * x**16)
    case(17)
      pnx = (3.338470459_8 * x - 169.149169922_8 * x**3 + 2486.492797852_8 * x**5 - &
        16339.809814453_8 * x**7 + 56735.450744629_8 * x**9 - 111407.794189453_8 * x**11 + &
        124262.539672852_8 * x**13 - 73374.071044922_8 * x**15 + 17804.002532959_8 * x**17)
    case(18)
      pnx = (-0.185470581_8 + 31.71546936_8 * x**2 - 888.03314209_8 * x**4 + &
        9531.555725098_8 * x**6 - 51061.905670166_8 * x**8 + 153185.717010498_8 * x**10 - &
        269235.502624512_8 * x**12 + 275152.766418457_8 * x**14 - 151334.021530151_8 * x**16 + &
        34618.893814087_8 * x**18)
    case(19)
      pnx = (-3.52394104_8 * x + 222.008285522_8 * x**3 - 4084.952453613_8 * x**5 + &
        34041.270446777_8 * x**7 - 153185.717010498_8 * x**9 + 403853.253936768_8 * x**11 - &
        642023.121643066_8 * x**13 + 605336.086120605_8 * x**15 - 311570.044326782_8 * x**17 + &
        67415.740585327_8 * x**19)
    case(20)
      pnx = (0.176197052_8 - 37.00138092_8 * x**2 + 1276.547641754_8 * x**4 - &
        17020.635223389_8 * x**6 + 114889.287757874_8 * x**8 - 444238.579330444_8 * x**10 + &
        1043287.572669983_8 * x**12 - 1513340.215301514_8 * x**14 + 1324172.688388824_8 * x**16 - &
        640449.535560608_8 * x**18 + 131460.694141388_8 * x**20)
    case default
      pnx = ONE ! correct for case(0), incorrect for the rest
    end select
  
  end function calc_pn

end module legendre
