module legendre

  use ace_header
  use constants

  implicit none

  contains

!===============================================================================
! CALC_INT_PN_TABLELIN calculates the Legendre expansion of a tabular angular
! distribution with linear-linear interpolation. This function returns
! all n orders. It is intended to operate on one set of points from the table.
! Since this function is called repeatedly,
! neither n or x is checked to see if they are in the applicable range.
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=10. If x is outside of the range, the return
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

      real(8) :: integrals(n)       ! The Legendres evaluated at x (0:n)

      integer :: l     ! Index of the Legendre order
      real(8) :: values

      ! The values below for integrals(l+1) are the integral of (f(x)*P_l(x))
      ! where f(x) is a straight line between flow and fhigh at xlow and xhigh
      integrals = ZERO

      ! Sometimes xlow and xhigh can be very near each other (on the order of
      ! machine precision) and this is a perfectly valid situation, however,
      ! these cases lead to division by zero errors when the actual result should
      ! approach zero. In that case, just return zero (which integrals has
      ! alrady been set to.
      if (xhigh-xlow < FP_PRECISION) return

      do l = 0, n - 1
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
    case(0)
      pnx = ONE
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
