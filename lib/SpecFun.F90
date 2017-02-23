
module Constants
  use iso_fortran_env, only: error_unit, dp => real64
  implicit none
  real (dp), parameter :: Pi    = 3.141592653589793_dp , ExpEuler = 1.781072417990198_dp,&
                          Zeta3 = 1.2020569031595942   , Zeta2 = 1.6449340668482262_dp,  &
                          Euler = 0.5772156649015329_dp, l2  = 0.6931471805599453_dp,    &
                          Pio2  = 1.5707963267948966_dp, Pi2 = 9.869604401089358_dp ,    &
                          sr2   = 1.4142135623730951_dp, prec = 1e-10_dp

    !**************************************************************
    !>
    !  Machine constants (replaces the old SLATEC [D1MACH](http://www.netlib.org/slatec/src/d1mach.f) function)
    !
    !  The traditional D1MACH constants are:
    !  * _dpD1MACH( 1) = B**(EMIN-1)_dp,           the smallest positive magnitude.
    !  * _dpD1MACH( 2) = B**EMAX*(1 - B**(-T))_dp, the largest magnitude.
    !  * _dpD1MACH( 3) = B**(-T)_dp,               the smallest relative spacing.
    !  * _dpD1MACH( 4) = B**(1-T)_dp,              the largest relative spacing.
    !  * _dpD1MACH( 5) = LOG10(B)_dp

    !**************************************************************

  real (dp), dimension(5), parameter :: d1mach = [  tiny(1.0_dp), huge(1.0_dp), &
           real(radix(1.0_dp),dp)**(-digits(1.0_dp)), epsilon(1.0_dp), &
           log10(real(radix(1.0_dp),dp)) ]

end module Constants

!ccccccccccccccc

module Legendre
  use Constants, only: dp
contains
  function LegendreList(n,x) result(list)
    integer      , intent(in) :: n
    real (dp)    , intent(in) :: x
    real (dp), dimension(0:n) :: list
    integer                   :: i

    list(:1) = [ 1._dp, x ]

    do i = 2, n
      list(i) = ( (2 * i - 1) * x * list(i - 1) - (i - 1) * list(i - 2) )/i
    end do

  end function LegendreList

!ccccccccccccccc

  function QLegendreList(n,x) result(list)
    integer      , intent(in) :: n
    real (dp)    , intent(in) :: x
    real (dp), dimension(0:n) :: list
    integer                   :: i

    list(:1) = [ 0._dp, 1._dp ]; list(2:) = 0

    do i = 2, n
      list(i) = (  (2 * i - 1) * ( x * list(i - 1) - (-1)**i ) - &
                   (i - 1) * list(i - 2)  )/i
    end do

  end function QLegendreList

end module Legendre

!ccccccccccccccc

subroutine dfZero(F, B, C, R, RE, AE, IFLAG)
  use constants, only: dp, d1mach
!! DFZERO finds a zero of a function in a given interval.
!
!***PURPOSE  Search for a zero of a function F(X) in a given interval
!            (B,C).  It is designed primarily for problems where F(B)
!            and F(C) have opposite signs.
!***LIBRARY   SLATEC
!***CATEGORY  F1B
!***TYPE      DOUBLE PRECISION (FZERO-S, DFZERO-D)
!***KEYWORDS  BISECTION, NONLINEAR, ROOTS, ZEROS
!***AUTHOR  Shampine, L. F., (SNLA)
!           Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     DFZERO searches for a zero of a DOUBLE PRECISION function F(X)
!     between the given DOUBLE PRECISION values B and C until the width
!     of the interval (B,C) has collapsed to within a tolerance
!     specified by the stopping criterion,
!        ABS(B-C)  <=  2.*(RW*ABS(B)+AE).
!     The method used is an efficient combination of bisection and the
!     secant rule and is due to T. J. Dekker.
!
!     Description Of Arguments
!
!   F     :EXT   - Name of the DOUBLE PRECISION external function.  This
!                  name must be in an EXTERNAL statement in the calling
!                  program.  F must be a function of one DOUBLE
!                  PRECISION argument.
!
!   B     :INOUT - One end of the DOUBLE PRECISION interval (B,C).  The
!                  value returned for B usually is the better
!                  approximation to a zero of F.
!
!   C     :INOUT - The other end of the DOUBLE PRECISION interval (B,C)
!
!   R     :IN    - A (better) DOUBLE PRECISION guess of a zero of F
!                  which could help in speeding up convergence.  If F(B)
!                  and F(R) have opposite signs, a root will be found in
!                  the interval (B,R);  if not, but F(R) and F(C) have
!                  opposite signs, a root will be found in the interval
!                  (R,C);  otherwise, the interval (B,C) will be
!                  searched for a possible root.  When no better guess
!                  is known, it is recommended that R be set to B or C,
!                  since if R is not interior to the interval (B,C), it
!                  will be ignored.
!
!   RE    :IN    - Relative error used for RW in the stopping criterion.
!                  If the requested RE is less than machine precision,
!                  then RW is set to approximately machine precision.
!
!   AE    :IN    - Absolute error used in the stopping criterion.  If
!                  the given interval (B,C) contains the origin, then a
!                  nonzero value should be chosen for AE.
!
!   IFLAG :OUT   - A status code.  User must check IFLAG after each
!                  call.  Control returns to the user from DFZERO in all
!                  cases.
!
!                1  B is within the requested tolerance of a zero.
!                   The interval (B,C) collapsed to the requested
!                   tolerance, the function changes sign in (B,C), and
!                   F(X) decreased in magnitude as (B,C) collapsed.
!
!                2  F(B) = 0.  However, the interval (B,C) may not have
!                   collapsed to the requested tolerance.
!
!                3  B may be near a singular point of F(X).
!                   The interval (B,C) collapsed to the requested tol-
!                   erance and the function changes sign in (B,C), but
!                   F(X) increased in magnitude as (B,C) collapsed, i.e.
!                     ABS(F(B out))  >  MAX(ABS(F(B in)),ABS(F(C in)))
!
!                4  No change in sign of F(X) was found although the
!                   interval (B,C) collapsed to the requested tolerance.
!                   The user must examine this case and decide whether
!                   B is near a local minimum of F(X), or B is near a
!                   zero of even multiplicity, or neither of these.
!
!                5  Too many ( >  500) function evaluations used.
!
!***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
!                 code, Report SC-TM-70-631, Sandia Laboratories,
!                 September 1970.
!               T. J. Dekker, Finding a zero by means of successive
!                 linear interpolation, Constructive Aspects of the
!                 Fundamental Theorem of Algebra, edited by B. Dejon
!                 and P. Henrici, Wiley-Interscience, 1969.
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   700901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DFZERO
  real (dp), intent(inout) :: B, C
  real (dp), intent(in   ) :: R, RE, AE
  integer        , intent(out  ) :: IFLAG

  real (dp), external      :: F

  real (dp) :: A, ACBS, ACMB, AW, CMB, ER, FA, FB, FC, FX, FZ, P, Q, RW, T, TOL, Z
  INTEGER         :: IC, KOUNT
!
!***FIRST EXECUTABLE STATEMENT  DFZERO
!
!   ER is two times the computer unit roundoff value which is defined
!   here by the function D1MACH.
!
  ER = 2 * D1MACH(4)
!
!   Initialize.
!
  Z = R
  if (R  <=  MIN(B,C)  .OR.  R  >=  MAX(B,C)) Z = C
  RW = MAX(RE,ER)
  AW = MAX(AE,0._dp)
  IC = 0
  T = Z
  FZ = F(T)
  FC = FZ
  T = B
  FB = F(T)
  KOUNT = 2
  if (SIGN(1.0_dp,FZ)  ==  SIGN(1.0_dp,FB)) go to 1
  C = Z
  go to 2
    1 if (Z  ==  C) go to 2
  T = C
  FC = F(T)
  KOUNT = 3
  if (SIGN(1.0_dp,FZ)  ==  SIGN(1.0_dp,FC)) go to 2
  B = Z
  FB = FZ
    2 A = C
  FA = FC
  ACBS = ABS(B-C)
  FX = MAX(ABS(FB),ABS(FC))
!
    3 if (ABS(FC)  >=  ABS(FB)) go to 4
!
!   Perform interchange.
!
  A = B
  FA = FB
  B = C
  FB = FC
  C = A
  FC = FA
!
    4 CMB = 0.5_dp*(C-B)
  ACMB = ABS(CMB)
  TOL = RW*ABS(B) + AW
!
!   Test stopping criterion and function count.
!
  if (ACMB  <=  TOL) go to 10
  if (FB  ==  0._dp) go to 11
  if (KOUNT  >=  500) go to 14
!
!   Calculate new iterate implicitly as B+P/Q, where we arrange
!   P  >=  0.  The implicit form is used to prevent overflow.
!
  P = (B-A)*FB
  Q = FA - FB
  if (P  >=  0._dp) go to 5
  P = -P
  Q = -Q
!
!   Update A and check for satisfactory reduction in the size of the
!   bracketing interval.  If not, perform bisection.
!
    5 A = B
  FA = FB
  IC = IC + 1
  if (IC  <  4) go to 6
  if (8.0_dp*ACMB  >=  ACBS) go to 8
  IC = 0
  ACBS = ACMB
!
!   Test for too small a change.
!
    6 if (P  >  ABS(Q)*TOL) go to 7
!
!   Increment by TOLerance.
!
  B = B + SIGN(TOL,CMB)
  go to 9
!
!   Root ought to be between B and (C+B)/2.
!
    7 if (P  >=  CMB*Q) go to 8
!
!   Use secant rule.
!
  B = B + P/Q
  go to 9
!
!   Use bisection (C+B)/2.
!
    8 B = B + CMB
!
!   Have completed computation for new iterate B.
!
    9 T = B
  FB = F(T)
  KOUNT = KOUNT + 1
!
!   Decide whether next step is interpolation or extrapolation.
!
  if (SIGN(1.0_dp,FB) /=  SIGN(1.0_dp,FC)) go to 3
  C = A
  FC = FA
  go to 3
!
!   Finished.  Process results for proper setting of IFLAG.
!
   10 if (SIGN(1.0_dp,FB)  ==  SIGN(1.0_dp,FC)) go to 13
  if (ABS(FB)  >  FX) go to 12
  IFLAG = 1
  return
   11 IFLAG = 2
  return
   12 IFLAG = 3
  return
   13 IFLAG = 4
  return
   14 IFLAG = 5
  return
end subroutine dfZero

!ccccccccccccccc

module carlson_elliptic_module
  use constants, only: dp, d1mach; use iso_fortran_env, only: error_unit

  implicit none

  private

!*******************************************************************************
!> author: Jacob Williams
!  license: BSD
!  date: 2/7/2016
!
!  Carlson symmetric forms of elliptic integrals.
!
!  These routines are refactored versions of the ones from [SLATEC](http://www.netlib.org/slatec/).
!  They have been converted into modern Fortran, and the documentation has been
!  converted to FORD syntax.

  public :: drc, drd, drf, drj

  contains

!*******************************************************************************
!>
!  Compute an approximation of the Carlson elliptic integral:
!  $$ R_C(x,y) = \frac{1}{2} \int_{0}^{\infty} (t+x)^{-1/2}(t+y)^{-1} dt $$
!  where \(x\ge0\) and \(y>0\).
!
!  The duplication theorem is iterated until the variables are nearly equal,
!  and the function is then expanded in Taylor series to fifth order.
!  Logarithmic, inverse circular, and inverse hyperbolic functions can be
!  expressed in terms of DRC.
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### DRC special cases
!
!  $$
!    \begin{array}{rcll}
!     R_C(x,x+z) + R_C(y,y+z) &=& R_C(0,z)             & x>0, y>0, ~\mathrm{and}~ z>0 ~\mathrm{and}~ x y = z^2 \\
!     R_C(0,1/4)              &=& R_C(1/16,1/8) = \pi  & \\
!     R_C(9/4,2)              &=& \ln(2)               &
!    \end{array}
!  $$
!
!### Special functions via DRC
!
!  $$
!  \begin{array}{rll}
!     \ln(x)        &= (x-1) R_C \left( \left( \frac{1+x}{2} \right)^2, x \right) & x>0 \\
!     \sin^{-1}(x)  &=  x R_C ( (1-x)^2 ,1 ) & -1 \le x \le 1 \\
!     \cos^{-1}(x)  &= \sqrt{1-x^2} R_C(x^2,1)  & 0 \le x \le 1 \\
!     \tan^{-1}(x)  &= x R_C(1,1+x^2) & -\infty \lt x \lt \infty \\
!     \cot^{-1}(x)  &= R_C(x^2, x^2+1) & 0 \le x \lt \infty \\
!     \sinh^{-1}(x) &= x R_C(1+x^2, 1) & -\infty \lt x \lt \infty \\
!     \cosh^{-1}(x) &= \sqrt{x^2-1} R_C(x^2,1)  & x \ge 1 \\
!     \tanh^{-1}(x) &= x R_C(1,1-x^2) & -1 < x < 1 \\
!     \coth^{-1}(x) &= R_C(x^2,x^2-1) & x > 1 \\
!  \end{array}
!  $$
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drc.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(dp) function drc(x,y,ier) ! NAG: S21BAF

    implicit none

    real(dp), intent(in ) :: x    !! nonnegative variable
    real(dp), intent(in ) :: y    !! positive variable
    integer , intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! _dpIER = 0_dp: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! _dpIER > 0_dp: Abnormal termination of the routine:
                                !! _dpIER = 1_dp: _dpx<0 or y<=0_dp
                                !! _dpIER = 2_dp: _dpx+y<LOLIM_dp
                                !! _dpIER = 3_dp: _dpmax(x,y) > UPLIM_dp

    character(len=16) :: xern3 , xern4 , xern5
    real (dp) :: lamda, mu , s , sn , xn , yn

    real (dp), parameter :: errtol = (d1mach(3)/16)**(1.0_dp/6)
        !! Determines the accuracy of the answer.
        !!
        !! The value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! machine precision.
        !!
        !! Relative error due to truncation is less than
        !! _dp16 * ERRTOL ** 6 / (1 - 2 * ERRTOL)_dp.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than):
        !! (1.0e-3, 2.0e-17),
        !! (3.0e-3, 2.0e-14),
        !! (1.0e-2, 2.0e-11),
        !! (3.0e-2, 2.0e-8),
        !! (1.0e-1, 2.0e-5)
        !!
        !! The accuracy of the computed approximation to the inte-
        !! gral can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the trunca-
        !! tion error there will be round-off error, but in prac-
        !! tice the total error from both sources is usually less
        !! than the amount given in the table.
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.
    real(dp), parameter :: lolim  = 5 * d1mach(1) !! Lower limit of valid arguments
    real(dp), parameter :: uplim  = d1mach(2)/5 !! Upper limit of valid arguments
    real(dp), parameter :: c1     = 1._dp/7
    real(dp), parameter :: c2     = 9._dp/22

    !initialize:
    drc = 0

    ! check for errors:
    if ( x < 0 .or. y <= 0 ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write(error_unit,'(a)') &
            'drc: x<0 .or. y<=0 where x = '//xern3//' and y = '//xern4
        return
    endif

    if ( max(x,y) > uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') uplim
        write(error_unit,'(a)') &
            'drc: max(x,y)>uplim where x = '//&
            xern3//' y = '//xern4//' and uplim = '//xern5
        return
    endif

    if ( x+y<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') lolim
        write(error_unit,'(a)') &
            'drc: x+y<lolim where x = '//xern3//&
            ' y = '//xern4//' and lolim = '//xern5
        return
    endif

    ier = 0;  xn = x;  yn = y

    do
        mu = (xn + yn + yn)/3
        sn = (yn + mu)/mu - 2
        if ( abs(sn) < errtol ) exit
        lamda = 2 * sqrt(xn)*sqrt(yn) + yn
        xn = (xn+lamda)/4;  yn = (yn+lamda)/4
    end do

    s = sn * sn * (0.3_dp + sn*(c1+sn*(0.375_dp + sn*c2)));  drc = (1 + s)/sqrt(mu)

    end function drc
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 2nd kind:
!  $$ R_D(x,y,z) = \frac{3}{2} \int_{0}^{\infty}
!                       (t+x)^{-1/2}
!                       (t+y)^{-1/2}
!                       (t+z)^{-3/2} dt $$
!  Where \(x\ge0\), \(y\ge0\), \(x+y>0\), and \(z>0\).
!
!  If \(x=0\) or \(y=0\), the integral is complete.
!
!  The duplication theorem is iterated until the variables are
!  nearly equal, and the function is then expanded in Taylor
!  series to fifth order.
!
!### DRD Special Comments
!
!  $$
!    \begin{array}{rl}
!      R_D(x,y,z) + R_D(y,z,x) + R_D(z,x,y) = \frac{3}{\sqrt{x y z}}  &  x>0, y>0, z>0
!    \end{array}
!  $$
!
!### Special functions via DRD and DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      E(\phi,k) = \sin \phi  R_F(\cos^2 \phi,1-k^2 \sin^2 \phi,1)
!      -\frac{k^2}{3} \sin^3 \phi R_D(\cos^2 \phi,1-k^2 \sin^2 \phi,1)
!    $$
!    When \( \phi = \pi /2 \) the integral is complete:
!    $$
!      \begin{array}{rcl}
!       E(k) &=& R_F(0,1-k^2 ,1) - \frac{k^2}{3} R_D(0,1-k^2 ,1) \\
!            &=& \int_{0}^{\pi/2} \sqrt{1-k^2 \sin^2 \phi}  d \phi
!      \end{array}
!    $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      \mathrm{EL2}(x,k_c,a,b) = ax R_F(1,1+k_c^2 x^2 ,1+x^2 )
!        + \frac{1}{3}(b-a) x^3 R_D(1,1+k_c^2 x^2 ,1+x^2 )
!    $$
!
!  * Legendre form of alternative ELLIPTIC INTEGRAL of 2nd kind:
!
!    $$
!      \begin{array}{rcl}
!      D(q,k) &=& \int_{0}^{q} \sin^2 p (1-k^2 \sin^2 p)^{-1/2} dp \\
!             &=& \frac{1}{3} (\sin^3 q) R_D(\cos^2 q,1-k^2 \sin^2 q,1)
!      \end{array}
!    $$
!
!  * Lemniscate constant B:
!
!    $$
!      \begin{array}{rcl}
!      B &=& \int_{0}^{1} s^2 (1-s^4)^{-1/2} ds \\
!        &=& \frac{1}{3} R_D (0,2,1)
!      \end{array}
!    $$
!
!  * Heuman's LAMBDA function:
!
!    $$
!      \begin{array}{rcl}
!      \frac{\pi}{2} \Lambda_0(a,b) &=&
!      \sin b \left(R_F(0,\cos^2 a,1)-\frac{1}{3} \sin^2 a
!      R_D(0,\cos^2 a,1) \right) R_F(\cos^2 b,1-\cos^2 a \sin^2 b,1) \\
!      & &-\frac{1}{3} \cos^2 a \sin^3 b R_F(0,\cos^2 a,1)
!      R_D(\cos^2 b,1-\cos^2 a \sin^2 b,1)
!      \end{array}
!    $$
!
!  * Jacobi ZETA function:
!
!    $$
!      \begin{array}{rcl}
!      Z(b,k) &=& \frac{k^2}{3} \sin b R_F(\cos^2 b,1-k^2 \sin^2 b,1)
!               R_D(0,1-k^2 ,1)/R_F(0,1-k^2 ,1) \\
!                & & -\frac{k^2}{3} \sin^3 b R_D(\cos^2 b,1-k^2 \sin^2 b,1)
!      \end{array}
!    $$
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Modify calls to XERMSG to put in standard form.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drd.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(dp) function drd(x,y,z,ier)

    implicit none

    real(dp), intent(in) :: x    !! nonnegative variable (\(x+y>0\))
    real(dp), intent(in) :: y    !! nonnegative variable (\(x+y>0\))
    real(dp), intent(in) :: z    !! positive variable
    integer, intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! _dpIER = 0_dp: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! _dpIER > 0_dp: Abnormal termination of the routine:
                                !! _dpIER = 1_dp: _dpmin(x,y) < 0_dp
                                !! _dpIER = 2_dp: _dpmin(x + y, z ) < LOLIM_dp
                                !! _dpIER = 3_dp: _dpmax(x,y,z) > UPLIM_dp

    character(len=16) :: xern3 , xern4 , xern5 , xern6
    real(dp) :: epslon, ea , eb , ec , ed , ef , lamda
    real(dp) :: mu , power4 , sigma , s1 , s2 , xn , xndev
    real(dp) :: xnroot , yn , yndev , ynroot , zn , zndev , znroot

    real(dp), parameter :: errtol = (d1mach(3)/3)**(1.0_dp/6)
        !! Determines the accuracy of the answer.
        !! The value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! machine precision.
        !!
        !! Relative error due to truncation is less than
        !! _dp3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2_dp.
        !!
        !! The accuracy of the computed approximation to the integral
        !! can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the truncation
        !! error there will be round-off error, but in practice the
        !! total error from both sources is usually less than the
        !! amount given in the table.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than):
        !! (1.0e-3, 4.0e-18),
        !! (3.0e-3, 3.0e-15),
        !! (1.0e-2, 4.0e-12),
        !! (3.0e-2, 3.0e-9),
        !! (1.0e-1, 4.0e-6)
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.

    real(dp), parameter :: lolim  = 2/(d1mach(2))**(2.0_dp/3) !! Lower limit of valid arguments
    real(dp), parameter :: tuplim = (errtol/10)**(1.0_dp/3)/&
                                    d1mach(1)**(1.0_dp/3)
    real(dp), parameter :: uplim  = tuplim**2  !! Upper limit of valid arguments
    real(dp), parameter :: c1     = 3.0_dp/14
    real(dp), parameter :: c2     = 1.0_dp/6
    real(dp), parameter :: c3     = 9.0_dp/22
    real(dp), parameter :: c4     = 3.0_dp/26

    ! initialize:
    drd = 0

    ! check for errors:
    if ( min(x,y)<0 ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write(error_unit,'(a)') 'drd: min(x,y)<0 where x = '//xern3// &
                                ' and y = '//xern4
        return
    endif

    if ( max(x,y,z)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') uplim
        write(error_unit,'(a)') 'drd: max(x,y,z)>uplim where x = '// &
                                xern3//' y = '//xern4//' z = '//xern5// &
                                ' and uplim = '//xern6
        return
    endif

    if ( min(x+y,z)<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') lolim
        write(error_unit,'(a)') 'drd: min(x+y,z)<lolim where x = '// &
                                xern3//' y = '//xern4//' z = '//xern5// &
                                ' and lolim = '//xern6
        return
    endif

    ier = 0; xn = x;  yn = y;  zn = z; sigma = 0;  power4 = 1

    do
        mu     = (xn + yn + 3 * zn)/5
        xndev  = (mu-xn)/mu
        yndev  = (mu-yn)/mu
        zndev  = (mu-zn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev))
        if ( epslon<errtol ) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
        sigma  = sigma + power4/(znroot*(zn+lamda))
        power4 = power4/4
        xn     = (xn+lamda)/4
        yn     = (yn+lamda)/4
        zn     = (zn+lamda)/4
    end do

    ea  = xndev * yndev
    eb  = zndev * zndev
    ec  = ea - eb
    ed  = ea - 6 * eb
    ef  = ed + ec + ec
    s1  = ed * ( c3 * ed/4 - c1 - 1.50_dp * c4 * zndev * ef)
    s2  = zndev * (c2 * ef + zndev * ( zndev * c4 * ea - c3 * ec))
    drd = 3 * sigma + power4*(1 + s1 + s2)/(mu*sqrt(mu))

    end function drd
!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 1st kind:
!  $$ R_F(x,y,z) = \frac{1}{2} \int_{0}^{\infty}
!                       (t+x)^{-1/2}
!                       (t+y)^{-1/2}
!                       (t+z)^{-1/2} dt $$
!  Where \(x\ge0\), \(y\ge0\), \(z\ge0\), and at most one of
!  them is \(=0\).
!
!  If \(x=0\), \(y=0\), or \(z=0\), the integral is complete.
!
!  The duplication theorem is iterated until the variables are
!  nearly equal, and the function is then expanded in Taylor
!  series to fifth order.
!
!### DRF Special Comments
!
!  $$
!    \begin{array}{rl}
!    R_F(x,x+z,x+w) + R_F(y,y+z,y+w) = R_F(0,z,w)
!    & x>0, y>0, z>0, x y = z w
!    \end{array}
!  $$
!
!### Special functions via DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 1st kind:
!
!    $$
!    \begin{array}{rl}
!      F(\phi,k) &= \sin \phi R_F( \cos^2 \phi,1-k^2 \sin^2 \phi,1) \\
!           K(k) &= R_F(0,1-k^2 ,1) = \int_{0}^{\pi/2} (1-k^2 sin^2 \phi )^{-1/2} d \phi
!    \end{array}
!    $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 1st kind:
!
!    $$
!      \mathrm{EL1}(x,k_c) = x R_F(1,1+k_c^2 x^2 ,1+x^2 )
!    $$
!
!  * Lemniscate constant A:
!
!    $$
!      A = \int_{0}^{1} (1-s^4 )^{-1/2}    ds = R_F(0,1,2) = R_F(0,2,1)
!    $$
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC))
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drf.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(dp) function drf(x,y,z,ier) ! NAG: S21BBF

    implicit none

    real(dp), intent(in) :: x    !! nonnegative variable
    real(dp), intent(in) :: y    !! nonnegative variable
    real(dp), intent(in) :: z    !! nonnegative variable
    integer, intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! _dpIER = 0_dp: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! _dpIER > 0_dp: Abnormal termination of the routine:
                                !! _dpIER = 1_dp: _dpmin(x,y,z) < 0_dp
                                !! _dpIER = 2_dp:_dp min(x+y,x+z,y+z) < LOLIM_dp
                                !! _dpIER = 3_dp: _dpmax(x,y,z) > UPLIM_dp

    character(len=16) :: xern3 , xern4 , xern5 , xern6
    real(dp) :: epslon, e2 , e3 , lamda, mu , s , xn , xndev, &
                xnroot , yn , yndev , ynroot , zn , zndev , znroot

    real(dp), parameter :: errtol = (4 * d1mach(3))**(1.0_dp/6)
        !! Determines the accuracy of the answer.
        !! The value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! machine precision.
        !!
        !! Relative error due to truncation is less than
        !! _dpERRTOL ** 6 / (4 * (1-ERRTOL)_dp.
        !!
        !! The accuracy of the computed approximation to the integral
        !! can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the truncation
        !! error there will be round-off error, but in practice the
        !! total error from both sources is usually less than the
        !! amount given in the table.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than):
        !! (1.0e-3, 3.0e-19),
        !! (3.0e-3, 2.0e-16),
        !! (1.0e-2, 3.0e-13),
        !! (3.0e-2, 2.0e-10),
        !! (1.0e-1, 3.0e-7)
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.

    real(dp), parameter :: lolim  = 5 * d1mach(1) !! Lower limit of valid arguments
    real(dp), parameter :: uplim  = d1mach(2)/5 !! Upper limit of valid arguments
    real(dp), parameter :: c1     = 1.0_dp/24
    real(dp), parameter :: c2     = 3.0_dp/44
    real(dp), parameter :: c3     = 1.0_dp/14

    ! initialize:
    drf = 0

    ! check for errors:
    if ( min(x,y,z) < 0 ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write(error_unit,'(a)') 'drf: min(x,y,z)<0 where x = '// &
                xern3//' y = '//xern4//' and z = '//xern5
        return
    endif

    if ( max(x,y,z)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') uplim
        write(error_unit,'(a)') 'drf: max(x,y,z)>uplim where x = '// &
                xern3//' y = '//xern4//' z = '//xern5// &
                ' and uplim = '//xern6
        return
    endif

    if ( min(x+y,x+z,y+z)<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') lolim
        write(error_unit,'(a)') 'drf: min(x+y,x+z,y+z)<lolim where x = '//xern3// &
                ' y = '//xern4//' z = '//xern5//' and lolim = '//xern6
        return
    endif

    ier = 0;  xn  = x;  yn  = y;  zn  = z

    do
        mu     = (xn+yn+zn)/3
        xndev  = 2 - (mu+xn)/mu
        yndev  = 2 - (mu+yn)/mu
        zndev  = 2 - (mu+zn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev))
        if ( epslon<errtol ) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
        xn     = (xn+lamda)/4
        yn     = (yn+lamda)/4
        zn     = (zn+lamda)/4
    end do
    e2  = xndev * yndev - zndev * zndev;  e3  = xndev * yndev * zndev
    s   = 1 + (c1*e2 - 0.10_dp - c2 * e3)*e2 + c3*e3;  drf = s/sqrt(mu)

    end function drf

!*******************************************************************************

!*******************************************************************************
!>
!  Compute an approximation for the incomplete or
!  complete elliptic integral of the 3rd kind:
!  $$
!    R_J(x,y,z,p) = \frac{3}{2} \int_{0}^{\infty}
!                   (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} (t+p)^{-1} dt
!  $$
!  where \(x\ge0\), \(y\ge0\), and \(z\ge0\),
!  and at most one of them \(=0\), and \(p>0\).
!
!  If \(x=0\) or \(y=0\) or \(z=0\), then the
!  integral is COMPLETE.
!
!  The duplication theorem is iterated
!  until the variables are nearly equal, and the function is
!  then expanded in Taylor series to fifth order.
!
!### DRJ Special Comments
!
!  $$
!  R_J(x,x+z,x+w,x+p) + R_J(y,y+z,y+w,y+p) +
!    (a-b) R_J(a,b,b,a) + \frac{3}{\sqrt{a}} = R_J(0,z,w,p)
!  $$
!  where:
!  $$
!      x>0, y>0, z>0, w>0, p>0
!  $$
!  and:
!  $$
!     \begin{array}{rl}
!         x y &= z w           \\
!           a &= p^2 (x+y+z+w) \\
!           b &= p (p+x) (p+y) \\
!       b - a &= p (p-z) (p-w)
!     \end{array}
!  $$
!  The sum of the third and
!  fourth terms on the left side is \( 3 R_C(a,b) \) .
!
!### Special functions via DRJ and DRF
!
!  * Legendre form of ELLIPTIC INTEGRAL of 3rd kind:
!
!  $$
!      \begin{array}{rcl}
!      P(\phi,k,n) &=& \int_{0}^{\phi} (1+n \sin^2 \theta )^{-1}
!                    (1-k^2 \sin^2 \theta )^{-1/2} d \theta \\
!                  &=& \sin \phi R_F(\cos^2 \phi, 1-k^2 sin^2 \phi,1)
!                    -\frac{n}{3} \sin^3 \phi R_J(\cos^2 \phi,1-k^2 \sin^2 \phi, 1,1+n \sin^2 \phi)
!      \end{array}
!  $$
!
!  * Bulirsch form of ELLIPTIC INTEGRAL of 3rd kind:
!
!  $$
!   \begin{array}{rcl}
!      \mathrm{EL3}(x,k_c,p) &=& x R_F(1,1+k_c^2 x^2 ,1+x^2 ) +
!      \frac{1}{3}(1-p) x^3 R_J(1,1+k_c^2 x^2 ,1+x^2 ,1+p x^2 ) \\
!      \mathrm{CEL}(k_c,p,a,b) &=& a R_F(0,k_c^2 ,1) +
!      \frac{1}{3}(b-pa) R_J(0,k_c^2 ,1,p)
!   \end{array}
!  $$
!
!  * Heuman's LAMBDA function:
!
!  $$
!   \begin{array}{rcl}
!   L(a,b,p) &=&  ( \cos^2 a \sin b \cos b /(1-\cos^2 a \sin^2 b )^{1/2} )
!          (\sin p R_F(\cos^2 p ,1-\sin^2 a \sin^2 p ,1) \\
!          & ~ & +
!          (\sin^2 a \sin^3 p /(3(1-\cos^2 a \sin^2 b )))
!          R_J(\cos^2 p ,1-\sin^2 a \sin^2 p ,1,1-
!          \sin^2 a \sin^2 p /(1-\cos^2 a \sin^2 b ))) \\
!   \frac{\pi}{2} \Lambda_0(a,b) &=& L(a,b,\pi/2) \\
!         &=& \cos^2 a  \sin b \cos b (1-\cos^2 a \sin^2 b )^{-1/2}
!         R_F(0,\cos^2 a ,1) \\
!         &~& + (1/3) \sin^2 a \cos^2 a
!         \sin b \cos b (1-\cos^2 a \sin^2 b )^{-3/2}
!         R_J(0,\cos^2 a ,1,\cos^2 a \cos^2 b /(1-\cos^2 a \sin^2 b ))
!   \end{array}
!  $$
!
!  * Jacobi ZETA function:
!
!  $$
!      Z(b,k) = \frac{k^2}{3} \sin b \cos b (1-k^2 \sin^2 b)^{1/2}
!               \frac{R_J(0,1-k^2 ,1,1-k^2 \sin^2 b)}{R_F (0,1-k^2 ,1)}
!  $$
!
!### Authors
!  * Carlson, B. C. Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Notis, E. M., Ames Laboratory-DOE, Iowa State University, Ames, IA  50011
!  * Pexton, R. L., Lawrence Livermore National Laboratory, Livermore, CA  94550
!
!### References
!  * B. C. Carlson and E. M. Notis, [Algorithms for incomplete
!    elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
!    ACM Transactions on Mathematical
!    Software 7, 3 (September 1981), pp. 398-403.
!  * B. C. Carlson, Computing elliptic integrals by
!    duplication, Numerische Mathematik 33, (1979),
!    pp. 1-16.
!  * B. C. Carlson, Elliptic integrals of the first kind,
!    SIAM Journal of Mathematical Analysis 8, (1977),
!    pp. 231-242.
!
!### History
!  * 790801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891009  Removed unreferenced statement labels.  (WRB)
!  * 891009  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900510  Changed calls to XERMSG to standard form, and some editorial changes.  (RWC)).
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jan 2016, Refactored [SLATEC routine](http://www.netlib.org/slatec/src/drj.f) into modern Fortran. (Jacob Williams)
!
!@warning Changes in the program may improve speed at the expense of robustness.

    real(dp) function drj(x,y,z,p,ier)! NAG : S21BDF

    implicit none

    real(dp), intent(in) :: x    !! nonnegative variable
    real(dp), intent(in) :: y    !! nonnegative variable
    real(dp), intent(in) :: z    !! nonnegative variable
    real(dp), intent(in) :: p    !! positive variable
    integer, intent(out) :: ier  !! indicates normal or abnormal termination:
                                !! _dpIER = 0_dp: Normal and reliable termination of the
                                !!  routine. It is assumed that the requested
                                !!  accuracy has been achieved.
                                !! _dpIER = 1_dp: _dpmin(x,y,z) < 0.0_dp_dp
                                !! _dpIER = 2_dp: _dpmin(x+y,x+z,y+z,p) < LOLIM_dp
                                !! _dpIER = 3_dp: _dpmax(x,y,z,p) > UPLIM_dp

    character(len=16) xern3 , xern4 , xern5 , xern6 , xern7
    real(dp) :: alfa, beta, ea, eb, ec, e2, e3, epslon, lamda, mu, pn, pndev, power4, &
    sigma, s1, s2, s3, xn, xndev, xnroot, yn, yndev, ynroot, zn, zndev, znroot

    real(dp), parameter :: errtol = (d1mach(3)/3)**(1.0_dp/6)
        !! Determines the accuracy of the answer
        !!
        !! the value assigned by the routine will result
        !! in solution precision within 1-2 decimals of
        !! "machine precision".
        !!
        !! Relative error due to truncation of the series for DRJ
        !! is less than _dp3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2_dp.
        !!
        !! The accuracy of the computed approximation to the integral
        !! can be controlled by choosing the value of ERRTOL.
        !! Truncation of a Taylor series after terms of fifth order
        !! introduces an error less than the amount shown in the
        !! second column of the following table for each value of
        !! ERRTOL in the first column.  In addition to the truncation
        !! error there will be round-off error, but in practice the
        !! total error from both sources is usually less than the
        !! amount given in the table.
        !!
        !! Sample choices:
        !! (ERRTOL, Relative truncation error less than}
        !! (1.0e-3, 4.0e-18),
        !! (3.0e-3, 3.0e-15),
        !! (1.0e-2, 4.0e-12),
        !! (3.0e-2, 3.0e-9),
        !! (1.0e-1, 4.0e-6)
        !!
        !! Decreasing ERRTOL by a factor of 10 yields six more
        !! decimal digits of accuracy at the expense of one or
        !! two more iterations of the duplication theorem.

    ! LOLIM and UPLIM determine the valid range of X, Y, Z, and P
    real(dp), parameter :: lolim  = (5 * d1mach(1) )**(1.0_dp/3)
        !! not less than the cube root of the value
        !! of LOLIM used in the routine for [[DRC]].
    real(dp), parameter :: uplim  = 0.3_dp * (d1mach(2)/5)**(1.0_dp/3)
        !! not greater than 0.3 times the cube root of
        !! the value of UPLIM used in the routine for [[DRC]].

    real(dp), parameter :: c1 = 3.0_dp/14
    real(dp), parameter :: c2 = 1.0_dp/3
    real(dp), parameter :: c3 = 3.0_dp/22
    real(dp), parameter :: c4 = 3.0_dp/26

    !initialize:
    drj = 0

    ! check for errors:
    if ( min(x,y,z) < 0 ) then
        ier = 1
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write(error_unit,'(a)') 'drj: min(x,y,z)<0 where x = '// &
            xern3//' y = '//xern4//' and z = '//xern5
        return
    endif

    if ( max(x,y,z,p)>uplim ) then
        ier = 3
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') p
        write (xern7,'(1pe15.6)') uplim
        write(error_unit,'(a)') 'drj: max(x,y,z,p)>uplim where x = '// &
            xern3//' y = '//xern4//' z = '//xern5//' p = '// &
            xern6//' and uplim = '//xern7
        return
    endif

    if ( min(x+y,x+z,y+z,p)<lolim ) then
        ier = 2
        write (xern3,'(1pe15.6)') x
        write (xern4,'(1pe15.6)') y
        write (xern5,'(1pe15.6)') z
        write (xern6,'(1pe15.6)') p
        write (xern7,'(1pe15.6)') lolim
        write(error_unit,'(a)') 'drj: min(x+y,x+z,y+z,p)<lolim where x = '//xern3// &
            ' y = '//xern4//' z = '//xern5//' p = '//xern6// &
            ' and lolim = '//xern7
        return
    endif

    ier = 0;  xn = x;  yn = y;  zn = z;  pn = p;  sigma  = 0;  power4 = 1

    do
        mu     = (xn+yn+zn+pn+pn)/5
        xndev  = (mu-xn)/mu
        yndev  = (mu-yn)/mu
        zndev  = (mu-zn)/mu
        pndev  = (mu-pn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev),abs(pndev))
        if ( epslon<errtol ) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda  = xnroot*(ynroot+znroot) + ynroot*znroot
        alfa   = pn*(xnroot+ynroot+znroot) + xnroot*ynroot*znroot
        alfa   = alfa*alfa
        beta   = pn*(pn+lamda)*(pn+lamda)
        sigma  = sigma + power4*drc(alfa,beta,ier)
        power4 = power4/4
        xn     = (xn+lamda)/4
        yn     = (yn+lamda)/4
        zn     = (zn+lamda)/4
        pn     = (pn+lamda)/4
    end do

    ea  = xndev * (yndev + zndev) + yndev * zndev
    eb  = xndev * yndev * zndev
    ec  = pndev * pndev
    e2  = ea - 3 * ec
    e3  = eb + 2 * pndev*(ea-ec)
    s1  = 1 + e2 * (0.750_dp * c3 * e2 - c1 - 1.50_dp * c4 * e3)
    s2  = eb * ( c2/2 + pndev*(pndev * c4 - 2 * c3 ) )
    s3  = pndev * ea * (c2 - pndev*c3) - c2*pndev*ec
    drj = 3 * sigma + power4 * (s1 + s2 + s3)/( mu * sqrt(mu) )

    end function drj
!*******************************************************************************

end module carlson_elliptic_module

!ccccccccccccccc

module Chaplin
  use constants, only: dp, zeta2, zeta3, Pi2
  implicit none

  real (dp), private, parameter :: border = 0.3_dp, eps = 1e-14_dp

  private :: bsli2_inside, bsli2_outside, bsli3_inside

  contains

!ccccccccccccccc

  recursive function cli2(z) result(ris)
    complex (dp), intent(in) :: z

    complex (dp) :: ris, zlocal
    real(dp) :: zabs, arg

    zabs = abs(z);  zlocal = z

    if (zabs > 1 + eps) then
      ris = - cli2(1/z) - zeta2 - log(-z)**2/2
    else if (zabs <= border) then
      ris = bsli2_inside(z)
    else
      if (zabs > 1) then
        arg = atan2( dimag(zlocal), dreal(zlocal) )
        zlocal = dcmplx( cos(arg), sin(arg) )
      endif
      ris = bsli2_outside(zlocal)
    endif

  end function cli2

!ccccccccccccccc

  recursive function cli3(z) result(ris)
    complex (dp), intent(in) :: z

    complex (dp) :: ris, zlocal
    real(dp) :: zabs, arg

    zabs = abs(z);  zlocal = z

    if (zabs > 1 + eps) then
      ris = cli3(1/z) - log(-z)**3/6 - zeta2 * log(-z)
    else if (zabs <= border) then
      ris = bsli3_inside(z)
    else
      if (zabs > 1) then
        arg = atan2(dimag( zlocal), dreal(zlocal) ); zlocal = dcmplx( cos(arg), sin(arg) )
      endif
      ris = bsli3_outside(zlocal)
    endif

  end function cli3

!ccccccccccccccc

  complex (dp) function bsli2_inside(z)
    complex (dp), intent(in) :: z

    integer, parameter :: Nmax = 11 ! this is half the order we want (coz odd bernoulli numbers are zero except BernoulliB(1)=-0.5_dp)
    integer :: i
    complex (dp) :: ris, zb
    real (dp), dimension(11), parameter :: bern =  [1._dp,0.8333333333333333e-1_dp,      &
      -0.1388888888888889e-2_dp, 0.3306878306878307e-4_dp, -0.8267195767195767e-6_dp,    &
       0.208767569878681e-7_dp , -0.5284190138687493e-9_dp, 0.1338253653068468e-10_dp,   &
      -0.3389680296322583e-12_dp, 0.8586062056277845e-14_dp, -0.2174868698558062e-15_dp ]

    zb = dcmplx(1_dp,0_dp) - z;  zb = - log(zb)
    ris = - zb**2/4          !accounting for BernoulliB(1) = -0.5_dp

    do i = 1, Nmax
      ris = ris + zb**(2*i-1) * bern(i)/(2 * i - 1)
    enddo

    bsli2_inside = ris

  end function bsli2_inside

!ccccccccccccccc

  complex (dp) function bsli2_outside(z)
    complex (dp), intent(in) :: z

    integer, parameter :: Nmax = 29 ! this is half the order we want (coz odd bernoulli numbers are zero except BernoulliB(1)=-0.5_dp)
    integer            :: i
    complex (dp) :: ris, zb

    real (dp), dimension(29), parameter :: zeta = [                                       &
    -0.01388888888888889_dp, 0.00006944444444444444_dp, -7.873519778281683e-7_dp,         &
     1.148221634332745e-8_dp,  -1.897886998897100e-10_dp, 3.387301370953521e-12_dp,       &
     -6.372636443183180e-14_dp, 1.246205991295067e-15_dp, -2.510544460899955e-17_dp,      &
     5.178258806090624e-19_dp, -1.088735736830085e-20_dp, 2.325744114302087e-22_dp,       &
     -5.035195213147390e-24_dp, 1.102649929438122e-25_dp, -2.438658550900734e-27_dp,      &
     5.440142678856252e-29_dp, -1.222834013121735e-30_dp, 2.767263468967951e-32_dp,       &
     -6.300090591832014e-34_dp, 1.442086838841848e-35_dp, -3.317093999159543e-37_dp,      &
     7.663913557920658e-39_dp, -1.777871473383066e-40_dp, 4.139605898234137e-42_dp,       &
     -9.671557036081102e-44_dp, 2.266718701676613e-45_dp, -5.327956311328254e-47_dp,      &
     1.255724838956433e-48_dp, -2.967000542247094e-50_dp]

    zb = log(z); ris = dcmplx(zeta2, 0_dp) + zb * ( 1 - log(-zb) ) + zb**2 * zeta2

    do i = 1, Nmax
      ris = ris + zb**(2 * i + 1) * zeta(i)
    enddo

    bsli2_outside=ris

  end function bsli2_outside

!ccccccccccccccc

  complex (dp) function bsli3_inside(z)
   complex (dp), intent(in) :: z

   integer, parameter :: Nmax = 21
   integer            :: i
   complex (dp) :: ris, zb

   real (dp), parameter, dimension(21) :: bern = [ 1._dp, -0.75_dp,                   &
   0.2361111111111111_dp, -0.03472222222222222_dp, 0.0006481481481481481_dp,          &
   0.0004861111111111111_dp, -0.00002393550012597632_dp, -0.00001062925170068027_dp,  &
   7.794784580498866e-7_dp, 2.526087595532040e-7_dp, -2.359163915200471e-8_dp,        &
   -6.168132746415575e-9_dp, 6.8244567489810783e-10_dp, 1.524285616929085e-10_dp,     &
   -1.916909414174054e-11_dp, -3.791718683693992e-12_dp, 5.277408409541286e-13_dp,    &
   9.471165533842511e-13_dp, -1.432311114490360e-14_dp, -2.372464515550457e-15_dp,    &
   3.846565792753191e-16_dp]

   zb = dcmplx(1_dp,0_dp) - z;  zb = - log(zb);  ris = dcmplx(0_dp, 0_dp)

   do i = 1, Nmax
     ris = ris + zb**i * bern(i)/i
   enddo

   bsli3_inside=ris

  end function bsli3_inside

!ccccccccccccccc

 complex (dp) function bsli3_outside(z)
   complex (dp), intent(in) :: z

   integer, parameter :: Nmax = 30
   integer            :: i
   complex (dp) :: ris, zb

   real (dp), parameter, dimension(29) :: zeta = &
   [-0.003472222222222222_dp, 0.00001157407407407407_dp, -9.841899722852104e-8_dp,  &
     1.148221634332745e-9_dp, -1.581572499080917e-11_dp,  2.419500979252515e-13_dp, &
   -3.982897776989488e-15_dp,  6.923366618305929e-17_dp, -1.255272230449977e-18_dp, &
    2.353754002768465e-20_dp, -4.536398903458687e-22_dp,  8.945169670392643e-24_dp, &
   -1.798284004695496e-25_dp,  3.675499764793738e-27_dp, -7.620807971564795e-29_dp, &
    1.600041964369486e-30_dp, -3.396761147560376e-32_dp,  7.282272286757765e-34_dp, &
   -1.575022647958003e-35_dp,  3.433540092480589e-37_dp, -7.538849998089870e-39_dp, &
    1.666068164765360e-40_dp, -3.703898902881387e-42_dp,  8.279211796468275e-44_dp, &
   -1.859914814630981e-45_dp,  4.197627225327060e-47_dp, -9.514207698800454e-49_dp, &
    2.165042825786954e-50_dp, -4.945000903745158e-52_dp]

   zb = log(z)
   ris = dcmplx(zeta2, 0_dp) + zb * zeta2 + zb**3 * zeta3 + zb**2 * ( 1.5_dp -log(-zb) )/2

   do i = 2, Nmax
     ris = ris + zb**(2*i) * zeta(i - 1)
   enddo

   bsli3_outside=ris

  end function bsli3_outside

end module Chaplin

!ccccccccccccccc

module DeriGamma
use constants
Private :: I1MACH, J4SAVE

contains

recursive subroutine DPSifN(X, N, KODE, M, ANS, NZ, IERR)
!
! Modified by Vicent Mateu on 26/01/2017; made it recursive, now accepts negative non-integer arguments
!
!! DPSifN computes derivatives of the Psi function.
!
!***LIBRARY   SLATEC
!***CATEGORY  C7C
!***TYPE      DOUBLE PRECISION (PSifN-S, DPSifN-D)
!***KEYWORDS  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION,
!             PSI FUNCTION
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!         The following definitions are used in DPSifN:
!
!      Definition 1
!         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of
!                  the log GAMMA function.
!      Definition 2
!                     K   K
!         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X).
!   ___________________________________________________________________
!      DPSifN computes a sequence of SCALED derivatives of
!      the PSI function; i.e. for fixed X and M it computes
!      the M-member sequence
!
!                    ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X)
!                       for K = N,...,N+M-1
!
!      where PSI(K,X) is as defined above.   For KODE=1, DPSifN returns
!      the scaled derivatives as described.  KODE=2 is operative only
!      when K=0 and in that case DPSifN returns -PSI(X) + LN(X).  That
!      is, the logarithmic behavior for large X is removed when KODE=2
!      and K=0.  When sums or differences of PSI functions are computed
!      the logarithmic terms can be combined analytically and computed
!      separately to help retain significant digits.
!
!         Note that call DPSifN(X,0,1,1,ANS) results in
!                   ANS = -PSI(X)
!
!     Input      X is DOUBLE PRECISION
!           X      - Argument, X  >  0._dp
!           N      - First member of the sequence, 0  <=  N  <=  100
!                    N=0 gives ANS(1) = -PSI(X)       for KODE=1
!                                       -PSI(X)+LN(X) for KODE=2
!           KODE   - Selection parameter
!                    KODE=1 returns scaled derivatives of the PSI
!                    function.
!                    KODE=2 returns scaled derivatives of the PSI
!                    function EXCEPT when N=0. In this case,
!                    ANS(1) = -PSI(X) + LN(X) is returned.
!           M      - Number of members of the sequence, M >= 1
!
!    Output     ANS is DOUBLE PRECISION
!           ANS    - A vector of length at least M whose first M
!                    components contain the sequence of derivatives
!                    scaled according to KODE.
!           NZ     - Underflow flag
!                    NZ == 0, A normal return
!                    NZ.ne.0, Underflow, last NZ components of ANS are
!                             set to zero, ANS(M-K+1)=0.0, K=1,...,NZ
!           IERR   - Error flag
!                    IERR=0, A normal return, computation completed
!                    IERR=1, Input error,     no computation
!                    IERR=2, Overflow,        X too small or N+M-1 too
!                            large or both
!                    IERR=3, Error,           N too large. Dimensioned
!                            array TRMR(NMAX) is not large enough for N
!
!         The nominal computational accuracy is the maximum of unit
!         roundoff (=D1MACH(4)) and 1.d-18 since critical constants
!         are given to only 18 digits.
!
!         PSifN is the single precision version of DPSifN.
!
! *Long Description:
!
!         The basic method of evaluation is the asymptotic expansion
!         for large X >= XMIN followed by backward recursion on a two
!         term recursion relation
!
!                  W(X+1) + X**(-N-1) = W(X).
!
!         This is supplemented by a series
!
!                  SUM( (X+K)**(-N-1) , K=0,1,2,... )
!
!         which converges rapidly for large N. Both XMIN and the
!         number of terms of the series are calculated from the unit
!         roundoff of the machine environment.
!
!***REFERENCES  Handbook of Mathematical Functions, National Bureau
!                 of Standards Applied Mathematics Series 55, edited
!                 by M. Abramowitz and I. A. Stegun, equations 6.3.5,
!                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
!               D. E. Amos, A portable Fortran subroutine for
!                 derivatives of the Psi function, Algorithm 610, ACM
!                 Transactions on Mathematical Software 9, 4 (1983),
!                 pp. 494-502.
!***ROUTINES CALLED  D1MACH, I1MACH
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPSifN
  implicit none
  INTEGER :: I, IERR, J, K, KODE, M, MM, MX, N, NMAX, NN, NP, NX, NZ, FN
  real(dp) :: ANS, ARG, B, DEN, ELIM, EPS, FLN, XM, XMIN, XQ, YINT, &
   FX, RLN, RXSQ, R1M4, R1M5, S, SLOPE, T, TA, TK, TOL, TOLS, TRM, &
   TRMR, TSS, TST, TT, T1, T2, WDTOL, X, XDMLN, XDMY, XINC, XLN

  DIMENSION B(22), TRM(22), TRMR(100), ANS(*)
  SAVE NMAX, B
  DATA NMAX /100/
!-----------------------------------------------------------------------
!             BERNOULLI NUMBERS
!-----------------------------------------------------------------------
  DATA B  /1_dp, -5e-1_dp,1.66666666666666667e-1_dp, &
   -3.33333333333333333e-2_dp, 2.38095238095238095e-2_dp, &
   -3.33333333333333333e-2_dp, 7.57575757575757576e-2_dp, &
   -2.53113553113553114e-1_dp, 1.16666666666666667_dp,    &
   -7.09215686274509804_dp,    5.49711779448621554e1_dp,  &
   -5.29124242424242424e2_dp,  6.19212318840579710e3_dp,  &
   -8.65802531135531136e4_dp,  1.42551716666666667e6_dp,  &
   -2.72982310678160920e7_dp,  6.01580873900642368e8_dp,  &
   -1.51163157670921569e10_dp, 4.29614643061166667e11_dp, &
   -1.37116552050883328e13_dp, 4.88332318973593167e14_dp, &
   -1.92965793419400681e16_dp/

   if (x < 0) then

     call DPSifN(1 + X, N, KODE, M, ANS(:M), NZ, IERR)
     ANS(:M) = ANS(:M) + [ ( 1/x**(i+1), i = N, N + M - 1 ) ]; return

   end if
!
!***FIRST EXECUTABLE STATEMENT  DPSifN
  IERR = 0
  NZ=0
  if (X <= 0._dp) IERR=1
  if (N < 0) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (M < 1) IERR=1
  if (IERR /= 0) RETURN
  MM=M
  NX = MIN(-I1MACH(15),I1MACH(16))
  R1M5 = D1MACH(5)
  R1M4 = D1MACH(4)*0.5_dp
  WDTOL = MAX(R1M4,0.5e-18_dp)
!-----------------------------------------------------------------------
!     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  ELIM = 2.302_dp*(NX*R1M5-3._dp)
  XLN = LOG(X)
   41 CONTINUE
  NN = N + MM - 1
  FN = NN
  T = (FN+1)*XLN
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X
!-----------------------------------------------------------------------
  if (ABS(T) > ELIM) go to 290
  if (X < WDTOL) go to 260
!-----------------------------------------------------------------------
!     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1
!-----------------------------------------------------------------------
  RLN = R1M5*I1MACH(14)
  RLN = MIN(RLN,18.06_dp)
  FLN = MAX(RLN,3._dp) - 3._dp
  YINT = 3.50_dp + 0.40_dp*FLN
  SLOPE = 0.21_dp + FLN*(0.0006038_dp*FLN+0.008677_dp)
  XM = YINT + SLOPE*FN
  MX = INT(XM) + 1
  XMIN = MX
  if (N == 0) go to 50
  XM = -2.302_dp*RLN - MIN(0._dp,XLN)
  ARG = XM/N
  ARG = MIN(0._dp,ARG)
  EPS = EXP(ARG)
  XM = 1._dp - EPS
  if (ABS(ARG) < 1e-3_dp) XM = -ARG
  FLN = X*XM/EPS
  XM = XMIN - X
  if (XM > 7._dp .AND. FLN < 15._dp) go to 200
   50 CONTINUE
  XDMY = X
  XDMLN = XLN
  XINC = 0._dp
  if (X >= XMIN) go to 60
  NX = INT(X)
  XINC = XMIN - NX
  XDMY = X + XINC
  XDMLN = LOG(XDMY)
   60 CONTINUE
!-----------------------------------------------------------------------
!     GENERATE W(N+MM-1,X) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
  T = FN*XDMLN
  T1 = XDMLN + XDMLN
  T2 = T + XDMLN
  TK = MAX(ABS(T),ABS(T1),ABS(T2))
  if (TK > ELIM) go to 380
  TSS = EXP(-T)
  TT = 0.5_dp/XDMY
  T1 = TT
  TST = WDTOL*TT
  if (NN /= 0) T1 = TT + 1._dp/FN
  RXSQ = 1._dp/(XDMY*XDMY)
  TA = 0.5_dp*RXSQ
  T = (FN+1)*TA
  S = T*B(3)
  if (ABS(S) < TST) go to 80
  TK = 2._dp
  DO 70 K=4,22
    T = T*((TK+FN+1)/(TK+1._dp))*((TK+FN)/(TK+2._dp))*RXSQ
    TRM(K) = T*B(K)
    if (ABS(TRM(K)) < TST) go to 80
    S = S + TRM(K)
    TK = TK + 2._dp
   70 CONTINUE
   80 CONTINUE
  S = (S+T1)*TSS
  if (XINC == 0._dp) go to 100
!-----------------------------------------------------------------------
!     BACKWARD RECUR FROM XDMY TO X
!-----------------------------------------------------------------------
  NX = INT(XINC)
  NP = NN + 1
  if (NX > NMAX) go to 390
  if (NN == 0) go to 160
  XM = XINC - 1._dp
  FX = X + XM
!-----------------------------------------------------------------------
!     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL
!-----------------------------------------------------------------------
  DO 90 I=1,NX
    TRMR(I) = FX**(-NP)
    S = S + TRMR(I)
    XM = XM - 1._dp
    FX = X + XM
   90 CONTINUE
  100 CONTINUE
  ANS(MM) = S
  if (FN == 0) go to 180
!-----------------------------------------------------------------------
!     GENERATE LOWER DERIVATIVES, J < N+MM-1
!-----------------------------------------------------------------------
  if (MM == 1) RETURN
  DO 150 J=2,MM
    FN = FN - 1
    TSS = TSS*XDMY
    T1 = TT
    if (FN /= 0) T1 = TT + 1._dp/FN
    T = (FN+1)*TA
    S = T*B(3)
    if (ABS(S) < TST) go to 120
    TK = 4 + FN
    DO 110 K=4,22
      TRM(K) = TRM(K)*(FN+1)/TK
      if (ABS(TRM(K)) < TST) go to 120
      S = S + TRM(K)
      TK = TK + 2._dp
  110   CONTINUE
  120   CONTINUE
    S = (S+T1)*TSS
    if (XINC == 0._dp) go to 140
    if (FN == 0) go to 160
    XM = XINC - 1._dp
    FX = X + XM
    DO 130 I=1,NX
      TRMR(I) = TRMR(I)*FX
      S = S + TRMR(I)
      XM = XM - 1._dp
      FX = X + XM
  130   CONTINUE
  140   CONTINUE
    MX = MM - J + 1
    ANS(MX) = S
    if (FN == 0) go to 180
  150 CONTINUE
  return
!-----------------------------------------------------------------------
!     RECURSION FOR N = 0
!-----------------------------------------------------------------------
  160 CONTINUE
  DO 170 I=1,NX
    S = S + 1._dp/(X+NX-I)
  170 CONTINUE
  180 CONTINUE
  if (KODE == 2) go to 190
  ANS(1) = S - XDMLN
  return
  190 CONTINUE
  if (XDMY == X) RETURN
  XQ = XDMY/X
  ANS(1) = S - LOG(XQ)
  return
!-----------------------------------------------------------------------
!     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,...
!-----------------------------------------------------------------------
  200 CONTINUE
  NN = INT(FLN) + 1
  NP = N + 1
  T1 = (N+1)*XLN
  T = EXP(-T1)
  S = T
  DEN = X
  DO 210 I=1,NN
    DEN = DEN + 1._dp
    TRM(I) = DEN**(-NP)
    S = S + TRM(I)
  210 CONTINUE
  ANS(1) = S
  if (N /= 0) go to 220
  if (KODE == 2) ANS(1) = S + XLN
  220 CONTINUE
  if (MM == 1) RETURN
!-----------------------------------------------------------------------
!     GENERATE HIGHER DERIVATIVES, J > N
!-----------------------------------------------------------------------
  TOL = WDTOL/5._dp
  DO 250 J=2,MM
    T = T/X
    S = T
    TOLS = T*TOL
    DEN = X
    DO 230 I=1,NN
      DEN = DEN + 1._dp
      TRM(I) = TRM(I)/DEN
      S = S + TRM(I)
      if (TRM(I) < TOLS) go to 240
  230   CONTINUE
  240   CONTINUE
    ANS(J) = S
  250 CONTINUE
  return
!-----------------------------------------------------------------------
!     SMALL X < UNIT ROUND OFF
!-----------------------------------------------------------------------
  260 CONTINUE
  ANS(1) = X**(-N-1)
  if (MM == 1) go to 280
  K = 1
  DO 270 I=2,MM
    ANS(K+1) = ANS(K)/X
    K = K + 1
  270 CONTINUE
  280 CONTINUE
  if (N /= 0) RETURN
  if (KODE == 2) ANS(1) = ANS(1) + XLN
  return
  290 CONTINUE
  if (T > 0._dp) go to 380
  NZ=0
  IERR=2
  return
  380 CONTINUE
  NZ=NZ+1
  ANS(MM)=0._dp
  MM=MM-1
  if (MM == 0) RETURN
  go to 41
  390 CONTINUE
  NZ=0
  IERR=3
end

!ccccccccccccccc

FUNCTION I1MACH (I)
!
!! I1MACH returns integer machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0  <=  X(I)  <  A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0  <=  X(I)  <  B for I=1,...,T,
!                0  <  X(1), and EMIN  <=  E  <=  EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   891012  Added VAX G-floating constants.  (WRB)
!   891012  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
!           (RWC)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added Convex -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
!           options.  (DWL, RWC and WRB).
!***END PROLOGUE  I1MACH
!
  integer i1mach
  INTEGER IMACH(16),OUTPUT
  SAVE IMACH
  EQUIVALENCE (IMACH(4),OUTPUT)
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        129 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1025 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA IMACH( 1) /          7 /
!     DATA IMACH( 2) /          2 /
!     DATA IMACH( 3) /          2 /
!     DATA IMACH( 4) /          2 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         33 /
!     DATA IMACH( 9) / Z1FFFFFFFF /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -256 /
!     DATA IMACH(13) /        255 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /       -256 /
!     DATA IMACH(16) /        255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /        -50 /
!     DATA IMACH(16) /         76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /     -32754 /
!     DATA IMACH(16) /      32780 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -4095 /
!     DATA IMACH(13) /       4094 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -4095 /
!     DATA IMACH(16) /       4094 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /    6LOUTPUT/
!     DATA IMACH( 5) /         60 /
!     DATA IMACH( 6) /         10 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         48 /
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /       -929 /
!     DATA IMACH(13) /       1070 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /       -929 /
!     DATA IMACH(16) /       1069 /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / Z'7FFFFFFF' /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16383 /
!     DATA IMACH(16) /      16383 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -pd8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 46 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         46 /
!     DATA IMACH( 9) / 1777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 64 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 777777777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     DATA IMACH( 1) /         11 /
!     DATA IMACH( 2) /         12 /
!     DATA IMACH( 3) /          8 /
!     DATA IMACH( 4) /         10 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          6 /
      DATA IMACH( 4) /          6 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        128 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         24 /
!     DATA IMACH( 6) /          3 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         23 /
!     DATA IMACH( 9) /    8388607 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         38 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /         43 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         63 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         39 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         55 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          7 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1015 /
!     DATA IMACH(16) /       1017 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) /  Z7FFFFFFF /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         54 /
!     DATA IMACH(15) /       -101 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         62 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1021 /
!     DATA IMACH(13) /       1024 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16381 /
!     DATA IMACH(16) /      16384 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          1 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /      -1024 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA IMACH( 1) /          1 /
!     DATA IMACH( 2) /          1 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!***FIRST EXECUTABLE STATEMENT  I1MACH
!
  if ( I < 1 .OR. I > 16 ) then
    WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
    STOP
  end if

  I1MACH = IMACH(I)

end

!ccccccccccccccc

function J4SAVE (IWHICH, IVALUE, ISET)
!
!! J4SAVE saves or recalls global variables needed by error handling routines.
!
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - if ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  if ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
  LOGICAL ISET
  INTEGER IPARAM(9)
  SAVE IPARAM
  DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
  DATA IPARAM(5)/1/
  DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
  J4SAVE = IPARAM(IWHICH)
  if (ISET) IPARAM(IWHICH) = IVALUE
end

!ccccccccccccccc

end module DeriGamma

!ccccccccccccccc

module Poly
  use constants;  implicit none

  private :: PolyLogSum, PolySum, BernoulliPoly, BernoulliPolyComp, Harmonic, Zeta, &
           factorial, bernoa, bernoulli, Pi, ComplexPolyLogSum, ComplexPolySum

 contains

 recursive function PolyLog(n, z) result (pl)
   integer, intent(in) :: n
   real (dp), intent(in) :: z

   real (dp) :: pl

   if (n < 0) then

     pl = 0

   else if (n == 0) then

     pl = z/(1 - z)

   else if (n == 1) then

     pl = - log(1 - z)

   else if ( abs(z - 1) < 1e-8_dp) then

     pl = Zeta(n)

   else if ( abs(z + 1) < 1e-8_dp) then

     pl = ( 2._dp**(1 - n) - 1 ) * Zeta(n)

   else if ( abs(z) <= 0.5 ) then

     pl = PolySum(n,z)

   else if ( abs(z) > 1 ) then

     pl = BernoulliPoly(n,z) - (-1)**n * PolyLog(n, 1/z)

   else if ( z < 0 ) then

     pl = 2 * PolyLog(n,z**2)/2**n - PolyLog(n, - z)

   else

     pl = PolyLogSum( n, log(z) )

   end if

 end function PolyLog

!ccccccccccccccc

 recursive function ComplexPolyLog(n, z) result (pl)
   implicit none
   integer, intent(in) :: n
   complex (dp), intent(in) :: z

   complex (dp) :: pl

   if ( abs( aimag(z) ) < 1e-8_dp .and. real(z) <= 1 ) then

      pl = (1._dp, 0._dp) * PolyLog( n, real(z) )

   else if (n < 0) then

     pl = 0

   else if (n == 0) then

     pl = z/(1 - z)

   else if (n == 1) then

     pl = - log(1 - z)

   else if ( abs(z) <= 0.5 ) then

     pl = ComplexPolySum(n,z)

   else if ( abs(z) > 1 ) then

     pl = ComplexBernoulliPoly(n,z) - (-1)**n * ComplexPolyLog(n, 1._dp/z)

   else

     pl = ComplexPolyLogSum( n, log(z) )

   end if

 end function ComplexPolyLog

!ccccccccccccccc

 real(dp) function PolyLogSum(n, lz)
   implicit none
   integer, intent(in) :: n
   real(dp), intent(in) :: lz

   real(dp) :: ff, lzpow, test
   integer :: i

   PolyLogSum = lz**(n-1)/factorial(n - 1) * ( harmonic(n - 1) - log(-lz) ) + &
   Zeta(n)

   ff = 1;   lzpow = 1;  i = 0

   enless_loop: do

   i = i + 1;   ff = i * ff;  lzpow = lz * lzpow

   if ( n - i == 1) cycle enless_loop

   test = Zeta(n - i) * lzpow/ff

   PolyLogSum = PolyLogSum + test
   if ( i > 20 .and. abs(test) < 1e-10_dp ) exit enless_loop

   end do enless_loop

end function PolyLogSum

!ccccccccccccccc

 complex (dp) function ComplexPolyLogSum(n, lz)
   implicit none
   integer, intent(in) :: n
   complex(dp), intent(in) :: lz

   complex(dp) :: ff, lzpow, test
   integer :: i

   ComplexPolyLogSum = lz**(n-1)/factorial(n - 1) * ( harmonic(n - 1) - log(-lz) ) + &
   Zeta(n)

   ff = 1;   lzpow = 1;  i = 0

   enless_loop: do

     i = i + 1;  ff = i * ff;  lzpow = lz * lzpow

     if ( n - i == 1 ) cycle enless_loop

     test = Zeta(n - i) * lzpow/ff

     ComplexPolyLogSum = ComplexPolyLogSum + test
     if ( i > 20 .and. 100 * max(  ( aimag(test) ), abs( real(test) )  ) < 1e-10_dp ) exit enless_loop

   end do enless_loop

end function ComplexPolyLogSum

!ccccccccccccccc

 real(dp) function PolySum(n,z)
   implicit none
   integer, intent(in) :: n
   real(dp), intent(in) :: z

   real(dp) :: test, zpow
   integer :: i

   if (n > 0) then

     zpow = z;  i = 1;  PolySum = z

     enless_loop: do

       zpow = zpow * z;  i = i + 1;  test = zpow/i**n;  PolySum = PolySum + test
       if (Abs(test) < 1e-10_dp) exit enless_loop

     end do enless_loop
   end if

 end function PolySum

!ccccccccccccccc

 complex (dp) function ComplexPolySum(n,z)
   implicit none
   integer          , intent(in) :: n
   complex(dp), intent(in) :: z

   complex (dp)            :: test, zpow
   integer :: i

   if (n > 0) then

     zpow = z;  i = 1;  ComplexPolySum = z

     enless_loop: do

       zpow = zpow * z;  i = i + 1;  test = zpow/i**n
       ComplexPolySum = ComplexPolySum + test
       if (max(  ( aimag(test) ), abs( real(test) )  ) < 1e-10_dp) exit enless_loop

     end do enless_loop
   end if

 end function ComplexPolySum

!ccccccccccccccc

complex (dp) function BernoulliPolyComp(n,z)
  implicit none
  integer           , intent(in) :: n
  complex (dp), intent(in) :: z

  real (dp) :: bn(0:n)
  integer :: i

  call bernoulli(n, bn)

  BernoulliPolyComp = 0

  do i = 0, n

    BernoulliPolyComp = BernoulliPolyComp + binom(n,i) * bn(n - i) * z**i

  end do

end function BernoulliPolyComp

!ccccccccccccccc

  real (dp) function BernoulliPoly(n,z)
    implicit none
    integer, intent(in) :: n
    real (dp), intent(in) :: z

    complex (dp) :: CI, Cphase, lz

    CI = (0._dp, 1._dp);  Cphase = 2 * Pi * CI;  lz = log( z * (1._dp,0._dp) )/Cphase

    BernoulliPoly = - Real( (2*Pi*CI)**n * BernoulliPolyComp(n,lz) )/factorial(n)

  end function BernoulliPoly

!ccccccccccccccc

  complex (dp) function ComplexBernoulliPoly(n,z)
    implicit none
    integer, intent(in) :: n
    complex (dp), intent(in) :: z

    complex (dp) :: CI, Cphase, lz

    CI = (0._dp, 1._dp);  Cphase = 2 * Pi * CI;  lz = log(z)

    ComplexBernoulliPoly = - (2*Pi*CI)**n * BernoulliPolyComp(n,lz/Cphase)/factorial(n)

    if ( aimag(z) < 0 .or. ( abs( aimag(z) ) < 1e-8_dp .and. real(z) > 1 ) ) then

      ComplexBernoulliPoly = ComplexBernoulliPoly - Cphase * lz**(n - 1)/factorial(n - 1)

    end if

  end function ComplexBernoulliPoly

!ccccccccccccccc

  recursive function harmonic(i) result(prod)
    implicit none
    integer, intent(in) :: i
   real (dp) :: prod

    select case(i)
      case(:-1)
       prod = 0._dp
    case(1)
       prod = 1._dp
    case default
       prod = harmonic(i - 1) + 1._dp/i
    end select
end function Harmonic

!ccccccccccccccc

  real(dp) function Zeta(n)
   implicit none
   integer, intent(in) :: n

   integer             :: i
   real(dp)      :: test

   if (n < 0) then

     Zeta = - bernoa(-n + 1)/(-n + 1)

   else if (n == 0) then

     Zeta = - 0.5

   else if (n == 1) then

     Zeta = 1.d100

   else if ( mod(n,2) == 0 ) then

     Zeta = (-1)**(n/2 + 1) * bernoa(n) * (2*Pi)**n/2/factorial(n)

   else

     Zeta = 1;   i = 1

     enless_loop: do

       i = i + 1

       test = 1/(1._dp*i)**n;  Zeta = Zeta + test
       if (Abs(1000 * test) < 1e-10_dp) exit enless_loop

     end do enless_loop
   end if

  end function Zeta

!ccccccccccccccc

  recursive function factorial(i) result(prod)
    integer, intent(in) :: i
    real (dp) :: prod

    select case(i)
    case(:-1)
      prod = 0._dp
    case(0:1)
      prod = 1._dp
    case default
      prod = i * factorial(i - 1)
    end select

  end function factorial

!ccccccccccccccc

  real (dp) function binom(i,j)
    integer, intent(in) :: i, j

    integer :: k

    if (j < 0) then

      binom = 0._dp

    else if (j == 0) then

      binom = 1._dp

    else

      binom = i

      do k = 1, j - 1

        binom = binom * (i - k)

      end do

      binom = binom / factorial(j)

    end if
  end function binom

!ccccccccccccccc

  real (dp) function bernoa(n)
    implicit none
    integer, intent(in) :: n

    real (dp) :: bn(0:n)

    call bernoulli(n, bn)

    bernoa = bn(n)

  end function bernoa

!ccccccccccccccc

  subroutine bernoulli(n, bn)

!*****************************************************************************80
!
!! BERNOA computes the Bernoulli number Bn.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    11 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index.
!
!    Output, real (dp) BN, the value of the N-th Bernoulli number.
!
  implicit none

   integer ( kind = 4 ), intent(in) :: n
   real (dp), intent(out) :: bn(0:n)

   real (dp) :: r, s
   integer   :: j, k, m

   bn(0) = 1._dp; bn(1) = -0.5_dp

   do m = 2, n
     s = - ( 1._dp / ( m + 1 ) - 0.5 )
     do k = 2, m - 1
       r = 1._dp
       do j = 2, k
         r = r * ( j + m - k ) / j
       end do
     s = s - r * bn(k)
    end do
    bn(m) = s
   end do

   do m = 3, n, 2
     bn(m) = 0
   end do

  end subroutine bernoulli

end module Poly

!ccccccccccccccc
!C     ALGORITHM 490 COLLECTED ALGORITHMS FROM ACM.
!C     ALGORITHM APPEARED IN COMM. ACM, VOL. 18, NO. 4N,
!C     P. 200.

real (dp) function DILOG(X)
 use constants, only: dp

!C REAL PART OF THE DILOGARITHM FUNCTION FOR A REAL
!C ARGUMENT. REF. NO. 1=L. LEWIN, *DILOGARITHMS +
!C ASSOCIATED FUNCTIONS*
!C                      (MAC-DONALD, LONDON, 1958).
!C NUMERICAL CONSTANTS USED ARE C(N)=(N(N+1)(N+2))**2
!C FOR N=1 TO 30, (PI**2)/3=3.289868...,
!C (PI**2)/6=1.644394..., AND ZERO OF DILOG ON THE
!C POSITIVE REAL AXIS, X0=12.59517...

  implicit none
  real (dp), intent(in) :: x
  real (dp) :: A, B, BY, C1, C2, C3, C4, DX, DY, TEST, W, X0, Y, Z
  integer :: N

  real (dp), dimension(30) :: C = [36._dp, 576._dp, 36.d2, 144.d2, 441.d2, 112896._dp, &
    254016._dp, 5184.d2, 9801.d2, 17424.d2, 2944656._dp, 4769856._dp, 74529.d2, 112896.d2,   &
    166464.d2, 23970816._dp, 33802596._dp, 467856.d2, 636804.d2, 853776.d2, 112911876._dp,   &
    147476736._dp, 19044.d4, 24336.d4, 3080025.d2, 386358336._dp, 480661776._dp, 5934096.d2, &
    7273809.d2, 8856576.d2]

  if (X > 12.6_dp) go to 10;  if (X >= 12.59_dp) go to 100;  if (X >= 2._dp) go to 10
  if (X > 1._dp) go to 20;  if (X == 1._dp) go to 30;  if (X > .5_dp) go to 40
  if (X > 1e-2_dp) go to 50;  if (X < -1._dp) go to 60;  if (X < -1e-2_dp) go to 70

!C DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(1).
      DILOG = X*(1._dp+X*(.25_dp+X*(1._dp/9._dp+X* (625e-4_dp+X*(4e-2_dp+X*(1._dp/36._dp+X*(1._dp/ &
   49._dp+X/64._dp)))))))
      RETURN
!C DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(6),
!C AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   10 Y = 1._dp/X
      BY = -1._dp - Y*(4._dp+Y)
      DILOG = 3.28986813369645287_dp - 0.5_dp*log(X)**2 + (Y*(4._dp+5.75_dp*Y)+3._dp* &
   (1._dp+Y)*(1._dp-Y)*log(1._dp-Y))/BY
      if (DILOG+4._dp*Y == DILOG) RETURN
      go to 80
!C DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7) WITH
!C X=1/X + EQ(6), AND DESCRIPTION OF THIS ALGORITHM,
!C EQ(4).
   20 Y = 1._dp - 1._dp/X
      DX = log(X)
      BY = 1._dp + Y*(4._dp+Y)
      DILOG = 1.64493406684822643_dp + &
   DX*(.5_dp*DX-log(X-1._dp)) + &
   (Y*(4._dp+5.75_dp*Y)-3._dp*(1._dp+Y)*DX/X)/BY
      go to 80
!C DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(2).
   30 DILOG = 1.64493406684822643_dp
      RETURN
!C DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7),
!C AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   40 Y = 1._dp - X
      DX = log(X)
      BY = -1._dp - Y*(4._dp+Y)
      DILOG = 1.64493406684822643_dp - DX*log(Y) + (Y*(4._dp+5.75_dp*Y)+3._dp*(1._dp+Y)*DX*X)/BY
      go to 80
!C DILOG COMPUTED FROM DESCRIPTION OF THIS ALGORITHM,
!C EQ(4)
   50 Y = X
      BY = 1._dp + Y*(4._dp+Y)
      DILOG = (Y*(4._dp+5.75_dp*Y)+3._dp*(1._dp+Y)* &
   (1._dp-Y)*log(1._dp-Y))/BY
      go to 80
!C DILOG COMPUTED FROM REF. NO. 1, P.245, EQ(12) WITH
!C X=-X, AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   60 Y = 1._dp/(1._dp-X)
      DX = log(-X)
      DY = log(Y)
      BY = 1._dp + Y*(4._dp+Y)
      DILOG = -1.64493406684822643_dp + 0.5_dp*DY*(DY+2._dp*DX) + (Y*(4._dp+5.75_dp*Y) &
   + 3._dp*(1._dp+Y)*(1._dp-Y)*(DX+DY))/BY
      if (DILOG+4._dp*Y == DILOG) RETURN
      go to 80
!C DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(8),
!C AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   70 Y = X/(X-1._dp)
      DX = log(1._dp-X)
      BY = -1._dp - Y*(4._dp+Y)
      DILOG = (Y*(4._dp+5.75_dp*Y)-3._dp*(1._dp+Y)* (1._dp-Y)*DX)/BY - .5_dp*DX*DX
   80 B = 4._dp*Y*Y/BY
      DO 90 N=1,30
        B = B*Y
        A = B/C(N)
        TEST = DILOG
        DILOG = DILOG + A
        if (DILOG == TEST) RETURN
   90 CONTINUE
      RETURN
!C DILOG COMPUTED FROM TAYLOR SERIES ABOUT ZERO OF
!C DILOG, X0.
  100 X0 = 12.5951703698450184_dp
      Y = X/X0 - 1._dp
      Z = 1._dp/11.5951703698450184_dp
      W = Y*Z
      C1 = (3._dp*X0-2._dp)/6._dp
      C2 = ((11._dp*X0-15._dp)*X0+6._dp)/24._dp
      C3 = (((5.D1*X0-104._dp)*X0+84._dp)*X0-24._dp)/ 12.D1
      C4 = ((((274._dp*X0-77.D1)*X0+94.D1)*X0-54.D1)* X0+12.D1)/72.D1
      DILOG = Y*(1._dp-Y*(.5_dp-Y*(1._dp/3._dp-Y* (.25_dp-Y*(.2_dp-Y/6._dp)))))*log(Z) - &
      W*X0*Y*(.5_dp-W*(C1-W*(C2-W*(C3-W*C4))))

end function DILOG

!ccccccccccccccc

subroutine f90Elliptic3(phi, k, c, res)
  use constants, only: dp; implicit none

  real (dp), intent(in ) :: phi, k, c
  real (dp), intent(out) :: res

  real (dp)              :: Elliptic3

  res = Elliptic3(phi, k, c)

end subroutine f90Elliptic3

!ccccccccccccccc

real (dp) function Elliptic3(c, phi, k)
  use constants, only: dp; use carlson_elliptic_module
  implicit none

  real (dp), intent(in ) :: phi, k, c

  integer                      :: ier

  Elliptic3 = Sin(phi) * drf( Cos(phi)**2, 1 - k * Sin(phi)**2, 1._dp, ier ) + &
              1._dp/3 * Sin(phi)**3 * c * drj( Cos(phi)**2, 1._dp, 1 - k * Sin(phi)**2, &
                                              1 - c * Sin(phi)**2, ier )

end function Elliptic3

!ccccccccccccccc

subroutine f90Elliptic2(phi, k, res)
  use constants, only: dp; use carlson_elliptic_module
  implicit none

  real (dp), intent(in )               :: phi, k
  real (dp), intent(out), dimension(2) :: res

  integer                                    :: ier

  res(1) = Sin(phi) * drf( Cos(phi)**2, 1 - k * Sin(phi)**2, 1._dp, ier )

  res(2) = res(1) - 1._dp/3 * Sin(phi)**3 * k * drj( Cos(phi)**2, 1 - k * Sin(phi)**2, &
                                                    1._dp, 1._dp, ier )

end subroutine f90Elliptic2


!ccccccccccccccc

module hyper
 use adapt, only: dGauss; use constants, only: dp, Euler, Pi, Pi2, prec
 use DeriGamma; implicit none

 real (dp), parameter, dimension(39) :: ContList = [2._dp, 2._dp, 2.1666666666666665_dp, &
 2.3333333333333335_dp, 2.4833333333333334_dp, 2.6166666666666667_dp, 2.7357142857142858_dp, &
 2.842857142857143_dp, 2.9400793650793653_dp, 3.028968253968254_dp, 3.110786435786436_dp, &
 3.1865440115440116_dp, 3.2570568320568323_dp, 3.322990897990898_dp, 3.38489565989566_dp, &
 3.4432289932289932_dp, 3.4983760520525227_dp, 3.550663633751869_dp, 3.6003712360910503_dp,&
 3.6477396571436818_dp, 3.6929777523817773_dp, 3.7362677956718207_dp, 3.7777697719564056_dp,&
 3.8176248444201737_dp, 3.855958177753507_dp, 3.892881254676584_dp, 3.9284937902891195_dp, &
 3.962885324680654_dp, 3.9961365562077473_dp, 4.028320464253724_dp, 4.059503259952649_dp, &
 4.08974519543652_dp, 4.119101256042581_dp, 4.147621755151315_dp, 4.1753528475882895_dp, &
 4.202336974572417_dp, 4.228613250848693_dp, 4.254217802769035_dp, 4.279184064577402_dp]

 contains

!ccccccccccccccc

real (dp) function pFq(a, b, x)
  real (dp), dimension(:), intent(in) :: a, b
  real (dp)              , intent(in) :: x
  integer                             :: i
  real (dp)                           :: fac
  real (dp), dimension( size(a) )     :: tabA
  real (dp), dimension( size(b) )     :: tabB

  pFq = 1; tabA = a; tabB = b; fac = 1

  do i = 1, 100
    fac = fac * x/i * product(tabA) / product(tabB)
    pFq = pFq + fac; if ( abs(fac/pFq) <= prec ) return
    if (i > 3 .and. abs(fac) >= pFQ ) exit
    tabA = tabA + 1;  tabB = tabB + 1
  end do

  pFq = 0

end function pFq

!ccccccccccccccc

real (dp) function F32(w, x)
  real (dp), intent(in) :: w, x
  integer               :: i
  real (dp)             :: fac, fac1, fac2

  F32 = 1; fac = 1; fac1 = 1; fac2 = 1 - w

  do i = 1, 100
    fac = fac * (-x)/i * fac1**3 / fac2;  fac1 = fac1 + 1;  fac = fac / fac1
    if (i > 2 .and. abs(fac) >= F32 ) exit; F32 = F32 + fac
    if ( abs(fac) <= prec ) return        ; fac2 = fac2 + 1
  end do

  F32 = 0

end function F32

!ccccccccccccccc

real (dp) function F32Inf(w, x)
  real (dp)  , intent(in) :: w, x
  real (dp), dimension(2) :: PG
  real (dp)               :: fact, fact2, lg, corr
  integer                 :: i, NZ, IERR

  call DPSIFN(-w, 0, 1, 2, PG, NZ, IERR)

  fact = w/x; lg = log(x) - Euler

  F32Inf = - fact * ( Pi2 + 2 * lg**2 + 2 * PG(1) * ( Pg(1) + 2 * Lg ) - 2 * PG(2) )/4

  PG(1) = PG(1) + lg; fact2 = 1

  do i = 1, 39
    fact = - fact * (w + i)/x; PG(1) = PG(1) - 1/(i + w); fact2 = fact2 * max(i - 1,1)
    corr = fact * ( ContList(i) + PG(1) )/fact2/i**2
    if (i > 3 .and. abs(corr) >= F32Inf ) exit
    F32Inf = F32Inf + corr; if ( abs(corr/F32Inf) <= prec ) return
  end do

  F32Inf = 0

end function F32Inf

!ccccccccccccccc

real (dp) function H3F2Exact(w, x)
 real (dp), intent(in) :: w, x

   if ( abs(w) <= tiny(1._dp) ) then
     H3F2Exact = log(1 + x)/x; return
   end if

   H3F2Exact = F32(w,x)
   if ( abs(H3F2Exact) > 0 ) return

   H3F2Exact = F32Inf(w, x)
   if ( abs(H3F2Exact) > 0 ) return

   H3F2Exact = dGauss(HyperInt, 0._dp, 1._dp, prec)

   contains

!ccccccccccccccc

  real (dp) function HyperInt(y)
    real (dp), intent(in) :: y
    HyperInt = Hyper2F1(1._dp, 1._dp, 1 - w, - x * y)
  end function HyperInt

end function H3F2Exact

!ccccccccccccccc

real (dp) function hyper2f1(a, b, c, x)
 implicit none

 real (dp), intent(in) :: a, b, c, x

 call f90hiper2f1 (a, b, c, x, hyper2f1)

end function hyper2f1

!ccccccccccccccc

subroutine f90hiper2f1 ( ain, bin, cin, xin, hf )
!*****************************************************************************80
!
!! HYGFX evaluates the hypergeometric function F(A,B,C,X).
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program that the copyright
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!    Modified by Vicent Mateu 19 September 2016
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real (dp) A, B, C, X, the arguments of the function.
!    C must not be equal to a nonpositive integer.
!    X < 1.
!
!    Output, real HF, the value of the function.
!
  implicit none

  real (dp), intent(in)  :: ain, bin, cin, xin
  real (dp), intent(out) :: hf

  real (dp) :: a, a0, aa, b, bb, c, c0, c1, eps, f0, f1, g0, g1, g2, g3, ga, gabc, x, x1, &
               gam, gb, gbm, gc, gca, gcab, gcb, gm, hw, pa, pb, r, r0, r1, rm, rp, sm,&
               sp0, sp

  integer  :: j, k, m, nm
  logical  :: l0, l1, l2, l3, l4, l5

  a = ain;  b = bin;  c = cin;  x = xin

  l0 = ( c == aint(c) ) .and. ( c < 0 )
  l1 = ( 1 - x < 1.e-15_dp ) .and. (c - a - b <= 0)
  l2 = ( a == aint(a) ) .and. ( a < 0 )
  l3 = ( b == aint(b) ) .and. ( b < 0 )
  l4 = ( c - a == aint(c - a) ) .and. ( c - a <= 0 )
  l5 = ( c - b == aint(c - b) ) .and. ( c - b <= 0 )

  if ( l0 .or. l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYGFX - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    return
  end if

  if ( 0.95 < x ) then
    eps = 1e-8_dp
  else
    eps = 1e-15_dp
  end if

  if ( x == 0 .or. a == 0 .or. b == 0 ) then

    hf = 1
    return

  else if ( 1 - x == eps .and. 0 < c - a - b ) then

 gc = gamma(c)
 gcab = gamma(c - a - b)
 gca = gamma(c - a)
 gcb = gamma(c - b)
    hf = gc * gcab /( gca *gcb )
    return

  else if ( 1 + x <= eps .and. abs ( c - a + b - 1 ) <= eps ) then

    g0 = sqrt (pi) * 2._dp**( - a )
    g1 = gamma(c)
    g2 = gamma(1 + a / 2 - b)
    g3 = gamma(0.5 + a/2)
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1;  r = 1

    do k = 1, nm
      r = r * ( a + k - 1 ) * ( b + k - 1 ) / ( k * ( c + k - 1 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1;  r  = 1
    do k = 1, nm
      r = r * ( c - a + k - 1 ) * ( c - b + k - 1 ) / ( k * ( c + k - 1 ) ) * x
      hf = hf + r
    end do
    hf = ( 1 - x )**( c - a - b ) * hf
    return

  end if

  aa = a;  bb = b;  x1 = x
!
!  WARNING: ALTERATION OF INPUT ARGUMENTS A AND B, WHICH MIGHT BE CONSTANTS.
!
  if ( x < 0 ) then
    x = x / ( x - 1 )
    if ( a < c .and. b < a .and. 0 < b ) then
      a = bb; b = aa
    end if
    b = c - b
  end if

  if ( 0.75 <= x ) then

    gm = 0

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1e-15_dp ) then

      m = int ( c - a - b )
 ga = gamma(a)

gb = gamma(b)

gc = gamma(c)
gam = gamma(a + m)
gbm = gamma(b + m)
      call psi ( a, pa )
      call psi ( b, pb )

      if ( m /= 0 ) then
        gm = 1
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1;  r0 = 1; r1 = 1; sp0 = 0; sp = 0

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1 ) * ( b + k - 1 ) / ( k * ( k - m ) ) * ( 1 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1 / ( a + k - 1 ) &
            + 1 / ( b + k - 1 ) - 1 / real ( k, dp )
        end do

        f1 = pa + pb + sp0 + 2 * Euler + log ( 1 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1 - a ) / ( k * ( a + k - 1 ) ) + ( 1 - b ) / ( k * ( b + k - 1 ) )

          sm = 0
          do j = 1, m
            sm = sm + ( 1 - a ) / ( ( j + k ) * ( a + j + k - 1 ) ) + 1 / (b + j + k - 1)
          end do

          rp = pa + pb + 2 * Euler + sp + sm + log ( 1 - x )

          r1 = r1 * ( a + m + k - 1 ) * ( b + m + k - 1 ) &
            / ( k * ( m + k ) ) * ( 1 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1 ) * ( b - m + k - 1 ) / ( k * ( k - m ) ) * (1 - x)
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1 / real ( k, dp )
        end do

        f1 = pa + pb - sp0 + 2 * Euler + log ( 1 - x )

        do k = 1, 250

          sp = sp + ( 1 - a ) / ( k * ( a + k - 1 ) ) + ( 1 - b ) / ( k * ( b + k - 1 ) )

          sm = 0
          do j = 1, m
            sm = sm + 1 / real ( j + k, dp )
          end do

          rp = pa + pb + 2 * Euler + sp - sm + log ( 1 - x )

          r1 = r1 * ( a + k - 1 ) * ( b + k - 1 ) / ( k * ( m + k ) ) * ( 1 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

ga = gamma(a)
gb = gamma(b)
gc = gamma(c)
gca = gamma(c - a)
gcb = gamma(c - b)
gcab = gamma(c - a - b)
gabc = gamma(a + b - c)
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1 - x )**( c - a - b )
      hf = 0;  r0 = c0;  r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1 ) * ( b + k - 1 ) / ( k * ( a + b - c + k ) ) * ( 1 - x )

        r1 = r1 * ( c - a + k - 1 ) * ( c - b + k - 1 ) &
          / ( k * ( c - a - b + k ) ) * ( 1 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1

    if ( a < c .and. c < 2 * a .and. b < c .and. c < 2 * b ) then

      a0 = ( 1 - x )**( c - a - b );  a = c - a;  b = c - b

    end if

    hf = 1;  r = 1

    do k = 1, 250

      r = r * ( a + k - 1 ) * ( b + k - 1 ) / ( k * ( c + k - 1 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( x1 < 0 ) then

    x = x1;  c0 = 1 / ( 1 - x )**aa;  hf = c0 * hf

  end if

  a = aa; b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYGFX - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

end subroutine f90hiper2f1

!ccccccccccccccc

pure subroutine psi (x, ps)
  use constants

!*****************************************************************************80
!
!! PSI computes the PSI function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program that the copyright
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!    Modified by Vicent Mateu 19 September 2016
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real (dp) X, the argument.
!
!    Output, real (dp) PS, the value of the PSI function.
!
  implicit none

  real (dp), intent(in ) :: x
  real (dp), intent(out) :: ps

  real (dp), parameter   :: a1 = -0.83333333333333333e-1_dp
  real (dp), parameter   :: a2 =  0.83333333333333333e-2_dp
  real (dp), parameter   :: a3 = -0.39682539682539683e-2_dp
  real (dp), parameter   :: a4 =  0.41666666666666667e-2_dp
  real (dp), parameter   :: a5 = -0.75757575757575758e-2_dp
  real (dp), parameter   :: a6 =  0.21092796092796093e-1_dp
  real (dp), parameter   :: a7 = -0.83333333333333333e-1_dp
  real (dp), parameter   :: a8 =  0.4432598039215686_dp

  integer (kind = 4)           :: k, n
  real (dp)              :: x2, xa, s

  xa = abs(x); s = 0

  if ( x == aint(x) .and. x <= 0 ) then

    ps = 1.d+300;    return

  else if ( xa == aint(xa) ) then

    n = int(xa)
    do k = 1, n - 1
      s = s + 1/real( k, dp )
    end do

    ps = s - Euler

  else if ( xa + 0.5 == aint ( xa + 0.5 ) ) then

    n = int ( xa - 0.5 )

    do k = 1, n
      s = s + 1._dp / real ( 2 * k - 1, dp )
    end do

    ps = 2 * s - Euler - 1.386294361119891_dp

  else

    if ( xa < 10 ) then

      n = 10 - int ( xa )
      do k = 0, n - 1
        s = s + 1/( xa + real ( k, dp ) )
      end do

      xa = xa + real ( n, dp )

    end if

    x2 = 1._dp / ( xa * xa )

    ps = log ( xa ) - 0.5 / xa + x2 * ((((((( a8 * x2 + a7 ) * x2 + a6 ) &
      * x2 + a5 ) * x2 + a4 ) * x2 + a3 ) * x2 + a2 ) * x2 + a1 )

    ps = ps - s

  end if

  if ( x < 0 ) then
    ps = ps - pi * cos (pi * x) / sin (pi * x) - 1 / x
  end if

end subroutine psi

end module hyper

!ccccccccccccccc
