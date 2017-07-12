
module Adapt
  use iso_fortran_env, only: error_unit, dp => real64
  implicit none
  real (dp)    , parameter :: R1 = 1, HF = R1/2

  contains

! ccccc

  Subroutine DAdapt(F, A, B, NSEG, RELTOL, ABSTOL, RES, ERR)
    real (dp)    , intent(in ) :: A, B, RELTOL, ABSTOL
    integer      , intent(in ) :: NSEG

    real (dp)    , intent(out) :: RES, ERR

    real (dp)     , external :: f
    real (dp)                :: Bige, bin, root, TE, TERSS, TVALS, XHIB, XLOB, XNEW
    integer                  :: i, IBIG, iter, NSEGD

    integer              , parameter :: NDIM = 100

    real (dp), dimension(NDIM), save :: XLO, XHI, TVAL, TERS
    integer                   , save :: NTER = 0

!     RES = Estimated Integral of F from A to B,
!     ERR = Estimated absolute error on RES.
!     NSEG  specifies how the adaptation is to be done:
!         = 0   means use previous binning,
!        =1   means fully automatic, adapt until tolerance attained.
!        =n>1 means first split interval into n equal segments,
!             then adapt as necessary to attain tolerance.
!     The specified tolerances are:
!            relative: RELTOL ;  absolute: ABSTOL.
!        It stops when one OR the other is satisfied, or number of
!        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
!        can be set to zero, in which case only the other is used.

    if (NSEG <= 0)  then

      if (NTER == 0) then
        NSEGD = 1 ; go to 2
      end if

      TVALS = 0;  TERSS = 0

      do i = 1,NTER
         CALL DGS56P( F, XLO(i), XHI(i), TVAL(i), TE )
         TERS(i) = TE**2;  TVALS=TVALS+TVAL(i);  TERSS=TERSS+TERS(i)
      end do

      ROOT = SQRT(2*TERSS);  go to 9

    end if

    NSEGD = MIN(NSEG, NDIM)
2   XHIB = A
    BIN = (B-A)/NSEGD

    do i = 1, NSEGD
      XLO(i) = XHIB  ;  XLOB = XLO(i);  XHI(i)=XHIB+BIN;  if (i == NSEGD) XHI(i) = B
      XHIB = XHI(i); CALL DGS56P(F,XLOB,XHIB,TVAL(i),TE);  TERS(i)=TE**2
    end do

    NTER = NSEGD

    do ITER = 1, NDIM
      TVALS = TVAL(1);  TERSS = TERS(1)

      do i = 2, NTER
        TVALS = TVALS + TVAL(i);  TERSS = TERSS + TERS(i)
      end do

      ROOT= SQRT(2 * TERSS)
      if ( ROOT <= ABSTOL .or. ROOT <= RELTOL * ABS(TVALS) ) go to 9
      if (NTER == NDIM) go to 9;  BIGE = TERS(1);  IBIG = 1

      do i = 2, NTER
        if (TERS(i) .GT. BIGE) then
          BIGE = TERS(i);  IBIG=I
        end if
      end do

      NTER = NTER + 1;  XHI(NTER) = XHI(IBIG); XNEW = HF * ( XLO(IBIG) + XHI(IBIG) )
      XHI(IBIG) = XNEW;  XLO(NTER) = XNEW
      CALL DGS56P( F, XLO(IBIG), XHI(IBIG), TVAL(IBIG), TE ) ;  TERS(IBIG) = TE**2
      CALL DGS56P( F, XLO(NTER), XHI(NTER), TVAL(NTER), TE ) ;  TERS(NTER) = TE**2
    end do

9   RES = TVALS
    ERR = ROOT

  end Subroutine Dadapt

! ccccc

  Subroutine DGS56P(F, A, B, RES, ERR)
    real (dp)    , intent(in ) :: A, B
    real ( dp )  , external    :: f

    real (dp)    , intent(out) :: RES, ERR
    real (dp)                  :: RANG, E5, E6
    integer                          :: i

    real (dp)   , dimension(5) :: X5, W5
    real (dp)   , dimension(6) :: X6, W6

    x5 = [ 4.6910077030668004e-2_dp, 2.3076534494715846e-1_dp, 0.5_dp, &
           7.6923465505284154e-1_dp, 9.5308992296933200e-1_dp ]

    w5 = [ 1.1846344252809454e-1_dp, 2.3931433524968324e-1_dp, 2.8444444444444444e-1_dp, &
           2.3931433524968324e-1_dp, 1.1846344252809454e-1_dp ]

    x6 = [ 3.3765242898423989e-2_dp, 1.6939530676686775e-1_dp, 3.8069040695840155e-1_dp, &
           6.1930959304159845e-1_dp, 8.3060469323313225e-1_dp, 9.6623475710157601e-1_dp]

    w6 = [ 8.5662246189585178e-2_dp, 1.8038078652406930e-1_dp, 2.3395696728634552e-1_dp, &
           2.3395696728634552e-1_dp, 1.8038078652406930e-1_dp, 8.5662246189585178e-2_dp ]

    RANG = B-A; E5 = 0; E6 = 0

    do i = 1,5
      E5 = E5 + W5(i) * F( A + RANG * X5(i) ); E6 = E6 + W6(i) * F( A + RANG * X6(i) )
    end do

    E6  = E6 + W6(6) * F( A + RANG * X6(6) ); RES = HF * (E6 + E5) * RANG
    ERR = ABS( (E6 - E5) * RANG )

  end Subroutine DGS56P

! ccccc

  real (dp) FUNCTION dGauss(f, A, B, EPS)
    real (dp)   , intent(in) :: A, B, EPS
    real (dp)   ,   external :: f

    real (dp), dimension(12) :: W, X
    real (dp)                :: AA, BB, C1, C2, U, S8, S16, CONST
    logical                        :: MFLAG
    integer                        :: i

!     ******************************************************************
!
!     ADAPTIVE real (dp) GAUSSIAN QUADRATURE.
!
!     dGauss IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
!     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
!     EPS.
!
!     ******************************************************************

      DATA W / 0.1012285362903762591525313543_dp , 0.2223810344533744705443559944_dp , &
               0.3137066458778872873379622020_dp , 0.3626837833783619829651504493_dp , &
               0.2715245941175409485178057246e-1_dp, 0.6225352393864789286284383699e-1_dp, &
               0.9515851168249278480992510760e-1_dp, 0.1246289712555338720524762822_dp , &
               0.1495959888165767320815017305_dp , 0.1691565193950025381893120790_dp , &
               0.1826034150449235888667636680_dp , 0.1894506104550684962853967232_dp/

      DATA X / 0.9602898564975362316835608686_dp, 0.7966664774136267395915539365_dp, &
               0.5255324099163289858177390492_dp, 0.1834346424956498049394761424_dp, &
               0.9894009349916499325961541735_dp, 0.9445750230732325760779884155_dp, &
               0.8656312023878317438804678977_dp, 0.7554044083550030338951011948_dp, &
               0.6178762444026437484466717640_dp, 0.4580167776572273863424194430_dp, &
               0.2816035507792589132304605015_dp, 0.9501250983763744018531933543e-1_dp/

      DATA W / 0.101228536290376259_dp , 0.222381034453374471_dp , 0.313706645877887287_dp , &
               0.362683783378361983_dp , 0.271524594117540949e-1_dp, 0.622535239386478929e-1_dp, &
               0.951585116824927848e-1_dp, 0.124628971255533872_dp , 0.149595988816576732_dp , &
               0.169156519395002538_dp,  0.182603415044923589_dp , 0.189450610455068496_dp/

      DATA X / 0.960289856497536232_dp, 0.796666477413626740_dp, 0.525532409916328986_dp, &
               0.183434642495649805_dp, 0.989400934991649933_dp, 0.944575023073232576_dp, &
               0.865631202387831744_dp, 0.755404408355003034_dp, 0.617876244402643748_dp, &
               0.458016777657227386_dp, 0.281603550779258913_dp, 0.950125098376374402e-1_dp/

!C
!C     ******************************************************************
!C
!C  START.
    dGauss = 0;  if(B == A) RETURN;  CONST = 0.005/(B-A);  BB = A
!C
!C  COMPUTATIONAL LOOP.
1   AA = BB
    BB = B
2   C1 = (BB+AA)/2
    C2 = (BB-AA)/2;  S8 = 0

    do i = 1, 4
      U = C2*X(I);  S8 = S8 + W(I) * ( F(C1 + U) + F(C1 - U) )
    end do

    S8 = C2 * S8;  S16 = 0

    do I = 5, 12
            U=C2*X(I)
            S16=S16+W(I)*(F(C1+U)+F(C1-U))
    end do

    S16 = C2*S16;  if( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) go to 5;  BB = C1
    if ( 1 + ABS(CONST * C2) /= 1) go to 2;  dGauss = 0

    if (MFLAG) then
        WRITE(*,6)
    end if
!       if(.NOT. RFLAG) CALL ABEND
    return

5     dGauss = dGauss + S16; if(BB /= B) go to 1

6 FORMAT( 4X, 'FUNCTION dGauss ... TOO HIGH ACCURACY REQUIRED')

  end function dGauss

end module Adapt
