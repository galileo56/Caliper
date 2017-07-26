
!ccccccccccccccc

subroutine f90Pi0(zr, zi, res)
  use SigmaClass, only: Pi0; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = Pi0( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi0

!ccccccccccccccc

subroutine f90Pi0Der(i, zr, zi, res)
  use SigmaClass, only: Pi0Der; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  integer  , intent(in) :: i
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = Pi0Der( i, cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi0Der

!ccccccccccccccc

subroutine f90Pi1Der(i, zr, zi, res)
  use SigmaClass, only: Pi1Der; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  integer  , intent(in) :: i
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = Pi1Der( i, cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi1Der

!ccccccccccccccc

subroutine f90Pi1(zr, zi, res)
  use SigmaClass, only: Pi1; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = Pi1( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi1

!ccccccccccccccc

subroutine f90DeltaBottomCharm(r1, r2, res)
  use constants, only: dp; use AnomDimClass, only: DeltaBottomCharm; implicit none
  real (dp), intent(in ) :: r1, r2
  real (dp), intent(out) :: res

  res = DeltaBottomCharm(r1,r2)

end subroutine f90DeltaBottomCharm

!ccccccccccccccc

subroutine f90gammaRBottomCharm(r1, r2, res)
  use constants, only: dp; use AnomDimClass, only: GammaRBottomCharm; implicit none
  real (dp), intent(in ) :: r1, r2
  real (dp), intent(out) :: res

  res = GammaRBottomCharm(r1,r2)

end subroutine f90gammaRBottomCharm

!ccccccccccccccc

subroutine f90DeltaCharm2(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm2; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharm2(r)

end subroutine f90DeltaCharm2

!ccccccccccccccc

subroutine f90p2(r, res)
  use constants, only: dp; use AnomDimClass, only: p; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = p(r)

end subroutine f90p2

!ccccccccccccccc

subroutine f90p2Int(r, res)
  use constants, only: dp; use AnomDimClass, only: pInt; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = pInt(r)

end subroutine f90p2Int

!ccccccccccccccc

subroutine f90p2Double(r1, r2, res)
  use constants, only: dp; use AnomDimClass, only: p2Int; implicit none
  real (dp), intent(in ) :: r1, r2
  real (dp), intent(out) :: res

  res = p2Int(r1, r2)

end subroutine f90p2Double

!ccccccccccccccc

subroutine f90DeltaCharmNh(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNh; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNh(r)

end subroutine f90DeltaCharmNh

!ccccccccccccccc

subroutine f90DeltaCharmGlue(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmGlue; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmGlue(r)

end subroutine f90DeltaCharmGlue

!ccccccccccccccc

subroutine f90DeltaCharmGlueDer(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmGlueDer; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmGlueDer(r)

end subroutine f90DeltaCharmGlueDer

!ccccccccccccccc

subroutine f90DeltaCharmNhDer(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNhDer; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNhDer(r)

end subroutine f90DeltaCharmNhDer

!ccccccccccccccc

subroutine f90DeltaCharmNL(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNL; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNL(r)

end subroutine f90DeltaCharmNL

!ccccccccccccccc

subroutine f90DeltaCharmNLDer(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNLDer; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNLDer(r)

end subroutine f90DeltaCharmNLDer

!ccccccccccccccc

subroutine f90DeltaCharm3(nl, nh, r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm3; implicit none
  real (dp), intent(in ) :: r
  integer  , intent(in)  :: nl, nh
  real (dp), intent(out) :: res

  res = deltaCharm3(nl, nh, r)

end subroutine f90DeltaCharm3

!ccccccccccccccc

subroutine f90DeltaCharm3Der(nl, nh, r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm3Der; implicit none
  real (dp), intent(in ) :: r
  integer  , intent(in)  :: nl, nh
  real (dp), intent(out) :: res

  res = deltaCharm3Der(nl, nh, r)

end subroutine f90DeltaCharm3Der

!ccccccccccccccc

subroutine f90GammaRCharm3(nl, nh, r, res)
  use constants, only: dp; use AnomDimClass, only: GammaRCharm3; implicit none
  real (dp), intent(in ) :: r
  integer  , intent(in)  :: nl, nh
  real (dp), intent(out) :: res

  res = GammaRCharm3(nl, nh, r)

end subroutine f90GammaRCharm3

!ccccccccccccccc

subroutine f90GammaRCharm2(r, res)
  use constants, only: dp; use AnomDimClass, only: gammaRcharm2; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = gammaRCharm2(r)

end subroutine f90GammaRCharm2

!ccccccccccccccc

subroutine f90DeltaCharm2Der(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm2Der; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharm2Der(r)

end subroutine f90DeltaCharm2Der

!ccccccccccccccc

subroutine f90LegendreList(n, k, x, res)
  use constants, only: dp; use Legendre; implicit none
  integer                  , intent(in)  :: n, k
  real (dp)                , intent(in)  :: x
  real (dp), dimension(0:n), intent(out) :: res

  res = LegendreList(n,k,x)

end subroutine f90LegendreList

!ccccccccccccccc

subroutine f90QLegendreList(n, x, res)
  use constants, only: dp; use Legendre; implicit none
  integer                  , intent(in)  :: n
  real (dp)                , intent(in)  :: x
  real (dp), dimension(0:n), intent(out) :: res

  res = QLegendreList(n,x)

end subroutine f90QLegendreList

!ccccccccccccccc

subroutine f90MCtop(shape, mt, Q, n, k, x, res)
  use constants, only: dp; use hyper; use MCtopClass; implicit none
  character (len = *), intent(in)  :: shape
  real (dp)          , intent(in)  :: mt, Q, x
  integer            , intent(in)  :: n, k
  real (dp)          , intent(out) :: res
  type (MCScales)                  :: MC

  MC = MCScales( shape(:6), n );     call MC%setMass(mt, Q)
  res = MC%Distribution(k,x)

end subroutine f90MCtop

!ccccccccccccccc

subroutine f90BreitUnstable(shape, mt, Q, gamma, n, k, x, res)
  use constants, only: dp; use hyper; use MCtopClass; implicit none
  character (len = *), intent(in)  :: shape
  real (dp)          , intent(in)  :: mt, Q, x, gamma
  integer            , intent(in)  :: n, k
  real (dp)          , intent(out) :: res
  type (MCtop)                     :: MC

  MC = MCtop( shape(:6), mt, Q, n );  res = MC%BreitUnstable(gamma, k,x)

end subroutine f90BreitUnstable

!ccccccccccccccc

subroutine f90ModelUnstable(shape, mt, Q, c, clen, lambda, n, k, p, res)
  use constants, only: dp; use hyper; use MCtopClass; use ModelClass; implicit none
  character (len = *)       , intent(in)  :: shape
  real (dp)                 , intent(in)  :: mt, Q, p, lambda
  real (dp), dimension(clen), intent(in)  :: c
  integer                   , intent(in)  :: n, clen, k
  real (dp)                 , intent(out) :: res
  type (MCtop)                            :: MC
  type (Model)                            :: Mod

  MC = MCtop( shape(:6), mt, Q, n ); Mod = Model(lambda, c, [0, 0], 'sum')
  res = Mod%ModelUnstable(MC, k, p)

end subroutine f90ModelUnstable

!ccccccccccccccc

subroutine f90BreitModelUnstable(shape, mt, Q, gamma, c, clen, lambda, n, k, p, res)
  use constants, only: dp; use hyper; use MCtopClass; use ModelClass; implicit none
  character (len = *)       , intent(in)  :: shape
  real (dp)                 , intent(in)  :: mt, Q, p, lambda, gamma
  real (dp), dimension(clen), intent(in)  :: c
  integer                   , intent(in)  :: n, clen, k
  real (dp)                 , intent(out) :: res
  type (MCtop)                            :: MC
  type (Model)                            :: Mod

  MC = MCtop( shape(:6), mt, Q, n ); Mod = Model(lambda, c, [0, 0], 'sum')
  res = Mod%BreitModUns(MC, gamma, k, p)

end subroutine f90BreitModelUnstable

!ccccccccccccccc

subroutine f90ModelUnstablediff(shape, mt, Q, c, clen, lambda, n, k, p, p2, res)
  use constants, only: dp; use hyper; use MCtopClass; use ModelClass; implicit none
  character (len = *)       , intent(in)  :: shape
  real (dp)                 , intent(in)  :: mt, Q, p, p2, lambda
  real (dp), dimension(clen), intent(in)  :: c
  integer                   , intent(in)  :: n, clen, k
  real (dp)                 , intent(out) :: res
  type (MCtop)                            :: MC
  type (Model)                            :: Mod

  MC = MCtop( shape(:6), mt, Q, n ); Mod = Model(lambda, c, [0, 0], 'sum')
  res = Mod%ModelUnstable(MC, k, p, p2)

end subroutine f90ModelUnstablediff

!ccccccccccccccc

subroutine f90DeltaMCtop(shape, mt, Q, res)
  use constants, only: dp; use hyper; use MCtopClass; implicit none
  character (len = *), intent(in)  :: shape
  real (dp)          , intent(in)  :: mt, Q
  real (dp)          , intent(out) :: res
  type (MCtop)                     :: MC

  MC = MCtop( shape(:6), mt, Q );  res = MC%delta()

end subroutine f90DeltaMCtop

!ccccccccccccccc

subroutine f90pfq(A, IP, B, IQ, Z, res)
  use constants, only: dp; use hyper; implicit none
  real (dp), dimension(IP), intent(in)  :: A
  real (dp), dimension(IQ), intent(in)  :: B
  integer                 , intent(in)  :: IP, IQ
  real (dp)               , intent(in)  :: Z
  real (dp)               , intent(out) :: res

  res = pFq(a, b, z)

end subroutine f90pfq

!ccccccccccccccc

subroutine f90PolyGamma(n, w, res)
  use constants, only: dp; use DeriGamma; implicit none
  integer      , intent(in) :: n
  real (dp)    , intent(in) :: w
  real (dp), dimension(0:n) :: res
  integer                   ::  NZ, IERR

  call DPSIFN(w, 0, 1, n + 1, res, NZ, IERR)

end subroutine f90PolyGamma

!ccccccccccccccc

subroutine f90MassiveProfList(terms, hard, shape, EShape, setup, gap, space, cum, scheme, &
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia
  use AnomDimClass; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, &
  cum, abs, current, Eshape, terms
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, ns, &
  runMass, orderMass, n, clen
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, &
  muC, Q, G3, lambda, R0, mu0, delta0, h, s3, muR0, j3, Rat0, n0, n1, t2, cnt, eS,&
  slope, eH, eJ, width, ts, xi, xiB, muLambda2, beta, deltaLambda, del0, delta1, &
  muM, mass, gammaZ, sin2ThetaW

  real (dp), dimension(n), intent(in)   :: tauList
  real (dp), dimension(n), intent(out)  :: res
  type (AnomDim), dimension(3:6)        :: AnDim
  type (ElectroWeak)                    :: EW
  type (Alpha)                          :: alphaAll
  type (MatricesElementsMass)           :: MatEl
  type (SingularMassScales)             :: Sing
  type (Model)                          :: Mod
  type (MassiveScales)                  :: nonSing
  type (ProfilesPythia)                 :: Prof
  type (CumulantMass)                   :: Cumul
  character (len = 5)                   :: alphaScheme
  integer                               :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, &
  delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), &
  current(:6), orderMass, matEl, EW )

  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListDist(terms(:7), cum(:4), Mod, setup(:15), gap(:12), space(:6), &
  order, R0, muR0, del0, h, tauList)

end subroutine f90MassiveProfList

!ccccccccccccccc

subroutine f90MassiveProfPieceList(terms, hard, shape, EShape, setup, gap, space, cum, scheme,  &
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,&
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,&
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, &
  cum, abs, current, setup, Eshape
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, &
  runMass, orderMass, n, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC,&
  G3, lambda, R0, mu0, delta0, s3, muR0, j3, muC, Q, Rat0, n0, n1, t2, cnt, eS,  &
   slope, eH, eJ, width, ts, xi, xiB, muLambda2, beta, deltaLambda, del0, delta1,&
   muM, mass, gammaZ, sin2ThetaW, h
  real (dp), dimension(n), intent(in)   :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out)  :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListDistPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, tauList)

end subroutine f90MassiveProfPieceList

!ccccccccccccccc

subroutine f90MassivePieceBin(terms, hard, shape, EShape, setup, gap, space, cum, scheme,  &
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,&
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,&
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, Eshape, &
                                        cum, abs, current, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, &
                                        runMass, orderMass, n
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q,&
                                        G3, lambda, R0, mu0, delta0, h, s3, muR0, j3, &
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW
  real (dp), dimension(2,n), intent(in) :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out) :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListBinPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, tauList)

end subroutine f90MassivePieceBin

!ccccccccccccccc

subroutine f90MassiveBinList(terms, hard, shape, EShape, setup, gap, space, cum, scheme, &
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *)        , intent(in) :: shape, scheme, setup, space, gap, &
  cum, abs, current, Eshape, terms, hard
  integer                    , intent(in) :: orderAlpha, order, runAlpha, run, &
  runMass, orderMass, n, clen, nf, ns
  real (dp)                  , intent(in) :: mZ, amZ, muLambda1, mT, muT, mB,   &
  muC, Q, G3, lambda, R0, mu0, delta0, h, s3, muR0, j3, Rat0, n0, n1, t2, cnt,  &
  slope, eH, eJ, width, ts, xiB, muLambda2, beta, deltaLambda, del0, delta1,&
  muM, mass, gammaZ, sin2ThetaW, muB, mC, eS, xi
  real (dp), dimension(2,n), intent(in) :: tauList
  real (dp), dimension(n) , intent(out) :: res
  integer                               :: i
  type (AnomDim), dimension(3:6)        :: AnDim
  type (ElectroWeak)                    :: EW
  type (Alpha)                          :: alphaAll
  type (MatricesElementsMass)           :: MatEl
  type (SingularMassScales)             :: Sing
  type (Model)                          :: Mod
  type (MassiveScales)                  :: nonSing
  type (ProfilesPythia)                 :: Prof
  type (CumulantMass)                   :: Cumul
  character (len = 5)                   :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
  ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
  orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListBin(terms(:7), cum(:9), Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, del0, h, tauList)

end subroutine f90MassiveBinList

!ccccccccccccccc

subroutine f90MassiveProf(terms, hard, shape, EShape, setup, gap, space, cum, scheme, &
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, terms, &
                                        cum, abs, current, Eshape
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (Model)                       :: Mod
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%Bin(terms(:7), cum(:4), Mod, setup(:15), gap(:12), space(:6), order, R0, &
  muR0, del0, h, 0, tau)

end subroutine f90MassiveProf

!ccccccccccccccc

subroutine f90MassOrigin(shape, EShape, gap, scheme, orderAlpha, runAlpha, orderMass,  &
 runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, &
 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, &
 mass, muM, R0, muR0, del0, h, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ElectroWeakClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia; use AnomDimClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, gap, Eshape
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        R0, mu0, delta0, h, muR0, Rat0, n0, n1, t2, cnt,  &
                                        eS, slope, eH, ts, muLambda2, beta, deltaLambda,  &
                                        del0, delta1, muM, mass, eJ
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (MatricesElementsMass)        :: MatEl
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, 0._dp, 0._dp)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, 0, EShape(:1), shape(:6) )
  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), 'abs', 'all', &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, 'expand' )
  Cumul     = CumulantMass(Prof, Sing, 0._dp, 0._dp, 0._dp)

  res = Cumul%FindOrigin( gap(:12), order, R0, muR0, del0, h )

end subroutine f90MassOrigin

!ccccccccccccccc

subroutine f90MassiveProfPiece(terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, &
 current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, &
 mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,    &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, &
                                        cum, abs, current, Eshape, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%BinPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, 0, tau)

end subroutine f90MassiveProfPiece

!ccccccccccccccc

subroutine f90MassiveProfDiffPiece(terms, hard, shape, EShape, setup, gap, space, cum, scheme,  &
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,&
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,&
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, &
                                        cum, abs, current, Eshape, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW, tau2
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%BinPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, 0, tau, tau2)

end subroutine f90MassiveProfDiffPiece

!ccccccccccccccc

subroutine f90MassiveProfDiff(terms, hard, shape, EShape, setup, gap, space, cum, scheme,&
 abs, current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, terms, &
                                        cum, abs, current, Eshape
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW, tau2
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (Model)                       :: Mod
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%Bin(terms(:7), cum(:9), Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, del0, h, 0, tau, tau2)

end subroutine f90MassiveProfDiff

!ccccccccccccccc

subroutine f90MassiveMoment(terms, hard, shape, EShape, setup, gap, space, scheme, abs,  &
 current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,     &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, pow, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, terms, &
                                        abs, current, Eshape
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, &
                                        runMass, orderMass, pow
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW, tau2
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (Model)                       :: Mod
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%Bin(terms(:7), 'Integrate', Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, del0, h, pow, tau, tau2)

end subroutine f90MassiveMoment

!ccccccccccccccc

subroutine f90MasslessProfList(terms, hard, shape, setup, gap, space, cum, orderAlpha, &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda,   &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda,    &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ
  real (dp), dimension(n), intent(in ) :: tauList
  real (dp), dimension(n), intent(out) :: res
  type (Alpha)                        :: alphaAll
  type (MatricesElements)             :: MatEl
  type (SingularScales)               :: Sing
  type (Model)                        :: Mod
  type (ProfilesMassless)             :: Prof
  type (CumulantMassless)             :: Cumul
  type (AnomDim), dimension(3:6)      :: AnDim
  integer                             :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, &
  eS, eJ, eR, ns)

  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListDist( terms(:7), cum(:4), Mod, setup(:15), gap(:12), space(:6), &
  order, R0, muR0, delta0, h, tauList )

end subroutine f90MasslessProfList

!ccccccccccccccc

subroutine f90MasslessProfPieceList(terms, hard, shape, setup, gap, space, cum, orderAlpha, &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda,     &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ
  real (dp), dimension(n), intent(in ) :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out) :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListDistPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, delta0, h, tauList )

end subroutine f90MasslessProfPieceList

!ccccccccccccccc

subroutine f90MasslessPieceBin(terms, hard, shape, setup, gap, space, cum, orderAlpha, &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda,   &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda,    &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS, &
  eR, tR, slope, eH, eJ, muC, j3
  real (dp), dimension(2,n), intent(in ) :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out) :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListBinPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, delta0, h, tauList )

end subroutine f90MasslessPieceBin

!ccccccccccccccc

subroutine f90MasslessBinList(terms, hard, shape, setup, gap, space, cum, orderAlpha, &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda,  &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda,   &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC,&
  Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS, &
  eR, tR, slope, eH, eJ, muC, j3
  real (dp), dimension(2,n), intent(in ) :: tauList
  real (dp), dimension(n)  , intent(out) :: res
  type (Alpha)                           :: alphaAll
  type (MatricesElements)                :: MatEl
  type (SingularScales)                  :: Sing
  type (Model)                           :: Mod
  type (ProfilesMassless)                :: Prof
  type (CumulantMassless)                :: Cumul
  type (AnomDim), dimension(3:6)         :: AnDim
  integer                                :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, &
  eS, eJ, eR, ns)

  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListBin( terms(:7), cum(:9), Mod, setup(:8), gap(:12), space(:6), order, &
                       R0, muR0, delta0, h, tauList )

end subroutine f90MasslessBinList

!ccccccccccccccc

subroutine f90MasslessProf(terms, hard, shape, setup, gap, space, cum, orderAlpha,   &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda,  &
 R0, muR0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS, &
  eR, tR, slope, eH, eJ, ts, muC, j3
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, &
  eS, eJ, eR, ns)

  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, &
  muT, mB, muB, mC, muC)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%Bin( terms(:7), cum(:4), Mod, setup(:8), gap(:12), space(:6), order, R0, &
  muR0, delta0, h, 0, tau )

end subroutine f90MasslessProf

!ccccccccccccccc

subroutine f90MasslessProfPiece(terms, hard, shape, setup, gap, space, cum, orderAlpha,   &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda,  &
 R0, muR0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ, &
                                        ts
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  integer                            :: i, j, k
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%BinPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), &
  order, R0, muR0, delta0, h, 0, tau )

end subroutine f90MasslessProfPiece

!ccccccccccccccc

subroutine f90MasslessProfDiffPiece(terms, hard, shape, setup, gap, space, cum, orderAlpha, &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda, R0, &
 muR0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ, &
                                        ts, tau2
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  integer                            :: i, j, k
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%BinPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), &
  order, R0, muR0, delta0, h, 0, tau, tau2 )

end subroutine f90MasslessProfDiffPiece

!ccccccccccccccc

subroutine f90MasslessProfDiff(terms, hard, shape, setup, gap, space, cum, orderAlpha,  &
 runAlpha, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda, R0,    &
 muR0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS,   &
  eR, tR, slope, eH, eJ, ts, tau2, muC, j3
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%Bin(terms(:7), cum(:9), Mod, setup(:15), gap(:12), space(:6), &
  order, R0, muR0, delta0, h, 0, tau, tau2)

end subroutine f90MasslessProfDiff

!ccccccccccccccc

subroutine f90MasslessMoment(terms, hard, shape, setup, gap, space, orderAlpha, runAlpha,&
 order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, &
 n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda, R0, muR0, delta0,  &
 h, tau, tau2, pow, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen, ns, pow
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ, &
                                        ts, tau2
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%Bin(terms(:7), 'Integrate', Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, delta0, h, pow, tau, tau2)

end subroutine f90MasslessMoment

!ccccccccccccccc

subroutine f90FindOrigin(shape, gap, orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, &
 muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, &
 eR, R0, muR0, delta0, h, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, gap
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, &
                                        Q, R0, mu0, delta0, h, muR0, ts, eH, slope, n0, &
                                        Rat0, n1, t2, cnt, eS, eR, tR
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim = AnomDim('MSbar', i, 0._dp)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, 0._dp, eR, 0)
  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl    = MatricesElements(alphaAll, nf, 0._dp, 0._dp, 0._dp, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), 'expand' )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%FindOrigin( gap(:12), order, R0, muR0, delta0, h)

end subroutine f90FindOrigin

!ccccccccccccccc

subroutine f90Profiles(Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, &
                       tau, res)
  use constants, only: dp; use ProfilesClass, only: profilesmassless; implicit none
  integer                 , intent(in) :: ns
  real (dp)               , intent(in) :: mu0, R0, n0, n1, t2, cnt, eS, eR, tR, &
  Q, eH, eJ, tau, ts, slope
  real (dp), dimension(5), intent(out) :: res

  real (dp)                            :: muS, muJ, muH, muNS, R
  type (ProfilesMassless)              :: Prof

  Prof = ProfilesMassless(Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)

  call Prof%Scales(tau, muNS, muH, muJ, muS, R)

  res = [ muH, muJ, muS, muNS, R ]

end subroutine f90Profiles

!ccccccccccccccc

subroutine f90ProfilesMass(Q, beta, mu0, deltaLambda, R0, n0, delta0, n1, delta1, &
  t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, def, EShape, tau, res)

  use constants, only: dp; use ProfilesClass, only: profilesPythia; implicit none

  character (len = *)    , intent(in ) :: EShape, def
  integer                , intent(in ) :: ns
  real (dp)              , intent(in ) :: Q, mu0, R0, n0, delta0, n1, t2, cnt, ts, slope, &
                                          eJ, mass, beta, delta1, deltaLambda, muM,   &
                                          eH, eS, tau
  real (dp), dimension(6), intent(out) :: res
  real (dp)                            :: muS, muJ, muH, muNS, R
  type (ProfilesPythia)                :: Prof

  Prof = ProfilesPythia(Q, beta, mu0, deltaLambda, R0, n0, delta0, n1, delta1, t2, ts, &
                        slope, cnt, eH, eS, eJ, mass, muM, ns, def, EShape)

  call Prof%Scales(tau, muNS, muH, muJ, muS, R)

  res = [ muH, muJ, muS, mum, muNS, R ]

end subroutine f90ProfilesMass

!ccccccccccccccc

subroutine f90SingularMass(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,      &
  current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, &
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,   &
  muM, mu, width, c, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass     ; use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass    ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, space, gap, Eshape, abs, hard,&
                                        current
  integer            , intent(in   ) :: orderAlpha, order, runAlpha, run, nf, clen,      &
                                        runMass, orderMass
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, xiB, gammaZ, sin2ThetaW, width, mC, muC, j3, &
                                        Rmass
  real (dp)          , intent(inout) :: muH, muJ, muS, R, muM
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC )

  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), &
  current(:6), orderMass, matEl, EW )

  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model( lambda, c, [0,0], 'sum' )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                               mu0, delta0, h, tau )
  else

    res = Sing%SingleSingWidth( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, setup(:15), gap(:12), space(:6), cum(:4), &
                                           order, R0, mu0, delta0, h, tau)

end subroutine f90SingularMass

!ccccccccccccccc

subroutine f90SingularMassDiff(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, &
  current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3,&
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,  &
  muM, mu, width, c, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, space, gap, Eshape, abs, hard,&
                                        current
  integer            , intent(in   ) :: orderAlpha, order, runAlpha, run, nf, clen,      &
                                        runMass, orderMass
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, xiB, gammaZ, sin2ThetaW, width, mC, muC, j3, &
                                        tau2, Rmass
  real (dp)          , intent(inout) :: muH, muJ, muS, R, muM
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model( lambda, c, [0,0], 'sum' )

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                               mu0, delta0, h, tau, tau2 )
  else

    res = Sing%SingleSingWidth( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau, tau2 )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, setup(:15), gap(:12), space(:6), cum(:4), order,&
                                               R0, mu0, delta0, h, tau, tau2)
end subroutine f90SingularMassDiff

!ccccccccccccccc

subroutine f90SingularMassPiece(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,        &
  current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, &
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,   &
  muM, mu, width, piece, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, abs, current,  &
                                        hard, scheme, setup
  integer              , intent(in ) :: orderAlpha, order, runAlpha, run, nf, runMass, &
                                        orderMass
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, gammaZ, sin2ThetaW, width, mC, muC, j3, muH, &
                                        muJ, muS, R, muM, xiB, Rmass
  real (dp)            , intent(out) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
  muT, mB, muB, mC, muC )

  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
  orderMass, matEl, EW )

  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], piece, 'piece')

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0,   &
                           delta0, h, tau )
  else

    res = Sing%SingleSingWidth( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0, delta0, &
                                h, tau )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, 'Model', gap(:12), space(:6), cum(:4), order,&
  R0, mu0, delta0, h, tau)

end subroutine f90SingularMassPiece

!ccccccccccccccc

subroutine f90SingularMassList(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,&
  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,   &
  width, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, abs, current,  &
                                        hard, scheme, setup
  integer              , intent(in ) :: orderAlpha, order, runAlpha, run, nf, runMass, &
                                        orderMass, clen
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, gammaZ, sin2ThetaW, width, mC, muC, j3, muH, &
                                        muJ, muS, R, muM, xiB, Rmass

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  integer                                              :: i, j, k
  type (ElectroWeak)                                   :: EW
  type (Alpha)                                         :: alphaAll
  type (MatricesElementsMass)                          :: MatEl
  type (SingularMassScales)                            :: Sing
  type (MassiveScales)                                 :: nonSing
  type (AnomDim), dimension(3:6)                       :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)                                  :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, gap(:12), space(:6), cum(:4), order, R0,   &
    mu0, delta0, h, tau )
  else

    res = Sing%SingleSingWidth( Mod, setup, gap(:12), space(:6), cum(:4), order, R0,   &
    mu0, delta0, h, tau )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, gap(:12), space(:6), cum(:4), order, &
  R0, mu0, delta0, h, tau)

end subroutine f90SingularMassList

!ccccccccccccccc

subroutine f90SingularMassDiffList(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,&
  current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, &
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,   &
  muM, mu, width, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, current,  &
  hard, scheme, setup, abs
  integer              , intent(in ) :: orderAlpha, order, runAlpha, run, nf, runMass, &
                                        orderMass, clen
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, gammaZ, sin2ThetaW, width, mC, muC, j3, muH, &
                                        muJ, muS, R, muM, tau2, xiB, Rmass

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  integer                                              :: i, j, k
  type (ElectroWeak)                                   :: EW
  type (Alpha)                                         :: alphaAll
  type (MatricesElementsMass)                          :: MatEl
  type (SingularMassScales)                            :: Sing
  type (MassiveScales)                                 :: nonSing
  type (AnomDim), dimension(3:6)                       :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)                                  :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau, tau2 )
  else

    res = Sing%SingleSingWidth( Mod, setup, gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau, tau2 )
  end if
  if (muJ >= muM) res = res + Sing%NonDist(Mod, gap(:12), space(:6), cum(:4), order, &
                                               R0, mu0, delta0, h, tau, tau2)
end subroutine f90SingularMassDiffList

!ccccccccccccccc

subroutine f90SingularMassDiffPiece(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,   &
  current, xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3,&
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,  &
  muM, mu, width, piece, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, current,  &
  hard, scheme, setup, abs
  integer              , intent(in ) :: orderAlpha, order, runAlpha, nf, runMass, &
  orderMass, run
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, &
  mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, mB, muB, xi, gammaZ, &
  sin2ThetaW, width, mC, muC, j3, muH, muJ, muS, R, muM, tau2, xiB, Rmass
  real (dp)            , intent(out) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,    &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], piece, 'piece')

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0, &
                           delta0, h, tau, tau2 )
  else

    res = Sing%SingleSingWidth( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0, delta0, &
                                h, tau, tau2 )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, 'Model', gap(:12), space(:6), cum(:4), order, &
                                               R0, mu0, delta0, h, tau)
end subroutine f90SingularMassDiffPiece

!ccccccccccccccc

subroutine f90MassNonDist(hard, shape, Eshape, setup, gap, space, cum, scheme,  &
  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c, &
  clen, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass     ; use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass    ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, gap, hard, Eshape, space
  integer            , intent(in   ) :: orderAlpha, order, runAlpha, run, clen, runMass, &
                                        orderMass, nf
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, &
                                        mC, muC, muT
  real (dp)          , intent(inout) :: muH, muS, R, muM, muJ, Rmass
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW       = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl    = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
                                    muLambda1, muLambda2 )
  nonSing  = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                              orderMass, matEl, EW )
  Sing     = SingularMassScales( nonSing, run, hard(:6) )
  Mod      = Model( lambda, c, [0,0], 'sum' )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res = Sing%NonDist(Mod, setup(:15), gap(:12), space(:3), cum(:4), order, R0, &
                     mu0, delta0, h, tau)

end subroutine f90MassNonDist

!ccccccccccccccc

subroutine f90MassNonDistDiff(hard, shape, Eshape, setup, gap, space, cum, scheme, &
  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,  &
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c,  &
  clen, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, gap, hard, Eshape, space
  integer            , intent(in   ) :: orderAlpha, order, runAlpha, run, clen, runMass, &
                                        orderMass, nf
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, h, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, tau2,&
                                        mC, muC, muT
  real (dp)          , intent(inout) :: muH, muS, R, muM, muJ, Rmass
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model( lambda, c, [0,0], 'sum' )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res = Sing%NonDist(Mod, setup(:15), gap(:12), space(:3), cum(:4), order, R0, &
                     mu0, delta0, h, tau, tau2)

end subroutine f90MassNonDistDiff

!ccccccccccccccc

subroutine f90MassNonDistPiece(hard, shape, Eshape, gap, space, cum, scheme, &
  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, &
  mu, c, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer, dimension(2), intent(in ) :: c
  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlpha, runAlpha, run, nf, runMass, &
                                        orderMass, order
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, &
                                        mC, muC, muT, muH, muS, R, muM, muJ, Rmass
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,     &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], c, 'piece')

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, 'Model', gap(:12), space(:3), cum(:4), order, &
  R0, mu0, delta0, h, tau )

end subroutine f90MassNonDistPiece

!ccccccccccccccc

subroutine f90MassNonDistList(hard, shape, Eshape, gap, space, cum, scheme, &
  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, &
  mu, clen, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlpha, runAlpha, run, nf, runMass, &
                                        orderMass, clen, order
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, &
                                        mC, muC, muT, muH, muS, R, muM, muJ, Rmass
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out  ) :: res

  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (AnomDim), dimension(3:6)     :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT, &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, &
                            delta0, h, tau )

end subroutine f90MassNonDistList

!ccccccccccccccc

subroutine f90MassNonDistDiffList(hard, shape, Eshape, gap, space, cum, scheme,  &
  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,&
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,   &
  clen, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlpha, runAlpha, run, nf, runMass, &
                                        orderMass, clen, order
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muJ, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, tau2,   &
                                        mC, muC, muT, muH, muS, R, muM, muB, Rmass, h
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out  ) :: res

  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,     &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, delta0, h, tau, tau2 )

end subroutine f90MassNonDistDiffList

!ccccccccccccccc

subroutine f90MassNonDistDiffPiece(hard, shape, Eshape, gap, space, cum, scheme, &
  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,&
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c,&
  lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer, dimension(2), intent(in ) :: c
  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlpha, order, runAlpha, run, nf, runMass, &
                                        orderMass
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muJ, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, tau2,&
                                        mC, muC, muT, muH, muS, R, muM, muB, Rmass
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,     &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], c, 'piece')

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, 'Model', gap(:12), space(:3), cum(:4), order, R0, mu0, delta0, h, &
                            tau, tau2 )

end subroutine f90MassNonDistDiffPiece

!ccccccccccccccc

subroutine f90NSMass(shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h, &
  gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass;  use ModelClass
  use constants, only: dp; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in ) :: shape, current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nf, clen, runMass, &
                                      orderMass
  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, mu,&
                                      muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muJ,  &
                                      delta0, h, R, width, mB, muB, muM, muLambda2
  real (dp)          , intent(out) :: res
  type (MassiveScales)             :: nonSing
  type (ElectroWeak)               :: EW
  type (MatricesElementsMass)      :: MatEl
  type (Alpha)                     :: alphaAll
  type (Model)                     :: Mod
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, c, [0,0], 'sum')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, &
                    muT, mB, muB, mC, muC)
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:10), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass( Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, delta0, &
                        h, t )

end subroutine f90NSMass

!ccccccccccccccc

subroutine f90NSMassDiff(shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h, &
  gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use ModelClass; use AnomDimClass
  use constants, only: dp; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in ) :: shape, current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nf, clen, runMass, &
                                      orderMass
  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, t2, &
                                      muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muB, h,&
                                      delta0, R, width, mB, muJ, muM, muLambda2, mu
  real (dp)          , intent(out) :: res
  type (MassiveScales)             :: nonSing
  type (ElectroWeak)               :: EW
  type (MatricesElementsMass)      :: MatEl
  type (Alpha)                     :: alphaAll
  type (Model)                     :: Mod
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, c, [0,0], 'sum')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, &
  muT, mB, muB, mC, muC)
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
  muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), &
  orderMass, matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass( Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, &
   delta0, h, t, t2 )

end subroutine f90NSMassDiff

!ccccccccccccccc

subroutine f90NSMassPiece(shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha,   &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, piece, width, lambda, R0, mu0, delta0, h,   &
  gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme
  integer              , intent(in ) :: orderAlpha, runAlpha, order, run, nf, runMass,   &
                                        orderMass
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda ,&
                                        muJ, delta0, R, h, mT, mB, width, muM,      &
                                        muLambda2
  real (dp)           , intent(out) :: res
  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  character (len = 5)               :: alphaScheme
  integer                           :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, [1._dp], piece, 'piece')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, 'Model', gap(:12), cum(:4), order, &
                       run, R0, mu0, delta0, h, t)

end subroutine f90NSMassPiece

!ccccccccccccccc

subroutine f90NSMassList(shape, setup, gap, cum, scheme, abs, current, orderAlpha, &
runAlpha, orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC,  &
muLambda1, muLambda2, Q, mu, muM, muJ, muS, R, Rmass, clen, width, lambda, R0, mu0,&
delta0, h, gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme, setup
  integer              , intent(in ) :: orderAlpha, runAlpha, order, run, nf, runMass,   &
                                        orderMass, clen
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda , &
                                        muJ, delta0, R, h, mT, mB, width, muM,       &
                                        muLambda2

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (AnomDim), dimension(3:6)    :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)               :: alphaScheme
  integer                           :: i, j, k

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, delta0, h, t)

end subroutine f90NSMassList

!ccccccccccccccc

subroutine f90NSMassDiffList(shape, setup, gap, cum, scheme, abs, current, orderAlpha, &
runAlpha, orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, &
muLambda1, muLambda2, Q, mu, muM, muJ, muS, R, Rmass, clen, width, lambda, R0, mu0,&
delta0, h, gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme, setup
  integer              , intent(in ) :: orderAlpha, runAlpha, order, run, nf, &
                                        orderMass, runMass, clen
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda, &
                                        muJ, delta0, R, h, mT, mB, width, muM, t2, &
                                        muLambda2

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (AnomDim), dimension(3:6)    :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)               :: alphaScheme
  integer                           :: i, j, k

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, delta0, h, t, t2)

end subroutine f90NSMassDiffList

!ccccccccccccccc

subroutine f90NSMassDiffPiece(shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha,   &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, piece, width, lambda, R0, mu0, delta0, h,   &
  gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme
  integer              , intent(in ) :: orderAlpha, runAlpha, order, run, nf, runMass,   &
                                        orderMass
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda, &
                                        muJ, delta0, R, h, mT, mB, width, muM, t2,  &
                                        muLambda2
  real (dp)           , intent(out) :: res
  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  character (len = 5)               :: alphaScheme
  integer                           :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, [1._dp], piece, 'piece')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlpha, runAlpha, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, 'Model', gap(:12), cum(:4), order, &
                       run, R0, mu0, delta0, h, t, t2)

end subroutine f90NSMassDiffPiece

!ccccccccccccccc

subroutine f90HJMNSMass(setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha,     &
 orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, &
 muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h,  &
 gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass; use MatrixElementsClass
  use ModelClass         ; use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in ) :: current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nf, clen, &
  orderMass, runMass

  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, sin2ThetaW, &
  mu, muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muB, delta0, R, width, mB, &
  muJ, muM, muLambda2, h, gammaZ

  real (dp)             , intent(out) :: res
  integer                             :: i, j
  type (AnomDim), dimension(3:6)      :: AnDim
  type (MassiveNS)                    :: nonSing
  type (ElectroWeak)                  :: EW
  type (Alpha)                        :: alphaAll
  type (MatrixElementsMass)           :: MatEl
  type (Model), dimension(clen, clen) :: Mod
  character (len = 5)                 :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, mB, muB, mC, muC)
  MatEl     = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp,  &
  0._dp, Q, mu, muJ, muM, muS, R, Rmass, muLambda1, muLambda2)

  nonSing = MassiveNS('HJM', 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
  matEl, EW, width, mu)

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = nonSing%HJMNSMass(c, Mod, setup(:15), gap(:12), cum(:4), &
  order, run, R0, mu0, delta0, h, t)

end subroutine f90HJMNSMass

!ccccccccccccccc

subroutine f90HJMNSMassDiff(setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, &
 orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,     &
 muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h, &
 gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in ) :: current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlpha, runAlpha, order, run, clen, &
  orderMass, runMass, nf

  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, gammaZ, mu, &
  muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muB, h, delta0, R, width, mB, &
  muJ, muM, muLambda2, t2, sin2ThetaW

  real (dp)             , intent(out) :: res
  type (AnomDim), dimension(3:6)      :: AnDim
  integer                             :: i, j
  type (MassiveNS)                    :: nonSing
  type (ElectroWeak)                  :: EW
  type (Alpha)                        :: alphaAll
  type (MatrixElementsMass)           :: MatEl
  type (Model), dimension(clen, clen) :: Mod
  character (len = 5)                 :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, mB, muB, mC, muC)
  MatEl     = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, &
  0._dp, Q, mu, muJ, muM, muS, R, Rmass, muLambda1, muLambda2)

  nonSing = MassiveNS('HJM', 'Q', scheme(:7), abs(:3), current(:6), orderMass,  &
  matEl, EW,width, mu)

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = nonSing%HJMNSMass(c, Mod, setup(:15), gap(:12), cum(:4), &
  order, run, R0, mu0, delta0, h, t, t2)

end subroutine f90HJMNSMassDiff

!ccccccccccccccc

subroutine f90FOMass(shape, current, m, Q, mZ, gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ;  use ElectroWeakClass; use MatrixElementsClass
  use AlphaClass         ;  use constants, only: dp; use AnomDimClass
  implicit none

  real (dp)          , intent(in ) :: t, m, Q, mZ, gammaZ, sin2ThetaW
  character (len = *), intent(in ) :: shape, current
  real (dp)          , intent(out) :: res
  type (MassiveNS)                 :: nonSing
  type (ElectroWeak)               :: EW
  type (MatrixElementsMass)        :: MatEl
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, m, m, m, m, m, m)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, Q, &
  Q, Q, tiny(1._dp), Q, Q, Q, 0._dp, 0._dp)

  nonSing   = MassiveNS( shape(:6), 'Q', 'pole', 'noAbs', current(:6), 0, matEl, &
  EW, 0._dp, Q)

  res = nonSing%FOMass(t)

end subroutine f90FOMass

!ccccccccccccccc

subroutine f90DiffDeltaGapMass(gap, order, R0, R, mu0, mu, muM, muLambda1, muLambda2, &
  orderAlpha, runAlpha, runMass, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, res)

  use MatrixElementsClass; use AlphaClass  ; use GapMassClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: gap
  integer            , intent(in ) :: order, nf, orderAlpha, runAlpha, runMass
  real (dp)          , intent(in ) :: R0, R, mu0, mu, muLambda1, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, muLambda2, muM
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaAll
  type (MatrixElementsMass)        :: MatEl
  type (AnomDim), dimension(3:6)   :: AnDim
  type (GapMass)                   :: MassGap
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, mB, muB, mC, muC)

  MatEl    = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp,  &
  0._dp, mu, mu, mu, muM, mu, mu, mu, muLambda1, muLambda2)

  MassGap = GapMass(order, gap(:12), MatEl, R0, mu0)

  res = MassGap%DeltaGapMass(R, mu)

end subroutine f90DiffDeltaGapMass

!ccccccccccccccc

subroutine f90ThrustNS1loop(t, res)
  use NonSingularClass; use constants, only: dp; implicit none
  real (dp), intent(in )               :: t
  real (dp), intent(out), dimension(3) :: res

  res = ThrustNS1loop(t)

end subroutine f90ThrustNS1loop

!ccccccccccccccc

subroutine f90ThrustNS2loop(er, t, res)
  use NonSingularClass; use AlphaClass;  use RunningClass; use AnomDimClass
  use constants, only: dp; implicit none

  real (dp), intent(in )               :: t, er
  real (dp), intent(out), dimension(2) :: res
  type (Alpha)                         :: alphaAll
  type (Running)                       :: alphaMass
  type (NonSingular)                   :: nonSing
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = Running(5, 0, alphaAll, 0._dp)
  nonSing   = NonSingular('thrust', alphaMass, 0._dp, 0._dp, 0)
  res       = nonSing%NS2loop(er, t)

end subroutine f90ThrustNS2loop

!ccccccccccccccc

subroutine f90EWfactors(nf, Q, mZ, gammaZ, sin2ThetaW, res)
  use ElectroWeakClass;  use constants, only: dp; implicit none
  integer                , intent(in ) :: nf
  real (dp)              , intent(in ) :: Q, mZ, gammaZ, sin2ThetaW
  real (dp), dimension(2), intent(out) :: res
  type (ElectroWeak)                   :: EW

  EW = ElectroWeak(mZ, gammaZ, sin2ThetaW); res = EW%EWfactors(nf, Q)

end subroutine f90EWfactors

!ccccccccccccccc

subroutine f90PolyLog(n, z, res)
  use Poly;  use constants, only: dp; implicit none
  integer  , intent(in ) :: n
  real (dp), intent(in ) :: z
  real (dp), intent(out) :: res

  res = polylog(n,z)

end subroutine f90PolyLog

!ccccccccccccccc

subroutine f90DiLog(z, res)
  use constants, only: dp; implicit none
  real (dp), intent(in ) :: z
  real (dp), intent(out) :: res
  real (dp)              :: dilog

  res = dilog(z)

end subroutine f90DiLog

!ccccccccccccccc

subroutine f90ComplexPolyLog(n, zre, zim, res)
  use Poly; use constants, only: dp; implicit none
  integer  , intent(in ) :: n
  real (dp), intent(in ) :: zre, zim
  real (dp), intent(out) :: res(2)
  complex (dp)           :: pl

  pl  = ComplexPolyLog( n, cmplx(zre, zim, dp) )
  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90ComplexPolyLog

!ccccccccccccccc

subroutine f90NGLSoft(nf, zre, zim, res)
  use MatrixElementsClass;  use constants, only: dp; implicit none
  integer  , intent(in )     :: nf
  real (dp), intent(in )     :: zre, zim
  real (dp), intent(out)     :: res(2)
  complex (dp)               :: pl

  pl = NGLSoft( nf, cmplx(zre, zim, dp) );  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90NGLSoft

!ccccccccccccccc

subroutine f90CLi2(zre, zim, res)
  use Chaplin;  use constants, only: dp; implicit none
  real (dp), intent(in ) :: zre, zim
  real (dp), intent(out) :: res(2)
  complex (dp)           :: pl

  pl = Cli2( cmplx(zre, zim, dp) );  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90CLi2

!ccccccccccccccc

subroutine f90CLi3(zre, zim, res)
  use Chaplin;  use constants, only: dp; implicit none
  real (dp), intent(in ) :: zre, zim
  real (dp), intent(out) :: res(2)
  complex (dp)           :: pl

  pl = Cli3( cmplx(zre, zim, dp) );  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90CLi3

!ccccccccccccccc

subroutine f90Model(c, clen, lambda, k, l, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, l, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%ShapeFun(k, l)

end subroutine f90Model

!ccccccccccccccc

subroutine f90ModelDiff(c, clen, lambda, k, l, l2, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer                   , intent(in)  :: k, clen
  real (dp)                 , intent(in)  :: lambda, l, l2
  real (dp), dimension(clen), intent(in)  :: c
  real (dp)                 , intent(out) :: res
  type (Model)                            :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%ShapeFun(k, l, l2)

end subroutine f90ModelDiff

!ccccccccccccccc

subroutine f90BreitModel(c, clen, lambda, width, k, l, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, l, width, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%BreitModel(width, k, l)

end subroutine f90BreitModel

!ccccccccccccccc

subroutine f90BreitModelDiff(c, clen, lambda, width, k, l, l2, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, l, l2, width, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%BreitModel(width, k, l, l2)

end subroutine f90BreitModelDiff

!ccccccccccccccc

subroutine f90Taylor(c, clen, lambda, k, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%Taylor(k)

end subroutine f90Taylor

!ccccccccccccccc

subroutine f90TaylorPiece(c, lambda, k, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in) :: k, c(2)
  real (dp), intent(in) :: lambda
  real(dp), intent(out) :: res
  type (Model)          :: Mod

  Mod = Model(lambda, [0._dp, 0._dp], c, 'piece');  res = Mod%Taylor(k)

end subroutine f90TaylorPiece

!ccccccccccccccc

subroutine f90MomentModel(c, clen, lambda, k, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%MomentModel(k)

end subroutine f90MomentModel

!ccccccccccccccc

subroutine f90ModelPiece(c, lambda, k, l, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in) :: k, c(2)
  real (dp), intent(in) :: lambda, l
  real(dp), intent(out) :: res
  type (Model)          :: Mod

  Mod = Model(lambda, [0._dp, 0._dp], c, 'piece');  res = Mod%ShapeFun(k, l)

end subroutine f90ModelPiece

!ccccccccccccccc

subroutine f90MSbarDeltaPiece(nl, nh, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                    , intent(in )  :: nl, nh
  real (dp), dimension(0:4,4), intent(out)  :: beta

  beta = MSbarDeltaPiece(nl, nh)

end subroutine f90MSbarDeltaPiece

!ccccccccccccccc

subroutine f90AlphaMatchingLog(str, direction, nf, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)          , intent(in )  :: str, direction
  integer                      , intent(in )  :: nf
  real (dp), dimension(0:4,5), intent(out)    :: beta
  type (AnomDim)                              :: run

  run = AnomDim(str, nf, 0._dp)
  beta = run%MatchingAlphaLog(direction)

end subroutine f90AlphaMatchingLog

!ccccccccccccccc

subroutine f90AlphaMatchingInverse(str, nf, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)      , intent(in )  :: str
  integer                  , intent(in )  :: nf
  real (dp), dimension(5), intent(out)    :: beta
  type (AnomDim)                          :: run

  run = AnomDim(str, nf, 0._dp)
  beta = getInverse( run%Matchingalpha() )

end subroutine f90AlphaMatchingInverse

!ccccccccccccccc

subroutine f90AnomDim(str, nf, G4, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in )  :: str
  integer                , intent(in )  :: nf
  real (dp)              , intent(in )  :: G4
  real (dp), dimension(5), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, G4);  beta = run%betaQCD(str)

end subroutine f90AnomDim

!ccccccccccccccc

subroutine f90cCoef(nf, order, n, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                  , intent(in ) :: nf, order, n
  real (dp), dimension(0:n), intent(out) :: beta
  type (AnomDim)                         :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%cCoeff(order, n)

end subroutine f90cCoef

!ccccccccccccccc

subroutine f90PScoef(nf, lg, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                , intent(in )  :: nf
  real (dp)              , intent(in )  :: lg
  real (dp), dimension(4), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%PScoef(lg)

end subroutine f90PScoef

!ccccccccccccccc

subroutine f90N12(str, order, nf, lambda, err, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: nf, order
  real (dp)          , intent(in ) :: lambda, err
  real (dp)          , intent(out) :: res
  type (AnomDim)                   :: run

  run = AnomDim('MSbar', nf, 0._dp, err);  res = run%N12(order, str, lambda)

end subroutine f90N12

!ccccccccccccccc

subroutine f90N12Generic(a, order, nf, lambda, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  real (dp), dimension(4), intent(in ) :: a
  integer                , intent(in ) :: nf, order
  real (dp)              , intent(in ) :: lambda
  real (dp)              , intent(out) :: res
  type (AnomDim)                       :: run

  run = AnomDim('MSbar', nf, 0._dp, 0._dp);  res = run%N12Generic(order, a, lambda)

end subroutine f90N12Generic

!ccccccccccccccc

subroutine f90Scoef(str, nf, beta)
  use AlphaClass;  use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass; implicit none

  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp), dimension(3), intent(out) :: beta
  real (dp)                            :: mu = 20._dp
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, &
  0._dp, 0._dp, 0._dp, 0._dp, 0._dp)

  MatEl     = MatrixElementsMass(alphaAll, nf, 0, 0._dp, 0._dp, 0._dp, mu, &
  0._dp, mu, mu, tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)

  beta = MatEl%Scoef(str)

end subroutine f90Scoef

!ccccccccccccccc

subroutine f90sCoefLambda(str, nf, lambda, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp)              , intent(in ) :: lambda
  real (dp), dimension(4), intent(out) :: beta
  type (AnomDim)                       :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%sCoefLambda( str, lambda )

end subroutine f90sCoefLambda

!ccccccccccccccc

subroutine f90sCoefGamma(gama, n, nf, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                , intent(in ) :: nf, n
  real (dp), dimension(n), intent(in ) :: gama
  real (dp), dimension(n), intent(out) :: beta
  type (AnomDim)                       :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%sCoef( gama )

end subroutine f90sCoefGamma

!ccccccccccccccc

subroutine f90wTilde(order, nf, gamma, a0, a1, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer  , intent(in )                 :: order, nf
  real (dp), intent(in )                 :: a0, a1
  real (dp), intent(in ), dimension(0:3) :: gamma
  real (dp), intent(out)                 :: res
  type (AnomDim)                         :: run

  run = AnomDim('MSbar', nf, 0._dp);  res = run%wTilde(order, gamma, a0, a1)

end subroutine f90wTilde

!ccccccccccccccc

subroutine f90kTilde(order, nf, gamma, a0, a1, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                 , intent(in) :: order, nf
  real (dp)               , intent(in ) :: a0, a1
  real (dp)               , intent(out) :: res
  real (dp), intent(in), dimension(0:3) :: gamma
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  res = run%kTilde(order, gamma, a0, a1)

end subroutine f90kTilde

!ccccccccccccccc

subroutine f90InteCorre(b, x0, x1, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  real (dp), intent(in ) :: b, x0, x1
  real (dp), intent(out) :: res

  res = InteCorre(b, x0, x1)

end subroutine f90InteCorre

!ccccccccccccccc

subroutine f90alphaQCD(str, method, order, run, nf, mZ, amZ, mT, muT, mB, muB, &
mC, muC, mu, res)
  use AlphaClass;  use constants, only: dp; use AnomDimClass; implicit none
  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaOb
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaOb = alpha( Andim, order, run, mZ, amZ, mT, muT, mB, &
  muB, mC, muC, method(:9) )

  res     = alphaOb%alphaQCD(nf, mu)

end subroutine f90alphaQCD

!ccccccccccccccc

subroutine f90alphaComplex(str, method, order, run, nf, mZ, amZ, mT, muT, mB, muB, &
mC, muC, muR, muI, res)
  use AlphaClass;  use constants, only: dp; use AnomDimClass; implicit none
  character (len = *)    , intent(in ) :: str, method
  integer                , intent(in ) :: order, run, nf
  real (dp)              , intent(in ) :: mZ, amZ, muR, muI, mT, muT, mB, muB, mC, muC
  real (dp), dimension(2), intent(out) :: res
  type (Alpha)                         :: alphaOb
  type (AnomDim), dimension(3:6)       :: AnDim
  complex (dp)                         :: alp
  integer                              :: i

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaOb = alpha( AnDim, order, run, mZ, amZ, mT, muT, mB, &
  muB, mC, muC, method(:9) )

  alp = alphaOb%alphaQCD( nf, cmplx(muR, muI, dp) )

  res = [ RealPart(alp), ImagPart(alp) ]

end subroutine f90alphaComplex

!ccccccccccccccc

subroutine f90MSbarMass(orderAlpha, runAlpha, run, nf, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass; implicit none

  integer          , intent(in ) :: orderAlpha, runAlpha, run, nf
  real (dp)        , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)        , intent(out) :: res
  type (Running)                 :: alphaMass
  type (Alpha)                   :: alphaAll
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  alphaMass = Running(nf, run, alphaAll, muC);  res = alphaMass%MSbarMass(mu)

end subroutine f90MSbarMass

!ccccccccccccccc

subroutine f90PoleMass(orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer            , intent(in ) :: orderAlpha, runAlpha, run, order, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  alphaMass = Running(nf, run, alphaAll, muC);  res = alphaMass%PoleMass(order, mu)

end subroutine f90PoleMass

!ccccccccccccccc

subroutine f90MSbarMassLow(orderAlpha, runAlpha, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer          , intent(in ) :: orderAlpha, runAlpha, run, nf
  real (dp)        , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)        , intent(out) :: res
  type (Running)                 :: alphaMass
  type (Alpha)                   :: alphaAll
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  alphaMass = Running(nf, run, alphaAll, muC); res = alphaMass%MSbarMassLow(mu)

end subroutine f90MSbarMassLow

!ccccccccccccccc

subroutine f90OptimalR(type, n, method, orderAlpha, runAlpha, order, run, nf, &
mZ, amZ, mT, muT, mB, muB, mC, muC, lambda, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in ) :: method, type
  integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, lambda, n
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  alphaMass = Running(nf, run, alphaAll, mu)
  res = alphaMass%OptimalR( type, n, order, lambda, method(:8) )

end subroutine f90OptimalR

!ccccccccccccccc

  subroutine f90UpsilonDeltaCharm(n, l, alp, mb, mc, res)
    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none
    real (dp)        , intent(in ) :: alp, mb, mc
    integer          , intent(in ) :: n, l
    real (dp)        , intent(out) :: res
    type (NRQCD)                   :: Upsilon
    type (VFNSMSR)                 :: MSR
    type (Alpha)                   :: alphaAll
    type (Running), dimension(2)   :: alphaMass
    type (AnomDim), dimension(3:6) :: AnDim
    integer                        :: i

    do i = 3, 6
      AnDim(i) = AnomDim('MSbar', i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, &
    0._dp, 0._dp, mc, mc, mc, mc)

    alphaMass = Running(4, 0, alphaAll, 100._dp)
    MSR       = VFNSMSR(alphaMass)
    Upsilon   = NRQCD( 'up', 'MSbar', 'no', MSR, n, l, 0, 0 )

    res = Upsilon%DeltaCharm(alp, mb)

  end subroutine f90UpsilonDeltaCharm

!ccccccccccccccc

  subroutine f90UpsilonDeltaCharmBin(n, l, alp, mb, mc, res)
    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none
    real (dp)        , intent(in ) :: alp, mb, mc
    integer          , intent(in ) :: n, l
    real (dp)        , intent(out) :: res
    type (NRQCD)                   :: Upsilon
    type (VFNSMSR)                 :: MSR
    type (Alpha)                   :: alphaAll
    type (Running), dimension(2)   :: alphaMass
    type (AnomDim), dimension(3:6) :: AnDim
    integer                        :: i

    do i = 3, 6
      AnDim(i) = AnomDim('MSbar', i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, &
    0._dp, 0._dp, mc, mc, mc, mc)

    alphaMass = Running(4, 0, alphaAll, 100._dp)
    MSR       = VFNSMSR(alphaMass)
    Upsilon   = NRQCD( 'up', 'MSbar', 'no', MSR, n, l, 0, 0 )

    res = Upsilon%DeltaCharmBin(alp, mb)

  end subroutine f90UpsilonDeltaCharmBin

!ccccccccccccccc

  subroutine f90DeltaCharmExact(charm, type, scheme, average, n, l, j, s, nl, &
    mH, mL, mu, alp, res)
    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *), intent(in ) :: type, scheme, average, charm
    integer            , intent(in ) :: nl, n, l, j, s
    real (dp)          , intent(in ) :: mH, mL, alp, mu
    real (dp)          , intent(out) :: res
    real (dp)                        :: mB, mC, mT
    type (NRQCD)                     :: Upsilon
    type (Alpha)                     :: alphaAll
    type (VFNSMSR)                   :: MSR
    type (Running), dimension(2)     :: alphaMass
    type (AnomDim), dimension(3:6)   :: AnDim
    integer                          :: i

    if (nl == 5) then
      mT = mH; mB = mL
    else if (nl == 4) then
      mB = mH; mC = mL
    else
      res = 0; return
    end if

    do i = 3, 6
      AnDim(i) = AnomDim(scheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, 4, 4, 91.187_dp, 0.118_dp, &
    mT, mT, mB, mB, mC, mC)

    alphaMass = [ Running(nl - 1, 0, alphaAll, mL), &
    Running(nl, 0, alphaAll, mH) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm, scheme, average(:3), MSR, n, l, j, s )

    res = Upsilon%DeltaCharmExact(type, mu, alp, mH)

  end subroutine f90DeltaCharmExact

!ccccccccccccccc

  subroutine f90NRQCD(n, l, j, s, charm, scheme, average, method, orderAlpha, &
  runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1,   &
  lambda2, lam, mu, R, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *), intent(in ) :: method, scheme, charm, average
    integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nl, n,&
    l, j, s
    real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, &
    lambda1, lambda2, lam, R
    real (dp)          , intent(out) :: res(5)
    character (len = 5)              :: alphaScheme
    type (NRQCD)                     :: Upsilon
    type (Alpha)                     :: alphaAll
    type (VFNSMSR)                   :: MSR
    type (Running), dimension(2)     :: alphaMass
    type (AnomDim), dimension(3:6)   :: AnDim
    integer                          :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%En( charm(:6), order, mu, R, lam, method(:8) )

  end subroutine f90NRQCD

!ccccccccccccccc

  subroutine f90FindMass(ord, n, l, j, s, iter, charm, scheme, average, method, &
  orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass,&
  lambda1, lambda2, lam, mu, R, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *), intent(in ) :: method, scheme, charm, iter, average
    integer            , intent(in ) :: ord, orderAlpha, runAlpha, order, run, &
    nl, n, l, j, s
    real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
    lambda1, lam, R, mass, lambda2
    real (dp)          , intent(out) :: res
    character (len = 5)              :: alphaScheme
    type (NRQCD)                     :: Upsilon
    type (VFNSMSR)                   :: MSR
    type (Alpha)                     :: alphaAll
    type (Running), dimension(2)     :: alphaMass
    type (AnomDim), dimension(3:6)   :: AnDim
    integer                          :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%MassFitter( iter(:10), charm(:6), ord, order, mu, R, mass, &
    lam, method(:8) )

  end subroutine f90FindMass

!ccccccccccccccc

  subroutine f90MassError(ord, n, l, j, s, iter, charm, scheme, average, method,&
  orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass,&
  lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, x, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *)    , intent(in ) :: method, scheme, charm, average, iter
    integer                , intent(in ) :: ord, orderAlpha, runAlpha, order, run, &
    nl, n, l, j, s
    real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, x, &
    lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2, muC
    real (dp), dimension(2), intent(out) :: res
    character (len = 5)                  :: alphaScheme
    type (NRQCD)                         :: Upsilon
    type (VFNSMSR)                       :: MSR
    type (Alpha)                         :: alphaAll
    type (Running), dimension(2)         :: alphaMass
    type (AnomDim), dimension(3:6)       :: AnDim
    integer                              :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%MassError( iter(:10), charm(:6), ord, order, mu0, mu1, &
    deltaMu, R0, R1, deltaR, x, mass, lam, method(:8) )

  end subroutine f90MassError

!ccccccccccccccc

  subroutine f90MassList(ord, n, l, j, s, iter, charm, scheme, average, method, &
  orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass,&
  lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *)    , intent(in ) :: method, scheme, charm, average, iter
    integer                , intent(in ) :: ord, orderAlpha, runAlpha, order, run, &
    nl, n, l, j, s
    real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, &
    lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2, muC

    real (dp), dimension( 3, 0:Floor( (mu1 - mu0)/deltaMu ), &
    0:Floor( (R1 - R0)/deltaR ) ), intent(out) :: res

    character (len = 5)                  :: alphaScheme
    type (NRQCD)                         :: Upsilon
    type (VFNSMSR)                       :: MSR
    type (Alpha)                         :: alphaAll
    type (Running), dimension(2)         :: alphaMass
    type (AnomDim), dimension(3:6)       :: AnDim
    integer                              :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%MassList( iter(:10), charm(:6), ord, order, mu0, mu1, &
    deltaMu, R0, R1, deltaR, mass, lam, method(:8) )

  end subroutine f90MassList

!ccccccccccccccc

  subroutine f90NRQCDList(n, l, j, s, iter, charm, scheme, average, method, &
  orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass,&
  lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *)    , intent(in ) :: method, scheme, charm, average, iter
    integer                , intent(in ) :: orderAlpha, runAlpha, order, run, &
    nl, n, l, j, s
    real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, &
    lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2, muC

    real (dp), dimension( 7, 0:Floor( (mu1 - mu0)/deltaMu ), &
    0:Floor( (R1 - R0)/deltaR ) ), intent(out) :: res

    character (len = 5)                  :: alphaScheme
    type (NRQCD)                         :: Upsilon
    type (VFNSMSR)                       :: MSR
    type (Alpha)                         :: alphaAll
    type (Running), dimension(2)         :: alphaMass
    type (AnomDim), dimension(3:6)       :: AnDim
    integer                              :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%NRQCDList( iter(:10), charm(:6), order, mu0, mu1, &
    deltaMu, R0, R1, deltaR, mass, lam, method(:8) )

  end subroutine f90NRQCDList

!ccccccccccccccc

  subroutine f90NRQCDError(n, l, j, s, iter, charm, scheme, average, method, &
  orderAlpha, runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
  mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, x, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *), intent(in ) :: method, scheme, charm, average, iter
    integer            , intent(in ) :: orderAlpha, runAlpha, order, run, &
    nl, n, l, j, s
    real (dp)          , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, muC, x, &
    lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2
    real (dp), dimension(2,0:4), intent(out) :: res
    character (len = 5)              :: alphaScheme
    type (NRQCD)                     :: Upsilon
    type (VFNSMSR)                   :: MSR
    type (Alpha)                     :: alphaAll
    type (Running), dimension(2)     :: alphaMass
    type (AnomDim), dimension(3:6)   :: AnDim
    integer                          :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%EnError( iter(:10), charm(:6), order, mu0, mu1, &
    deltaMu, R0, R1, deltaR, x, mass, lam, method(:8) )

  end subroutine f90NRQCDError

!ccccccccccccccc

  subroutine f90MassIter(n, l, j, s, charm, scheme, average, method, orderAlpha, &
  runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, &
  lambda2, lam, mu, R, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *), intent(in ) :: method, scheme, average, charm
    integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nl, &
    n, l, j, s
    real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
    lambda1, lambda2, lam, R, mass
    real (dp)          , intent(out) :: res(5)
    character (len = 5)              :: alphaScheme
    type (VFNSMSR)                   :: MSR
    type (NRQCD)                     :: Upsilon
    type (Alpha)                     :: alphaAll
    type (Running), dimension(2)     :: alphaMass
    type (AnomDim), dimension(3:6)   :: AnDim
    integer                          :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%MassIter( charm(:6), order, mu, R, mass, lam, method(:8) )

  end subroutine f90MassIter

!ccccccccccccccc

  subroutine f90MassExpand(n, l, j, s, charm, scheme, average, method, orderAlpha, &
  runAlpha, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass, lambda1, &
  lambda2, lam, mu, R, res)

    use RunningClass;  use AlphaClass;  use constants, only: dp
    use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

    character (len = *), intent(in ) :: method, scheme, average, charm
    integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nl, &
    n, l, j, s
    real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
    lambda1, lambda2, lam, R, mass
    real (dp)          , intent(out) :: res(5)
    character (len = 5)              :: alphaScheme
    type (VFNSMSR)                   :: MSR
    type (NRQCD)                     :: Upsilon
    type (Alpha)                     :: alphaAll
    type (Running), dimension(2)     :: alphaMass
    type (AnomDim), dimension(3:6)   :: AnDim
    integer                          :: i

    alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

    do i = 3, 6
      AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
    end do

    alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
    mT, muT, mB, muB, mC, muC)

    alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
    Running(nl, run, alphaAll, lambda1) ]

    MSR     = VFNSMSR(alphaMass)
    Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

    res = Upsilon%EnExpand( charm(:6), order, mu, R, mass, lam, method(:8) )

  end subroutine f90MassExpand

!ccccccccccccccc

subroutine f90mmfromMSR(type, orderAlpha, runAlpha, order, run, nf, mZ, amZ, &
mT, muT, mB, muB, mC, muC, mu, R, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer            , intent(in ) :: orderAlpha, runAlpha, order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, R
  character (len = *), intent(in ) :: type
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  type (Alpha)                     :: alphaAll
  real (dp)                        :: mass
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  if (nf == 5) mass = mT; if (nf == 4) mass = mB; if (nf == 3) mass = mC

  alphaMass = Running(nf, run, alphaAll, mu)
  res = alphaMass%mmFromMSR(type, mass, order, R)

end subroutine f90mmfromMSR

!ccccccccccccccc

subroutine f90MSRMass(type, method, orderAlpha, runAlpha, order, run, nf, mZ, &
amZ, mT, muT, mB, muB, mC, muC, lambda, mu, R, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in) :: method, type
  integer           , intent(in ) :: orderAlpha, runAlpha, order, run, nf
  real (dp)         , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, R, lambda
  real (dp)         , intent(out) :: res
  type (Running)                  :: alphaMass
  type (Alpha)                    :: alphaAll
  type (AnomDim), dimension(3:6)  :: AnDim
  integer                         :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  alphaMass = Running(nf, run, alphaAll, mu)
  res       = alphaMass%MSRMass(type, order, R, lambda, method)

end subroutine f90MSRMass

!ccccccccccccccc

subroutine f90MSRVFNS(up, type, method, orderAlpha, runAlpha, order, run, nf, mZ, &
amZ, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, R, res)
  use RunningClass;  use AlphaClass   ;  use constants, only: dp
  use AnomDimClass;  use VFNSMSRClass ;  implicit none

  character (len = *), intent(in) :: method, type, up
  integer           , intent(in ) :: orderAlpha, runAlpha, order, run, nf
  real (dp)         , intent(in ) :: mZ, amZ, mu1, mu2, mT, muT, mB, muB, mC, &
  muC, R, lambda
  real (dp)         , intent(out) :: res
  type (Running), dimension(2)    :: alphaMass
  type (AnomDim), dimension(3:6)  :: AnDim
  type (Alpha)                    :: alphaAll
  type (VFNSMSR)                  :: MSRVFNS
  integer                         :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  alphaMass(1) = Running(nf - 1, run, alphaAll, mu2)
  alphaMass(2) = Running(nf    , run, alphaAll, mu1)
  MSRVFNS      = VFNSMSR( alphaMass )
  res          = MSRVFNS%MSRMass( up(:4), type(:4), order, R, lambda, method(:8) )

end subroutine f90MSRVFNS

!ccccccccccccccc

subroutine f90JetMass(orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, muLambda, R, mu, res)
 use AlphaClass;    use MatrixElementsClass;  use constants, only: dp
 use AnomDimClass;  implicit none

  integer  , intent(in )      :: orderAlpha, order, runAlpha, run, nf
  real (dp), intent(in )      :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, R
  real (dp), intent(out)      :: res
  type (Alpha)                :: alphaAll
  type (MatricesElementsMass) :: MatEl
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                     :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  MatEl     = MatricesElementsMass(alphaAll, nf, run, 0._dp, 0._dp, 0._dp, 0._dp,  &
  muLambda, muLambda)

  res = MatEl%JetMass(order, run, R, mu)

end subroutine f90JetMass

!ccccccccccccccc

subroutine f90mmFromJetMass(orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, &
muT, mB, muB, mC, muC, muLambda, R, mu, res)
 use AlphaClass;   use MatrixElementsClass;  use constants, only: dp
 use AnomDimClass; implicit none

  integer  , intent(in )         :: orderAlpha, order, runAlpha, run, nf
  real (dp), intent(in )         :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, R
  real (dp), intent(out)         :: res
  type (Alpha)                   :: alphaAll
  type (MatricesElementsMass)    :: MatEl
  real (dp)                      :: mass
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  MatEl     = MatricesElementsMass(alphaAll, nf, run, 0._dp, 0._dp, 0._dp, 0._dp,  &
  muLambda, muLambda)

  if (nf == 6) mass = mT; if (nf == 5) mass = mB; if (nf == 4) mass = mC

  res = MatEl%mmFromJetMass(order, run, R, mu, mass)

end subroutine f90mmFromJetMass

!ccccccccccccccc

subroutine f90Singular(hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, &
 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, &
 mu, c, clen, lambda, R0, mu0, delta0, h, tau, res)
  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muC, j3
  real (dp)          , intent(inout) :: muH, muJ, muS, R
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  integer                            :: i
  type (AnomDim), dimension(3:6)     :: AnDim

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')

    call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
    call Sing%SetRunning(muJ, muS, R, mu)

  res       = Sing%SingleSing(Mod, setup(:15), gap(:12), space(:6), cum(:4), &
  order, R0, mu0, delta0, h, tau)

end subroutine f90Singular

!ccccccccccccccc

subroutine f90SingularDiff(hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, &
 nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, &
 c, clen, lambda, R0, mu0, delta0, h, tau1, tau2, res)
  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, hard, gap
  integer            , intent(in)    :: orderAlpha, order, runAlpha, run, nf, clen
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC,&
  mu,  Q, G3, lambda, tau1, tau2, R0, mu0, delta0, h, s3, j3, muC
  real (dp)          , intent(inout) :: muH, muJ, muS, R
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 3, 6
    AnDim  = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')

  call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
  call Sing%SetRunning(muJ, muS, R, mu)

  res       = Sing%SingleSing(Mod, setup(:15), gap(:12), space(:6), cum(:4), &
  order, R0, mu0, delta0, h, tau1, tau2)

end subroutine f90SingularDiff

!ccccccccccccccc

subroutine f90SingularHJM(hard, setup, gap, space, cum, orderAlpha, runAlpha, order, run, &
  isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
  muH, muJ, muS, R, mu, c, clen, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass; use MatrixElementsClass; use ModelClass; use SingularClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in) :: cum, setup, space, hard, gap
  integer            , intent(in) :: orderAlpha, order, runAlpha, run, nf, clen, isoft
  real (dp)      ,  intent(inout) :: muH, muJ, muS, R
  real (dp)          , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC,&
  muC, j3, Q, G3, lambda, tau, R0, mu0, delta0, h, s31, s32, s3

  real (dp)          ,   intent(out)  :: res
  type (Alpha)                        :: alphaAll
  type (MatrixElements)               :: MatEl
  type (SingularMassless)             :: Sing
  type (Model), dimension(clen, clen) :: Mod
  type (AnomDim), dimension(3:6)      :: AnDim
  integer                             :: i, j

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = Sing%HJMSing(c, Mod, 'sum', setup(:15), gap(:12), space(:3), cum(:4), order, &
                     isoft, s31, s32, R0, mu0, delta0, h, tau)

end subroutine f90SingularHJM

!ccccccccccccccc

subroutine f90SingularHJM1D(hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, &
 j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, &
 muS, R, mu, c, clen, lambda, R0, mu0, delta0, h, tau, res)
  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), intent(in), dimension (clen) :: c
  character (len = *), intent(in ) :: cum, hard, gap
  integer            , intent(in ) :: orderAlpha, order, runAlpha, run, nf, clen, isoft
  real (dp)          , intent(in ) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC , &
                                     s3, muH, Q, muJ, muS, R, G3, lambda, tau, R0, mu0 , &
                                     delta0, h, s31, s32, j3
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaAll
  type (MatrixElements)            :: MatEl
  type (SingularMassless)          :: Sing
  type (Model)                     :: Mod
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')

  res       = Sing%HJMSing1D(Mod, gap(:12), cum(:4), order, isoft, s31, s32, R0, mu0, &
                             delta0, h, tau)

end subroutine f90SingularHJM1D

!ccccccccccccccc

subroutine f90SingularHJM1DPiece(hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, &
 j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS , &
 R, mu, c, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in ) :: cum, hard, gap
  integer              , intent(in ) :: orderAlpha, order, runAlpha, run, nf, isoft
  integer, dimension(2), intent(in ) :: c
  real (dp)            , intent(in ) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC , &
                                       s3, muH, Q, muJ, muS, R, G3, lambda, tau, mu0, &
                                       delta0, h, s31, s32, j3, muC, R0
  real (dp)           , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatrixElements)             :: MatEl
  type (SingularMassless)           :: Sing
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  integer                           :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))
  Mod       = Model(lambda, [1._dp], c, 'piece')

  res       = Sing%HJMSing1D(Mod, gap(:12), cum(:4), order, isoft, s31, s32, R0, mu0, &
                             delta0, h, tau)

end subroutine f90SingularHJM1DPiece

!ccccccccccccccc

subroutine f90SingularDouble(hard, setup, gap, space, cum1, cum2, orderAlpha, runAlpha, order, &
run, isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,   &
 muH, muJ, muS, R, mu, c, clen, lambda, R0, mu0, delta0, h, rho1, rho2, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in) :: hard, cum1, cum2, setup, space, gap
  integer            , intent(in) :: orderAlpha, order, runAlpha, run, nf, clen, isoft
  real (dp)          , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC , &
                                     j3, s3, h, Q, G3, lambda, rho1, rho2, R0, mu0, s31,&
                                     s32, delta0
  real (dp)       , intent(inout) :: muH, muJ, muS, R
  real (dp)       , intent(out)   :: res
  type (Alpha)                        :: alphaAll
  type (MatrixElements)               :: MatEl
  type (SingularMassless)             :: Sing
  type (AnomDim), dimension(3:6)      :: AnDim
  type (Model), dimension(clen, clen) :: Mod
  integer                             :: i, j

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = Sing%DoubleSing(c, Mod, 'sum', setup(:15), gap(:12), space(:3), cum1(:4), &
                        cum2(:4), order, isoft, s31, s32, R0, mu0, delta0, h, rho1, rho2)

end subroutine f90SingularDouble

!ccccccccccccccc

subroutine f90SingularHJMPiece(hard, gap, space, cum, orderAlpha, runAlpha, order, run , &
 isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
 muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau, res)

 use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
 use constants, only: dp; use AnomDimClass; implicit none

 character (len = *)  , intent(in) :: cum, space, hard, gap
 integer              , intent(in) :: orderAlpha, order, runAlpha, run, nf, isoft
 integer, dimension(4), intent(in) :: c
 real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                     j3, s3, muH, Q, muJ, muS, R, G3, lambda, tau, R0 , &
                                     mu0, delta0, h, s31, s32

 real (dp)     , intent(out) :: res

  type (Alpha)                     :: alphaAll
  type (MatrixElements)            :: MatEl
  type (SingularMassless)          :: Sing
  type (AnomDim), dimension(3:6)   :: AnDim
  type (Model)   , dimension(2,2)  :: ModList
  real (dp)      , dimension(2,2)  :: clist = 1
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  ModList(1,1) = Model(lambda, [1._dp], c(:2), 'piece'); ModList(2,2) = ModList(1,1)
  ModList(1,2) = Model(lambda, [1._dp], c(3:), 'piece'); ModList(2,1) = ModList(1,2)

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  res = Sing%HJMSing(clist, ModList, 'piece', 'Model', &
  gap(:12), space(:3), cum(:4), order, isoft, s31, s32, R0, mu0, delta0, h, tau)

end subroutine f90SingularHJMPiece

!ccccccccccccccc

subroutine f90SingularDoublePiece(hard, gap, space, cum1, cum2, orderAlpha, runAlpha, order , &
 run, isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
 muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, rho1, rho2, res)

 use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
 use constants, only: dp; use AnomDimClass; implicit none

 integer, dimension(4), intent(in ) :: c
 character (len = *)  , intent(in ) :: cum1, cum2, space, hard, gap
 integer              , intent(in ) :: orderAlpha, order, runAlpha, run, nf, isoft
 real (dp)            , intent(in ) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                       j3, s3, muH, muJ, muS, R, G3, lambda, rho1, rho2, &
                                       R0, mu0, Q, delta0, h, s31, s32
 real (dp)            , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatrixElements)             :: MatEl
  type (SingularMassless)           :: Sing
  type (AnomDim) , dimension(3:6)   :: AnDim
  type (Model)   , dimension(2,2)   :: ModList
  real (dp)      , dimension(2,2)   :: clist = 1
  integer                           :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  ModList(1,1) = Model(lambda, [1._dp], c(:2), 'piece'); ModList(2,2) = ModList(1,1)
  ModList(1,2) = Model(lambda, [1._dp], c(3:), 'piece'); ModList(2,1) = ModList(1,2)

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  res = Sing%DoubleSing(clist, modList, 'piece'  , &
                        'Model', gap(:12), space(:3), cum1(:4), cum2(:4), order, isoft, &
                         s31, s32, R0, mu0, delta0, h, rho1, rho2)

end subroutine f90SingularDoublePiece

!ccccccccccccccc

subroutine f90SingularPiece(hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, &
j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, piece, &
lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in) :: shape, cum, space, hard, gap
  integer              , intent(in) :: orderAlpha, order, runAlpha, run, nf
  integer, dimension(2), intent(in) :: piece
  real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, &
  Q, muJ, muS, R, G3, lambda, tau, R0, mu0, delta0, j3, s3, muH, h, mC, muC
  real (dp)           , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatricesElements)           :: MatEl
  type (SingularScales)             :: Sing
  type (Model)                      :: Mod
  type (AnomDim),  dimension(3:6)   :: AnDim
  integer                           :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Mod       = Model(lambda, [1._dp], piece, 'piece')
  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
 mB, muB, mC, muC)

  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )

  call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
  call Sing%SetRunning(muJ, muS, R, mu)

  res = Sing%SingleSing(Mod, 'Model', gap(:12), space(:3), cum(:4), order, R0, &
  mu0, delta0, h, tau)

end subroutine f90SingularPiece

!ccccccccccccccc

subroutine f90SingularList(hard, shape, gap, space, cum, orderAlpha, runAlpha, order  , &
run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, &
mu, clen, lambda, R0, mu0, delta0, h, tau, res)

 use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
 use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in) :: shape, cum, space, hard, gap
  integer              , intent(in) :: orderAlpha, order, runAlpha, run, nf, clen
  real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                      Q, muJ, muS, R, G3, lambda, tau, R0, mu0, delta0,  &
                                      j3, s3, muH, h
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  integer                          :: i, j, k
  type (AnomDim), dimension(3:6)   :: AnDim
  type (Alpha)                     :: alphaAll
  type (MatrixElements)            :: MatEl
  type (SingularMassless)          :: Sing
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless( MatEl, run, mu, shape(:6), hard(:6) )

  res = Sing%SingleSing(Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, delta0, h, tau)

end subroutine f90SingularList

!ccccccccccccccc

subroutine f90SingularDiffPiece(hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run,&
nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu,   &
piece, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in) :: shape, cum, space, hard, gap
  integer              , intent(in) :: orderAlpha, order, runAlpha, run, nf
  integer, dimension(2), intent(in) :: piece
  real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                      Q, muJ, muS, R, G3, lambda, tau, R0, mu0, delta0, &
                                      j3, s3, muH, h, tau2
  real (dp)           , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatricesElements)           :: MatEl
  type (SingularScales)             :: Sing
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  integer                           :: i

 do i = 3, 6
   AnDim(i) = AnomDim('MSbar', i, G3)
 end do

 Mod       = Model(lambda, [1._dp], piece, 'piece')
 alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                   mB, muB, mC, muC)
 MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
 Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )

 call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
 call Sing%SetRunning(muJ, muS, R, mu)

 res = Sing%SingleSing(Mod, 'Model', gap(:12), space(:3), cum(:4), order, R0, mu0, &
                       delta0, h, tau, tau2)

end subroutine f90SingularDiffPiece

!ccccccccccccccc

subroutine f90SingularDiffList(hard, shape, gap, space, cum, orderAlpha, runAlpha, order,&
 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, &
 mu, clen, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in) :: shape, cum, space, hard, gap
  integer            , intent(in) :: orderAlpha, order, runAlpha, run, nf, clen
  real (dp)          , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, &
                                     Q, muJ, muS, R, G3, lambda, tau, R0, delta0, &
                                     j3, s3, muH, h, tau2, mu0, muC
  integer                         :: i, j, k
  type (Alpha)                    :: alphaAll
  type (MatricesElements)         :: MatEl
  type (SingularScales)           :: Sing
  type (AnomDim), dimension(3:6)  :: AnDim
  real (dp)   , dimension( (clen + 2) * (clen + 1)/2 ) , intent(out) :: res
  type (Model), dimension( (clen + 2) * (clen + 1)/2 )               :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i

    Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1

    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC)
  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )

   call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
   call Sing%SetRunning(muJ, muS, R, mu)

  res = Sing%SingleSing(Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, &
                        delta0, h, tau, tau2)

end subroutine f90SingularDiffList

!ccccccccccccccc

subroutine f90Rhad(str, orderAlpha, runAlpha, order, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, mu, Q, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: order, runAlpha, orderAlpha, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, Q
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, mB, muB, mC, muC)
  alphaMass = Running(nf, runAlpha, alphaAll, muC)
  EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)

  MatEl = Sigma(alphaMass, EW);  res = MatEl%Rhad(order, mu, Q)

end subroutine f90Rhad

!ccccccccccccccc

subroutine f90RhadCoefs(nf, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  integer                , intent(in ) :: nf
  real (dp), dimension(4), intent(out) :: res
  type (Running)                       :: alphaMass
  type (Alpha)                         :: alphaAll
  type (Sigma)                         :: MatEl
  type (ElectroWeak)                   :: EW
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, &
  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = Running(nf, 0, alphaAll, 0._dp)
  EW = ElectroWeak(0._dp, 2.4952_dp, 0.23119_dp)

  MatEl = Sigma(alphaMass, EW);  res = MatEl%RhadCoefs()

end subroutine f90RhadCoefs

!ccccccccccccccc

subroutine f90RhadMass(str, curr, orderAlpha, runAlpha, runMass, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, mT, muT, mB, muB, mC, muC, mu, Q, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlpha, orderAlpha, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, gammaZ, &
  sin2ThetaW, Q
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)

  res       = MatEl%RhadMass(curr(:6), order, mu, Q)

end subroutine f90RhadMass

!ccccccccccccccc

subroutine f90lambdaQCD(str, order, runAlpha, run, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, mu, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: order, runAlpha, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim( str(:5), i, 0._dp )
  end do

  alphaAll  = Alpha(AnDim, order, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)
  alphaMass = Running(nf, run, alphaAll, mu)
  res       = alphaMass%lambdaQCD(run)

end subroutine f90lambdaQCD

!ccccccccccccccc

subroutine f90delta(str, nf, mu, R, res)
  use AlphaClass; use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass; implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp)              , intent(in ) :: mu, R
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, &
  0._dp, 0._dp, 0._dp, 0._dp, 0._dp)

  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, &
  mu, mu, mu, tiny(1._dp), mu, R, R, 0._dp, 0._dp)

  res       = MatEl%delta( str(:12) )

end subroutine f90delta

!ccccccccccccccc

subroutine f90PSdelta(orderAlpha, runAlpha, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, R, lg, res)
  use AlphaClass;    use RunningClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer                , intent(in ) :: nf, orderAlpha, runAlpha
  real (dp)              , intent(in ) :: mu, R, mZ, amZ, mT, muT, mB, muB, mC, muC, lg
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (Running)                       :: alphaMass
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  alphaMass = Running(nf, 0, alphaAll, 0._dp)

  res       = alphaMass%PSdelta( R, mu, lg )

end subroutine f90PSdelta

!ccccccccccccccc

subroutine f90deltaGap(str, orderAlpha, runAlpha, runMass, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, R, res)
  use AlphaClass;    use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf, orderAlpha, runAlpha, runMass
  real (dp)              , intent(in ) :: mu, R, mZ, amZ, mT, muT, mB, muB, mC, muC
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  real (dp)                            :: muM
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  muM = tiny(1._dp)

  if ( str(:3) == 'MSR' .or. str(:8) == 'MSbarLow' .or. str(:7) == 'JetMass' ) muM = 2 * mu

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, &
  mT, muT, mB, muB, mC, muC)

  MatEl     = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, &
  0._dp, mu, mu, mu, muM, mu, R, R, 0._dp, 0._dp)

  res       = MatEl%delta( str(:12) )

end subroutine f90deltaGap

!ccccccccccccccc

subroutine f90deltaMSbar(orderAlpha, runAlpha, run, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, mu, res)
  use RunningClass; use AlphaClass;  use MatrixElementsClass
  use constants, only: dp; use AnomDimClass; implicit none
  integer                , intent(in ) :: orderAlpha, runAlpha, run, nf
  real (dp)              , intent(in ) :: mZ, amZ, mu, mT, muT, muB, mB, mC, muC
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlpha, runAlpha, mZ, amZ, mT, muT, &
  mB, muB, mC, muC)

  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, mu, mu, &
  tiny(1._dp), mu, mu, mu, muC, muC)

  res       = MatEl%delta('MSbar')

end subroutine f90deltaMSbar

!ccccccccccccccc

subroutine f90CoefMat(str, nf, s3, res)
  use AlphaClass; use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  character (len = *)        , intent(in ) :: str
  integer                    , intent(in ) :: nf
  real (dp)                  , intent(in ) :: s3
  real (dp), dimension(0:4,3), intent(out) :: res
  real (dp)                                :: mu = 20._dp
  type (Alpha)                             :: alphaAll
  type (MatrixElementsMass)                :: MatEl
  integer                                  :: i
  type (AnomDim), dimension(3:6)           :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                     0._dp, 0._dp, 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, s3, s3, s3, s3, mu, mu, mu, tiny(1._dp),&
                             mu, mu, mu, 0._dp, 0._dp)
  res = MatEl%CoefMat(str)

end subroutine f90CoefMat

!ccccccccccccccc

subroutine f90NGLintegral(nf, pow, w1, w2, res)
  use AlphaClass; use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer     , intent(in ) :: nf, pow
  real (dp)   , intent(in ) :: w1, w2
  real (dp)   , intent(out) :: res
  real (dp)                 :: mu = 20._dp
  type (MatrixElementsMass) :: MatEl
  type (Alpha)              :: alphaAll
  integer                   :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                     0._dp, 0._dp, 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, mu, mu, &
                         tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)

 res = MatEl%NGLIntegral(pow, w1, w2)

end subroutine f90NGLintegral

!ccccccccccccccc

subroutine f90NGLDoubleintegral(nf, pow, w1, w2, ratio, res)
  use MatrixElementsClass;  use constants, only: dp; implicit none
  integer        , intent(in ) :: nf, pow
  real (dp)      , intent(in ) :: w1, w2, ratio
  real (dp)      , intent(out) :: res

  res = NGLDoubleIntegral(nf, pow, w1, w2, ratio)

end subroutine f90NGLDoubleintegral

!ccccccccccccccc

subroutine f90NGLfunction(nf, x, res)
  use AlphaClass;  use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer     , intent(in ) :: nf
  real (dp)   , intent(in ) :: x
  real (dp)   , intent(out) :: res
  real (dp)                 :: mu = 20._dp
  type (Alpha)              :: alphaAll
  type (MatrixElementsMass) :: MatEl
  integer                   :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                     0._dp, 0._dp, 0._dp)
  MatEl    = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, mu, mu, &
                         tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)
  res      = MatEl%NGLfunction(x)

end subroutine f90NGLfunction

!ccccccccccccccc

subroutine f90DiffDeltaGap(str, scheme, order, R0, R, mu0, mu, muLambda, orderAlpha, &
                           runAlpha, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, res)
  use MatrixElementsClass; use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  character (len = *), intent(in ) :: str, scheme
  integer            , intent(in ) :: order, nf, orderAlpha, runAlpha
  real (dp)          , intent(in ) :: R0, R, mu0, mu, muLambda, mZ, amZ, mT, muT, mB, &
                                      muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaAll
  type (MatrixElementsMass)        :: MatEl
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(scheme(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlpha, runAlpha, mZ, amZ, mT, muT, mB, muB, mC, muC)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, &
  mu, mu, tiny(1._dp), mu, mu, mu, muLambda, muLambda)

  res       = MatEl%DiffDeltaGap( str(:12), order, R0, R, mu0, mu )

end subroutine f90DiffDeltaGap

!ccccccccccccccc

subroutine f90kernels(n, width, w, mu, p, res)
  use KernelsClass;  use constants, only: dp; implicit none
  integer  , intent(in ) :: n
  real (dp), intent(in ) :: width, w, mu, p
  real (dp), intent(out) :: res(0:n)
  type (Kernels)         :: ker

  ker = Kernels(n, width, w, mu, p);  res = ker%KernelList()

end subroutine f90kernels

!ccccccccccccccc

subroutine f90GammaDerList(n, w, res)
  use KernelsClass;  use constants, only: dp; implicit none
  integer   , intent(in) :: n
  real (dp), intent(in ) :: w
  real (dp), intent(out) :: res(0:n)
  type (Kernels)         :: ker

 ker = Kernels(n, 0._dp, w, 3._dp, 3._dp)
 res = ker%DerList()

end subroutine f90GammaDerList

!ccccccccccccccc

subroutine f90NGLkernels(n, n1, n2, width, w, mu, p, res)
  use KernelsClass;  use constants, only: dp; implicit none
  integer                , intent(in ) :: n, n1, n2
  real (dp), dimension(2), intent(in ) :: w, p
  real (dp)              , intent(in ) :: mu, width
  real (dp)              , intent(out) :: res(n)
  type (Kernels)                       :: ker1, ker2
  real (dp), dimension(0:2*n+n1)       :: kerList1
  real (dp), dimension(0:2*n+n2)       :: kerList2

  ker1     = Kernels( 2*n + n1, width, w(1), mu, p(1) )
  ker2     = Kernels( 2*n + n2, width, w(2), mu, p(2) )

  kerList1 = ker1%KernelList(); kerList2 = ker2%KernelList()

  res      = KernelsCombine(kerList1, kerList2, n)

end subroutine f90NGLkernels

!ccccccccccccccc

subroutine f90hyper2f1 ( a, b, c, x, res )
  use hyper; use constants, only: dp; implicit none
  real (dp), intent(in)  :: a, b, c, x
  real (dp), intent(out) :: res

  res = hyper2f1(a, b, c, x)

end subroutine f90hyper2f1

!ccccccccccccccc

subroutine f90hyperf32exact(w2, x2, res)
  use hyper; use constants, only: dp; implicit none
  real (dp), intent(in ) :: w2, x2
  real (dp), intent(out) :: res

  res = H3F2Exact(w2, x2)

end subroutine f90hyperf32exact

!ccccccccccccccc

subroutine f90GammaR(str, nf, beta)
  use MatrixElementsClass;  use AlphaClass; use constants
  use AnomDimClass;  implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp), dimension(3), intent(out) :: beta
  real (dp)                            :: mu = 20._dp
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                    0._dp, 0._dp, 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, nf, 0, 0._dp, 0._dp, 0._dp, mu, 0._dp, mu, mu, &
                            tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)

  beta = MatEl%GammaR( str(:7) )

end subroutine f90GammaR

!ccccccccccccccc
