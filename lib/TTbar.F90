
! ccccccccccc

module RNRQCDClass
  use constants, only: dp, Pi, Pi2, l2, Zeta3, Euler, Pi3, prec; use ToppikClass
  use DeriGamma, only: DiGam, trigam ; use RunningClass; use AnomDimClass
  use AlphaClass; use NRQCDClass; use VFNSMSRClass; use SigmaClass
  use ElectroWeakClass; use adapt; use QuadPack, only: qags
  implicit none ;  private

  real (dp)   , parameter :: FPi = 4 * Pi
  complex (dp), parameter :: cI = (0,1)

  public :: qFromV, vC, VcsLL, XiNNLLSoftMixLogc1, vStar, vRootStar, SwitchOff

! ccccccccccc

  type, public                :: RNRQCD
    integer                   :: nl, ordMass, nf
    real (dp)                 :: a1, a2, mass, Qt, gt, lambda, m1S, MSRm, QSwitch
    real (dp), dimension(0:4) :: beta
    type (ElectroWeak)        :: EW
    type (Running)            :: run
    type (Alpha)              :: AlphaOb
    type (AnomDim)            :: andim
    character (len = 4)       :: scheme
    character (len = 8)       :: method
    real (dp), dimension(0:4) :: m1Spiece

  contains

    procedure, pass (self), public :: VceffsNNLL, VrsLL, V2sLL, VssLL, Vk1sLL, &
    Vk2sLL, VkeffsLL, XiNLL, XiNNLLnonmix, XiNNLLmixUsoft, MLLc2, MNLLc1, &
    MNLLplusNNLLnonmixc1, MNNLLAllc1InclSoftMixLog, A1pole, Xsec, Delta1S, &
    MNLLc1Square, MNNLLc1Square, Rexp, RQCD, Rmatched, Switch, RmatchedList, &
    SigmaMatchedRadiative, SigmaMatched, SigmaMatchedRadiativeCone, &
    SigmaMatchedRadiativeCum, SigmaMatchedRadiativeConeCum

    procedure, pass (self), private :: xc01, xc11, xc12, xc22

  end type RNRQCD

!ccccccccccccccc

  interface RNRQCD
    module procedure RNRQCDIn
  end interface RNRQCD

contains

  type (RNRQCD) function RNRQCDIn(MSR, EW, scheme, method, gt, ordMass, ord1S, &
  R1S, lambda)
    type (VFNSMSR)     , intent(in) :: MSR
    type (ElectroWeak) , intent(in) :: EW
    real (dp)          , intent(in) :: gt, lambda, R1S
    integer            , intent(in) :: ordMass, ord1S
    character (len = *), intent(in) :: scheme, method
    type (Running)                  :: run
    type (AnomDim)                  :: adim
    type (NRQCD)                    :: Upsilon
    integer                         :: nl
    real (dp), dimension(4)         :: A

    run = MSR%RunArray(2); RNRQCDIn%gt = gt; RNRQCDIn%ordMass = ordMass
    RNRQCDIn%method = method(:8); RNRQCDIn%lambda = lambda

    RNRQCDIn%scheme = scheme(:4); if ( scheme(:2) == 'S1' ) RNRQCDIn%scheme = 'MSRn'

    nl = run%numFlav(); adim = run%adim(); RNRQCDIn%run = run; RNRQCDIn%nl = nl
    RNRQCDIn%beta = adim%betaQCD('beta') ; RNRQCDIn%AlphaOb = run%alphaAll()
    RNRQCDIn%mass = run%scales('mH')     ; RNRQCDIn%andim = adim

    RNRQCDIn%a1 = 31._dp/3 - 10._dp * nl/9
    RNRQCDIn%a2 = 456.74883699902244_dp - 66.35417150661816_dp * nl + &
    1.2345679012345678_dp * nl**2

    Upsilon = NRQCD('up', 'MSRn', 'no', MSR, 1, 0, 1, 1)

    if (nl == 5 .or. nl == 3 .or. nl == 1) RNRQCDIn%Qt =   2._dp/3
    if (nl == 4 .or. nl == 2 .or. nl == 0) RNRQCDIn%Qt = - 1._dp/3

    RNRQCDIn%nf = nl + 1;  RNRQCDIn%EW = EW

    if ( scheme(:4) == 'pole' ) then

      RNRQCDIn%m1S = RNRQCDIn%mass

    else

      RNRQCDIn%m1Spiece = Upsilon%En(ordMass, R1S, R1S, lambda, method, 'ttbar')/2
      RNRQCDIn%m1S = sum( RNRQCDIn%m1Spiece(:ord1S) )

      RNRQCDIn%MSRm = run%MSRMass( scheme(:4), ordMass, RNRQCDIn%mass, &
      lambda, method)

      A = powList( RNRQCDIn%MSRm/RNRQCDIn%m1S, 4 )

      RNRQCDIn%QSwitch = 2 * RNRQCDIn%m1S + Sqrt(  ( 6.25e-6_dp - A(1)/2000 + &
      3 * A(2)/200 - A(3)/5 + A(4) ) * RNRQCDIn%m1S**2 - gt**2  )

    end if

  end function RNRQCDIn

! ccccccccccc

  function RmatchedList(self, order, h, hnu, v1, v2, Q0, Q1, DeltaQ) result(list)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q0, Q1, DeltaQ, hnu, v1, v2
    real (dp)                          :: nu, Q
    integer                            :: i, imax
    real(dp), dimension( 2,0:Floor( (Q1 - Q0)/deltaQ ) ) :: list

    imax = Floor( (Q1 - Q0)/deltaQ )

    do i = 0, imax
      Q = Q0 + i * deltaQ; list(1,i) = Q; nu = hnu * vStar(Q, self%m1S, self%gt)
      list(2,i) = self%Rmatched(order, h, nu, v1, v2, Q)
    end do

  end function RmatchedList

! ccccccccccc

  real (dp) function Rmatched(self, order, h, nu, v1, v2, Q)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu, v1, v2

    Rmatched = self%RQCD(order, h, Q) + ( self%Xsec(order, Q, h, nu) &
    - self%Rexp(order, h, nu, Q) ) * SwitchOff(Q, self%m1S, self%gt, v1, v2)

  end function Rmatched

! ccccccccccc

  real (dp) function SigmaMatched(self, order, h, nu, v1, v2, Q)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu, v1, v2
    real (dp)        , dimension(2)    :: EWfactors

    EWfactors    = self%EW%EWfactors(self%nf, Q)
    SigmaMatched = FPi * ( self%run%alphaQED(Q)/Q )**2 * EWfactors(1) * &
    self%Rmatched(order, h, nu, v1, v2, Q)

  end function SigmaMatched

! ccccccccccc

  real (dp) function SigmaMatchedRadiative(self, order, h, nu, v1, v2, Q, x, theta)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu, v1, v2, x, theta

    SigmaMatchedRadiative = self%run%alphaQED(Q)/Pi2 * g(x, theta) * &
    self%SigmaMatched( order, h, nu, v1, v2, Q * sqrt(1 - 2 * x) )

  end function SigmaMatchedRadiative

! ccccccccccc

  real (dp) function SigmaMatchedRadiativeCum(self, order, h, nu, v1, v2, Q, &
  x0, x1, theta)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu, v1, v2, x0, x1, theta
    real (dp)                          :: abserr

    call DAdapt(integrand, x0, x1, 1, prec, prec, SigmaMatchedRadiativeCum, &
    abserr)

  contains

    real (dp) function integrand(x)
      real (dp), intent(in) :: x

      integrand = self%SigmaMatchedRadiative(order, h, nu, v1, v2, Q, x, theta)

    end function integrand

  end function SigmaMatchedRadiativeCum

! ccccccccccc

  real (dp) function SigmaMatchedRadiativeConeCum(self, order, h, nu, v1, v2, Q, &
  x0, x1, theta, deltaTheta)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu, v1, v2, x0, x1, theta, deltaTheta
    real (dp)                          :: abserr

    call DAdapt(integrand, x0, x1, 1, prec, prec, SigmaMatchedRadiativeConeCum,&
    abserr)

  contains

    real (dp) function integrand(x)
      real (dp), intent(in) :: x

      integrand = self%SigmaMatchedRadiativeCone(order, h, nu, v1, v2, Q, x, &
      theta, deltaTheta)

    end function integrand

  end function SigmaMatchedRadiativeConeCum

! ccccccccccc

  real (dp) function SigmaMatchedRadiativeCone(self, order, h, nu, v1, v2, Q, &
  x, theta, deltaTheta)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu, v1, v2, x, theta, deltaTheta

    SigmaMatchedRadiativeCone = self%run%alphaQED(Q) * gInt(x, theta, deltaTheta) * &
    self%SigmaMatched( order, h, nu, v1, v2, Q * sqrt(1 - 2 * x) )/Pi2

  end function SigmaMatchedRadiativeCone

! ccccccccccc

  real (dp) function Switch(self)
    class (RNRQCD), intent(in) :: self

    Switch = self%Qswitch

  end function Switch

! ccccccccccc

  real (dp) function Rexp(self, order, h, nu, Q)
    class (RNRQCD)     , intent(inout) :: self
    Integer            , intent(in)    :: order
    real (dp)          , intent(in)    :: h, Q, nu
    real (dp)                          :: mu, m, lnu, R
    complex (dp)                       :: v
    integer                            :: i, j, n
    real (dp), dimension( 0:min(order,3) )   :: rQ

    n = min(order,3); mu = h * self%m1S

    Rexp = 0; if (order < 1) return

    if ( self%scheme(:4) == 'pole' ) then
      m = self%mass
    else
      R = self%m1S * vStar(q, self%m1S, self%gt)
      m = self%run%MSRMass(self%scheme(:4), self%ordMass, R, self%lambda, self%method)
    end if

    v = vC(q, m, m, self%gt); rQ = 0; lnu = log(nu)

    do i = 0, n
      do j = 1, n
        rQ(i) = rQ(i) + EXPterms3(v, i, j - i, m, mu, R, self%scheme)
      end do
    end do

    do i = 0, order
      do j = n + 1, 3 * (n - 1)
        rQ(i) = rQ(i) + higherOrderLogs3(v, i, j - i, m, mu, lnu, R, self%scheme)
      end do
    end do

    Rexp = 4 * m**2/q**2 * dot_product( PowList0( self%run%alphaQCD(mu)/Pi, n ), rQ )

  end function Rexp

!ccccccccccccccc

  real (dp) function RQCD(self, order, h, Q)
    class (RNRQCD)     , intent(inout)       :: self
    integer            , intent(in)          :: order
    real (dp)          , intent(in)          :: h, Q
    complex (dp)                             :: z
    real (dp)                                :: h2, mu, m, R, lr
    integer                                  :: n
    real (dp), dimension(5)                  :: b
    real (dp), dimension(0:4)                :: delta
    real (dp), dimension( 0:min(order,3) )   :: rQ
    real (dp), dimension( 0:min(order,3) - 1 )                   :: rQDer1
    real (dp), dimension( 0:min(order,3) - 2 )                   :: rQDer2
    real (dp), dimension( 0:min(order,3) - 1, min(order,3) )     :: Rcoef
    real (dp), dimension( 0:min(order,3) - 2, min(order,3) - 1 ) :: RcoefDer1

    if ( self%scheme(:4) == 'pole' ) then
      m = self%mass; delta = 0; h2 = h; mu = h * m
    else if ( self%scheme(:3) == 'MSR' ) then

      if (Q < self%QSwitch) then
        R = self%m1S * vStar(q, self%m1S, self%gt)
        m = self%run%MSRMass(self%scheme(:4), self%ordMass, R, self%lambda, self%method)
      else
        R = self%mass; m = self%MSRm
      end if

      delta = R * self%andim%betaQCD( self%scheme(:4) // 'delta')/m
      if (order > 1) lr = log(R/m); mu = h * self%m1S
      h2 = mu/R

    end if

    z = ( Q + (0,1) * self%gt )**2/4/m**2 ; n = min(order,3)

    Rcoef = 0; RQCD = 0; if (n < 0) return; rQ = PiCoef(z,n)

    b = getInverse( alphaMatch( 'pole', self%nl + 1) )
    call alphaReExpand( rQ(1:), b(:n) ); if (order > 1) Rcoef(0,:) = rQ(1:)

    if (order > 1) call self%andim%expandAlpha(Rcoef)

    if ( self%scheme(:3) == 'MSR' ) then

      if (order > 1) rQ(2:) = matmul( PowList0( lr, n - 1 ), Rcoef(:,2:) )
      if (order > 2) rQ(3) = rQ(3) - Rcoef(1,2) * delta(1)

      delta(1:3) = [ 0._dp, 3 * delta(1)**2, 6 * delta(1) * delta(2) - &
      4 * delta(1)**3 ] - 2 * delta(1:3)

      if (order > 0) then

        rQDer1 = PiCoefDer1(z,n)

        if (order > 2) then
          RcoefDer1 = 0; RcoefDer1(0,:) = rQDer1(1:)
          call self%andim%expandAlpha(RcoefDer1)
          rQDer1(2:) = matmul( PowList0( lr, n - 2 ), RcoefDer1(:,2:) )
        end if

        rQ(1:) = rQ(1:) + delta(1) * rQDer1

        if (order > 1) then
          rQDer2 = PiCoefDer2(z,n)
          rQ(2:) = rQ(2:) + delta(2) * rQDer1(:n-2) + delta(1)**2 * rQDer2/2
        end if

        if (order > 2) rQ(3) = rQ(3) + PiCoefDer3(z) * delta(1)**3/6 + &
         delta(3) * rQDer1(0) + delta(1) * delta(2) * rQDer2(0)

      end if

    end if

    if (order > 1) then
      Rcoef(0,:) = rQ(1:);  call self%andim%expandAlpha(Rcoef)
      rQ(2:) = matmul( PowList0( log(h2), n - 1 ), Rcoef(:,2:) )
    end if

    RQCD = 64 * Pi * m**2/3/Q**2 * &
    dot_product( PowList0( self%run%alphaQCD(mu)/Pi, n ), rQ )

  end function RQCD

! ccccccccccc

  function Delta1S(self) result(res)
    class (RNRQCD)  , intent(inout) :: self
    real (dp), dimension(0:4)       :: res

    res = self%m1Spiece

  end function Delta1S

! ccccccccccc

  real (dp) function Xsec(self, order, Q, h, nu)
    class (RNRQCD)  , intent(inout) :: self
    real (dp)          , intent(in) :: h, nu, Q
    integer            , intent(in) :: order
    complex (dp)                    :: vt
    real (dp)                       :: ac, ah, asLL, asNLL, auLL, En, c1run, L,&
    mu1, mu2, asNNLL, c2run, Vcc, V22, Vss, Vrr, VkkCACF, VkkCF2, Vkk, Vkkk1I, &
    Vkkk2T, VcsNNLL, DelmNNLL, rCoul, r2, rd, rr, rk, rkin, r1S, pre, inM, mP, &
    DeltaLL, aS, mPLL, DeltaNLL, DeltaNNLL, asPi

    Xsec = 0; inM = self%m1S

    if ( self%scheme(:4) == 'pole' ) then
      mP = inM; mPLL = inM
    end if

    mu1 = h * inM;  mu2 = nu * mu1;  ah = self%run%AlphaQCD(mu1)
    asLL = self%AlphaOb%alphaGenericFlavor('inverse', 1, self%nl, mu1, ah, mu2)
    ac = - VcsLL(asLL)/FPi

    if ( order == 2 .or. order == 3 ) then

      asNLL = self%AlphaOb%alphaGenericFlavor('inverse', 2, self%nl, &
      mu1, ah, mu2)

      auLL = self%AlphaOb%alphaGenericFlavor('inverse', 1, self%nl, &
      mu1, ah, mu2 * nu)

    end if

    if ( self%scheme(:2) == 'S1' .or. self%scheme(:3) == 'MSR' ) then

      if ( order == 1 ) then
        DeltaLL = ac**2/8;  mP = inM * (1 + DeltaLL)
      else if ( order == 2 ) then
        as = - VcsLL(asNLL)/FPi; L = Log(h * nu/as)
        DeltaLL = as**2/8; as = 3 * as/4
        DeltaNLL = DeltaLL * as * ( self%beta(0) * (L + 1) + self%a1/2 )/Pi
        mpLL = 1 + DeltaLL; mP = inM * (mpLL + DeltaNLL); mpLL = inM * mpLL
      end if

    end if

    if ( order == 1 ) then

      vt = vC(q, mp, mp, self%gt); ac = 4 * asLL/3

      Xsec =  18 * (self%Qt * mP/q)**2 *  ImagPart(  cI * vt - ac * (  &
      Euler - 0.5_dp + l2 + Log( - cI * mP * vt/h/nu/inM ) + &
      DiGam( 1 - cI * ac/(2 * vt) )  )  )

    else if ( order == 2 ) then

      c1run = self%MNLLc1Square(ah, asLL, auLL);  En = q - 2 * mP

      Xsec = (2 * mpLL * self%Qt/Q)**2 * c1run * &
      self%A1pole(2, En, mPLL, asNLL, 0._dp, mu2)

    else if ( order == 3 ) then

      asNNLL = self%AlphaOb%alphaGenericFlavor('inverse', 3, self%nl, &
      mu1, ah, mu2)

      Vcc = - ac

      c1run = self%MNNLLc1Square(ah, asLL, auLL, nu, h, 1._dp)
      c2run = self%MLLc2(ah, auLL)
      V22 = self%V2sLL(ah, asLL, auLL)/FPi
      Vss = self%VssLL(ah, asLL)/FPi
      Vrr = self%VrsLL(asLL, auLL)/FPi
      Vkk = - 28 * asLL**2/9;  VkkCACF = - 4 * asLL**2;  VkkCF2 = 8 * asLL**2/9
      Vkkk1I = asLL**2 * self%Vk1sLL( asLL, auLL)
      Vkkk2T = asLL**2 * self%Vk2sLL(asLL, auLL)
      VcsNNLL = self%VceffsNNLL(asNNLL, asLL, auLL)

      DelmNNLL = (V22/2 + Vss + 3 * Vrr/8) * Vcc**3/2 - (Vkk + 6 * Vkkk1I + &
      4 * Vkkk2T) * Vcc**2/8 + 5 * Vcc**4/128

      if ( self%scheme(:2) == 'S1' .or. self%scheme(:3) == 'MSR' ) then

        as = - VcsLL(asNNLL)/FPi; DeltaLL = as**2/8; mpLL = inM * (1 + DeltaLL)

        as = - VcsNNLL/FPi; L = Log(h * nu/as); asPi = asNNLL/Pi

        DeltaLL = as**2/8
        DeltaNLL  = DeltaLL * asPi * ( self%beta(0) * (L + 1) + self%a1/2 )
        DeltaNNLL = DeltaLL * asPi**2 * ( self%beta(0)**2 * (3 * L**2/4 + L  + &
        Zeta3/2 + Pi2/24 + 0.25_dp) + self%beta(0) * self%a1/2 * (3 * L/2 + 1) &
        + self%beta(1) * (L + 1)/4 + self%a1**2/16 + self%a2/8  )

        mp = inM * (1 + DeltaLL + DeltaNLL + DeltaLL**2 + DeltaNNLL)

      end if

      vt = vC(q, mpLL, inM, self%gt); En = q - 2 * mp

      rCoul = c1run * mpLL**2 * self%A1pole(3, En, mpLL, &
      asNNLL, VcsNNLL, h * inM * nu)/18/Pi

      r2 = 2 * c2run * inM**2/FPi * ImagPart(vt**2 * ggg(ac, vt, l2 - 0.5_dp) )

      rd = FPi * ( V22 + 2 * Vss) * ImagPart( dgd(ac, vt) )
      rr = FPi * Vrr * ImagPart( dgr(ac, vt) )

      rk = (  VkkCACF * ImagPart( dgkCACF(ac, vt) ) + VkkCF2 * &
      ImagPart( dgkCF2(ac, vt) ) + Vkkk1I * ImagPart( dgkk1I(ac, vt) ) + &
      Vkkk2T * ImagPart( dgkk2T(ac, vt) )  )

      rkin = ImagPart( dgkin(ac, vt) )
      r1S = 0

      if ( self%scheme(:2) == 'S1' .or. self%scheme(:3) == 'MSR' ) then
        r1S = DelmNNLL * inM**2/FPi * ImagPart( dggg(ac, vt) )
      end if

      Pre = 72 * Pi/Q**2

      Xsec = self%Qt**2 * pre * (rCoul + rr + r1S + r2 + rkin + rd + rk)

    end if

  contains

! ccccccccccc

    complex (dp) function ggg(a, v, X)
      real (dp)   , intent(in) :: a, X
      complex (dp), intent(in) :: v

      ggg = cI * vt - a * (  Log( - cI * v/nu/h ) + X + Euler + &
      DiGam( 1 - cI * a/2/v )  )

    end function ggg

! ccccccccccc

    complex (dp) function dggga(a, v, X)
      real (dp)   , intent(in) :: a, X
      complex (dp), intent(in) :: v

      dggga = - Euler - X - Log( - cI * v/nu/h ) - DiGam( 1 - cI * a/2/v ) &
      + cI * a * TriGam( 1 - cI * a/2/v )/2/v

    end function dggga

! ccccccccccc

    complex (dp) function dgggv(a, v, X)
      real (dp)   , intent(in) :: a, X
      complex (dp), intent(in) :: v

      dgggv = cI - a * ( 1/v + cI * a * TriGam( 1 - cI * a/2/v )/2/v**2 )

    end function dgggv

! ccccccccccc

  complex (dp) function dgk(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgk = - inM**2 * ( ggg( a, v, 2 * l2 - 2 )**2 + v**2 )/2/FPi/a

    end function dgk

! ccccccccccc

  complex (dp) function dgkCACF(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgkCACF = - inM**2 * ( ggg( a, v, l2 - 1.25_dp )**2 + v**2 )/2/FPi/a

    end function dgkCACF

! ccccccccccc

  complex (dp) function dgkk2T(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgkk2T = - inM**2 * ( ggg( a, v, l2 - 21._dp/16 )**2 + v**2 )/2/Pi/a

    end function dgkk2T

! ccccccccccc

  complex (dp) function dgkCF2(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgkCF2 = - inM**2 * ( ggg( a, v, l2 - 1 )**2 + v**2 )/2/FPi/a

    end function dgkCF2

! ccccccccccc

  complex (dp) function dgd(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgd = - inM**2 * ggg( a, v, l2 - 0.5_dp )**2/FPi**2

    end function dgd

! ccccccccccc

  complex (dp) function dgr(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgr = - inM**2/FPi**2 * ( ggg(a, v, l2 - 1.5_dp)**2 + v**2 + &
    v**2 * dggga(a, v, l2 - 0.5_dp) )

    end function dgr

! ccccccccccc

  complex (dp) function dgkin(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgkin = inM**2/4/FPi * (  a * ggg(a, v, l2 - 1.5_dp)**2 + a * v**2 + &
    2 * v**2 * ( ggg(a, v, l2 - 0.5_dp) + a * dggga(a, v, l2 - 0.5_dp) + &
    v/4 *  dgggv(a, v, l2 - 0.5_dp) )  )

    end function dgkin

! ccccccccccc

  complex (dp) function dggg(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dggg = - dgggv(a, v, 1._dp)/v

    end function dggg

! ccccccccccc

  complex (dp) function dgkk1I(a, v)
    real (dp)   , intent(in) :: a
    complex (dp), intent(in) :: v

    dgkk1I = - 3 * inM**2 * ( ggg( a, v, l2 - 17._dp/12 )**2 + v**2 )/FPi/a

    end function dgkk1I

  end function Xsec

! ccccccccccc

  real (dp) function A1pole(self, order, En, mass, asoft, VcsNNLL, musoft)
    class (RNRQCD), intent(in) :: self
    integer       , intent(in) :: order
    real (dp)     , intent(in) :: En, asoft, musoft, VcsNNLL, mass
    type (Toppik)              :: TP
    real (dp), dimension(2)    :: res
    real (dp)                  :: x0, x1, x2, aSuPi

    x0 = 0; x1 = 0; x2 = 0

    if ( order == 1 ) then
      x0 = 1
    else if ( order == 2 ) then
      aSuPi = asoft/FPi; x0 = self%xc01(aSuPi); x1 = self%xc11(aSuPi)
    else if ( order == 3 ) then
      aSuPi = asoft/FPi; x1 = self%xc12(aSuPi); x2 = self%xc22(aSuPi)
      x0 = self%a1 * asuPi + self%a2 * asuPi**2 - 3 * VcsNNLL/16/asoft/Pi
    end if

    TP = Toppik(self%nl, En, mass, self%gt, asoft, musoft, 175000000._dp, &
    175000000._dp, x0, x1, x2, 0._dp, 0._dp, 0._dp, &
    0._dp, 0._dp, 0._dp, 0._dp, 0, 0, 0._dp, 0)

    res = TP%CrossSec();  A1pole = 9 * res(1)/4

  end function A1pole

! ccccccccccc

  real (dp) function VcsLL(as)
    real (dp), intent(in) :: as

    VcsLL = - 16 * Pi * as/3

  end function VcsLL

! ccccccccccc

  real (dp) function VceffsNNLL(self, asNNLL, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: asNNLL, as, au

    VceffsNNLL = VcsLL(asNNLL) + 24 * as**3 * Pi * Log(au/as)/self%beta(0)

  end function VceffsNNLL

! ccccccccccc

  real (dp) function VrsLL(self, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: as, au

    VrsLL = VcsLL(as) * ( 1 - 8 * Log(au/as)/self%beta(0) )

  end function VrsLL

! ccccccccccc

  real (dp) function V2sLL(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: as, au, ah
    real (dp)                  :: rat1

    rat1 = as/ah

    V2sLL = - 4 * Pi * ah/3 * (  (1 - rat1) * (425._dp/117 - 319/self%beta(0)/39) &
    - (15 - self%beta(0))/(6 - self%beta(0)) * ( 1 - rat1**(1 - 6/self%beta(0)) ) - &
    616 * (33 - 3 * self%beta(0)) * ( 1 - rat1**(1 - 13/self%beta(0)/2) )/117/&
    (39 - 6 * self%beta(0)) ) + 4 * VcsLL(as) * log(au/as)/self%beta(0)/9

  end function V2sLL

! ccccccccccc

  real (dp) function VssLL(self, ah, as)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as

    VssLL = - 8 * ah * Pi/( 6 - self%beta(0) ) * (  1 - (as/ah)**( 1 - 6/self%beta(0) ) &
    * ( 21 - 2 * self%beta(0) )/9  )

  end function VssLL

! ccccccccccc

  real (dp) function Vk1sLL(self, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: as, au

    Vk1sLL = 16 * log(au/as)/self%beta(0)/9

  end function Vk1sLL

! ccccccccccc

  real (dp) function Vk2sLL(self, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: as, au

    Vk2sLL = 112 * log(au/as)/self%beta(0)/9

  end function Vk2sLL

! ccccccccccc

  real (dp) function VkeffsLL(self, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: as, au

    VkeffsLL = as**2/9 * ( 544 * log(au/as)/self%beta(0) - 28 )

  end function VkeffsLL

! ccccccccccc

  real (dp) function XiNLL(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au
    real (dp)                  :: rat, rat2

    rat = as/ah; rat2 = as/au

    XiNLL = ah * Pi * ( - 2 * ( 1276._dp/3 + 2278 * self%beta(0)/9 ) * &
    (1 - rat)/117/self%beta(0)**2 - 8 * ( 1 - rat**(1 - 6/self%beta(0)) ) * &
    (39 - 5 * self%beta(0))/27/(6 - self%beta(0))**2 + 9856 * &
    ( 1 - rat**(1 - 13/self%beta(0)/2) ) * (33 - 3 * self%beta(0))/ &
    351/(39 - 6 * self%beta(0))**2 + 16 * (152 * self%beta(0) + 2 * self%beta(0)**2 &
    - 957) * Log(rat)/9/(6 - self%beta(0))/self%beta(0)**2/(39 - 6 * self%beta(0)) &
    - 7072 * ( rat - 1 + rat2 * Log(rat2) )/81/self%beta(0)**2 )

  end function

! ccccccccccc

  real (dp) function XiNNLLnonmix(self, ah, as, au, hh, ss)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au, hh, ss
    real (dp)                  :: rat, rat2

    rat = as/ah; rat2 = au/as

    XiNNLLnonmix = ah**2 * (  8 * ( 36 * (3 + 2 * l2) + 32 * (19 - 18 * l2)/9 &
    + 9 * (1 + 18 * l2) )/27/self%beta(0) + 3536 * Log(hh)/81/self%beta(0) ) * &
    ( 1 - rat + 2 * Log(rat2) ) + ss * ( 8 * ah**2 * (self%beta(0) - 12) * (957 &
    - 152 * self%beta(0) - 2 * self%beta(0)**2) * (1 - rat)/27/(6 - self%beta(0))&
    /self%beta(0)**2/(39 - 6 * self%beta(0)) + ah**2 * (1 - rat**2) * (61248 + &
    450491 * self%beta(0)/3 - 50537 * self%beta(0)**2/3)/8424/self%beta(0)**2 &
    - 1232 * ah**2 * (3465 - 513 * self%beta(0) + 18 * self%beta(0)**2) * &
    ( 1 - rat**(2 - 13/self%beta(0)/2) )/3159/(39 - 6 * self%beta(0))/(39 - &
    12 * self%beta(0)) + ah**2 * (1116 - 198 * self%beta(0) + 5 * self%beta(0)**2) &
    * ( 1 - rat**(2 - 6/self%beta(0)) )/81/(6 - self%beta(0))/(3 - self%beta(0)) + &
    2 * ah**2 * (2521 * self%beta(0)/9 - 7072._dp/3) * ( 2.5_dp - 2 * rat - &
    rat**2/2 + (4 - rat**2) * Log(rat2) )/27/self%beta(0)**2)

  end function XiNNLLnonmix

! ccccccccccc

  real (dp) function XiNNLLmixUsoft(self, ah, as)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as
    real (dp)                  :: rat, DiLog

    rat = as/ah

    XiNNLLmixUsoft = - 1768 * ah**2 * ( 3 * (47 + 6 * Pi2) - 5 * self%nl ) * &
    ( 3 - 2 * rat - rat**2 - 4 * Log(2 - rat) )/729/self%beta(0)**2 - 1768 * &
    ah**2 * self%beta(1) * ( Pi2/6 - 1.75_dp - 2 * DiLog(rat/2) - Log(rat/2)**2 &
    + rat**2 * ( 0.75_dp - Log(rat)/2 ) + rat * (  1 - Log( rat/(2 - rat) )  ) &
    + Log( rat/(2 - rat) )**2)/81/self%beta(0)**3

  end function XiNNLLmixUsoft

! ccccccccccc

  real (dp) function MLLc2(self, ah, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, au

    MLLc2 = - 1._dp/6 - 32 * log(au/ah)/9/self%beta(0)

  end function MLLc2

! ccccccccccc

  real (dp) function MNLLc1(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au

    MNLLc1 = (1 - 8 * ah/3/Pi) * exp( self%xiNLL(ah, as, au) )

  end function MNLLc1

! ccccccccccc

  real (dp) function MNLLc1Square(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au

    MNLLc1Square = (1 - 16 * ah/3/Pi) * exp( 2 * self%xiNLL(ah, as, au) )

  end function MNLLc1Square

! ccccccccccc

  real (dp) function MNLLplusNNLLnonmixc1(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au

    MNLLplusNNLLnonmixc1 = 1 - 8 * ah/3/Pi + ah**2 * &
    (0.04127900074317465_dp * self%nl - 3.55919055282743_dp - 14._dp/9 * l2)

    MNLLplusNNLLnonmixc1 = MNLLplusNNLLnonmixc1 * &
    exp( self%xiNLL(ah, as, au) + self%xiNNLLnonmix(ah, as, au, 1._dp, 1._dp) )

  end function MNLLplusNNLLnonmixc1

! ccccccccccc

  real (dp) function XiNNLLSoftMixLogc1(ah, nu, hh)
    real (dp), intent(in) :: ah, nu, hh

    XiNNLLSoftMixLogc1 = ah**3/Pi * ( 1361._dp/405 - 140 * Log(hh)/81 ) * Log(nu)

  end function XiNNLLSoftMixLogc1

! ccccccccccc

  real (dp) function MNNLLAllc1InclSoftMixLog(self, ah, as, au, nu, hh, ss)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au, nu, hh, ss

    MNNLLAllc1InclSoftMixLog = 1 - 0.8488263631567752_dp * ah + &
    ah**2 * ( 0.04467257170808474_dp * self%nl - 4.654387355189673_dp - &
    (2.5925925925925926_dp + 0.13509491152311703_dp * self%beta(0)) * Log(hh) )

    MNNLLAllc1InclSoftMixLog = MNNLLAllc1InclSoftMixLog * &
    exp( self%xiNLL(ah, as, au) + self%xiNNLLmixUsoft(ah, as) + ss * &
    xiNNLLSoftMixLogc1(ah, nu, hh) + self%xiNNLLnonmix(ah, as, au, hh, ss) )

  end function MNNLLAllc1InclSoftMixLog

! ccccccccccc

  real (dp) function MNNLLc1Square(self, ah, as, au, nu, hh, ss)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au, nu, hh, ss

    MNNLLc1Square = 1 - 1.6976527263135504_dp * ah + ah**2 *  ( &
    0.0825580014863493_dp * self%nl - 8.554332805940287_dp - &
    ( 5.185185185185185_dp + 0.27018982304623407_dp * self%beta(0) ) * Log(hh) )

    MNNLLc1Square = MNNLLc1Square * &
    exp(  2 * ( self%xiNLL(ah, as, au) + self%xiNNLLmixUsoft(ah, as) + ss * &
    xiNNLLSoftMixLogc1(ah, nu, hh) + self%xiNNLLnonmix(ah, as, au, hh, ss) )  )

  end function MNNLLc1Square

! ccccccccccc

  real (dp) function qFromV(v, m, gt)
    real (dp), intent(in) :: v, m, gt

    qFromV = (8 * m**2 * v**2 + 4 * m**2 * v**4 - gt**2)/4/m/v**2

  end function qFromV

! ccccccccccc

  complex (dp) function vC(q, m1, m2, gt)
    real (dp), intent(in) :: q, m1, m2, gt

    vC = sqrt(  ( q - 2 * m1 + cI * gt ) /m2  )

  end function vC

! ccccccccccc

  real (dp) function vStar(q, m, gt)
    real (dp), intent(in) :: q, m, gt

    vStar = 0.05_dp + Abs( vC(q, m, m, gt) )

  end function vStar

! ccccccccccc

  real (dp) function vRootStar(q, m, gt)
    real (dp), intent(in) :: q, m, gt

    vRootStar = sqrt( vStar(q, m, gt) )

  end function vRootStar

! ccccccccccc

  real (dp) function SwitchOff(q, m, gt, v0, v1)
    real (dp), intent(in) :: q, m, gt, v0, v1
    real (dp)             :: v

    SwitchOff = 0; v = realpart( vC(q, m, m, gt) )

    if (v < v0) then
      SwitchOff = 1
    else if ( v > v0 .and. v <= (v0 + v1)/2 ) then
      SwitchOff = 1 - 2 * (v - v0)**2/(v1 - v0)**2
    else if ( v <= v1 .and. v > (v0 + v1)/2 ) then
      SwitchOff = 2 * (v - v1)**2/(v1 - v0)**2
    end if

  end function SwitchOff

! ccccccccccc

  real (dp) function xc01(self, asuPi)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: asuPi
    xc01 = 1 + self%a1 * asuPi
  end function xc01

! ccccccccccc

  real (dp) function xc11(self, asuPi)
    class (RNRQCD), intent(in) :: self
    real (dp), intent(in)      :: asuPi
    xc11 = self%beta(0) * asuPi
  end function xc11

! ccccccccccc

  real (dp) function xc12(self, asuPi)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: asuPi

    xc12 = self%beta(0) * aSuPi + aSuPi**2 * (self%beta(1) + 2 * self%beta(0) * self%a1)

  end function xc12

! ccccccccccc

  real (dp) function xc22(self, asuPi)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: asuPi

    xc22 = asuPi**2 * self%beta(0)**2

  end function xc22

! ccccccccccc

  real (dp) function higherOrderLogs3(v, i, j, m, mu, lnu, R, scheme)
    real (dp)          , intent(in) :: m, mu, lnu, R
    character (len = *), intent(in) :: scheme
    complex (dp)       , intent(in) :: v
    integer            , intent(in) :: i, j
    complex (dp)                    :: higherOrderLogs

    higherOrderLogs = 0

    if (i == 1 .and. j == 3) then
      higherOrderLogs = 128 * cI * v**3 * lnu/9
    else if (i == 2 .and. j == 1) then
      higherOrderLogs = - 280 * cI * v * Pi2 * Log(v)/27
    else if (i == 2 .and. j == 2) then
      higherOrderLogs = 170 * cI * v**2 * Pi2 * lnu/27
    else if (i == 2 .and. j == 3) then
      higherOrderLogs = - 1472 * cI * v**3 * lnu**2/27
    else if (i == 3 .and. j == 1) then
      higherOrderLogs = - cI * lnu * v * Pi2 * ( 64343 + 140280 * l2 &
      - 8940 * lnu + 2955 * Pi2 + 96600 * Log(mu/m) + 17880 * Log(v) )/1215
    else if (i == 3 .and. j == 2) then
      higherOrderLogs = - 2461 * cI * v**2 * Pi2 * lnu**2/27
    else if (i == 3 .and. j == 3) then
      higherOrderLogs = 67712 * cI * v**3 * lnu**3/243
    end if

    higherOrderLogs3 = ImagPart(higherOrderLogs)

    if ( scheme(:4) == 'pole' ) return

    higherOrderLogs = 0

    if (i == 2 .and. j == 2) then
      higherOrderLogs = - 512 * cI * R * v * lnu/9/m
    else if (i == 3 .and. j == 3) then
      higherOrderLogs = - 8 * cI * (  lnu * ( 2 * R * (85 * m * Pi2 * v - 192 * R) + &
      805 * m**2 * Pi2 * v**2 * Log(R/mu) ) + 23 * m * Pi2 * v * ( 6 * R - &
      35 * m * v * Log(R/mu) ) * Log(v)  )/81/m**2/v
    else if (i == 3 .and. j == 2) then
      higherOrderLogs = 5888 * cI * lnu**2 * R * v/27/m
    end if

    higherOrderLogs3 = higherOrderLogs3 + ImagPart(higherOrderLogs)

  end function higherOrderLogs3

! ccccccccccc

  real (dp) function EXPterms3(v, i, j, m, mu, R, scheme)
    real (dp)          , intent(in) :: m, mu, R
    character (len = *), intent(in) :: scheme
    complex (dp)       , intent(in) :: v
    integer            , intent(in) :: i, j
    complex (dp)                    :: l1, l3, EXPterms

    EXPterms = 0

    if (i == 0 .and. j == 1) then
      EXPterms = 2 * cI * v
    else if (i == 0 .and. j == 3) then
      EXPterms = 7 * v**3 * cI/12
    else if (i == 1 .and. j == 0) then
      EXPterms = - 8 * Pi * Log(- cI * v)/3
    else if (i == 1 .and. j == 1) then
      EXPterms = - 32 * cI * v/3
    else if (i == 1 .and. j == 2) then
      EXPterms = v**2 * Pi * ( 33 - 40 * Log(-2 * cI * m * v/mu) )/9
    else if (i == 2 .and. j == -1) then
      EXPterms = 8 * cI * Pi2**2/v/27
    else if (i == 2 .and. j == 0) then
      l1 = Log(- cI * m * v/mu)
      EXPterms = 2 * Pi * ( l1 * (69 * l1 + 138 * l2 - 43) + &
      192 * log(- cI * v) )/27
    else if (i == 2 .and. j == 1) then
      EXPterms = - 14.513280990412623_dp * cI * v * Pi2 * ( Log(-cI * v) - &
      0.31755295480049456_dp - 0.28545651550321627_dp * Log(- cI * m * v/mu) )
    else if (i == 3 .and. j == -2) then
      EXPterms = - 32 * Pi3 * Zeta3/27/v**2
    else if (i == 3 .and. j == -1) then
      EXPterms = - 4 * cI * Pi2 * (53 * Pi2 + 828 * Zeta3 + &
      138 * Pi2 * Log(-2 * cI * m * v/mu) )/243/v
    else if (i == 3 .and. j == 0) then
      l1 = Log(- cI * m * v/mu); l3 = Log(- cI * v)
      EXPterms = - 1.3234297814023874_dp * Pi3 * ( l1**3 + l1**2 * &
      (3 * l2 - 4.565197036493422_dp) + (-4.748801815703692_dp - &
      20.890672104680654_dp * l2 - 14.621887456729755_dp * l3) * l3 + l1 * &
      (20.107979523788522_dp + l2 * (5.491493383742911_dp + 3 * l2) + &
      14.621887456729755_dp * l3) )
    end if

    EXPterms3 = ImagPart(EXPterms)

    if ( scheme(:4) == 'pole' ) return

    EXPterms = 0

    if (i == 1 .and. j == 0) then
      EXPterms = 8 * ( Pi * Log(R/m) - cI * R/m/v )/3
    else if (i == 1 .and. j == 2) then
      EXPterms = 5 * v * (3 * cI * R + 8 * m * Pi * v * Log(R/mu) )/9/m
    else if (i == 2 .and. j == -1) then
      EXPterms =  16 * R * (2 * m * Pi * v - cI * R)/9/m**2/v**3
    else if (i == 2 .and. j == 0) then
      EXPterms = (   552 * m * Pi3 * v + 4416 * cI * R * Log(R/mu) + &
      32 * m * Pi * v * (  192 * Log(m/mu) - Log(R/mu) * (149 + 138 * l2 + &
      69 * Log(R/mu) )  ) - 3 * cI * R * ( 291 + 16 * Pi2 * (2 * l2 - 1) -  &
      48 * Zeta3 )   )/432/m/v
    else if (i == 2 .and. j == 1) then
      EXPterms = (  2 * R * (  4 * m * Pi * v * (40 * l2 - 19 - 8 * cI * Pi) - &
      51 * cI * R + 64 * m * Pi * v * Log(v * m/R) )  )/27/m**2/v
    else if (i == 3 .and. j == -2) then
      EXPterms = 32 * R * ( 12 * m * Pi * R * v - 6 * cI * R**2 &
      + cI * m**2 * Pi2**2 * v**2)/81/m**3/v**5
    else if (i == 3 .and. j == -1) then
      EXPterms = R * ( -12.20650397291053_dp * cI * R - (7.776708753396424_dp, &
      - 67.25952628890525_dp) * m * v + 13.62962962962963_dp * cI * R * Log(R/mu)&
      - 42.818744315594216_dp * m * v * Log((m*R*v)/mu**2) )/m**2/v**3
    else if (i == 3 .and. j == 0) then
      EXPterms = ( -2.8641975308641974_dp * cI * R**3 + 4.964492094561649_dp * &
      m * R**2 * v + (214.36438198725799_dp,115.57430402546646_dp) * m**2 * R *&
      v**2 + 251.67998748274093_dp * m**3 * v**3 + m**2 * v**2 * &
      (-600.0037412368929_dp * m * v * Log(m/mu)**2 + ( 83.07628673312446_dp * cI * R &
      - 52.34229117492672_dp * m * v) * Log(R/mu) + (-39.18518518518519_dp * cI * R - &
      102.00205692336024_dp * m * v) * Log(R/mu)**2 + 41.03462996911112_dp * m * v * &
      Log(R/mu)**3 + Log(m/mu) * ( -54.518518518518526_dp * cI * R + 789.059507252147_dp&
      * m * v + 600.0037412368929_dp * m * v * Log(R/mu) ) + &
     136.46860406444546_dp * cI * R * Log(v)) )/m**3/v**3
    end if

    EXPterms3 = EXPterms3 + ImagPart(EXPterms)

  end function EXPterms3

! ccccccccccc

end module RNRQCDClass
