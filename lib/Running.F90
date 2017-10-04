
module RunningClass
  use AnomDimClass;  use AlphaClass; use Constants, only: dp, Pi, ExpEuler, d1mach, prec
  use QuadPack, only: qags         ; implicit none   ;  private

  public                        :: DeltaComputer, AddAlpha

!ccccccccccccccc

  type, public                  :: Running
    private
    character (len = 5)         :: str
    integer                     :: nf, runMass
    type (AnomDim)              :: andim
    type (Alpha)                :: AlphaOb
    real (dp), dimension(0:4)   :: beta, gamma
    real (dp), dimension(4)     :: bHat, sCoefP, sCoefN, gammaRn, gammaRp
    real (dp)                   :: muLambda, mH, mL
    real (dp), dimension(0:4,4) :: tab

    contains

    procedure, pass(self), public :: MSbarMass, MSRMass, lambdaQCD, DiffDeltaMu,  &
    DiffDelta, orders, adim, DiffDeltaHadron, MSbarDeltaMu, MSbarMassLow, &
    AlphaAll, DeltaGapMatching, DiffRMass, mmFromMSR, PoleMass, gammaR, sCoef,  &
    SetMTop, SetMBottom, SetMCharm, SetLambda, SetAlpha, scales, numFlav,DiffR, &
    PSdelta, OptimalR, MSRMatching, MSREvol, AlphaQED, scheme, OptimalR2

    procedure, pass(self), private :: alphaQCDReal, alphaQCDComplex, &
    wTildeReal, wTildeComplex, kTildeReal, kTildeComplex, RunningMass

    generic, public :: alphaQCD => alphaQCDReal, alphaQCDComplex
    generic, public :: wTilde   => wTildeReal  , wTildeComplex
    generic, public :: kTilde   => kTildeReal  , kTildeComplex

  end type Running

!ccccccccccccccc

  interface Running
    module procedure  InitRun
  end interface Running

  contains

!ccccccccccccccc

   type (Running) function InitRun(nf, runMass, AlphaOb, muLambda)
    integer     , intent(in)  :: runMass, nf
    real (dp)   , intent(in)  :: muLambda
    type (Alpha), intent(in)  :: AlphaOb
    real (dp), dimension(0:4) :: sCoef

    InitRun%AlphaOb = AlphaOb          ;  InitRun%andim    = AlphaOb%adim(nf)
    InitRun%str     = AlphaOb%scheme() ;  InitRun%muLambda = muLambda
    InitRun%nf      = nf               ;  InitRun%runMass  = runMass
    InitRun%beta    = InitRun%andim%betaQCD('beta')
    InitRun%gamma   = InitRun%andim%betaQCD('mass')

    if (nf == 6) then
      InitRun%mH = 0                   ; InitRun%mL = AlphaOb%scales('mT')
    else if (nf == 5) then
      InitRun%mH = AlphaOb%scales('mT'); InitRun%mL = AlphaOb%scales('mB')
    else if (nf == 4) then
      InitRun%mH = AlphaOb%scales('mB'); InitRun%mL = AlphaOb%scales('mC')
    else if (nf == 3) then
      InitRun%mH = AlphaOb%scales('mC'); InitRun%mL = 0
    else
      InitRun%mH = 0; InitRun%mL = 0
    end if

    InitRun%tab = MSbarDeltaPiece(nf - 1, 1)
    sCoef = InitRun%andim%betaQCD('sCoefMSRp')
    InitRun%sCoefP = sCoef(1:); sCoef = 0;
    sCoef = InitRun%andim%betaQCD('sCoefMSRn')
    InitRun%sCoefN = sCoef(1:)
    sCoef = 0; sCoef = InitRun%andim%betaQCD('bHat')
    InitRun%bHat = sCoef(1:)
    sCoef = 0; sCoef = InitRun%andim%betaQCD('MSRpdelta')
    InitRun%gammaRp = sCoef(1:)
    sCoef = 0; sCoef = InitRun%andim%betaQCD('MSRndelta')
    InitRun%gammaRn = sCoef(1:)

   end function InitRun

!ccccccccccccccc

  function gammaR(self, type) result(res)
    class (Running)    , intent(in) :: self
    character (len = *), intent(in) :: type
    real (dp), dimension(4)         :: res

    if ( type(:4) == 'MSRn' ) then
      res = self%gammaRn
    else if ( type(:4) == 'MSRp' ) then
      res = self%gammaRp
    end if

  end function gammaR

!ccccccccccccccc

  function sCoef(self, type) result(res)
    class (Running)    , intent(in) :: self
    character (len = *), intent(in) :: type
    real (dp), dimension(4)         :: res

    if ( type(:4) == 'MSRn' ) then
      res = self%sCoefN
    else if ( type(:4) == 'MSRp' ) then
      res = self%sCoefP
    end if

  end function sCoef

!ccccccccccccccc

  subroutine SetLambda(self, lambda)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: lambda

    self%muLambda = lambda

  end subroutine SetLambda

!ccccccccccccccc

  subroutine SetAlpha(self, alpha)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: alpha

    call self%AlphaOb%SetAlpha(alpha, 0._dp)

  end subroutine SetAlpha

!ccccccccccccccc

  subroutine SetMTop(self, mT, muT)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: mT, muT

    call self%AlphaOb%SetMTop(mT, muT)
    if (self%nf == 6) self%mL = mT; if (self%nf == 5) self%mH = mT

  end subroutine SetMTop

!ccccccccccccccc

  subroutine SetMBottom(self, mB, muB)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: mB, muB

    call self%AlphaOb%SetMBottom(mB, muB)
    if (self%nf == 5) self%mL = mB; if (self%nf == 4) self%mH = mB

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, mC, muC)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: mC, muC

    call self%AlphaOb%SetMCharm(mC, muC)
    if (self%nf == 4) self%mL = mC; if (self%nf == 3) self%mH = mC

  end subroutine SetMCharm

!ccccccccccccccc

   character (len = 5) function scheme(self)
    class (Running), intent(in) :: self
    scheme = self%str
   end function scheme

!ccccccccccccccc

  function PSdelta(self, R, mu, lg) result(res)
    class (Running), intent(in) :: self
    real (dp)      , intent(in) :: R, mu, lg
    real (dp), dimension(4)     :: res
    real (dp), dimension(0:4,4) :: coef

    coef(0,:) = R * self%andim%PScoef(lg)
    call self%andim%expandAlpha(coef)
    call AddAlpha(  coef, powList( self%alphaQCD(mu)/Pi, 4 )  )

    res = DeltaComputer(coef, powList( log(mu/R), 4 ), 0)

  end function PSdelta

!ccccccccccccccc

   real (dp) function wTildeReal(self, order, gamma, mu0, mu1)
    class (Running)          , intent(in) :: self
    integer                  , intent(in) :: order
    real (dp)                , intent(in) :: mu0, mu1
    real (dp), dimension(0:3), intent(in) :: gamma

    wTildeReal = self%andim%wTilde( order, gamma, self%alphaQCD(mu0), self%alphaQCD(mu1) )

   end function wTildeReal

!ccccccccccccccc

   real (dp) function wTildeComplex(self, order, gamma, mu0, mu1)
     class (Running)          , intent(in) :: self
     integer                  , intent(in) :: order
     real (dp)                , intent(in) :: mu1
     complex (dp)             , intent(in) :: mu0
     real (dp), dimension(0:3), intent(in) :: gamma

     wTildeComplex = self%andim%wTilde( order, gamma, self%alphaQCD(mu0), self%alphaQCD(mu1) )

   end function wTildeComplex

!ccccccccccccccc

   real (dp) function kTildeReal(self, order, gamma, mu0, mu1)
     class (Running)          , intent(in) :: self
     integer                  , intent(in) :: order
     real (dp)                , intent(in) :: mu0, mu1
     real (dp), dimension(0:3), intent(in) :: gamma

     kTildeReal = self%andim%kTilde( order, gamma, self%alphaQCD(mu0), self%alphaQCD(mu1) )

   end function kTildeReal

!ccccccccccccccc

   real (dp) function kTildeComplex(self, order, gamma, mu0, mu1)
     class (Running)          , intent(in) :: self
     integer                  , intent(in) :: order
     complex (dp)             , intent(in) :: mu0
     real (dp)                , intent(in) :: mu1
     real (dp), dimension(0:3), intent(in) :: gamma

     kTildeComplex = self%andim%kTilde( order, gamma, self%alphaQCD(mu0), self%alphaQCD(mu1) )

   end function kTildeComplex

!ccccccccccccccc

   pure real (dp) function alphaQCDReal(self, mu)
     real (dp)      , intent(in) :: mu
     class (Running), intent(in) :: self

     alphaQCDReal = self%AlphaOb%alphaQCD(self%nf, mu)

   end function alphaQCDReal

!ccccccccccccccc

   pure real (dp) function alphaQED(self, mu)
     real (dp)      , intent(in) :: mu
     class (Running), intent(in) :: self

     alphaQED = self%AlphaOb%alphaQED(self%nf, mu)

   end function alphaQED

!ccccccccccccccc

   complex (dp) function alphaQCDComplex(self, mu)
     complex (dp)   , intent(in) :: mu
     class (Running), intent(in) :: self

     alphaQCDComplex = self%AlphaOb%alphaQCD(self%nf, mu)

   end function alphaQCDComplex

!ccccccccccccccc

   real (dp) function MSbarMass(self, mu)
     real (dp)      , intent(in) :: mu
     class (Running), intent(in) :: self

     if ( self%str == 'MSbar') then
       MSbarMass = self%RunningMass(self%mL, mu)
     else
       MSbarMass = self%mL
     end if

   end function MSbarMass

 !ccccccccccccccc

   integer function numFlav(self)
     class (Running), intent(in) :: self
     numFlav = self%nf
   end function numFlav

!ccccccccccccccc

   real (dp) function PoleMass(self, order, mu)
     real (dp)      , intent(in)   :: mu
     integer        , intent(in)   :: order
     class (Running), intent(in)   :: self

     real (dp)    , dimension(0:3) :: delta

     delta(0) = 1; delta(1:) = self%MSbarDeltaMu(mu)

     PoleMass = self%MSbarMass(mu) * sum(  delta( :min(3,order) )  )

   end function PoleMass

!ccccccccccccccc

   real (dp) function MSbarMassLow(self, mu)
     real (dp)      , intent(in) :: mu
     class (Running), intent(in) :: self

     real (dp)                   :: aMatch

     aMatch = self%AlphaOb%alphaQCD(self%nf + 1, self%mH)/Pi

     MSbarMassLow = (1 + 0.10356715567659536_dp * aMatch**2) * self%RunningMass(self%mH, mu)

   end function MSbarMassLow

!ccccccccccccccc

   pure real (dp) function RunningMass(self, m, mu)
     real (dp)      , intent(in) :: m, mu
     class (Running), intent(in) :: self
     real (dp)                   :: a0, a1

     a0 = self%alphaQCD(m);  a1 = self%alphaQCD(mu)

     RunningMass = m * exp( self%andim%wTilde( self%runMass, self%gamma, a0, a1) )

   end function RunningMass

!ccccccccccccccc

   function MSbarDeltaMu(self, mu) result(delta)
     real (dp)      , intent(in) :: mu
     class (Running), intent(in) :: self
     real (dp), dimension(4)     :: alphaList, delta
     real (dp), dimension(0:4)   :: logList

     alphaList = PowList( self%alphaQCD(mu)/Pi, 4 )
     LogList = PowList0(  log( mu/self%MSbarMass(mu) ), 4  )

     delta = alphaList * matmul( logList, self%tab )

   end function MSbarDeltaMu

!ccccccccccccccc

  pure subroutine AddAlpha(Mat, alphaList)
    real (dp), dimension(:  ), intent(in   ) :: alphaList
    real (dp), dimension(:,:), intent(inout) :: Mat
    integer                                  :: i

    do i = 1, Min( size(alphaList), size(Mat,2) )
      Mat(:,i) = alphaList(i) * Mat(:,i)
    enddo

  end subroutine AddAlpha

!ccccccccccccccc

  real (dp) function MSRMass(self, type, order, R, lambda, method)
    class (Running)                    , intent(in) :: self
    real (dp)                          , intent(in) :: R
    real (dp), optional                , intent(in) :: lambda
    character (len = *)                , intent(in) :: type
    character (len = *), optional      , intent(in) :: method
    integer                            , intent(in) :: order
    real (dp)                        , dimension(4) :: a
    character (len = 8)                             :: met
    integer                                         :: i
    real (dp)                                       :: alphaM, lan, matching

    if ( .not. present(lambda) ) then
      lan = 1
    else
      lan = lambda
    end if

    if ( .not. present(method) ) then
      met = 'analytic'
    else
      met = method
    end if

    matching = 1

     if ( self%runMass > 0 .and. type(:4) == 'MSRn' ) then
      alphaM = self%AlphaOb%alphaQCD(self%nf + 1, self%mH)/Pi
      a = self%MSRMatching(type); i = min(order, 4)
      matching = 1 + dot_product( a(:i), PowList(alphaM, i) )
     end if

     MSRMass = self%mH * matching + self%MSREvol(type, order, R, lan, met)

  end function MSRMass

!ccccccccccccccc

  real (dp) function MSREvol(self, type, order, R, lambda, method)
    class (Running)                    , intent(in) :: self
    real (dp)                          , intent(in) :: R, lambda
    character (len = *)                , intent(in) :: type, method
    integer                            , intent(in) :: order
    real (dp)                        , dimension(4) :: gammaR
    real (dp), dimension(self%runMass,self%runMass) :: gammaLog
    real (dp), dimension(self%runMass)              :: gamm, alphaList, lglist
    integer                                         :: neval, ier
    real (dp)                                       :: abserr, lan, delta, &
    mu, Rstep, delta1, delta2

    gammaR = self%gammaR(type)

    if ( abs(lambda - 1) >= tiny(1._dp) ) then

      gammaLog(1,:) = gammaR(:self%runMass) ; lan = lambda
      call self%andim%expandAlpha( gammaLog )

      gamm = lambda * (   gammaLog(1,:) + &
      matmul(  powList( - log(lambda), self%runMass-1 ), gammaLog(2:,:)  )   )

    else

      lan = 1

    end if

    if ( method(:8) == 'analytic' ) then

      if ( abs(lambda - 1) <= tiny(1._dp) ) then
        MSREvol = self%DiffR( self%sCoef(type), self%runMass, self%mH, R )
      else
        MSREvol = self%DiffR( self%andim%sCoefRecursive(gamm), &
        self%runMass, self%mH/lambda, R/lambda )
       end if

       MSREvol = self%lambdaQCD(self%runMass) * MSREvol

    else if ( method(:7) == 'numeric' ) then

      if ( abs(lambda - 1) <= tiny(1._dp) ) then
        gamm = gammaR(:self%runMass)
      end if

      gamm = self%andim%GammaRComputer( gamm )

      call qags( Integrand, R/lan, self%mH/lan, prec, prec, MSREvol, abserr, neval, ier )

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( ( self%mH - R )/delta )  )
      delta = (self%mH - R)/neval; MSREvol = 0; Rstep = self%mH + delta
      gammaLog(1,:) = gammaR(:self%runMass)
      call self%andim%expandAlpha( gammaLog )

      do ier = 1, neval

        Rstep = Rstep - delta; mu = (Rstep - delta/2)/lan
        alphaList = powList( self%alphaQCD(mu)/Pi, self%runMass )

        lgList = powList( log(mu/Rstep), self%runMass )
        delta1 = Rstep * dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )

        lgList = powList( log(mu/(Rstep - delta) ), self%runMass )
        delta2 = (Rstep - delta) * dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )

        MSREvol = MSREvol + delta1 - delta2

      end do

    end if

    contains

!ccccccccccccccc

    real (dp) function Integrand(R)
      real (dp)             , intent(in) :: R
      real (dp), dimension(self%runMass) :: aPi

      aPi = powList( self%alphaQCD(R)/4/Pi, self%runMass )

      Integrand = dot_product( gamm, aPi )

    end function Integrand

  end function MSREvol

!ccccccccccccccc

  real (dp) function mmFromMSR(self, type, mass, order, R)
    class (Running) , intent(inout) :: self
    character (len = *), intent(in) :: type
    real (dp)       , intent(in   ) :: R, mass
    integer         , intent(in   ) :: order
    integer                         :: IFLAG
    real (dp)                       :: a, b, c, rat, massOr

     a = mass/2; b = 2 * mass; massOr = self%mH

    if (self%nf == 5) rat = self%alphaOb%scales('muT')/self%mH
    if (self%nf == 4) rat = self%alphaOb%scales('muB')/self%mH
    if (self%nf == 3) rat = self%alphaOb%scales('muC')/self%mH

    call DFZERO(MSR, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG); mmFromMSR = a

    if (self%nf == 5) call self%SetMTop(   massOr, rat * massOr)
    if (self%nf == 4) call self%SetMBottom(massOr, rat * massOr)
    if (self%nf == 3) call self%SetMCharm( massOr, rat * massOr)

   contains

   real (dp) function MSR(mm)
     real (dp), intent(in)   :: mm

     if (self%nf == 5) call self%SetMTop(   mm, rat * mm)
     if (self%nf == 4) call self%SetMBottom(mm, rat * mm)
     if (self%nf == 3) call self%SetMCharm (mm, rat * mm)

     MSR = self%MSRMass(type, order, R) - mass

   end function MSR

  end function mmFromMSR

!ccccccccccccccc

  real (dp) function OptimalR(self, type, n, order, lambda, method)
    class (Running)                    , intent(in) :: self
    character (len = *)                , intent(in) :: type
    real (dp)          , optional      , intent(in) :: lambda
    character (len = *), optional      , intent(in) :: method
    integer                            , intent(in) :: order
    real (dp)                          , intent(in) :: n
    integer                                         :: IFLAG
    real (dp)                                       :: a, b, c

    a = 0.5_dp; b = self%mH

    call DFZERO(root, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG); OptimalR = a

  contains

    real (dp) function root(x)
      real (dp), intent(in) :: x

      root = n * x - 4 * self%MSRMass(type, order, x, lambda, method) &
      * self%alphaQCD(x)/3

    end function root

  end function OptimalR

!ccccccccccccccc

  real (dp) function OptimalR2(self, n, mass)
    class (Running)                    , intent(in) :: self
    real (dp)                          , intent(in) :: n, mass
    integer                                         :: IFLAG
    real (dp)                                       :: a, b, c

    a = 0.5_dp; b = self%mH

    call DFZERO(root, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG); OptimalR2 = a

  contains

    real (dp) function root(x)
      real (dp), intent(in) :: x

      root = n * x - 4 * Mass * self%alphaQCD(x)/3

    end function root

  end function OptimalR2

!ccccccccccccccc

  function MSRMatching(self, type) result(a)
    class (Running)    , intent(in) :: self
    character (len = *), intent(in) :: type
    real (dp) , dimension(4)        :: a, c
    real (dp) , dimension(5)        :: b

    a = 0

    if ( type(:4) == 'MSRp' ) then
      return
    else if ( type(:4) == 'MSRn' ) then
      c = MSbarDelta(self%nf, 1); a = MSbarDelta(self%nf, 0)
      b = self%andim%MatchingAlphaUp()
      call alphaReExpand( a, b(:4) );  a = c - a
    else if ( type(:5) == 'charm' ) then
      a = self%andim%MSRDelta('charm') - MSbarDelta(self%nf - 1, 1)
    end if

  end function MSRMatching

!ccccccccccccccc

   real (dp) function lambdaQCD(self, run)
    integer        , intent(in) :: run
    class (Running), intent(in) :: self
    real (dp)                   :: t

    t = - 2 * pi/self%beta(0)/self%alphaQCD(self%muLambda)

    lambdaQCD = self%muLambda * self%andim%Gfun(run, t)

   end function lambdaQCD

!ccccccccccccccc

  real (dp) function DiffR(self, sCoef, order, r0, r1)
    class (Running)        , intent(in) :: self
    integer                , intent(in) :: order
    real (dp)              , intent(in) :: r0, r1
    real (dp), dimension(3), intent(in) ::sCoef

    DiffR = 0; if ( abs(r0 - r1) <= d1mach(1) ) return

    DiffR = self%andim%DeltaR( sCoef, order, self%alphaQCD(r0), self%alphaQCD(r1) )

  end function DiffR

!ccccccccccccccc

  real (dp) function DiffRMass(self, order, m, r0, r1)
    integer        , intent(in) :: order
    real (dp)      , intent(in) :: r0, r1, m
    class (Running), intent(in) :: self
    real (dp)                   :: a0, a1

    DiffRMass = 0; if ( abs(r0 - r1) <= d1mach(1) .or. order < 2) return

    a0 = self%alphaQCD(r0); a1 = self%alphaQCD(r1)

    DiffRMass = self%andim%DeltaRMass( order, self%lambdaQCD(order), m, a0, a1 )

  end function DiffRMass

!ccccccccccccccc

  real (dp) function DiffDeltaHadron(self, gamma, order, r0, r1, mu0, mu1)
    integer                , intent(in) :: order
    real (dp)              , intent(in) :: r0, r1, mu0, mu1
    class (Running)        , intent(in) :: self
    real (dp), dimension(3), intent(in) :: gamma
    real (dp)                           :: a, b, delta, a0, a1, abserr
    real (dp), dimension(3)             :: gammaBet(3)
    integer                             :: neval, ier

    gammaBet = gamma * PowList(1/( 2 * self%beta(0) ), 3)

    a = 6/self%beta(0) * log( self%alphaQCD(mu1)/self%alphaQCD(r1) )
    b = 6/self%beta(0) * log( self%alphaQCD(mu1)/self%alphaQCD(r0) )

    a0 = self%alphaQCD(r0);  a1 = self%alphaQCD(r1)
    delta = 6/self%beta(0) * log( self%beta(0)/2/Pi * self%alphaQCD(mu1) )

    call qags( inteHadron, 0._dp, 1._dp, prec, prec, DiffDeltaHadron, abserr, neval, ier )

    DiffDeltaHadron = self%lambdaQCD(order) * DiffDeltaHadron +            &
    0.8862269254527579_dp * ( dgamma(1 + a)/dgamma(1.5_dp + a) *              &
    self%DiffDeltaMu(order, R1, R1, mu1) + dgamma(1 + b)/dgamma(1.5_dp + b) * &
    self%DiffDeltaMu(order, R0, mu0, R0) )

    contains

      real (dp) function inteHadron(r)
        real (dp), intent(in) :: r
        inteHadron = (1 - r**2)**delta * &
        self%andim%DeltaRHadron( self%andim%sCoefHadron(gammaBet,r), r, order, a0, a1)
      end

  end function DiffDeltaHadron

!ccccccccccccccc

  real (dp) function DiffDeltaMu(self, order, R, mu0, mu1)
    integer        , intent(in) :: order
    real (dp)      , intent(in) :: R, mu0, mu1
    class (Running), intent(in) :: self

    DiffDeltaMu = 0; if ( abs(mu0 - mu1) <= d1mach(1) ) return

    DiffDeltaMu = self%andim%DeltaMu( order, R, self%alphaQCD(mu0), self%alphaQCD(mu1) )

  end function DiffDeltaMu

!ccccccccccccccc

  real (dp) function DeltaGapMatching(self)
    class (Running), intent(in) :: self

    DeltaGapMatching = ExpEuler * self%mL * ( self%alphaQCD(self%mL)/Pi )**2 * deltaMass(0._dp, 1._dp)

  end function DeltaGapMatching

!ccccccccccccccc

  real (dp) function DiffDelta(self, sCoef, a, order, r0, r1, mu0, mu1)
    class (Running)        , intent(in) :: self
    integer                , intent(in) :: order
    real (dp)              , intent(in) :: a, r0, r1, mu0, mu1
    real (dp), dimension(3), intent(in) :: sCoef

    DiffDelta = a * ( self%DiffDeltaMu(order, R1, R1, mu1)   + &
                      self%DiffDeltaMu(order, R0, mu0, R0) ) + &
                      self%lambdaQCD(order) * self%DiffR(sCoef, order, r0, r1)

  end function DiffDelta

!ccccccccccccccc

  pure real (dp) function scales(self, str)
    class (Running)    , intent(in) :: self
    character (len = *), intent(in) :: str

    scales = 0

    if ( str(:2) == 'mZ'       ) scales = self%AlphaOb%scales('mZ')
    if ( str(:2) == 'mC'       ) scales = self%AlphaOb%scales('mC')
    if ( str(:2) == 'mB'       ) scales = self%AlphaOb%scales('mB')
    if ( str(:2) == 'mT'       ) scales = self%AlphaOb%scales('mT')
    if ( str(:2) == 'mH'       ) scales = self%mH
    if ( str(:2) == 'mL'       ) scales = self%mL
    if ( str(:3) == 'amZ'      ) scales = self%AlphaOb%scales('amZ')
    if ( str(:7) == 'muLambda' ) scales = self%muLambda

  end function scales

!ccccccccccccccc

  pure integer function orders(self, str)
    class (Running)    , intent(in) :: self
    character (len = *), intent(in) :: str

    orders = 0

    if ( str(:5) == 'order'    ) orders = self%AlphaOb%orders('order')
    if ( str(:8) == 'runAlpha' ) orders = self%AlphaOb%orders('run'  )
    if ( str(:7) == 'runMass'  ) orders = self%runMass

  end function orders

!ccccccccccccccc

  pure type(AnomDim) function adim(self)
    class (Running), intent(in) :: self
    adim = self%andim
  end function

!ccccccccccccccc

  pure type(Alpha) function AlphaAll(self)
    class (Running), intent(in) :: self
    AlphaAll = self%AlphaOb
  end function

!ccccccccccccccc

  pure function DeltaComputer(Mat, logList, pow) result(delta)
    real (dp), dimension(:  ), intent(in) :: logList
    real (dp), dimension(:,:), intent(in) :: Mat
    integer                  , intent(in) :: pow
    real (dp), dimension( size(Mat,2) )   :: delta
    integer                               :: i

    delta = Mat(1,:)

    do i = 1, min( size(logList), size(Mat,1) - 1 )
      if (pow == 1) delta = delta + (i + 1) * logList(i) * Mat(i + 1, :)
      if (pow == 0) delta = delta + logList(i) * Mat(i + 1, :)
    enddo

  end function DeltaComputer

!ccccccccccccccc

end module RunningClass
