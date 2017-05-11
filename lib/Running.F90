
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
    real (dp), dimension(4)     :: bHat, sCoef, sCoefNatural, gammaR, gammaRNatural
    real (dp)                   :: muLambda, mH, mL
    real (dp), dimension(0:4,4) :: tab

    contains

    procedure, pass(self), public :: MSbarMass, MSRMass, lambdaQCD, DiffDeltaMu,  &
    DiffDelta, orders, adim, DiffDeltaHadron, scheme, MSbarDeltaMu, MSbarMassLow, &
    AlphaAll, DeltaGapMatching, DiffRMass, MSRNaturalMass, mmFromMSR, PoleMass,   &
    SetMTop, SetMBottom, SetMCharm, SetLambda, mmFromMSRNatural, SetAlpha, scales,&
    DiffR, PSdelta, OptimalR, OptimalRNatural, numFlav

    procedure, pass(self), private :: MSRMatching, alphaQCDReal, alphaQCDComplex, &
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
    else if (nf < 4) then
      InitRun%mH = AlphaOb%scales('mC'); InitRun%mL = 0
    end if

    InitRun%tab = MSbarDeltaPiece(nf - 1, 1)
    sCoef = InitRun%andim%betaQCD('sCoefMSR')
    InitRun%sCoef = sCoef(1:); sCoef = 0;
    sCoef = InitRun%andim%betaQCD('sCoefMSRNatural')
    InitRun%sCoefNatural = sCoef(1:)
    sCoef = 0; sCoef = InitRun%andim%betaQCD('bHat')
    InitRun%bHat = sCoef(1:)
    sCoef = 0; sCoef = InitRun%andim%betaQCD('MSRdelta')
    InitRun%gammaR = sCoef(1:)
    sCoef = 0; sCoef = InitRun%andim%betaQCD('MSRNaturaldelta')
    InitRun%gammaRNatural = sCoef(1:)

   end function InitRun

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

    call self%AlphaOb%SetAlpha(alpha)

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
     LogList(0) = 1; LogList(1:) = PowList(  log( mu/self%MSbarMass(mu) ), 4  )

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

  real (dp) function OptimalR(self, n, lambda, method)
    class (Running)                    , intent(in) :: self
    real (dp)          , optional      , intent(in) :: lambda
    character (len = *), optional      , intent(in) :: method
    real (dp)                          , intent(in) :: n
    integer                                         :: IFLAG
    real (dp)                                       :: a, b, c

    a = 0.5_dp; b = self%mH

    call DFZERO(root, a, b, c, 1e-10_dp, 1e-10_dp, IFLAG); OptimalR = a

  contains

    real (dp) function root(x)
      real (dp), intent(in) :: x

      root = n * x - 4 * self%MSRmass(x, lambda, method) * self%alphaQCD(x)/3

    end function root

  end function OptimalR

!ccccccccccccccc

   real (dp) function MSRmass(self, R, lambda, method)
     class (Running)                    , intent(in) :: self
     real (dp)                          , intent(in) :: R
     real (dp)          , optional      , intent(in) :: lambda
     character (len = *), optional      , intent(in) :: method
     real (dp), dimension(self%runMass,self%runMass) :: gammaLog
     real (dp), dimension(self%runMass)              :: gamm, alphaList, lglist
     integer                                         :: neval, ier
     real (dp)                                       :: corr, abserr, lan, delta, &
     mu, Rstep, delta1, delta2

     if ( present(lambda) ) then

       gammaLog(1,:) = self%gammaR(:self%runMass) ; lan = lambda
       call self%andim%expandAlpha( gammaLog )

       gamm = lambda * (  gammaLog(1,:) + &
       matmul( powList(-log(lambda), self%runMass-1), gammaLog(2:,:) )  )

     else

       lan = 1

     end if

     if (  ( .not. present(method) ) .or. method(:8) == 'analytic' ) then

       if ( .not. present(lambda) ) then
         corr = self%DiffR( self%sCoef, self%runMass, self%mH, R )
       else
        corr = self%DiffR( self%andim%sCoefRecursive(gamm), &
        self%runMass, self%mH/lambda, R/lambda )
       end if

       corr = self%lambdaQCD(self%runMass) * corr

    else if ( method(:7) == 'numeric' ) then

      gammaLog(1,:) = self%gammaR(:self%runMass)

      if ( .not. present(lambda) ) then
        gamm = self%gammaR(:self%runMass)
      end if

        gamm = self%andim%GammaRComputer( gamm )

      call qags( Integrand, R/lan, self%mH/lan, prec, prec, corr, abserr, neval, ier )

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( ( self%mH - R )/delta )  )
      delta = (self%mH - R)/neval; corr = 0; Rstep = self%mH + delta
      gammaLog(1,:) = self%gammaR(:self%runMass)
      call self%andim%expandAlpha( gammaLog )

      do ier = 1, neval

        Rstep = Rstep - delta; mu = (Rstep - delta/2)/lan
        alphaList = powList( self%alphaQCD(mu)/Pi, self%runMass )

        lgList = powList( log(mu/Rstep), self%runMass )
        delta1 = Rstep * dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )

        lgList = powList( log(mu/(Rstep - delta) ), self%runMass )
        delta2 = (Rstep - delta) * dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )

        corr = corr + delta1 - delta2

      end do

    end if

    MSRmass = self%mH + corr

    contains

!ccccccccccccccc

   real (dp) function Integrand(R)
     real (dp)             , intent(in) :: R
     real (dp), dimension(self%runMass) :: aPi

     aPi = powList( self%alphaQCD(R)/4/Pi, self%runMass )

     Integrand = dot_product( gamm, aPi )

   end function Integrand

   end function MSRmass

!ccccccccccccccc

  real (dp) function mmFromMSR(self, mass, R)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: R, mass
    integer                        :: IFLAG
    real (dp)                      :: a, b, c, rat, massOr

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
     real (dp), intent(in) :: mm

     if (self%nf == 5) call self%SetMTop(   mm, rat * mm)
     if (self%nf == 4) call self%SetMBottom(mm, rat * mm)
     if (self%nf == 3) call self%SetMCharm (mm, rat * mm)

     MSR = self%MSRMass(R) - mass


   end function MSR

  end function mmFromMSR

!ccccccccccccccc

  real (dp) function mmFromMSRNatural(self, mass, order, R)
    class (Running), intent(inout) :: self
    real (dp)      , intent(in   ) :: R, mass
    integer        , intent(in   ) :: order
    integer                        :: IFLAG
    real (dp)                      :: a, b, c, rat, massOr

     a = mass/2; b = 2 * mass; massOr = self%mH

    if (self%nf == 5) rat = self%alphaOb%scales('muT')/self%mH
    if (self%nf == 4) rat = self%alphaOb%scales('muB')/self%mH
    if (self%nf == 3) rat = self%alphaOb%scales('muC')/self%mH

    call DFZERO(MSR, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG); mmFromMSRNatural = a

    if (self%nf == 5) call self%SetMTop(   massOr, rat * massOr)
    if (self%nf == 4) call self%SetMBottom(massOr, rat * massOr)
    if (self%nf == 3) call self%SetMCharm( massOr, rat * massOr)

   contains

   real (dp) function MSR(mm)
     real (dp), intent(in)   :: mm

     if (self%nf == 5) call self%SetMTop(   mm, rat * mm)
     if (self%nf == 4) call self%SetMBottom(mm, rat * mm)
     if (self%nf == 3) call self%SetMCharm (mm, rat * mm)

     MSR = self%MSRNaturalMass(order, R) - mass

   end function MSR

  end function mmFromMSRNatural

!ccccccccccccccc

  real (dp) function OptimalRNatural(self, n, order, lambda, method)
    class (Running)                    , intent(in) :: self
    real (dp)          , optional      , intent(in) :: lambda
    character (len = *), optional      , intent(in) :: method
    integer                            , intent(in) :: order
    real (dp)                          , intent(in) :: n
    integer                                         :: IFLAG
    real (dp)                                       :: a, b, c

    a = 0.5_dp; b = self%mH

    call DFZERO(root, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG); OptimalRNatural = a

  contains

    real (dp) function root(x)
      real (dp), intent(in) :: x

      root = n * x - 4 * self%MSRNaturalMass(order, x, lambda, method) &
      * self%alphaQCD(x)/3

    end function root

  end function OptimalRNatural

!ccccccccccccccc

   real (dp) function MSRNaturalMass(self, order, R, lambda, method)
     class (Running)                    , intent(in) :: self
     real (dp)                          , intent(in) :: R
     real (dp), optional                , intent(in) :: lambda
     character (len = *), optional      , intent(in) :: method
     integer                            , intent(in) :: order
     real (dp)                        , dimension(4) :: a
     real (dp), dimension(self%runMass,self%runMass) :: gammaLog
     real (dp), dimension(self%runMass)              :: gamm, alphaList, lglist
     integer                                         :: neval, ier, i
     real (dp)                                       :: corr, abserr, lan, delta, &
     alphaM, matching, mu, Rstep, delta1, delta2

     if ( present(lambda) ) then

       gammaLog(1,:) = self%gammaRNatural(:self%runMass) ; lan = lambda
       call self%andim%expandAlpha( gammaLog )

       gamm = lambda * (  gammaLog(1,:) + &
       matmul( powList(-log(lambda), self%runMass-1), gammaLog(2:,:) )  )

     else

       lan = 1

     end if

     if (  ( .not. present(method) ) .or. method(:8) == 'analytic' ) then

       if ( .not. present(lambda) ) then
         corr = self%DiffR( self%sCoefNatural, self%runMass, self%mH, R )
       else
        corr = self%DiffR( self%andim%sCoefRecursive(gamm), &
        self%runMass, self%mH/lambda, R/lambda )
       end if

       corr = self%lambdaQCD(self%runMass) * corr

    else if ( method(:7) == 'numeric' ) then

      gammaLog(1,:) = self%gammaRNatural(:self%runMass)

      if ( .not. present(lambda) ) then
        gamm = self%gammaRNatural(:self%runMass)
      end if

        gamm = self%andim%GammaRComputer( gamm )

      call qags( Integrand, R/lan, self%mH/lan, prec, prec, corr, abserr, neval, ier )

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( ( self%mH - R )/delta )  )
      delta = (self%mH - R)/neval; corr = 0; Rstep = self%mH + delta
      gammaLog(1,:) = self%gammaRNatural(:self%runMass)
      call self%andim%expandAlpha( gammaLog )

      do ier = 1, neval

        Rstep = Rstep - delta; mu = (Rstep - delta/2)/lan
        alphaList = powList( self%alphaQCD(mu)/Pi, self%runMass )

        lgList = powList( log(mu/Rstep), self%runMass )
        delta1 = Rstep * dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )

        lgList = powList( log(mu/(Rstep - delta) ), self%runMass )
        delta2 = (Rstep - delta) * dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )

        corr = corr + delta1 - delta2

      end do

    end if

     if (self%runMass > 0) alphaM = self%AlphaOb%alphaQCD(self%nf + 1, self%mH)/Pi ! note: perform matching with nl active flavors?
     a = self%MSRMatching(); i = min(order, 4)

     matching = 1 + dot_product( a(:i), PowList(alphaM, i) )

     MSRNaturalMass = self%mH * matching + corr

    contains

!ccccccccccccccc

   real (dp) function Integrand(R)
     real (dp)             , intent(in) :: R
     real (dp), dimension(self%runMass) :: aPi

     aPi = powList( self%alphaQCD(R)/4/Pi, self%runMass )

     Integrand = dot_product( gamm, aPi )

   end function Integrand

   end function MSRNaturalMass

!ccccccccccccccc

  function MSRMatching(self) result(a)
    class (Running), intent(in) :: self
    real (dp) , dimension(4)    :: a, c
    real (dp) , dimension(5)    :: b

    c = MSbarDelta(self%nf, 1); a = MSbarDelta(self%nf, 0)
    b = self%andim%MatchingAlphaUp()
    call alphaReExpand( a, b(:4) );  a = c - a

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
