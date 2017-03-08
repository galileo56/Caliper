
module RunningClass
  use AnomDimClass;  use AlphaClass; use Constants, only: dp, Pi, ExpEuler, d1mach
  use QuadPack, only: qags         ; implicit none   ;  private

!ccccccccccccccc

  type, public                  :: Running
    private
    character (len = 5)         :: str
    integer                     :: nf, runMass
    type (AnomDim)              :: andim
    type (Alpha)                :: AlphaOb
    real (dp), dimension(0:4)   :: beta, gamma
    real (dp), dimension(4)     :: bHat, sCoef, sCoefNatural
    real (dp)                   :: muLambda, mH, mL
    real (dp), dimension(0:4,4) :: tab

    contains

    procedure, pass(self), public :: MSbarMass, MSRMass, lambdaQCD, DiffDeltaMu,  &
    DiffDelta, orders, adim, DiffDeltaHadron, scheme, MSbarDeltaMu, MSbarMassLow, &
    AlphaAll, DeltaGapMatching, DiffRMass, MSRNaturalMass, mmFromMSR, PoleMass,   &
    SetMTop, SetMBottom, SetMCharm, SetLambda, mmFromMSRNatural, SetAlpha, scales,&
    DiffR

    procedure, pass(self), private :: MSRMatching, alphaQCDReal, alphaQCDComplex, &
    wTildeReal, wTildeComplex, kTildeReal, kTildeComplex, RunningMass, sCoefLambda

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

   real (dp) function MSRmass(self, R, lambda)
     class (Running)    , intent(in) :: self
     real (dp)          , intent(in) :: R
     real (dp), optional, intent(in) :: lambda
     real (dp)                       :: corr

     if ( .not. present(lambda) ) then
       corr = self%DiffR( self%sCoef, self%runMass, self%mH, R )
     else
       corr = self%DiffR( self%sCoefLambda('Practical', lambda), &
       self%runMass, self%mH/lambda, R/lambda )
     end if

     MSRmass = self%mH + self%lambdaQCD(self%runMass) * corr

   end function MSRmass

!ccccccccccccccc

  function sCoefLambda(self, type, lambda) result(res)
     class (Running)      , intent(in) :: self
     real (dp)            , intent(in) :: lambda
     character (len = *)  , intent(in) :: type
     real (dp)          , dimension(4) :: res, scoef
     real (dp)                         :: lg

     scoef = 0; lg = log(lambda)

     if ( type(:7) == 'Natural' ) then
       scoef = self%sCoefNatural
     else if ( type(:9) == 'Practical' ) then
       scoef = self%sCoef
     end if

     res = scoef

     res(2) = res(2) - lg * res(1)
     res(3) = res(3) - 2 * lg * res(2) + res(1) * (  lg**2 - ( 2 * self%bHat(1) &
     + self%bHat(2) ) * lg  )

     res(4) = res(4) - 2 * lg * res(3) + res(2) * (  lg**2 - ( 3 * self%bHat(1) &
     + self%bHat(2) ) * lg  ) + res(1) * (  (self%bHat(3) - 3 * self%bHat(1)**2 +&
     3 * self%bHat(2) - self%bHat(1) * self%bHat(2)) * lg + ( 9 * self%bHat(1)/2 &
     + 2 * self%bHat(2) ) * lg**2 - lg**3  )

     res = lambda * res

  end function sCoefLambda

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

   real (dp) function MSRNaturalMass(self, order, R, lambda)
     class (Running)    , intent(in) :: self
     real (dp)          , intent(in) :: R
     real (dp), optional, intent(in) :: lambda
     integer            , intent(in) :: order
     real (dp)        , dimension(3) :: a
     real (dp)                       :: alphaM, matching, corr
     integer                         :: i

     if ( .not. present(lambda) ) then
       corr = self%DiffR( self%sCoefNatural, self%runMass, self%mH, R )
     else
       corr = self%DiffR( self%sCoefLambda('Natural', lambda), &
       self%runMass, self%mH/lambda, R/lambda )
     end if

     if (self%runMass > 0) alphaM = self%AlphaOb%alphaQCD(self%nf + 1, self%mH)/Pi
     a = self%MSRMatching(); i = min(order, 4)

     matching = 1 + dot_product( a(:i), PowList(alphaM, i) )

     MSRNaturalMass = self%mH * matching + self%lambdaQCD(self%runMass) * corr

   end function MSRNaturalMass

!ccccccccccccccc

  function MSRMatching(self) result(a)
    class (Running), intent(in)    :: self
    real (dp) , dimension(4)       :: a, b, c
    real (dp) , dimension(0:3,0:3) :: tab

    c   = MSbarDelta(self%nf, 1); a = MSbarDelta(self%nf, 0)
    tab = self%andim%alphaMatching(self%nf + 1)

    b = tab(:,0); call alphaReExpand(a,b);  a = c - a

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
    real (dp)                           :: a, b, delta, a0, a1, abserr, tny
    real (dp), dimension(3)             :: gammaBet(3)
    integer                             :: neval, ier

    gammaBet = gamma * PowList(1/( 2 * self%beta(0) ), 3); tny = tiny(1._dp)

    a = 6/self%beta(0) * log( self%alphaQCD(mu1)/self%alphaQCD(r1) )
    b = 6/self%beta(0) * log( self%alphaQCD(mu1)/self%alphaQCD(r0) )

    a0 = self%alphaQCD(r0);  a1 = self%alphaQCD(r1)
    delta = 6/self%beta(0) * log( self%beta(0)/2/Pi * self%alphaQCD(mu1) )

    call qags( inteHadron, 0._dp, 1._dp, tny, tny, DiffDeltaHadron, abserr, neval, ier )

    DiffDeltaHadron = self%lambdaQCD(order) * DiffDeltaHadron +            &
    0.8862269254527579_dp * ( dgamma(1 + a)/dgamma(1.5 + a) *               &
    self%DiffDeltaMu(order, R1, R1, mu1) + dgamma(1 + b)/dgamma(1.5 + b) * &
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

end module RunningClass
