
module SingularClass
  use MatrixElementsClass; use AnomDimClass ; use RunningClass; use KernelsClass
  use ModelClass         ; use GapMassClass ; use MassiveNSClass
  use Hyper              ; use MCtopClass   ; use QuadPack, only: qagi, qags
  use constants, only: dp, Pi, Pi2, Zeta3, ExpEuler, d1mach, prec ; implicit none
  private

!ccccccccccccccc

  type, public                    :: SingularMassless
    private
    class (MatricesElements), allocatable :: MatEl
    type (Running)                :: alpha
    character (len = 8)           :: hard
    integer                       :: run, run1, nf
    Type (AnomDim)                :: andim
    character (len = 6)           :: shape, EShape
    real (dp), dimension(0:3)     :: beta, cusp, noncusp, HardExp, HardMassExp, ExpHard
    real (dp), dimension(3)       :: alphaList
    real (dp), dimension(0:4,3)   :: MatExp, MatHard
    real (dp), dimension(6)       :: NGL2loop, NGLog2loop
    real (dp), dimension(0:6,0:3) :: MatExpanded, MatExpandedHJM, MatAdded, MatHJMFO, MatFO
    real (dp)                     :: Prefact, muS, alphaS, alphaJ, w, R, tmin, compfact, &
                                     ESFac, alphaR, muSEuler, A0, Q, muH, k, PiResum, QoMuH
    contains

    final                         :: delete_object
    procedure, pass(self), public :: Prefactor, ExpMat, s2rho, s2theta, SetMat, SetHard,  &
                                     HJMSing, s2log, Resum, wTilde, HJMSing1D, SetRunning,&
                                     DoubleSing, runOrd, setAlpha, SetMTop, SetMBottom,   &
                                     SetMCharm, setGammaShift, ESmin
    procedure, pass(self)         :: NGLIntComputer, IntMassJet, SingleSingMod, SingleSingList

    generic, public               :: SingleSing => SingleSingMod, SingleSingList

  end type SingularMassless

!ccccccccccccccc

  type, extends (SingularMassless), public :: SingularScales
    private
    real (dp), dimension(0:4,3)  , private :: MatJet, MatSoft
    real (dp), dimension(0:3)    , private :: Soft, Jet

  end type SingularScales

!ccccccccccccccc

  type, extends (SingularMassless), public :: SingularMass
    private
    real (dp), dimension(3)        :: deltaM
    real (dp)                      :: width, lm, LQ, alphaM, mass, muM, muJ, m2, m
    class (MassiveNS), allocatable :: MassNS
    type (MCtop)                   :: MC

    contains

    procedure, pass(self)          :: SoftMatchingExp, JetNonDist, NonDistMod, SetMass, &
                                      SingleSingWidthMod, SingleSingWidthList, NonDistList

    procedure, pass(self), public  :: NSMass

    generic, public                :: NonDist         => NonDistMod        , NonDistList
    generic, public                :: SingleSingWidth => SingleSingWidthMod, SingleSingWidthList

  end type SingularMass

!ccccccccccccccc

  type, extends (SingularMass), public   :: SingularMassScales
    private
    real (dp), dimension(0:4,3), private :: MatMassHard, MatJet, MatSoft, &
                                            MatSoftNl, MatBJet
    real (dp), dimension(6)    , private :: NGLog2loopNl
    real (dp)                  , private :: Hmat, QoMpole
    type (Running)             , private :: alphaNl
    real (dp), dimension(0:3)  , private :: cuspNl, Soft, SoftNl, Jet, bJet, Hm

  contains

    procedure, pass(self)                :: SetEverything
    procedure, pass(self), public        :: setHardMass, NSScales

  end type SingularMassScales

!ccccccccccccccc

  interface SingularMassless
    module procedure InSing
  end interface SingularMassless

!ccccccccccccccc

  interface SingularScales
    module procedure InSingScales
  end interface SingularScales

!ccccccccccccccc

  interface SingularMass
    module procedure InSingMass
  end interface SingularMass

!ccccccccccccccc

  interface SingularMassScales
    module procedure InMassScales
  end interface SingularMassScales

  contains ! subroutines start here

!ccccccccccccccc

 subroutine delete_object(this)
   type (SingularMassless) :: this
     if ( allocated(this%MatEl ) ) deallocate(this%MatEl)
  end subroutine delete_object

!ccccccccccccccc

  type (SingularScales) function InSingScales(MatEl, run, shape, hard)
    type (MatricesElements)      , intent(in) :: MatEl
    integer                      , intent(in) :: run
    character (len = *)          , intent(in) :: shape
    character (len = *), optional, intent(in) :: hard

    type (AnomDim)                            :: andim
    type (Running)                            :: alpha

    InSingScales%HardMassExp = [1, 0, 0, 0]; InSingScales%run = run
    InSingScales%run1 = run - 1

    allocate( MatricesElements :: InSingScales%MatEl ); InSingScales%Shape = Shape

    select type (selector => InSingScales%MatEl)
      type is (MatricesElements);  selector = MatEl
    end select

    andim = MatEl%adim('nf'); alpha = MatEl%run() ;  InSingScales%andim = andim
    InSingScales%beta = andim%betaQCD('beta')     ;  InSingScales%alpha = alpha

    InSingScales%cusp = andim%betaQCD('cusp') ; InSingScales%Jet  = andim%betaQCD('jet')
    InSingScales%Soft = andim%betaQCD('soft') ; InSingScales%tmin = 0
    InSingScales%NonCusp = andim%betaQCD('hard'); InSingScales%A0 = 1

    InSingScales%NGLog2loop = MatEl%NonGlobalList(); InSingScales%nf = andim%numFlav()

    InSingScales%hard = 'expand'; if ( present(hard) ) InSingScales%hard = hard

    InSingScales%MatSoft = MatEl%CoefMat(shape); InSingScales%ESFac = shapeFactor(shape)
    InSingScales%MatJet  = 2 * MatEl%CoefMat('jet')

  end function InSingScales

!ccccccccccccccc

  type (SingularMassless) function InSing(MatEl, run, mu, shape, hard)
    type (MatrixElements)        , intent(in) :: MatEl
    integer                      , intent(in) :: run
    real (dp)                    , intent(in) :: mu
    character (len = *)          , intent(in) :: shape
    character (len = *), optional, intent(in) :: hard

    real (dp)                                 :: Q, muS, muJ, muH, k, wJet, QoMuH
    real (dp), dimension(0:3)                 :: cusp, NonCusp, soft, jet
    complex (dp)                              :: muHC
    integer                                   :: i, j, run1
    type (AnomDim)                            :: andim
    type (Running)                            :: alpha

    InSing%HardMassExp = [1, 0, 0, 0];  allocate( MatrixElements :: InSing%MatEl )

    select type (selector => InSing%MatEl)
      type is (MatrixElements);  selector = MatEl
    end select

    InSing%Shape  = Shape;    andim = MatEl%adim('nf');  InSing%MatAdded = 0
    InSing%alphaS = MatEl%alphaScale('soft'); alpha = MatEl%run() ; InSing%tmin = 0
    InSing%alphaJ = MatEl%alphaScale('jet') ; InSing%andim = andim; InSing%MatExp = 0
    InSing%alphaR = MatEl%alphaScale('R')   ; InSing%beta  = andim%betaQCD('beta')

    muS = MatEl%scales('muS'); muJ = MatEl%scales('muJ'); InSing%R = MatEl%scales('R')
    Q   = MatEl%scales('Q')  ; muH = abs( MatEl%scales('muH') ); InSing%run = run
    QoMuH = Q/muH; soft = andim%betaQCD('soft'); jet = andim%betaQCD('jet')
    InSing%run1 = run - 1

    InSing%muS = muS ;  InSing%Q = Q ;  InSing%NGL2loop = MatEl%NonGlobalList()

    InSing%nf = andim%numFlav(); InSing%A0 = 1; InSing%preFact = 1; run1 = run - 1

    if ( .not. present(hard) .or. hard(:6) == 'expand' ) then
      InSing%MatExp  = MatEl%CoefMat('hard');  InSing%HardExp = [1, 0, 0, 0]
    else
      InSing%HardExp = expandExpVec( MatEl%delta('hard') )
    end if

    InSing%MatExp = InSing%MatExp + 2 * MatEl%CoefMat('jet') + MatEl%CoefMat(shape)

    InSing%ESFac  = shapeFactor(shape); cusp = andim%betaQCD('cusp')
    NonCusp   = andim%betaQCD('hard')

    wJet = 2 * alpha%wTilde( run, cusp, muJ, mu )

    InSing%w = - 2 * alpha%wTilde( run, cusp, muS, mu) + wJet

    k = - 2 * alpha%kTilde(run, cusp, muS, mu) - 2 * alpha%kTilde(run , cusp, muH, mu) &
        + 4 * alpha%kTilde(run, cusp, muJ, mu) +     alpha%wTilde(run1, soft, muS, mu) &
        +     alpha%wTilde(run1, jet, muJ, mu) +     alpha%wTilde(run1, NonCusp, muH, mu)

    if (MatEl%scales('muH') < 0) then

      muHC = dcmplx(0._dp, - muH)

      k = k + alpha%wTilde(run1, NonCusp, muHC, muH) - 2 * alpha%kTilde(run, cusp, muHC, muH)

      InSing%preFact = InSing%preFact * QoMuH**( 2 * alpha%wTilde(run, cusp, muHC, muH ) )

    end if

    InSing%preFact =  InSing%preFact * (muJ**2/Q/muS)**wJet * exp(k) * QoMuH**( &
      2 * alpha%wTilde(run, cusp, muH, mu) )

    InSing%MatExpanded    = expandExp(InSing%MatExp  )
    InSing%MatExpandedHJM = expandExp(InSing%MatExp/2)

    do i = 0, 3
      do j = 0, 2 * i
        InSing%MatAdded(j,i) = Sum(  InSing%MatExpanded( j, ceiling(j/2.):i )  )
      end do
    end do

    InSing%muSEuler = muS * ExpEuler;   InSing%compfact = 1/InSing%muSEuler

    InSing%MatFO    = matmul( PosToMomMatrix(), InSing%MatAdded       )
    InSing%MatHJMFO = matmul( PosToMomMatrix(), InSing%MatExpandedHJM )

  end function InSing

!ccccccccccccccc

  type (SingularMass) function InSingMass(MassNS, run, xiJ, xiB, mu, hard)
    type (MassiveNS)             , intent(in) :: MassNS
    real (dp)                    , intent(in) :: mu, xiJ, xiB
    integer                      , intent(in) :: run
    character (len = *), optional, intent(in) :: hard
    character (len = 6)                       :: shape, EShape
    real (dp), dimension(0:3)                 :: cusp, cuspNl, jet, bJet, NonCusp, Hmass, &
                                                 Soft, SoftNl
    type (MatrixElementsMass)                 :: MatEl
    type (AnomDim)                            :: andim, andimNl
    type (Running)                            :: alpha, alphaNl
    integer                                   :: i, j, run1
    complex (dp)                              :: muHC
    real    (dp)                              :: Q, muS, muJ, muH, muM, k, wJet, wSoft, &
                                                 mPole, kSoft, wSoftnf, wHard, kHard  , &
                                                 kJet, wSoftNl, mass, wJetNf, wJetNl  , &
                                                 QoMu, QoM, JetArg, alphaM, m2, A0, wB

    mPole   = MassNS%MassVar('mPole')   ; InSingMass%deltaM = MassNS%MassDeltaShift()
    shape   = MassNS%EShape('shape')    ; InSingMass%shape  = shape; run1 = run - 1
    EShape  = MassNS%EShape('EShape')   ; InSingMass%mass   = mPole
    InSingMass%EShape = EShape          ; MatEl = MassNS%matElementNS()
    InSingMass%nf = MassNS%numFlav()    ; andimNl = MatEl%adim('nl')
    InSingMass%R = MatEl%scales('R')    ; A0 = MassNS%MassVar('A0'); InSingMass%A0 = A0

    allocate( MatrixElementsMass :: InSingMass%MatEl )

    select type (selector => InSingMass%MatEl)
      type is (MatrixElementsMass);  selector = MatEl
    end select

    andim = MatEl%adim('nf'); InSingMass%HardMassExp = [1, 0, 0, 0]
    InSingMass%alphaS = MatEl%alphaScale('soft');  alpha = MatEl%run()
    InSingMass%alphaJ = MatEl%alphaScale('jet') ;  InSingMass%andim = andim
    InSingMass%alphaR = MatEl%alphaScale('R')   ;  alphaNl = MatEl%runNl()
    InSingMass%beta   = andimNl%betaQCD('beta') ;  InSingMass%MatAdded = 0

    muS = MatEl%scales('muS')  ; muJ = MatEl%scales('muJ'); muH = abs( MatEl%scales('muH') )
    Q   = MatEl%scales('Q')    ; InSingMass%run = run     ; muM = MatEl%scales('muM')
    InSingMass%MatExp = 0      ; mass    = MassNS%MassVar('mass')     ; QoMu = Q/muH
    JetArg = muJ**2/Q/muS      ; alphaM  = matEl%alphaScale('massNf') ; QoM  = Q/mPole
    InSingMass%alphaM = alphaM ; InSingMass%muJ = muJ; InSingMass%run1 = run1

    InSingMass%MC = MCtop(shape, mass, Q)

    if (muJ < muM) then
      InSingMass%MatExp = MatEl%CoefMat('bJet')
      call MatEl%AddLogs( InSingMass%MatExp, log(muJ * mass/Q/muS) )
      InSingMass%MatExp = 2 * ( InSingMass%MatExp + MassNS%JetMassExp(InSingMass%hard, xiJ, xiB)/A0 )
    else
      InSingMass%MatExp = 2 * ( MatEl%CoefMat('jet') + MassNS%JetMassExp(InSingMass%hard, xiJ)/A0 )
    end if

    InSingMass%muS    = muS   ;  InSingMass%Q = Q ;  wSoftnf = 0; InSingMass%muM = muM
    andimNl = MatEl%adim('nl');  InSingMass%NGL2loop = MatEl%NonGlobalList()

    if ( .not. present(hard) .or. hard(:6) == 'expand') then
      InSingMass%MatExp = InSingMass%MatExp + MatEl%CoefMat('hard')
      if (muJ < muM) then
        InSingMass%MatExp = InSingMass%MatExp + MatEl%HardMassExp(mPole)
      end if
      InSingMass%HardExp = [1, 0, 0, 0]
    else
      InSingMass%HardExp = expandExpVec( MatEl%delta('hard') )
      InSingMass%HardExp(1) = InSingMass%HardExp(1) + (1 - xiJ) * MassNS%MassVar('Hcorr') &
                                                               * MassNS%MassVar('alphaH')
      if (muJ < muM) then
        InSingMass%HardMassExp = expandExpVec( MatEl%HardMassVec(mPole) )
        InSingMass%HardMassExp(1) = InSingMass%HardMassExp(1) + xiJ * (1 - xiB) * &
                                    MassNS%MassVar('Hcorr') * MassNS%MassVar('alphaM')
      end if
    end if

    InSingMass%MatExp = InSingMass%MatExp + MatEl%CoefMat(shape)

    InSingMass%ESFac  = shapeFactor(shape)

!   First take care of the Jet sector

    cusp    = andim%betaQCD('cusp'); cuspNl = andimNl%betaQCD('cusp')
    jet     = andim%betaQCD('jet') ; bJet   = andimNl%betaQCD('bJet')
    NonCusp = andim%betaQCD('hard'); Hmass  = andimNl%betaQCD('Hm')
    Soft    = andim%betaQCD('soft'); SoftNl = andimNl%betaQCD('soft')

    jet_scenario: if (muJ > muM) then  ! SCET scenarios III & IV
      ordering_if: if (mu > muM) then  ! regular mu choice

        wJet  = 2 * alpha%wTilde(run, cusp, muJ, mu)
        wHard = 2 * alpha%wTilde(run, cusp, muH, mu)

        kJet = 4 * alpha%kTilde(run, cusp, muJ, mu) + alpha%wTilde(run1, jet , muJ, mu)

        kHard = alpha%wTilde(run1, NonCusp, muH, mu) - 2 * alpha%kTilde(run, cusp, muH, mu)

        InSingMass%preFact = QoMu**wHard * JetArg**wJet

      else  ! anti-natural mu choice

        wJet  = 2 *   alpha%wTilde(run, cusp  , muJ, muM)
        wB    = 2 * alphaNl%wTilde(run, cuspNl, muM, mu )
        wHard = 2 *   alpha%wTilde(run, cusp  , muH, muM)

        kJet = 4 * alpha%kTilde(run, cusp, muJ, muM) + alpha%wTilde(run1, jet, muJ, muM) &
         + 2 * alphaNl%kTilde(run, cuspNl, muM, mu) + alphaNl%wTilde(run1, bJet, muM, mu)

        kHard = alpha%wTilde(run1, NonCusp, muH, muM)      &
          - 2 * alpha%kTilde(run, cusp, muH, muM) + alphaNl%wTilde(run1, Hmass, muM, mu )

        InSingMass%preFact = QoMu**wHard * JetArg**wJet * (muM/QoM/muS)**wB *  &
                 QoM**(  2 * alphaNl%wTilde(run, cuspNl, muM, mu  )  )

        wJet = wJet + wB

      end if ordering_if

    else ! Scenario II

      ordering_if_2: if (mu < muM) then ! sensible ordering

        wJet  =   2 * alphaNl%wTilde(run, cuspNl, muJ, mu)
        wSoft = - 2 * alphaNl%wTilde(run, cuspNl, muS, mu)

        kJet  = 2 * alphaNl%kTilde(run, cuspNl, muJ, mu) + alphaNl%wTilde(run1, bJet, muJ, mu)

        kSoft = - 2 * alphaNl%kTilde(run , cuspNl, muS, mu) +  &
                      alphaNl%wTilde(run1, softNl, muS, mu)

        kHard = alpha%wTilde(run1, NonCusp, muH, muM) - &
         2 * alpha%kTilde(run , cusp, muH, muM) + alphaNl%wTilde(run1, Hmass, muM, mu)


        InSingMass%preFact = (muJ/QoM/muS)**wJet &
      * QoMu**( 2 * alpha%wTilde(run, cusp, muH, muM) )  &
      *  QoM**( 2 * alphaNl%wTilde(run, cuspNl, muM, mu) )

      else  ! nonsensical ordering

        wJetNf = 2 *   alpha%wTilde(run,   cusp, muM, mu )
        wJetNl = 2 * alphaNl%wTilde(run, cuspNl, muJ, muM)

        wJet =  wJetNf + wJetNl

        kJet  = 2 * alphaNl%kTilde(run , cuspNl, muJ, muM)      &
                +   alphaNl%wTilde(run1, bJet  , muJ, muM)      &
                + 4 * alpha%kTilde(run , cusp  , muM, mu )      &
                +     alpha%wTilde(run1, jet   , muM, mu )

        wSoftnf = - 2 *   alpha%wTilde(run,   cusp, muM, mu )
        wSoft   = - 2 * alphaNl%wTilde(run, cuspNl, muS, muM) + wSoftnf

        kSoft = - 2 * alphaNl%kTilde(run , cuspNl, muS, muM)  + &
                      alphaNl%wTilde(run1, SoftNl, muS, muM)  - &
                  2 *   alpha%kTilde(run , cusp  , muM, mu )  + &
                        alpha%wTilde(run1, soft  , muM, mu )

        kHard =       alpha%wTilde(run1, NonCusp, muH, mu)    &
                - 2 * alpha%kTilde(run , cusp, muH, mu)

        InSingMass%preFact = (muM/muS)**wSoftnf * (muJ/QoM/muS)**wJetNl * &
                               (muM**2/Q/muS)**wJetNf * QoMu**(  2      * &
                               alpha%wTilde(run, cusp, muH, mu )  )
      end if ordering_if_2
    end if  jet_scenario

!   We take care of the soft sector next

    soft_scenario: if (muM < muS) then ! scenario IV

      ordering_if_3: if (mu > muM) then ! sensible ordering

        wSoft = - 2 * alpha%wTilde(run , cusp, muS, mu)
        kSoft = - 2 * alpha%kTilde(run , cusp, muS, mu) + alpha%wTilde(run1, soft, muS, mu)

      else          ! not sensible ordering

        wSoftNl = - 2 * alphaNl%wTilde(run, cuspNl, muM, mu )
        wSoft   = - 2 *   alpha%wTilde(run,   cusp, muS, muM) + wSoftNl

        kSoft = - 2 * alpha%kTilde(run , cusp, muS, muM) + &
        alpha%wTilde(run1, soft, muS, muM) - 2 * alphaNl%kTilde(run, cuspNl, muM, mu ) + &
                    alphaNl%wTilde(run1, softNl, muM, mu )

        InSingMass%preFact = InSingMass%preFact * (muM/muS)**wSoftnl

      end if ordering_if_3

    else if (muJ >= muM) then ! scenario III

       InSingMass%MatExp = InSingMass%MatExp + InSingMass%SoftMatchingExp()

      ordering_if_4: if (mu > muM) then ! sensible ordering

        wSoftnf = - 2 *   alpha%wTilde(run,   cusp, muM, mu )
        wSoft   = - 2 * alphaNl%wTilde(run, cuspNl, muS, muM) + wSoftnf

        kSoft = - 2 * ( alphaNl%kTilde(run, cuspNl, muS, muM)    + &
                          alpha%kTilde(run,   cusp, muM, mu ) ) + &
        alphaNl%wTilde(run1, SoftNl, muS, muM) + alpha%wTilde(run1, Soft, muM, mu )

        InSingMass%preFact = InSingMass%preFact * (muM/muS)**wSoftnf

      else          ! nonsensical ordering

        wSoft = - 2 * alphaNl%wTilde(run, cuspNl, muS, mu)

        kSoft = - 2 * alphaNl%kTilde(run , cuspNl, muS, mu)  + &
                      alphaNl%wTilde(run1, SoftNl, muS, mu)

      end if ordering_if_4
    end if   soft_scenario

    m2 = MassNS%MassVar('m2') ; InSingMass%width = MassNS%MassVar('width')
    InSingMass%m2 = m2        ; InSingMass%Lm    = 2 * MassNS%MassVar('lm')

    InSingMass%tmin = MassNS%MassVar('tmin');  InSingMass%w = wSoft + wJet

    k = kSoft + kHard + kJet

    if ( MatEl%scales('muH') < 0 ) then

      muHC = dcmplx( 0._dp, - muH )

      k = k + alpha%wTilde(run1, NonCusp, muHC, muH) &
        - 2 * alpha%kTilde(run , cusp, muHC, muH)

     InSingMass%preFact = InSingMass%preFact * QoMu**( &
     2 * alpha%wTilde(run, cusp, muHC, muH )  )

    end if

    InSingMass%preFact = InSingMass%preFact * exp(k)

    InSingMass%MatExpanded    = InSingMass%A0 * expandExp(InSingMass%MatExp  )
    InSingMass%MatExpandedHJM = InSingMass%A0 * expandExp(InSingMass%MatExp/2)

    do i = 0, 3
      do j = 0, 2 * i
        InSingMass%MatAdded(j,i) = Sum(  InSingMass%MatExpanded( j, ceiling(j/2.):i )  )
      end do
    end do

    InSingMass%muSEuler = muS * ExpEuler; InSingMass%compfact = 1/InSingMass%muSEuler

    InSingMass%MatFO    = matmul( PosToMomMatrix(), InSingMass%MatAdded       )
    InSingMass%MatHJMFO = matmul( PosToMomMatrix(), InSingMass%MatExpandedHJM )

  end function InSingMass

!ccccccccccccccc

  type (SingularMassScales) function InMassScales(MassNS, run, hard)
    type (MassiveScales)         , intent(in) :: MassNS
    integer                      , intent(in) :: run
    character (len = *), optional, intent(in) :: hard

    character (len = 6)                       :: shape, EShape
    type (MatricesElementsMass)               :: MatEl
    type (AnomDim)                            :: andim, andimNl
    type (Running)                            :: alpha, alphaNl

    shape   = MassNS%EShape('shape')    ; InMassScales%shape = shape
    EShape  = MassNS%EShape('EShape')   ; InMassScales%mass  = MassNS%MassVar('mPole')
    InMassScales%EShape = EShape        ; MatEl   = MassNS%matElementScales()
    InMassScales%run1   = run - 1       ; andimNl = MatEl%adim('nl')

    allocate ( MassiveScales :: InMassScales%MassNS )

    select type (selector => InMassScales%MassNS)
    type is (MassiveScales)
      selector = MassNS
    end select

    allocate( MatricesElementsMass :: InMassScales%MatEl )

    select type (selector => InMassScales%MatEl)
     type is (MatricesElementsMass);  selector = MatEl
    end select

    andim = MatEl%adim('nf');  InMassScales%nf = MassNS%numFlav(); InMassScales%run = run
    alpha = MatEl%run();  InMassScales%andim = andim ;    alphaNl = MatEl%runNl()
    InMassScales%beta    = andimNl%betaQCD('cusp'); InMassScales%alpha = alpha
    InMassScales%alphaNl = alphaNl;  InMassScales%HardMassExp = [1, 0, 0, 0]

    InMassScales%cusp    = andim%betaQCD('cusp')
    InMassScales%SoftNl  = andimNl%betaQCD('soft')
    InMassScales%cuspNl  = andimNl%betaQCD('cusp')
    InMassScales%Soft    = andim%betaQCD('soft')
    InMassScales%Jet     = andim%betaQCD('jet')
    InMassScales%bJet    = andimNl%betaQCD('bJet')
    InMassScales%Hm      = andimNl%betaQCD('Hm')
    InMassScales%NonCusp = andim%betaQCD('hard')

    InMassScales%hard = 'expand'; if ( present(hard) ) InMassScales%hard = hard

    InMassScales%MatJet    = 2 * MatEl%CoefMat('jet')
    InMassScales%MatBJet   = 2 * MatEl%CoefMat('bJet')
    InMassScales%MatSoft   = MatEl%CoefMat(shape)
    InMassScales%MatSoftNl = MatEl%CoefMat('Nl'//shape)
    InMassScales%ESFac     = shapeFactor(shape)

    andimNl = MatEl%adim('nl') ;  InMassScales%NGLog2loop = MatEl%NonGlobalList()
    InMassScales%NGLog2loopNl = MatEl%NonGlobalListNl()

  end function InMassScales

!ccccccccccccccc

  subroutine SetMass(self)
    class (SingularMass), intent(inout) :: self

    call self%massNS%setMasses(); self%mass = self%MassNS%MassVar('mPole')

  end subroutine SetMass

!ccccccccccccccc

  subroutine SetAlpha(self, alpha)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: alpha

    call self%alpha%setAlpha(alpha)

    select type (self)
    type is (SingularMass); call self%massNS%setAlpha(alpha)
    end select

  end subroutine SetAlpha

!ccccccccccccccc

  subroutine SetMtop(self, m, mu)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: m, mu

    call self%alpha%SetMtop(m, mu)

    select type (self)
    type is (SingularMass); call self%massNS%SetMtop(m, mu)
    end select

  end subroutine SetMtop

!ccccccccccccccc

  subroutine SetMBottom(self, m, mu)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: m, mu

    call self%alpha%SetMBottom(m, mu)

    select type (self)
    type is (SingularMass); call self%massNS%SetMBottom(m, mu)
    end select

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, m, mu)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: m, mu

    call self%alpha%SetMCharm(m, mu)

    select type (self)
    type is (SingularMass); call self%massNS%SetMCharm(m, mu)
    end select

  end subroutine SetMCharm

!ccccccccccccccc

  subroutine setHard(self, Q, muH)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: Q, muH

    complex (dp)                            :: muHC

    self%Q = Q; self%muH = abs(muH);  self%QoMuH = Q/self%muH
    call self%MatEl%setHard(Q, muH)

      if ( self%hard(:6) == 'expand' ) then
        self%MatHard = self%MatEl%CoefMat('hard');  self%ExpHard = [1, 0, 0, 0]
      else
        self%ExpHard = expandExpVec( self%MatEl%delta('hard') );  self%MatHard = 0
      end if

      if (muH < 0) then

        muHC = dcmplx(0._dp, muH)

        self%k = self%alpha%wTilde( self%run1, self%NonCusp, muHC, self%muH ) &
           - 2 * self%alpha%kTilde( self%run , self%cusp   , muHC, self%muH )

        self%PiResum = self%QoMuH**( 2 * self%alpha%wTilde(self%run, self%cusp, muHC, self%muH )  )

      else
        self%PiResum = 1; self%k = 0
      end if

  end subroutine setHard

!ccccccccccccccc

  subroutine setHardMass(self, muM)
    class (SingularMassScales), intent(inout) :: self
    real (dp)                 , intent(in)    :: muM

    select type (selector => self%MatEl)
    class is (MatricesElementsMass)
      call selector%SetHardMass(muM)
    end select

    self%muM     = muM             ; self%Lm = 2 * self%MassNS%MassVar('lm')
    self%QoMpole = self%Q/self%mass; self%alphaM = self%matEl%alphaScale('massNf')

    select type (selector => self%MatEl)
    class is (MatricesElementsMass)

      if ( self%hard(:6) == 'expand' ) then
        self%MatMassHard = selector%HardMassExp(self%mass)
      else
        self%MatMassHard(:3,1) = expandExpVec( selector%HardMassVec(self%mass) )
        self%Hmat = 0; self%MatHard = 0
      end if

    end select

  end subroutine setHardMass

!ccccccccccccccc

  real (dp) function NonDistMod(self, Mod, setup, gap, space, cum, order, R0, &
                                mu0, delta0, h, tau, tau2)
    class (SingularMass) , intent(in) :: self
    type (Model)         , intent(in) :: Mod
    character (len = *)  , intent(in) :: setup, space, gap, cum
    real (dp)            , intent(in) :: R0, mu0, delta0, tau, h
    real (dp), optional  , intent(in) :: tau2
    integer              , intent(in) :: order

    real (dp)                         :: p, p2, Qm2, shift, w
    real (dp), dimension(1)           :: res
    real (dp), dimension(order,order) :: delta
    logical                           :: dobsing

    NonDistMod = 0

    if (present (tau2) ) then
      if (tau >= tau2) return
    end if

    if ( setup(:5) == 'Model' ) then
      if ( present(tau2) ) then
        res = self%NonDistList([Mod], gap, space, cum, order, R0, mu0, delta0, h, tau, tau2)
      else
        res = self%NonDistList([Mod], gap, space, cum, order, R0, mu0, delta0, h, tau)
      end if
      NonDistMod = res(1); return
    end if

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    p = (self%Q * tau - self%Q * self%tmin)/self%ESFac - shift; p2 = p
    if ( present(tau2) ) p2 = (self%Q * tau2 - self%Q * self%tmin)/self%ESFac - shift

    if ( p <= 0 .and. present(tau2) ) p = p2

    if ( p2 <= 0 .or. order <= 0 ) return

    dobsing = .not. present(tau2) .or. ( present(tau2) .and. abs(p2 - p) <= d1mach(1) )

    Qm2 = self%Q * self%m2; w = self%w + cumConst(cum)

    if ( setup(:2) == 'FO' ) then

      if ( dobsing ) then
        NonDistMod = 2 * self%alphaJ * self%JetNonDist(cum, self%Q * p)
      else
        NonDistMod = 2 * self%alphaJ * ( self%JetNonDist(cum, self%Q * p2) - &
                                      self%JetNonDist(cum, self%Q * p ) )
      end if

    else

      if ( dobsing ) then
        NonDistMod = 2 * self%muSEuler**w * self%alphaJ/3 * &
                   NonDist(p)/gamma(1 - w)/self%m2/self%Q
      else
        NonDistMod = 2 * self%muSEuler**w * self%alphaJ/3/self%m2/self%Q * &
                   ( NonDist(p2) - NonDist(p) )/gamma(1 - w)
      end if

    end if

    NonDistMod = Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) ) * &
    NonDistMod * self%Prefact * self%compFact**cumConst(cum) * &
                 (self%Q/self%ESFac)**( 1 + cumConst(cum) )

    contains

!ccccccccccccccc

  real (dp) function NonDist(x)
    real (dp), intent(in) :: x

    NonDist = self%IntMassJet(x/Qm2, w)/x**w

  end function NonDist

  end function NonDistMod

!ccccccccccccccc

  function NonDistList(self, ModList, gap, space, cum, order, R0, mu0, delta0, &
                       h, tau, tau2) result(resList)
    class (SingularMass)      , intent(in) :: self
    type (Model), dimension(:), intent(in) :: ModList
    character (len = *)       , intent(in) :: gap, cum, space
    real (dp)                 , intent(in) :: R0, mu0, delta0, tau, h
    real (dp), optional       , intent(in) :: tau2
    integer                   , intent(in) :: order
    real (dp)                              :: p, p2, Qm2, shift, w, absErr, gam, &
                                              newTerm, wProd
    real (dp), dimension(order,order)      :: delta
    real (dp), dimension( size(ModList) )  :: resList, Omega
    integer                                :: neval, ier, i, l, fac
    logical                                :: dobsing

    resList = 0; w = self%w + cumConst(cum); wProd = 1

    if ( present (tau2) ) then
      if (tau >= tau2) return
    end if

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    p = (self%Q * tau - self%Q * self%tmin)/self%ESFac - shift; p2 = p
    if ( present(tau2) ) p2 = (self%Q * tau2 - self%Q * self%tmin)/self%ESFac - shift

    if ( p <= 0 .and. present(tau2) ) p = p2

    if ( p2 <= 0 .or. order <= 0 ) return

    dobsing = .not. present(tau2) .or. ( present(tau2) .and. abs(p2 - p) <= d1mach(1) )

    Qm2 = self%Q * self%m2; fac = 1; gam = self%muSEuler**w / gamma(1 - w)

    if ( space(:3) == 'OPE' ) then

      if ( dobsing ) then
        resList = NonDist(p)
      else
        resList = NonDist(p, p2)
      end if

      resList = resList * [  ( ModList(i)%MomentModel(0), i = 1, size(modList) )  ]

      do l = 1, 30
        wProd = - wProd * w; w = w + 1; fac = fac * l
        Omega = [  ( ModList(i)%MomentModel(l), i = 1, size(modList) )  ]

        if ( dobsing ) then
          newTerm = NonDist(p)/fac
        else
          newTerm = NonDist(p, p2)/fac
        end if

        newTerm = newTerm * wProd * (-1)**l

        if ( abs( newTerm * maxval(Omega) ) > abs( maxval(resList) ) ) then
          resList = 0; w = self%w + cumConst(cum);  go to 10
        end if

        resList = resList + newTerm * Omega
        if ( abs( maxval(newTerm * Omega/resList) ) < prec ) go to 20

      end do

      go to 20

    end if

  10 continue

      do i = 1, size(modList)
        call qags( IntegrandNonDist, 0._dp, p2, prec, prec, resList(i), abserr, neval, ier )
      end do

  20 continue

    resList = 2 * gam * self%alphaJ/3 * resList * self%Prefact * &
    self%compFact**cumConst(cum) * (self%Q/self%ESFac)**( 1 + cumConst(cum) )/&
    self%m2/self%Q * Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) )

    contains

!ccccccccccccccc

  real (dp) function IntegrandNonDist(x)
    real (dp), intent(in) :: x

    if ( dobsing ) then
      IntegrandNonDist = NonDist(x) * ModList(i)%ShapeFun(0, p - x)
    else
      IntegrandNonDist = NonDist(x) * ModList(i)%ShapeFun(0, p - x, p2 - x)
    end if

  end function IntegrandNonDist

!ccccccccccccccc

    real (dp) function IntegrandNonDistCum(x)
      real (dp), intent(in) :: x

      if ( dobsing ) then
        IntegrandNonDistCum = NonDist(x) * ModList(i)%ShapeFun(-1, p - x)
      else
        IntegrandNonDistCum = NonDist(x) * ModList(i)%ShapeFun(-1, p - x, p2 - x)
      end if

    end function IntegrandNonDistCum

!ccccccccccccccc

  real (dp) function NonDist(x, x2)
    real (dp)          , intent(in) :: x
    real (dp), optional, intent(in) :: x2

      NonDist = self%IntMassJet(x/Qm2, w)/x**w

    if ( present(x2) ) NonDist = self%IntMassJet(x2/Qm2, w)/x2**w - NonDist

  end function NonDist

  end function NonDistList

!ccccccccccccccc

  real (dp) function ESmin(self, order, gap, delta0, R0, mu0, h)
    class (SingularMassless)         , intent(in ) :: self
    integer                          , intent(in ) :: order
    character (len = *)              , intent(in ) :: gap
    real (dp)                        , intent(in ) :: R0, mu0, delta0, h
    real (dp)                                      :: shift
    real (dp)             , dimension(order,order) :: delta

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    ESmin = shift * self%ESFac/self%Q + self%tmin

  end

!ccccccccccccccc

  subroutine setGammaShift( self, order, gap, delta0, R0, mu0, h, shift, delta )
    class (SingularMassless)         , intent(in ) :: self
    integer                          , intent(in ) :: order
    character (len = *)              , intent(in ) :: gap
    real (dp)                        , intent(in ) :: R0, mu0, delta0, h
    real (dp)                        , intent(out) :: shift
    real (dp), dimension(order,order), intent(out) :: delta

    type (GapMass)                                 :: GapMassive

    delta = 0; shift = 0

    select type (self)
    class is (SingularMassless)

      call self%MatEl%setGammaShift( order, self%run, gap, self%shape, delta0, R0, mu0, &
                                     h, shift, delta )
    class is (SingularMass)

    select type (selector => self%MatEl)
      class is (MatricesElementsMass); GapMassive = GapMass(self%run, gap, selector, R0, mu0)
    end select

      call GapMassive%setGammaShift( order, self%shape, delta0, h, self%R, &
                                     self%muS, shift, delta(1,:) )

      delta(1,:) = delta(1,:) + self%Q * self%deltaM(:order)

      if (order > 1) delta(2 ,2) = delta(1,1)**2/2
      if (order > 2) delta(2:,3) = [ delta(1,1) * delta(1,2), delta(1,1)**3/6 ]

    end select

  end subroutine setGammaShift

!ccccccccccccccc

  subroutine SetRunning(self, muJ, muS, R, mu)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: muJ, muS, R, mu

    real (dp)                               :: wJet, k, kHard, kJet, kSoft, wB, wHard, &
                                               wSoftnf, wSoftnl, wJetNf, wJetNl, wSoft

    self%preFact = 1;  self%muS = muS;  self%R = R

    select type (self)
    class is (SingularScales)

    wJet = 2 * self%alpha%wTilde( self%run, self%cusp, muJ, mu )

    self%w = - 2 * self%alpha%wTilde( self%run, self%cusp, muS, mu) + wJet

    k = self%k - 2 * self%alpha%kTilde( self%run, self%cusp, muS, mu )  &
        - 2 * self%alpha%kTilde( self%run , self%cusp, self%muH, mu )       &
        + 4 * self%alpha%kTilde( self%run , self%cusp, muJ     , mu )       &
        +     self%alpha%wTilde( self%run1, self%soft, muS , mu )       &
        +     self%alpha%wTilde( self%run1, self%jet , muJ , mu )       &
        +     self%alpha%wTilde( self%run1, self%NonCusp, self%muH, mu )

    self%preFact = self%PiResum * (muJ**2/self%Q/muS)**wJet * exp(k) * &
    self%QoMuH**( 2 * self%alpha%wTilde(self%run, self%cusp, self%muH, mu )  )

    call self%MatEl%SetDelta(muS, R, self%alphaList)

    class is (SingularMassScales)

    select type (selector => self%MatEl)
    type is (MatricesElementsMass)
      select type (or => self%MassNS)
      type is (MassiveScales);  selector = or%MatElementScales()
      end select
    end select

    self%m2    = self%MassNS%MassVar('m2')   ; self%muJ = muJ
    self%width = self%MassNS%MassVar('width')
    self%tmin  = self%MassNS%MassVar('tmin') ; self%deltaM = self%MassNS%MassDeltaShift()

    self%MC = MCtop( self%shape, self%MassNS%MassVar('mass'), self%Q)

    if (muJ > self%muM) then  ! SCET scenarios, III and IV

     ! Jet and Hard Running, same for both scenarios III and IV

      if (mu > self%muM) then ! regular mu choice

        wJet = 2 * self%alpha%wTilde( self%run, self%cusp, muJ, mu )

        kJet = 4 * self%alpha%kTilde( self%run    , self%cusp, muJ, mu ) &
                 + self%alpha%wTilde( self%run1, self%jet , muJ, mu )

        wHard = 2 * self%alpha%wTilde( self%run, self%cusp, self%muH, mu )

        kHard = self%alpha%wTilde( self%run1, self%NonCusp, self%muH, mu ) &
          - 2 * self%alpha%kTilde( self%run , self%cusp   , self%muH, mu )

        self%preFact = self%QoMuH**wHard * (muJ**2/self%Q/muS)**wJet

      else  ! anti-natural mu choice

        wJet = 2 * self%alpha%wTilde( self%run, self%cusp, muJ, self%muM )
        wB   = 2 * self%alphaNl%wTilde( self%run, self%cuspNl, self%muM, mu  )

        kJet = 4  *   self%alpha%kTilde( self%run , self%cusp  , muJ, self%muM )   &
                  +   self%alpha%wTilde( self%run1, self%jet   , muJ, self%muM )   &
             +  2 * self%alphaNl%kTilde( self%run , self%cuspNl, self%muM, mu  )   &
                  + self%alphaNl%wTilde( self%run1, self%bjet, self%muM, mu  )

        wHard = 2 * self%alpha%wTilde( self%run, self%cusp, self%muH, self%muM )

        kHard = self%alpha%wTilde( self%run1, self%NonCusp, self%muH, self%muM ) &
          - 2 * self%alpha%kTilde( self%run , self%cusp, self%muH, self%muM )        &
          +   self%alphaNl%wTilde( self%run1, self%Hm, self%muM     , mu  )

        self%preFact = self%QoMuH**wHard * (muJ**2/self%Q/muS)**wJet *         &
                   (self%muM * self%mass/self%Q/muS)**wB * (self%Q/self%mass)**(  2 *   &
                    self%alphaNl%wTilde(self%run, self%cuspNl, self%muM, mu  )  )

        wJet = wJet + wB

      end if

    ! soft sector, different for Scenarios III and IV

    if (self%muM < muS) then ! scenario IV

      if (mu > self%muM) then ! sensible ordering

      wSoft = - 2 * self%alpha%wTilde( self%run , self%cusp, muS, mu )

      kSoft = - 2 * self%alpha%kTilde( self%run , self%cusp, muS, mu ) + &
                    self%alpha%wTilde( self%run1, self%soft, muS, mu )

      else          ! not sensible ordering

      wSoftNl = - 2 * self%alphaNl%wTilde( self%run, self%cuspNl, self%muM, mu  )
      wSoft   = - 2 *   self%alpha%wTilde( self%run, self%cusp  , muS, self%muM ) + wSoftNl

      kSoft = - 2 *   self%alpha%kTilde( self%run ,   self%cusp, muS, self%muM ) + &
                      self%alpha%wTilde( self%run1,   self%soft, muS, self%muM ) - &
                2 * self%alphaNl%kTilde( self%run , self%cuspNl, self%muM, mu  ) + &
                    self%alphaNl%wTilde( self%run1, self%softNl, self%muM, mu  )

        self%preFact = self%preFact * (self%muM/muS)**wSoftnl

      end if

    else ! scenario III

      if (mu > self%muM) then ! sensible ordering

        wSoftnf = - 2 *   self%alpha%wTilde( self%run,   self%cusp, self%muM, mu  )
        wSoft   = - 2 * self%alphaNl%wTilde( self%run, self%cuspNl, muS, self%muM ) + wSoftnf

        kSoft = - 2 * (  self%alphaNl%kTilde( self%run , self%cuspNl, muS, self%muM )    + &
                           self%alpha%kTilde( self%run , self%cusp  , self%muM, mu  )  ) + &
                         self%alphaNl%wTilde( self%run1, self%softNl, muS, self%muM )    + &
                           self%alpha%wTilde( self%run1, self%Soft  , self%muM, mu  )

        self%preFact = self%preFact * (self%muM/muS)**wSoftnf

      else          ! nonsensical ordering

        wSoft = - 2 * self%alphaNl%wTilde( self%run, self%cuspNl, muS, mu )

        kSoft = - 2 * self%alphaNl%kTilde( self%run , self%cuspNl, muS, mu )  + &
                      self%alphaNl%wTilde( self%run1, self%SoftNl, muS, mu )
      end if
    end if

    else ! Scenario II, bHQET

    if (mu < self%muM) then ! sensible ordering

      wJet = 2 * self%alphaNl%wTilde( self%run, self%cuspNl, muJ, mu )

      kJet  = 2 * self%alphaNl%kTilde( self%run , self%cuspNl, muJ, mu  )     &
              +   self%alphaNl%wTilde( self%run1, self%bJet  , muJ, mu  )

      wSoft = - 2 * self%alphaNl%wTilde( self%run , self%cuspNl, muS, mu )

      kSoft = - 2 * self%alphaNl%kTilde( self%run , self%cuspNl, muS, mu ) +  &
                    self%alphaNl%wTilde( self%run1, self%softNl, muS, mu )

      kHard =       self%alpha%wTilde( self%run1, self%NonCusp, self%muH, self%muM )  &
              - 2 * self%alpha%kTilde( self%run , self%cusp, self%muH, self%muM )  &
              +   self%alphaNl%wTilde( self%run1, self%Hm, self%muM     , mu  )

      self%preFact = (muJ/self%QoMpole/muS)**wJet * &
    self%QoMuH**(      2 *   self%alpha%wTilde(self%run, self%cusp  , self%muH, self%muM )  )  &
    * self%QoMpole**(  2 * self%alphaNl%wTilde(self%run, self%cuspNl, self%muM, mu  )  )

    else  ! nonsensical ordering

      wJetNf = 2 *   self%alpha%wTilde( self%run, self%cusp, self%muM, mu )
      wJetNl = 2 * self%alphaNl%wTilde( self%run, self%cuspNl, muJ, self%muM )

      wJet =  wJetNf + wJetNl

      kJet  = 2 * self%alphaNl%kTilde( self%run , self%cuspNl, muJ, self%muM )      &
              +   self%alphaNl%wTilde( self%run1, self%bJet, muJ, self%muM )        &
              + 4 * self%alpha%kTilde( self%run , self%cusp  , self%muM, mu  )      &
              +     self%alpha%wTilde( self%run1, self%jet   , self%muM, mu  )

      wSoftnf = - 2 *   self%alpha%wTilde( self%run,   self%cusp, self%muM, mu  )
      wSoft   = - 2 * self%alphaNl%wTilde( self%run, self%cuspNl, muS, self%muM ) + wSoftnf

      kSoft = - 2 * self%alphaNl%kTilde( self%run , self%cuspNl, muS, self%muM )  + &
                    self%alphaNl%wTilde( self%run1, self%softNl, muS, self%muM )  - &
                2 *   self%alpha%kTilde( self%run ,   self%cusp, self%muM, mu  )  + &
                      self%alpha%wTilde( self%run1,   self%soft, self%muM, mu  )

      kHard =       self%alpha%wTilde( self%run1, self%NonCusp, self%muH, mu )  &
              - 2 * self%alpha%kTilde( self%run , self%cusp, self%muH, mu )

      self%preFact = (self%muM/muS)**wSoftnf * (muJ/self%QoMpole/muS)**wJetNl * &
                     (self%muM**2/self%Q/muS)**wJetNf * self%QoMuH**(  2      * &
                      self%alpha%wTilde(self%run, self%cusp, self%muH, mu )  )
    end if

    end if

    self%w = wSoft + wJet

    self%preFact = self%preFact * self%PiResum * exp(self%k + kHard + kJet + kSoft)

    end select

  end subroutine SetRunning

!ccccccccccccccc

  type (MassiveNS) function NSMass(self)
    class (SingularMass), intent(in) :: self

    select type (selector => self%MassNS)
    type is (MassiveNS);  NSMass = selector
    end select

  end function NSMass

!ccccccccccccccc

  type (MassiveScales) function NSScales(self)
    class (SingularMassScales), intent(in) :: self

    select type (selector => self%MassNS)
    type is (MassiveScales);  NSScales = selector
    end select

  end function NSScales

!ccccccccccccccc

  subroutine SetMat(self, muJ, muS, xiJ, xiB)
    class (SingularMassless), intent(inout) :: self
    real (dp)               , intent(in)    :: muJ, muS
    real (dp), optional     , intent(in)    :: xiJ, xiB
    real (dp), dimension(0:4,3)             :: MatJet, MatSoft
    real (dp)                               :: lJet, XJ, xB
    integer                                 :: i, j

    self%MatExp = 0; xJ = 0; xB = 0
    if ( present(xiJ) ) XJ = xiJ; if ( present(xiB) ) XB = xiB

    self%HardExp = self%ExpHard

    select type (self)
    class is (SingularScales)

    self%alphaList = powList( self%alpha%alphaQCD(muS)/Pi, 3 )

    self%NGLog2loop = self%alphaList(2) * self%NGL2loop

    MatSoft = self%MatSoft; MatJet = self%MatJet

    call self%MatEl%AddLogs(MatJet, log(muJ**2/self%Q/muS))
    call AddAlpha(MatSoft, self%alphaList)
    call AddAlpha(MatJet, powList( self%alpha%alphaQCD(muJ)/Pi, 3 ))

    self%MatExp = self%MatHard + MatSoft + MatJet

    class is (SingularMassScales)

    if ( self%hard(:6) == 'expand' ) then
      if (muJ < self%muM) then
        self%MatExp = self%MatExp + self%MatMassHard
      end if
    else
      self%HardExp(1) = self%HardExp(1) + (1 - xJ) * self%MassNS%MassVar('Hcorr') * &
                                                     self%MassNS%MassVar('alphaH')
      if (muJ < self%muM) then
        self%HardMassExp    = self%MatMassHard(:3,1)
        self%HardMassExp(1) = self%HardMassExp(1) + xB * self%MassNS%MassVar('Hcorr') &
                                           * (1 - xJ) * self%MassNS%MassVar('alphaM')
      end if
    end if

    if ( muS >= self%muM ) then
      self%alphaList = powList( self%alpha%alphaQCD(muS)/Pi, 3 )
      MatSoft = self%MatSoft; self%NGL2loop = self%alphaList(2) * self%NGLog2loop
    else
      self%alphaList = powList( self%alphaNl%alphaQCD(muS)/Pi, 3 )
      MatSoft = self%MatSoftNl; self%NGL2loop = self%alphaList(2) * self%NGLog2loopNl
    end if

    if ( muJ >= self%muM ) then
      self%alphaJ = self%alpha%alphaQCD(muJ)/Pi  ;  MatJet = self%MatJet
      lJet = log(muJ**2/self%Q/muS)
    else
      self%alphaJ = self%alphaNl%alphaQCD(muJ)/Pi;  MatJet = self%MatBJet
      lJet = log(muJ * self%MassNS%MassVar('mass')/self%Q/muS)
    end if

    call self%MatEl%AddLogs( MatJet, lJet )

    call AddAlpha(  MatSoft, self%alphaList  )
    call AddAlpha(  MatJet , powList( self%alphaJ, 3 )  )

    self%MatExp = self%MatExp + self%MatHard + MatSoft + MatJet

    self%A0 = self%MassNS%MassVar('A0')

    if (muJ >= self%muM) then
      self%MatExp = self%MatExp + 2 * self%MassNS%JetMassExp(self%hard, xJ)/self%A0
      if (muS < self%muM) self%MatExp = self%MatExp + self%SoftMatchingExp()
    else
      self%MatExp = self%MatExp + 2 * self%MassNS%JetMassExp(self%hard, xJ, xB)/self%A0
    end if

    end select

    self%MatExpanded    = self%A0 * expandExp(self%MatExp  ) ; self%MatAdded = 0
    self%MatExpandedHJM = self%A0 * expandExp(self%MatExp/2)

    do i = 0, 3
      do j = 0, 2 * i
        self%MatAdded(j,i) = Sum(  self%MatExpanded( j, ceiling(j/2.):i )  )
      end do
    end do

    self%muSEuler = muS * ExpEuler;   self%compfact = 1/self%muSEuler

    self%MatFO    = matmul( PosToMomMatrix(), self%MatAdded       )
    self%MatHJMFO = matmul( PosToMomMatrix(), self%MatExpandedHJM )

  end subroutine SetMat

!ccccccccccccccc

  subroutine SetEverything(self, mu, Q, muH, muJ, muS, R, Rmass, muM, width)
    class (SingularMassScales), intent(inout) :: self
    real (dp)                 , intent(in)    :: mu, muJ, muS, R, Rmass, width, muM, Q, muH

    call self%MassNS%SetEverything(mu, Q, muH, muJ, muS, R, Rmass, muM, width)

  end

!ccccccccccccccc

  integer function runOrd(self)
    class (SingularMassless) , intent(in) :: self
    runOrd = self%run
  end function runOrd

!ccccccccccccccc

  real (dp) function SingleSingMod(self, Mod, setup, gap, space, cum, order, R0, &
                                    mu0, delta0, h, tau, tau2)
    class (SingularMassless) , intent(in)            :: self
    type (Model)       , intent(in)                  :: Mod
    character (len = *), intent(in)                  :: setup, space, gap, cum
    real (dp)          , intent(in)                  :: R0, mu0, delta0, tau, h
    real (dp), optional, intent(in)                  :: tau2
    integer            , intent(in)                  :: order
    real (dp), dimension(order,order)                :: delta
    real (dp), dimension(0:2 * order)                :: kerSingle
    real (dp), dimension( order, 0:2 * (order - 1) ) :: deltaAdd
    real (dp), dimension(0:order,0:2 * order)        :: kerMatrix
    real (dp), dimension(2 * order)                  :: logList
    logical                                          :: dobsing
    integer                                          :: i, j, k, neval, ier
    type (Kernels)                                   :: ker
    real (dp), dimension(1)                          :: res
    real (dp)                                        :: w, p, shift, p2, &
    pshift, pshift2, abserr, result

    SingleSingMod = 0

    if (present (tau2) .or. order < 0 ) then
      if (tau >= tau2)  return
    end if

    if ( setup(:5) == 'Model' ) then
      if ( present(tau2) ) then
        res = self%SingleSingList([Mod], gap, space, cum, order, R0, mu0, &
                                      delta0, h, tau, tau2)
      else
        res = self%SingleSingList([Mod], gap, space, cum, order, R0, mu0, &
                                      delta0, h, tau)
      end if
      SingleSingMod = res(1); return
    end if

    deltaAdd = 0; kerSingle = 0; kermatrix = 0; w = self%w + cumConst(cum)

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    pshift = self%Q * tau/self%ESFac - shift  ;  pshift2 = pshift
    p = pshift - self%Q * self%tmin/self%ESFac;  p2 = p

    if ( present(tau2) ) then
      pshift2 = self%Q * tau2/self%ESFac - shift
      p2 = pshift2 - self%Q * self%tmin/self%ESFac
    end if

    if ( p <= 0 .and. present(tau2) ) p = p2

! Result = 0 if p is negative

    if ( p2 <= 0 .and. setup(:15) /= 'NoModelUnstable') return

    dobsing = .not. present(tau2) .or. ( present(tau2) .and. abs(p2 - p) <= d1mach(1) )

! FO will only work if all scales are set equal and there is no gap

    if ( setup(:2) == 'FO' ) then
      if ( cum(:3) == 'cum' ) then

        logList = powList(log(p/self%muS), 2*order)

        if ( present(tau2) .and. p2 > p ) logList = powList(log(p2/self%muS), 2*order) - logList

        logList = [ ( logList(i)/i, i = 1, 2 * order ) ]

        if ( dobsing ) SingleSingMod = self%MatFO(0, order)

        SingleSingMod = SingleSingMod + dot_product( self%MatFO(1:2*order, order), logList )

      else

        logList(1) = 1/tau; logList(2:) = powList(log(p/self%muS), 2*order - 1)/tau

        if ( present(tau2) .and. p2 > p ) then
          logList(1) = 1/tau2 - logList(1)
          logList(2:) = powList(log(p2/self%muS), 2*order - 1)/tau2 - logList(2:)
        end if

        SingleSingMod = dot_product( self%MatFO(1:2*order, order), logList )

      end if
      return
    end if

    ker = Kernels(n = 2 * order, w = w); kerSingle = ker%DerList()
    kerMatrix(0,:) = KernelMatrixList( kerSingle, self%MatAdded(:2*order,order) )

    if ( sum( abs(delta) ) > 0 ) then

      do j = 1, order
        do k = 0, 2 * (order - j)

          deltaAdd(j,k) = sum(   self%MatAdded( k, order - j: ceiling(k/2.): -1 ) * &
                                 delta( j,j:order - ceiling(k/2.) )   )
        end do
      end do

      do i = 1, order

        kerSingle( :2 * (order - i) ) = DerAdd1(-1, kerSingle( :2 * (order - i) ), w + i - 1  )
        kerMatrix( i,:2 * (order - i) ) = KernelMatrixList( kerSingle( :2 * (order - i) ),&
         deltaAdd(i,:2 * (order - i)) )

      end do
    end if

    SingleSingMod = NoMod(p)
    if ( present(tau2) .and. p2 > p ) SingleSingMod = NoMod(p2) - SingleSingMod

    select type (self)
    class is (SingularMass)

      if ( setup(8:15) == 'Unstable' ) then

        if ( self%shape(:6) == 'thrust' ) SingleSingMod = NoMod(p) * self%MC%Delta()

        if ( present(tau2) .and. self%shape(:6) == 'thrust' ) SingleSingMod = &
        self%MC%Delta() * NoMod(p2) - SingleSingMod

        call qags( UnstableInt, pshift - self%MC%maxP(), pshift2, prec, &
        prec, result, abserr, neval, ier )

        SingleSingMod = result + SingleSingMod

      end if

    end select

    SingleSingMod = SingleSingMod * self%Prefact * self%compFact**cumConst(cum) * &
           (self%Q/self%ESFac)**( 1 + cumConst(cum) )  * &
            Sum( self%HardMassExp(:order) ) * Sum( self%HardExp(:order) )

    return ! the rest will not be evaluated

    select type (self)
    class is (SingularMass)

      if ( self%muJ < self%muM ) return

      if ( .not. present(tau2) ) then

        SingleSingMod = SingleSingMod + self%NonDistMod(Mod, setup, gap, space, &
                                           cum, order, R0, mu0, delta0, h, tau)
       else

        SingleSingMod = SingleSingMod + self%NonDistMod(Mod, setup, gap, space, &
                                     cum, order, R0, mu0, delta0, h, tau, tau2)
      end if
    end select

  contains

!ccccccccccccccc

  real (dp) function UnstableInt(x)
    real (dp), intent(in) :: x

    select type (self)
    class is (SingularMass)

      if (.not. present(tau2) ) then
        UnstableInt = NoMod(x) * self%MC%Qdist(pshift - x)
      else
        UnstableInt = NoMod(x) * self%MC%Qdist( pshift - x, pshift2 - x )
      end if

    end select

  end function UnstableInt

!ccccccccccccccc

  recursive real (dp) function NoMod(x, x2) result(res)
    real (dp)          , intent(in) :: x
    real (dp), optional, intent(in) :: x2

    res = 0

    if ( present(x2) ) then
      if (x2 > x) res = NoMod(x2) - NoMod(x); return
    else
      if (x < 0) return
    end if

    res = LogMatrixSum( kerMatrix(0,:), w, self%muS, x )

    if (  sum( abs(delta) ) > 0  ) then
      do j = 1, order

        res = res + (-self%compFact)**j * LogMatrixSum(  kerMatrix( j,:2 * (order - j) ), &
                         w + j, self%muS, x  )
      end do
    end if

  end function NoMod

  end function SingleSingMod

!ccccccccccccccc

  function SingleSingList(self, ModList, gap, space, cum, order, R0, mu0, delta0, h, &
                          tau, tau2) result(resList)
    class (SingularMassless)  , intent(in)           :: self
    type (Model), dimension(:), intent(in)           :: ModList
    character (len = *)       , intent(in)           :: space, gap, cum
    real (dp)                 , intent(in)           :: R0, mu0, delta0, tau, h
    real (dp)       , optional, intent(in)           :: tau2
    integer                   , intent(in)           :: order
    real (dp), dimension(order,order)                :: delta
    real (dp), dimension(0:2 * order)                :: kerSingle, kerSingleOPE
    real (dp), dimension( order, 0:2 * (order - 1) ) :: deltaAdd
    real (dp), dimension(0:order,0:2 * order)        :: kerMatrix
    complex (dp), dimension(3)                       :: sumtot2
    logical                                          :: dobsing
    integer                                          :: neval, ier, i, j, k, l, fac
    type (Kernels)                                   :: ker
    real (dp)                                        :: w, abserr, p, shift, p2, &
                                                        pmid, newTerm
    real (dp), dimension( size(ModList) )            :: resList, Omega
    character (len = 6)                              :: space2

    deltaAdd = 0; kerSingle = 0; kermatrix = 0; resList = 0; space2 = space
    w = self%w + cumConst(cum);  shift = 0; l = 0

    if (present (tau2) .or. order < 0 ) then
      if (tau >= tau2) return
    end if

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    p = (self%Q * tau - self%Q * self%tmin)/self%ESFac - shift; p2 = p
    if ( present(tau2) ) p2 = (self%Q * tau2 - self%Q * self%tmin)/self%ESFac - shift

    if ( p <= 0 .and. present(tau2) ) p = p2

    pmid = (p + p2)/2

! Result = 0 if p is negative

    if ( p2 <= 0 ) return

    dobsing = .not. present(tau2) .or. ( present(tau2) .and. abs(p2 - p) <= d1mach(1) )

    if ( space(:3) == 'mom' .and. cum(:6) == 'cumMod' ) w = self%w

    if ( space(:3) /= 'pos' ) then

      ker = Kernels(n = 2 * order, w = w); kerSingle = ker%DerList(); kerSingleOPE = kerSingle
      kerMatrix(0,:) = KernelMatrixList( kerSingle, self%MatAdded(:2*order,order) )

    end if

    if ( sum( abs(delta) ) > 0 .and. space(:6) /= 'posExp' ) then

      do j = 1, order
        do k = 0, 2 * (order - j)

          deltaAdd(j,k) = sum(   self%MatAdded( k, order - j: ceiling(k/2.): -1 ) * &
                                 delta( j,j:order - ceiling(k/2.) )   )
        end do
      end do

      if ( space(:3) /= 'pos' ) then

        do i = 1, order

          kerSingle( :2 * (order - i) ) = DerAdd1(-1, kerSingle( :2 * (order - i) ), w + i - 1  )
          kerMatrix( i,:2 * (order - i) ) = KernelMatrixList(  kerSingle( :2 * (order - i) ),&
           deltaAdd( i,:2 * (order - i) )  )

        end do
      end if
    end if

    if ( space(:6) == 'posExp' ) then
      sumtot2 = 0; sumtot2(:order) = self%MatExp(0,:order)
    end if

    if ( space(:3) == 'OPE' ) then

      if (dobsing) then
        resList = NoMod(p)
      else
        resList = NoMod(p2,p)
      end if

      resList = resList * [  ( ModList(i)%MomentModel(0), i = 1, size(modList) )  ]
      kerMatrix = 0; fac = 1

      do l = 1, 30

        kerSingle = DerAdd1(-1, kerSingleOPE, w); kerSingleOPE = kerSingle; w = w + 1
        kerMatrix(0,:) = KernelMatrixList( kerSingle, self%MatAdded(:2*order,order) )

        if (  sum( abs(delta) ) > 0  ) then
          do i = 1, order

            kerSingle( :2 * (order - i) ) = DerAdd1(-1, kerSingle( :2 * (order - i) ), w + i - 1  )
            kerMatrix( i,:2 * (order - i) ) = KernelMatrixList(  kerSingle( :2 * (order - i) ),&
            deltaAdd( i,:2 * (order - i) )  )

          end do
        end if

        Omega = [  ( ModList(i)%MomentModel(l), i = 1, size(modList) )  ]; fac = fac * l

        if (dobsing) then
          newTerm = (-self%compFact)**l * NoMod(p)/fac
        else
          newTerm = (-self%compFact)**l * ( NoMod(p2,p) )/fac
        end if

        if ( abs( newTerm * maxval(Omega) ) > abs( maxval(resList) ) ) then
          space2 = 'posExp'; resList = 0; w = self%w + cumConst(cum)
          sumtot2 = 0; sumtot2(:order) = self%MatExp(0,:order);  exit
        end if

        resList = resList + newTerm * Omega; kerMatrix = 0
        if ( abs( maxval(newTerm * Omega/resList) ) < prec ) exit

      end do

    else if ( space(:6) == 'Taylor' ) then

      fac = 1

      do l = 0, 80

        Omega = [  ( ModList(i)%Taylor(l), i = 1, size(modList) )  ];  fac = fac * max(l,1)
        kerSingle = DerSubtract1(-1, kerSingleOPE, w); kerSingleOPE = kerSingle; w = w - 1
        if (  maxval( abs(Omega) ) <= d1mach(1)  ) cycle
        kerMatrix(0,:) = KernelMatrixList( kerSingle, self%MatAdded(:2*order,order) )

        if (  sum( abs(delta) ) > 0  ) then
          do i = 1, order

            kerSingle( :2 * (order - i) ) = DerAdd1(-1, kerSingle( :2 * (order - i) ), w + i - 1  )
            kerMatrix( i,:2 * (order - i) ) = KernelMatrixList(  kerSingle( :2 * (order - i) ),&
            deltaAdd( i,:2 * (order - i) )  )

          end do
        end if

        if (dobsing) then
          newTerm = self%compFact**(-l - 1) * NoMod(p) * fac
        else
          newTerm = self%compFact**(-l - 1) * ( NoMod(p2,p) ) * fac
        end if

        resList = resList + newTerm * Omega; kerMatrix = 0
        if ( abs( newTerm * maxval(Omega) ) < prec ) exit

      end do

    end if

    do i = 1, size(modList)

      if ( space2(:3) == 'mom' .and. cum(:6) == 'cumMod' ) then

        call qags( ModCumInt, 0._dp, p2, prec, prec, resList(i), abserr, neval, ier )

        resList(i) = resList(i)/self%muSEuler

      else if ( space2(:3) == 'mom' ) then

        call qags( ModInt, 0._dp, p2, prec, prec, resList(i), abserr, neval, ier )

      else if ( space2(:6) == 'posExp' ) then

        call qagi( FourierExp, 0._dp, 1, prec, prec, resList(i), abserr, neval, ier )

        resList(i) = self%A0 * resList(i)/Pi/pmid

      else if ( space2(:3) == 'pos' ) then

        call qagi( Fourier, 0._dp, 1, prec, prec, resList(i), abserr, neval, ier )

        resList(i) = resList(i)/Pi/pmid

      end if
    end do

    resList = resList * self%Prefact * self%compFact**cumConst(cum) * &
                 (self%Q/self%ESFac)**( 1 + cumConst(cum) )         * &
                Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) )

    return ! the rest will not be evaluated

    select type (self)
    class is (SingularMass)

      if ( self%muJ < self%muM ) return

      if ( .not. present(tau2) ) then

        resList = resList + self%NonDistList(ModList, gap, space, cum, order, R0, &
                                                 mu0, delta0, h, tau)
      else

        resList = resList + self%NonDistList(ModList, gap, space, cum, order, R0, &
                                                 mu0, delta0, h, tau, tau2)
      end if
    end select

  contains

!ccccccccccccccc

  recursive real (dp) function NoMod(x, x2) result(res)
    real (dp)          , intent(in) :: x
    real (dp), optional, intent(in) :: x2
    integer                         :: j

    res = 0

    if ( present(x2) ) then
      if (x2 > x) res = NoMod(x2) - NoMod(x); return
    end if

    res = LogMatrixSum( kerMatrix(0,:), w, self%muS, x )

    if (  sum( abs(delta) ) > 0  ) then
      do j = 1, order

        res = res + (-self%compFact)**j * LogMatrixSum(  kerMatrix( j,:2 * (order - j) ), &
                         w + j, self%muS, x  )
      end do
    end if

  end function NoMod

!ccccccccccccccc

  real (dp) function ModInt(x)
    real (dp), intent(in) :: x

    if (dobsing) then
      ModInt = NoMod(x) * ModList(i)%ShapeFun(0, p2 - x)
    else
      ModInt = NoMod(x) * ModList(i)%ShapeFun(0, p - x, p2 - x)
    end if

  end function ModInt

!ccccccccccccccc

  real (dp) function ModCumInt(x)
    real (dp), intent(in) :: x

    if ( dobsing ) then
      ModCumInt = NoMod(x) * ModList(i)%ShapeFun(-1, p2 - x)
    else
      ModCumInt = NoMod(x) * ModList(i)%ShapeFun(-1, p - x, p2 - x)
    end if

  end function ModCumInt

!ccccccccccccccc

  real (dp) function FourierExp(x)
    real (dp),      intent(in) :: x
    complex (dp)               :: y, z, t, tmuS, CI, C11, core, logTab(order + 1)
    complex (dp), dimension(3) :: sumtot
    integer                    :: j

    CI = (0, 1);  C11 = (1, 1);  z = C11 * x - CI/10; y = z/pmid
    t = CI * y; tmuS = t * self%muSEuler; logTab = powList( log(tmuS), order + 1 )

    if (order == 0) then
      core = 1
    else

      sumtot = sumtot2

      if (  abs( sum(delta) ) > 0  ) sumtot(:order) = sumtot(:order) - t * delta(1,:)

      do j = 1, order + 1
        sumtot = sumtot + self%MatExp(j,:) * logTab(j)
      enddo

      core = expandExp( sumtot(:order) )

    end if

    core = core * C11 * tmuS**w * ModList(i)%SoftFourier(y)

    if ( dobsing ) then
      FourierExp = real(  core * exp(CI * z)  )
    else
      FourierExp = real(  core * ( exp( p2 * t) - exp(p * t) )  )
    end if

  end function FourierExp

!ccccccccccccccc

  real (dp) function Fourier(x)
    real (dp), intent(in) :: x
    complex (dp)          :: y, z, t, tmuS, lgc, CI, C11, core, logTab(0:2*order)
    integer               :: j

    CI = (0, 1);  C11 = (1, 1);  z = C11 * x - CI/10; y = z/pmid
    t = CI * y ;  tmuS = t * self%muSEuler; lgc = log(tmuS)

    logTab(0) = 1; logTab(1:) = powList(lgc, 2*order)

    core = sum( self%MatAdded(:2 * order, order) * logTab )

    if (  abs( sum(delta) ) > 0  ) then

      do j = 1, order
        core = core + (-t)**j * (   sum(  deltaAdd( j, :2 * (order - j) ) &
                                   * logTab( :2 * (order - j) )  )   )
      end do
    end if

    core = core * C11 * tmuS**w * ModList(i)%SoftFourier(y)

    if ( dobsing ) then
      Fourier = real(  core * exp(CI * z)  )
    else
      Fourier = real(  core * ( exp(p2 * t) - exp(p * t) )  )
    end if
  end function Fourier

  end function SingleSingList

!ccccccccccccccc

  real (dp) function SingleSingWidthMod(self, Mod, setup, gap, space, cum, order, &
    R0, mu0, delta0, h, tau, tau2)
    class (SingularMass), intent(in)                 :: self
    type (Model)        , intent(in)                 :: Mod
    character (len = *) , intent(in)                 :: setup, gap, cum, space
    real (dp)           , intent(in)                 :: R0, mu0, delta0, tau, h
    real (dp)           , intent(in), optional       :: tau2
    integer             , intent(in)                 :: order
    type (Kernels)                                   :: ker
    real (dp), dimension(order,order)                :: delta
    real (dp), dimension( order, 0:2 * (order - 1) ) :: deltaAdd
    real (dp), dimension( 0:order, 3, 0:2*order )    :: derInvList
    integer                                          :: i, j, neval, ier
    real (dp), dimension(1)                          :: res
    real (dp)                                        :: w, p, shift, p2, pshift, &
    pshift2, result, abserr

    deltaAdd = 0; SingleSingWidthMod = 0; derInvList = 0

    if (present (tau2) .or. order < 0 ) then
      if (tau >= tau2) return
    end if

    if ( self%width <= d1mach(1) ) then
      if ( .not. present(tau2) ) SingleSingWidthMod = self%SingleSingMod(Mod, setup, gap, &
                                 'posExp', cum, order, R0, mu0, delta0, h, tau)
      if (       present(tau2) ) SingleSingWidthMod = self%SingleSingMod(Mod, setup, gap, &
                                 'posExp', cum, order, R0, mu0, delta0, h, tau, tau2)
      return
    end if

    if ( setup(:5) == 'Model' ) then
      if ( present(tau2) ) then
        res = self%SingleSingWidthList([Mod], setup, gap, space, cum, order, R0, mu0, &
                                        delta0, h, tau, tau2)
      else
        res = self%SingleSingWidthList([Mod], setup, gap, space, cum, order, R0, mu0, &
                                     delta0, h, tau)
      end if
      SingleSingWidthMod = res(1); return
    end if

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    pshift = self%Q * tau/self%ESFac - shift; pshift2 = pshift
    p = pshift - self%Q * self%tmin/self%ESFac; p2 = p

    if ( present(tau2) ) then
      pshift2 = self%Q * tau2 - shift
      p2 = pshift2 - self%Q * self%tmin/self%ESFac
    end if

    w = self%w + cumConst(cum)

    if ( sum( abs(delta) ) > 0 ) then
      do i = 1, order
        do j = 0, 2 * (order - i)

          deltaAdd(i,j) = sum(   self%MatAdded( j, order - i: ceiling(j/2.): -1 ) * &
                                 delta( i,i:order - ceiling(j/2.) )   )
        end do
      end do
    end if

    ker = Kernels(n = 2 * order, width = self%width, w = w)
    derInvList(0,1,:) = ker%DerList()

    if (  sum( abs(delta) ) > 0  ) then
      do i = 1, order
        derInvList(i,1,:) = DerAdd1(1, derInvList( i - 1,1,:2 * (order - i) ), w + i + 1  )
      end do
    end if

    SingleSingWidthMod = NoMod(p)
    if ( present(tau2) ) SingleSingWidthMod = NoMod(p2) - SingleSingWidthMod

    if ( setup(8:15) == 'Unstable' ) then

      if ( self%shape(:6) == 'thrust' ) SingleSingWidthMod = NoMod(p) * self%MC%Delta()

      if ( present(tau2) .and. self%shape(:6) == 'thrust' ) SingleSingWidthMod = &
      self%MC%Delta() * NoMod(p2) - SingleSingWidthMod

      call qags( UnstableInt, pshift - self%MC%maxP(), pshift2, prec, &
      prec, result, abserr, neval, ier )

      SingleSingWidthMod = result + SingleSingWidthMod

    end if

    SingleSingWidthMod = SingleSingWidthMod * self%Prefact * self%compFact**cumConst(cum) * &
                        Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) ) * &
                        (self%Q/self%ESFac)**( 1 + cumConst(cum) )

    return ! the rest will not be evaluated

    select type (self)
    class is (SingularMass)
      if ( self%muJ < self%muM ) return

      if ( .not. present(tau2) ) then

        SingleSingWidthMod = SingleSingWidthMod + self%NonDist(Mod, setup, space, &
                                            gap, cum, order, R0, mu0, delta0, h, tau)
      else

        SingleSingWidthMod = SingleSingWidthMod + self%NonDist(Mod, setup, space, &
                                     gap, cum, order, R0, mu0, delta0, h, tau, tau2)
      end if
    end select

  contains

!ccccccccccccccc

  real (dp) function UnstableInt(x)
    real (dp), intent(in) :: x

    if (.not. present(tau2) ) then
      UnstableInt = NoMod(x) * self%MC%QDist( pshift - x )
    else
      UnstableInt = NoMod(x) * self%MC%QDist( pshift - x, pshift2 - x )
    end if

  end function UnstableInt

!ccccccccccccccc

  recursive real (dp) function NoMod(x, x2) result(res)
    real (dp)            , intent(in) :: x
    real (dp)  , optional, intent(in) :: x2
    real (dp), dimension(0:2 * order) :: kerMatrix
    real (dp), dimension(0:2 * order) :: kerSingle
    real (dp)                         :: y

    res = 0

    if ( present(x2) ) then
      if (x2 > x) res = NoMod(x2) - NoMod(x); return
    end if

    kermatrix = 0; res = 0; kerSingle = 0; y = 0.5_dp + Atan(x/self%width)/Pi

    derInvList(0,2:,:) = KernelWidth(2*order, w, y)
    kerSingle          = SumKernelWidth(DerInvList(0,:,:), y)
    kerMatrix          = KernelMatrixList( kerSingle, self%MatAdded(:2*order,order) )
    res              = LogMatrixSum( kerMatrix, w, self%muS, x, self%width )

    if (  sum( abs(delta) ) > 0  ) then
      do i = 1, order

        derInvList( i,2:,:2 * (order - i) ) = KernelWidth(2 * (order - i), w + i, y)
        kerSingle( :2 * (order - i) ) = SumKernelWidth(DerInvList(i,:,:2 * (order - i)), y)
        kerMatrix( :2 * (order - i) ) = KernelMatrixList( kerSingle( :2 * (order - i) ),&
        deltaAdd(i,:2 * (order - i)) )

        res = res + (-self%compFact)**i * LogMatrixSum(  kerMatrix( :2 * (order - i) ), &
        w + i, self%muS, x, self%width  )

      end do
    end if

  end function NoMod

  end function SingleSingWidthMod

!ccccccccccccccc

  function SingleSingWidthList(self, ModList, setup, gap, space, cum, order, R0, &
  mu0, delta0, h, tau, tau2) result(resList)
    class (SingularMass), intent(in)                 :: self
    type (Model)        , intent(in), dimension(:)   :: ModList
    character (len = *) , intent(in)                 :: space, gap, cum, setup
    real (dp)           , intent(in)                 :: R0, mu0, delta0, tau, h
    real (dp)           , intent(in), optional       :: tau2
    integer             , intent(in)                 :: order
    type (Kernels)                                   :: ker, kerOr
    real (dp), dimension(order,order)                :: delta
    real (dp), dimension( order, 0:2 * (order - 1) ) :: deltaAdd
    integer                                          :: neval, ier, i, j, l, fac
    real (dp), dimension( size(modList) )            :: resList, Omega
    real (dp), dimension( 0:order, 3, 0:2*order )    :: derInvList, derInvListOr
    real (dp)                                        :: w, abserr, p, shift, p2, &
    pshift, pshift2, newTerm, result

    delta = 0; deltaAdd = 0; shift = 0; derInvList = 0; fac = 1

    if (present (tau2) .or. order < 0 ) then
      if (tau >= tau2) return
    end if

    if ( self%width <= d1mach(1) ) then
      if ( .not. present(tau2) ) resList = self%SingleSingList(ModList, gap, &
                                 'posExp', cum, order, R0, mu0, delta0, h, tau)
      if (       present(tau2) ) resList = self%SingleSingList(ModList, gap, &
                                 'posExp', cum, order, R0, mu0, delta0, h, tau, tau2)
      return
    end if

    call self%setGammaShift( order, gap, delta0, R0, mu0, h, shift, delta )

    pshift = self%Q * tau/self%ESFac - shift; pshift2 = pshift
    p = pshift - self%Q * self%tmin/self%ESFac; p2 = p

    if ( present(tau2) ) then
      pshift2 = self%Q * tau2/self%ESFac - shift
      p2 = pshift2 - self%Q * self%tmin/self%ESFac
    end if

    w = self%w + cumConst(cum); if ( cum(:6) == 'cumMod' ) w = self%w

    if ( sum( abs(delta) ) > 0 ) then
      do i = 1, order
        do j = 0, 2 * (order - i)

          deltaAdd(i,j) = sum(   self%MatAdded( j, order - i: ceiling(j/2.): -1 ) * &
                                 delta( i,i:order - ceiling(j/2.) )   )
        end do
      end do
    end if

    ker = Kernels(n = 2 * order, width = self%width, w = w)
    derInvList(0,1,:) = ker%DerList()

    if (  sum( abs(delta) ) > 0  ) then
      do i = 1, order
        derInvList(i,1,:) = DerAdd1(1, derInvList( i - 1,1,:2 * (order - i) ), w + i + 1  )
      end do
    end if

    kerOr = ker; derInvListOr = derInvList

    if ( setup(6:13) == 'Unstable' ) then

      do i = 1, size(modList)

        call qagi( UnstableInt, 0._dp, 1, prec, prec, resList(i), abserr, neval, ier )

        if ( self%shape(:6) == 'thrust' ) then
          call qagi( ModInt, 0._dp, 1, prec, prec, result, abserr, neval, ier )
          resList(i) = resList(i) + result * self%MC%Delta()
        end if

      end do

    else

      if ( space(:3) == 'OPE' ) then

        resList = NoMod(p); if ( present(tau2) ) resList = NoMod(p2) - resList
        resList = resList * [  ( ModList(i)%MomentModel(0), i = 1, size(modList) )  ]

        do l = 1, 30

          Omega = [  ( ModList(i)%MomentModel(l), i = 1, size(modList) )  ]
          w = w + 1; ker = Kernels(n = 2 * order, width = self%width, w = w)
          derInvList(0,1,:) = ker%DerList(); fac = fac * l

          if (  sum( abs(delta) ) > 0  ) then
            do i = 1, order
              derInvList(i,1,:) = DerAdd1(1, derInvList( i - 1,1,:2 * (order - i) ), w + i + 1  )
            end do
          end if

          newTerm = NoMod(p); if ( present(tau2) ) newTerm = NoMod(p2) - newTerm
          newTerm = (-self%compFact)**l * newTerm/fac

          if ( abs( newTerm * maxval(Omega) ) > abs( maxval(resList) ) ) then
            resList = 0; w = self%w + cumConst(cum); ker = kerOr
            if ( cum(:6) == 'cumMod' ) w = self%w
            derInvList = derInvListOr;  go to 10
          end if

          resList = resList + newTerm * Omega
          if ( abs( maxval(newTerm * Omega/resList) ) < prec ) go to 20

        end do

        go to 20

      end if

      10 continue

      do i = 1, size(modList)

        if ( cum(:6) == 'cumMod' ) then
          call qagi( ModCumInt, p2, -1, prec, prec, resList(i), abserr, neval, ier )
          resList(i) = resList(i)/self%muSEuler
        else
          call qagi( ModInt, 0._dp, 1, prec, prec, resList(i), abserr, neval, ier )
        end if

      end do

      20 continue

    end if

    resList = resList * self%Prefact * self%compFact**cumConst(cum) * &
    Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) ) * &
    (self%Q/self%ESFac)**( 1 + cumConst(cum) )

    return ! the rest will not be evaluated

    select type (self)
    class is (SingularMass)

      if (self%muJ < self%muM) return

      if ( .not. present(tau2) ) then
        resList = resList + self%NonDistList(ModList, gap, 'space', cum, order, &
        R0, mu0, delta0, h, tau)
      else
        resList = resList + self%NonDistList(ModList, gap, 'space', cum, order, &
        R0, mu0, delta0, h, tau, tau2)
      end if
    end select

  contains

!ccccccccccccccc

  real (dp) function UnstableInt(x, x2)
    real (dp)          , intent(in) :: x
    real (dp), optional, intent(in) :: x2

    if ( .not. present(tau2) ) then
      UnstableInt = NoMod(pshift - x) * ModList(i)%ModelUnstable(self%MC, 0, x)
    else
      UnstableInt = NoMod(pshift - x, pshift2 - x) * ModList(i)%ModelUnstable(self%MC, 0, x) ! the smart way does not work...
    end if

  end function UnstableInt

!ccccccccccccccc

  real (dp) function ModInt(x)
    real (dp), intent(in) :: x

    if ( .not. present(tau2) ) then
      ModInt = NoMod(p2 - x) * ModList(i)%ShapeFun(0, x)
    else
      ModInt = NoMod(p - x, p2 - x) * ModList(i)%ShapeFun(0, x) ! the smart way does not work...
    end if

  end function ModInt

!ccccccccccccccc

  real (dp) function ModCumInt(x)
    real (dp), intent(in) :: x

    if ( .not. present(tau2) ) then
      ModCumInt = NoMod(x) * ModList(i)%ShapeFun(-1, p2 - x)
    else
      ModCumInt = NoMod(x) * ModList(i)%ShapeFun(-1, p - x, p2 - x)
    end if

  end function ModCumInt

!ccccccccccccccc

  recursive real (dp) function NoMod(x, x2) result(res)
    real (dp)            , intent(in) :: x
    real (dp)  , optional, intent(in) :: x2
    real (dp), dimension(0:2 * order) :: kerMatrix
    real (dp), dimension(0:2 * order) :: kerSingle
    integer                           :: ii
    real (dp)                         :: y

    res = 0

    if ( present(x2) ) then
      if (x2 > x) res = NoMod(x2) - NoMod(x); return
    end if

    kermatrix = 0; res = 0; kerSingle = 0; y = 0.5_dp + Atan(x/self%width)/Pi

    derInvList(0,2:,:) = KernelWidth(2*order, w, y)
    kerSingle          = SumKernelWidth(DerInvList(0,:,:), y)
    kerMatrix          = KernelMatrixList( kerSingle, self%MatAdded(:2*order,order) )
    res                = LogMatrixSum( kerMatrix, w, self%muS, x, self%width )

    if (  sum( abs(delta) ) > 0  ) then
      do ii = 1, order

        derInvList( ii,2:,:2 * (order - ii) ) = KernelWidth(2 * (order - ii), w + ii, y)
        kerSingle( :2 * (order - ii) ) = SumKernelWidth(DerInvList(ii,:,:2 * (order - ii)), y)
        kerMatrix( :2 * (order - ii) ) = KernelMatrixList( kerSingle( :2 * (order - ii) ),&
         deltaAdd(ii,:2 * (order - ii)) )

        res = res + (-self%compFact)**ii * LogMatrixSum(  kerMatrix( :2 * (order - ii) ), &
                         w + ii, self%muS, x, self%width  )

      end do
    end if

  end function NoMod

  end function SingleSingWidthList

!ccccccccccccccc

  real (dp) function HJMSing(self, c, modList, piece, setup, gap, space, cum, &
    order, isoft, s31, s32, R0, mu0, delta0, h, rho)
    class (SingularMassless)             , intent(in) :: self
    character (len = *)         , intent(in) :: setup, space, gap, cum, piece
    real (dp)                   , intent(in) :: R0, mu0, delta0, rho, h, s31, s32
    integer                     , intent(in) :: order, isoft

    real (dp)   , dimension(:,:), intent(in) :: c
    type (Model), dimension(:,:), intent(in) :: modList

    real (dp) :: srho, slog

    HJMSing = 2**( 1 + cumConst(cum) ) * self%DoubleSing(c, modList, piece, setup, &
             gap, space, cum, 'cum', order, isoft, s31, s32, R0, mu0, delta0, h, rho, rho)

    if ( setup(:2) == 'FO' ) then

      srho = self%s2rho(isoft)
      slog = self%s2log(isoft)

      if (order > 1 .and. cum(:3) == 'cum')  HJMSing = HJMSing + 4 * self%alphaS**2 * srho

      if (order > 2 .and. cum(:4) == 'diff') HJMSing = HJMSing - self%alphaS**3 *          &
      ( 32 * (slog + srho * log(rho) )/3 + (52 - 8 * self%nf/3._dp) * srho  )/rho

     if (order > 2 .and. cum(:3) == 'cum') HJMSing = HJMSing + self%alphaS**3 *            &
     ( - 8 * srho/3 + 4 * srho * (11 - 2 * self%nf/3._dp) * log(self%muS/self%Q)           &
       - 0.10280837917801414_dp * s31 + 0.6088068189625152_dp * s32 + (8 * self%nf/3._dp - &
      52) * slog - 16._dp/3 * self%s2Theta(isoft) + 8.772981689857206_dp * srho - ( 32 *      &
      slog/3 + (52 - 8 * self%nf/3._dp) * srho ) * log(rho) - 16 * srho/3 * log(rho)**2  )

   end if

  end function HJMSing

!ccccccccccccccc

  real (dp) function HJMSing1D(self, mod, gap, cum, order, isoft, s31, s32, & ! TODO: can be optimized
                                     R0, mu0, delta0, h, rho)
    class (SingularMassless)    , intent(in) :: self
    character (len = *), intent(in) :: gap, cum
    real (dp)          , intent(in) :: R0, mu0, delta0, rho, h, s31, s32
    integer            , intent(in) :: order, isoft
    type (Model)       , intent(in) :: mod

    real (dp)                       :: p, abserr
    integer                         :: neval, ier, cnt1, cnt2

    p = self%Q * rho; cnt1 = 1 + cumConst(cum); cnt2 = 1 - cumConst(cum)

    call qags( ModHJMInt, 0._dp, p, prec, prec, HJMSing1D, abserr, neval, ier )

    contains

!ccccccccccccccc

    real (dp) function ModHJMInt(x)
      real (dp)      , intent(in)     :: x

      real (dp)      , dimension(1,1) :: c
      type (Model)   , dimension(1,1) :: modList

     c(1,1) = 1; modList(1,1) = mod

      ModHJMInt = self%HJMSing(c, modList, gap, 'NoModel', gap, 'mom', cum,  &
      order, isoft, s31, s32, R0, mu0, delta0, h, x/self%Q) * mod%ShapeFun(0, p - x)

    end function ModHJMInt

  end function HJMSing1D

!ccccccccccccccc

  real (dp) function DoubleSing(self, c, modList, piece, setup, gap, space, &
    cum1, cum2, order, isoft, s31, s32, R0, mu0, delta0, h, rho1, rho2)
    class (SingularMassless), intent(in)     :: self
    character (len = *), intent(in)          :: setup, space, gap, cum1, cum2, piece
    real (dp)          , intent(in)          :: R0, mu0, delta0, rho1, rho2, h, s31, s32
    integer            , intent(in)          :: order, isoft
    real (dp)   , intent(in), dimension(:,:) :: c
    type (Model), intent(in), dimension(:,:) :: modList

    type (GapMass)  :: GapMassive
    integer         :: i, j, k, l, cnt, iNGLMax, ialpha, iNGL2, iNGL3, clen
    type (Kernels)  :: ker
    logical         :: same
    type (Model)    :: ShapeFun
    real (dp)       :: p1, p2, p, w1, w2, w, shift, NGL0, NGL1, NGL2, SumTerm, lg, facw, &
                       Gammaw1w2, rat, lp1, lp2, ratinv

    real (dp), dimension( 0:order, order, 0:2*(order - 1) ) :: deltaAdd
    real (dp), dimension( 0:order, 0:order, 0:2*order )     :: kerMatrix1, &
                                                               kerMatrix2, kerMatrix

    real (dp), dimension(2*order)            :: logList
    real (dp), dimension(0:2*isoft)          :: lgList, NGLdelta1, NGLdelta2
    real (dp), dimension(0:2*isoft, 0:2)     :: SoftNGL1, SoftNGL2, NG1, NG2
    real (dp), dimension(0:2)                :: l1List
    real (dp), dimension(0:4)                :: NG3loop
    real (dp), dimension(order,order)        :: delt
    real (dp), dimension(0:order)            :: SumKer1, SumKer2, sumkerProd
    real (dp), dimension(:,:)  , allocatable :: list1, list2
    real (dp), dimension(:)    , allocatable :: list
    real (dp), dimension(:)    , allocatable :: NGLPert
    real (dp), dimension(:,:,:), allocatable :: resul1, resul2, resulNGL1, resulNGL2, &
                                                resulNGL1Delta, resulNGL2Delta

    lgList = 0  ; NGLdelta1 = 0; NGLdelta2 = 0; l1List = 0     ; NG3loop    = 0; delt = 0
    SumKer1 = 0 ; SumKer2   = 0; shift = 0    ; sumkerProd = 0 ; kerMatrix  = 0
    SoftNGL1 = 0; SoftNGL2  = 0; NG1   = 0    ; NG2   = 0      ; DoubleSing = 0
    clen = size(c, 1)          ; iNGL2 = 0    ; iNGL3 = 0      ; kerMatrix2 = 0

    if ( order == 3 .and. abs(s31) > tiny(1._dp)  ) iNGL3 = 1; kerMatrix1 = 0; deltaAdd = 0
    if ( order == 3 .and. abs(s32) > tiny(1._dp)  ) iNGL3 = 2; if (order > 1) iNGL2 = isoft

    w1 = self%w/2 + cumConst(cum1);  w2 = self%w/2 + cumConst(cum2); NGL1 = 0;  NGL0 = 0

    iNGLMax = max(iNGL2, iNGL3); NGL2 = 0; cnt = cumConst(cum1) + cumConst(cum2)

    allocate( NGLPert(iNGLMax) ); NGLPert = 0; same = ( cum1 /= cum2 )

    select type (self)
    class is (SingularMassless)

      call self%MatEl%setGammaShift( order, self%run, gap, self%shape, delta0, R0, mu0, &
                                     h, shift, delt )
      do i = 1, order
        delt(i,:) = delt(i,:)/2**i
      end do

      p1 = self%Q * rho1 - shift/2; p2 = self%Q * rho2 - shift/2

    class is (SingularMass)

    select type (selector => self%MatEl)
      class is (MatricesElementsMass); GapMassive = GapMass(self%run, gap, selector, R0, mu0)
    end select

      call GapMassive%setGammaShift( order, self%shape, delta0, h, self%R, &
                                     self%muS, shift, delt(1,:) )

      delt(1,:) = delt(1,:)/2 + self%Q * self%deltaM(:order)

      if (order > 1) delt(2 ,2) = delt(1,1)**2/2
      if (order > 2) delt(2:,3) = [ delt(1,1) * delt(1,2), delt(1,1)**3/6 ]

      p1 = self%Q * rho1 - shift/2 - self%Q * self%tmin
      p2 = self%Q * rho2 - shift/2 - self%Q * self%tmin

    end select

    if ( gap(:6) == 'thrust' .or. gap(:6) == 'Cparam' ) then

      do i = 1, order
        do j = 1, i                 ! subtraction
          do k = 0, 2 * (i - j)     ! log power

             deltaAdd(j,i,k) = sum(   self%MatExpandedHJM( k, i - j: ceiling(k/2.): -1 ) &
                                       * delt( j, j:i - ceiling(k/2.) )   )

          end do
        end do
      end do
    end if

    if ( p1 <= 0 .or. p2 <= 0 ) return

    if ( setup(:2) == 'FO') then   ! only works if all scales are set equal. NGL only for HJM

      lg = log(p1/self%muS)

      if ( cum1(:3) == 'cum' ) then
        logList = [ ( lg**i/i, i = 1, 2*order ) ]; SumKer1 = self%MatHJMFO(0,:order)
        SumKer1 = SumKer1 + matmul( logList, self%MatHJMFO(1:2*order, :order) )
      else
        logList(1) = 1; logList(2:) = powList(lg, 2*order - 1)
        SumKer1 = matmul( logList, self%MatHJMFO(1:2*order, :order) )/rho1
      end if

      lg = log(p2/self%muS)

      if ( cum2(:3) == 'cum' ) then
        logList = [ ( lg**i/i, i = 1, 2*order ) ]; SumKer2 = self%MatHJMFO(0,:order)
        SumKer2 = SumKer2 + matmul( logList, self%MatHJMFO(1:2*order, :order) )
      else
        logList(1) = 1; logList(2:) = powList(lg, 2*order - 1)
        SumKer2 = matmul( logList, self%MatHJMFO(1:2*order, :order) )/rho2
      end if

      DoubleSing = Sum( expandProd( sumker1, sumker2) )

      return
    end if

    allocate(   list(     0:2*max(order, iNGL2 + order/3, iNGL3) )  ); list  = 0
    allocate(  list1( 0:3,0:2*max(order, iNGL2 + order/3, iNGL3) )  ); list1 = 0
    allocate(  list2( 0:3,0:2*max(order, iNGL2 + order/3, iNGL3) )  ); list2 = 0

    kernel_if: if (  .not. ( setup(:5) == 'Model' .and. space(:3) == 'pos' )  ) then

    ker = Kernels( n = 2 * max(order, iNGL2 + order/3, iNGL3 ), w = w1 )
    list1(0,:) = ker%DerList(); list2(0,:) = list1(0,:)
    if (same) ker = Kernels( n = 2*max(order, iNGL2 + order/3, iNGL3), w = w2 )
    if (same) list2(0,:) = ker%DerList()

    do i = 0, order
      kerMatrix1(0,i,:2*i) = KernelMatrixList( list1(0,:2*i), self%MatExpandedHJM(:2*i,i) )
      if (same) kerMatrix2(0,i,:2*i) = KernelMatrixList( list2(0,:2*i), self%MatExpandedHJM(:2*i,i) )
    end do

      if (isoft > 0 .and. order > 1) then

        SoftNGL1(:,0) = NGLSum( KernelsSum( list1(0,:2*isoft), list2(0,:2*isoft) ), &
        self%NGL2loop(:isoft) ); SoftNGL2(:,0) = SoftNGL1(:,0)

        if (order > 2) then

          NG3loop = NGLSum( KernelsSum( list1(0,:4), list2(0,:4) ), &
                            self%alphaS**3/32 * [s31, s32] )

          SoftNGL1(:,1) = NGLSum( KernelsSum( list1(0,:2*isoft + 1), list2(0,:2*isoft) ), &
                                self%NGL2loop(:isoft) )

          SoftNGL1(:,2) = NGLSum( KernelsSum( list1(0,:2*isoft + 2), list2(0,:2*isoft) ), &
                                self%NGL2loop(:isoft) )

          NG1 = AddBinomial(SoftNGL1, self%MatExpandedHJM(:2,1) + self%beta(0)/2 * self%alphaS * &
                                                                [0,1,0])

          SoftNGL2 =  SoftNGL1; NG2 = NG1

          if (same) then

            SoftNGL2(:,1) = NGLSum( KernelsSum( list2(0,:2*isoft + 1), list1(0,:2*isoft) ), &
                                  self%NGL2loop(:isoft) )

            SoftNGL2(:,2) = NGLSum( KernelsSum( list2(0,:2*isoft + 2), list1(0,:2*isoft) ), &
                                  self%NGL2loop(:isoft) )

            NG2 = AddBinomial(SoftNGL2, self%MatExpandedHJM(:2,1) + self%beta(0)/2 * self%alphaS * &
                                                                  [0,1,0])

          end if
        end if
      end if

      if ( gap(:6) == 'thrust' .or. gap(:6) == 'Cparam' ) then

        do i = 1, order
          list1(i,:) = DerAdd1(-1, list1(i - 1,:), w1 + i - 1 )
          list2(i,:) = DerAdd1(-1, list2(i - 1,:), w2 + i - 1 )
        end do

        if ( order > 2 .and. isoft > 0) then

          NGLdelta1 = NGLSum( KernelsSum( list1(1,:2*isoft), list2(0,:2*isoft) ), &
                            self%NGL2loop(:isoft) ); NGLdelta2 = NGLdelta1

          if (same) NGLdelta2 = NGLSum( KernelsSum( list1(0,:2*isoft), list2(1,:2*isoft) ) &
                                    , self%NGL2loop(:isoft) )
        end if

        order_do: do i = 1, order
          subtraction_do: do j = 1, i

            kerMatrix1( j, i,:2*(i - j) ) = KernelMatrixList( list1( j,:2*(i - j) ), &
                                            deltaAdd( j,i,:2*(i - j) ) )

            if (same) kerMatrix2( j, i,:2*(i - j) ) = KernelMatrixList( list2( j,:2*(i - j) ), &
                                                      deltaAdd( j,i,:2*(i - j) ) )

          end do subtraction_do
        end do order_do
      end if
    end if Kernel_if

    if (.not. same) list2 = list1

    same = ( same .or. abs(rho1 - rho2) > 0 )

    model_if: if ( setup(:7) == 'NoModel' ) then

      do i = 0, order
        sumker1(i) = LogMatrixSum( kerMatrix1(0,i,:2*i), w1, self%muS, p1 )
        if (same) sumker2(i) = LogMatrixSum( kerMatrix2(0,i,:2*i), w2, self%muS, p2 )

        if (  sum( abs(delt) ) > 0 .and. i > 0 ) then

          do j = 1, i

            sumker1(i) = sumker1(i) + (-self%compFact)**j * LogMatrixSum(  &
            kerMatrix1( j, i,:2 * (i - j) ), w1 + j, self%muS, p1  )

            if (same) sumker2(i) = sumker2(i) + (-self%compFact)**j * LogMatrixSum(  &
            kerMatrix2( j, i,:2 * (i - j) ), w2 + j, self%muS, p2  )

          end do
        end if
      end do

      if ( .not. same ) sumker2 = sumker1
      sumker1 = expandProd(sumker1, sumker2);  DoubleSing = Sum(sumker1)

      if (isoft > 0 .and. order > 1) then

        facw = self%muSEuler**(w1 + w2)/p1**(1+w1)/p2**(1+w2)

        lg = log(p2/p1); lgList(0) = 1; lgList(1:) = PowList(lg, 2*isoft)

        DoubleSing = DoubleSing + Sum( lgList * SoftNGL1(:,0) ) * facw

        if (order > 2) then

          DoubleSing = DoubleSing + Sum( lgList(:4) * NG3loop ) * facw

          if (  abs( delt(1,1) ) > 0 ) DoubleSing = DoubleSing - delt(1,1) * facw * &
                      (  Sum( lgList * NGLdelta1 )/p1 + Sum( lgList * NGLdelta2 )/p2  )

          lg = log(self%muSEuler/p1); l1List = [1._dp, lg, lg**2]

          DoubleSing = DoubleSing + Sum( lgList * matmul(NG1, l1list) ) * facw

          lg = log(self%muSEuler/p2); l1List = [1._dp, lg, lg**2]
          lgList = [ (lgList(i) * (-1)**i, i = 0, 2*isoft) ]

          DoubleSing = DoubleSing + Sum( lgList * matmul(NG2, l1list) ) * facw

        end if
      end if

      exactNGL: if (iNGL2 == -1) then ! compute the NGL integral exactly

        rat = p2/p1
        Gammaw1w2 = self%alphaS**2 * self%muSEuler**( w1 + w2 ) / &
        p1**(1 + w1 ) / p2**(1 + w2 ) / gamma( - 1 - w1 - w2 )

        NGL0 = NGLDoubleIntegral(self%nf, 0, w1, w2, rat) * Gammaw1w2

        DoubleSing = DoubleSing + NGL0

        three_loop: if (order == 3) then

          ratInv = p1/p2;  lp1 = log(self%muSEuler/p1);  lp2 = log(self%muSEuler/p2)

          NGL1 = NGLDoubleIntegral(self%nf, 1, w1, w2, rat)
          NGL2 = NGLDoubleIntegral(self%nf, 1, w2, w1, ratInv)

         DoubleSing = DoubleSing + NGL0 * (  2 * self%MatExpandedHJM(0,1) + ( self%beta(0)/2 * self%alphaS &
         + self%MatExpandedHJM(1,1) ) *  (lp1 + lp2) + 0*self%MatExpandedHJM(2,1) * (lp1**2 + lp2**2)  ) &
         + 0*self%MatExpandedHJM(2,1) * Gammaw1w2 *     &
         ( NGLDoubleIntegral(self%nf, 2, w1, w2, rat) + 2 * lp1 * NGL1 +    &
         NGLDoubleIntegral(self%nf, 2, w2, w1, ratInv) + 2 * lp2 * NGL2 ) + &
         ( self%beta(0)/2 * self%alphaS + self%MatExpandedHJM(1,1) ) * (NGL1 + NGL2) * Gammaw1w2

          delta_if: if (  abs( delt(1,1) ) > 0  ) then

             DoubleSing = DoubleSing - ( NGLDoubleIntegral(self%nf, 0, w1 + 1, &
             w2, rat) + NGLDoubleIntegral(self%nf, 0, w1, &
             w2 + 1, rat) ) * Gammaw1w2 * delt(1,1) * self%compfact

          end if delta_if
        end if three_loop
      end if exactNGL

    else if ( setup(:5) == 'Model') then

      allocate(  resul1( 0:3, clen, clen ), resul2( 0:3, clen, clen )  )
      allocate(  resulNGL1     ( 0:2*iNGLMax + 2*(order/3), clen, clen )  )
      allocate(  resulNGL2     ( 0:2*iNGLMax + 2*(order/3), clen, clen )  )
      allocate(  resulNGL1Delta( 0:2*iNGL2,                 clen, clen )  )
      allocate(  resulNGL2Delta( 0:2*iNGL2,                 clen, clen )  )

      resulNGL1 = 0; resulNGL2 = 0; resulNGL1Delta = 0; resulNGL2Delta = 0; resul1 = 0
      resul2 = 0
      if ( piece(:5) == 'piece' ) same = .false.

      iloop: do i = 1, clen
        jloop: do j = 1, i

          if ( piece(:5) == 'piece' .and. j > 0  ) exit iloop
          if ( piece(:5) == 'piece' .and. i == 1 ) then
            w1 = w2; p1 = p2; kerMatrix1 = kerMatrix2; list1 = list2
          end if

          ShapeFun  = modList(i,j); w = w1; p = p1
          kerMatrix = kerMatrix1 ; resul1(:,i,j) = VectComputer(space)

          if (same) then
             w = w2; p = p2; kerMatrix = kerMatrix2; resul2(:,i,j) = VectComputer(space)
          end if

          if (iNGLMax > 0) then

            list = list1(0,:);   w = w1; p = p1; resulNGL1(:,i,j) = NGLIntComp( space, 2*isoft + 2*(order/3) )
            if (same) then
              list = list2(0,:); w = w2; p = p2; resulNGL2(:,i,j) = NGLIntComp( space, 2*isoft + 2*(order/3) )
            end if

            if (order == 3 .and. iNGL2 > 0 .and. abs( delt(1,1) ) > 0  ) then

              list = list1(1,:); w = w1 + 1; p = p1; resulNGL1Delta(:,i,j) = NGLIntComp(space, 2*isoft)

              if (same) then
                list = list2(1,:); w = w2 + 1; p = p2; resulNGL2Delta(:,i,j) = NGLIntComp(space, 2*isoft)
              end if
            end if
          end if
        end do jloop
      end do   iloop

      call symmetrize( resul1 ); if (same) call symmetrize( resul2 )
      if (iNGLMax > 0) call symmetrize( resulNGL1 )
      if (same .and. iNGLMax > 0) call symmetrize( resulNGL2 )

      if ( iNGLMax > 0 .and. order == 3 .and. iNGL2 > 0 .and. abs( delt(1,1) ) > 0 ) then
        call symmetrize( resulNGL1Delta )
        if (same) call symmetrize( resulNGL2Delta )
      end if

      if ( .not. same ) resul2 = resul1; if ( .not. same ) resulNGL2 = resulNGL1
      if ( .not. same ) resulNGL2Delta = resulNGL1Delta

      if ( piece(:5) == 'piece' ) then
        resul2(:,0,0)    = resul1(:,1,0);   resulNGL2Delta(:,0,0) = resulNGL1Delta(:,1,0)
        resulNGL2(:,0,0) = resulNGL1(:,1,0)
      end if

      i2loop: do i = 1, clen
        do j = 1, clen
          do k = 1, clen
            do l = 1, clen

              if ( piece(:5) == 'piece' .and. l > 0 ) exit i2loop

              sumkerProd = expandProd( resul1(:,i,j), resul2(:,k,l) )
              SumTerm = Sum( sumkerProd )

              if (iNGLMax > 0) then

                NGLPert = KernelsCombine( resulNGL1(:2*iNGLMax,i,j), resulNGL2(:2*iNGLMax,k,l), iNGLMax )
                if (isoft > 0) NGL0 = Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) )
                if (isoft > 0) SumTerm = SumTerm + NGL0

                if (order == 3) then

                  if (iNGL3 >= 1) SumTerm = SumTerm + self%alphaS**3 * s31 * NGLPert(1)/32
                  if (iNGL3 == 2) SumTerm = SumTerm + self%alphaS**3 * s32 * NGLPert(2)/32

                  if (isoft > 0) then

                    NGLPert = KernelsCombine( resulNGL1(1:iNGL2+1,i,j), resulNGL2(:iNGL2,k,l), iNGL2 ) + &
                              KernelsCombine( resulNGL2(1:iNGL2+1,i,j), resulNGL1(:iNGL2,k,l), iNGL2 )
                    NGL1 = Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) )

                    NGLPert = KernelsCombine( resulNGL1(2:iNGL2+2,i,j), resulNGL2(:iNGL2,k,l), iNGL2 ) + &
                              KernelsCombine( resulNGL2(2:iNGL2+2,i,j), resulNGL1(:iNGL2,k,l), iNGL2 )
                    NGL2 = Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) )

                    SumTerm = SumTerm + 2 * NGL0 * self%MatExpandedHJM(0,1) + &
                    NGL1 * ( self%beta(0)/2 * self%alphaS + self%MatExpandedHJM(1,1) ) + &
                    NGL2 * self%MatExpandedHJM(2,1)

                    if (  abs( delt(1,1) ) > 0  ) then

                      NGLPert = KernelsCombine( resulNGL1Delta(:,i,j), resulNGL2(:iNGL2,k,l), iNGL2 ) + &
                                KernelsCombine( resulNGL2Delta(:,i,j), resulNGL1(:iNGL2,k,l), iNGL2 )

                      SumTerm = SumTerm - Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) ) * delt(1,1) * self%compfact

                    end if
                  end if
                end if
              end if

              DoubleSing = DoubleSing + c(i,k) * c(j,l) * SumTerm
              if ( piece(:5) == 'piece' ) DoubleSing = SumTerm

            end do
          end do
        end do
      end do i2loop

      deallocate(resul1, resul2, resulNGL1, resulNGL2, resulNGL1Delta, resulNGL2Delta)

    end if model_if

    DoubleSing = Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) ) * &
                 self%Prefact * DoubleSing * self%Q**(2 - cnt) * self%muSEuler**cnt

    deallocate(list1, list2, list, NGLPert)

    contains

!ccccccccccccccc

  function NGLIntComp(space, dim) result(Vec)
    character (len = *), intent(in) :: space
    integer            , intent(in) :: dim

    real (dp), dimension(0:dim)     :: Vec
    real (dp)                       :: abserr
    integer                         :: ier, neval

    Vec = 0

    do ialpha = 0, dim

      space_if: if ( space(:3) == 'mom' ) then

        call qags( MomLog, 0._dp, p, prec, prec, Vec(ialpha), abserr, neval, ier )

      else if ( space(:3) == 'pos' ) then

        call qagi( FourierL, 0._dp, 1, prec, prec, Vec(ialpha), abserr, neval, ier )
        Vec(ialpha) = Vec(ialpha)/Pi/p

      end if space_if

    end do

    if ( space(:3) == 'mom' )  Vec = KerList( list(:dim), Vec)

  end function NGLIntComp

!ccccccccccccccc

  real (dp) function MomLog(x)
    real (dp), intent(in) :: x
    real (dp)             :: arg

    arg = self%muSEuler/x

    MomLog = log(arg)**ialpha * arg**w * ShapeFun%ShapeFun(0, p - x)/x

  end function MomLog

!ccccccccccccccc

  real (dp) function FourierL(x)
    real (dp), intent(in) :: x

    complex (dp)          :: y, lgc, CI, C11, t, z, z2

    CI = (0, 1);  C11 = (1, 1); y = ( - CI/10 + C11 * x ); z2 = y/p
    t = CI * z2; z = self%muSEuler * t; lgc = log(z)**ialpha

    FourierL = real(  lgc * exp(CI * y) * C11 * z**w * ShapeFun%SoftFourier(z2)  )

  end function FourierL

!ccccccccccccccc

  function VectComputer(space) result(Vec)
    character (len = *), intent(in) :: space

    real (dp), dimension(0:3)       :: Vec
    real (dp)                       :: abserr
    integer                         :: ier, neval

    Vec = 0

    do ialpha = 0, order

      space_if: if ( space(:3) == 'mom') then

        call qags( ModVectInt, 0._dp, p, prec, prec, Vec(ialpha), abserr, neval, ier )

      else if ( space(:3) == 'pos') then

        call qagi( FourierDouble, 0._dp, 1, prec, prec, Vec(ialpha), abserr, neval, ier )
        Vec(ialpha) = Vec(ialpha)/Pi/p

      end if space_if
    end do

  end function VectComputer

!ccccccccccccccc

  real (dp) function ModVectInt(x)
    real (dp), intent(in) :: x

    integer               :: j

    ModVectInt = LogMatrixSum( kerMatrix(0,ialpha,:2*ialpha), w, self%muS, x )

    if (  sum( abs(delt) ) > 0 .and. ialpha > 0 ) then
      do j = 1, ialpha

        ModVectInt = ModVectInt + (-self%compFact)**j * LogMatrixSum(  &
        kerMatrix( j, ialpha,:2 * (ialpha - j) ), w + j, self%muS, x  )

      end do
    end if

    ModVectInt = ModVectInt * ShapeFun%ShapeFun(0, p - x)

  end function ModVectInt

!ccccccccccccccc

  real (dp) function FourierDouble(x)
    real (dp), intent(in)                   :: x

    complex (dp)                            :: y, t, z, z2, lgc, CI, C11
    complex (dp), dimension(:), allocatable :: sumtot
    integer                                 :: j

     if ( allocated(sumtot) ) deallocate(sumtot); allocate( sumtot(ialpha) ); sumtot = 0

    CI = (0, 1);  C11 = (1, 1); y = ( - CI/10 + C11 * x ); z2 = y/p
    t = CI * z2; z = self%muSEuler * t; lgc = log(z)

    sumtot = self%MatExp(0,:ialpha)/2
    if (   sum(  abs( delt(1,:ialpha) )  ) > 0   ) sumtot = sumtot - t * delt(1,:ialpha)

    do j = 1, ialpha + 1

      sumtot = sumtot + self%MatExp(j,:) * lgc**j/2

    enddo

    FourierDouble = real(  expandExpOrder(ialpha, sumtot) * exp(CI * y) * &
     C11 * z**w * ShapeFun%SoftFourier(z2)  )

     if ( allocated(sumtot) ) deallocate(sumtot)

  end function FourierDouble

  end function DoubleSing

!ccccccccccccccc

  real (dp) function DoubleSingWidth(self, c, modList, piece, setup, gap, &
    cum1, cum2, order, isoft, s31, s32, R0, mu0, delta0, h, rho1, rho2)
    class (SingularMass), intent(in)          :: self
    character (len = *) , intent(in)          :: setup, gap, cum1, cum2, piece
    real (dp)           , intent(in)          :: R0, mu0, delta0, rho1, rho2, h, s31, s32
    integer             , intent(in)          :: order, isoft

    real (dp)    , intent(in), dimension(:,:) :: c
    type (Model) , intent(in), dimension(:,:) :: modList

    type (GapMass)  :: GapMassive
    integer         :: i, j, k, l, cnt, iNGLMax, ialpha, iNGL2, iNGL3, clen
    type (Kernels)  :: ker
    logical         :: same
    type (Model)    :: ShapeFun
    real (dp)       :: p1, p2, p, w1, w2, w, shift, NGL0, NGL1, NGL2, SumTerm, lg, facw

    real (dp), dimension( 0:order, order, 0:2*(order - 1) ) :: deltaAdd
    real (dp), dimension( 0:order, 0:order, 0:2*order )     :: kerMatrix1, kerMatrix2

    real (dp), dimension(0:2*isoft)          :: lgList, NGLdelta1, NGLdelta2
    real (dp), dimension(0:2*isoft, 0:2)     :: SoftNGL1, SoftNGL2, NG1, NG2
    real (dp), dimension(0:2)                :: l1List
    real (dp), dimension(0:4)                :: NG3loop
    real (dp), dimension(order,order)        :: delt
    real (dp), dimension(0:order)            :: SumKer1, SumKer2, sumkerProd
    real (dp), dimension(:,:)  , allocatable :: list1, list2
    real (dp), dimension(:)    , allocatable :: list
    real (dp), dimension(:)    , allocatable :: NGLPert
    real (dp), dimension(:,:,:), allocatable :: resul1, resul2, resulNGL1, resulNGL2, &
                                                resulNGL1Delta, resulNGL2Delta

    lgList = 0  ; NGLdelta1 = 0; NGLdelta2 = 0; l1List = 0     ; NG3loop = 0
    SumKer1 = 0 ; SumKer2   = 0; shift = 0    ; sumkerProd = 0 ; delt    = 0
    SoftNGL1 = 0; SoftNGL2  = 0; NG1   = 0    ; NG2   = 0      ; DoubleSingWidth = 0
    clen = size(c, 1)          ; iNGL2 = 0    ; iNGL3 = 0      ; kerMatrix2 = 0

    select type (selector => self%MatEl)
      class is (MatricesElementsMass)
      GapMassive = GapMass(self%run, gap, selector, R0, mu0)
    end select

    if ( self%width <= d1mach(1) ) then
      DoubleSingWidth = self%DoubleSing(c, modList, piece, setup, gap, 'posExp', &
      cum1, cum2, order, isoft, s31, s32, R0, mu0, delta0, h, rho1, rho2);  return
    end if

    if ( order == 3 .and. abs(s31) > tiny(1._dp)  ) iNGL3 = 1; kerMatrix1 = 0; deltaAdd = 0
    if ( order == 3 .and. abs(s32) > tiny(1._dp)  ) iNGL3 = 2; if (order > 1) iNGL2 = isoft

    w1 = self%w/2 + cumConst(cum1);  w2 = self%w/2 + cumConst(cum2); NGL1 = 0;  NGL0 = 0

    iNGLMax = max(iNGL2, iNGL3); NGL2 = 0; cnt = - cumConst(cum1) - cumConst(cum2)

    allocate( NGLPert(iNGLMax) ); NGLPert = 0; same = ( cum1 /= cum2 )

   call GapMassive%setGammaShift( order, self%shape, delta0, h, self%R, &
                                     self%muS, shift, delt(1,:) )

    delt(1,:) = delt(1,:)/2 + self%Q * self%deltaM(:order)

    if (order > 1) delt(2 ,2) = delt(1,1)**2/2
    if (order > 2) delt(2:,3) = [ delt(1,1) * delt(1,2), delt(1,1)**3/6 ]

    p1 = self%Q * rho1 - shift/2 - self%Q * self%tmin
    p2 = self%Q * rho2 - shift/2 - self%Q * self%tmin

    if ( gap(:6) == 'thrust' .or. gap(:6) == 'Cparam' ) then

      do i = 1, order
        do j = 1, i                 ! subtraction
          do k = 0, 2 * (i - j)     ! log power

            deltaAdd(j,i,k) = sum(   self%MatExpandedHJM( k, i - j: ceiling(k/2.): -1 ) * &
                                   delt( j, j:i - ceiling(k/2.) )   )
          end do
        end do
      end do
    end if

    allocate(   list(     0:2*max(order, iNGL2 + order/3, iNGL3) )  ); list  = 0
    allocate(  list1( 0:3,0:2*max(order, iNGL2 + order/3, iNGL3) )  ); list1 = 0
    allocate(  list2( 0:3,0:2*max(order, iNGL2 + order/3, iNGL3) )  ); list2 = 0

    ker = Kernels( n = 2 * max(order, iNGL2 + order/3, iNGL3 ), w = w1 )
    list1(0,:) = ker%DerList(); list2(0,:) = list1(0,:)
    if (same) ker = Kernels( n = 2*max(order, iNGL2 + order/3, iNGL3), w = w2 )
    if (same) list2(0,:) = ker%DerList()

    do i = 0, order
      kerMatrix1(0,i,:2*i) = KernelMatrixList( list1(0,:2*i), self%MatExpandedHJM(:2*i,i) )
      if (same) kerMatrix2(0,i,:2*i) = KernelMatrixList( list2(0,:2*i), self%MatExpandedHJM(:2*i,i) )
    end do

      if (isoft > 0 .and. order > 1) then

        SoftNGL1(:,0) = NGLSum( KernelsSum( list1(0,:2*isoft), list2(0,:2*isoft) ), &
        self%NGL2loop(:isoft) ); SoftNGL2(:,0) = SoftNGL1(:,0)

        if (order > 2) then

          NG3loop = NGLSum( KernelsSum( list1(0,:4), list2(0,:4) ), &
                            self%alphaS**3/32 * [s31, s32] )

          SoftNGL1(:,1) = NGLSum( KernelsSum( list1(0,:2*isoft + 1), list2(0,:2*isoft) ), &
                                self%NGL2loop(:isoft) )

          SoftNGL1(:,2) = NGLSum( KernelsSum( list1(0,:2*isoft + 2), list2(0,:2*isoft) ), &
                                self%NGL2loop(:isoft) )

          NG1 = AddBinomial(SoftNGL1, self%MatExpandedHJM(:2,1) + self%beta(0)/2 * self%alphaS * &
                                                                [0,1,0])

          SoftNGL2 =  SoftNGL1; NG2 = NG1

          if (same) then

            SoftNGL2(:,1) = NGLSum( KernelsSum( list2(0,:2*isoft + 1), list1(0,:2*isoft) ), &
                                  self%NGL2loop(:isoft) )

            SoftNGL2(:,2) = NGLSum( KernelsSum( list2(0,:2*isoft + 2), list1(0,:2*isoft) ), &
                                  self%NGL2loop(:isoft) )

            NG2 = AddBinomial(SoftNGL2, self%MatExpandedHJM(:2,1) + self%beta(0)/2 * self%alphaS * &
                                                                  [0,1,0])

          end if
        end if
      end if

      if ( gap(:6) == 'thrust' .or. gap(:6) == 'Cparam' ) then

        do i = 1, order
          list1(i,:) = DerAdd1(-1, list1(i - 1,:), w1 + i - 1 )
          list2(i,:) = DerAdd1(-1, list2(i - 1,:), w2 + i - 1 )
        end do

        if ( order > 2 .and. isoft > 0) then

          NGLdelta1 = NGLSum( KernelsSum( list1(1,:2*isoft), list2(0,:2*isoft) ), &
                            self%NGL2loop(:isoft) ); NGLdelta2 = NGLdelta1

          if (same) NGLdelta2 = NGLSum( KernelsSum( list1(0,:2*isoft), list2(1,:2*isoft) ) &
                                    , self%NGL2loop(:isoft) )

        end if

        order_do: do i = 1, order
          subtraction_do: do j = 1, i

            kerMatrix1( j, i,:2*(i - j) ) = KernelMatrixList( list1( j,:2*(i - j) ), &
                                            deltaAdd( j,i,:2*(i - j) ) )

            if (same) kerMatrix2( j, i,:2*(i - j) ) = KernelMatrixList( list2( j,:2*(i - j) ), &
                                                      deltaAdd( j,i,:2*(i - j) ) )

          end do subtraction_do
        end do order_do
      end if

    if (.not. same) list2 = list1

    same = ( same .or. abs(rho1 - rho2) > 0 )

    model_if: if ( setup(:7) == 'NoModel' ) then
    ker = Kernels(n = 2 * order, width = self%width, w = w1, p = p1); list1(0,:) = ker%DerList()
    if (same) ker = Kernels(n = 2 * order, width = self%width, w = w2, p = p2); list2(0,:) = ker%DerList()

      do i = 0, order
        kerMatrix1(0,i,:2*i) = KernelMatrixList( list1(0,:2*i), self%MatExpandedHJM(:2*i,i) )
        sumker1(i) = LogMatrixSum( kerMatrix1(0,i,:2*i), w1, self%muS, p1 )
        if (same) kerMatrix2(0,i,:2*i) = KernelMatrixList( list2(0,:2*i), self%MatExpandedHJM(:2*i,i) )
        if (same) sumker2(i) = LogMatrixSum( kerMatrix2(0,i,:2*i), w2, self%muS, p2 )

        if (  sum( abs(delt) ) > 0 .and. i > 0 ) then

          do j = 1, i

            sumker1(i) = sumker1(i) + (-self%compFact)**j * LogMatrixSum(  kerMatrix1( j, i,:2 * (i - j) ), &
                         w1 + j, self%muS, p1  )

            if (same) sumker2(i) = sumker2(i) + (-self%compFact)**j * LogMatrixSum(  kerMatrix2( j, i,:2 * (i - j) ), &
                         w2 + j, self%muS, p2  )

          end do
        end if
      end do

      if ( .not. same ) sumker2 = sumker1
      sumker1 = expandProd(sumker1, sumker2); DoubleSingWidth = Sum(sumker1)

      if (isoft > 0 .and. order > 1) then

        facw = self%muSEuler**(w1 + w2)/p1**(1+w1)/p2**(1+w2)

        lg = log(p2/p1); lgList(0) = 1; lgList(1:) = PowList(lg, 2*isoft)

        DoubleSingWidth = DoubleSingWidth + Sum( lgList * SoftNGL1(:,0) ) * facw

        if (order > 2) then

          DoubleSingWidth = DoubleSingWidth + Sum( lgList(:4) * NG3loop ) * facw

          if (  abs( delt(1,1) ) > 0 ) DoubleSingWidth = DoubleSingWidth - delt(1,1) * facw * &
                      (  Sum( lgList * NGLdelta1 )/p1 + Sum( lgList * NGLdelta2 )/p2  )

          lg = log(self%muSEuler/p1); l1List = [1._dp, lg, lg**2]

          DoubleSingWidth = DoubleSingWidth + Sum( lgList * matmul(NG1, l1list) ) * facw

          lg = log(self%muSEuler/p2); l1List = [1._dp, lg, lg**2]
          lgList = [ (lgList(i) * (-1)**i, i = 0, 2*isoft) ]

          DoubleSingWidth = DoubleSingWidth + Sum( lgList * matmul(NG2, l1list) ) * facw

        end if
      end if

    else if ( setup(:5) == 'Model') then

      allocate(  resul1( 0:3, clen, clen ), resul2( 0:3, clen, clen )  )
      allocate(  resulNGL1     ( 0:2*iNGLMax + 2*(order/3),  clen, clen )  )
      allocate(  resulNGL2     ( 0:2*iNGLMax + 2*(order/3),  clen, clen )  )
      allocate(  resulNGL1Delta( 0:2*iNGL2,                  clen, clen )  )
      allocate(  resulNGL2Delta( 0:2*iNGL2,                  clen, clen )  )

      resulNGL1 = 0; resulNGL2 = 0; resulNGL1Delta = 0; resulNGL2Delta = 0; resul1 = 0
      resul2 = 0
      if ( piece(:5) == 'piece' ) same = .false.

      iloop: do i = 1, clen
        jloop: do j = 1, i

          if ( piece(:5) == 'piece' .and. j > 0  ) exit iloop
          if ( piece(:5) == 'piece' .and. i == 1 ) then
            w1 = w2; p1 = p2
          end if

          ShapeFun  = modList(i,j); w = w1; p = p1; resul1(:,i,j) = VectComputer()

          if (same) then
             w = w2; p = p2; resul2(:,i,j) = VectComputer()
          end if

          if (iNGLMax > 0) then

            list = list1(0,:);   w = w1; p = p1; resulNGL1(:,i,j) = NGLIntComp( 2*isoft + 2*(order/3) )
            if (same) then
              list = list2(0,:); w = w2; p = p2; resulNGL2(:,i,j) = NGLIntComp( 2*isoft + 2*(order/3) )
            end if

            if (order == 3 .and. iNGL2 > 0 .and. abs( delt(1,1) ) > 0  ) then

              list = list1(1,:); w = w1 + 1; p = p1; resulNGL1Delta(:,i,j) = NGLIntComp(2*isoft)

              if (same) then
                list = list2(1,:); w = w2 + 1; p = p2; resulNGL2Delta(:,i,j) = NGLIntComp(2*isoft)
              end if
            end if
          end if
        end do jloop
      end do   iloop

      call symmetrize( resul1 ); if (same) call symmetrize( resul2 )
      if (iNGLMax > 0) call symmetrize( resulNGL1 )
      if (same .and. iNGLMax > 0) call symmetrize( resulNGL2 )

      if ( iNGLMax > 0 .and. order == 3 .and. iNGL2 > 0 .and. abs( delt(1,1) ) > 0 ) then
        call symmetrize( resulNGL1Delta )
        if (same) call symmetrize( resulNGL2Delta )
      end if

      if ( .not. same ) resul2 = resul1; if ( .not. same ) resulNGL2 = resulNGL1
      if ( .not. same ) resulNGL2Delta = resulNGL1Delta

      if ( piece(:5) == 'piece' ) then
        resul2(:,0,0)    = resul1(:,1,0);   resulNGL2Delta(:,0,0) = resulNGL1Delta(:,1,0)
        resulNGL2(:,0,0) = resulNGL1(:,1,0)
      end if

      i2loop: do i = 1, clen
        do j = 1, clen
          do k = 1, clen
            do l = 1, clen

              if ( piece(:5) == 'piece' .and. l > 0 ) exit i2loop

              sumkerProd = expandProd( resul1(:,i,j), resul2(:,k,l) )
              SumTerm = Sum( sumkerProd )

              if (iNGLMax > 0) then

                NGLPert = KernelsCombine( resulNGL1(:2*iNGLMax,i,j), resulNGL2(:2*iNGLMax,k,l), iNGLMax )
                if (isoft > 0) NGL0 = Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) )
                if (isoft > 0) SumTerm = SumTerm + NGL0

                if (order == 3) then

                  if (iNGL3 >= 1) SumTerm = SumTerm + self%alphaS**3 * s31 * NGLPert(1)/32
                  if (iNGL3 == 2) SumTerm = SumTerm + self%alphaS**3 * s32 * NGLPert(2)/32

                  if (isoft > 0) then

                    NGLPert = KernelsCombine( resulNGL1(1:iNGL2+1,i,j), resulNGL2(:iNGL2,k,l), iNGL2 ) + &
                              KernelsCombine( resulNGL2(1:iNGL2+1,i,j), resulNGL1(:iNGL2,k,l), iNGL2 )
                    NGL1 = Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) )

                    NGLPert = KernelsCombine( resulNGL1(2:iNGL2+2,i,j), resulNGL2(:iNGL2,k,l), iNGL2 ) + &
                              KernelsCombine( resulNGL2(2:iNGL2+2,i,j), resulNGL1(:iNGL2,k,l), iNGL2 )
                    NGL2 = Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) )

                    SumTerm = SumTerm + 2 * NGL0 * self%MatExpandedHJM(0,1) + &
                    NGL1 * ( self%beta(0)/2 * self%alphaS + self%MatExpandedHJM(1,1) ) + &
                    NGL2 * self%MatExpandedHJM(2,1)

                    if (  abs( delt(1,1) ) > 0  ) then

                      NGLPert = KernelsCombine( resulNGL1Delta(:,i,j), resulNGL2(:iNGL2,k,l), iNGL2 ) + &
                                KernelsCombine( resulNGL2Delta(:,i,j), resulNGL1(:iNGL2,k,l), iNGL2 )

                      SumTerm = SumTerm - Sum( self%NGL2loop(:isoft) * NGLPert(:isoft) ) * delt(1,1) * self%compfact

                    end if
                  end if
                end if
              end if

              DoubleSingWidth = DoubleSingWidth + c(i,k) * c(j,l) * SumTerm
              if ( piece(:5) == 'piece' ) DoubleSingWidth = SumTerm

            end do
          end do
        end do
      end do i2loop

      deallocate(resul1, resul2, resulNGL1, resulNGL2, resulNGL1Delta, resulNGL2Delta)

    end if model_if

    DoubleSingWidth = Sum( self%HardExp(:order) ) * Sum( self%HardMassExp(:order) ) * &
                      self%Prefact * DoubleSingWidth * self%Q**(2 - cnt) * self%muSEuler**cnt

    deallocate(list1, list2, list, NGLPert)

    contains

!ccccccccccccccc

  function NGLIntComp(dim) result(Vec)
    integer            , intent(in) :: dim

    real (dp), dimension(0:dim)     :: Vec
    real (dp)                       :: abserr
    integer                         :: ier, neval

    Vec = 0

    do ialpha = 0, dim
      call qagi( MomLog, 0._dp, 1, prec, prec, Vec(ialpha), abserr, neval, ier )
    end do

    Vec = KerList( list(:dim), Vec)

  end function NGLIntComp

!ccccccccccccccc

  real (dp) function MomLog(x)
    real (dp), intent(in) :: x
    real (dp)             :: arg

    arg = self%muSEuler/x

    MomLog = log(arg)**ialpha * arg**w * ShapeFun%ShapeFun(0, p - x)/x

  end function MomLog

!ccccccccccccccc

  function VectComputer() result(Vec)

    real (dp), dimension(0:3)       :: Vec
    real (dp)                       :: abserr
    integer                         :: ier, neval

    Vec = 0

    do ialpha = 0, order
      call qagi( ModVectInt, p, -1, prec, prec, Vec(ialpha), abserr, neval, ier )
    end do

  end function VectComputer

!ccccccccccccccc

  real (dp) function ModVectInt(x)
    real (dp), intent(in) :: x

    integer               :: j

    ker = Kernels(n = 2 * order, width = self%width, w = w, p = x); list1(0,:) = ker%DerList()
    kerMatrix1(0,ialpha,:2*ialpha) = KernelMatrixList( list1(0,:2*ialpha), self%MatExpandedHJM(:2*ialpha,ialpha) )

    ModVectInt = LogMatrixSum( kerMatrix1(0,ialpha,:2*ialpha), w, self%muS, x )

    if (  sum( abs(delt) ) > 0 .and. ialpha > 0  ) then
      do j = 1, ialpha

        ModVectInt = ModVectInt + (-self%compFact)**j * LogMatrixSum(  &
        kerMatrix1( j, ialpha,:2 * (ialpha - j) ), w + j, self%muS, x  )

      end do
    end if

    ModVectInt = ModVectInt * ShapeFun%ShapeFun(0, p - x)

  end function ModVectInt

  end function DoubleSingWidth

!ccccccccccccccc

  function NGLIntComputer(self, space, isoft, pow, wSoft, mod, p) result(Vec)
    character (len = *)             :: space
    class (SingularMassless), intent(in)     :: self
    integer        , intent(in)     :: isoft, pow
    real (dp)      , intent(in)     :: wsoft, p
    type(Model)    , intent(in)     :: Mod

    real (dp), dimension(0:2*isoft) :: Vec
    real (dp)                       :: abserr, l
    integer                         :: ier, neval, i, powLog

    Vec = 0; l = p

    do i = 0, 2 * isoft

      powLog = i + pow

      space_if: if ( space(:3) == 'mom') then

        call qags( MomentumLog, 0._dp, p, prec, prec, Vec(i), abserr, neval, ier )

      else if ( space(:3) == 'pos') then

        call qagi( FourierLog, 0._dp, 1, prec, prec, Vec(i), abserr, neval, ier )
        Vec(i) = Vec(i)/Pi

      end if space_if

    end do

    contains

!ccccccccccccccc

  real (dp) function FourierLog(x)
    real (dp), intent(in) :: x

    complex (dp)          :: y, lgc, CI, C11

    CI = (0, 1);  C11 = (1, 1)

    y = ( - CI/10 + C11 * x )/l; lgc = log( CI * self%muSEuler * y )**powLog

    FourierLog = real(  lgc * exp( CI * l * y ) * &
     C11 * ( CI * self%muSEuler * y )**wSoft * mod%SoftFourier(y)  )/l

  end function FourierLog

!ccccccccccccccc

  real (dp) function MomentumLog(x)
    real (dp), intent(in)            :: x

    type (kernels)                   :: ker
    real (dp), dimension( 0:powLog ) :: kerlista

    ker = Kernels(n = powLog, w = wSoft, mu = self%muS, p = x); kerlista = ker%KernelList()

    MomentumLog = kerlista(powlog) * mod%ShapeFun(0, l - x)

  end function MomentumLog

  end function NGLIntComputer

!ccccccccccccccc

 function ExpMat(self) result(mat)
    class (SingularMassless), intent(in) :: self
    real (dp), dimension(0:4,3) :: mat
    mat = self%MatExp
  end function

!ccccccccccccccc

  real (dp) function Prefactor(self)
    class (SingularMassless), intent(in) :: self
    Prefactor = self%preFact
  end function Prefactor

!ccccccccccccccc

  integer function Resum(self)
    class (SingularMassless), intent(in) :: self
    Resum = self%run
  end function Resum

!ccccccccccccccc

  real (dp) function wTilde(self)
    class (SingularMassless), intent(in) :: self
    wTilde = self%w
  end function wTilde

!ccccccccccccccc

  integer function cumConst(cum)
    character (len = *), intent(in) :: cum

    if ( cum(:3) == 'cum' ) then
      cumConst = -1
    else if ( cum(:4) == 'diff' ) then
      cumConst = 0
    else
      read (cum(:1),*) cumConst
    end if

  end function cumConst

!ccccccccccccccc

  subroutine symmetrize(mat)
    real (dp), intent(inout), dimension(:,:,:) :: mat
    integer                                    :: i, j, len

    len = size(mat,2)

    do i = 1, len
      do j = i + 1, len
        mat(:,i,j) = mat(:,j,i)
      end do
    end do
  end subroutine symmetrize

!ccccccccccccccc

  real (dp) function s2rho(self, isoft)
    class (SingularMassless), intent(in) :: self
    integer, intent(in)         :: isoft

      select case(isoft)
      case(-1)
       s2rho = 5.817594867359059_dp + 0.2895745232754116_dp  * self%nf
    case(1)
       s2rho = 6.234083201979812_dp + 0.2741556778080381_dp  * self%nf
    case(2)
       s2rho = 5.873308790742764_dp + 0.28918794494291505_dp * self%nf
    case(3)
       s2rho = 5.825832781183491_dp + 0.2898944331804042_dp  * self%nf
    case(4)
       s2rho = 5.8188363048966805_dp + 0.2897145802589072_dp * self%nf
    case(5)
       s2rho = 5.817775727311131_dp + 0.28961899917968487_dp * self%nf
    case(6)
       s2rho = 5.817618283362493_dp + 0.28958722156747796_dp * self%nf
    case default
       s2rho = 0
    end select

    s2rho = s2rho/4

  end function s2rho

!ccccccccccccccc

  real (dp) function s2log(self, isoft)
    class (SingularMassless), intent(in) :: self
    integer, intent(in)         :: isoft

      select case(isoft)
      case(-1)
       s2log = - 4.054190364656132_dp - 0.2182114365222843_dp  * self%nf
    case(1)
       s2log = - 4.5556371521743655_dp - 0.2003428171932647_dp * self%nf
    case(2)
       s2log = - 4.137691998772623_dp - 0.2177571985850039_dp  * self%nf
    case(3)
       s2log = - 4.068393349247342_dp - 0.21878842848865399_dp * self%nf
    case(4)
       s2log = - 4.056556964310233_dp - 0.21848415983631592_dp * self%nf
    case(5)
       s2log = - 4.054562606594437_dp - 0.21830442488299617_dp * self%nf
    case(6)
       s2log = - 4.054241525672442_dp - 0.21823961969286576_dp * self%nf
    case default
       s2log = 0
    end select

    s2log = s2log/4

  end function s2log


!ccccccccccccccc

  real (dp) function s2theta(self, isoft)
    class (SingularMassless), intent(in) :: self
    integer, intent(in)         :: isoft

    select case(isoft)
      case(-1)
       s2theta = 1.5905128313389865_dp + 0.10549437091963328_dp * self%nf
    case(1)
       s2theta = 2.0509311668995327_dp + 0.09019360280921593_dp * self%nf
    case(2)
       s2theta = 1.6913305514803658_dp + 0.10517696178501458_dp * self%nf
    case(3)
       s2theta = 1.6108837253562018_dp + 0.106374087173767_dp   * self%nf
    case(4)
       s2theta = 1.5943476166811446_dp + 0.10594900641404088_dp * self%nf
    case(5)
       s2theta = 1.591174171979759_dp + 0.10566301011019416_dp  * self%nf
    case(6)
       s2theta = 1.5906104121095412_dp + 0.10554922394890451_dp * self%nf
    case default
       s2theta = 0
    end select

    s2theta = s2theta/4

  end function s2theta

!ccccccccccccccc

  function SoftMatchingExp(self) result(SoftMatExp) ! Issue with mass scheme
    class (SingularMass) , intent(in) :: self

    real (dp)                         :: LQ, Lm
    integer                           :: nl
    real (dp), dimension(0:4,3)       :: SoftMatExp

    nl = self%nf - 1;  Lm = Log(self%mass/self%muM); LQ = log(self%muS/self%muM)

    SoftMatExp = 0

    SoftMatExp(0,2) = - (56._dp/9 + 40 * Lm/3 + 8 * Lm**2)/9 * self%alphaM**2
    SoftMatExp(0,1) = - SoftMatExp(0,2) * LQ

    SoftMatExp(1,2) = self%alphaM**2 * ( - 20 * Lm**2 - 8 * Lm**3 + Lm * (Pi2 -   &
     112._dp/3) - 82._dp/3 + 5 * Pi2/12 + 7 * Zeta3 )/27 + self%alphaM**3 * LQ *  &
     ( (16 * nl - 232) * Lm**3/9 + Lm**2 * (238._dp/9 - 4 * Pi2) + Lm *           &
     (625/27._dp + 136 * nl/27._dp - 20 * Pi2/3 + 52 * Zeta3) )/9

  end function SoftMatchingExp

!ccccccccccccccc

  real (dp) function FQCD(self)
    use Poly
    class (SingularMass), intent(in) :: self
    real (dp)                        :: root, mx, lmx, DiLog, Li2, Li3, fac1, fac2, mx2

    root = Sqrt(1 + 4 * self%m2); fac1 = 19 + 46 * self%m2; fac2 = 1 - 6 * self%m2**2
    mx  = (root - 1)/2/self%m; lmx = log(mx); mx2 = mx**2
    Li2 = DiLog(mx2);  Li3 = PolyLog(3, mx2)

    FQCD = 5 * self%Lm/324 * (53 + 132 * self%m2) + Lmx**2 * root * fac1/27 +        &
           2 * Lmx * fac2 * (2 * Lmx**2 - Pi2)/27 + 3355/1944._dp + 119 * self%m2/27 &
           - root * fac1 * (Pi2 - 6 * Li2)/162 + 2 * fac2 * (Li3 - Zeta3)/9

  end function FQCD

!ccccccccccccccc

  real (dp) function FSCET(self)
    class (SingularMass), intent(in) :: self

    FSCET = (  875 + 48 * self%Lm**2 * (7 - 3 * self%LQ) + 48 * Pi2 - 360 * Zeta3  + &
               4 * self%LQ * (  65 + 12 * (2 - 3 * self%LQ) * self%LQ + 3 * Pi2 )  + &
               4 * self%Lm * ( 121 + 36 * (3 - 2 * self%LQ) * self%LQ + 6 * Pi2 )  )/1296

  end function FSCET

!ccccccccccccccc

  real (dp) function FMS(self)
    class (SingularMass), intent(in) :: self

    FMS = (  3 * Pi2 - 4 * self%Lm**3 + 6 * self%Lm**2 * (3 - 4 * self%LQ) - 8 * Zeta3 + &
          12 * self%Lm * (     Pi2 -  8 + ( 6 -  4 * self%LQ) * self%LQ ) +              &
           2 * self%LQ * ( 5 * Pi2 - 48 + (27 - 14 * self%LQ) * self%LQ )  )/216

  end function FMS

!ccccccccccccccc

  real (dp) function JetNonDist(self, cum, s)
    class (SingularMass), intent(in) :: self
    character (len = *) , intent(in) :: cum
    real (dp)           , intent(in) :: s

    real (dp)                        :: dilog

    JetNonDist = 0; if (s < 0) return

    if ( self%EShape(:1) == 'Q' ) then

      if ( cum(:4) == 'diff' ) then

        JetNonDist = ( s/(self%m2 + s)**2 - 4/s * log(1 + s/self%m2) )/3

      else if ( cum(:3) == 'cum' ) then

        JetNonDist = ( log(1 + s/self%m2) -s/(self%m2 + s) + 4 * dilog(-s/self%m2) )/3

      end if

    else if ( self%EShape(:1) == 'P' ) then

      if ( cum(:4) == 'diff' ) then

        JetNonDist = ( (s - 7 * self%m2)/(s - self%m2)**2 - 2 * (2 * s - 5 * self%m2) &
        * s * log(s/self%m2)/(s - self%m**3)**2 )/3

      else if ( cum(:3) == 'cum' ) then

        JetNonDist = (  3 * s/(s - self%m2) - 4 * dilog(-s/self%m2) +  ( s * (s - 4 * &
        self%m2) - 4 * (s - self%m2)**2 * log(1 - s/self%m2) ) * log(s/self%m2)  )/3

      end if
    end if

  end function JetNonDist

!ccccccccccccccc

  real (dp) function IntMassJet(self, p, w)
    class (SingularMassless), intent(in) :: self
    real (dp)       , intent(in) :: p, w

    real (dp)                    :: abserr
    integer                      :: neval, ier

    IntMassJet = 0; if (p <= 0) return

    if ( self%EShape(:1) == 'Q' ) then

      IntMassJet = p/(1 + p)**2/(1 - w) * (1 - w**2 + (w - p) * w * &
      Hyper2F1(1._dp, 1._dp, 2 - w, -p) ) - 4 * H3F2Exact(w, p)

    else

      call qags( ThrustJet, 0._dp, 1._dp, prec, prec, IntMassJet, abserr, neval, ier )

      IntMassJet = p * IntMassJet - 7

    end if

    contains

    real (dp) function ThrustJet(y)
      real (dp), intent(in) :: y

      real (dp)             :: x

       x = p * y

       if (self%EShape(:1) == 'Q') then
         ThrustJet = ( - x *  (4 + 7 * x + 5 * x**2) +    &
         4 * (1 + x)**3 * log(1 + x) )/x**2/(1 + x)**3
       else
         ThrustJet =  ( 28 * x - 23 - 5 * x**2 +  &
         2 * (2 * x**2 - 5 - 6 * x) * log(x) )/(x - 1)**4
       end if

       ThrustJet = (1 - y)**(- w) * ThrustJet

    end function ThrustJet

  end function IntMassJet

!ccccccccccccccc

  pure integer function shapeFactor(str)
    character (len = *), intent(in) :: str

    shapeFactor = 1; if ( str(:6) == 'Cparam' ) shapeFactor = 6

  end function shapeFactor

!ccccccccccccccc

end module SingularClass
