
module MatrixElementsClass
  use AnomDimClass;  use RunningClass;  use AlphaClass   ; use Chaplin, only: CLi2, CLi3
  use Constants, only: dp, Pi, Pi2, Pio2, Zeta3, ExpEuler, prec; use QuadPack, only: qags
  implicit none;  private

  public :: expandExp, expandExpOrder, cuspConst, expandProd, expandExpVec, gapCons, &
  posToMomMatrix, NGLSoft, NGLDoubleIntegral, factList

  interface expandExpOrder
    module procedure :: expandExpComplexOrder, expandExpRealOrder
  end interface expandExpOrder

  interface expandExp
    module procedure :: expandExpComplex, expandExpReal, expandExpMat
  end interface expandExp

!ccccccccccccccc

  type, public                  :: MatricesElements
    private
    type (AnomDim)              :: andim
    integer                     :: nf
    type (Running)              :: alphaMass
    type (Alpha)                :: AlphaAll
    real (dp)                   :: muH, Q, R, muS, muJ
    real (dp), dimension(6)     :: NGlist, alphaList
    real (dp), dimension(0:4  ) :: ListFact
    real (dp), dimension(0:4,3) :: softExpT, softExpC, jetExp, HardExp, coefCusp, coefHard
    real (dp), dimension(0:4)   :: beta
    real (dp), dimension(3)     :: deltaT, sCoefT, CparamGammaR, deltaTHadron, &
    ThrustGammaR, sCoefC, deltaCHadron, deltaC

    contains

    procedure, pass(self), private :: SoftNGLThetaApprox, SetHardFunc
    procedure, pass(self), public  :: CoefMat, delta, sCoef, DiffDeltaGap, adim,  &
    DiffDeltaGapHadron, NonGlobalList, NGLIntegral, run, alphaScale, NGLfunction, &
    scales, SetDelta, AddLogs, SetScales, SetAlpha, SetMTop, SetMBottom, SetHard, &
    setGammaShift, SetMCharm, GammaR

  end type MatricesElements

!ccccccccccccccc

  type, extends (MatricesElements), public  :: MatrixElements

  end type MatrixElements

!ccccccccccccccc

  type, extends (MatricesElements), public  :: MatricesElementsMass
    private
    character (len = 5)             :: scheme
    type (AnomDim)                  :: andimNl
    type (Running)                  :: alphaMassNl
    real (dp)                       :: Rmass, mm, muM
    integer                         :: nl
    real (dp), dimension(0:4,3)     :: BJetExp, softExpTNl, softExpCNl, HmExp
    real (dp), dimension(0:3)       :: jetMatching
    real (dp), dimension(0:4)       :: betaNl
    real (dp), dimension(6)         :: NGlistNl
    real (dp), dimension(3)         :: sCoefTNl, sCoefCNl, sCoefbJet, alphaMMlist, &
    CparamGammaRNl, ThrustGammaRNl, JetGammaR

    contains

    procedure, pass(self), public   :: MassDeltaSet, DiffDeltaGapNl, setMass, &
    runNl, NonGlobalListNl, SetHardMass, HardMassVect, HardMassExp, HardMassVec,&
    JetMass, mmFromJetMass, sCoefNl

  end type MatricesElementsMass

!ccccccccccccccc

  type, extends (MatricesElementsMass), public  :: MatrixElementsMass
    private
    type (AnomDim)                              :: andimS
    type (Running)                              :: alphaMassS
    real (dp)                                   :: MSMass, MSLow, alphaB
    real (dp), dimension(3)                     :: DeltaMJet
    real (dp), dimension(4)                     :: deltaMSRNatural, DeltaMSR, &
    deltaMSLow, deltaMS

    contains

    procedure, pass(self), public :: MassDeltaScheme

  end type MatrixElementsMass

!ccccccccccccccc

  interface MatricesElementsMass
    module procedure InMatMass
  end interface MatricesElementsMass

!ccccccccccccccc

  interface MatricesElements
    module procedure InMatrices
  end interface MatricesElements

!ccccccccccccccc

  interface MatrixElements
    module procedure InMatEl
  end interface MatrixElements

!ccccccccccccccc

  interface MatrixElementsMass
    module procedure InMatElMass
  end interface MatrixElementsMass

  contains

!ccccccccccccccc

  type (MatricesElements) function InMatrices(AlphaAll, nf, s3T, s3C, j3, muLambda)
    type (Alpha)   , intent(in) :: AlphaAll
    integer        , intent(in) :: nf
    real (dp)      , intent(in) :: s3T, s3C, j3, muLambda
    integer                     :: i
    type (Running)              :: alphaMass
    real (dp), dimension(0:4)   :: ListFact, beta, betaList
    real (dp), dimension(0:4,3) :: coefCusp, coefSoft, coefJet, coefHard
    integer  , dimension(3)     :: TwoList
    type (AnomDim)              :: andim

    alphaMass = Running(nf, 0, AlphaAll, muLambda) ; coefCusp = 0; ListFact = factList(4)

    andim     = alphaMass%adim() ; InMatrices%nf = nf ; InMatrices%ListFact = ListFact
    InMatrices%andim = andim     ; beta = andim%betaQCD('beta');  InMatrices%beta = beta

  ! Expand Cusp anomalous dimension for nf and nl flavors

    call andim%kTildeExpand( andim%betaQCD('cusp'), coefCusp )
    InMatrices%coefCusp = coefCusp

! Assign all constants

    coefJet = 0; coefHard = 0;  coefSoft = 0;  InMatrices%alphaMass = alphaMass

  ! Expand anomalous dimension

    call andim%wTildeExpand( andim%betaQCD('hard') , coefHard )
    call andim%wTildeExpand( andim%betaQCD('soft') , coefSoft )
    call andim%wTildeExpand( andim%betaQCD('jet')  , coefJet  )

    InMatrices%coefHard = coefHard

! Initialize everything to zero

    InMatrices%softExpT = 0 ;  InMatrices%softExpC = 0;  InMatrices%HardExp = 0
    InMatrices%jetExp   = 0 ;  InMatrices%AlphaAll = AlphaAll

! Non-Global 2-loop soft function

    InMatrices%NGlist = NGLOas2(nf)

! Assign non-log matrix element terms

    InMatrices%softExpT(0,:) = GenMat('SoftThrust', nf, s3T)
    InMatrices%softExpC(0,:) = GenMat('SoftCparam', nf, s3C)
    InMatrices%jetExp(0,:)   = GenMat('jet' , nf, j3)
    InMatrices%HardExp(0,:)  = GenMat('hard', nf, 0._dp)

! Compute full matrix elements:
  ! Expand alphaS

    call andim%expandAlpha(InMatrices%softExpT)
    call andim%expandAlpha(InMatrices%softExpC)
    call andim%expandAlpha(InMatrices%jetExp)
    call andim%expandAlpha(InMatrices%HardExp)

  ! Add alphaQCD expansion + anomalous dimension expansion

    InMatrices%softExpT = InMatrices%softExpT - 2 * coefCusp + coefSoft
    InMatrices%softExpC = InMatrices%softExpC - 2 * coefCusp + coefSoft
    InMatrices%jetExp   = InMatrices%jetExp   + 2 * coefCusp + coefJet/2

    betaList = andim%betaQCD('betaList')

    InMatrices%ThrustGammaR = ExpEuler * andim%GammaRComputer( InMatrices%softExpT(1,:) )
    InMatrices%CparamGammaR = ExpEuler * andim%GammaRComputer( InMatrices%softExpC(1,:) )

    InMatrices%sCoefT =  andim%sCoef( InMatrices%ThrustGammaR * betaList(1:3) )
    InMatrices%sCoefC =  andim%sCoef( InMatrices%CparamGammaR * betaList(1:3) )

 ! rescale jet logs because of their mass dimension

    TwoList = powList(2,4);

    do i = 1, 4
      InMatrices%jetExp(i,:) = InMatrices%jetExp(i,:)/TwoList(i)
    end do

   end function InMatrices

!ccccccccccccccc

  type (MatricesElementsMass) function InMatMass(AlphaAll, nf, runMass, s3T, s3C, j3, b3, &
                                                  muLambda, muLambdaNl)
    type (Alpha)   , intent(in)   :: AlphaAll
    integer        , intent(in)   :: nf, runMass
    real (dp)      , intent(in)   :: s3T, s3C, j3, b3, muLambda, muLambdaNl
    integer                       :: i, nl
    type (Running)                :: alphaMass, alphaMassNl
    real (dp), dimension(0:4)     :: ListFact, beta, betaList, betaListNl
    integer  , dimension(3)       :: TwoList
    type (AnomDim)                :: andim, andimNl
    real (dp)                     :: mm
    real (dp), dimension(0:3,0:3) :: tab
    real (dp), dimension(3)       :: a, b
    real (dp), dimension(4)       :: c
    real (dp), dimension(0:4,3)   :: coefCusp, coefCuspNl, coefSoft, coefJet, &
    coefSoftNl, coefBJet, coefHard

    alphaMass = Running(nf, runMass, AlphaAll, muLambda); ListFact = factList(4)
    nl = nf - 1; alphaMassNl = Running(nl, runMass, AlphaAll, muLambdaNl); coefCusp = 0
    coefCuspNl = 0; coefSoftNl = 0; mm  = alphaMass%scales('mL'); InMatMass%mm = mm
    InMatMass%scheme = andim%scheme();  InMatMass%nl = nl

    andim   = alphaMass%adim()  ; InMatMass%nf = nf      ; InMatMass%ListFact = ListFact
    andimNl = alphaMassNl%adim(); InMatMass%andim = andim; InMatMass%andimNl  = andimNl

    InMatMass%alphaMMlist = powList(alphaMass%alphaQCD(mm)/Pi,3)
    InMatMass%alphaMassNl = alphaMassNl

  ! Expand Cusp anomalous dimension for nf and nl flavors

    call   andim%kTildeExpand( andim%betaQCD('cusp')  , coefCusp   )
    call andimNl%kTildeExpand( andimNl%betaQCD('cusp'), coefCuspNl )
    InMatMass%coefCusp = coefCusp

    betaListNl = andimNl%betaQCD('betaList'); InMatMass%betaNl = andimNl%betaQCD('beta')

! compute Logs

    beta = andim%betaQCD('beta'); InMatMass%beta = beta

    coefJet = 0; coefHard = 0;  coefSoft = 0

  ! Expand anomalous dimension

    call andim%wTildeExpand( andim%betaQCD('hard')    , coefHard )
    call andim%wTildeExpand( andim%betaQCD('soft')    , coefSoft )
    call andimNl%wTildeExpand( andimNl%betaQCD('soft'), coefSoftNl )
    call andim%wTildeExpand( andim%betaQCD('jet')     , coefJet  )
    call andimNl%wTildeExpand( andimNl%betaQCD('bJet'), coefBJet )

    InMatMass%coefHard = coefHard;  InMatMass%alphaMass = alphaMass

! Initialize everything to zero

    InMatMass%softExpT = 0 ;  InMatMass%softExpC = 0;  InMatMass%HardExp = 0
    InMatMass%jetExp   = 0 ;  InMatMass%AlphaAll = AlphaAll

! Non-Global 2-loop soft function

    InMatMass%NGlist = NGLOas2(nf);  InMatMass%NGlistNl = NGLOas2(nl)

! Assign non-log matrix element terms

    InMatMass%softExpTNl(0,:) = GenMat('SoftThrust', nl, s3T)
    InMatMass%softExpT(0,:)   = GenMat('SoftThrust', nf, s3T)
    InMatMass%softExpC(0,:)   = GenMat('SoftCparam', nf, s3C)
    InMatMass%softExpCNl(0,:) = GenMat('SoftCparam', nl, s3C)
    InMatMass%jetExp(0,:)     = GenMat('jet' , nf, j3)
    InMatMass%HardExp(0,:)    = GenMat('hard', nf, 0._dp)
    InMatMass%BJetExp(0,:)    = GenMat('bJet', nl, b3)

! Compute full matrix elements:
  ! Expand alphaS

    call andim%expandAlpha(InMatMass%softExpT)
    call andim%expandAlpha(InMatMass%softExpC)
    call andimNl%expandAlpha(InMatMass%softExpTNl)
    call andimNl%expandAlpha(InMatMass%softExpCNl)
    call andim%expandAlpha(InMatMass%jetExp)
    call andim%expandAlpha(InMatMass%HardExp)
    call andimNl%expandAlpha(InMatMass%BJetExp)

  ! Add alphaQCD expansion + anomalous dimension expansion

    InMatMass%softExpT = InMatMass%softExpT - 2 * coefCusp + coefSoft
    InMatMass%softExpC = InMatMass%softExpC - 2 * coefCusp + coefSoft
    InMatMass%jetExp   = InMatMass%jetExp   + 2 * coefCusp + coefJet/2
    InMatMass%BJetExp  = InMatMass%BJetExp  +   coefCuspNl + coefBJet/2

    InMatMass%softExpTNl = InMatMass%softExpTNl - 2 * coefCuspNl + coefSoftNl
    InMatMass%softExpCNl = InMatMass%softExpCNl - 2 * coefCuspNl + coefSoftNl

    betaList = andim%betaQCD('betaList')

    InMatMass%ThrustGammaR = ExpEuler * andim%GammaRComputer( InMatMass%softExpT(1,:) )
    InMatMass%CparamGammaR = ExpEuler * andim%GammaRComputer( InMatMass%softExpC(1,:) )
    InMatMass%JetGammaR    = ExpEuler * andimNl%GammaRComputer( InMatMass%BJetExp(1,:) )/2
    InMatMass%ThrustGammaRNl = ExpEuler * andimNl%GammaRComputer( InMatMass%softExpTNl(1,:) )
    InMatMass%CparamGammaRNl = ExpEuler * andimNl%GammaRComputer( InMatMass%softExpCNl(1,:) )

    InMatMass%sCoefT    =  andim%sCoef( InMatMass%ThrustGammaR * betaList(1:3) )
    InMatMass%sCoefC    =  andim%sCoef( InMatMass%CparamGammaR * betaList(1:3) )
    InMatMass%sCoefbJet = andimNl%sCoef( andimNl%GammaRComputer( InMatMass%BJetExp(1,:) ) * betaListNl(1:3) ) * ExpEuler/2
    InMatMass%sCoefTNl  =  andimNl%sCoef( InMatMass%ThrustGammaRNl * betaListNl(1:3) )
    InMatMass%sCoefCNl  =  andimNl%sCoef( InMatMass%CparamGammaRNl * betaListNl(1:3) )

 ! rescale jet logs because of their mass dimension

    TwoList = powList(2,4);

    do i = 1, 4
      InMatMass%jetExp(i,:) = InMatMass%jetExp(i,:)/TwoList(i)
    end do

     InMatMass%alphaList = 0; InMatMass%jetMatching(0) = 1

     a = ExpEuler * DeltaComputer(InMatMass%BJetExp(1:,:), [0,0,0] * 0._dp, 1)/2
     c = MSbarDelta(nf - 1, 1); tab = andim%alphaMatching(nf)
     b = tab(:2,0); call alphaReExpand(a,b);  a = c(:3) - a
     a = a * InMatMass%alphaMMlist

     InMatMass%jetMatching(1:) = a

   end function InMatMass

!ccccccccccccccc

  type (MatrixElements) function InMatEl(AlphaAll, nf, s3T, s3C, j3, Q, muH, muJ, &
  muS, R, muLambda)
    type (Alpha)        , intent(in) :: AlphaAll
    integer             , intent(in) :: nf
    real (dp)           , intent(in) :: s3T, s3C, j3, muLambda
    real (dp), optional , intent(in) :: muJ, muS, R, Q, muH
    integer                          :: i
    type (Running)                   :: alphaMass
    real (dp), dimension(0:4)        :: beta, betaList
    real (dp), dimension(0:4)        :: ListFact
    real (dp), dimension(3)          :: alphaSList, alphaJList
    integer  , dimension(3)          :: TwoList
    type (AnomDim)                   :: andim
    real (dp), dimension(0:4,3)      :: coefCusp, coefSoft, coefJet, coefHard
    real (dp)                        :: alphaJ, alphaS, alphaR

    alphaMass = Running(nf, 0, AlphaAll, muLambda) ; coefCusp = 0

    InMatEl%AlphaAll = AlphaAll  ;  InMatEl%nf  = nf  ; ListFact = factList(4)
    InMatEl%ListFact = ListFact  ;  andim = alphaMass%adim()

  ! Expand Cusp anomalous dimension for nf and nl flavors

    call andim%kTildeExpand( andim%betaQCD('cusp'), coefCusp )

    InMatEl%coefCusp = coefCusp

! Assign all constants

    if ( present(muS) ) then
      alphaS = alphaMass%alphaQCD(muS)/Pi; alphaSList = powList(alphaS, 3)
      InMatEl%muS = muS
    else
      alphaS = 1; alphaSlist = 1; InMatEl%muS = 0
    end if

    if ( present(muJ) ) then
      alphaJ = alphaMass%alphaQCD(muJ)/Pi; alphaJList = powList(alphaJ, 3)
      InMatEl%muJ = muJ
    else
      alphaJ = 1; alphaJlist = 1; InMatEl%muJ = 0
    end if

    if ( present(R  ) ) then
      alphaR = alphaMass%alphaQCD(R)/Pi; InMatEl%R = R
    else
      alphaR = 1; InMatEl%R = 0
    end if

    beta = andim%betaQCD('beta');  InMatEl%beta = beta

  ! Expand anomalous dimension

    coefJet  = 0; coefHard = 0;  coefSoft = 0

    call andim%wTildeExpand( andim%betaQCD('hard') , coefHard )
    call andim%wTildeExpand( andim%betaQCD('soft') , coefSoft )
    call andim%wTildeExpand( andim%betaQCD('jet')  , coefJet  )

    InMatEl%coefHard = coefHard

    InMatEl%alphaList(2:) = [ alphaJ, 0._dp, 0._dp, alphaS, alphaR ]

    InMatEl%alphaMass = alphaMass ; InMatEl%andim = andim

! Initialize everything to zero

    InMatEl%softExpT = 0; InMatEl%softExpC = 0;  InMatEl%HardExp = 0; InMatEl%jetExp = 0

! Non-Global 2-loop soft function

    InMatEl%NGlist = alphaSList(2) * NGLOas2(nf)

! Assign non-log matrix element terms

    InMatEl%softExpT(0,:) = GenMat('SoftThrust', nf, s3T)
    InMatEl%softExpC(0,:) = GenMat('SoftCparam', nf, s3C)
    InMatEl%jetExp(0,:)   = GenMat('jet' , nf, j3)
    InMatEl%HardExp(0,:)  = GenMat('hard', nf, 0._dp)

! Compute full matrix elements:
  ! Expand alphaS

    call andim%expandAlpha(InMatEl%softExpT)
    call andim%expandAlpha(InMatEl%softExpC)
    call andim%expandAlpha(InMatEl%jetExp)
    call andim%expandAlpha(InMatEl%HardExp)

    if ( present(Q) .and. present(muH) ) call InMatEl%SetHard(Q, muH)

  ! Add alphaQCD expansion + anomalous dimension expansion

    InMatEl%softExpT = InMatEl%softExpT - 2 * coefCusp + coefSoft
    InMatEl%softExpC = InMatEl%softExpC - 2 * coefCusp + coefSoft
    InMatEl%jetExp   = InMatEl%jetExp   + 2 * coefCusp + coefJet/2

    betaList = andim%betaQCD('betaList')

    InMatEl%ThrustGammaR = ExpEuler * andim%GammaRComputer( InMatEl%softExpT(1,:) )
    InMatEl%CparamGammaR = ExpEuler * andim%GammaRComputer( InMatEl%softExpC(1,:) )

    InMatEl%sCoefT =  andim%sCoef( InMatEl%ThrustGammaR * betaList(1:3) )
    InMatEl%sCoefC =  andim%sCoef( InMatEl%CparamGammaR * betaList(1:3) )

    TwoList = powList(2,4)

 ! rescale jet logs because of their mass dimension

     log_rescale: do i = 1, 4
       InMatEl%jetExp(i,:) = InMatEl%jetExp(i,:)/TwoList(i)
     end do log_rescale

 ! compute Soft subtractions

    if ( present(muS) .and. present(R) ) call InMatEl%SetDelta(muS, R, alphaSList)

 ! multiply orders by power of alphaQCD

  if ( present(muS) ) call AddAlpha(InMatEl%softExpT, alphaSList)
  if ( present(muS) ) call AddAlpha(InMatEl%softExpC, alphaSList)
  if ( present(muJ) ) call AddAlpha(InMatEl%jetExp  , alphaJList)

 ! compute logs

    if ( present(muS) .and. present(muJ) ) call InMatEl%AddLogs( InMatEl%jetExp, log(muJ**2/Q/muS) )

   end function InMatEl

!ccccccccccccccc

  type (MatrixElementsMass) function InMatElMass(AlphaAll, nf, runMass, s3T, &
  s3C, j3, b3, Q, muH, muJ, muM, muS, R, Rmass, muLambda, muLambdaNl)
    type (Alpha)       , intent(in) :: AlphaAll
    integer            , intent(in) :: nf, runMass
    real (dp)          , intent(in) :: s3T, s3C, b3, j3, muLambda, muLambdaNl
    real (dp), optional, intent(in) :: Q, muH, muJ, muS, R, Rmass, muM
    character (len = 5)             :: str
    integer                         :: i, nl, nS
    type (Running)                  :: alphaMass, alphaMassNl
    real (dp), dimension(0:4)       :: betaS, betaList, betaListNl
    real (dp), dimension(3)         :: lgRList
    real (dp), dimension(0:4)       :: ListFact
    real (dp), dimension(3)         :: alphaSList, a, b
    real (dp), dimension(4)         :: c, alphaJList, lgRmassList, lgMSLowList, lgMSmassList
    integer  , dimension(3)         :: TwoList
    type (AnomDim)                  :: andim, andimNl, andimH, andimS
    real (dp)                       :: alphaJ, alphaS, mm, alphaR, EuR
    real (dp), dimension(0:3,0:3)   :: tab
    real (dp), dimension(0:4,4)     :: coefMSR, coefMassLow, coefMSRNatural, coefMass
    real (dp), dimension(0:4,3)     :: coefCusp, coefSoft, coefJet, coefHard, &
    coefHm, coefCuspnl, coefCuspS, coefSoftNl, coefBjet

! Initialize everything to zero

    InMatElMass%softExpT = 0 ;  InMatElMass%softExpC = 0  ;  InMatElMass%softExpCNl = 0
    InMatElMass%HardExp  = 0 ;  InMatElMass%softExpTNl = 0;  coefCuspNl = 0; alphaS = 0

    alphaMass   = Running(nf, runMass, AlphaAll, muLambda)  ; nl = nf - 1;  alphaJ = 0
    alphaMassNl = Running(nl, runMass, AlphaAll, muLambdaNl); coefCusp = 0; alphaR = 0
    InMatElMass%alphaMassNl = alphaMassNl;  InMatElMass%alphaMass = alphaMass

    andim   = alphaMass%adim()   ; mm  = alphaMass%scales('mL') ; InMatElMass%nf      = nf
    andimNl = alphaMassNl%adim() ; str = andim%scheme()         ; InMatElMass%BJetExp = 0

    if ( present(R) ) then
      InMatElMass%R = R; euR = ExpEuler * R
    end if

    if ( present(Rmass) ) InMatElMass%Rmass = Rmass

    InMatElMass%mm = mm; InMatElMass%AlphaAll = AlphaAll;  InMatElMass%HmExp = 0
    ListFact = factList(4)     ; TwoList = powList(2,4) ;  InMatElMass%ListFact = ListFact
    InMatElMass%alphaMMlist = 0; InMatElMass%muJ = 0    ;  InMatElMass%jetExp = 0

  ! Expand Cusp anomalous dimension for nf and nl flavors

    coefHard = 0;  betaListNl = andimNl%betaQCD('betaList')

    call   andim%kTildeExpand(   andim%betaQCD('cusp'), coefCusp   )
    call andimNl%kTildeExpand( andimNl%betaQCD('cusp'), coefCuspNl )
    call   andim%wTildeExpand(   andim%betaQCD('hard'), coefHard   )

    InMatElMass%coefHard = coefHard;  InMatElMass%coefCusp = coefCusp
    InMatElMass%HardExp(0,:)  = GenMat('hard', nF, 0._dp)
    call andimH%expandAlpha(InMatElMass%HardExp)

    if ( present(Q) .and. present(muH) ) call InMatElMass%SetHard(Q, muH)
    if ( present(muM) ) InMatElMass%muM = muM

    ! Jet Sector

    if ( .not. present(muJ) .or. muJ < muM) then
      InMatElMass%alphaMMlist = powList(alphaMass%alphaQCD(mm)/Pi,3)
      coefBJet = 0; call andimNl%wTildeExpand( andimNl%betaQCD('bJet'), coefBJet )
      InMatElMass%BJetExp(0,:) = GenMat('bJet', nl, b3)
      coefHm = 0; call andim%kTildeExpand( andim%betaQCD('cusp') , coefHm )

      call andimNl%expandAlpha(InMatElMass%BJetExp)
      call InMatElMass%SetHardMass(muM)

      InMatElMass%BJetExp   = InMatElMass%BJetExp + coefCuspNl + coefBJet/2
      InMatElMass%JetGammaR = ExpEuler * andimNl%GammaRComputer( InMatElMass%BJetExp(1,:) )/2
      InMatElMass%sCoefbJet = andimNl%sCoef( andimNl%GammaRComputer( InMatElMass%BJetExp(1,:) ) * betaListNl(1:) ) * ExpEuler/2
      InMatElMass%jetMatching(0) = 1

      a = ExpEuler * DeltaComputer(InMatElMass%BJetExp(1:,:), [0,0,0] * 0._dp, 1)/2
      c = MSbarDelta(nf - 1, 1); tab = andim%alphaMatching(nf)
      b = tab(:2,0); call alphaReExpand(a,b);  a = c(:3) - a
      a = a * InMatElMass%alphaMMlist

      InMatElMass%jetMatching(1:) = a

    end if

    if ( .not. present(muJ) .or. muJ >= muM) then
      call  andim%wTildeExpand( andim%betaQCD('jet'), coefJet  )
      InMatElMass%jetExp(0,:) = GenMat('jet', nf, j3)
      call andim%expandAlpha(InMatElMass%jetExp)
      InMatElMass%jetExp = InMatElMass%jetExp + 2 * coefCusp + coefJet/2

 ! rescale jet logs because of their mass dimension

      log_rescale: do i = 1, 4
        InMatElMass%jetExp(i,:) = InMatElMass%jetExp(i,:)/TwoList(i)
      end do log_rescale
    end if

    if ( present(muJ) ) then

      InMatElMass%muJ = muJ

      if (muJ < muM) then

        alphaJ  = alphaAll%alphaQCD(nl, muJ)/Pi; alphaJList = powList(alphaJ, 4)
        InMatElMass%MSLow = alphaMassNl%MSbarMassLow(muJ); InMatElMass%alphaB = alphaJ
        coefMassLow = InMatElMass%MSLow * MSbarDeltaPiece(nl, 0)

        if ( present(Rmass) ) then
          coefMSR = 0       ; coefMSR(0,:)        = Rmass * andimNl%MSRDelta()
          coefMSRNatural = 0; coefMSRNatural(0,:) = Rmass * MSbarDelta(nl, 0)
          call andimNl%expandAlpha(coefMSR); call andimNl%expandAlpha(coefMSRNatural)
          lgRMassList = powList( log(muJ/Rmass), 4 )
          call AddAlpha(coefMSR       , alphaJList)
          call AddAlpha(coefMSRNatural, alphaJList)
        end if

        call AddAlpha( InMatElMass%BJetExp, alphaJList(:3) )
        call AddAlpha(coefMassLow         , alphaJList     )

        lgMSLowList            = powList( log(muJ/InMatElMass%MSLow), 4 )
        InMatElMass%deltaMSLow = DeltaComputer(coefMassLow, lgMSLowList, 0)

        if ( present(Rmass) ) then
          InMatElMass%deltamJet = ExpEuler * Rmass * &
          DeltaComputer( InMatElMass%BJetExp(1:,:), lgRMassList, 1 )/2

          InMatElMass%deltaMSR        = DeltaComputer(coefMSR       , lgRMassList, 0)
          InMatElMass%deltaMSRNatural = DeltaComputer(coefMSRNatural, lgRMassList, 0)
        end if
      else

        InMatElMass%MSmass = alphaMass%MSbarMass(muJ); coefJet = 0
        lgMSmassList = powList( log(muJ/InMatElMass%MSmass), 4)
        alphaJ = alphaMass%alphaQCD(muJ)/Pi;  alphaJList = powList(alphaJ, 4)
        coefMass = InMatElMass%MSmass * MSbarDeltaPiece(nl, 1)

        call AddAlpha( InMatElMass%jetExp, alphaJList(:3) )
        call AddAlpha( coefMass          , alphaJList     )

        InMatElMass%deltaMS = coefMass(0,:) + matmul( lgMSmassList, coefMass(1:,:) )

        call InMatElMass%AddLogs(InMatElMass%jetExp, log(muJ**2/Q/muS) )

      end if
    end if

! compute Logs

    if ( present(muS) .and. present(R) ) lgRList = powList( log(muS/R), 3 )

! Assign all constants

    coefSoft = 0;  coefSoftNl = 0

    if ( present(muS) ) then

      InMatElMass%muS = muS

      if (muS >= muM) then
        if ( present(R) ) alphaR = alphaMass%alphaQCD(R)/Pi
        alphaS = alphaMass%alphaQCD(muS)/Pi ; InMatElMass%alphaMassS = alphaMass
        andimS = andim; coefCuspS = coefCusp; nS = nf

      else
        if ( present(R) ) alphaR = alphaAll%alphaQCD(nl,R)/Pi
        alphaS = alphaAll%alphaQCD(nl,muS)/Pi    ; InMatElMass%alphaMassS = alphaMassNl
        andimS = andimNl; coefCuspS = coefCuspNl ; nS     = nl
      end if

      InMatElMass%andimS = andimS          ; betaS      = andimS%betaQCD('beta')
      InMatElMass%beta   = betaS           ; alphaSList = powList(alphaS, 3)

    else
      alphaSList = 1; InMatElMass%beta = andim%betaQCD('beta'); andimS = andim; nS = nF
    end if

    InMatElMass%betaNl = andimNl%betaQCD('beta')

  ! Expand anomalous dimension

    call  andimS%wTildeExpand( andimS%betaQCD('soft'), coefSoft )

    InMatElMass%alphaList(2) = alphaJ; InMatElMass%alphaList(5:) = [alphaS, alphaR]

    InMatElMass%andim = andim;  InMatElMass%andimNl   = andimNl

! Non-Global 2-loop soft function

    InMatElMass%NGlist = alphaSList(2) * NGLOas2(nS)

! Assign non-log matrix element terms

    InMatElMass%softExpT(0,:) = GenMat('SoftThrust', nS, s3T)
    InMatElMass%softExpC(0,:) = GenMat('SoftCparam', nS, s3C)

    if (  muS >= muM .or. ( .not. present(muS) )  ) then
      InMatElMass%softExpTNl(0,:) = GenMat('SoftThrust', nl, s3T)
      InMatElMass%softExpCNl(0,:) = GenMat('SoftCparam', nl, s3C)
      call andimNl%expandAlpha(InMatElMass%softExpTNl)
      call andimNl%expandAlpha(InMatElMass%softExpCNl)
    end if

! Compute full matrix elements:
  ! Expand alphaS

    call andimS%expandAlpha(InMatElMass%softExpT)
    call andimS%expandAlpha(InMatElMass%softExpC)

  ! Add alphaQCD expansion + anomalous dimension expansion

    InMatElMass%softExpT = InMatElMass%softExpT - 2 * coefCuspS  + coefSoft
    InMatElMass%softExpC = InMatElMass%softExpC - 2 * coefCuspS  + coefSoft

    betaList = andimS%betaQCD('betaList')

    InMatElMass%ThrustGammaR = ExpEuler * andimS%GammaRComputer( InMatElMass%softExpT(1,:) )
    InMatElMass%CparamGammaR = ExpEuler * andimS%GammaRComputer( InMatElMass%softExpC(1,:) )

    InMatElMass%sCoefT =  andimS%sCoef( InMatElMass%ThrustGammaR * betaList(1:3) )
    InMatElMass%sCoefC =  andimS%sCoef( InMatElMass%CparamGammaR * betaList(1:3) )

    if (  muS >= muM .or. (.not. present(muS) )  ) then

      call andimNl%wTildeExpand(andimNl%betaQCD('soft'), coefSoftNl  )

      InMatElMass%softExpTNl = InMatElMass%softExpTNl - 2 * coefCuspNl + coefSoftNl
      InMatElMass%softExpCNl = InMatElMass%softExpCNl - 2 * coefCuspNl + coefSoftNl

      InMatElMass%ThrustGammaRNl = ExpEuler * andimNl%GammaRComputer( InMatElMass%softExpTNl(1,:) )
      InMatElMass%CparamGammaRNl = ExpEuler * andimNl%GammaRComputer( InMatElMass%softExpCNl(1,:) )

      InMatElMass%sCoefTNl = andimNl%sCoef( InMatElMass%ThrustGammaRNl * betaListNl(1:3) )
      InMatElMass%sCoefCNl = andimNl%sCoef( InMatElMass%CparamGammaRNl * betaListNl(1:3) )

    else if ( muS < muM .and. present(muS) ) then

      InMatElMass%softExpTNl     = InMatElMass%softExpT
      InMatElMass%softExpCNl     = InMatElMass%softExpC
      InMatElMass%ThrustGammaRNl = InMatElMass%ThrustGammaR
      InMatElMass%CparamGammaRNl = InMatElMass%CparamGammaR
      InMatElMass%sCoefTNl       = InMatElMass%sCoefT
      InMatElMass%sCoefCNl       = InMatElMass%sCoefC

    end if

 ! compute Soft subtraction

   if ( present(muS) .and. present(R) ) call InMatElMass%SetDelta(muS, R, alphaSList)

 ! multiply orders by power of alphaQCD

   call AddAlpha(InMatElMass%softExpT, alphaSList)
   call AddAlpha(InMatElMass%softExpC, alphaSList)

   end function InMatElMass

!ccccccccccccccc

  subroutine SetAlpha(self, alpha)
    class (MatricesElements), intent(inout) :: self
    real (dp)               , intent(in)    :: alpha

    call self%AlphaAll%SetAlpha(alpha); call self%AlphaMass%SetAlpha(alpha)

    select type (self)
    class is (MatricesElementsMass)
      call self%alphaMassNl%SetAlpha(alpha)
    type is (MatrixElementsMass)
      call self%alphaMassS%SetAlpha(alpha)
    end select

  end subroutine SetAlpha

!ccccccccccccccc

  subroutine SetMass(self, m, mu)
    class (MatricesElementsMass), intent(inout) :: self
    real (dp)                   , intent(in)    :: m, mu

    if (self%nf == 6) call self%SetMTop   (m, mu)
    if (self%nf == 5) call self%SetMBottom(m, mu)
    if (self%nf == 4) call self%SetMCharm (m, mu)

  end subroutine SetMass

!ccccccccccccccc

  subroutine SetMtop(self, mT, muT)
    class (MatricesElements), intent(inout) :: self
    real (dp)               , intent(in)    :: mT, muT

    call self%AlphaAll%SetMtop(mT, muT); call self%AlphaMass%SetMtop(mT, muT)

    select type (self)
    class is (MatricesElementsMass)
      call self%alphaMassNl%SetMtop(mT, muT)

      if (self%nf == 6) then
        self%mm = mT; self%JetMatching(1:) = self%JetMatching(1:)/self%alphaMMlist
        self%alphaMMlist = powList(self%alphaMass%alphaQCD(self%mm)/Pi,3)
        self%JetMatching(1:) = self%JetMatching(1:) * self%alphaMMlist
      end if

    type is (MatrixElementsMass)
      call self%alphaMassS%SetMtop(mT, muT)
    end select

  end subroutine SetMtop

!ccccccccccccccc

  subroutine SetMBottom(self, mB, muB)
    class (MatricesElements), intent(inout) :: self
    real (dp)               , intent(in)    :: mB, muB

    call self%AlphaAll%SetMBottom(mB, muB); call self%AlphaMass%SetMBottom(mB, muB)

    select type (self)
    class is (MatricesElementsMass)
      call self%alphaMassNl%SetMBottom(mB, muB)

      if (self%nf == 5) then
        self%mm = mB; self%JetMatching(1:) = self%JetMatching(1:)/self%alphaMMlist
        self%alphaMMlist = powList(self%alphaMass%alphaQCD(self%mm)/Pi,3)
        self%JetMatching(1:) = self%JetMatching(1:) * self%alphaMMlist
      end if

    type is (MatrixElementsMass)
      call self%alphaMassS%SetMBottom(mB, muB)
    end select

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, mC, muC)
    class (MatricesElements), intent(inout) :: self
    real (dp)               , intent(in)    :: mC, muC

    call self%AlphaAll%SetMCharm(mC, muC); call self%AlphaMass%SetMCharm(mC, muC)

    select type (self)
    class is (MatricesElementsMass)
      call self%alphaMassNl%SetMCharm(mC, muC)

      if (self%nf == 4) then
        self%mm = mC; self%JetMatching(1:) = self%JetMatching(1:)/self%alphaMMlist
        self%alphaMMlist = powList(self%alphaMass%alphaQCD(self%mm)/Pi,3)
        self%JetMatching(1:) = self%JetMatching(1:) * self%alphaMMlist
      end if

    type is (MatrixElementsMass)
      call self%alphaMassS%SetMCharm(mC, muC)
    end select

  end subroutine SetMCharm

!ccccccccccccccc

  subroutine SetHard(self, Q, muH)
    class (MatricesElements), intent(inout) :: self
    real (dp)               , intent(in)    :: Q, muH
    complex (dp)                            :: alphaHComplex
    real (dp)                               :: alphaH, alphaIm
    real (dp), dimension(3)                 :: alphaHList

    self%Q = Q; self%muH = muH

    if (muH > 0) then
      alphaH = self%alphaMass%alphaQCD(muH)/Pi
    else
      alphaHComplex = self%alphaMass%alphaQCD( dcmplx(0._dp, muH) )/Pi
      alphaH = real(alphaHComplex); alphaIm = ( imag(alphaHComplex)/alphaH )**2
    end if

    alphaHList = powList(alphaH, 3)

    self%HardExp(0,:) = 2 * self%SetHardFunc( log( Q/abs(muH) ), alphaIm, self%coefCusp, &
                                              self%coefHard)
    self%HardExp(1:,:) = 0

    call AddAlpha(self%HardExp , alphaHList)

    select type (self)
    class is (MatrixElements)
      self%alphaList(1) = alphaH
    class is (MatricesElementsMass)
      self%alphaList(1) = alphaH
    end select

  end subroutine SetHard

!ccccccccccccccc

  subroutine SetHardMass(self, muM)
    class (MatricesElementsMass), intent(inout) :: self
    real (dp)                   , intent(in)    :: muM
    real (dp)                                   :: LQ, Lm

    self%muM = muM; LQ = 2 * log(self%Q/self%mm); Lm = 2 * log(self%mm/muM); self%HmExp = 0
    self%alphaList(3:4) = [self%alphaMass%alphaQCD(muM), self%alphaMassNl%alphaQCD(muM)]/Pi

    self%HmExp(0,1) = self%alphaList(3) * ( 3.7632893778988175_dp + 2 * Lm * (Lm - 1)/3 )
    self%HmExp(0,2) = 22.02886430312693_dp - 28 * LQ/81 - &
    Lm**2 * (- 0.5473424096795991_dp + LQ/9 + 13._dp * self%nf/54) + &
    Lm**3/54 * (2 * self%nf - 37) - 2.4985480650201737_dp * self%nf - &
    Lm * (9.51716816196516_dp + 10 * LQ/27 - 1.2063904494634092_dp * self%nf)

    if ( self%scheme == 'MSbar') self%HmExp(0,2) = self%HmExp(0,2) + 16 * (2 * Lm - 1)/9

    self%HmExp(0,2) = self%alphaList(3)**2 * self%HmExp(0,2)

  end subroutine SetHardMass

!ccccccccccccccc

  type (AnomDim) function adim(self, str)
    class (MatricesElements)     , intent(in) :: self
    character (len = *), optional, intent(in) :: str

    select type (self)
     class is (MatricesElements)
       adim = self%andim
     class is (MatricesElementsMass)
       if ( present(str) .and. str(:2) == 'nl' ) then
         adim = self%andimNl
       else
         adim = self%andim
       end if
    end select
  end function

!ccccccccccccccc

  type (Running) function run(self)
    class (MatricesElements), intent(in) :: self
    run = self%alphaMass
  end function

!ccccccccccccccc

  type (Running) function runNl(self)
    class (MatricesElementsMass), intent(in) :: self
    runNl = self%alphaMassNl
  end function

!ccccccccccccccc

  subroutine SetDelta(self, muS, R, alphaList)
    class (MatricesElements), intent(inout) :: self
    real (dp)               , intent(in   ) :: muS, R
    real (dp), dimension(3) , intent(in   ) :: alphaList

    real (dp)                               :: a, Pre, euR, alphaR, beta0
    real (dp), dimension(3)                 :: lgRList

    lgRList = powList( log(muS/R), 3 ) ; euR = ExpEuler * R ; self%R = R ; self%muS = muS

    select type (self)
    class is (MatricesElements)

      alphaR      = self%alphaMass%alphaQCD(R); beta0 = self%beta(0)
      self%deltaT = euR * DeltaComputer( self%softExpT(1:,:), lgRList, 1 ) * alphaList
      self%deltaC = euR * DeltaComputer( self%softExpC(1:,:), lgRList, 1 ) * alphaList

    class is (MatricesElementsMass)

      if (muS >= self%muM) then
        alphaR      = self%alphaMass%alphaQCD(R); beta0 = self%beta(0)
        self%deltaT = euR * DeltaComputer( self%softExpT(1:,:), lgRList, 1 ) * alphaList
        self%deltaC = euR * DeltaComputer( self%softExpC(1:,:), lgRList, 1 ) * alphaList
      else
        alphaR      = self%alphaMassNl%alphaQCD(R); beta0 = self%betaNl(0)
        self%deltaT = euR * DeltaComputer( self%softExpTNl(1:,:), lgRList, 1 ) * alphaList
        self%deltaC = euR * DeltaComputer( self%softExpCNl(1:,:), lgRList, 1 ) * alphaList
      end if

    end select

    a = 6/beta0 * log( Pi * alphaList(1)/alphaR )
    Pre = sqrt(pi)/2 * gamma(1 + a)/gamma(1.5_dp + a)

    self%deltaTHadron = Pre * self%deltaT;  self%deltaCHadron = Pre * self%deltaC

  end subroutine SetDelta

!ccccccccccccccc

  real (dp) function alphaScale(self,str)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str

    alphaScale = 0

      if ( str(:4) == 'hard'   ) alphaScale = self%alphaList(1)
      if ( str(:3) == 'jet'    ) alphaScale = self%alphaList(2)
      if ( str(:6) == 'massNf' ) alphaScale = self%alphaList(3)
      if ( str(:6) == 'massNl' ) alphaScale = self%alphaList(4)
      if ( str(:4) == 'soft'   ) alphaScale = self%alphaList(5)
      if ( str(:1) == 'R'      ) alphaScale = self%alphaList(6)

  end function alphaScale

!ccccccccccccccc

  pure function expandExpMat(a) result(b)
    real (dp), dimension(0:4,3), intent(in) :: a
    real (dp), dimension(0:6,0:3)           :: b

    b        = 0;  b(4,2) = a(2,1)**2/2; b(1,1:2) = [ a(1,1), a(0,1) * a(1,1) + a(1,2) ]
    b(0,:3)  = [1._dp, a(0,1), a(0,1)**2/2 + a(0,2), a(0,1)**3/6 + a(0,1) * a(0,2) + a(0,3)]
    b(1,3)   = a(0,1)**2 * a(1,1)/2 + a(0,2) * a(1,1) + a(0,1) * a(1,2) + a(1,3)
    b(2,1:2) = [ a(2,1), a(1,1)**2/2 + a(0,1) * a(2,1) + a(2,2) ]
    b(2,3)   = a(0,1) * a(1,1)**2/2 + a(1,1) * a(1,2) + a(0,1)**2 * a(2,1)/2 +  &
    a(0,2) * a(2,1) + a(0,1) * a(2,2) + a(2,3)

    b(3,2)   = a(1,1) * a(2,1) + a(3,2)
    b(3,3)   = a(1,1)**3/6 + a(0,1) * a(1,1) * a(2,1) + a(1,2) * a(2,1) + &
    a(1,1) * a(2,2) + a(0,1) * a(3,2) + a(3,3)

    b(4,3)   = a(1,1)**2 * a(2,1)/2 + a(0,1) * a(2,1)**2/2 + a(2,1) * a(2,2) + &
    a(1,1) * a(3,2) + a(4,3)

    b(5,3)   = a(1,1) * a(2,1)**2/2 + a(2,1) * a(3,2);  b(6,3) = a(2,1)**3/6

  end function expandExpMat

!ccccccccccccccc

  pure function expandExpVec(a) result(b)
    real (dp), dimension(3), intent(in) :: a
    real (dp), dimension(0:3)           :: b

    b = [ 1._dp, a(1), a(1)**2/2 + a(2), a(1)**3/6 + a(1) * a(2) + a(3) ]

  end function expandExpVec

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  pure complex (dp) function expandExpComplex(a)
    complex (dp), dimension(:), intent(in) :: a
    integer                                :: i, order

    order = size(a); expandExpComplex = 0

    do i = 0, min(order,3)
      expandExpComplex = expandExpComplex + expandExpComplexOrder(i, a)
    end do

 end function expandExpComplex

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  pure complex (dp) function expandExpReal(a)
    real (dp), dimension(:), intent(in) :: a
    integer                             :: i, order

    order = size(a); expandExpReal = 0

    do i = 0, min(order,3)
      expandExpReal = expandExpReal + expandExpRealOrder(i, a)
    end do

 end function expandExpReal

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  pure complex (dp) function expandExpComplexOrder(order, a)
    integer                       , intent(in) :: order
    complex (dp), dimension(order), intent(in) :: a

    expandExpComplexOrder = 0
    if (order == 0) expandExpComplexOrder = (1._dp, 0._dp)
    if (order == 1) expandExpComplexOrder = a(1)
    if (order == 2) expandExpComplexOrder = a(1)**2/2 + a(2)
    if (order == 3) expandExpComplexOrder = a(1)**3/6 + a(1) * a(2) + a(3)

  end function expandExpComplexOrder

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  pure real (dp) function expandExpRealOrder(order, a)
    integer                    , intent(in) :: order
    real (dp), dimension(order), intent(in) :: a

    expandExpRealOrder = 0
    if (order == 0) expandExpRealOrder = 1
    if (order == 1) expandExpRealOrder = a(1)
    if (order == 2) expandExpRealOrder = a(1)**2/2 + a(2)
    if (order == 3) expandExpRealOrder = a(1)**3/6 + a(1) * a(2) + a(3)

  end function expandExpRealOrder

!ccccccccccccccc

 pure function expandProd(a, b) result(c)
   real (dp), dimension(:), intent(in)               :: a, b
   real (dp), dimension(  min( size(a), size(b) )  ) :: c
   integer                                           :: i, j

   c = 0

   do i = 1, min( size(a), size(b) )
     do j = 1, i
       c(i) = c(i) + a(1 + i - j) * b(j)
     end do
   end do

  end function expandProd

!ccccccccccccccc

  function CoefMat(self, str) result(bet)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str
    real (dp)         , dimension(0:4,3) :: bet

    bet = 0

    if ( str(:4)  == 'hard'  ) bet = self%HardExp
    if ( str(:3)  == 'jet'   ) bet = self%JetExp
    if ( str(:6) == 'thrust' ) bet = self%softExpT
    if ( str(:6) == 'Cparam' ) bet = self%softExpC

    select type (self)
     class is (MatricesElementsMass)
       if ( str(:4) == 'bJet'     ) bet = self%BjetExp
       if ( str(:8) == 'Nlthrust' ) bet = self%softExpTNl
       if ( str(:8) == 'NlCparam' ) bet = self%softExpCNl
       if ( str(:2) == 'Hm'       ) bet = self%HmExp
    end select

  end function CoefMat

!ccccccccccccccc

  function delta(self, str) result(bet)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str
    real (dp)             , dimension(4) :: bet

    bet = 0

    if ( str(:4) == 'hard' ) bet(:3) = self%HardExp(0,:)

    if ( str(:12) == 'thrustHadron' .or. str(:9) == 'SJMHadron' ) then
      bet(:3) = self%deltaTHadron
    else if ( str(:6) == 'thrust' .or. str(:3) == 'SJM')   then
      bet(:3) = self%deltaT
    end if

    if ( str(:12) == 'CparamHadron' ) then
      bet(:3) = self%deltaCHadron
    else if ( str(:6) == 'Cparam' )   then
      bet(:3) = self%deltaC
    end if

    if ( str(:9) == 'HJMHadron' )     then
      bet(:3) = self%deltaTHadron/2
    else if ( str(:3) == 'HJM' )      then
      bet(:3) = self%deltaT/2
    end if

    select type (self)
    class is (MatrixElementsMass)
       if ( str(:7 ) == 'JetMass'    ) bet(:3) = self%deltamJet
       if ( str(:5 ) == 'MSbar'      ) bet = self%deltaMS
       if ( str(:3 ) == 'MSR'        ) bet = self%deltaMSR
       if ( str(:10) == 'MSRNatural' ) bet = self%deltaMSRNatural
       if ( str(:8 ) == 'MSbarLow'   ) bet = self%deltaMSLow
       if ( str(:8 ) == 'hardMass'   ) bet(:3) = self%HmExp(0,:)
    end select

  end function delta

!ccccccccccccccc

  function NonGlobalList(self) result(bet)
    class (MatricesElements), intent(in) :: self
    real (dp)             , dimension(6) :: bet
    bet = self%NGlist
  end function NonGlobalList

!ccccccccccccccc

  function NonGlobalListNl(self) result(bet)
    class (MatricesElementsMass), intent(in) :: self
    real (dp)                 , dimension(6) :: bet
    bet = self%NGlistNl
  end function NonGlobalListNl

!ccccccccccccccc

  function sCoef(self, str) result(bet)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str
    real (dp)             , dimension(3) :: bet
    real (dp)           , dimension(0:4) :: aux

    bet = 0

    if ( str(:6) == 'thrust' .or. str(:3) == 'SJM' ) bet = self%sCoefT
    if ( str(:3) == 'HJM'    ) bet = self%sCoefT/2
    if ( str(:6) == 'Cparam' ) bet = self%sCoefC

    select type (self)
     class is (MatricesElementsMass)
       if ( str(:4) == 'bJet'   ) bet = self%sCoefBJet

       if ( str(:3) == 'MSR' .and. str(:10) == 'MSRNatural' ) then
         aux = self%andim%betaQCD('MSRsCoef')       ;  bet = aux(1:3)
       else if ( str(:3) == 'MSR' ) then
         aux = self%andim%betaQCD('sCoefMSRNatural');  bet = aux(1:3)
       end if

    end select

  end function sCoef

!ccccccccccccccc

  function sCoefNl(self, str) result(bet)
    class (MatricesElementsMass), intent(in) :: self
    character (len = *)         , intent(in) :: str
    real (dp)                 , dimension(3) :: bet

    bet = 0

    if ( str(:6) == 'thrust' .or. str(:3) == 'SJM' ) bet = self%sCoefTNl
    if ( str(:3) == 'HJM'    ) bet = self%sCoefTNl/2
    if ( str(:6) == 'Cparam' ) bet = self%sCoefCNl
    if ( str(:4) == 'bJet'   ) bet = self%sCoefBJet

  end function sCoefNl

!ccccccccccccccc

  function GammaR(self, str) result(bet)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str
    real (dp)             , dimension(3) :: bet

    bet = 0

    if ( str(:6) == 'thrust' .or. str(:3) == 'SJM' ) bet = self%ThrustGammaR
    if ( str(:3) == 'HJM'     ) bet = self%ThrustGammaR/2
    if ( str(:6) == 'Cparam'  ) bet = self%CparamGammaR

    select type (self)
     class is (MatricesElementsMass)
       if ( str(:7) == 'JetMass' ) bet = self%JetGammaR
    end select
  end function GammaR

!ccccccccccccccc

  real (dp) function JetMass(self, order, run, R, mu)
    class (MatricesElementsMass), intent(in) :: self
    integer                   , intent(in) :: run, order
    real (dp)                 , intent(in) :: R, mu

    JetMass = self%mm * sum(  self%jetMatching( :min(order,3) )  ) + &
    self%DiffDeltaGapNl('bJet', run, self%mm, R, self%mm, mu)

  end function JetMass

!ccccccccccccccc

  subroutine setGammaShift(self, order, run, gap, shape, delta0, R0, mu0, h, shift, delta)
    class (MatricesElements) , intent(in)                           :: self
    real (dp)                , intent(in)                           :: h, delta0, R0, mu0
    character (len = *)      , intent(in)                           :: gap, shape
    integer                  , intent(in)                           :: order, run
    real (dp)                , intent(out), dimension(order, order) :: delta
    real (dp)                , intent(out)                          :: shift
    real (dp),                              dimension(3)            :: del

    del = 0;  delta = 0;  shift = 0

    if ( gap(:6) == 'thrust' .or. gap(:6) == 'Cparam' .or. gap(2:3) == 'JM' ) then

      del = self%delta(gap);  delta(1,:) = del(:order)

      shift = self%DiffDeltaGap(gap, run, R0, self%R, mu0, self%muS) &
      + gapCons(gap) * delta0 + (1 - h) * Sum( delta(1,:) )

      delta = h * delta

      if (  ( shape(:6) == 'thrust' .or. shape(2:3) == 'JM' ) .and. gap(:6) == 'Cparam'  ) then
        shift = 4 * shift/Pi;  delta = 4 * delta/Pi
      else if (  ( gap(:6) == 'thrust' .or. gap(2:3) == 'JM' ) .and. shape(:6) == 'Cparam'  ) then
        shift = Pi * shift/4;  delta = Pi * delta/4
      end if
    end if

    if (order > 1) delta(2 ,2) = delta(1,1)**2/2
    if (order > 2) delta(2:,3) = [ delta(1,1) * delta(1,2), delta(1,1)**3/6 ]

  end subroutine setGammaShift

!ccccccccccccccc

  real (dp) function DiffDeltaGap(self, str, order, r0, r1, mu0, mu1)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str
    integer                 , intent(in) :: order
    real (dp)               , intent(in) :: r0, r1, mu0, mu1

    if ( str(7:12) == 'Hadron' ) then

      DiffDeltaGap = self%DiffDeltaGapHadron( str(:6), order, r0, r1, mu0, mu1 )

    else

      select type (self)
      class is (MatricesElements)
        DiffDeltaGap = self%alphaMass%DiffDelta( self%sCoef(str), &
        cuspConst(str), order, r0, r1, mu0, mu1 )

       class is (MatrixElementsMass)

        DiffDeltaGap = self%alphaMassS%DiffDelta( self%sCoef(str), &
        cuspConst(str), order, r0, r1, mu0, mu1 )

      end select

    end if

  end function DiffDeltaGap

!ccccccccccccccc

  real (dp) function DiffDeltaGapNl(self, str, order, r0, r1, mu0, mu1)
    class (MatricesElementsMass), intent(in) :: self
    character (len = *)         , intent(in) :: str
    integer                     , intent(in) :: order
    real (dp)                   , intent(in) :: r0, r1, mu0, mu1

    DiffDeltaGapNl = self%alphaMassNl%DiffDelta( self%sCoefNl(str), &
    cuspConst(str), order, r0, r1, mu0, mu1 )

  end function DiffDeltaGapNl

!ccccccccccccccc

  real (dp) function DiffDeltaGapHadron(self, str, order, r0, r1, mu0, mu1)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str
    integer                 , intent(in) :: order
    real (dp)               , intent(in) :: r0, r1, mu0, mu1

    select type (self)
     class is (MatricesElements)

      DiffDeltaGapHadron = self%alphaMass%DiffDeltaHadron(   &
      self%gammaR( str(:6) ), order,r0, r1, mu0, mu1  )

    class is (MatrixElementsMass)

      DiffDeltaGapHadron = self%alphaMassS%DiffDeltaHadron(   &
      self%gammaR( str(:6) ), order,r0, r1, mu0, mu1  )

    end select

  end function DiffDeltaGapHadron

!ccccccccccccccc

  subroutine SetScales(self, str, mu)
    class (MatricesElements), intent(inout) :: self
    character (len = *)     , intent(in)    :: str
    real (dp)               , intent(in)    :: mu

    if ( str(:3) == 'muH' ) self%muH = mu; if ( str(:1) == 'Q' ) self%Q = mu

    select type (self)
    class is (MatricesElements)
      if ( str(:3) == 'muJ'   ) self%muJ   = mu
      if ( str(:3) == 'muS'   ) self%muS   = mu
      if ( str(:1) == 'R'     ) self%R     = mu
    class is (MatricesElementsMass)
      if ( str(:3) == 'muJ'   ) self%muJ   = mu
      if ( str(:3) == 'muS'   ) self%muS   = mu
      if ( str(:1) == 'R'     ) self%R     = mu
      if ( str(:3) == 'muM'   ) self%muM   = mu
      if ( str(:5) == 'Rmass' ) self%Rmass = mu
    end select

  end subroutine SetScales

!ccccccccccccccc

  real (dp) function scales(self, str)
    class (MatricesElements), intent(in) :: self
    character (len = *)     , intent(in) :: str

    scales = 0;

    if ( str(:3) == 'muH' ) scales = self%muH; if ( str(:1) == 'Q' ) scales = self%Q

    select type (self)
    class is (MatricesElements)
      if ( str(:3) == 'muJ'   ) scales = self%muJ
      if ( str(:3) == 'muS'   ) scales = self%muS
      if ( str(:1) == 'R'     ) scales = self%R
    class is (MatricesElementsMass)
      if ( str(:3) == 'muJ'   ) scales = self%muJ
      if ( str(:3) == 'muS'   ) scales = self%muS
      if ( str(:1) == 'R'     ) scales = self%R
      if ( str(:2) == 'mm'    ) scales = self%mm
      if ( str(:3) == 'muM'   ) scales = self%muM
      if ( str(:5) == 'Rmass' ) scales = self%Rmass
    end select

  end function scales

!ccccccccccccccc

  subroutine MassDeltaScheme(self, scheme, massOrd, m, deltaM)
    class (MatrixElementsMass), intent(in ) :: self
    character (len = *)       , intent(in ) :: scheme
    integer                                 :: massOrd
    real (dp)                 , intent(out) :: m
    real (dp), dimension(3)   , intent(out) :: deltaM

    deltaM = 0

    if ( scheme(:4) == 'pole' ) m = self%alphaMass%scales('mL')

    if ( self%muM < self%muJ ) then
      if ( scheme(:4) /= 'pole' ) then
        m = self%MSmass; deltaM = self%delta('MSbar')
      end if
    else

      deltaM = self%delta(scheme)

      if ( scheme(:3) == 'MSR' .and. scheme(:10) /= 'MSRNatural') then
        m = self%alphaMassNl%MSRmass(self%Rmass)
      else if ( scheme(:10) == 'MSRNatural' ) then
        m = self%alphaMassNl%MSRNaturalMass(self%alphaMass%orders('runMass'), self%Rmass)
      else if ( scheme(:8) == 'MSbarLow' ) then
        m = self%MSLow
      else if ( scheme(:7) == 'JetMass'  ) then
        m = self%JetMass( massOrd, self%alphaMass%orders('runMass'), self%Rmass, self%muJ )
      end if
    end if

  end subroutine MassDeltaScheme

!ccccccccccccccc

  subroutine MassDeltaSet(self, scheme, massOrd, muJ, Rmass, m, deltaM, alphaJ)
    class (MatricesElementsMass), intent(in ) :: self
    character (len = *)         , intent(in ) :: scheme
    integer                     , intent(in ) :: massOrd
    real (dp)                   , intent(in ) :: muJ, Rmass
    real (dp)                   , intent(out) :: m, alphaJ
    real (dp), dimension(3)     , intent(out) :: deltaM
    real (dp), dimension(0:4,3)               :: coef
    real (dp), dimension(3)                   :: lgList, alphaList

    deltaM = 0; coef = 0

    if ( scheme(:4) == 'pole' ) then;  m = self%alphaMass%scales('mL'); return; end if

    if ( self%muM < muJ ) then

      m      = self%alphaMass%MSbarMass(muJ)   ; coef = m * MSbarDeltaPiece(self%nl, 1)
      alphaJ = self%alphaMass%alphaQCD(muJ)/Pi ; call AddAlpha( coef, powList(alphaJ, 3) )

      deltaM = coef(0,:) + matmul(  powList( log(muJ/m), 3 ), coef(1:3,:)  )

    else

      alphaJ = self%alphaAll%alphaQCD(self%nl, muJ)/Pi; alphaList = powList(alphaJ, 3)

      if ( scheme(:8) /= 'MSbarLow' ) lgList = powList( log(muJ/Rmass), 3 )

      if ( scheme(:3) == 'MSR' .and. scheme(:10) /= 'MSRNatural') then

        coef(0,:) = Rmass * self%andimNl%MSRDelta(); call self%andimNl%expandAlpha(coef)
        call AddAlpha(coef, alphaList);  deltaM = DeltaComputer(coef, lgList, 0)
        m = self%alphaMassNl%MSRmass(Rmass)

      else if ( scheme(:10) == 'MSRNatural' ) then

        coef(0,:) = Rmass * MSbarDelta(self%nl, 0);  call self%andimNl%expandAlpha(coef)
        call AddAlpha(coef, alphaList);          deltaM = DeltaComputer(coef, lgList, 0)
        m = self%alphaMassNl%MSRNaturalMass(massOrd, Rmass)

      else if ( scheme(:8) == 'MSbarLow' ) then

        m = self%alphaMassNl%MSbarMassLow(muJ);  coef = m * MSbarDeltaPiece(self%nl, 0)
        call AddAlpha(coef, alphaList)
        deltaM = DeltaComputer(  coef, powList( log(muJ/m), 3 ), 0  )

      else if ( scheme(:7) == 'JetMass'  ) then

        deltaM = alphaList * ExpEuler * Rmass * DeltaComputer( self%BJetExp(1:,:), lgList, 1 )/2
        m = self%JetMass(massOrd, self%alphaMass%orders('runMass'), Rmass, muJ)

      end if
    end if

  end subroutine MassDeltaSet

!ccccccccccccccc

  real (dp) function mmFromJetMass(self, order, run, R, mu, mass)
    class (MatricesElementsMass), intent(inout) :: self
    real (dp)                   , intent(in   ) :: mu, R, mass
    integer                     , intent(in   ) :: order, run
    integer                                     :: IFLAG
    real (dp)                                   :: a0, b0, c0, rat, massOr

    a0 = mass/2; b0 = 2 * mass; massOr = self%mm; rat = 1

    if (self%nf == 6) rat = self%alphaAll%scales('muT')/self%mm
    if (self%nf == 5) rat = self%alphaAll%scales('muB')/self%mm
    if (self%nf == 4) rat = self%alphaAll%scales('muC')/self%mm

    call DFZero(JMass, a0, b0, c0, 1e-9_dp, 1e-9_dp, IFLAG); mmFromJetMass = a0
    call self%SetMass(massOr, rat * massOr)

    contains

!ccccccccccccccc

    real (dp) function JMass(mm)
      real (dp), intent(in) :: mm
      call self%SetMass(mm, rat * mm);   JMass = self%JetMass(order, run, R, mu) - mass
    end function JMass

  end function mmFromJetMass

!ccccccccccccccc

  pure function factList(n) result(list)
    integer, intent(in)       :: n
    integer                   :: i
    real (dp), dimension(0:n) :: list

    list(0) = 1

    do i = 1, n
      list(i) = list(i - 1) * i
    end do

  end function factList

!ccccccccccccccc

  real (dp) function NGLDoubleIntegral(pow, nf, w1, w2, r)
    real (dp), intent(in) :: w1, w2, r
    integer  , intent(in) :: pow, nf
    integer               :: neval, ier
    real (dp)             :: abserr

    call qags( NGLDoubleintegrand, 0._dp, Pi, prec, prec, NGLDoubleIntegral, &
    abserr, neval, ier )

    NGLDoubleIntegral = NGLDoubleIntegral/Pi

    contains

!ccccccccccccccc

  real (dp) function NGLDoubleintegrand(x)
    real (dp), intent(in) :: x
    complex (dp)          :: expTheta, expTheta1

    expTheta = Exp( (0,1) * x );  expTheta1 = 1 + expTheta

    NGLDoubleintegrand = Real( NGLSoft(nf, r * expTheta) * expTheta1**( - 1 - w1 ) * &
    (1 + 1/expTheta)**( - 1 - w2 ) * log(expTheta1)**pow  )

  end function NGLDoubleintegrand

  end function NGLDoubleIntegral

!ccccccccccccccc

  real (dp) function NGLIntegral(self, pow, w1, w2)
    class (MatricesElements), intent(in) :: self
    real (dp)               , intent(in) :: w1, w2
    integer                 , intent(in) :: pow
    integer                              :: neval, ier
    real (dp)                            :: abserr

    call qags( NGLintegrand, 0._dp, Pi, prec, prec, NGLIntegral, abserr, neval, ier )

    NGLIntegral = NGLIntegral/Pi

    contains

!ccccccccccccccc

  real (dp) function NGLintegrand(x)
    real (dp), intent(in) :: x
    real (dp)             :: cTheta, wdiff, ltheta, cdiff, sdiff

    cTheta = 2 * Cos(x/2) ;   wdiff = w1 - w2 ;  cdiff = Cos(wdiff/2 * x)

    if (pow > 0) ltheta = Log(cTheta);  if (pow > 0) sdiff  = Sin(wdiff/2 * x)

    NGLintegrand = cTheta**( - 2 - w1 - w2 ) * self%NGLfunction(x)

    if (pow == 0) NGLintegrand = NGLintegrand *   cdiff
    if (pow == 1) NGLintegrand = NGLintegrand * ( cdiff * ltheta + x * sdiff/2 )
    if (pow == 2) NGLintegrand = NGLintegrand * ( cdiff * (ltheta**2 - x**2/4) + &
                                                      x * ltheta * sdiff )

  end function NGLintegrand

  end function NGLIntegral

!ccccccccccccccc

  real (dp) function NGLfunction(self, y)
    class (MatricesElements), intent(in) :: self
    real (dp)               , intent(in) :: y
    real    (dp)                         :: x, ImPoly2, ImPoly22, ImPoly3, RePoly3
    complex (dp)                         :: expIx

    x = abs(y)

    if (x < 1e-8_dp) then
      NGLfunction = self%SoftNGLThetaApprox(6,x)
    else

      expIx = Exp( (0,1) * x )

      ImPoly2  = aimag( CLi2(expIx) )     ;  RePoly3  =  real( CLi3(expIx ) )
      ImPoly22 = aimag( CLi2(1 - expIx) ) ;  ImPoly3  = aimag( CLi3(1 - expIx) )

      NGLfunction = (  - 1880.0581346047973_dp + 107.39746270332105_dp * self%nf + &
      118.4352528130723_dp*x**2 + 15 * x**4 - (24 - 8 * self%nf) * x/Tan(x/2) +    &
      (6 - 2 * self%nf) * x**2/Sin(x/2)**2 + 32 * (33 - 2 * self%nf) * x * ImPoly2 &
       - 144 * ImPoly22**2 - 288 * x * ImPoly3 + (264 - 16 * self%nf) * x**2 *     &
      Log( 2 * Sin(x/2) ) + (1584 - 96 * self%nf) * RePoly3  )/72

   end if

  end function NGLfunction

!ccccccccccccccc

  complex (dp) function NGLSoft(nf, z)
    complex (dp) , intent(in) :: z
    integer      , intent(in) :: nf
    complex (dp)              :: z1, z2, lz, lz1, CL, CL2, CL3, CL4

    z1 = 1 - z    ;  z2 = - z1/z     ;  lz1 = log(z1);  lz = log(z)
    CL2 = CLi2(z) ;  CL3 = Cli3(1/z) ;  CL4 = Cli3(z);  CL = CLi2(z2)

    NGLSoft = 1._dp/3 - nf/9._dp - 22 * Zeta3 + 36 * nf * Zeta3/27 - 2 * nf &
    * CL4/3 + lz/3/z1 - nf * lz/9/z1 + z * lz/3/z1 - nf * z * lz/9/z1 - Pi2 *  &
    lz**2/3 + z * lz**2/3/z1**2 - nf * z * lz**2/9/z1**2 - 11 * lz1 * lz**2/6 +     &
    nf * lz1 * lz**2/9 - 11 * Log(z2) * lz**2/6 + nf * Log(z2) * lz**2/9 +     &
    lz**4/12 - Pi2 * CLi2(z1)/3 + CLi2(z1)**2 + 2 * (33 - 2 * nf) * lz * CLi2(1/z)  &
    /9 - Pi2 * CL/3  - 22 * lz/3 * CL2 + 4 * nf * lz * CL2/9 + 2 * lz * CLi3(z1) +  &
    11 * CL3 - 2 * nf * CL3/3 - 2 * lz * CLi3(z2) + 11 * CL4 + CL**2

  end function NGLSoft

!ccccccccccccccc

  real (dp) function SoftNGLThetaApprox(self, isoft, x)
    class (MatricesElements), intent(in) :: self
    integer                 , intent(in) :: isoft
    real (dp)               , intent(in) :: x
    integer                              :: i

    SoftNGLThetaApprox = 0

    do i = 1, isoft
      SoftNGLThetaApprox = SoftNGLThetaApprox + self%NGList(i) * x**(2*i) * (-1)**i
    end do

  end function SoftNGLThetaApprox

!ccccccccccccccc

  pure real(dp) function cuspConst(str)
    character (len = *), intent(in) :: str

    cuspConst = 0

    if ( str(:6) == 'thrust' .or. str(:6) == 'Cparam' .or. str(2:3) == 'JM' ) cuspConst = 1
    if ( str(:4) == 'bJet' ) cuspConst = - 0.25_dp

  end function cuspConst

!ccccccccccccccc

  pure integer function gapCons(str)
    character (len = *), intent(in) :: str
    gapCons = 0
    if ( str(:6) == 'thrust' .or. str(:3) == 'SJM') gapCons = 2
    if ( str(:6) == 'Cparam' .or. str(:3) == 'HJM') gapCons = 1
  end function gapCons

!ccccccccccccccc

  pure function posToMomMatrix() result(PosToMom)
    real (dp), dimension(0:6,0:6) :: PosToMom

    PosToMom = 0;   PosToMom(0,2) = - 1.6449340668482262_dp;   PosToMom(0,0) = 1
    PosToMom(0,3) =  - 2.4041138063191885_dp;    PosToMom(1,1) =  - 1
    PosToMom(0,4) =    1.623484850566707_dp ;    PosToMom(2,2) =    2
    PosToMom(0,5) =   14.659820882505041_dp ;    PosToMom(3,3) =  - 3
    PosToMom(0,6) =   29.184858319032745_dp ;    PosToMom(4,4) =    4
    PosToMom(1,3) =    4.934802200544679_dp ;    PosToMom(5,5) =  - 5
    PosToMom(1,4) =    9.616455225276754_dp ;    PosToMom(6,6) =    6
    PosToMom(1,5) = -  8.117424252833533_dp ;    PosToMom(1,6) = - 87.95892529503024_dp
    PosToMom(2,4) = - 19.739208802178716_dp ;    PosToMom(2,5) = - 48.08227612638377_dp
    PosToMom(2,6) =   48.70454551700121_dp  ;    PosToMom(3,5) =   49.34802200544679_dp
    PosToMom(3,6) =  144.2468283791513_dp   ;    PosToMom(4,6) = - 98.69604401089359_dp

  end function posToMomMatrix

!ccccccccccccccc

  pure subroutine AddLogs(self, Mat, lg)
    class (MatricesElements) , intent(in   ) :: self
    real (dp)                , intent(in   ) :: lg
    real (dp), dimension(:,:), intent(inout) :: Mat
    integer                                  :: i, k, n
    real (dp), dimension( size(Mat,1) - 1 )  :: lgList

    n = size(Mat,1); lgList = powList(lg,n - 1)

    non_dist_logs: do i = 0, n - 1
      sum_logs: do k = i + 1, n - 1

       Mat(i + 1,:) = Mat(i + 1,:) + Mat(k + 1,:) * lgList(k - i) * self%ListFact(k)/&
       self%ListFact(i)/self%ListFact(k - i)

      enddo sum_logs
    enddo non_dist_logs

  end subroutine AddLogs

!ccccccccccccccc

  function HardMassVect(self, mass) result(HmExp)
    class (MatricesElementsMass), intent(in) :: self
    real (dp)                   , intent(in) :: mass
    real (dp)                 , dimension(3) :: HmExp
    real (dp)                                :: Lm, LQ

    HmExp = 0; Lm = 2 * log(mass/self%muM); LQ = 2 * log(self%Q/mass)

    HmExp(1) = 3.7632893778988175_dp + 2 * Lm * (Lm - 1)/3

    HmExp(2) = 22.02886430312693_dp - 28 * LQ/81 - Lm**2 * (LQ/9 - 0.5473424096795991_dp &
    + 13._dp * self%nf/54) + Lm**3/54 * (2 * self%nf - 37) &
    - 2.4985480650201737_dp * self%nf - Lm * (9.51716816196516_dp + 10 * LQ/27 &
    - 1.2063904494634092_dp * self%nf)

    HmExp(3) = - lQ * (4 * Lm**2 + 40 * Lm/3 + 112._dp/9)/36

  end function HardMassVect

!ccccccccccccccc

  function HardMassVec(self, mass) result(HmExp)
    class (MatricesElementsMass), intent(in) :: self
    real (dp)                   , intent(in) :: mass
    real (dp)                 , dimension(3) :: HmExp
    real (dp)                 , dimension(2) :: alphaList

    HmExp = self%HardMassVect(mass); alphaList = powList( self%alphaList(3), 2)

    HmExp(:2) = alphaList * HmExp(:2)
    HmExp(1)  = HmExp(1) + alphaList(2) * HmExp(3); HmExp(3) = 0

  end function HardMassVec

!ccccccccccccccc

  function HardMassExp(self, mass) result(HmExp)
    class (MatricesElementsMass), intent(in) :: self
    real (dp)                   , intent(in) :: mass
    real (dp)             , dimension(0:4,3) :: HmExp

    HmExp = 0; HmExp(0,:) = self%HardMassVec(mass)

  end function HardMassExp

!ccccccccccccccc

  pure function GenMat(ES, nf, s3) result(res)
    character (len = *)  , intent(in) :: ES
    integer              , intent(in) :: nf
    real (dp)            , intent(in) :: s3
    real (dp)          , dimension(3) :: res

    res = 0

    if ( ES(:10) == 'SoftThrust' ) then

      res = [ - 3.2898681336964524_dp, - 14.12474197259218_dp + 1.8079382571738185_dp * nf, &
             s3/32 ]

    else if ( ES(:10) == 'SoftCparam' ) then

      res = [ - 1.096622711232151_dp, (14.60597_dp * nf -  115.9508_dp)/8, s3/32 ]

    else if ( ES(:3) == 'jet' ) then

      res = [ 0.1400879108690316_dp, - 0.03381606438755153_dp &
            - 0.44946310960521874_dp * nf, 0.004279026412919151_dp + s3/64 + &
              0.06296434803729367_dp * nf ]

    else if ( ES(:4) == 'bJet' ) then

      res = [ 1.8816446889494087_dp, 15.768694405652397_dp - 1.2102726427482482_dp * nf, &
              s3/32 ]

    else if ( ES(:4) == 'hard' ) then

      res = [ - 2.118355311050591_dp, - 12.39778859389983_dp + 1.4232390700560242_dp * nf,&
       - 128.52222623172312_dp + 23.31580959050584_dp * nf - 0.7731256344257622_dp * nf**2 ]

    end if

  end function GenMat

!ccccccccccccccc

  pure function NGLOas2(nf) result(res)
    integer     , intent(in) :: nf
    real (dp) , dimension(6) :: res

    res = [ - 1.8949340668482262_dp - nf/12._dp                 , &
    - 0.01851851851851849_dp    + 0.0007716049382716084_dp  * nf, &
      0.00034567901234566767_dp - 5.1440329218114395e-6_dp  * nf, &
    - 6.636252384548125e-6_dp   - 1.7059292852904945e-7_dp  * nf, &
      1.2457657960415758e-7_dp  + 1.1227055980202183e-8_dp  * nf, &
    - 2.2144753634989467e-9_dp  - 4.4695740912623766e-10_dp * nf ]

  end function NGLOas2

!ccccccccccccccc

  function SetHardFunc(self, lgH, alphaIm, coefCusp, coefHard) result(HardExp)
    class (MatricesElements)   , intent(in) :: self
    real (dp), dimension(0:4,3), intent(in) :: coefCusp, coefHard
    real (dp)                  , intent(in) :: lgH, alphaIm
    real (dp), dimension(0:4)               :: lgList, PiList
    real (dp)                               :: cumul
    integer                                 :: i, j
    real (dp), dimension(3)                 :: HardExp

    HardExp = self%HardExp(0,:)

    if (self%muH > 0) then

      lgList(0) = 1; lgList(1:) = powList(lgH , 4)
      PiList(0) = 1; PiList(1:) = powList(pio2, 4)

      log_hard : do i = 1, 4

         cumul = 0

         sum_pi: do j = 0, i/2

           cumul = cumul + self%ListFact(i)/self%ListFact(2 * j)/self%ListFact(i - 2 * j) *  &
           (-1)**(i - j) * lgList(i - 2 * j) * PiList(2 * j)

         end do sum_pi

         HardExp = HardExp + cumul * ( self%HardExp(i,:) + coefHard(i,:)/2 - coefCusp(i,:) )

       enddo log_hard

     else

       order_do: do i = 1, 3

         cumul = 0

         alpha_do: do j = 1, i/2

           cumul = cumul + self%ListFact(i)/self%ListFact(2 * j)/self%ListFact(i - 2 * j) &
           * (-1)**j * alphaIm**j

         end do alpha_do

           lgList = powList(-lgH, 5)

           log_loop: do j = 1, 4

           HardExp(i) = HardExp(i) + ( self%HardExp(j,i) + coefHard(j,i)/2 - &
           coefCusp(j,i) ) * lgList(j)

         end do log_loop

         HardExp(i) = HardExp(i) * cumul

       end do order_do
     end if
   end function SetHardFunc

!ccccccccccccccc

end module MatrixElementsClass
