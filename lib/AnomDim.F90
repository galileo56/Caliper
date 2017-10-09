
module AnomDimClass
  use Constants, only: dp, d1mach, Pi, Pi2, ExpEuler, prec
  use QuadPack, only: qags, qagi; use adapt, only: dGauss
  implicit none ;  private

  real (dp), parameter :: al = 0.634_dp, bet = 1.035_dp, gam = - 23.6_dp, &
  ep = 1.19_dp, del = - 0.481_dp, fourPi = 4 * Pi

  public :: inteCorre, alphaReExpand, deltaMass, MSbarDelta, deltaCharmGlue, &
  PowList, getInverse, MSbarDeltaPiece, AlphaExpand, alphaMatchingLog, pint, &
  deltaCharm2, deltaCharm2Der, gammaRcharm2, gammaRcharm3, deltaCharm3, p,   &
  deltaCharm3Der, deltaCharmNh, deltaCharmNhDer, P2int, DeltaBottomCharm,    &
  GammaRBottomCharm, deltaCharmNL, deltaCharmNLDer, deltaCharmGlueDer, &
  alphaMatch, PowList0, factList

  interface PowList
    module procedure   :: PowListDP, PowListInt, PowListComp
  end interface PowList

  interface PowList0
    module procedure   :: PowList0DP, PowList0Int, PowList0Comp
  end interface PowList0

!ccccccccccccccc

  type, public :: AnomDim
    private
    integer                       :: nf
    character (len = 5)           :: str
    real (dp)                     :: G4, err, beta0QED
    real (dp), dimension(5)       :: AlphaMatch, AlphaMatchInv, AlphaMatchUp
    real (dp), dimension(0:4,5)   :: AlphaMatchLog, AlphaMatchInvLog
    real (dp), dimension(0:3,0:3) :: gammaHm
    real (dp), dimension(0:4)     :: bHat, bCoef, cCoef
    real (dp), dimension(0:4)     :: beta, beta2List
    real (dp), dimension(0:4)     :: gammaMass
    real (dp), dimension(0:3)     :: gammaHard, gammaB, gammaJet, gammaSoft, &
    cusp, gtilde, gl
    real (dp), dimension(4)       :: sCoefMSRp, sCoefMSRn, betaList, gammaRp, &
    gammaRn, sCoefMSRInc1, gammaRInc1, sCoefMSRInc2, gammaRInc2, gammaRInc3,  &
    sCoefMSRInc3, aRS, aRSn, aRSp, sCoefRS, gammaRS

  contains

    procedure, pass(self), public :: expandAlpha, wTildeExpand, kTildeExpand, &
    sCoef, DeltaMu, betaQCD, numFlav, DeltaR, DeltaRHadron, Gfun, DeltaRMass, &
    bHQETgamma, wTildeHm, GammaRComputer, sCoefRecursive, N12Generic, PScoef, &
    sCoefHadron, scheme, MSRDelta, sCoefLambda, P12, sCoefGeneric, cCoeff,    &
    P12Generic, MatchingAlpha, MatchingAlphaLog, MatchingAlphaUp, betaQED,    &
    alphaMatching, N12, N12Ratio, N12Residue

    procedure, pass(self), private ::  wTildeReal, wTildeComplex, kTildeReal, &
    kTildeComplex, rootReal, rootComplex

    generic, public                :: wTilde => wTildeReal, wTildeComplex
    generic, public                :: kTilde => kTildeReal, kTildeComplex
    generic, public                :: root   => rootReal, rootComplex

  end type AnomDim

!ccccccccccccccc

  interface AnomDim
    module procedure InAdim
  end interface AnomDim

  contains

!ccccccccccccccc

   type (AnomDim) function InAdim(str, nf, G4, err)
    character (len = *), intent(in) :: str
    integer            , intent(in) :: nf         ! number of active flavors
    real (dp)          , intent(in) :: G4         ! cusp anomalous dimension
    real (dp), optional, intent(in) :: err        ! 4-loop MS-bar to pole conversion error
    real (dp), dimension(0:4)       :: beta2List, beta
    real (dp), dimension(4)         :: betaList
    real (dp), dimension(-3:6)      :: poch
    integer                         :: n, k, l
    integer, dimension(0:4)         :: signlist
    real (dp)                       :: suma

    InAdim%err = 0; if ( present(err) ) InAdim%err = err

    beta = betaFun(nf);    InAdim%bCoef = beta/beta(0); InAdim%beta = beta

    InAdim%str = str ; InAdim%nf   = nf   ; InAdim%gammaHm = 0
    InAdim%G4  = G4  ; InAdim%GammaMass = massAnomDim(nf)

    InAdim%cusp = [ 16._dp/3, 66.47322097196788_dp - 5.925925925925926_dp * nf, &
    1174.8982718294073_dp - 183.18743044213898_dp * nf - 0.7901234567901234_dp * nf**2, G4 ]

    InAdim%gammaJet = [ 8._dp, 11.643665868860921_dp - 17.79927174385542_dp * nf, &
    407.55606459953486_dp - 170.69272520622684_dp * nf + 1.4510585125043_dp * nf**2, 0._dp ]

    InAdim%gammaB = [ 16._dp/3, 9.7609731286894_dp - 8.532462893504388_dp * nf, &
    154.56312636466896_dp - 108.96471855555455_dp * nf - 1.1952605707205706_dp * nf**2, 0._dp ]

    InAdim%gammaHard = [ - 8._dp, - 74.8217346132464_dp + 15.19273477627696_dp * nf, &
    - 1499.402316365377_dp + 259.2839279841602_dp * nf - 1.8561956264347472_dp * nf**2, 0._dp ]

    InAdim%gammaSoft    = - InAdim%gammaHard - InAdim%gammaJet
    InAdim%gammaHm(:,0) = - InAdim%gammaSoft - InAdim%gammaB

    if (str(:5) == 'MSbar') then

      InAdim%gammaHm(:1,1) = [ 512._dp/9, 2142.6359241459822_dp - 177.726619506625_dp * nf ]
      InAdim%gammaHm(0,2)  =  709.0476903676573_dp -  63.20987654320987_dp * nf

    end if

    betaList = PowList( 1/beta(0)/2, 4 ); InAdim%betaList = betaList
    beta2List = PowList0( beta(0)/2, 4 ); InAdim%beta2List = beta2List

    InAdim%bHat(0) = 1; InAdim%cCoef(0) = 1

    do n = 0, 3
      InAdim%cCoef(n+1) = - sum( InAdim%cCoef(:n) * InAdim%bCoef(n+1:1:-1) )
      InAdim%bHat(n+1)  = - (-1)**n * InAdim%cCoef(n+1) * betaList(n+1)
    end do

    InAdim%gl(0) = 1; InAdim%gTilde(0) = 1

    do n = 0, 2
      InAdim%gl(n + 1) = - sum( powList(-1,n + 1) * InAdim%bHat(2:n+2) * InAdim%gl(n:0:-1) )/(n+1)
      InAdim%gTilde(n + 1) = sum( powList(-1,n + 1) * InAdim%bHat(2:n+2) * InAdim%gTilde(n:0:-1) )/(n+1)
    end do

    poch(0:-3:-1) = 1/PochHammerList( - InAdim%bHat(1), 3 ) * powList0(-1,3)
    poch(0:3) = PochHammerList( 1 + InAdim%bHat(1), 3 )

    do n = 1, 4
      InAdim%aRS(n) =  Pi * sum( InAdim%gl * poch(n-1:n-4:-1) ) * beta2List(n - 1)
    end do

    InAdim%gammaRS       = InAdim%GammaRComputer( InAdim%aRS )
    InAdim%gammaRp       = InAdim%GammaRComputer( InAdim%MSRDelta('MSRp') )
    InAdim%gammaRn       = InAdim%GammaRComputer( MSbarDelta(nf    , 0, InAdim%err) )
    InAdim%gammaRInc1    = InAdim%GammaRComputer( MSbarDelta(nf    , 1, InAdim%err) )
    InAdim%gammaRInc2    = InAdim%GammaRComputer( MSbarDelta(nf - 1, 1, InAdim%err) )
    InAdim%gammaRInc3    = InAdim%GammaRComputer( MSbarDelta(nf + 1, 1, InAdim%err) )

    InAdim%sCoefRS       = InAdim%sCoefRecursive( InAdim%aRS )
    InAdim%sCoefMSRp     = InAdim%sCoefRecursive( InAdim%MSRDelta('MSRp') )
    InAdim%sCoefMSRn     = InAdim%sCoefRecursive( MSbarDelta(nf   ,  0, InAdim%err) )
    InAdim%sCoefMSRInc1  = InAdim%sCoefRecursive( MSbarDelta(nf   ,  1, InAdim%err) )
    InAdim%sCoefMSRInc2  = InAdim%sCoefRecursive( MSbarDelta(nf - 1, 1, InAdim%err) )
    InAdim%sCoefMSRInc3  = InAdim%sCoefRecursive( MSbarDelta(nf + 1, 1, InAdim%err) )

    InAdim%AlphaMatch    = InAdim%alphaMatching(nf)
    InAdim%AlphaMatchInv = getInverse( InAdim%AlphaMatch )
    InAdim%AlphaMatchUp  = InAdim%alphaMatching(nf + 1)
    InAdim%AlphaMatchLog = alphaMatchingLog(InAdim%AlphaMatch, nf, nf - 1)
    InAdim%AlphaMatchInvLog = alphaMatchingLog(InAdim%AlphaMatchInv, nf - 1, nf)

    InAdim%beta0QED = 4 ! ( electron + muon + up + down + strange )

    if (nf >= 4) InAdim%beta0QED = InAdim%beta0QED + 4._dp/3 + 1 ! add tau and charm
    if (nf >= 5) InAdim%beta0QED = InAdim%beta0QED + 1._dp/3     ! add bottom
    if (nf == 6) InAdim%beta0QED = InAdim%beta0QED + 4._dp/3     ! add top

    InAdim%beta0QED = 4 * InAdim%beta0QED/3

    InAdim%aRSn = 0; InAdim%aRSp = 0; signlist = PowList0(-1,4)

    do n = 1, 4
      do k = 0, 3
        poch(0:4 + k - n) = PochHammerList( - InAdim%bHat(1) - k, 4 + k - n )
        l = max(n - k, 0)
        suma = sum( InAdim%gl(l:) * signlist(l:)/poch(1+k+l-n:4+k-n) )
        InAdim%aRSn(n) = InAdim%aRSn(n) + InAdim%sCoefMSRn(k+1) * suma * (-1)**k
        InAdim%aRSp(n) = InAdim%aRSp(n) + InAdim%sCoefMSRp(k+1) * suma * (-1)**k
      end do
    end do

    InAdim%aRSn = InAdim%MSRDelta('MSRn') - InAdim%aRSn * beta2List(1:) * signlist(1:)
    InAdim%aRSp = InAdim%MSRDelta('MSRp') - InAdim%aRSp * beta2List(1:) * signlist(1:)

   end function InAdim

!ccccccccccccccc

  pure function cCoeff(self, order, m) result(res)
    class (AnomDim)              , intent(in) :: self
    integer                      , intent(in) :: order, m
    real (dp)  , dimension( 0:m )             :: res
    real (dp)  , dimension( min(order, 4) )   :: bCoef, ListPow
    integer                                   :: n, ord

    ord = min(order, 4); ListPow = powList(1/fourPi, ord); res(0) = 1

    res(1:ord) = self%cCoef(1:ord) * ListPow; bCoef = self%bCoef(1:ord) * ListPow

    do n = ord, m - 1
      res(n + 1) = - sum( res(n + 1 - ord:n) * bCoef(ord:1:-1) )
    end do

  end function cCoeff

!ccccccccccccccc

   pure character (len = 5) function scheme(self)
    class (AnomDim), intent(in) :: self
    scheme = self%str
   end function scheme

!ccccccccccccccc

   pure real (dp) function rootReal(self, order, aMu, a0, aLLInv)
    class (AnomDim), intent(in) :: self
    integer        , intent(in) :: order
    real (dp)      , intent(in) :: aMu, a0, aLLInv
    real(dp)                    :: aSuPi, aSuPi0, radical, raiz

    aSuPi = aMu/fourPi; aSuPi0 = a0/fourPi

    rootReal = aLLInv - self%bCoef(1)/fourPi * log(aMu/a0)

    if ( order == 1 ) then

      rootReal = rootReal + self%bCoef(1)/fourPi * &
      log( (1 + self%bCoef(1) * aSuPi)/(1 + self%bCoef(1) * aSuPi0)  )

    else if ( order == 2 ) then

      radical = self%bCoef(1)**2 - 4 * self%bCoef(2); raiz = sqrt( abs(radical) )

      rootReal = rootReal + self%bCoef(1)/fourPi/2 * &
      log( (1 + self%bCoef(1) * aSuPi + self%bCoef(2) * aSuPi**2)/&
      (1 + self%bCoef(1) * aSuPi0 + self%bCoef(2) * aSuPi0**2)  )

      if (radical < 0) then

        rootReal = rootReal + ( self%bCoef(1)**2 - 2 * self%bCoef(2) )/fourPi/raiz * &
        (  atan( (self%bCoef(1) + 2 * self%bCoef(2) * aSuPi)/raiz ) - &
        atan( (self%bCoef(1) + 2 * self%bCoef(2) * aSuPi0)/raiz )  )

      else

        rootReal = rootReal - ( self%bCoef(1)**2 - 2 * self%bCoef(2) )/fourPi/raiz * &
        (  atanh( (self%bCoef(1) + 2 * self%bCoef(2) * aSuPi)/raiz ) - &
        atanh( (self%bCoef(1) + 2 * self%bCoef(2) * aSuPi0)/raiz )  )

      end if

    else if ( order == 3 ) then
      select case ( self%nf )
      case (1)
        rootReal = rootReal - 0.6404675122670069_dp * (  &
        ATan( 0.5243427048805052_dp * (2 * aMu - 0.32585902383203375_dp) )   -  &
        ATan( 0.5243427048805052_dp * (2 * a0  - 0.32585902383203375_dp) )  ) + &
        0.44444989158991516_dp * log( (aMu + 0.9651061041938452_dp)/&
        (a0 + 0.9651061041938452_dp) ) + 0.12175509250042996_dp * &
        Log( (0.9358509616463481_dp - 0.32585902383203375_dp * aMu + aMu**2)/&
        (0.9358509616463481_dp - 0.32585902383203375_dp * a0 + a0**2) )
      case (2)
        rootReal = rootReal - 0.5909374157641947_dp * (  &
        ATan( 0.48547873662772983_dp * (2 * aMu - 0.37047299521078925_dp) )   -  &
        ATan( 0.48547873662772983_dp * (2 * a0  - 0.37047299521078925_dp) )  ) + &
        0.4177502777889539_dp * log( (aMu + 1.0315084572808646_dp)/&
        (a0 + 1.0315084572808646_dp) ) + 0.10669069654635029_dp * &
        Log( (1.0950296899338223_dp - 0.37047299521078925_dp * aMu + aMu**2)/&
        (1.0950296899338223_dp - 0.37047299521078925_dp * a0 + a0**2) )
      case (3)
        rootReal = rootReal - 0.5342835587527397_dp * (  &
        ATan( 0.44205618349190395_dp * (2 * aMu - 0.4428447212033426_dp) )   -  &
        ATan( 0.44205618349190395_dp * (2 * a0  - 0.4428447212033426_dp) )  ) + &
        0.3907037552748204_dp * log( (aMu + 1.1120254691513045_dp)/&
        (a0 + 1.1120254691513045_dp) ) + 0.08759024341484811_dp * &
        Log( (1.3283651814941193_dp - 0.4428447212033426_dp * aMu + aMu**2)/&
        (1.3283651814941193_dp - 0.4428447212033426_dp * a0 + a0**2) )
      case (4)
        rootReal = rootReal - 0.46716448667125304_dp * (  &
        ATan( 0.3928365460051735_dp * (2 * aMu - 0.5735173459270638_dp) )   -  &
        ATan( 0.3928365460051735_dp * (2 * a0  - 0.5735173459270638_dp) )  ) + &
        0.36500373291627086_dp * log( (aMu + 1.2090182171482846_dp)/&
        (a0 + 1.2090182171482846_dp) ) + 0.06259674590338339_dp * &
        Log( (1.7022351111458442_dp - 0.5735173459270638_dp * aMu + aMu**2)/&
        (1.7022351111458442_dp - 0.5735173459270638_dp * a0 + a0**2) )
      case (5)
        rootReal = rootReal - 0.3836359569088666_dp * (  &
        ATan( 0.3365914582110142_dp * (2 * aMu - 0.8495407965950523_dp) )   -  &
        ATan( 0.3365914582110142_dp * (2 * a0  - 0.8495407965950523_dp) )  ) + &
        0.3441107587583152_dp * log( (aMu + 1.320588432330028_dp)/&
        (a0 + 1.320588432330028_dp) ) + 0.02861824451931911_dp * &
        Log( (2.3870817866590763_dp - 0.8495407965950523_dp * aMu + aMu**2)/&
        (2.3870817866590763_dp - 0.8495407965950523_dp * a0 + a0**2) )
      case (6)
        rootReal = rootReal - 0.27090635085432485_dp * (  &
        ATan( 0.27520884668510553_dp * (2 * aMu - 1.5929881907058037_dp) )   -  &
        ATan( 0.27520884668510553_dp * (2 * a0  - 1.5929881907058037_dp) )  ) + &
        0.3341470572175008_dp * log( (aMu + 1.4277939462375213_dp)/&
        (a0 + 1.4277939462375213_dp) ) - 0.01928679573770475_dp * &
        Log( (3.9351725745955957_dp - 1.5929881907058037_dp * aMu + aMu**2)/&
        (3.9351725745955957_dp - 1.5929881907058037_dp * a0 + a0**2) )
      case default
        rootReal = 0
      end select

    else if ( order > 3 ) then
      select case ( self%nf )
      case (1)
        rootReal = rootReal - 0.5489848787403228_dp * (  &
        ATan( 0.6354705274484186_dp * (2 * aMu - 0.7677302261302864_dp) )   -  &
        ATan( 0.6354705274484186_dp * (2 * a0  - 0.7677302261302864_dp) )  ) - &
        0.23200408345167492_dp * (  &
        ATan( 0.8864954614918331_dp * (2 * aMu + 1.542321166011278_dp) )    -  &
        ATan( 0.8864954614918331_dp * (2 * a0  + 1.542321166011278_dp) )  ) + &
        0.04040335124693508_dp * &
        log( (0.7664358581916029_dp - 0.7677302261302864_dp * aMu + aMu**2)/&
        (0.7664358581916029_dp - 0.7677302261302864_dp * a0 + a0**2) ) + &
        0.30357668704845114_dp * &
        Log( (0.9651061041938452_dp + 1.542321166011278_dp * aMu + aMu**2)/&
        (0.9651061041938452_dp + 1.542321166011278_dp * a0 + a0**2) )
      case (2)
        rootReal = rootReal - 0.5088583673947228_dp * (  &
        ATan( 0.5807020947660965_dp * (2 * aMu - 0.8268729922000535_dp) )   -  &
        ATan( 0.5807020947660965_dp * (2 * a0  - 0.8268729922000535_dp) )  ) - &
        0.1812488725394936_dp * (  &
        ATan( 0.8171726383562911_dp * (2 * aMu + 1.7419567282921797_dp) )    -  &
        ATan( 0.8171726383562911_dp * (2 * a0  + 1.7419567282921797_dp) )  ) + &
        0.03843876906477947_dp * &
        log( (0.9122966896897218_dp - 0.8268729922000535_dp * aMu + aMu**2)/&
        (0.9122966896897218_dp - 0.8268729922000535_dp * a0 + a0**2) ) + &
        0.2771270663760476_dp * &
        Log( (1.1329830828367182_dp + 1.7419567282921797_dp * aMu + aMu**2)/&
        (1.1329830828367182_dp + 1.7419567282921797_dp * a0 + a0**2) )
      case (3)
        rootReal = rootReal - 0.4648829323834203_dp * (  &
        ATan( 0.5186296543090279_dp * (2 * aMu - 0.9045061008589117_dp) )   -  &
        ATan( 0.5186296543090279_dp * (2 * a0  - 0.9045061008589117_dp) )  ) - &
        0.10665133580718698_dp * (  &
        ATan( 0.7478910571108373_dp * (2 * aMu + 2.069827558808925_dp) )    -  &
        ATan( 0.7478910571108373_dp * (2 * a0  + 2.069827558808925_dp) )  ) + &
        0.03534462037703361_dp * &
        log( (1.1339812941648535_dp - 0.9045061008589117_dp * aMu + aMu**2)/&
        (1.1339812941648535_dp - 0.9045061008589117_dp * a0 + a0**2) ) + &
        0.24759750067522504_dp * &
        Log( (1.5180010453345267_dp + 2.069827558808925_dp * aMu + aMu**2)/&
        (1.5180010453345267_dp + 2.069827558808925_dp * a0 + a0**2) )
      case (4)
        rootReal = rootReal - 0.4162619666563603_dp * (  &
        ATan( 0.4467310336890205_dp * (2 * aMu - 1.0104662655416288_dp) )   -  &
        ATan( 0.4467310336890205_dp * (2 * a0  - 1.0104662655416288_dp) )  ) + &
        0.03691444757280972_dp * (  &
        ATan( 0.7190586675946767_dp * (2 * aMu + 2.742110429973893_dp) )    -  &
        ATan( 0.7190586675946767_dp * (2 * a0  + 2.742110429973893_dp) )  ) + &
        0.029930532499675307_dp * &
        log( (1.507962493506846_dp - 1.0104662655416288_dp * aMu + aMu**2)/&
        (1.507962493506846_dp - 1.0104662655416288_dp * a0 + a0**2) ) + &
        0.21516807986184325_dp * &
        Log( (2.3633089675134706_dp + 2.742110429973893_dp * aMu + aMu**2)/&
        (2.3633089675134706_dp + 2.742110429973893_dp * a0 + a0**2) )
      case (5)
        rootReal = rootReal - 0.36091735121516355_dp * (  &
        ATan( 0.3626874145171533_dp * (2 * aMu - 1.163068663777_dp) )   -  &
        ATan( 0.3626874145171533_dp * (2 * a0  - 1.163068663777_dp) )  ) + &
        0.3943399946498003_dp * log( (aMu + 1.5721862037478178_dp)/&
        (a0 + 1.5721862037478178_dp) ) - &
        0.027723547491953493_dp * log( (aMu + 3.511042056691402_dp)/&
        (a0 + 3.511042056691402_dp) ) + 0.017365400319553415_dp * &
        Log( (2.2387135279454564_dp - 1.163068663777_dp * aMu + aMu**2)/&
        (2.2387135279454564_dp - 1.163068663777_dp * a0 + a0**2) )
      case (6)
        rootReal = rootReal - 0.2706517070577726_dp * (  &
        ATan( 0.27659246828604006_dp * (2 * aMu - 1.6112137089692589_dp) )   -  &
        ATan( 0.27659246828604006_dp * (2 * a0  - 1.6112137089692589_dp) )  ) + &
        0.33463363709028904_dp * log( (aMu + 1.4322403491695925_dp)/&
        (a0 + 1.4322403491695925_dp) ) - &
        3.242619128589197e-8_dp * log( (aMu + 114.63875643838061_dp)/&
        (a0 + 114.63875643838061_dp) ) - 0.019530069461003236_dp * &
        Log( (3.916831300483894_dp - 1.6112137089692589_dp * aMu + aMu**2)/&
        (3.916831300483894_dp - 1.6112137089692589_dp * a0 + a0**2) )
      case default
        rootReal = 0
      end select
    end if

   end function rootReal

!ccccccccccccccc

   pure complex (dp) function rootComplex(self, order, aMu, a0, aLLInv)
    class (AnomDim), intent(in) :: self
    integer        , intent(in) :: order
    complex (dp)   , intent(in) :: aMu, aLLInv
    real (dp)      , intent(in) :: a0
    real (dp)                   :: aSuPi0
    complex (dp)                :: aSuPi, raiz

    aSuPi = aMu/fourPi; aSuPi0 = a0/fourPi

    rootComplex = aLLInv - self%bCoef(1)/fourPi * log(aMu/a0)

    if ( order == 1 ) then

      rootComplex = rootComplex + self%bCoef(1)/fourPi * &
      log( (1 + self%bCoef(1) * aSuPi)/(1 + self%bCoef(1) * aSuPi0)  )

    else if ( order == 2 ) then

      raiz = sqrt( self%bCoef(1)**2 - 4 * self%bCoef(2) )

      rootComplex = rootComplex + self%bCoef(1)/fourPi/2 * &
      log( (1 + self%bCoef(1) * aSuPi + self%bCoef(2) * aSuPi**2)/&
      (1 + self%bCoef(1) * aSuPi0 + self%bCoef(2) * aSuPi0**2)  )+&
      ( self%bCoef(1)**2 - 2 * self%bCoef(2) )/fourPi/raiz * &
      (  atan( (self%bCoef(1) + 2 * self%bCoef(2) * aSuPi)/raiz ) - &
      atan( (self%bCoef(1) + 2 * self%bCoef(2) * aSuPi0)/raiz )  )

    else if ( order == 3 ) then
      select case ( self%nf )
      case (1)
        rootComplex = rootComplex - 0.6404675122670069_dp * (  &
        ATan( 0.5243427048805052_dp * (2 * aMu - 0.32585902383203375_dp) )   -  &
        ATan( 0.5243427048805052_dp * (2 * a0  - 0.32585902383203375_dp) )  ) + &
        0.44444989158991516_dp * log( (aMu + 0.9651061041938452_dp)/&
        (a0 + 0.9651061041938452_dp) ) + 0.12175509250042996_dp * &
        Log( (0.9358509616463481_dp - 0.32585902383203375_dp * aMu + aMu**2)/&
        (0.9358509616463481_dp - 0.32585902383203375_dp * a0 + a0**2) )
      case (2)
        rootComplex = rootComplex - 0.5909374157641947_dp * (  &
        ATan( 0.48547873662772983_dp * (2 * aMu - 0.37047299521078925_dp) )   -  &
        ATan( 0.48547873662772983_dp * (2 * a0  - 0.37047299521078925_dp) )  ) + &
        0.4177502777889539_dp * log( (aMu + 1.0315084572808646_dp)/&
        (a0 + 1.0315084572808646_dp) ) + 0.10669069654635029_dp * &
        Log( (1.0950296899338223_dp - 0.37047299521078925_dp * aMu + aMu**2)/&
        (1.0950296899338223_dp - 0.37047299521078925_dp * a0 + a0**2) )
      case (3)
        rootComplex = rootComplex - 0.5342835587527397_dp * (  &
        ATan( 0.44205618349190395_dp * (2 * aMu - 0.4428447212033426_dp) )   -  &
        ATan( 0.44205618349190395_dp * (2 * a0  - 0.4428447212033426_dp) )  ) + &
        0.3907037552748204_dp * log( (aMu + 1.1120254691513045_dp)/&
        (a0 + 1.1120254691513045_dp) ) + 0.08759024341484811_dp * &
        Log( (1.3283651814941193_dp - 0.4428447212033426_dp * aMu + aMu**2)/&
        (1.3283651814941193_dp - 0.4428447212033426_dp * a0 + a0**2) )
      case (4)
        rootComplex = rootComplex - 0.46716448667125304_dp * (  &
        ATan( 0.3928365460051735_dp * (2 * aMu - 0.5735173459270638_dp) )   -  &
        ATan( 0.3928365460051735_dp * (2 * a0  - 0.5735173459270638_dp) )  ) + &
        0.36500373291627086_dp * log( (aMu + 1.2090182171482846_dp)/&
        (a0 + 1.2090182171482846_dp) ) + 0.06259674590338339_dp * &
        Log( (1.7022351111458442_dp - 0.5735173459270638_dp * aMu + aMu**2)/&
        (1.7022351111458442_dp - 0.5735173459270638_dp * a0 + a0**2) )
      case (5)
        rootComplex = rootComplex - 0.3836359569088666_dp * (  &
        ATan( 0.3365914582110142_dp * (2 * aMu - 0.8495407965950523_dp) )   -  &
        ATan( 0.3365914582110142_dp * (2 * a0  - 0.8495407965950523_dp) )  ) + &
        0.3441107587583152_dp * log( (aMu + 1.320588432330028_dp)/&
        (a0 + 1.320588432330028_dp) ) + 0.02861824451931911_dp * &
        Log( (2.3870817866590763_dp - 0.8495407965950523_dp * aMu + aMu**2)/&
        (2.3870817866590763_dp - 0.8495407965950523_dp * a0 + a0**2) )
      case (6)
        rootComplex = rootComplex - 0.27090635085432485_dp * (  &
        ATan( 0.27520884668510553_dp * (2 * aMu - 1.5929881907058037_dp) )   -  &
        ATan( 0.27520884668510553_dp * (2 * a0  - 1.5929881907058037_dp) )  ) + &
        0.3341470572175008_dp * log( (aMu + 1.4277939462375213_dp)/&
        (a0 + 1.4277939462375213_dp) ) - 0.01928679573770475_dp * &
        Log( (3.9351725745955957_dp - 1.5929881907058037_dp * aMu + aMu**2)/&
        (3.9351725745955957_dp - 1.5929881907058037_dp * a0 + a0**2) )
      case default
        rootComplex = 0
      end select

    else if ( order > 3 ) then
      select case ( self%nf )
      case (1)
        rootComplex = rootComplex - 0.5489848787403228_dp * (  &
        ATan( 0.6354705274484186_dp * (2 * aMu - 0.7677302261302864_dp) )   -  &
        ATan( 0.6354705274484186_dp * (2 * a0  - 0.7677302261302864_dp) )  ) - &
        0.23200408345167492_dp * (  &
        ATan( 0.8864954614918331_dp * (2 * aMu + 1.542321166011278_dp) )    -  &
        ATan( 0.8864954614918331_dp * (2 * a0  + 1.542321166011278_dp) )  ) + &
        0.04040335124693508_dp * &
        log( (0.7664358581916029_dp - 0.7677302261302864_dp * aMu + aMu**2)/&
        (0.7664358581916029_dp - 0.7677302261302864_dp * a0 + a0**2) ) + &
        0.30357668704845114_dp * &
        Log( (0.9651061041938452_dp + 1.542321166011278_dp * aMu + aMu**2)/&
        (0.9651061041938452_dp + 1.542321166011278_dp * a0 + a0**2) )
      case (2)
        rootComplex = rootComplex - 0.5088583673947228_dp * (  &
        ATan( 0.5807020947660965_dp * (2 * aMu - 0.8268729922000535_dp) )   -  &
        ATan( 0.5807020947660965_dp * (2 * a0  - 0.8268729922000535_dp) )  ) - &
        0.1812488725394936_dp * (  &
        ATan( 0.8171726383562911_dp * (2 * aMu + 1.7419567282921797_dp) )    -  &
        ATan( 0.8171726383562911_dp * (2 * a0  + 1.7419567282921797_dp) )  ) + &
        0.03843876906477947_dp * &
        log( (0.9122966896897218_dp - 0.8268729922000535_dp * aMu + aMu**2)/&
        (0.9122966896897218_dp - 0.8268729922000535_dp * a0 + a0**2) ) + &
        0.2771270663760476_dp * &
        Log( (1.1329830828367182_dp + 1.7419567282921797_dp * aMu + aMu**2)/&
        (1.1329830828367182_dp + 1.7419567282921797_dp * a0 + a0**2) )
      case (3)
        rootComplex = rootComplex - 0.4648829323834203_dp * (  &
        ATan( 0.5186296543090279_dp * (2 * aMu - 0.9045061008589117_dp) )   -  &
        ATan( 0.5186296543090279_dp * (2 * a0  - 0.9045061008589117_dp) )  ) - &
        0.10665133580718698_dp * (  &
        ATan( 0.7478910571108373_dp * (2 * aMu + 2.069827558808925_dp) )    -  &
        ATan( 0.7478910571108373_dp * (2 * a0  + 2.069827558808925_dp) )  ) + &
        0.03534462037703361_dp * &
        log( (1.1339812941648535_dp - 0.9045061008589117_dp * aMu + aMu**2)/&
        (1.1339812941648535_dp - 0.9045061008589117_dp * a0 + a0**2) ) + &
        0.24759750067522504_dp * &
        Log( (1.5180010453345267_dp + 2.069827558808925_dp * aMu + aMu**2)/&
        (1.5180010453345267_dp + 2.069827558808925_dp * a0 + a0**2) )
      case (4)
        rootComplex = rootComplex - 0.4162619666563603_dp * (  &
        ATan( 0.4467310336890205_dp * (2 * aMu - 1.0104662655416288_dp) )   -  &
        ATan( 0.4467310336890205_dp * (2 * a0  - 1.0104662655416288_dp) )  ) + &
        0.03691444757280972_dp * (  &
        ATan( 0.7190586675946767_dp * (2 * aMu + 2.742110429973893_dp) )    -  &
        ATan( 0.7190586675946767_dp * (2 * a0  + 2.742110429973893_dp) )  ) + &
        0.029930532499675307_dp * &
        log( (1.507962493506846_dp - 1.0104662655416288_dp * aMu + aMu**2)/&
        (1.507962493506846_dp - 1.0104662655416288_dp * a0 + a0**2) ) + &
        0.21516807986184325_dp * &
        Log( (2.3633089675134706_dp + 2.742110429973893_dp * aMu + aMu**2)/&
        (2.3633089675134706_dp + 2.742110429973893_dp * a0 + a0**2) )
      case (5)
        rootComplex = rootComplex - 0.36091735121516355_dp * (  &
        ATan( 0.3626874145171533_dp * (2 * aMu - 1.163068663777_dp) )   -  &
        ATan( 0.3626874145171533_dp * (2 * a0  - 1.163068663777_dp) )  ) + &
        0.3943399946498003_dp * log( (aMu + 1.5721862037478178_dp)/&
        (a0 + 1.5721862037478178_dp) ) - &
        0.027723547491953493_dp * log( (aMu + 3.511042056691402_dp)/&
        (a0 + 3.511042056691402_dp) ) + 0.017365400319553415_dp * &
        Log( (2.2387135279454564_dp - 1.163068663777_dp * aMu + aMu**2)/&
        (2.2387135279454564_dp - 1.163068663777_dp * a0 + a0**2) )
      case (6)
        rootComplex = rootComplex - 0.2706517070577726_dp * (  &
        ATan( 0.27659246828604006_dp * (2 * aMu - 1.6112137089692589_dp) )   -  &
        ATan( 0.27659246828604006_dp * (2 * a0  - 1.6112137089692589_dp) )  ) + &
        0.33463363709028904_dp * log( (aMu + 1.4322403491695925_dp)/&
        (a0 + 1.4322403491695925_dp) ) - &
        3.242619128589197e-8_dp * log( (aMu + 114.63875643838061_dp)/&
        (a0 + 114.63875643838061_dp) ) - 0.019530069461003236_dp * &
        Log( (3.916831300483894_dp - 1.6112137089692589_dp * aMu + aMu**2)/&
        (3.916831300483894_dp - 1.6112137089692589_dp * a0 + a0**2) )
      case default
        rootComplex = 0
      end select
    end if

   end function rootComplex

!ccccccccccccccc

   pure real (dp) function Gfun(self, order, t)
     class (AnomDim)    , intent(in) :: self
     real (dp)          , intent(in) :: t
     integer            , intent(in) :: order
     real (dp), dimension( min(3,order - 1) ) :: listT
     integer                         :: i, ord

     ord = min(3,order - 1)

     listT = 1 / powList(t,ord) / [ (i, i = 1, ord) ]

     if (order >= 0) Gfun = t - dot_product( self%bHat(2:ord + 1), listT )

     Gfun = exp(Gfun);     if (order > 0) Gfun = (-t)**self%bHat(1) * Gfun

   end function Gfun

!ccccccccccccccc

   pure real (dp) function gamma4(self)
    class (AnomDim), intent(in) :: self
    gamma4 = self%G4
   end

!ccccccccccccccc

   pure integer function numFlav(self)
    class (AnomDim), intent(in) :: self
    numFlav = self%nf
   end

!ccccccccccccccc

   pure real (dp) function wTildeReal(self, order, gam, a0, a1)
    class(AnomDim)               , intent(in) :: self
    integer                      , intent(in) :: order
    real (dp)                    , intent(in) :: a0, a1
    real (dp), dimension(0:order), intent(in) :: gam
    real (dp)                                 :: a0Pi, a1Pi
    integer                                   :: i

    wTildeReal = 0; if (order < 0) return

    if (order >= 0) wTildeReal = gam(0) * log(a1/a0)

    if (order > 0) then
      a0Pi = a0/fourPi; a1Pi = a1/fourPi
    end if

    do i = 1, order

      wTildeReal = wTildeReal + (a1Pi**i - a0Pi**i)/i * &
      dot_product( gam(:i), self%cCoef(i:0:-1) )

    end do

    wTildeReal = - wTildeReal/self%beta(0)

   end function wTildeReal

!ccccccccccccccc

   pure real (dp) function wTildeComplex(self, order, gam, a0, a1)
    class(AnomDim)               , intent(in) :: self
    integer                      , intent(in) :: order
    real (dp)                    , intent(in) :: a1
    complex (dp)                 , intent(in) :: a0
    real (dp), dimension(0:order), intent(in) :: gam
    complex (dp)                              :: wTilde, a0Pi
    real (dp)                                 :: a1Pi
    integer                                   :: i

    wTilde = 0; if (order < 0) return

    if (order >= 0) wTilde = gam(0) * log(a1/a0)

    if (order == 0) then
      a0Pi = a0/fourPi; a1Pi = a1/fourPi
    end if

    do i = 1, order

      wTilde = wTilde + (a1Pi**i - a0Pi**i)/i * &
      dot_product( gam(:i), self%cCoef(i:0:-1) )

    end do

    wTildeComplex = - realpart(wTilde)/self%beta(0)

   end function wTildeComplex

!ccccccccccccccc

   pure real (dp) function wTildeHm(self, order, gamma, a0, a1)
    class (AnomDim)              , intent(in) :: self
    integer                      , intent(in) :: order
    real (dp)                    , intent(in) :: a0, a1
    real (dp), dimension(0:3,0:3), intent(in) :: gamma
    real (dp), dimension(0:3)                 :: gammaAux

    gammaAux = gamma(:,0)

    if (order > 0) gammaAux(0) = gammaAux(0) + gamma(0,1)

    if (order > 1) then
     gammaAux(0) = gammaAux(0) + gamma(0,2); gammaAux(1) = gammaAux(1) + gamma(1,1)
   end if

    if (order > 2) then
      gammaAux(0) = gammaAux(0) + gamma(0,3); gammaAux(1) = gammaAux(1) + gamma(1,2)
      gammaAux(2) = gammaAux(2) + gamma(2,1)
    end if

   wTildeHm = self%wTilde(order, gammaAux, a0, a1)

   end function wTildeHm

!ccccccccccccccc

   pure real (dp) function kTildeReal(self, order, gam, a0, a1)
    class(AnomDim)               , intent(in) :: self
    integer                      , intent(in) :: order
    real (dp)                    , intent(in) :: a0, a1
    real (dp), dimension(0:order), intent(in) :: gam
    real (dp), dimension(0:order)             :: d
    real (dp)                                 :: lg, a0Pi, a1Pi, a0j, a1j, cd, suma
    integer                                   :: i, j

    kTildeReal = 0; if (order < 0) return

    d(0) = gam(0)

    do i = 1, order
      d(i) = dot_product( gam(:i), self%cCoef(i:0:-1) )/i
    end do

    lg = log(a1/a0); a0Pi = a0/fourPi; a1Pi = a1/fourPi

    kTildeReal = d(0) * (1/a1Pi - 1/a0Pi + lg/a0Pi)

    if (order > 0) kTildeReal = kTildeReal + d(1) * (a1/a0 - 1 - lg) + &
    d(0) * self%cCoef(1) * lg**2/2

    do j = 2, order

      a1j = a1Pi**(j-1); a0j = a0Pi**(j - 1); cd = d(0) * self%cCoef(j)/(j-1)

      suma = 0

      do i = 2, j - 1
        suma = suma + d(j - i) * self%cCoef(i) * &
        ( a1j - a1Pi**(j - i) * a0Pi**(i - 1) )/(i - 1)
      end do

      kTildeReal = kTildeReal + d(j) * a1j * (a1/a0 - 1) + suma              &
      + ( a1j - a0j ) * ( cd - sum( self%cCoef(:j-1) * d(j:1:-1) ) )/(j - 1) &
      + ( d(j-1) * self%cCoef(1) * a1j - cd * a0j ) * lg

    end do

    kTildeReal = kTildeReal/self%beta(0)**2/2

   end function kTildeReal

!ccccccccccccccc

   pure real (dp) function kTildeComplex(self, order, gam, a0, a1)
    class(AnomDim)               , intent(in) :: self
    integer                      , intent(in) :: order
    real (dp)                    , intent(in) :: a1
    complex (dp)                 , intent(in) :: a0
    real (dp), dimension(0:order), intent(in) :: gam
    real (dp), dimension(0:order)             :: d
    complex (dp)                              :: a0Pi, a0j, lg, kTilde, suma
    real (dp)                                 :: cd, a1Pi, a1j
    integer                                   :: i, j

    kTilde = 0; if (order < 0) return

    d(0) = gam(0)

    do i = 1, order
      d(i) = dot_product( gam(:i), self%cCoef(i:0:-1) )/i
    end do

    lg = log(a1/a0); a0Pi = a0/fourPi; a1Pi = a1/fourPi

    kTilde = d(0) * (1/a1Pi - 1/a0Pi + lg/a0Pi)

    if (order > 0) kTilde = kTilde + d(1) * (a1/a0 - 1 - lg) + &
    d(0) * self%cCoef(1) * lg**2/2

    do j = 2, order

      a1j = a1Pi**(j-1); a0j = a0Pi**(j - 1); cd = d(0) * self%cCoef(j)/(j-1)

      suma = 0

      do i = 2, j - 1
        suma = suma + d(j - i) * self%cCoef(i) * &
        ( a1j - a1Pi**(j - i) * a0Pi**(i - 1) )/(i - 1)
      end do

      kTilde = kTilde + d(j) * a1j * (a1/a0 - 1) + suma                      &
      + ( a1j - a0j ) * ( cd - sum( self%cCoef(:j-1) * d(j:1:-1) ) )/(j - 1) &
      + ( d(j-1) * self%cCoef(1) * a1j - cd * a0j ) * lg

    end do

    kTildeComplex = realpart(kTilde)/self%beta(0)**2/2

   end function kTildeComplex

!ccccccccccccccc

   pure subroutine expandAlpha(self, coef)
    class (AnomDim)          , intent(in   ) :: self
    real (dp), dimension(:,:), intent(inout) :: coef

    call AlphaExpand(self%beta, coef)

   end subroutine expandAlpha

!ccccccccccccccc

   pure subroutine AlphaExpand(beta, coef)
    real (dp), dimension(0:4), intent(in)    :: beta
    real (dp), dimension(:,:), intent(inout) :: coef
    integer                                  :: n, j, i

    coef(2:,:) = 0

    do n = 2, size(coef,2)
      do j = 1, n - 1

        do i = j, n - 1
          coef(j + 1,n) = coef(j + 1,n) + i * beta(n - i - 1) * &
          coef(j,i)/2**( 2 * (n - i) )
        end do

      coef(j + 1,n) = 2 * coef(j + 1,n)/j

      end do
    end do

   end subroutine AlphaExpand

!ccccccccccccccc

   pure subroutine wTildeExpand(self, gamma, coef)
    class (AnomDim), intent(in )             :: self
    real (dp), intent(in ), dimension(0:3)   :: gamma
    real (dp), intent(out), dimension(0:4,3) :: coef
    integer                                  :: n, i, j

    coef = 0; coef(1,:) = 2 * gamma(:2)/powList(4,3)

    do n = 2, 3
      do j = 2, n

        do i = j - 1, n - 1
          coef(j,n) = coef(j,n) + i * self%beta(n - i - 1) * &
          coef(j-1,i)/2**( 2 * (n - i) )
        end do

      coef(j,n) = 2 * coef(j,n)/j

      end do
    end do

   end subroutine wTildeExpand

!ccccccccccccccc

   pure subroutine kTildeExpand(self, gamma, coef)
    class (AnomDim)            , intent(in ) :: self
    real (dp), dimension(0:3)  , intent(in ) :: gamma
    real (dp), dimension(0:4,3), intent(out) :: coef
    integer                                  :: n, i, j

    coef = 0; coef(2,:) = gamma(:2)/powList(4,3)

    do n = 2, 3
      do j = 3, n + 1

        do i = j - 2, n - 1
          coef(j,n) = coef(j,n) + i * self%beta(n - i - 1) * &
          coef(j-1,i)/2**( 2 * (n - i) )
        end do

      coef(j,n) = 2 * coef(j,n)/j

      end do
    end do

   end subroutine kTildeExpand

!ccccccccccccccc

   real (dp) function DeltaMu(self, order, R, alpha0, alpha1)
     class (AnomDim), intent(in) :: self
     integer        , intent(in) :: order
     real (dp)      , intent(in) :: R, alpha0, alpha1

     DeltaMu = 2 * ExpEuler * R * self%wTilde(order - 1, self%cusp, alpha0, alpha1)

   end

!ccccccccccccccc

  pure real (dp) function betaQED(self)
    class (AnomDim), intent(in) :: self
    betaQED = self%beta0QED
  end function betaQED

!ccccccccccccccc

  pure function betaQCD(self, str) result(bet)
    class (AnomDim)    , intent(in) :: self
    character (len = *), intent(in) :: str
    real (dp)      , dimension(0:4) :: bet

    bet = 0

    if ( str( :4) == 'beta'            ) bet     = self%beta
    if ( str( :4) == 'cusp'            ) bet(:3) = self%cusp
    if ( str( :4) == 'mass'            ) bet     = self%gammaMass
    if ( str( :4) == 'hard'            ) bet(:3) = self%gammaHard
    if ( str( :3) == 'jet'             ) bet(:3) = self%gammaJet
    if ( str( :4) == 'soft'            ) bet(:3) = self%gammaSoft
    if ( str( :4) == 'bJet'            ) bet(:3) = self%gammaB
    if ( str( :2) == 'Hm'              ) bet(:3) = self%gammaHm(:,0)
    if ( str( :4) == 'bHat'            ) bet     = self%bHat
    if ( str( :5) == 'bCoef'           ) bet     = self%bCoef
    if ( str( :5) == 'cCoef'           ) bet     = self%cCoef
    if ( str( :2) == 'gl'              ) bet(:3) = self%gl
    if ( str( :6) == 'gTilde'          ) bet(:3) = self%gTilde
    if ( str( :9) == 'MSRpdelta'       ) bet(1:) = self%MSRDelta('MSRp')
    if ( str( :9) == 'MSRndelta'       ) bet(1:) = self%MSRDelta('MSRn')
    if ( str( :9) == 'sCoefMSRp'       ) bet(1:) = self%sCoefMSRp
    if ( str(:12) == 'sCoefMSRInc1'    ) bet(1:) = self%sCoefMSRInc1
    if ( str(:12) == 'sCoefMSRInc2'    ) bet(1:) = self%sCoefMSRInc2
    if ( str(:12) == 'sCoefMSRInc3'    ) bet(1:) = self%sCoefMSRInc3
    if ( str( :9) == 'sCoefMSRn'       ) bet(1:) = self%sCoefMSRn
    if ( str( :7) == 'sCoefRS'         ) bet(1:) = self%sCoefRS
    if ( str( :8) == 'betaList'        ) bet(1:) = self%betaList
    if ( str( :9) == 'gammaRInc1'      ) bet(1:) = self%gammaRInc1
    if ( str(:10) == 'gammaRInc2'      ) bet(1:) = self%gammaRInc2
    if ( str(:10) == 'gammaRInc3'      ) bet(1:) = self%gammaRInc3
    if ( str( :7) == 'gammaRp'         ) bet(1:) = self%gammaRp
    if ( str( :7) == 'gammaRn'         ) bet(1:) = self%gammaRn
    if ( str( :7) == 'gammaRS'         ) bet(1:) = self%gammaRS
    if ( str( :6) == 'betQED'          ) bet(1:) = self%beta0QED
    if ( str( :7) == 'RSdelta'         ) bet(1:) = self%aRS
    if ( str( :8) == 'RSndelta'        ) bet(1:) = self%aRSn
    if ( str( :8) == 'RSpdelta'        ) bet(1:) = self%aRSp

  end function betaQCD

!ccccccccccccccc

  real (dp) function P12Generic(self, order, a, lambda)
     class (AnomDim)        , intent(in) :: self
     real (dp)              , intent(in) :: lambda
     integer                , intent(in) :: order
     real (dp), dimension(4), intent(in) :: a
     real (dp), dimension(4)             :: sCoef
     integer                             :: i

     sCoef = self%sCoefRecursive(a)
     sCoef = self%sCoefGeneric(sCoef, lambda)

     P12Generic = 0

     do i = 0, min(order, 3)
       P12Generic = P12Generic + sCoef(i + 1)/gamma( 1 + self%bHat(1) + i )
     end do

  end function P12Generic

!ccccccccccccccc

  real (dp) function N12Generic(self, order, a, lambda)
     class (AnomDim)        , intent(in) :: self
     real (dp)              , intent(in) :: lambda
     real (dp), dimension(4), intent(in) :: a
     integer                , intent(in) :: order

     N12Generic = self%beta(0) * gamma( 1 + self%bHat(1) )/2/Pi * self%P12Generic(order, a, lambda)

  end function N12Generic

!ccccccccccccccc

  real (dp) function P12(self, order, type, lambda)
     class (AnomDim)      , intent(in) :: self
     real (dp)            , intent(in) :: lambda
     character (len = *)  , intent(in) :: type
     integer              , intent(in) :: order
     real (dp)          , dimension(4) :: sCoef
     integer                           :: i

     sCoef = self%sCoefLambda(type, lambda)

     P12 = 0

     do i = 0, min(order, 3)
       P12 = P12 + sCoef(i + 1)/gamma( 1 + self%bHat(1) + i )
     end do

  end function P12

!ccccccccccccccc

  real (dp) function N12(self, order, type, lambda)
     class (AnomDim)      , intent(in) :: self
     real (dp)            , intent(in) :: lambda
     character (len = *)  , intent(in) :: type
     integer              , intent(in) :: order

     N12 = self%beta(0) * gamma( 1 + self%bHat(1) )/2/Pi * self%P12(order, type, lambda)

  end function N12

!ccccccccccccccc

  real (dp) function N12Ratio(self, order, type, lambda)
     class (AnomDim)            , intent(in) :: self
     real (dp)                  , intent(in) :: lambda
     character (len = *)        , intent(in) :: type
     integer                    , intent(in) :: order
     real (dp)              , dimension(0:3) :: sCoef
     real (dp), dimension( 0:min(order, 3) ) :: poch
     real (dp)                               :: b1
     integer                                 :: k, n

     N12Ratio = 0; n = min(order, 3); b1 = self%bHat(1)

     sCoef = self%sCoefLambda(type, lambda)

     do k = 0, n

       poch(:n - k) = PochHammerList(1 + b1 + k, n - k)

       N12Ratio = N12Ratio + sCoef(k) * sum( self%gl(:n-k) * poch(n-k:0:-1)  )

     end do

     poch = PochHammerList(1 + b1, n)

     N12Ratio =  N12Ratio * self%beta(0)/2/Pi/sum( self%gl(:n) * poch(n:0:-1) )

  end function N12Ratio

!ccccccccccccccc

  real (dp) function N12Residue(self, order, type, lambda)
     class (AnomDim)            , intent(in) :: self
     real (dp)                  , intent(in) :: lambda
     character (len = *)        , intent(in) :: type
     integer                    , intent(in) :: order
     real (dp)              , dimension(0:3) :: sCoef
     real (dp), dimension( 0:min(order, 3) ) :: poch, listfact
     real (dp), dimension(   min(order, 4) ) :: a
     real (dp)                               :: b1
     integer                                 :: k, n, i

     N12Residue = 0; n = min(order, 3); b1 = self%bHat(1); a = 0

     sCoef = self%sCoefLambda(type, lambda)

     do i = 1, n + 1
       do k = 0, i - 1

       poch(:n - 1 - k) = PochHammerList(1 + b1 + k, i - 1 - k)

       a(i) = a(i) + sCoef(k) * sum( self%gl(:i - 1 - k) * poch(i - 1 - k:0:-1)  )

       end do
     end do

     poch = PochHammerList(-b1, n); listfact = factList(n)

     do k = 0, n
       N12Residue =  N12Residue + a(k + 1) * poch(n - k)/listfact(k)/listfact(n - k)
     end do

     N12Residue = N12Residue * self%beta(0)/2/Pi

  end function N12Residue


  !ccccccccccccccc

  real (dp) function N12RS(self, order, m, type, lambda)
     class (AnomDim)            , intent(in) :: self
     real (dp)                  , intent(in) :: lambda
     character (len = *)        , intent(in) :: type
     integer                    , intent(in) :: order
     real (dp)              , dimension(0:3) :: sCoef
     real (dp), dimension(0:6)               :: poch
     integer  , dimension(0:4)               :: signlist
     real (dp)                               :: b1, suma, a
     integer                                 :: k, n, m, l

     N12RS = 0; if (m < 4) return
     n = min(order, 3); b1 = self%bHat(1); a = 0

     sCoef = self%sCoefLambda(type, lambda)

     do k = 0, m - 1

       poch(:n - 1 - k) = PochHammerList(1 + b1 + k, m - 1 - k)

       a = a + sCoef(k) * sum( self%gl(:m - 1 - k) * poch(m - 1 - k:0:-1)  )

     end do

     signlist = PowList0(-1,4)

     do k = 0, 3
       poch(0:4 + k - m) = PochHammerList( - b1 - k, 4 + k - n )
       l = max(n - k, 0)
       suma = sum( self%gl(l:) * signlist(l:)/poch(1+k+l-n:4+k-m) )
       N12RS = N12RS + sCoef(k) * suma * (-1)**k
     end do

     N12RS = a - N12RS * self%beta2List(m) * signlist(m)

  end function N12RS

!ccccccccccccccc

  function sCoefLambda(self, type, lambda) result(res)
     class (AnomDim)      , intent(in) :: self
     real (dp)            , intent(in) :: lambda
     character (len = *)  , intent(in) :: type
     real (dp)          , dimension(4) :: res, sCoef
     real (dp)                         :: lg

     sCoef = 0; lg = log(lambda)

     if ( type(:4) == 'Inc1' ) then
       sCoef = self%sCoefMSRInc1
     else if ( type(:4) == 'Inc2' ) then
       sCoef = self%sCoefMSRInc2
     else if ( type(:4) == 'Inc3' ) then
       sCoef = self%sCoefMSRInc3
     else if ( type(:7) == 'Natural' ) then
       sCoef = self%sCoefMSRn
     else if ( type(:9) == 'Practical' ) then
       sCoef = self%sCoefMSRp
     end if

     res = self%sCoefGeneric(sCoef, lambda)

  end function sCoefLambda

!ccccccccccccccc

  function sCoefGeneric(self, sCoef, lambda) result(res)
     class (AnomDim)        , intent(in) :: self
     real (dp)              , intent(in) :: lambda
     real (dp), dimension(4), intent(in) :: sCoef
     real (dp)            , dimension(4) :: res
     real (dp)                           :: lg

     res = sCoef; lg = log(lambda)

     res(2) = res(2) - lg * sCoef(1)

     res(3) = res(3) - 2 * lg * sCoef(2) + sCoef(1) * (  lg**2 - ( 2 * self%bHat(1) &
     + self%bHat(2) ) * lg  )

     res(4) = res(4) - 3 * lg * sCoef(3) + sCoef(2) * (  3 * lg**2 - ( 3 * self%bHat(1) &
     + self%bHat(2) ) * lg  ) + sCoef(1) * (  (self%bHat(3) - 3 * self%bHat(1)**2 +&
     3 * self%bHat(2) - self%bHat(1) * self%bHat(2)) * lg + ( 9 * self%bHat(1)/2 &
     + 2 * self%bHat(2) ) * lg**2 - lg**3  )

     res = lambda * res

  end function sCoefGeneric

!ccccccccccccccc

  pure function sCoef(self, gamma) result(s)
    class (AnomDim)        , intent(in) :: self
    real (dp), dimension(:), intent(in) :: gamma
    real (dp), dimension( size(gamma) ) :: s
    integer                             :: i, j, k
    real (dp)                           :: sum1, sum2

    if ( size(gamma) > 0 ) s = gamma

    do i = 1, size(gamma)

      sum1 = 0

      do j = 0, i - 2

        sum2 = 0

        do k = 1, i - j - 1
          sum2 = sum2 + (-1)**(i - j - k) * self%gTilde(k) * self%bHat(i - j - k)
        end do

        sum1 = sum1 + gamma(j+1) * &
        ( self%gTilde(i - j) + (-1)**(i - j) * self%bHat(i - j) + sum2 )

      end do

      s(i + 1) = s(i + 1) + sum1 - sum( self%bHat(1:2) ) * gamma(i)

    end do

  end function sCoef

!ccccccccccccccc

  pure function sCoefRecursive(self, a) result(s)
    class (AnomDim)        , intent(in) :: self
    real (dp), dimension(:), intent(in) :: a
    real (dp), dimension( 0:size(a)-1 ) :: s
    real (dp), dimension(size(a) )      :: tilA
    integer                             :: k, n
    real (dp)                           :: suma

    tilA = a * self%betaList( :size(a) ) * powList(4,size(a))

    do k = 0, size(a) - 1

      suma = 0

      do n = 0, k - 1

        suma = suma + s(n) * dot_product( self%gl(k - n:0:-1), &
        PochHammerList(1 + self%bHat(1) + n, k - n) )

      end do

      s(k) = tilA(k+1) - suma

    end do

  end function sCoefRecursive

!ccccccccccccccc

  pure function PochHammerList(a,n) result(poch)
    real (dp), intent(in)     :: a
    integer  , intent(in)     :: n
    real (dp), dimension(0:n) :: poch
    integer                   :: j

    poch(0)  = 1
    poch(1:) = [ (a + j, j = 0, n - 1) ]

    do j = n, 1, -1
      poch(j) = Product( poch(1:j) )
    end do

  end function PochHammerList

!ccccccccccccccc

  pure function sCoefHadron(self, gamma, r) result(s)
    class (AnomDim), intent(in)         :: self
    real (dp), intent(in)               :: r
    real (dp), intent(in), dimension(3) :: gamma
    real (dp)            , dimension(3) :: s
    real (dp)                           :: gammaHadron

    gammaHadron = 6 * log(1 - r**2)/self%beta(0);  s = self%sCoef(gamma)

    s(2:3) = s(2:3) + gammaHadron * gamma(:2)   ;  s(3) = s(3) + gammaHadron * gamma(1)

  end function sCoefHadron

!ccccccccccccccc

  real (dp) function DeltaR(self, sCoef, order, a0, a1)
    integer                    , intent(in) :: order
    real (dp)                  , intent(in) :: a0, a1
    class (AnomDim)            , intent(in) :: self
    real (dp), dimension(order), intent(in) :: sCoef
    integer                                 :: i
    real (dp)                               :: t0, t1, Gcorre

    t0 = - 2 * pi/self%beta(0);  t1 = t0/a1;  t0 = t0/a0;  DeltaR = 0

    do i = 1, order
      if (  abs( sCoef(i) ) < d1mach(1)  ) cycle
      Gcorre = inteCorre( - self%bhat(1) - i, t0, t1 )
      DeltaR = DeltaR + sCoef(i) * Gcorre
    end do

  end function DeltaR

!ccccccccccccccc

  real (dp) function DeltaRMass(self, order, lambda, mass, a0, a1)
    integer        , intent(in) :: order
    real (dp)      , intent(in) :: a0, a1, mass, lambda
    class (AnomDim), intent(in) :: self
    integer                     :: neval, ier
    real (dp)                   :: a, t0, t1, abserr

    t0 = - 2 * pi/self%beta(0);  t1 = t0/a1;  t0 = t0/a0;   DeltaRMass = 0
    a  = - self%bhat(1) - 2

    call qags( correMass, t0, t1, prec, prec, DeltaRMass, abserr, neval, ier )

    DeltaRMass = DeltaRMass * self%betaList(2)

    contains

!ccccccccccccccc

    pure real (dp) function correMass(t)
      real (dp), intent(in) :: t
      real (dp)             :: x

      x = lambda/self%Gfun(order, t)/mass
      correMass = Exp(-t) * (-t)**a * GammaRMass(x)

    end function correMass

  end function DeltaRMass

!ccccccccccccccc

  real (dp) function DeltaRHadron(self, sCoef, r, order, a0, a1)
    integer                , intent(in) :: order
    real (dp)              , intent(in) :: a0, a1, r
    class (AnomDim)        , intent(in) :: self
    real (dp), dimension(3), intent(in) :: sCoef
    integer                             :: i
    real (dp)                           :: t0, t1, gammaHadron, Gcorre

    gammaHadron = 6 * log(1 - r**2)/self%beta(0)

    t0 = - 2 * pi/self%beta(0)/a0;  t1 = - 2 * pi/self%beta(0)/a1

    DeltaRHadron = 0

    do i = 2, order

      Gcorre = inteCorre( - self%bhat(1) + gammaHadron - i, t0, t1 )
      DeltaRHadron = DeltaRHadron + sCoef(i) * Gcorre

    end do

  end function DeltaRHadron

!ccccccccccccccc

  pure function bHQETgamma(self) result(tab)
    class (AnomDim), intent(in)   :: self
    real (dp), dimension(0:3,0:3) :: tab

    tab = self%gammaHm

  end function bHQETgamma

!ccccccccccccccc

  pure function MatchingAlpha(self) result(tab)
    class (AnomDim), intent(in) :: self
    real (dp)    , dimension(5) :: tab

    tab = self%alphaMatch

  end function MatchingAlpha

!ccccccccccccccc

  pure function MatchingAlphaLog(self, direction) result(tab)
    class (AnomDim)    , intent(in) :: self
    character (len = *), intent(in) :: direction
    real (dp)    , dimension(0:4,5) :: tab

    if ( direction(:4) == 'down' ) then
      tab = self%alphaMatchLog
    else if ( direction(:2) == 'up' ) then
      tab = self%AlphaMatchInvLog
    end if

end function MatchingAlphaLog

!ccccccccccccccc

  pure function MatchingAlphaUp(self) result(tab)
    class (AnomDim), intent(in)     :: self
    real (dp)  , dimension(5)       :: tab

    tab = self%alphaMatchUp

  end function MatchingAlphaUp

!ccccccccccccccc

  pure function alphaMatching(self, nf) result(tab)
    class (AnomDim), intent(in)     :: self
    integer        , intent(in)     :: nf
    real (dp)  , dimension(5)       :: tab

    tab = alphaMatch(self%str, nf)

  end function alphaMatching

!ccccccccccccccc

  pure function alphaMatch(str, nf) result(tab)
    character (len = *), intent(in) :: str
    integer        , intent(in)     :: nf
    real (dp)  , dimension(5)       :: tab, tab2
    real (dp)  , dimension(4)       :: b
    real (dp)  , dimension(0:4,0:4) :: d
    real (dp)  , dimension(0:4,5)   :: c
    integer                         :: i, j, k

    tab = 0; tab2 = 0; tab2(1) = 1
    tab2(3:) = - [7._dp/24, 5.586361025786356_dp - 0.26247081195432964_dp * nf,&
    95.80685705811597_dp - 10.171370332526523_dp * nf + 0.23954195789610408_dp * nf**2]

    if ( str(:5) == 'MSbar' ) then

      b = MSbarDelta(nf - 1, 1); d = 0; d(0,0) = 1

      do i = 1, 4

        do j = 1, i - 1
          d(1,i) = d(1,i) + j * d(1,j) * b(i - j)
        end do

        d(1,i) = b(i) - d(1,i)/i

      end do

      c = alphaMatchingLog(tab2, nf, nf - 1)

      do i = 1, 3
        do j = i, 3
          d(i + 1,j) = sum( d(1,1:j - i) * d(i,j - 1:i:-1) )
        end do
      end do

      do k = 1, 5
        do i = 1, k
          do j = 0, min(i - 1,k - i)
            tab(k) = tab(k) + (-1)**j * c(j,i) * d(j, k - i)
          end do
        end do
      end do

    else if ( str(:4) == 'pole' ) then

      tab = tab2

    end if

  end function alphaMatch

!ccccccccccccccc

  pure function alphaMatchingLog(d, nf, nl) result(tab)
    integer                , intent(in) :: nf, nl
    real (dp), dimension(5), intent(in) :: d
    real (dp), dimension(0:4,5)         :: tab  ! (log, alpha)
    real (dp), dimension(0:4,0:5)       :: b, c ! (log, alpha)
    real (dp), dimension(0:4,5,5)       :: ePow ! (log,alpha,power)
    integer                             :: n, i, j, k, l, m

    tab = 0; tab(0,1) = 1; ePow = 0

    b = 0; c = 0; c(0,1) = 1 ;  b(0,1:) = d
    call AlphaExpand( betaFun(nf), b(:,1:) )
    call AlphaExpand( betaFun(nl), c(:,1:) )

    do n = 2, 5

      ePow(:,n - 1,1) = tab(:,n - 1)

      do i = 2, n
        do l = 0, n - i

          do m = 1, n + 1 - i
            j = max(0,l + i - 1 + m - n); k = min(l, m - 1)
            ePow(l,n,i) = ePow(l,n,i) + sum( tab(j:k,m) * ePow(l-j:l-k:-1,n-m,i-1) )
          end do

        end do
      end do

      do l = 0, n - 1

        tab(l,n) = b(l,n)

        do i = 2, n
          j = max(0,l + i - n); k = min(l,i - 1)
          tab(l,n) = tab(l,n) - sum( c(j:k,i) * ePow(l - j:l - k:-1,n,i) )
        end do

      end do

    end do

  end function alphaMatchingLog

!ccccccccccccccc pole - MSR mass with mu = R

  pure function MSRDelta(self, type) result(coef)
    class (AnomDim)    , intent(in) :: self
    character (len = *), intent(in) :: type
    real (dp), dimension(4)         :: coef
    real (dp), dimension(5)         :: b

    if ( type(:4) == 'MSRn' ) then

      coef = MSbarDelta(self%nf, 0, self%err)

    else if ( type(:4) == 'MSRp' ) then

      coef = MSbarDelta(self%nf, 1, self%err)
      b = getInverse( self%alphaMatching(self%nf + 1) )
      call alphaReExpand( coef, b(:4) )

    else if ( type(:2) == 'RS' ) then

      coef = self%aRS

    else if ( type(:5) == 'charm' ) then

      coef = MSbarDelta(self%nf - 1, 2, self%err)
      b = getInverse( self%alphaMatching(self%nf + 1) )
      call alphaReExpand( coef, b(:4) )

    end if

  end function MSRDelta

!ccccccccccccccc

  pure function GammaRComputer(self, gamma) result(gammaMSR)
    class (AnomDim)        , intent(in) :: self
    real (dp), dimension(:), intent(in) :: gamma
    real (dp), dimension( size(gamma) ) :: gammaMSR
    integer  , dimension( size(gamma) ) :: List4
    integer                             :: i, j

    List4 = powList( 4, size(gamma) ) ; gammaMSR = List4 * gamma

    do i = 1, size(gamma)
      do j = 1, i - 1
        gammaMSR(i) = gammaMSR(i) - 2 * List4(j) * j * gamma(j) * self%beta(i - 1 - j)
      enddo
    end do

  end function GammaRComputer

!ccccccccccccccc ! Expand series in alpha when alpha is another series in alpha' (e.g. flavor matching)

  pure subroutine alphaReExpand(a, b)
    real (dp), dimension(:), intent(inout)   :: a
    real (dp), dimension(:), intent(in)      :: b
    real (dp), dimension( size(b), size(b) ) :: c
    integer                                  :: n, j

    c(1,:) = b

    do n = 1, size(a) - 1

      do j = n + 1, size(a)
        c(n+1,j) = sum( b(:j-n) * c(n,j-1:n:-1) )
      end do

      a(n+1) = sum( a(:n+1) * c(:n+1,n+1) )

    end do

  end

!ccccccccccccccc

  real (dp) function inteCorre(a, x0, x1)
   real (dp), intent(in) :: x0, x1, a
   integer               :: i
   real (dp)             :: den, fac, corr

     den = a + 1; fac = 1; inteCorre = ( (-x0)**den - (-x1)**den )/den; i = 0

     do
       i = i + 1; den = den + 1; fac = fac * i
       corr = ( (-x0)**den - (-x1)**den )/den/fac
       if ( i .gt. 5 .and. abs(corr) > abs(inteCorre) ) exit
       inteCorre = inteCorre + corr; if ( abs(corr) <= prec ) return
     end do

     inteCorre = dGauss(corre, x0, x1, prec)

   contains

     pure real (dp) function corre(x)
       real (dp), intent(in) :: x
       corre = Exp(-x) * (-x)**a
     end function corre

  end function inteCorre


!ccccccccccccccc

  pure real (dp) function GammaRMass(x)
    real (dp), intent(in) :: x
    real (dp)             :: lx

    lx = Log(x)

    GammaRMass = 27 * gam * ( ep * (2 * x**bet + al * bet) + x * (del * x**bet + del *   &
    al * bet + x * al * bet) ) + 8 * x * (   ep * (9 * (12 - 12 * lx + 12 * lx**2 + Pi2) &
    * x**bet + (56 - 60 * lx + 36 * lx**2 + 3 * Pi2) * al * bet) + x * (  x * ( (12 * lx &
    + 36 * lx**2 + 3 * Pi2 - 4) * x**bet + (56 - 60 * lx + 36 * lx**2 + 3 * Pi2) * al *  &
    bet) + del * ( (52 - 48 * lx + 72 * lx**2 + 6 * Pi2) * x**bet + (56 - 60 * lx + 36 * &
    lx**2 + 3 * Pi2) * al * bet )  )   )

    GammaRMass = 4 * exp(-x**(-bet) * al) * x**(1 - bet) * GammaRMass/81/( ep + x * (del + x) )**2

    GammaRMass = ExpEuler * ( GammaRMass + 1024._dp/81 - 128 * lx/27 - 128 * lx**2/9 )

  end function GammaRMass

!ccccccccccccccc

  pure real(dp) function DeltaMass(L, x)
    real (dp), intent(in) :: L, x
    real (dp)             :: lx

    lx = log(x); lx = 160 * lx/9 - 448._dp/27 - 32 * lx**2/3

    deltaMass = (  32 * L**2/3 + 160 * L/9 + 224._dp/27 + lx + ( exp(-x**(-bet) * al) * &
    (  8 * Pi2/9 - lx + gam/x) )/(1 + ep/x**2 + del/x)  )/12

  end function DeltaMass

!ccccccccccccccc

  pure function PowList0int(alpha, n) result(list)
    integer    , intent(in) :: alpha
    integer    , intent(in) :: n
    integer, dimension(0:n) :: list

    list(0) = 1; if (n > 0) list(1:) = PowList(alpha, n)

  end function PowList0Int

!ccccccccccccccc

  pure function PowList0DP(alpha, n) result(list)
    real (dp)    , intent(in) :: alpha
    integer      , intent(in) :: n
    real (dp), dimension(0:n) :: list

    list(0) = 1; if (n > 0) list(1:) = PowList(alpha, n)

  end function PowList0DP

!ccccccccccccccc

  pure function PowList0Comp(alpha, n) result(list)
    complex (dp)    , intent(in) :: alpha
    integer         , intent(in) :: n
    complex (dp), dimension(0:n) :: list

    list(0) = 1; if (n > 0) list(1:) = PowList(alpha, n)

  end function PowList0Comp

!ccccccccccccccc

  pure function PowListDP(alpha, n) result(list)
    real (dp)  , intent(in) :: alpha
    integer    , intent(in) :: n
    integer                 :: i
    real (dp), dimension(n) :: list
    real (dp)               :: a

    a = 1

    do i = 1, n
      a = a * alpha; list(i) = a
    end do

  end function PowListDP

!ccccccccccccccc

  pure function PowListComp(alpha, n) result(list)
    complex (dp)  , intent(in) :: alpha
    integer       , intent(in) :: n
    integer                    :: i
    complex (dp), dimension(n) :: list
    complex (dp)               :: a

    a = 1

    do i = 1, n
      a = a * alpha; list(i) = a
    end do

  end function PowListComp

!ccccccccccccccc

  pure function PScoef(self, lg) result(list)
    class (AnomDim), intent (in) :: self
    real (dp)      , intent (in) :: lg
    real (dp), dimension(4)      :: list

    list = [ 4._dp/3, 97._dp/9 - 22._dp * self%nf/27, 173.61795863880744_dp &
    - 23.78877355147744_dp * self%nf + 0.6460905349794239_dp * self%nf**2,  &
    3567.7037988167936_dp - 701.2304086375466_dp * self%nf +                &
    41.15507811732836 * self%nf**2 - 0.67466849565615_dp * self%nf**3 +     &
    88.82643960980423_dp * lg ]


  end function PScoef

!ccccccccccccccc

  pure function PowListInt(alpha, n) result(list)
    integer  , intent(in) :: alpha, n
    integer               :: i, a
    integer, dimension(n) :: list

    a = 1

    do i = 1, n
      a = a * alpha; list(i) = a
    end do

  end function PowListInt

!ccccccccccccccc pole - MSbar mass with mu = m(m)

  pure function MSbarDelta(nl, nh, err) result(coef)
    integer            , intent(in) :: nh, nl
    real (dp), optional, intent(in) :: err
    real (dp)        , dimension(4) :: coef

    coef = [ 4._dp/3, 13.33982910125516_dp + 0.10356715567659536_dp * nh          - &
    1.041366911171631_dp * nl, 188.67172035165487_dp + 1.8591544419385237_dp * nh + &
    0.06408045019609998_dp * nh**2 - 26.677375174269212_dp * nl                   + &
    0.022243482163948114_dp * nh * nl + 0.6526907490815437_dp * nl**2, &
    !3560.8519915203624_dp + 6.958064783286616_dp * nh  + 0.02484_dp * nh**3 - &
    !0.23497888797223965_dp * nh**2 - 744.8538175070678_dp * nl - &
    !0.9031405141668719_dp * nh * nl + 0.03617_dp * nh**2 * nl +  &
    !43.37904803138152_dp * nl**2 + 0.017202086604103457_dp * nh * nl**2 - &
    !0.6781410256045151_dp * nl**3 ]
    3560.8898310265354_dp + 6.959412931434764_dp * nh - &
    0.2350144435277955_dp * nh**2 + 0.07451_dp * nh**3/3 - &
    744.8539586181788_dp * nl - 0.9031271808335385_dp * nh * nl + &
    0.108515_dp * nh**2/3 * nl + 43.379048031381515_dp * nl**2 + &
    0.01720208660410344_dp * nh * nl**2 - 0.6781410256045151_dp * nl**3 ]

    if ( present(err) ) coef(4) = coef(4) + err * (1.63_dp + 0.12_dp * nh + &
    0.0027_dp * nh**2 + 0.04_dp * nl**2 + 0.0004_dp * nh * nl)

  end function MSbarDelta

!ccccccccccccccc pole - MSbar mass with m(mu)

  pure function MSbarDeltaPiece(nl, nh) result(coef)
    integer          , intent(in) :: nh, nl
    real (dp), dimension(0:4,4)   :: coef
    real (dp), dimension(0:4,0:4) :: b
    real (dp), dimension(0:4    ) :: beta, gam
    integer                       :: i, j, n, nf

    coef = 0; nf = nl + nh; b = 0; b(0,1:) = MSbarDelta(nl, nh); b(0,0) = 1
    beta = betaFun(nf); gam = massAnomDim(nf)

    do i = 1, 4
      do j = 0, i - 1

        b(j + 1,i) = 0

        do n = 1, i - j - 1
          b(j + 1,i) = b(j + 1,i) + (  ( (i - n) * beta(n - 1) - gam(n - 1) ) &
          * b(j, i - n) + (j + 1) * gam(n - 1) * b(j + 1, i - n)  )/4**n
        end do

        b(j + 1,i) = 2 * ( b(j + 1,i) + ( j * beta(i - j - 1) - gam(i - j - 1) ) &
        * b(j, j)/4**(i - j) )/(j + 1)

      end do
    end do

    coef = b(:,1:)

  end function MSbarDeltaPiece

!ccccccccccccccc

  pure function getInverse(c) result(d)
    real (dp), dimension(:)         , intent(in) :: c
    real (dp), dimension( size(c) )              :: d
    real (dp), dimension( 0:size(c) - 1, 0:size(c) ) :: e
    integer :: n, j

    d(1) = 1; e = 0; e(0,0) = 1

    do n = 2, size(c)

      do j = n - 1, size(c)
        e(n - 1,j) = sum( c(:j - n + 2) * e(n - 2,j - 1:n - 2:-1) )
      end do

      d(n) = - sum( d(:n-1) * e(1:n-1,n) )

    end do

  end function getInverse

!ccccccccccccccc

  pure function massAnomDim(nf) result(gam)
    integer    , intent(in) :: nf
    real (dp), dimension(5) :: gam

    gam = [ - 4._dp, - 16 * (101._dp/24 - 5._dp * nf/36), &
    - 64 * (1249._dp/64 - 2.284121493373736_dp * nf - 35._dp/1296 * nf**2 ), &
    - 256 * (98.9434142552029_dp - 19.1074619186354_dp * nf +  &
    0.005793222354358382_dp * nf**3 + 0.27616255142989465_dp * nf**2 ), &
    573139.8612330721_dp - 147134.94237385172_dp * nf + &
    7661.955304616065_dp * nf**2 + 110.91796496576201_dp * nf**3 - &
    0.08740751602244723_dp * nf**4    ]

  end function massAnomDim

!ccccccccccccccc

  real (dp) function deltaCharmGlue(r)
    real (dp), intent(in)   :: r
    real (dp)               :: lr
    real (dp), dimension(8) :: rPow

    rPow = PowList(r,8); lr = log(r); deltaCharmGlue = 0

    if ( r <= tiny(1._dp) ) return

    deltaCharmGlue = r * ( 146.55437276141436_dp - 68.5502265986958_dp * r + &
    72.20843181777957_dp * rPow(2) + 44.70786258835728_dp * rPow(3) + &
    17.872638638242837_dp * rpow(4) - 0.6287605934381015_dp * rPow(5) + &
    1.1859159626883737_dp * rPow(6) + 0.673242497280189_dp * rPow(7) + &
    1.9017595879357783_dp * rPow(8) )/(6.5849022304691935_dp + r) + &
   lr * r * (      -18.056564104234784_dp + lr * rPow(3) * ( -1.2733398103055087_dp + &
   0.4074074503956954_dp * lr - 0.852006762408542_dp * rpow(2) - &
   0.30539079178379064_dp * rPow(4) ) + r * (     -10.686984081790198_dp + &
   r * (    -18.66917624034727_dp + r * (   -0.2506821619225726_dp + r * &
   (  -2.8512190595320313_dp + r * ( 0.5518674634880775_dp + &
   (-0.5451781523380799_dp - 0.09796100275651486_dp * r) * r )  )   )    )     )      )

  end function deltaCharmGlue

!ccccccccccccccc

  real (dp) function deltaCharmGlueDer(r)
    real (dp), intent(in)   :: r
    real (dp)               :: lr
    real (dp), dimension(9) :: rPow

    rPow = PowList(r,9); lr = log(r); deltaCharmGlueDer = 0

    if ( r <= tiny(1._dp) ) return

    deltaCharmGlueDer = (  182.09667057261015_dp - 1603.9921461502795_dp * r + &
    389.59112788526295_dp * rpow(2) + 1054.578281805314_dp * rPow(3) + &
    576.9693355262091_dp * rPow(4) + 32.777403225777135_dp * rPow(5) + &
    32.29751350364612_dp * rPow(6) + 31.705680846164665_dp * rPow(7) + &
    115.5835006492528_dp * rPow(8) + 15.116115700729708_dp * rPow(9) + &
    lr**3 * rPow(3) * ( 70.66227578689724_dp + 21.461905830563058_dp * r + &
    1.6296298015827817_dp * rPow(2) ) + lr**2 * rPow(3) * ( -167.85612429694254_dp &
    - 50.98211588328542_dp * r - 225.53400814756782_dp * rPow(2) - &
    67.32457476190613_dp * rPow(3) - 111.048288577996_dp * rPow(4) - &
    32.17549609570938_dp * rPow(5) - 2.443126334270325_dp * rPow(6) ) + lr * &
    ( -782.9495455090408_dp - 1164.5967136957051_dp * r - 2728.0864913255296_dp &
    * rPow(2) - 912.8878394351315_dp * rPow(3) - 720.9101971207356_dp * rPow(4) &
    - 121.61007629328383_dp * rPow(5) - 158.565750813136_dp * rPow(6) - &
    109.11754563690172_dp * rPow(7) - 22.181139099099877_dp * rPow(8) - &
    1.3944696056197_dp * rPow(9) )  )/(6.5849022304691935_dp + r)**2

  end function deltaCharmGlueDer

!ccccccccccccccc

  real (dp) function deltaCharmNL(r)
    real (dp), intent(in)   :: r
    real (dp)               :: lr
    real (dp), dimension(8) :: rPow

    rPow = PowList(r,8); lr = log(r); deltaCharmNL = 0

    if ( r <= tiny(1._dp) ) return

    deltaCharmNL = - 0.0025170068027210884_dp * lr * r * ( -435.6852393273681_dp &
    - 435.6852393273681_dp * rPow(2) - 112.12030500101459_dp * rPow(3) - &
    95.64564564564564_dp * lr * rPow(3) + 29.42942942942943_dp * lr**2 * rPow(3) &
    + 10.594594594594595_dp * rPow(5) + rPow(7) ) + (  r * ( -1.5169370316230057_dp &
    - 0.7139583194060264_dp * r + 0.5740149072765712_dp * rPow(2) - &
    0.28690147708463154_dp * rPow(3) - 1.0821107015842863_dp * rPow(4) - &
    0.2198718445207236_dp * rPow(5) + 0.09197784283922195_dp * rPow(6) - &
    0.01027610419769106_dp * rPow(7) + 0.006725122803730151_dp * rPow(8) )  )&
    /(1.4606375667533493_dp + r)

  end function deltaCharmNL

!ccccccccccccccc

  real (dp) function deltaCharmNLDer(r)
    real (dp), intent(in)   :: r
    real (dp)               :: lr
    real (dp), dimension(9) :: rPow

    rPow = PowList(r,9); lr = log(r); deltaCharmNLDer = 0

    if ( r <= tiny(1._dp) ) return

    deltaCharmNLDer = (  0.12390777917265539_dp + 1.1178679723201213_dp * r + &
    5.237550598120457_dp * rPow(2) + 3.277409326369551_dp * rPow(3) - &
    6.842533472692735_dp * rPow(4) - 6.030046014500649_dp * rPow(5) - &
    0.2368358449111415_dp * rPow(6) + 0.3997531410979107_dp * rPow(7) + &
    0.009121104305927087_dp * rPow(8) + 0.05128397562712011_dp * rPow(9) + &
    lr**3 * rPow(3) * ( -0.6321369189366356_dp - 0.8655630025205033_dp * r - &
    0.2962962962962963_dp * rPow(2) ) + lr**2 * rPow(3) * ( 1.5803422973415886_dp &
    + 2.1639075063012583_dp * r + 0.7407407407407407_dp * rPow(2) ) + &
    lr * ( 2.339602993960532_dp + 3.2035366571611794_dp * r + &
    8.115431693113747_dp * rPow(2) + 13.046149089565898_dp * rPow(3) + &
    7.994031844504943_dp * rPow(4) + 1.2689578268985664_dp * rPow(5) - &
    0.46740402136107173_dp * rPow(6) - 0.20295950898079584_dp * rPow(7) - &
    0.05882295506925053_dp * rPow(8) - 0.020136054421768707_dp * rPow(9) )  )&
    /(1.4606375667533493_dp + r)**2

  end function deltaCharmNLDer

!ccccccccccccccc

  real (dp) function deltaCharmNH(r)
    real (dp), intent(in)   :: r
    real (dp)               :: lr
    real (dp), dimension(7) :: rPow

    rPow = PowList(r,7); lr = log(r); deltaCharmNH = 0

    if ( r <= tiny(1._dp) ) return

    deltaCharmNH = rPow(2) * (   ( 0.327412590226266_dp - 0.6530803188815482_dp * r &
    + 0.4620726956320225_dp * rPow(2) - 0.11513821076485932_dp * rPow(3) + &
    0.09025690005091372_dp * rPow(4) - 0.061259252872525434_dp * rPow(5) + &
    0.001927757595523288_dp * rPow(6) - 0.0020551485377330423_dp * rPow(7) )/&
    (1.4733594652016009_dp - r) + lr * (  3.5253490810187646e-7_dp - &
    0.1741826695338979906286_dp * rPow(2) - 0.05629628654108489566051_dp * rPow(4)&
     - 0.00120748327054723502409_dp * rPow(6) + lr * rPow(4) * &
     ( 0.02962963801364881794712_dp + 0.007142868560795129282328_dp * rPow(2) )  )   )

  end function deltaCharmNH

!ccccccccccccccc

  real (dp) function deltaCharmNHDer(r)
    real (dp), intent(in)   :: r
    real (dp)               :: lr
    real (dp), dimension(8) :: rPow

    rPow = PowList(r,8); lr = log(r); deltaCharmNHDer = 0

    if ( r <= tiny(1._dp) ) return

    deltaCharmNHDer = r * (  0.9647936429506729_dp - 3.2140798371305803_dp * r + &
    3.651244040544068_dp * rpow(2) - 1.7211505804162772_dp * rpow(3) + &
    0.9620480117953624_dp * rpow(4) - 0.9171934673489468_dp * rpow(5) + &
    0.33136027956270137_dp * rpow(6) - 0.03718794231137417_dp * rpow(7) + &
    0.015233705031317102_dp * rpow(8) + lr**2 * rPow(4) * ( 0.3859179960794236_dp &
    - 0.5238612914148798_dp * r + 0.3018230614378081_dp * rPow(2) - &
    0.16838420804381504_dp * rPow(3) + 0.057142948486361034_dp * rPow(4) ) + &
   lr * ( 1.5305571763431494e-6_dp - 2.077642574663505e-6_dp * r - &
   1.5124539694764731_dp * rPow(2) + 2.0530694786548094_dp * rPow(3) - &
   1.301335204121848_dp * rPow(4) + 0.8207155691004806_dp * rPow(5) - &
   0.26847665752978905_dp * rPow(6) - 0.013631141519220097_dp * rPow(7) + &
   0.004625870957212379_dp * rPow(8) )  )/(1.4733594652016009_dp - r)**2

 end function deltaCharmNHDer

!ccccccccccccccc

  real (dp) function deltaCharm3(nl, nh, r)
    real (dp), intent(in) :: r
    integer  , intent(in) :: nl, nh
    real (dp)             :: lr

    deltaCharm3 = 0; if ( r <= tiny(1._dp) ) return

    lr = log(r)

    deltaCharm3 = deltaCharmGlue(r) + nl * deltaCharmNl(r) + &
    nh * deltaCharmNh(r) + 2 * deltaCharm2(r)/3

  end function deltaCharm3

!ccccccccccccccc

  real (dp) function deltaCharm3Der(nl, nh, r)
    real (dp), intent(in) :: r
    integer  , intent(in) :: nl, nh
    real (dp)             :: lr

    deltaCharm3Der = 0; if ( r <= tiny(1._dp) ) return

    lr = log(r)

    deltaCharm3Der = deltaCharmGlueDer(r) + nl * deltaCharmNlDer(r) + &
    nh * deltaCharmNhDer(r) + 2 * deltaCharm2Der(r)/3

  end function deltaCharm3Der

!ccccccccccccccc

  real (dp) function deltaCharm2(r)
    real (dp)   , intent(in) :: r
    real (dp)                :: lr, r2, r3, dilog, r4, corr
    real (dp), dimension(10) :: a
    integer                  :: i

    deltaCharm2 = 0; if ( r <= tiny(1._dp) ) return; a = 0

    if (1000 * r <= 1) then

      lr = log(r); r2 = r**2; r4 = r**4

      deltaCharm2 = pi2/6 * r * ( 1 + r2 ) - r2 + (13 * lr/18 - 151._dp/216    - &
      pi2/18 - lr**2/3) * r4

      i = 2

      do
        r4 = r2 * r4; i = i + 1; corr = ( 2 * F(i) * lr + G(i) ) * r4
        deltaCharm2 = deltaCharm2 - 4 * corr/3
        if ( abs(corr) <= 1e-13_dp ) exit
      end do

    else if ( 1000 * abs(r - 1) <= 1 ) then

      a = powList(1 - r, 10)

      deltaCharm2 = pi2/6 - 0.5_dp + (4._dp/3 - 2 * pi2/9) * a(1) - a(5)/90  + &
      (pi2/12 - 1) * a(2) + (pi2/12 - 1) * a(3)/3 - a(6)/180 - a(8)/600 - &
      (pi2 - 9) * a(4)/36 - 2 * a(7)/675 - 779 * a(9)/793800 - 191 * a(10)/317520

    else

      lr = log(r); r2 = r**2; r3 = r * r2

      deltaCharm2 = ( lr**2 + pi2/6 - r2 * (1.5_dp + lr) + &
      (1 + r) * (1 + r3) * ( Dilog(-r) - lr**2/2 + lr * log(1 + r) + pi2/6 ) - &
      (1 - r) * (1 - r3) * ( Dilog(1 - r) + lr**2/2 + pi2/6 ) )/3

    end if

  end function deltaCharm2

!ccccccccccccccc

  real (dp) function deltaCharm2Der(r)
    real (dp)   , intent(in) :: r
    real (dp)                :: lr, r2, r3, dilog, corr
    real (dp), dimension(10) :: a
    integer                  :: i

    deltaCharm2Der = 0; if ( r <= tiny(1._dp) ) return; a = 0

    if (1000 * r <= 1) then

      lr = log(r); r2 = r**2; r3 = r * r2

      deltaCharm2Der = pi2/6 * ( 1 + 3 * r2 ) - 2 * r + (20 * lr/9 - 56._dp/27 - &
      2 * pi2/9 - 4 * lr**2/3) * r3

      i = 2

      do
        r3 = r2 * r3; i = i + 1
        corr = ( 2 * i * F(i) * lr + i * G(i) + F(i)  ) * r3
        deltaCharm2Der = deltaCharm2Der - 8 * corr/3
        if ( abs(corr) <= 1e-13_dp ) exit
      end do

    else if ( 1000 * abs(r - 1) <= 1 ) then

      a = powList(1 - r, 10)

      deltaCharm2Der = 2 * pi2/9 - 4._dp/3 + (2 - pi2/6) * a(1) + a(5)/30 + &
      (1 - pi2/12) * a(2) - (1 - pi2/9) * a(3) + 14 * a(6)/675  + a(4)/18 + &
      a(7)/75 + 779 * a(8)/88200 + 191 * a(9)/31752 + 3337 * a(10)/793800

    else

      lr = log(r); r2 = r**2; r3 = r * r2

      deltaCharm2Der = ( Pi**2/3 - 4 * r + Pi**2 * r2 - 4 * lr**2 * r**3 + &
      lr * ( log(1 + r) * (1 + 3 * r2 + 4 * r3) - 2 * r )                + &
      (1 + 3 * r2 - 4 * r3) * DiLog(1 - r) + (1 + 3 * r2 + 4 * r3) * DiLog(-r) )/3

    end if

  end function deltaCharm2Der

!ccccccccccccccc

  real (dp) function DeltaBottomFun(r, a)
    real (dp), intent(in) :: r, a
    real (dp)             :: lg, a2, a3, a4, r2

    lg = log(a); a2 = a**2; a3 = a2 * a; a4 = a2**2; r2 = r**2

    DeltaBottomFun = 0;  if ( r <= tiny(1._dp) ) return

    DeltaBottomFun = 64.94343831699491_dp * a * r + 31.76898209724912_dp * a2 * r &
    - 0.16416623684141785_dp * a3 * r - 13.407339285356375_dp * a2 * r2 - &
    0.7994215175732678_dp * a3 * r2 + 6.562504757079997_dp * a4 * r2 &
    - 70.5335117178178_dp * a * r * lg - 9.23584680826455_dp * a2 * r * lg + &
    10.068212741421073_dp * a2 * r2 * lg - 14.630200980712981_dp * a3 * r2 * lg

  end function DeltaBottomFun

!ccccccccccccccc

  real (dp) function DeltaBottomFunDer(r, a)
    real (dp), intent(in) :: r, a
    real (dp)             :: lg, a2, r2

    lg = log(a); a2 = a**2; r2 = r**2

    DeltaBottomFunDer = 0;  if ( r <= tiny(1._dp) ) return

    DeltaBottomFunDer = a2 * (13.407339285356375_dp + 0.7994215175732678_dp * a &
    + 6.562504757079997_dp * a2) * r2 - (10.068212741421073_dp - &
    14.630200980712981_dp * a) * a2 * r2 * lg

  end function DeltaBottomFunDer

!ccccccccccccccc

  real (dp) function DeltaBottomCharm(r1, r2)
    real (dp), intent(in) :: r1, r2

    DeltaBottomCharm = 0

    if ( r1 * r2 <= tiny(1.d0) ) return

    if (r1 >= r2) then
      DeltaBottomCharm = DeltaBottomFun(r1, r2/r1)
    else
      DeltaBottomCharm = DeltaBottomFun(r2, r1/r2)
    end if

    DeltaBottomCharm = DeltaBottomCharm/64

  end function DeltaBottomCharm

!ccccccccccccccc

  real (dp) function gammaRBottomCharm(r1, r2)
    real (dp), intent(in) :: r1, r2

    gammaRBottomCharm = 0

    if ( r1 * r2 <= tiny(1.d0) ) return

    if (r1 >= r2) then
      gammaRBottomCharm = DeltaBottomFunDer(r1, r2/r1)
    else
      gammaRBottomCharm = DeltaBottomFunDer(r2, r1/r2)
    end if

  end function gammaRBottomCharm

!ccccccccccccccc

  real (dp) function gammaRcharm2(r)
    real (dp), intent(in) :: r
    gammaRcharm2 = 16 * ( deltaCharm2(r) - r * deltaCharm2Der(r) )
  end function gammaRcharm2

!ccccccccccccccc

  real (dp) function gammaRcharm3(nl, nh, r)
    real (dp), intent(in) :: r
    integer  , intent(in) :: nl, nh
    real (dp)             :: bet

    bet = 11 - 2 * nl/3._dp

    gammaRcharm3 = 64 * ( deltaCharm3(nl, nh, r) - bet * deltaCharm2(r) &
    - r * deltaCharm3Der(nl, nh, r) )

  end function gammaRcharm3

!ccccccccccccccc

  pure real (dp) function F(n)
    integer, intent(in) :: n
    F = 3._dp *  (n - 1)/4/n/(n - 2)/(2 * n - 1)/(2 * n - 3)
  end function F

!ccccccccccccccc

  pure real (dp) function G(n)
    integer, intent(in) :: n
    G = - 3._dp *  (6 - 38 * n + 67 * n**2 - 48 * n**3 + 12 * n**4)/&
    4/(1 - 2 * n)**2/(3 - 2 * n)**2/ (n - 2)**2/ n**2
  end function G

!ccccccccccccccc

  pure function betaFun(nf) result(bet)
    integer    , intent(in) :: nf
    real (dp), dimension(5) :: bet

  bet = [ 11 - 2 * nf/3._dp, 102 - 38._dp * nf/3, 1428.5_dp - 5033 * nf/18._dp + &
  325 * nf**2/54._dp, 29242.964136194132_dp - 6946.289617003555_dp * nf + &
  405.08904045986293_dp * nf**2 + 1.4993141289437586_dp * nf**3, &
  537147.6740702358_dp - 186161.94951432804_dp * nf + &
  17567.757653436838_dp * nf**2 - 231.27767265113647_dp * nf**3 - &
  1.8424744081239024_dp * nf**4 ]

  end function betaFun

!ccccccccccccccc

  real (dp) function Pint(r)
    real (dp), intent(in) :: r
    real (dp)             :: abserr1, abserr2, P0, lr, r2, r4
    integer               :: neval1, neval2, ier1, ier2

    Pint = 0; if (r <= tiny(1._dp) ) return

    lr = log(r); r2 = r**2; r4 = r2**2

    call qags( h, 0._dp, 1._dp, prec, prec, P0  , abserr1, neval1, ier1 )
    call qagi( f, 1._dp, 1    , prec, prec, Pint, abserr2, neval2, ier2 )

    Pint = 32 * (P0 + Pint + r2 * (2 + 33 * r2/4 - 21 * r2 * lr) - &
    (17090689 - 8928/r4 + 125552/r2 + 7018410 * lr)/252000 )/27

  contains

    real (dp) function g(z)
      real (dp), intent(in) :: z
      real (dp)             :: t

      t = z/2
      g = ( t + (1 - t) * sqrt(1 + 4/z) ) * p(r2/z) * ( log(z) - 5._dp/3 )

    end function g

!ccccccccccccccc

    real (dp) function f(z)
      real (dp), intent(in) :: z
      real (dp)             :: lz

      lz = log(z)

      f = g(z) - r2 * (  6 * z - 8 - 3 * r2 + 6 * r2 * ( 2 * log(r) - lz )  ) * &
      (3 * lz - 5)/z**3

    end function f

!ccccccccccccccc

    real (dp) function h(z)
      real (dp), intent(in) :: z
      real (dp)             :: lz, sz, fac

      sz = sqrt(z); lz = log(z); fac = lz - 5._dp/3

      h = g(z) + fac * (  3 * sz**3 * ( 1/r4/70 + 1/r2/20 + 3 * (2 * lr - &
      fac)/64 ) - sz * (3 * lz/4 - 5._dp/4 + 2/r2/5 - 3 * lr/2) - &
      6 * (2 * lr - fac)/3/sz - z * (6 * lr - 3 * fac)/6 )

    end function h

  end function Pint

!ccccccccccccccc

  real (dp) function P2int(r1, r2)
    real (dp), intent(in) :: r1, r2
    integer               :: neval1, neval2, ier1, ier2
    real (dp)             :: abserr1, abserr2, r12, r22, r14, r24, lr1, lr2, P0


    r12 = r1**2; r22 = r2**2; r14 = r12**2; r24 = r22**2
    lr1 = log(r1); lr2 = log(r2)

    call qags( h, 0._dp, 1._dp, prec, prec, P0   , abserr1, neval1, ier1 )
    call qagi( f, 1._dp, 1    , prec, prec, P2int, abserr2, neval2, ier2 )

    P2int = 32 * (   P0 + P2int + 6 * r12 * r22 * (1 - 5 * r12 - 5 * r22 + &
    12 * r12 * lr1 + 12 * r22 * lr2) + (  112 * r12 * (72 + 1121 * r12) * r22 &
    - 8928 * r14 + 30 * lr2 * (3472 * r12 - 288 + 233947 * r14) * r24 + &
    ( 125552 * r12 - 8928 + 17090689 * r14) * r24 + 30 * lr1 * r14 * ( 3472 * r22 &
    + 233947 * r24 - 288 + 124110 * lr2 * r24 )  )/252000/r14/r24  )/27 !+ pInt(r1)

  contains

    real (dp) function g(z)
      real (dp), intent(in) :: z
      real (dp)             :: t

      t = z/2

      g = ( t + (1 - t) * sqrt(1 + 4/z) ) * p(r12/z) * p(r22/z)

    end function g

!ccccccccccccccc

    real (dp) function f(z)
      real (dp), intent(in) :: z
      real (dp)             :: lz

      lz = log(z)

      f = g(z) - 108 * r12 *  r22/z**3 + 18 * r12 * r22 * (8 + 3 * r12 + &
      3 * r22 - 12 * r12 * lr1 - 12 * r22 * lr2 + 6 * (r12 + r22) * lz)/z**4

    end function f

!ccccccccccccccc

    real (dp) function h(z)
      real (dp), intent(in) :: z
      real (dp)             :: lz, sz, fac, fac1, fac2

      sz = sqrt(z); lz = log(z); fac = 5._dp/3 - lz
      fac1 = fac + 2 * lr1;  fac2 = fac + 2 * lr2

      h = g(z) - ( 2 * fac1 * fac2/sz + (  ( 8 * r12 * fac1 - &
      fac2 * r22 * ( 15 * fac1 * r12 - 8) ) * sz  )/20/r12/r22 - &
      (  ( 480 * fac1 * r12 * r14 * r22 + (112 * ( 15 * fac1 * r12 - 8 ) * r14 &
      + 15 * fac2 * (32 * r12 + 7 * (16 + 15 * fac1 * r12) * r14) * r22) * r24 ) &
      * sz**3  )/11200/r12/r14/r22/r24 + fac1 * fac2 * z/2 )

    end function h

  end function P2int

!ccccccccccccccc

  pure real (dp) function p(z)
    real (dp), intent(in) :: z
    real (dp)             :: t, lz

    p = 0

    if ( z <= tiny(1._dp) ) return

      if (z > 20) then

        p = 5._dp/3 + 1/z**9/5819814 - 3/z**8/3695120 + 1/z**7/255255 - &
        1/z**6/51480 + 1/z**5/10010 - 1/z**4/1848 + 1/z**3/315 - 3/z**2/140 &
        + 1/z/5 + log(z)

      else if (100 * z < 1) then

        lz = log(z)
        p = 6 * z - z**9 * (1738808._dp/315 + 4576 * lz) - &
        z**7 * (17364._dp/35 + 432 * lz) - z**5 * (248._dp/5 + 48 * lz) + &
        z**3 * (-16._dp/3 - 8 * lz) + z**2 * (6 * lz - 3) + &
        z**4 * (33._dp/2 + 18 * lz) + z**6 * (463._dp/3 + 140 * lz) + &
        z**8 * (32749._dp/20 + 1386 * lz)

      else

      t = sqrt(1 + 4 * z)
      p = 2 + log(z) - (1 - 2 * z) * (  2 - t * log( (t + 1)/(t - 1) )  )

    end if

  end function p

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

end module AnomDimClass
