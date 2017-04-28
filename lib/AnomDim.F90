
module AnomDimClass
  use Constants, only: dp, d1mach, Pi, Pi2, ExpEuler, prec;  use QuadPack, only: qags
  use adapt, only: dGauss; implicit none ;  private

  real (dp), parameter :: al = 0.634_dp, bet = 1.035_dp, gam = - 23.6_dp, &
  ep = 1.19_dp, del = - 0.481_dp, fourPi = 4 * Pi

  public               :: inteCorre, alphaReExpand, deltaMass, MSbarDelta, &
  PowList, getInverse, MSbarDeltaPiece, AlphaExpand, alphaMatchingLog

  interface PowList
    module procedure   :: PowListDP, PowListInt, PowListComp
  end interface PowList

!ccccccccccccccc

  type, public :: AnomDim
    private
    integer                       :: nf
    character (len = 5)           :: str
    real (dp)                     :: G4, err
    real (dp), dimension(0:3,0:3) :: gammaHm
    real (dp), dimension(0:4)     :: bHat, bCoef, cCoef
    real (dp), dimension(0:4)     :: beta
    real (dp), dimension(0:4)     :: gammaMass
    real (dp), dimension(0:3)     :: gammaHard, gammaB, gammaJet, gammaSoft, &
    cusp, gtilde, gl
    real (dp), dimension(4)       :: sCoefMSR, sCoefMSRNatural, betaList, &
    gammaR,  gammaRNatural, sCoefMSRInc1, gammaRInc1, sCoefMSRInc2, gammaRInc2, &
    sCoefMSRInc3, gammaRInc3

    contains

    procedure, pass(self), public :: expandAlpha, wTildeExpand, kTildeExpand, &
    sCoef, DeltaMu, betaQCD, numFlav, DeltaR, DeltaRHadron, Gfun, DeltaRMass, &
    bHQETgamma, alphaMatching, wTildeHm, GammaRComputer, sCoefRecursive, PScoef,&
    sCoefHadron, scheme, MSRDelta, sCoefLambda, N12, P12, sCoefGeneric, cCoeff, &
    P12Generic, N12Generic

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

   type(AnomDim) function InAdim(str, nf, G4, err)
    character (len = *), intent(in) :: str
    integer            , intent(in) :: nf         ! number of active flavors
    real (dp)          , intent(in) :: G4         ! cusp anomalous dimension
    real (dp), optional, intent(in) :: err        ! 4-loop MS-bar to pole conversion error
    real (dp), dimension(4)         :: betaList
    real (dp), dimension(0:4)       :: beta
    integer                         :: n

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

    InAdim%gammaR        = InAdim%GammaRComputer( InAdim%MSRDelta() )
    InAdim%gammaRNatural = InAdim%GammaRComputer( MSbarDelta(nf    , 0, InAdim%err) )
    InAdim%gammaRInc1    = InAdim%GammaRComputer( MSbarDelta(nf    , 1, InAdim%err) )
    InAdim%gammaRInc2    = InAdim%GammaRComputer( MSbarDelta(nf - 1, 1, InAdim%err) )
    InAdim%gammaRInc3    = InAdim%GammaRComputer( MSbarDelta(nf + 1, 1, InAdim%err) )

    InAdim%sCoefMSR        = InAdim%sCoefRecursive( InAdim%MSRDelta() )
    InAdim%sCoefMSRNatural = InAdim%sCoefRecursive( MSbarDelta(nf   ,  0, InAdim%err) )
    InAdim%sCoefMSRInc1    = InAdim%sCoefRecursive( MSbarDelta(nf   ,  1, InAdim%err) )
    InAdim%sCoefMSRInc2    = InAdim%sCoefRecursive( MSbarDelta(nf - 1, 1, InAdim%err) )
    InAdim%sCoefMSRInc3    = InAdim%sCoefRecursive( MSbarDelta(nf + 1, 1, InAdim%err) )

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
    if ( str( :8) == 'MSRdelta'        ) bet(1:) = self%MSRDelta()
    if ( str(:15) == 'MSRNaturaldelta' ) bet(1:) = MSbarDelta(self%nf, 0, self%err)
    if ( str( :8) == 'sCoefMSR'        ) bet(1:) = self%sCoefMSR
    if ( str(:12) == 'sCoefMSRInc1'    ) bet(1:) = self%sCoefMSRInc1
    if ( str(:12) == 'sCoefMSRInc2'    ) bet(1:) = self%sCoefMSRInc2
    if ( str(:12) == 'sCoefMSRInc3'    ) bet(1:) = self%sCoefMSRInc3
    if ( str(:15) == 'sCoefMSRNatural' ) bet(1:) = self%sCoefMSRNatural
    if ( str( :8) == 'betaList'        ) bet(1:) = self%betaList
    if ( str( :9) == 'gammaRInc1'      ) bet(1:) = self%gammaRInc1
    if ( str(:10) == 'gammaRInc2'      ) bet(1:) = self%gammaRInc2
    if ( str(:10) == 'gammaRInc3'      ) bet(1:) = self%gammaRInc3
    if ( str( :6) == 'gammaR'          ) bet(1:) = self%gammaR
    if ( str(:13) == 'gammaRNatural'   ) bet(1:) = self%gammaRNatural

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
       sCoef = self%sCoefMSRNatural
     else if ( type(:9) == 'Practical' ) then
       sCoef = self%sCoefMSR
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

  pure function alphaMatching(self, nf) result(tab)
    class (AnomDim), intent(in)     :: self
    integer        , intent(in)     :: nf
    real (dp)  , dimension(5)       :: tab, tab2
    real (dp)  , dimension(4)       :: b
    real (dp)  , dimension(0:4,0:4) :: d
    real (dp)  , dimension(0:4,5)   :: c
    integer                         :: i, j, k

    tab = 0; tab2 = 0; tab2(1) = 1
    tab2(3:) = - [7._dp/24, 5.586361025786356_dp - 0.26247081195432964_dp * nf,&
    95.80685705811597_dp - 10.171370332526523_dp * nf + 0.23954195789610408_dp * nf**2]

    if ( self%str(:5) == 'MSbar' ) then

      b = MSbarDelta(nf - 1, 1); d = 0; d(0,0) = 1

      do i = 1, 4

        do j = 1, i - 1
          d(1,i) = d(1,i) + j * d(1,j) * b(i - j)
        end do

        d(1,i) = b(i) - d(1,i)/i

      end do

      c = alphaMatchingLog(tab2, nf)

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

    else if ( self%str(:4) == 'pole' ) then

      tab = tab2

    end if

  end function alphaMatching

!ccccccccccccccc

  pure function alphaMatchingLog(d, nf) result(tab)
    integer                , intent(in) :: nf
    real (dp), dimension(5), intent(in) :: d
    real (dp), dimension(0:4,5)         :: tab  ! (log, alpha)
    real (dp), dimension(0:4,0:5)       :: b, c ! (log, alpha)
    real (dp), dimension(0:4,5,5)       :: ePow ! (log,alpha,power)
    integer                             :: n, i, j, k, l, m

    tab = 0; tab(0,1) = 1; ePow = 0

    b = 0; c = 0; c(0,1) = 1 ;  b(0,1:) = d
    call AlphaExpand( betaFun(nf    ), b(:,1:) )
    call AlphaExpand( betaFun(nf - 1), c(:,1:) )

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

 pure function MSRDelta(self) result(coef)
    class (AnomDim), intent(in)   :: self
    real (dp), dimension(4)       :: coef, b

    coef = MSbarDelta(self%nf, 1, self%err)
    b = getInverse( self%alphaMatching(self%nf + 1) ); call alphaReExpand(coef, b)

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
    0.2350144435277955_dp * nh**2 + 0.024836666666666663_dp * nh**3 - &
    744.8539586181788_dp * nl - 0.9031271808335385_dp * nh * nl + &
    0.03617166666666667_dp * nh**2 * nl + 43.37905086931079_dp * nl**2 + &
    0.017201666666666664_dp * nh * nl**2 - 0.6781416666666666_dp * nl**3 ]

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

  pure function betaFun(nf) result(bet)
    integer    , intent(in) :: nf
    real (dp), dimension(5) :: bet

  bet = [ 11 - 2 * nf/3._dp, 102 - 38._dp * nf/3, 1428.5_dp - 5033 * nf/18._dp + &
  325 * nf**2/54._dp, 29242.964136194132_dp - 6946.289617003555_dp * nf + &
  405.08904045986293_dp * nf**2 + 1.4993141289437586_dp * nf**3, &
  537147.6740702358_dp - 186161.94951432804_dp * nf + &
  17567.757653436835_dp * nf**2 - 231.27767265113647_dp * nf**3 - &
  1.8424744081239026_dp * nf**4 ]

  end function betaFun

!ccccccccccccccc

end module AnomDimClass
