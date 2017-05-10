module NRQCDClass
  use Constants, only: dp, d1mach, Pi ;  use RunningClass; use AlphaClass
  use AnomDimClass; implicit none ; private

!ccccccccccccccc

  type, public :: NRQCD
    private
    real (dp), dimension(0:4,0:3) :: c
    character (len = 5)           :: scheme
    real (dp)                     :: mH
    type (Running)                :: alphaMass
    type (AnomDim)                :: Andim
    integer                       :: n, nl
    integer, dimension(0:3)       :: listFact

  contains

    procedure, pass(self), public  :: En
    procedure, pass(self), private :: Binomial

  end type NRQCD

!ccccccccccccccc

  interface NRQCD
    module procedure InNRQCD
  end interface NRQCD

  contains

!ccccccccccccccc

  type (NRQCD) function InNRQCD(scheme, alphaMass, n, l, j, s)
    integer            , intent(in) :: n, l, j, s
    character (len = *), intent(in) :: scheme
    type (Running)     , intent(in) :: alphaMass
    real (dp)     , dimension(0:4)  :: beta
    character (len = 5)             :: alphaScheme
    integer                         :: nl, i, jj, k

    InNRQCD%alphaMass  = alphaMass; InNRQCD%c = 0; InNRQCD%n = n; nl = alphaMass%numFlav()
    InNRQCD%Andim = InNRQCD%alphaMass%adim();  beta  = InNRQCD%Andim%betaQCD('beta')
    InNRQCD%mH = InNRQCD%alphaMass%scales('mH'); alphaScheme = alphaMass%scheme()
    InNRQCD%listFact = factList(3); InNRQCD%nl = nl

    if (  alphaScheme(:4) == 'pole') then
      InNRQCD%scheme = 'pole'
    else
      InNRQCD%scheme = scheme
    end if

    InNRQCD%c(:,0) = [ 1._dp, 31._dp/6 - 5._dp * nl/9, c2(nl, n, l, j, s), &
    c3(nl, n, l, j, s), c3log(nl, n, l, j, s) ]

    do k = 1, 3
      do jj = 0, k - 1

        do i = jj + 1, k - 1

          InNRQCD%c(k,jj + 1) = InNRQCD%c(k,jj + 1) + beta(k - 1 - i) * &
          ( (i + 2) * InNRQCD%c(i, jj) - (jj + 1) * InNRQCD%c(i, jj + 1) )/4**(k - i)

        end do

        InNRQCD%c(k,jj + 1) = 2 * (  InNRQCD%c(k,jj + 1) + &
        (jj + 2) * InNRQCD%c(jj,jj) * beta(k - 1 - jj)/4**(k - jj)  )/(jj + 1)

      end do
    end do

  end function InNRQCD

!ccccccccccccccc

  function En(self, order, mu, R, lambda, method) result(list)
    class (NRQCD)      , intent(in) :: self
    character (len = *), intent(in) :: method
    integer            , intent(in) :: order
    real (dp)          , intent(in) :: mu, R, lambda
    real (dp), dimension(0:4)       :: list, alphaList
    real (dp), dimension(0:3)       :: logList
    real (dp), dimension(0:4)       :: delta, lgmList
    real (dp), dimension(0:1,0:2)   :: deltaLog  ! (power, order)
    real (dp), dimension(0:4,4)     :: coefMSR
    real (dp), dimension(4,0:3)     :: c
    real (dp)                       :: alp, Rmass, mass, factor
    integer                         :: i, j, k

    list(0) = 2     ; alp           = self%alphaMass%alphaQCD(mu); coefMSR = 0
    alphaList(0) = 1; alphaList(1:) = PowList(alp/Pi,4); delta(0) = 1
    factor = - 4 * alp**2/9/self%n**2 ; logList(0) = 1

    if ( self%scheme(:4) == 'pole' ) then
      delta(1:) = 0; Rmass = 0; mass = self%mH
    else if ( self%scheme(:5) == 'MSbar' ) then
      coefMSR(0,:) = self%mH * self%andim%MSRDelta()
      Rmass = self%mH; mass = self%mH
    else if ( self%scheme(:4) == 'MSRp' ) then
      coefMSR(0,:) = R * self%andim%MSRDelta(); Rmass = R
      mass = self%alphaMass%MSRmass(R, lambda, method)
    else if ( self%scheme(:4) == 'MSRn' ) then
      coefMSR(0,:) = R * MSbarDelta(self%nl, 0); Rmass = R
      mass = self%alphaMass%MSRNaturalMass(order, R, lambda, method)
    end if

    logList(1:) = PowList( log(3 * self%n * mu / 4 / alp / mass) , 3)

    if ( self%scheme(:4) == 'pole' ) then

      list(1:) = matmul( self%c(:3,:), logList );  list(4) = list(4) + self%c(4,0)
      list(1:) = factor * alphaList(:3) * list(1:)

    else

      call self%andim%expandAlpha(coefMSR); lgmList = PowList( log(mu/Rmass), 4 )
      call AddAlpha( coefMSR, alphaList(1:) )        ; c = 0; deltaLog = 0
      delta(1:) = DeltaComputer(coefMSR, lgmList, 0)/mass; deltaLog(0,0) = 1
      deltaLog(1,1:) = delta(1:2) - [ 0._dp, delta(1)**2/2 ]

      do i = 1, 4
        do k = 0, i - 1

          do j = k, i - 1
            c(i,k) = c(i,k) + self%c(i - 1,j) * self%Binomial(i,j) * logList(j - k)
          end do

          c(i,k) = (-1)**k * c(i,k)

        end do

        c(i,:) = factor * alphaList(i - 1) * c(i,:)

      end do

      ! list(1:) = c(:,0)

      do j = 1, 4
        do i = 1, j
          do k = 0, min(i - 1, j - i)
            list(j) = list(j) + deltaLog(k,j - i) * c(i,k)
          end do
        end do
      end do

      do j = 4, 1, -1
        list(j) = sum( delta(0:j) * list(j:0:-1) )
      end do

    end if

    list = mass * list

  end function En

!ccccccccccccccc

  pure real (dp) function c2(nl, n, l, j, s)
    integer, intent(in) :: nl, n, l, j, s
    integer             :: nl2, ss, ss2, jj, jj2

    c2 = 0; nl2 = nl**2; ss = s * (1 + s); ss2 = ss**2; jj = j * (j + 1)
    jj2 = jj**2

    if (n == 1) then

      if (l == 0) then
        c2 = 220.34174177286914_dp - 19.798187714015196_dp * nl + &
        0.496190504426009_dp * nl2 - 11.69730891980961_dp * ss
      end if

    else if (n == 2) then

      if (l == 0) then
        c2 = 234.57814959979672_dp - 24.946605003852223_dp * nl + &
        0.6522031495725855_dp * nl2 - 5.848654459904805_dp * ss
      else if (l == 1) then
        c2 = 146.39649245990387_dp - 1.7545963379714415_dp * jj   + &
        0.2193245422464302_dp  * jj2 - 16.951844980548632_dp * nl + &
        0.40993769432096167_dp * nl2 + (1.169730891980961_dp      - &
        0.4386490844928604_dp * jj) * ss + 0.2193245422464302_dp * ss2
      end if

    else if (n == 3) then

      if (l == 0) then
        c2 = 243.61209732248133_dp - 27.345022293689244_dp * nl + &
        0.7248824613858287_dp * nl2 - 3.8991029732698697_dp * ss
      else if (l == 1) then
        c2 = 170.30854367434102_dp - 1.169730891980961_dp * jj    +  &
        0.14621636149762013_dp * jj2 - 20.255693524026842_dp * nl +  &
        0.5100543168506044_dp * nl2 + (0.7798205946539739_dp      -  &
        0.29243272299524026_dp * jj) * ss + 0.14621636149762013_dp * ss2
      else if (l == 2) then
        c2 = 132.27874873495742_dp - 0.2228058841868497_dp * jj     + &
        0.0069626838808390535_dp * jj2 - 16.128767492567004_dp * nl + &
        0.3849959522609123_dp * nl2 + (0.1671044131401373_dp        - &
        0.013925367761678107_dp * jj) * ss + 0.0069626838808390535_dp * ss2
      end if

    else if (n == 4) then
      if (l == 0) then
        c2 = 249.6485479478922_dp - 28.758871682291705_dp * nl + &
        0.7677263822525698_dp * nl2 - 2.9243272299524024_dp * ss
      else if (l == 1) then
        c2 = 186.45187350545456_dp - 0.8772981689857208_dp * jj  + &
        0.1096622711232151_dp * jj2 - 22.445631564883403_dp * nl + &
        0.5764160756644395_dp * nl2 + (0.5848654459904805_dp     - &
        0.2193245422464302_dp * jj) * ss + 0.1096622711232151_dp * ss2
      else if (l == 2) then
        c2 = 150.96237118757924_dp - 0.1671044131401373_dp * jj    + &
        0.0052220129106292906_dp * jj2 - 18.50593326997487_dp * nl + &
        0.45703127884902944_dp * nl2 + (0.12532830985510296_dp     - &
        0.010444025821258581_dp * jj) * ss + 0.0052220129106292906_dp * ss2
      else if (l == 3) then
        c2 = 126.83861491499124_dp - 0.059182812987131954_dp * jj   + &
        0.0008703354851048817_dp * jj2 - 15.738818473520091_dp * nl + &
        0.37317931532009674_dp * nl2 + (0.04525744522545385_dp      - &
        0.0017406709702097634_dp * jj) * ss + 0.0008703354851048817_dp * ss2
      end if

    end if

  end function c2

!ccccccccccccccc

  pure real (dp) function c3log(nl, n, l, j, s)
    integer, intent(in) :: nl, n, l, j, s
    integer             :: nl2, ss, ss2, jj, jj2

    c3log = 0; nl2 = nl**2; ss = s * (1 + s); ss2 = ss**2; jj = j * (j + 1)
    jj2 = jj**2

    if (n == 1) then

      if (l == 0) then
        c3log = 597.1110662659062_dp - 61.41087182900045_dp * ss
      end if

    else if (n == 2) then

      if (l == 0) then
        c3log = 329.5351247252613_dp - 30.705435914500224_dp * ss
      else if (l == 1) then
        c3log = 114.45085696226215_dp - 4.1671663026821735_dp * jj + &
        0.6579736267392905_dp * jj2 + (2.4125699647107317_dp     - &
        1.315947253478581_dp  * jj) * ss + 0.6579736267392905_dp * ss2
      end if

    else if (n == 3) then

      if (l == 0) then
        c3log = 236.44404123844322_dp - 20.470290609666815_dp * ss
      else if (l == 1) then
        c3log = 93.05452939644374_dp - 2.7781108684547826_dp * jj + &
        0.4386490844928604_dp * jj2 + (1.6083799764738214_dp    - &
        0.8772981689857208_dp * jj) * ss + 0.4386490844928604_dp * ss2
      else if (l == 2) then
        c3log = 72.13862701840323_dp - 0.522201291062929_dp * jj + &
        0.020888051642517162_dp * jj2 + (0.3550968779227917_dp - &
        0.041776103285034324_dp * jj) * ss + 0.020888051642517162_dp * ss2
      end if

    else if (n == 4) then
      if (l == 0) then
        c3log = 189.16741768754605_dp - 15.352717957250112_dp * ss
      else if (l == 1) then
        c3log = 81.62528380604644_dp - 2.0835831513410867_dp * jj + &
        0.32898681336964525_dp * jj2 + (1.2062849823553659_dp   - &
        0.6579736267392905_dp * jj) * ss + 0.32898681336964525_dp * ss2
      else if (l == 2) then
        c3log = 65.93835702251604_dp - 0.3916509682971967_dp * jj + &
        0.01566603873188787_dp * jj2 + (0.2663226584420938_dp   - &
        0.03133207746377574_dp * jj) * ss + 0.01566603873188787_dp * ss2
      else if (l == 3) then
        c3log = 59.17062829034048_dp - 0.1383833421316762_dp * jj  + &
        0.0026110064553146453_dp * jj2 + (0.09660723884664187_dp - &
        0.0052220129106292906_dp * jj) * ss + 0.0026110064553146453_dp * ss2
      end if

    end if

  end function c3log

!ccccccccccccccc

  pure real (dp) function c3(nl, n, l, j, s)
    integer, intent(in) :: nl, n, l, j, s
    integer             :: nl3, nl2, ss, ss2, jj, jj2

    c3 = 0; nl2 = nl**2; ss = s * (1 + s); ss2 = ss**2; jj = j * (j + 1)
    jj2 = jj**2; nl3 = nl**3

    if (n == 1) then

      if (l == 0) then
        c3 = 1928.7628660728874_dp - 418.00274061279185_dp * nl   + &
        27.350849570543684_dp * nl2 - 0.4478789305261393_dp * nl3 + &
        (218.58922102668993_dp - 11.52783363067168_dp * nl) * ss
      end if

    else if (n == 2) then

      if (l == 0) then
        c3 = 1555.6647767522752_dp - 427.2864081082992_dp * nl    + &
        29.077705021060194_dp * nl2 - 0.4700406103132774_dp * nl3 + &
        (189.25038089189587_dp - 10.715520511240458_dp * nl) * ss
      else if (l == 1) then
        c3 = 1919.4334788948734_dp + jj * (32.43347193767492_dp    - &
        1.224489010765641_dp * nl) + jj2 * (-3.9574555630864614_dp + &
        0.16829199733504055_dp * nl) - 412.57539959999093_dp * nl  + &
        25.345102512389868_dp * nl2 - 0.41382274015759046_dp * nl3 + &
        (-22.002104070692376_dp + jj * (7.914911126172923_dp       - &
        0.3365839946700811_dp * nl) + 0.9381729750917773_dp * nl) * ss + &
        (-3.9574555630864614_dp + 0.16829199733504055_dp * nl) * ss2
      end if

    else if (n == 3) then

      if (l == 0) then
        c3 = 1419.3458052875883_dp - 418.4769180394629_dp * nl   + &
        28.6078892939055_dp * nl2 - 0.45420057781924034_dp * nl3 + &
        (158.96019573759003_dp - 9.145048480340092_dp * nl) * ss
      else if (l == 1) then
        c3 = 1988.787252623633_dp + jj * (30.640977403131618_dp      - &
        1.3788285034575731_dp  * nl) + jj2 * (-3.7569974635218113_dp + &
        0.18250747692508693_dp * nl) - 444.90038297813294_dp * nl    + &
        27.73818853998447_dp * nl2 - 0.4544693728045063_dp * nl3     + &
        (-20.703548812349908_dp + jj * (7.513994927043625_dp         - &
        0.36501495385017385_dp * nl) + 1.0004503142481713_dp * nl) * ss + &
        (-3.7569974635218113_dp + 0.18250747692508693_dp * nl) * ss2
      else if (l == 2) then
        c3 = 1901.5017267203552_dp + jj * (4.079182167775893_dp      - &
        0.11282798880462161_dp * nl) + jj2 * (-0.1299105555452488_dp + &
        0.0040335703497889385_dp * nl) - 402.34691146347046_dp * nl  + &
        24.612549396122763_dp * nl2 - 0.4008723708453693_dp * nl3    + &
        (-3.0561439858025263_dp + jj * (0.25982111109049766_dp       - &
        0.008067140699577877_dp * nl) + 0.10222110919114273_dp * nl) * ss + &
        (-0.1299105555452488_dp + 0.0040335703497889385_dp * nl) * ss2
      end if

    else if (n == 4) then
      if (l == 0) then
        c3 = 1334.1657153067044_dp - 407.28551220371247_dp * nl  + &
        27.83594776731537_dp * nl2 - 0.4347159974378718_dp * nl3 + &
        (137.000933835361_dp - 7.940117353183779_dp * nl) * ss
      else if (l == 1) then
        c3 = 1995.6832394075916_dp + jj * (27.958016166615735_dp       - &
        1.3404873242204196_dp * nl) + jj2 * (-3.4373507564943924_dp    + &
        0.17517635102222015_dp * nl) - 457.0465561083301_dp * nl       + &
        28.689641476076723_dp * nl2 - 0.4683744068804726_dp * nl3      + &
        (-18.85267096658803_dp + jj * (6.874701512988785_dp            - &
        0.3503527020444403_dp * nl) + 0.9545817001042881_dp * nl) * ss + &
        (-3.4373507564943924_dp + 0.17517635102222015_dp * nl) * ss2
      else if (l == 2) then
        c3 = 1997.2069835963655_dp + jj * (3.958438263203731_dp          - &
        0.14118158483048285_dp * nl) + jj2 * (-0.12522904650665853_dp    + &
        0.004792696300685975_dp * nl) - 431.6105518485525_dp * nl        + &
        26.621747666930165_dp * nl2 - 0.4360214587004526_dp * nl3        + &
        (-2.968790587941932_dp + jj * (0.25045809301331706_dp            - &
        0.00958539260137195_dp * nl) + 0.11908627681361951_dp * nl) * ss + &
        (-0.12522904650665853_dp + 0.004792696300685975_dp * nl) * ss2
      else if (l == 3) then
        c3 = 1882.185489900684_dp + jj * (1.1136805133484151_dp             - &
        0.02558450978851483_dp * nl) + jj2 * (-0.017109601029735064_dp      + &
        0.00044023804726528275_dp * nl) - 396.62888749856256_dp * nl        + &
        24.241424872934203_dp * nl2 - 0.3942467526859063_dp * nl3           + &
        (-0.844278574298179_dp + jj * (0.03421920205947013_dp               - &
        0.0008804760945305655_dp * nl) + 0.024342937599636186_dp * nl) * ss + &
        (-0.017109601029735064_dp + 0.00044023804726528275_dp * nl) * ss2
      end if

    end if

  end function c3

!ccccccccccccccc

   real (dp) function Binomial(self, i, j)
     class (NRQCD), intent(in) :: self
     integer      , intent(in) :: i, j
     integer                   :: k

    select case(j)
    case(0)  ;  Binomial = 1
    case(:-1);  Binomial = 0
    case default

      Binomial = i

      do k = 1, j - 1
        Binomial = Binomial * (i - k)
      end do

      Binomial = Binomial/self%Listfact(j)

    end select

   end function Binomial

!ccccccccccccccc

  pure function factList(n) result(list)
    integer, intent(in)     :: n
    integer                 :: i
    integer, dimension(0:n) :: list

    list(0) = 1

    do i = 1, n
      list(i) = list(i - 1) * i
    end do

  end function factList

!ccccccccccccccc

end module NRQCDClass
