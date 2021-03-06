
module MCtopClass
  use Constants, only: dp, prec, Pi; use QuadPack, only: qags; use Legendre
  implicit none; private

  public                                 :: MCtop, MCScales, BreitWigner

!ccccccccccccccc

  type, public                           :: MCtop
    character (len = 6)                  :: shape
    real (dp)                            :: ESmin, ESmax, pmax, Dirac, Q
    integer                              :: n, sf
    real (dp), dimension(:), allocatable :: coefs

  contains

  !  final                          :: delete_object
   procedure, pass(self), public  :: Distribution, QDist, Delta, maxES, Qval, &
   maxP, BreitUnstable, setMass

  end type MCtop

!ccccccccccccccc

  type, extends(MCtop), public :: MCScales

  contains

  !  final                          :: delete_scales

  end type MCScales

!ccccccccccccccc

  interface MCtop
    module procedure InMCtop
  end interface MCtop

!ccccccccccccccc

  interface MCScales
    module procedure InMCScales
  end interface MCScales

contains

!ccccccccccccccc

  ! subroutine delete_object(self)
  !   type (MCtop) :: self
  !   if ( allocated(self%coefs) ) deallocate(self%coefs)
  ! end subroutine delete_object

!ccccccccccccccc

  ! subroutine delete_scales(self)
  !   type (MCScales) :: self
  !   if ( allocated(self%coefs) ) deallocate(self%coefs)
  ! end subroutine delete_scales

!ccccccccccccccc

  type (MCtop) function InMCtop(Eshape, mt, Q, n, ESmin)
    real (dp), optional, intent(in) :: ESmin
    integer  , optional, intent(in) :: n
    character (len = *), intent(in) :: Eshape
    real (dp)          , intent(in) :: mt, Q

    InMCtop%ESmin = 0; InMCtop%shape = Eshape

    if ( Eshape(:6) == 'thrust' ) then
      allocate( InMCtop%coefs(10) ); InMCtop%n = 0; InMCtop%sf = 1

    else if ( EShape(:6) == 'Cparam') then

      allocate( InMCtop%coefs(0:39) ); InMCtop%sf = 6

      if ( present(n) ) then
        InMCtop%n = min(39,n)
      else
        InMCtop%n = 39
      end if

      if ( present(ESmin) ) InMCtop%ESmin = ESmin

    end if

    call InMCtop%setMass(mt, Q)

  end function InMCtop

!ccccccccccccccc

  type (MCScales) function InMCScales(Eshape, n)
    integer  , optional, intent(in) :: n
    character (len = *), intent(in) :: Eshape

    InMCScales%ESmin = 0; InMCScales%shape = Eshape

    if ( Eshape(:6) == 'thrust' ) then
      allocate( InMCScales%coefs(10) ); InMCScales%n = 0; InMCScales%sf = 1

    else if ( EShape(:6) == 'Cparam') then

      allocate( InMCScales%coefs(0:39) ); InMCScales%sf = 6

      if ( present(n) ) then
        InMCScales%n = min(39,n)
      else
        InMCScales%n = 39
      end if

    end if

  end function InMCScales

!ccccccccccccccc

  subroutine setMass(self, mt, Q)
    class (MCtop), intent(inout) :: self
    real (dp)    , intent(in)    :: mt, Q
    real (dp)                    :: moQ, moQ2

    moQ = mt/Q; moQ2 = moQ**2; self%Q = Q

    if ( self%shape(:6) == 'thrust' ) then

      self%ESmax = 1 - sqrt(1 - 4 * moQ2)
      self%Dirac = 1 - 5.996687441937819_dp * moQ2 - 4.418930183289158_dp * &
      moQ2**2 + 11.76452036058563_dp * moQ2**3
      self%coefs = ThrustCoefs(moQ)
      self%coefs(10) = self%coefs(10) * self%ESmax/(1 - self%Dirac)

    else if ( self%Shape(:6) == 'Cparam') then
      self%ESmax = 12 * moQ2 * (1 - 3 * moQ2); self%Dirac = 0

      self%coefs = LagCoef(moQ)/self%ESmax

    end if

    self%pmax = Q * self%ESmax/self%sf

  end subroutine setMass

!ccccccccccccccc

  real (dp) function Qval(self)
    class (MCtop), intent(in) :: self
    Qval = self%Q
  end function Qval

!ccccccccccccccc

  real (dp) function QDist(self, p, p2)
    class (MCtop)      , intent(in) :: self
    real (dp)          , intent(in) :: p
    real (dp), optional, intent(in) :: p2

    if ( .not. present(p2) ) then
      QDist = self%sf * self%Distribution(0, self%sf * p/self%Q)/self%Q
    else
      QDist = self%sf * self%Distribution(0, self%sf * p/self%Q, self%sf * p2/self%Q)/self%Q
    end if

  end function QDist

!ccccccccccccccc

  real (dp) function BreitUnstable(self, Gamma, k, l, l2)
    class (MCtop)      , intent(in) :: self
    real (dp)          , intent(in) :: l, Gamma
    real (dp), optional, intent(in) :: l2
    integer            , intent(in) :: k
    real (dp)                       :: lint, abserr
    integer                         :: ier, neval

    BreitUnstable = 0; lint = l

    if ( present (l2) ) then
      if (l >= l2) return; lint = l
    end if

    call qags( integrand, l - self%pmax, lint, prec, &
    prec, BreitUnstable, abserr, neval, ier )

  contains

    real (dp) function integrand(x)
      real (dp), intent(in) :: x

      if (.not. present(l2) ) then
        integrand = BreitWigner(Gamma, k, x) * self%QDist( l - x )
      else
        integrand = BreitWigner(Gamma, k, x) * self%QDist( l - x, l2 - x )
      end if

    end function integrand

  end function BreitUnstable

!ccccccccccccccc

  recursive real (dp) function Distribution(self, k, x, x2) result(res)
    class (MCtop)       , intent(in) :: self
    real (dp)           , intent(in) :: x
    real (dp), optional , intent(in) :: x2
    integer             , intent(in) :: k
    real (dp)                        :: y

    res = 0

    if ( present(x2) ) then
      if (x2 > x) res = self%Distribution(k, x2) - self%Distribution(k, x); return
    end if

    if (x <= self%ESmin) return
    if (x >= self%ESmax .and. k /= -1) return

    y = (x - self%ESmin)/(self%ESmax - self%ESmin)

    if ( self%shape(:6) == 'Cparam' ) then

      if (y > 1 .and. k == -1) y = 1

      res = (2._dp)**k * dot_product( self%coefs(0:self%n), &
      LegendreList(self%n, k, 2 * y - 1 ) )

      if ( k == 0 .and. res < 0 .and. x > 0.8_dp * self%ESmax ) res = 0

    else if ( self%shape(:6) == 'thrust' ) then

      if ( 2 * y <= 1 ) then
        if (k == 0) then
          res = self%coefs(1) * y**3
        else if (k == 1) then
          res = 3 * self%coefs(1) * y**2
        else if (k == 2) then
          res = 6 * self%coefs(1) * y
        else if (k == 3) then
          res = 6 * self%coefs(1)
        else if (k > 3) then
          res = 0
        else if (k == -1) then
          res = self%coefs(1) * y**4/4
        end if
      else if ( y <= 0.61_dp ) then
        if (k == 0) then
          res = dot_product( self%coefs(8:9), [y, 1._dp] )
        else if (k == 1) then
          res = self%coefs(8)
        else if (k > 1) then
          res = 0
        else if (k == -1) then
          res =  dot_product( self%coefs(8:9), [y**2/2 - 0.125_dp, y - 0.5_dp] ) &
          + self%coefs(1)/64
        end if
      else if ( y <= 0.89_dp ) then
        if (k == 0) then
          res = dot_product( self%coefs(5:7), [y**2, y, 1._dp] )
        else if (k == 1) then
          res = dot_product( self%coefs(5:6), [2 * y, 1._dp] )
        else if (k == 2) then
          res = 2 * self%coefs(5)
        else if (k > 2) then
          res = 0
        else if (k == -1) then
          res = dot_product( self%coefs(8:9), [0.06105_dp, 0.11_dp] ) &
            + self%coefs(1)/64+ dot_product( self%coefs(5:7), &
            [ (y**3 - 0.226981_dp)/3, y**2/2 - 0.18605_dp, y - 0.61_dp] )
        end if
      else if (y < 1) then
        if (k == 0) then
          res = dot_product( self%coefs(3:4), [y, 1._dp] )
        else if (k == 1) then
          res = self%coefs(3)
        else if (k > 1) then
          res = 0
        else if (k == -1) then
          res = dot_product( self%coefs(8:9), [0.06105_dp, 0.11_dp] )        + &
          dot_product( self%coefs(3:4), [y**2/2 - 0.39605_dp, y - 0.89_dp] ) + &
          self%coefs(1)/64 + dot_product( self%coefs(5:7),                     &
          [ 0.477988_dp/3, 0.21_dp, 0.28_dp] )
        end if
      else
        if (k == -1) then
          res = dot_product( self%coefs(8:9), [0.06105_dp, 0.11_dp] )            + &
          dot_product( self%coefs(3:4), [0.10395_dp, 0.11_dp] ) + self%coefs(1)/64 &
          + dot_product( self%coefs(5:7), [ 0.477988_dp/3, 0.21_dp, 0.28_dp] )
        end if
      end if

      res = res/self%coefs(10)

    end if

    res = res / (self%ESmax - self%ESmin)**k

   end function Distribution

!ccccccccccccccc

  real (dp) function Delta(self)
    class (MCtop), intent(in) :: self

    delta = self%Dirac

   end function Delta

!ccccccccccccccc

  real (dp) function maxES(self)
    class (MCtop), intent(in) :: self

    maxES = self%ESMax

   end function maxES

!ccccccccccccccc

  real (dp) function maxp(self)
    class (MCtop), intent(in) :: self

    maxp = self%pMax

  end function maxp

!ccccccccccccccc

  function ThrustCoefs(m) result(res)
    real (dp)   , intent(in) :: m
    real (dp), dimension(11) :: res

    res(:6) = [ 0.3310188452495506_dp + 29.408029059401223_dp * m**2 + 420.78161872939825_dp * m**4, &
             10.582416684643045_dp + 59.168098208694204_dp * m**2 - 86.0283647010381_dp * m**4, &
             0.5661957300078854_dp + 24.96314852772734_dp * m - 199.1120591815817_dp * m**2, &
             9.343012917211_dp - 16.997770889484364_dp * m + 141.4278848881199_dp * m**2, &
             42.85501883992358_dp + 50.26811493972574_dp * m - 540.819475023059_dp * m**2, &
             -43.27240470704688_dp - 52.90240803297828_dp * m + 655.761339938332_dp * m**2 ]

    res(7)   = res(4) + 0.89_dp * ( res(3) - res(6) ) - res(5) * 0.89_dp**2
    res(8:9) = [ res(2), res(1) * 0.5_dp**3 - res(2) * 0.5_dp ]

    ! norm computation

    res(10) = res(1) * 0.5_dp**4/4 + res(9) * 0.11_dp + res(8) * ( 0.61_dp**2 - 0.5_dp**2)/2 &
    + res(7) * 0.28_dp + res(6) * ( 0.89_dp**2 - 0.61_dp**2)/2 + res(5) * ( 0.89_dp**3 - 0.61_dp**3)/3 &
    + res(4) * 0.11_dp + res(3) * ( 1 - 0.89_dp**2)/2

  end function ThrustCoefs

!ccccccccccccccc

  function LagCoef(m) result(res)
    real (dp)     , intent(in) :: m
    real (dp), dimension(0:39) :: res

    res = 0; if ( m < 0.07_dp .and. m > 0.3_dp) return

    res = [1._dp, 3.0246456232675523_dp - 2.101563195974475_dp * m - 13.429530904866745_dp * m**2 &
    - 244.72141333301047_dp * m**3 &
 + 4589.027290813364_dp * m**4 - 38050.32937871064_dp * m**5 + 190569.80718322814_dp * m**6 &
 - 582915.5096903557_dp * m**7 + 998671.1223194505_dp * m**8 - 731363.6896566917_dp * m**9, &
5.050610576623612_dp - 5.874145317417216_dp * m - 175.3737693027115_dp * m**2 + 909.2095465998219_dp * m**3 &
 + 1345.2853844483411_dp * m**4 - 46215.238544938176_dp * m**5 + 311433.85378123383_dp * m**6 &
 - 1.098277765491919e6_dp * m**7 + 2.050969642246135e6_dp * m**8 - 1.591545480177363e6_dp * m**9, &
7.0536314955730575_dp - 10.540403950602274_dp * m - 595.5071040939538_dp * m**2 + 5728.24506883116_dp * m**3 &
 - 30912.229831755412_dp * m**4 + 110947.8985922118_dp * m**5 - 243377.37430261893_dp * m**6 &
 + 217457.2520679997_dp * m**7 + 205294.69579931744_dp * m**8 - 445602.4656028675_dp * m**9, &
9.033798633722743 - 16.998864167190817_dp * m - 1295.3424213765895_dp * m**2 + 15619.953032063751_dp * m**3 &
 - 105348.88568686211_dp * m**4 + 498509.8399992952_dp * m**5 - 1.6808828999144225e6_dp * m**6 &
 + 3.8035079361195057e6_dp * m**7 - 5.115225628342381e6_dp * m**8 + 3.0621926368442522e6_dp * m**9, &
11.08429290559317 - 32.673064420743835_dp * m - 2083.843283634082_dp * m**2 + 28583.011018031197_dp * m**3 &
 - 202138.62555161747_dp * m**4 + 952256.1293819118_dp * m**5 - 3.125467296989716e6_dp * m**6 &
 + 6.89027373723319e6_dp * m**7 - 9.141996456692675e6_dp * m**8 + 5.482511511780976e6_dp * m**9, &
13.245382109210354 - 60.98540199441203_dp * m - 2895.962493366059_dp * m**2 + 45001.99384872943_dp * m**3 &
 - 333306.9567541629_dp * m**4 + 1.5558493248524466e6_dp * m**5 - 4.855251960933014e6_dp * m**6 &
 + 9.923197713062515e6_dp * m**7 - 1.2079605342000252e7_dp * m**8 + 6.654728014383942e6_dp * m**9, &
14.98606645124631 - 68.82406234642987_dp * m - 4665.2218551040505_dp * m**2 + 81019.3504885875_dp * m**3 &
 - 670290.5984985809_dp * m**4 + 3.4420226099618226e6_dp * m**5 - 1.1593280790700175e7_dp * m**6 &
 + 2.5047576465778865e7_dp * m**7 - 3.1550522961975873e7_dp * m**8 + 1.762350917949359e7_dp * m**9, &
16.094573469252467 - 45.56239942638033_dp * m - 7630.1217843517325_dp * m**2 + 141149.74980783224_dp * m**3 &
 - 1.2677321434520585e6_dp * m**4 + 7.001451442234771e6_dp * m**5 - 2.503100482328536e7_dp * m**6 &
 + 5.66803432322886e7_dp * m**7 - 7.399915871224932e7_dp * m**8 + 4.2436362665615775e7_dp * m**9, &
17.73400233106048 - 68.27162470291364_dp * m - 9604.073460802441_dp * m**2 + 192142.10793340328_dp * m**3 &
 - 1.817606623855565e6_dp * m**4 + 1.038936711400024e7_dp * m**5 - 3.796565973504382e7_dp * m**6 &
 + 8.712894644888964e7_dp * m**7 - 1.1463469050932011e8_dp * m**8 + 6.5998818040929265e7_dp * m**9, &
21.41906654519393 - 234.0391923791992_dp * m - 7902.8368803977455_dp * m**2 + 193621.77945061412_dp * m**3 &
 - 1.9515809466711509e6_dp * m**4 + 1.1450662691224426e7_dp * m**5 - 4.222263288613712e7_dp * m**6 &
 + 9.69324637669229e7_dp * m**7 - 1.2703624972086267e8_dp * m**8 + 7.271365452033658e7_dp * m**9, &
27.464091085168082 - 564.4461041376098_dp * m - 1863.6807838635505_dp * m**2 + 135382.64108443772_dp * m**3 &
 - 1.5821739845253306e6_dp * m**4 + 9.729213353866134e6_dp * m**5 - 3.6300545970591106e7_dp * m**6 &
 + 8.302335752540383e7_dp * m**7 - 1.0764684237834299e8_dp * m**8 + 6.079682709171351e7_dp * m**9, &
34.575261838377614 - 982.2512909088207_dp * m + 6591.219173947654_dp * m**2 + 45014.20670022755_dp * m**3 &
 - 964680.6059437819_dp * m**4 + 6.781304774688543e6_dp * m**5 - 2.636009234925626e7_dp * m**6 &
 + 6.0583313980111904e7_dp * m**7 - 7.76000132017184e7_dp * m**8 + 4.29405781535072e7_dp * m**9, &
40.83035175410649 - 1371.9804627977408_dp * m + 14555.534924637044_dp * m**2 - 35812.195389771616_dp * m**3 &
 - 479193.8418868497_dp * m**4 + 4.884702875179248e6_dp * m**5 - 2.1297326955107525e7_dp * m**6 &
 + 5.134853205112173e7_dp * m**7 - 6.702261889484295e7_dp * m**8 + 3.721584800066794e7_dp * m**9, &
45.5284308203888 - 1693.2619756410963_dp * m + 21116.000276432842_dp * m**2 - 95526.05768687108_dp * m**3 &
 - 222322.54471269998_dp * m**4 + 4.597456446955638e6_dp * m**5 - 2.328677787804186e7_dp * m**6 &
 + 6.0688340534439355e7_dp * m**7 - 8.343550055568105e7_dp * m**8 + 4.8158642517094456e7_dp * m**9, &
50.08720719665579 - 2034.099793643092_dp * m + 28704.447883826546_dp * m**2 - 172053.32175272636_dp * m**3 &
 + 166332.40662943083_dp * m**4 + 3.75165758300126e6_dp * m**5 - 2.399366916398135e7_dp * m**6 &
 + 6.866263482050174e7_dp * m**7 - 9.970795481081429e7_dp * m**8 + 5.969903710064894e7_dp * m**9, &
56.697421294633266 - 2528.9629404500592_dp * m + 40946.40371504433_dp * m**2 - 320771.6944429163_dp * m**3 &
 + 1.2061581287027872e6_dp * m**4 - 747422.5335216314_dp * m**5 - 1.164474376882483e7_dp * m**6 &
 + 4.7507549637287274e7_dp * m**7 - 7.878579975222978e7_dp * m**8 + 5.041263316521136e7_dp * m**9, &
65.40862576538392 - 3181.0421032615573_dp * m + 58008.26463382088_dp * m**2 - 545778.3340795825_dp * m**3 &
 + 2.9445290495554246e6_dp * m**4 - 9.181634387734788e6_dp * m**5 + 1.4627071487832597e7_dp * m**6 &
 - 3.931567182260689e6_dp * m**7 - 2.0740528390033588e7_dp * m**8 + 2.145843404364027e7_dp * m**9, &
73.78012828585094 - 3839.85367122378_dp * m + 76001.49518120868_dp * m**2 - 791245.1865479778_dp * m**3 &
 + 4.882584363081332e6_dp * m**4 - 1.8636080358183626e7_dp * m**5 + 4.362274991378952e7_dp * m**6 &
 - 5.8377583022248864e7_dp * m**7 + 3.630303932702363e7_dp * m**8 - 3.92404406080585e6_dp * m**9, &
80.78651201975016 - 4444.375575395826_dp * m + 93487.98540601635_dp * m**2 - 1.0396821474621181e6_dp * m**3 &
 + 6.899362586966885e6_dp * m**4 - 2.862101523532484e7_dp * m**5 + 7.42341862346831e7_dp * m**6 &
 - 1.1474555998435119e8_dp * m**7 + 9.269829831882124e7_dp * m**8 - 2.690428909946104e7_dp * m**9, &
88.55370059514114 - 5131.085700314765_dp * m + 114344.02882752768_dp * m**2 - 1.3544692437221785e6_dp * m**3 &
 + 9.64157609806654e6_dp * m**4 - 4.336246589549065e7_dp * m**5 + 1.2409333887802215e8_dp * m**6 &
 - 2.1840231339437953e8_dp * m**7 + 2.140096096516368e8_dp * m**8 - 8.82051622661604e7_dp * m**9, &
99.25874856504055 - 6037.519212433856_dp * m + 142400.46902943595_dp * m**2 - 1.7972535347020212e6_dp * m**3 &
 + 1.373134548906169e7_dp * m**4 - 6.689624311518229e7_dp * m**5 + 2.099414548977667e8_dp * m**6 &
 - 4.1215479983808345e8_dp * m**7 + 4.6160982473789173e8_dp * m**8 - 2.2555121324761388e8_dp * m**9, &
111.22343849143047 - 7057.7883518674935_dp * m + 174876.65256694023_dp * m**2 - 2.3285236853465238e6_dp * m**3 &
 + 1.8831425043508448e7_dp * m**4 - 9.740607027364239e7_dp * m**5 + 3.2551457162271464e8_dp * m**6 &
 - 6.824650525648329e8_dp * m**7 + 8.186372040669223e8_dp * m**8 - 4.2962144509474057e8_dp * m**9, &
118.40173153696084 - 7812.9529878720605_dp * m + 201655.41435957828_dp * m**2 - 2.7976860482022595e6_dp * m**3 &
 + 2.3562208320128188e7_dp * m**4 - 1.2678881464073764e8_dp * m**5 + 4.401435431626159e8_dp * m**6 &
 - 9.568670902700018e8_dp * m**7 + 1.1877107147181227e9_dp * m**8 - 6.435275293916061e8_dp * m**9, &
115.23780058412257 - 7954.098547423372_dp * m + 213401.56075989947_dp * m**2 - 3.0656847210683743e6_dp * m**3 &
 + 2.6652992787059106e7_dp * m**4 - 1.4763398125921527e8_dp * m**5 + 5.261297100928776e8_dp * m**6 &
 - 1.171089728692452e9_dp * m**7 + 1.4844693347257988e9_dp * m**8 - 8.193667886327482e8_dp * m**9, &
101.07640784409782 - 7428.426895860094_dp * m + 208435.16419515907_dp * m**2 - 3.104862494686558e6_dp * m**3 &
 + 2.7841258940538246e7_dp * m**4 - 1.5843670065714395e8_dp * m**5 + 5.782367580166419e8_dp * m**6 &
 - 1.3144959937357693e9_dp * m**7 + 1.697654743868332e9_dp * m**8 - 9.526442974725475e8_dp * m**9, &
77.99285047093665 - 6338.071621669256_dp * m + 188774.5971540118_dp * m**2 - 2.936203149088958e6_dp * m**3 &
 + 2.7251808025285892e7_dp * m**4 - 1.5962241949398154e8_dp * m**5 + 5.972102778005555e8_dp * m**6 &
 - 1.387410470632458e9_dp * m**7 + 1.8264524172120025e9_dp * m**8 - 1.0424803191139958e9_dp * m**9, &
46.812060301731634 - 4704.414147921791_dp * m + 154161.01765534276_dp * m**2 - 2.543931572120198e6_dp * m**3 &
 + 2.4637648870198518e7_dp * m**4 - 1.4918439690209216e8_dp * m**5 + 5.735316641460102e8_dp * m**6 &
 - 1.3632338571666594e9_dp * m**7 + 1.830188293678673e9_dp * m**8 - 1.0625768005994213e9_dp * m**9, &
7.847139312252604 - 2526.8511472074597_dp * m + 103980.69360599533_dp * m**2 - 1.9097785043032581e6_dp * m**3 &
 + 1.9749925135998152e7_dp * m**4 - 1.2520905046062276e8_dp * m**5 + 4.983813858268824e8_dp * m**6 &
 - 1.2177673740966372e9_dp * m**7 + 1.6723574047458694e9_dp * m**8 - 9.895910016077784e8_dp * m**9, &
-36.60844348195546 + 64.38963022310429_dp * m + 41282.822897670965_dp * m**2 - 1.0725075574993882e6_dp * m**3 &
 + 1.288503276288955e7_dp * m**4 - 8.911418200815153e7_dp * m**5 + 3.7600165709962046e8_dp * m**6 &
 - 9.586239780373292e8_dp * m**7 + 1.36027155172798e9_dp * m**8 - 8.262629817584653e8_dp * m**9, &
-83.65415075660984 + 2898.2558327002844_dp * m - 29739.586432201475_dp * m**2 - 88425.4522298853_dp * m**3 &
 + 4.50114243939258e6_dp * m**4 - 4.325052026927162e7_dp * m**5 + 2.14011580197559e8_dp * m**6 &
 - 6.008343895710565e8_dp * m**7 + 9.102272112463979e8_dp * m**8 - 5.799037097249833e8_dp * m**9, &
-131.45151770796767 + 5864.9109097763785_dp * m - 106397.38104899065_dp * m**2 + 1.0069247724963345e6_dp * m**3 &
 - 5.12222495824354e6_dp * m**4 + 1.1028894311098354e7_dp * m**5 + 1.6403335728605384e7_dp * m**6 &
 - 1.5114413255285746e8_dp * m**7 + 3.2773396794217116e8_dp * m**8 - 2.5173318128572628e8_dp * m**9, &
-177.01970592171443 + 8778.893094042838_dp * m - 183842.38347962237_dp * m**2 + 2.1433427375289723e6_dp * m**3 &
 - 1.5360872436360464e7_dp * m**4 + 7.016835088510215e7_dp * m**5 - 2.0377307680117235e8_dp * m**6 &
 + 3.605330889556964e8_dp * m**7 - 3.4816066470694965e8_dp * m**8 + 1.3604082855689755e8_dp * m**9, &
-215.78048538318092 + 11352.825876924011_dp * m - 254429.49287131111_dp * m**2 + 3.207491974753509e6_dp * m**3 &
 - 2.5178624444914337e7_dp * m**4 + 1.2808400346707451e8_dp * m**5 - 4.2348438074164295e8_dp * m**6 &
 + 8.797702668176477e8_dp * m**7 - 1.0443713913137125e9_dp * m**8 + 5.408018471360576e8_dp * m**9, &
-243.58444537069815 + 13318.522877293628_dp * m - 310823.8967104135_dp * m**2 + 4.0877750347270467e6_dp * m**3 &
 - 3.353054768043964e7_dp * m**4 + 1.785061108936398e8_dp * m**5 - 6.185295733576396e8_dp * m**6 &
 + 1.348414258784031e9_dp * m**7 - 1.6817205086301703e9_dp * m**8 + 9.15886398640447e8_dp * m**9, &
-257.39444091154286 + 14468.670181564157_dp * m - 347087.6302364655_dp * m**2 + 4.690172468461062e6_dp * m**3 &
 - 3.950564922252864e7_dp * m**4 + 2.1580492206352308e8_dp * m**5 - 7.666172459932783e8_dp * m**6 &
 + 1.7117139859169667e9_dp * m**7 - 2.1842549100865316e9_dp * m**8 + 1.2158166643726773e9_dp * m**9, &
-255.9203830507791 + 14704.52505926997_dp * m - 360168.1882953858_dp * m**2 + 4.963792794528394e6_dp * m**3 &
 - 4.259370541379616e7_dp * m**4 + 2.367583889392443e8_dp * m**5 - 8.548276507702368e8_dp * m**6 &
 + 1.9377415670979712e9_dp * m**7 - 2.5076084922317166e9_dp * m**8 + 1.4140578927492626e9_dp * m**9, &
-240.15330516183064 + 14070.28769389826_dp * m - 350829.8409094997_dp * m**2 + 4.914851726834618e6_dp * m**3 &
 - 4.2812846698245235e7_dp * m**4 + 2.4128669838462934e8_dp * m**5 - 8.822829720294756e8_dp * m**6 &
 + 2.0233060607547004e9_dp * m**7 - 2.6462405533555937e9_dp * m**8 + 1.5067437610278237e9_dp * m**9, &
-213.1902407380797 + 12738.93413893343_dp * m - 323182.79604288685_dp * m**2 + 4.598161669252364e6_dp * m**3 &
 - 4.061688939100559e7_dp * m**4 + 2.3181894019177043e8_dp * m**5 - 8.574207588201002e8_dp * m**6 &
 + 1.9868061404940634e9_dp * m**7 - 2.6230656524502583e9_dp * m**8 + 1.5063226330792553e9_dp * m**9, &
-178.28284526587058 + 10894.22349973113_dp * m - 281637.1554072727_dp * m**2 + 4.0729883935745903e6_dp * m**3 &
 - 3.64994470684883e7_dp * m**4 + 2.110103705069777e8_dp * m**5 - 7.895126105205859e8_dp * m**6 &
 + 1.8486109981118639e9_dp * m**7 - 2.463766069494559e9_dp * m**8 + 1.4270264140040352e9_dp * m**9 ]

  end function LagCoef

!ccccccccccccccc

  recursive pure real (dp) function BreitWigner(Gamma, i, l, l2) result(res)
    real (dp)          , intent(in)  :: Gamma, l
    real (dp), optional, intent(in)  :: l2
    integer            , intent(in)  :: i
    real (dp)                        :: Gammal

    res = 0

    if ( present(l2) ) then
      if (l2 > l) res = BreitWigner(Gamma, i, l2) - BreitWigner(Gamma, i, l)
      return
    end if

    Gammal = Gamma**2 + l**2

    select case(i)
    case default; res = 0
    case(0);  res = Gamma/Gammal/Pi
    case(1);  res = - 2 * Gamma * l/Gammal**2/Pi
    case(2);  res = - 2 * Gamma * (Gamma**2 - 3 * l**2)/Gammal**3/Pi
    case(3);  res = - 24 * Gamma * l * (Gamma**2 - l**2)/Gammal**4/Pi
    case(-1); res = 0.5_dp + ATan(l/Gamma)/Pi
    case(-2); res = ( l/2 + Gamma * log(Gammal)/Pi/2 + l/Pi * ATan(l/Gamma) )/Gammal
    case(-3); res = log(Gammal) * ( 0.5_dp + ATan(l/Gamma)/Pi )/2
    end select

  end function BreitWigner

!ccccccccccccccc

end module MCtopClass
