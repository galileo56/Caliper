
! ccccccccccc

module RNRQCDClass
  use constants, only: dp, Pi, Pi2, l2, Zeta3, Euler ; use ToppikClass
  use DeriGamma, only: DiGam ; use RunningClass; use AnomDimClass
  use AlphaClass; implicit none ;  private

  real (dp), parameter :: Qt = 2._dp/3

  public :: qFromV, vC, VcsLL, XiNNLLSoftMixLogc1, vStar, vRootStar, SwitchOff

! ccccccccccc

  type, public :: RNRQCD
    integer                   :: nl
    real (dp)                 :: a1, a2, mass
    real (dp), dimension(0:4) :: beta
    type (Running)            :: run
    type (Alpha)              :: AlphaOb

  contains

    procedure, pass(self), public :: VceffsNNLL, VrsLL, V2sLL, VssLL, Vk1sLL,  &
    Vk2sLL, VkeffsLL, XiNLL, XiNNLLnonmix, XiNNLLmixUsoft, MLLc2, MNLLc1, &
    MNLLplusNNLLnonmixc1, MNNLLAllc1InclSoftMixLog, A1pole, NRQCD

    procedure, pass(self), private :: xc01, xc11, xc12, xc22

  end type RNRQCD

!ccccccccccccccc

  interface RNRQCD
    module procedure RNRQCDIn
  end interface RNRQCD

contains

  type (RNRQCD) function RNRQCDIn(run)
    type (Running), intent(in) :: run
    type (AnomDim)             :: adim
    integer :: nl

    nl = run%numFlav(); adim = run%adim(); RNRQCDIn%run = run; RNRQCDIn%nl = nl
    RNRQCDIn%beta = adim%betaQCD('beta') ; RNRQCDIn%AlphaOb = run%alphaAll()
    RNRQCDIn%mass = run%scales('mH')

    RNRQCDIn%a1 = 31._dp/3 - 10._dp * nl/9
    RNRQCDIn%a2 = 456.74883699902244_dp - 66.35417150661816_dp * nl + &
    1.2345679012345678_dp * nl**2

  end function RNRQCDIn

! ccccccccccc

  real (dp) function NRQCD(self, order, scheme, q, gt, h, nu)
    class (RNRQCD)     , intent(in) :: self
    character (len = *), intent(in) :: order, scheme
    real (dp)          , intent(in) :: gt, h, nu, q
    real (dp)                       :: ac, ah, asLL
    complex (dp)                    :: v

    NRQCD = 0

    if ( scheme(:4) == 'pole') then

      if ( order(:2) == 'LL' ) then

        ah = self%run%AlphaQCD(h * self%mass); v = vC(q, self%mass, gt)

        asLL = self%AlphaOb%alphaGenericFlavor('inverse', 1, self%nl, &
        h * self%mass, ah, h * nu * self%mass)

        ac = 4 * asLL/3

        NRQCD =  18 * (Qt * self%mass/q)**2 *  ImagPart(  (0,1) * v - ac * (  &
        Euler - 0.5_dp + l2 + Log( - (0,1) * v/(h * nu) ) + &
        DiGam( 1 - (0,1) * ac/(2 * v) )  )  )

      end if

    end if

  end function NRQCD

! ccccccccccc

  real (dp) function A1pole(self, order, En, gamtop, asoft, VcsNNLL, musoft)
    class (RNRQCD)     , intent(in) :: self
    character (len = *), intent(in) :: order
    real (dp)          , intent(in) :: En, gamtop, asoft, musoft, VcsNNLL
    type (Toppik)                   :: TP
    real (dp), dimension(2)         :: res
    real (dp)                       :: x0, x1, x2, aSuPi

    x0 = 0; x1 = 0; x2 = 0

    if ( order(:2) == 'LL' ) then
      x0 = 1
    else if ( order(:3) == 'NLL' ) then
      aSuPi = asoft/4/Pi; x0 = self%xc01(aSuPi); x1 = self%xc11(aSuPi)
    else if ( order(:4) == 'NNLL' .or. order(:4) == 'N2LL' ) then
      aSuPi = asoft/4/Pi; x1 = self%xc12(aSuPi); x2 = self%xc22(aSuPi)
      x0 = self%a1 * asuPi + self%a2 * asuPi**2 - 3 * VcsNNLL/16/asoft/Pi
    end if

    TP = Toppik(self%nl, En, self%mass, gamtop, asoft, musoft, 175000000._dp, &
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

    VssLL = - 8 * ah * Pi/(6 - self%beta(0)) * ( 1 - (as/ah)**(1 - 6/self%beta(0)) &
    * (21 - 2 * self%beta(0))/9 )

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

    MLLc2 = -1._dp/6 - 32 * log(au/ah)/9/self%beta(0)

  end function MLLc2

! ccccccccccc

  real (dp) function MNLLc1(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au

    MNLLc1 = (1 - 8 * ah/3/Pi) * exp( self%xiNLL(ah, as, au) )

  end function MNLLc1

! ccccccccccc

  real (dp) function MNLLplusNNLLnonmixc1(self, ah, as, au)
    class (RNRQCD), intent(in) :: self
    real (dp)     , intent(in) :: ah, as, au

    ! MNLLplusNNLLnonmixc1 = 1 - 8 * ah/3/Pi + ah**2 * (   4 * (l2/2 - 0.625_dp) &
    ! + 16 * (l2/3 - 2/Pi**2 - 31._dp/24)/9 + (  8 * (11/Pi**2 - 2)/27 + 22 *    &
    ! self%nl/Pi2/27 + 16 * ( 4 * l2/3 + (9.75_dp - Zeta3)/Pi**2 - 35._dp/18 )/9 &
    ! - 4 * ( 8 * l2/3 - 179._dp/72 + (151._dp/36 + 13 * Zeta3/2)/Pi**2 )  )/2   )

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

  real (dp) function qFromV(v, m, gt)
    real (dp), intent(in) :: v, m, gt

    qFromV = (8 * m**2 * v**2 + 4 * m**2 * v**4 - gt**2)/4/m/v**2

  end function qFromV

! ccccccccccc

  complex (dp) function vC(q, m, gt)
    real (dp), intent(in) :: q, m, gt

    vC = sqrt(  ( q - 2 * m + (0,1) * gt ) /m  )

  end function vC

! ccccccccccc

  real (dp) function vStar(q, m, gt)
    real (dp), intent(in) :: q, m, gt

    vStar = 0.05_dp + Abs( vC(q, m, gt) )

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

    SwitchOff = 0; v = realpart( vC(q, m, gt) )

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

end module RNRQCDClass
