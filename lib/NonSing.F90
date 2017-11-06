
module NonSingularClass
  use RunningClass;  use ModelClass  ; use Constants, only: dp, Pi, Pi2, Zeta3
  use AnomDimClass;  implicit none   ; private

  public :: NonSingular, ThrustNS1loop

!ccccccccccccccc

  type NonSingular

    integer                  , private     :: n2, n3, nf, isoft
    character (len = 6)      , private     :: shape
    type (Running)           , private     :: alphaMass
    type (AnomDim)           , private     :: andim
    real (dp), dimension(0:3), private     :: beta
    real (dp)                , private     :: s2, s2rho, s3
    real (dp), dimension(:,:), allocatable, private :: interList2loop, interList3loop
    real (dp), dimension(:)  , allocatable, private :: interDer2loop, interDer3loop

  contains

    ! final                                  :: delete_object

    procedure :: HJMNS2loopFit, ThrustNS3loopFit, HJMNS3loopFit, Interpol2loop, &
    ThrustSing2loop, Interpol3loop, NS2loop, ThrustSing3loop

  end type NonSingular

!ccccccccccccccc

  interface NonSingular
    module procedure InitNonSingular
  end interface NonSingular

  contains

!ccccccccccccccc

  ! subroutine delete_object(this)
  !  type(NonSingular) :: this
  !    if ( allocated(this%interList2loop) ) deallocate(this%interList2loop)
  !    if ( allocated(this%interList3loop) ) deallocate(this%interList3loop)
  !    if ( allocated(this%interDer2loop ) ) deallocate(this%interDer2loop )
  !    if ( allocated(this%interDer3loop ) ) deallocate(this%interDer3loop )
  !
  ! end subroutine delete_object

!ccccccccccccccc

   type(NonSingular) function InitNonSingular(shape, alphaMass, s3, j3, isoft)
    character (len = *), intent(in) :: shape
    integer            , intent(in) :: isoft
    type(Running)      , intent(in) :: alphaMass
    real (dp)          , intent(in) :: s3, j3

    real (dp), dimension(102,2) :: interListThrust2loop
    real (dp), dimension(153,2) :: interListHJM2loop
    real (dp), dimension(18,2)  :: interListHJM3loop
    real (dp), dimension(27,2)  :: interListThrust3loop
    real (dp)                   :: s2rho, s2
    integer                     :: nf
    type(AnomDim)               :: andim

    andim = alphaMass%adim();         InitNonSingular%andim = andim; nf = andim%numFlav();
    InitNonSingular%s3 = s3 + j3;     InitNonSingular%alphaMass = alphaMass
    InitNonSingular%shape = shape;    InitNonSingular%s2 = s2;
    InitNonSingular%isoft = isoft;    InitNonSingular%nf = nf

    select case(isoft)
      case(-1)
       s2rho = 5.817594867359059_dp + 0.2895745232754116_dp * nf
    case(1)
       s2rho = 6.234083201979812_dp + 0.2741556778080381_dp * nf
    case(2)
       s2rho = 5.873308790742764_dp + 0.28918794494291505_dp * nf
    case(3)
       s2rho = 5.825832781183491_dp + 0.2898944331804042_dp * nf
    case(4)
       s2rho = 5.8188363048966805_dp + 0.2897145802589072_dp * nf
    case(5)
       s2rho = 5.817775727311131_dp + 0.28961899917968487_dp * nf
    case(6)
       s2rho = 5.817618283362493_dp + 0.28958722156747796_dp * nf
    case default
       s2rho = 0
    end select

    s2rho = s2rho/4 - 40.68040549378476_dp

    InitNonSingular%s2rho = s2rho

    interListThrust2loop(:,1) = [ 0.09, 0.0925, 0.095, 0.0975, 0.1, 0.105, 0.115,0.125, &
    0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195, 0.1975, 0.2025, 0.205, 0.2075,     &
    0.2125, 0.2175, 0.2225, 0.2275, 0.2325, 0.2375, 0.2425, 0.2475, 0.2525, 0.2575,     &
    0.2625, 0.2675, 0.2725, 0.2775, 0.2825, 0.2875, 0.2925, 0.2975, 0.3025, 0.3075,     &
    0.309, 0.311, 0.3125, 0.315, 0.317, 0.319, 0.321, 0.323, 0.325, 0.327, 0.329, 0.331,&
    0.333, 0.335, 0.337, 0.339, 0.341, 0.343, 0.345, 0.347, 0.349, 0.351, 0.353, 0.355, &
    0.357, 0.359, 0.361, 0.363, 0.365, 0.367, 0.369, 0.371, 0.373, 0.375, 0.377, 0.379, &
    0.381, 0.383, 0.385, 0.387, 0.389, 0.391, 0.393, 0.395, 0.397, 0.399, 0.401, 0.403, &
    0.405, 0.407, 0.409, 0.411, 0.413, 0.415, 0.417, 0.419, 0.421, 0.423, 0.425, 0.427, &
    0.429, 0.431,  0.433 ]

    interListThrust2loop(:,2) = [ -232.794, -230.591, -228.4, -226.219, -224.051,     &
    -219.417, -211.09, -202.831, -194.936, -187.417, -180.159, -173.596, -166.684,    &
    -160.596, -154.936, -153.665, -151.025, -149.615, -148.204, -145.959, -143.33,    &
    -141.248, -139.188, -137.112, -135.198, -133.697, -131.886, -130.66, -129.375,    &
    -128.213, -127.551, -126.797, -126.439, -126.478, -126.795, -127.77, -129.326,    &
    -131.845, -134.421, -135.257, -136.457, -137.349, -138.878, -139.966, -141.179,   &
    -142.331, -143.288, -144.388, -145.359, -146.153, -146.49, -149.242, -167.842,    &
    -177.966, -183.189, -185.961, -187.224, -187.601, -187.365, -186.653, -185.597,   &
    -184.271, -182.725, -181.022, -179.183, -177.232, -175.192, -173.085, -170.926,   &
    -168.714, -166.468, -164.206, -161.924, -159.624, -157.322, -155.014, -152.71,    &
    -150.409, -148.117, -145.836, -143.563, -141.308, -139.068, -136.843, -134.638,   &
    -132.45, -130.285, -128.14, -126.017, -123.917, -121.838, -119.783, -117.75,      &
    -115.739,-113.751, -111.784, -109.839, -107.92, -106.029, -104.164, -102.326,     &
    -100.514 ]

    interListHJM2loop(:,1) = [ 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0665, 0.0735, 0.0805, &
    0.0875, 0.0945, 0.1015, 0.1085, 0.1155, 0.1225, 0.1295, 0.1365, 0.1435, 0.1505,   &
    0.1575, 0.1645, 0.1715, 0.1785, 0.1855, 0.1925, 0.1995, 0.2065, 0.2135, 0.2205,   &
    0.2275, 0.2325, 0.2345, 0.2375, 0.2425, 0.2475, 0.2525, 0.2575, 0.2625, 0.2675,   &
    0.2725, 0.2775, 0.2825, 0.2875, 0.2925, 0.2975, 0.3025, 0.3075, 0.3125, 0.3175,   &
    0.3195, 0.3205, 0.3225, 0.3235, 0.3245, 0.3255, 0.3265, 0.3275, 0.3285, 0.3295,   &
    0.3305, 0.3315, 0.3325, 0.3335, 0.3345, 0.3355, 0.3365, 0.3375, 0.3385, 0.3395,   &
    0.3405, 0.3415, 0.3425, 0.3435, 0.3445, 0.3455, 0.3465, 0.3475, 0.3485, 0.3495,   &
    0.3505, 0.3515, 0.3525, 0.3535, 0.3545, 0.3555, 0.3565, 0.3575, 0.3585, 0.3595,   &
    0.3605, 0.3615, 0.3625, 0.3635, 0.3645, 0.3655, 0.3665, 0.3675, 0.3685, 0.3695,   &
    0.3705, 0.3715, 0.3725, 0.3735, 0.3745, 0.3755, 0.3765, 0.3775, 0.3785, 0.3795,   &
    0.3805, 0.3815, 0.3825, 0.3835, 0.3845, 0.3855, 0.3865, 0.3875, 0.3885, 0.3895,   &
    0.3905, 0.3915, 0.3925, 0.3935, 0.3945, 0.3955, 0.3965, 0.3975, 0.3985, 0.3995,   &
    0.4005, 0.4015, 0.4025, 0.4035, 0.4045, 0.4055, 0.4065, 0.4075, 0.4085, 0.4095,   &
    0.4105, 0.4115, 0.4125, 0.4135, 0.4145, 0.4155, 0.4165, 0.4175, 0.4185, 0.4195,   &
    0.4205, 0.4215, 0.4225, 0.4235, 0.4245 ]

    interListHJM2loop(:,2) = [ -420.262_dp + 0.0153806_dp * s2rho, -415.979_dp + 0.0234176_dp * s2rho, &
    -411.619_dp + 0.027652_dp * s2rho, -407.201_dp + 0.0285161_dp * s2rho, -402.742_dp + 0.0265489_dp * s2rho, &
    -389.161_dp, -378.498_dp, -366.477_dp, -355.161_dp, -344.874_dp, -335.065_dp, -325.875_dp, -317.187_dp, &
    -310.952_dp, -304.802_dp, -301.691_dp, -299.798_dp, -297.908_dp, -294.811_dp, -292.13_dp, -288.186_dp,  &
    -283.905_dp, -279.537_dp, -274.846_dp, -270.328_dp, -265.407_dp, -260.917_dp, -256.121_dp, -251.554_dp, &
    -248.507_dp, -247.161_dp, -245.238_dp, -242.018_dp, -239.051_dp, -236.052_dp, -232.753_dp, -229.76_dp,  &
    -226.463_dp, -223.233_dp, -219.838_dp, -216.076_dp, -212.252_dp, -207.995_dp, -203.323_dp, -198.008_dp, &
    -191.961_dp, -184.662_dp, -176.005_dp, -172.01_dp, -169.789_dp, -165.178_dp, -162.584_dp, -159.9_dp,    &
    -157.031_dp, -153.918_dp, -150.454_dp, -146.838_dp, -142.825_dp, -138.301_dp, -132.871_dp, -126.16_dp,  &
    -121.066_dp, -126.62_dp, -129.75_dp, -131.613_dp, -132.749_dp, -133.414_dp, -133.683_dp, -133.745_dp,   &
    -133.564_dp, -133.271_dp, -132.844_dp, -132.329_dp, -131.735_dp, -131.071_dp, -130.368_dp, -129.631_dp, &
    -128.826_dp, -128.016_dp, -127.167_dp, -126.292_dp, -125.407_dp, -124.486_dp, -123.569_dp, -122.628_dp, &
    -121.674_dp, -120.721_dp, -119.747_dp, -118.776_dp, -117.802_dp, -116.818_dp, -115.836_dp, -114.847_dp, &
    -113.858_dp, -112.867_dp, -111.877_dp, -110.891_dp, -109.896_dp, -108.907_dp, -107.924_dp, -106.937_dp, &
    -105.957_dp, -104.976_dp, -104.001_dp, -103.028_dp, -102.059_dp, -101.094_dp, -100.131_dp, -99.1743_dp, &
    -98.2215_dp, -97.272_dp, -96.3284_dp, -95.39_dp, -94.4544_dp, -93.5233_dp, -92.5993_dp, -91.6809_dp,    &
    -90.7635_dp, -89.8546_dp, -88.9509_dp, -88.0519_dp, -87.1585_dp, -86.2695_dp, -85.3882_dp, -84.5112_dp, &
    -83.639_dp, -82.7746_dp, -81.9142_dp, -81.0586_dp, -80.2094_dp, -79.3662_dp, -78.528_dp, -77.6958_dp,   &
    -76.8696_dp, -76.0488_dp, -75.2331_dp, -74.4236_dp, -73.6197_dp, -72.8214_dp, -72.0287_dp, -71.2418_dp, &
    -70.4602_dp, -69.6845_dp, -68.9143_dp, -68.1498_dp, -67.3908_dp, -66.6375_dp, -65.8898_dp, -65.1478_dp, &
    -64.4113_dp, -63.6803_dp, -62.9549_dp, -62.2349_dp ]

     interListHJM3loop(:,1) = [ 0.265, 0.27, 0.275, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, &
     0.41, 0.43, 0.45, 0.47, 0.49, 0.5, 0.505, 0.51, 0.515 ]

     interListHJM3loop(:,1) = [ - 4695.24_dp - 0.248188_dp * s3 + 6.28379_dp * s2 + &
     13.2304_dp * s2rho, -4562.34_dp - 0.255086_dp * s3 + 6.0601_dp * s2 + 12.8604_dp* s2rho, &
     - 4432.38_dp -  0.261031_dp * s3 + 5.84463_dp * s2 + 12.4979_dp * s2rho, -4202.26_dp + &
     6.72232_dp * s2 + 15.2315_dp * s2rho, -3212.76_dp + 6.28831_dp * s2 + 14.8226_dp * s2rho, &
     -1696.15_dp + 5.90696_dp * s2 + 14.4296_dp * s2rho, -1047.16_dp + 5.56923_dp * s2 + &
     14.0534_dp * s2rho, -781.843_dp + 5.26804_dp * s2 + 13.6943_dp * s2rho, -430.778_dp + &
     4.99777_dp * s2 + 13.352_dp * s2rho, -129.632_dp + 4.75387_dp * s2 + 13.0259_dp * s2rho, &
     124.271_dp + 4.53268_dp * s2 + 12.7154_dp * s2rho, 334.432_dp + 4.33116_dp * s2 + &
     12.4197_dp * s2rho, 509.709_dp + 4.1468_dp * s2 + 12.1379_dp * s2rho, 655.781_dp + &
     3.97749_dp * s2 + 11.8693_dp * s2rho, 720.287_dp + 3.8974_dp * s2 + 11.7391_dp * s2rho, &
     749.895_dp + 3.85881_dp * s2 + 11.6755_dp * s2rho, 778.156_dp + 3.82098_dp * s2 + &
     11.6125_dp * s2rho, 805.127_dp + 3.78389_dp * s2 + 11.5503_dp * s2rho ]

     interListThrust3loop(:,1) = [ 0.3, 0.305, 0.31, 0.3125, 0.315, 0.325, 0.335, 0.345,&
     0.355, 0.365, 0.375, 0.385, 0.395, 0.405, 0.415, 0.425, 0.435, 0.445, 0.455, 0.465, &
     0.475, 0.485, 0.495, 0.5, 0.505, 0.51, 0.515 ]

     interListThrust3loop(:,2) = [ -1699.87_dp + 0.0000538264_dp * s3, -1587.43_dp     &
     -0.00036209_dp * s3, -1474.77_dp - 0.000799743_dp * s3, -1418.39_dp             &
     -0.00102629_dp * s3, -1362._dp - 0.00125775_dp * s3, -1552.12_dp, -2179.12_dp,   &
     -3559.86_dp, -4026.42_dp, -4056.19_dp, -3951.16_dp, -3669.09_dp, -3358.12_dp, -3075.06_dp, &
     -2808.38_dp , -2561.97_dp, -2334.96_dp, -2124.33_dp, -1928.75_dp, -1747.01_dp, -1578.01_dp,&
     -1420.75_dp, -1274.35_dp, -1204.55_dp, -1137.57_dp, -1072.92_dp, -1010.49_dp ]

     if ( shape(:6) == 'thrust' ) then

       InitNonSingular%n2 = 102;  allocate( InitNonSingular%interList2loop(102,2)  )
       InitNonSingular%n3 =  27;  allocate( InitNonSingular%interList3loop( 27,2)  )
       allocate( InitNonSingular%interDer2loop(102), InitNonSingular%interDer3loop(27)  )

       InitNonSingular%interList2loop = interListThrust2loop
       InitNonSingular%interList3loop = interListThrust2loop

     else if ( shape(:3) == 'HJM' ) then

       InitNonSingular%n2 = 153;  allocate( InitNonSingular%interList2loop(153,2)  )
       InitNonSingular%n3 =  18;  allocate( InitNonSingular%interList3loop( 18,2)  )
       allocate( InitNonSingular%interDer2loop(153), InitNonSingular%interDer3loop(18)  )

       InitNonSingular%interList2loop = interListHJM2loop
       InitNonSingular%interList3loop = interListHJM3loop

     end if

     call spline_cubic_set( InitNonSingular%n2, InitNonSingular%interList2loop(:,1), &
     InitNonSingular%interList2loop(:,2), 0, 0._dp, 0, 0._dp, InitNonSingular%interDer2loop )

     call spline_cubic_set ( InitNonSingular%n3, InitNonSingular%interList3loop(:,1), &
     InitNonSingular%interList3loop(:,2), 0, 0._dp, 0, 0._dp, InitNonSingular%interDer3loop )

   end function InitNonSingular

!ccccccccccccccc

  function NS2loop(self, er, t) result(beta)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: er, t
    real (dp), dimension(2) :: beta

    beta = 0

    if (t <= 0) then

      beta = 0

    else if (t < 0.095) then

      beta = ThrustNS2loopFit(t) + er * ThrustNS2loopErr(t)

    else if (t < 0.433) then

      beta = self%Interpol2loop(t)

    else

      beta = - self%ThrustSing2loop(t)

    endif

  end function NS2loop

!ccccccccccccccc

  real (dp) function NS3loop(self, er, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in)    :: er, t

    if (t <= 0) then

      NS3loop = 0

    else if (t < 0.315) then

      NS3loop = self%ThrustNS3loopFit(t) + er * ThrustNS3loopErr(t)

    else if (t < 0.515) then

      NS3loop = self%Interpol3loop(t)

    else

      NS3loop = - self%ThrustSing3loop(t)

    endif
  end

!ccccccccccccccc

  function Interpol2loop(self, t) result(beta)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp), dimension(2) :: beta
    real (dp) :: aux

    if ( t > self%interList2loop(self%n2,1) .or. t < self%interList2loop(1,1) ) then

      beta = 0

    else

      call spline_cubic_val( self%n2, self%interList2loop(:,1), self%interList2loop(:,2), &
      self%interDer2loop, t, beta(1), beta(2), aux )

    end if

  end function Interpol2loop

!ccccccccccccccc

  real (dp) function Interpol3loop(self, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp) :: aux

    if ( t > self%interList3loop(1,1) .or. t < self%interList3loop(1,self%n2) ) then

      Interpol3loop = 0

    else

      call spline_cubic_val( self%n2, self%interList3loop(1,:), self%interList3loop(2,:), &
      self%interDer3loop, t, Interpol3loop, aux, aux )

    end if

  end function Interpol3loop

!ccccccccccccccc

  real (dp) function ThrustReg1loop(t)
    real (dp), intent(in) :: t
    real (dp) :: lt, lt3, lt15

    lt = log(t);  lt3 = 3 * lt;  lt15 = 15 * lt

    ThrustReg1loop = 4 * (lt - 1)/3 - 2 * t * (4 * lt - 5)/3 - 8 * t**2 * (11 + lt3)/9 - &
    8 * t**3 * (20 + lt3)/9 - 16 * t**4 * (166 + lt15)/90 - 16 * t**5 * (272 + lt15)/90

  end function ThrustReg1loop

!ccccccccccccccc

  real (dp) function ThrustFO1loop(t)
    real (dp), intent(in) :: t

    real (dp)             :: t2

    t2 = t**2

    ThrustFO1loop = ( 4/t/3 * (3 - 9 * t + 9 * t**3 + log(1/t - 2) * &
    (6 * t - 4 - 6 * t2) - 3 * t2) )/(t - 1)/2

  end function ThrustFO1loop

!ccccccccccccccc

  real (dp) function ThrustSing1loop(t)
    real (dp), intent(in) :: t

      ThrustSing1loop = - 2 * ( 1 + 4 * log(t)/3 )/t

  end function ThrustSing1loop

!ccccccccccccccc

  function ThrustSing2loop(self, t) result(beta)
    class(NonSingular), intent(in) :: self
    real (dp)   , intent(in) :: t

    real (dp), dimension(2) :: beta, o
    real (dp)               :: lg

    lg = log(t)

    beta(1) = - 12.673650220838139635048946729512 + 5._dp/6 * self%nf            &
    + ( - 14.5773382234326209783148442511447 + 11._dp/27 * self%nf ) * lg        &
    + (19 - 2._dp/3 * self%nf) * lg**2 + 32 * lg**3/9

    beta(1) = beta(1)/t


    o(1) = t**2;  o(2) = 1/o(1)

    beta(2) = ( -16.1332705288964373124827034189366 * o(2) - 128 * o(2) * lg**3/9 + &
              175.494538078915662815404630237026 * o(2) * lg - 20 * o(2) * lg**2 )/4

  end function ThrustSing2loop

!ccccccccccccccc

  real (dp) function ThrustSing3loop(self, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t

    integer :: nf2
    real (dp) :: lgt

    nf2 = self%nf**2;  lgt = log(t)

    ThrustSing3loop = ( - 143.525152363904551222617556049954 - 64 * lgt**5/27 &
    + 12.0822164629231987298396688856883 * self%nf - 6.30586823940581808756178361363709d-2 * nf2 &
    + (-122.613065878586624535273585934192 + 2.42690837988652807766243313380983 * self%nf &
    + 3.15710396734716702837886259658262d-1 * nf2) * lgt + (307.891498297039767351179762044922 &
    - 25.9452771711369534202162867586594 * self%nf + 17._dp * nf2/54) * lgt**2 &
    + (- 42.9438225359046032991727770422585 + 6.39506172839506170646473037777469 * self%nf &
    - 1.72839506172839496578319540276425d-1 * nf2) * lgt**3 + (40._dp * self%nf/27 - 100._dp/3) * lgt**4 )/t

  end function ThrustSing3loop

!ccccccccccccccc

  function ThrustNS1loop(t) result(beta)
    real (dp), intent(in) :: t
    real (dp), dimension(3) :: beta

    real (dp) :: lg, t1, t2, o(4)

    if (t <= 0) then

      beta = 0._dp

    elseif (t < 1.d-3) then

      beta(1) = ThrustReg1loop(t)

      lg = log(t)
      t1 = 8/t/3 + t**5 * (- 965.942857142857036478744703345001 - 32 * lg)             &
      + 4 * (1 - 4 * lg)/3 - 16 * t**2 * (7 + lg) - 80 * t**4 * (55 + 3 * lg)/9

      t2 = - 16 * t * (25 + 6 * lg)/9 - 112 * t**6/45 * (761 + 15 * lg)                 &
      - 1.37142857142857144125969171000179 * t**8 * (5323 + 35 * lg) -                  &
      16 * t**3/45 * (679 + 60 * lg) - 1.69312169312169302770598733332008d-2 * t** 7    &
      * (219323 + 2520 * lg) - 7.69600769600769574196874600602314d-3 * t**9 *           &
      (1871029 + 6930 * lg)

       beta(2) = (t1 + t2)/2

       o(1) = t**2;  o(2) = lg

       t1 = - 16/t/3 - 8/o(1)/3 + t**8 * ( -129648.415584415577228583060787059  &
     - 480 * o(2) ) - o(1) * ( 745.6 + 640 * o(2) ) - 16 * t * ( 15 + 2 *o(2) )  &
     - 16 * ( 31 + 6 * o(2) )/9

       t2 = - 22.4 * t**5 * ( 509 + 10 * o(2) ) - 80 * t**3/9 * ( 223 + 12 * o(2) )   &
       - 2.28571428571428558740308289998211   * t**4 * (   2127 +   70 * o(2) ) -     &
        1.37142857142857144125969171000179    * t**7 * (  42619 +  280 * o(2) ) -     &
        9.31216931216931165238293033326045d-1 * t**9 * ( 306017 +  630 * o(2) ) -     &
        1.18518518518518511939419113332406d-1 * t**6 * ( 219683 + 2520 * o(2) )

       beta(3) = (t1 + t2)/2

    elseif (t < 1._dp/3) then

      beta(1) = ThrustFO1loop(t) - ThrustSing1loop(t)

      t2 = t**2

      beta(2) = (    (   4 * (  4 * (1 - 2 * t)**2 * log(1/t - 2)                      &
      + (t - 1) * ( t * (11 * t - 6 + 18 * t**3 - 27 * t2) - 4 * log(t) * (1 - 3 * t + &
      2 * t2) )  )/3   )/( (t - 1)**2 * (2 * t - 1) * t2 )    )/2

        o(1) = t - 1;  o(2) = t**3;  o(3) = - 3 * t;  o(4) = t2

        beta(3) = (    8 * (   - 4 * (1 - 2*t)**2 * log(1/t - 2) *                       &
      ( 1 + o(3) + 3 * o(4) ) + o(1) * (  4 * log(t) * (1 + o(3) + 2 * o(4) )**2 +       &
      t * (7 - 36 * t - 12 * o(2) + 45 * o(4) )  )   )/3    )/( (1 - 2 * t)**2 * o(1)**3 &
      * o(2) )/2

    else

      beta(1) = - ThrustSing1loop(t);  beta(2) = 2 * ( 1 - 4 * log(t) )/t**2/3
      beta(3) = (  8 * ( 4 * log(t) - 3 )/3  )/t**3/2

    endif

  end function ThrustNS1loop

!ccccccccccccccc

  function ThrustNS2loopFit(t) result(beta)
    real (dp), intent(in) :: t
    real (dp), dimension(2) :: beta

    real (dp) :: lt

    if (t > 0.095 .or. t < 0) then

      beta = 0

    else

      lt = log(t)

      beta(1) = ( - 29.11314324370869 - 13.428738641123346 * lt           &
      - 29.166275841929746 * lt**2 - 4.586460484332804 * lt**3            &
      + 104.2046742967832 * t * lt**3 )

      beta(2) = ( - 13.4287386411524489737701060221298                               - &
      58.3325516838740476543989643687382 * lt + (- 13.759381452999948081128422927577 + &
      312.614022890382337038772675441578 * t) * lt**2 +                                &
      104.204674296794119747744389314903 * t  * lt**3 )/t

    end if

  end function ThrustNS2loopFit

!ccccccccccccccc

  function ThrustNS2loopErr(t) result(beta)
    real (dp), intent(in) :: t
    real (dp), dimension(2) :: beta

    real (dp), dimension(6) :: o
    integer :: i

    if (t > 0.095 .or. t < 0) then

       beta = 0

    else

      o(1) = log(t)

      do i = 2, 6
        o(i) = o(1)**i
      end do

!       beta(1) = sqrt( 19874.91438325717 + 55791.62706807897 * o(1) &
!     + 52502.12485612066 * o(1)**2 + (19631.49966977773 - 23750.91712106788 * t) * o(1)**3  &
!     + (3507.233229654627 - 33352.90487307808_dp * t) * o(1)**4 + (303.2527662844569         &
!     - 7973.174677553763 * t) * o(1)**5 + (10.26944327499955 - 539.4757158995853 * t        &
!     + 7108.384721165775 * t**2) * o(1)**6 )/4

       beta(1) = sqrt(19874.9143832571717460666604893049 + 55791.6270680789594393900188151747 * o(1)&
       +52502.1248561206554938962653977796 * o(1)**2 + (19631.4996697777321976730036112713 - &
       23750.9171210678804442295586341061 * t) * o(1)**3 + (3507.2332296546271379611425800249 - &
       33352.9048730780752762825613899622 * t) * o(1)**4 + (303.252766284456898304711103264708 &
       -7973.17467755376263482958165695891 * t) * o(1)**5 + (10.2694432749995501197304292873014- &
       539.475715899585317458786448696628 * t + 7108.38472116577513304491731105372 * t**2) * o(1)**6)

      beta(2) = (   ( 27895.8135340394797196950094075873 + &
      52502.1248561206554938962653977796 * o(1) + (29447.2495046665949658404315414373   -  &
      35626.3756816018251072364364517853 * t) * o(2)+(7014.4664593092542759222851600498 -  &
      78581.2683066900863337878035963513 * t) * o(3)+(758.131915711142223557317265658639   &
      -36609.3891304234420047691855870653 * t) * o(4)+21325.1541634973262873131716332864 * &
      (-0.257219237434926517593680728168692+t)*(-5.61658750226658920468025826266967d-3+t) *&
      o(5) + t * (-269.737857949792658729393224348314 + 7108.38472116577513304491731105372*&
      t)*o(6) )/(  t * sqrt( 19874.9143832571717460666604893049 + &
      55791.6270680789594393900188151747 * o(1) + 52502.1248561206554938962653977796 * o(2)&
      +(19631.4996697777321976730036112713 - 23750.9171210678804442295586341061 * t) * o(3)&
      +(3507.2332296546271379611425800249 - 33352.9048730780752762825613899622 * t) * o(4) &
      +(303.252766284456898304711103264708-7973.17467755376263482958165695891 * t) *       &
      o(5) + (10.2694432749995501197304292873014 - 539.475715899585317458786448696628 * t  &
      +7108.38472116577513304491731105372 * t**2) * o(6) )  )   )/4

    endif
  end function ThrustNS2loopErr

!ccccccccccccccc

  function HJMNS2loopErr(t) result(beta)
    real (dp), intent(in)   :: t
    real (dp), dimension(2) :: beta

    real (dp) :: lt, o(6)
    integer :: i

    lt = log(t); o(1) = lt

    do i = 2, 6
      o(i) = o(1)**i
    end do

    beta(1) = sqrt( 42613.18684465871 + 120607.69825299302 * lt +                      &
    114355.56461109986 * lt**2 + (42813.71975342122 - 51894.17520264695 * t) * lt**3 + &
    (7641.643373149764 - 73857.76290961262 * t) * lt**4 + (659.5648803927595 -         &
    17628.676876506346 * t) * lt**5 + (22.289966484133263 - 1188.281878192541 * t +    &
    16035.082280316727 * t**2) * lt**6 )/4

    beta(2) = ( 60303.8491264965159643907099962234 + &
    114355.564611099858041143306763843 * o(1) + (64220.5796301318265761892689624801      - &
    77841.2628039704213023242118651979 * t) * o(2) + (15283.2867462995269924874719436048 - &
    173662.613420548694875833461992443 * t) * o(3) + (1648.9122009818988345841717091389  - &
    81000.5736460721870173529168823734 * t) * o(4) + 48105.246840950179887386184418574   * &
    (-2.51815205035326972193843175773509d-1 + t) * (-5.52021855702979813429465139051899d-3 &
    + t) * o(5) + t * (-594.140939096270503938512774766423 + &
    16035.082280316725888980045056087 * t)* o(6) )/(  t * &
    sqrt( 42613.1868446587169785289006540552 + 120607.698252993023046997222991195 * o(1) + &
    114355.564611099858041143306763843 * o(2) + (42813.7197534212177174595126416534 -      &
    51894.1752026469504954775402438827 * t) * o(3) + (7641.64337314976371828834089683369   &
    -73857.762909612620916277592186816 * t) * o(4) + (659.56488039275953383366868365556    &
    -17628.67687650634618279354981496 * t) * o(5) + (22.2899664841332656450845206563827 -  &
    1188.28187819254105228594653453911 * t + 16035.082280316725888980045056087 * t**2) * o(6) )  )

  end function HJMNS2loopErr

!ccccccccccccccc

  real (dp) function ThrustNS3loopErr(t)
    real (dp), intent(in) :: t

    real (dp) :: lt, lt2

    if (t > 0.315) then

       ThrustNS3loopErr = 0._dp

    else

      lt  = log(t);  lt2 = lt**2

      ThrustNS3loopErr = (  66.3570537969770946062908478779718 * sqrt( lt2 * &
      (1.53609472457735374284482077200664 + 2.46331080820060854819075757404789 * lt &
      + lt2) * (3.8376549851267780510966076690238 + 3.81531269385931448567816914874129 * lt &
      + lt2) * (10.7197084488020610848479918786325 + 6.35650399375444674632262831437401 * lt &
      + lt2) * (47.4943333612261220366690395167097 + 13.727393246491039757728458425845 * lt &
      + lt2) )  )/8

    endif

  end function ThrustNS3loopErr

!ccccccccccccccc

  real (dp) function HJMNS3loopErr(t)
    real (dp), intent(in) :: t

    real (dp) :: lt

    lt = log(t)

    HJMNS3loopErr = (  11.804875053380671 * Sqrt( lt**6 * (4.094141071178425 +       &
    3.8643010319557316 * lt + lt**2) * (32.93500856923923 + 11.44458318375263 * lt &
    + lt**2) )  )/8

  end function HJMNS3loopErr

!ccccccccccccccc

  real (dp) function ThrustNS3loopFit(self, t)
    real (dp), intent(in) :: t
    class(NonSingular), intent(in) :: self

    real (dp) :: lt

!        o(1)=5.d-1 * 8998.054502925055_dp
!        self%s3=A3+j3+o(1)

    lt = log(t)

    ThrustNS3loopFit = (  (-22638.646534945290333951106731547 - 1349.7568075515132  &
    + 4.27111946801019470854043902363628d-1 * self%s3) * lt +                       &
    (-41034.460807046251318297436228022 + 905.899405165561 +                        &
    8.69059070422017576618145540123805d-1 * self%s3) * lt**2 +                      &
    (-24394.3073042593994159688008949161 + 361.01421334437197 +                     &
    6.20310699650405883431858455878682d-1 * self%s3) * lt**3 +                      &
    (-5867.95724199124091313706230721436 + 16.739842058772542 +                     &
    1.83196704512979069434663870197255d-1 * self%s3) * lt**4 +                      &
    (-457.708008404026944759834805154242 - 12.07645927763132 +                       &
    1.88791814905350663345018347172299d-2 * self%s3) * lt**5 )/8

  end function ThrustNS3loopFit

!ccccccccccccccc

  real (dp) function HJMNS3loopFit(self, t)
    real (dp), intent(in) :: t
    class(NonSingular), intent(in) :: self

    real (dp) :: lt

    lt = log(t)

    HJMNS3loopFit = ( 4332.502710938817 * lt**3 - 3.968041292810291 * self%s2 * lt**3  &
    - 11.92014658128103 * self%s2rho * lt**3 + 1.014458922948192 * self%s3 * lt**3 + &
    2018.001331545123 * lt**4 - 1.032416871289466 * self%s2 * lt**4 -                &
    5.140627417476885 * self%s2rho * lt**4 + 0.9613964981686627 * self%s3 * lt**4  &
    + 199.64534524153854 * lt**5 - 0.04871769071767448 * self%s2 * lt**5 -           &
    0.31497005317075477 * self%s2rho * lt**5 + 0.20880811120915 * self%s3 * lt**5 )/8

  end function HJMNS3loopFit

!ccccccccccccccc

  function HJMNS2loopFit(self, t) result(beta)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp), dimension(2) :: beta

    real (dp) :: lt
    lt = log(t)

    beta(1) = ( - 6847.412797589124 + 403.28576694938266 * self%s2rho -        &
    9270.035609868502 * lt + 546.371559245844 * self%s2rho * lt -              &
    2349.2073741763643 * lt**2 + 135.64141753046036 * self%s2rho * lt**2 -     &
    175.59127856914506 * lt**3 +  9.914222038931396 * self%s2rho * lt**3 +     &
    3705.7183041666844 * t * lt**3 - 210.34367126583035 * self%s2rho * t * lt**3 )/4

    beta(2) = ( - 9270.03560986850239089562819572166 + 546.371559245844018448678980348632 &
    * self%s2rho + (- 4698.41474835272876475755765568465 + 271.282835060920701764075602113735 &
    * self%s2rho) * lt + (- 526.773835707435189590341906296089 + self%s2rho * &
    (29.7426661167941874808207103342284 - 631.031013797491091565916576655582 * t) + &
    11117.1549125000534985474587301724 * t) * lt**2 + (3705.7183041666847955752928100992 &
     - 210.343671265830334249358202214353 * self%s2rho) * t * lt**3 )/t

  end function HJMNS2loopFit

!ccccccccccccccc

  real (dp) function CumThrustSing1loop(t)
    real (dp), intent(in) :: t
    real (dp) :: lt

    if (t > 0) then

      lt = log(t);  CumThrustSing1loop = 1.526578755797635 - 2 * lt - 4 * lt**2/3

    else

      CumThrustSing1loop = 0

    end if

  end function CumThrustSing1loop

!ccccccccccccccc

  real (dp) function CumThrustSing2loop(self, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp) :: lg

    if (t > 0.) then

       lg = log(t)

       CumThrustSing2loop = 23.1052507826620967534836381673813 - &
     5.8647508177380291982672133599408d-1 * self%nf + (-12.673650220838139635048946729512  &
     + 5 * self%nf/6) * lg + (-7.28866911171631048915742212557234 +                        &
     11 * self%nf/54._dp) * lg**2 + (19/3._dp - 2 * self%nf/9._dp) * lg**3 + 8 * lg**4/9

    else

       CumThrustSing2loop = 0

    end if

  end function CumThrustSing2loop

!ccccccccccccccc

  real (dp) function CumThrustSing3loop(self, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp) :: lg
    integer :: nf2

    if (t > 0) then

      nf2 = self%nf**2;  lg = log(t)

      CumThrustSing3loop = 447.055436980793885481944016646594 +1.5625d-2 * self%s3 &
      -58.3515328305320846169479409581982 * self%nf + 1.43961598822072578407471610262292 * &
      nf2+(-143.525152363904551222617556049954 + 12.0822164629231987298396688856883 * self%nf &
      -6.30586823940581808756178361363709d-2 * nf2) * lg + (-61.306532939293312267636792967096&
      +1.21345418994326403883121656690491 * self%nf + 1.57855198367358351418943129829131d-1 * &
      nf2) * lg**2 + (102.630499432346589117059920681641 - 8.64842572371231810279823548626155 &
      * self%nf + 1.04938271604938271330809129722184d-1 * nf2) * lg**3 + &
      (-10.7359556339761508247931942605646 + 1.59876543209876542661618259444367 * self%nf - &
      4.32098765432098730343568604439497d-2 * nf2) * lg**4 + (-20/3._dp + &
      2.96296296296296279848547783331014d-1 * self%nf) * lg**5 - &
      3.95061728395061706464730377774686d-1 * lg**6

    else

       CumThrustSing3loop = 0

    end if

  end function CumThrustSing3loop

!ccccccccccccccc

  double precision function CumThrustNS1loop(t)
    real (dp), intent(in) :: t
    real (dp)             :: dilog, lg

    if (t < 0) then

      CumThrustNS1loop = 0

    else if (t < 1.e-3_dp) then

      lg = log(t)

      CumThrustNS1loop = 2 * t * (lg - 2)/3 - t**4 * (77 + 12 * lg)/9 + t**8 * &
      (-33.6436507936507922522650915198028 - 2 * lg/3) - 5.92592592592592559697095566662028d-1 &
      * t**3 * (10 + 3 * lg) - 2 * t**2 * (4 * lg - 7)/3 - 16 * t**5 * (163 + 15 * lg)/225 &
      - 2.96296296296296279848547783331014d-2 * t**6 * (539 + 30 * lg) - &
      7.2562358276643994514643054571934d-3 * t**7 * (3137 + 105 * lg) - &
      1.88124632569077032684390360373072d-3 * t**9 * (27341 + 315 * lg)

        CumThrustNS1loop = CumThrustNS1loop/2._dp

    else if (t < 1._dp/3) then

      lg = log(t)

      CumThrustNS1loop = 2 * (  19.739208802178718205055929502123 + 9 * t * (4 + 3 * t) &
      -6._dp*log(1 - 2 * t) * ( -2.27411277760218766275102098006755d-1 + 6 * t + &
      4 * log(1 - t) ) + 36 * t * lg - 24 * dilog(1 - t) - 24 * dilog(2 * t)  &
      - 24 * dilog(2 * t - 1)  )/18

    else

        CumThrustNS1loop = 1 - CumThrustSing1loop(t)

    end if

  end function CumThrustNS1loop

!ccccccccccccccc

  real (dp) function CumThrustNS2loopFit(t)
    real (dp), intent(in) :: t

    real (dp) :: lg

    if (t < 0) then

      CumThrustNS2loopFit = 0

    else if (t .ge. 0.095_dp) then

      CumThrustNS2loopFit = - 22.78579241472016

    else

      lg = log(t)

      CumThrustNS2loopFit = -46.4981933804481251115703344112262 * t &
      -39.0767528612957448785891756415367 * t**2 + t * &
      (17.3850501367303422739496454596519 + 78.1535057225914897571783512830734 * t) * lg &
      + (-15.4068943889341203323795070900815 - 78.1535057225914897571783512830734 * t)   &
      * t * lg**2 + t * (-4.58646048433308806124841794371605 + &
      52.102337148394326504785567522049 * t) * lg**3

    end if

      CumThrustNS2loopFit = CumThrustNS2loopFit/4

  end function CumThrustNS2loopFit

!ccccccccccccccc

  real (dp) function CumThrustNS3loopFit(self, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp) :: lg

    if (t < 0) then

      CumThrustNS3loopFit = 0

    else if (t .ge. 0.315) then

      CumThrustNS3loopFit = - 1345.85610273526978808433796075406 - self%s3/4

    else

      lg = log(t)

      CumThrustNS3loopFit = t * (     - 81.6625163046992419424441322917119 - &
      2.79638874412130089552874778746627d-1 * self%s3 + lg * &
      (    81.6625163046992419424441322917119 + 2.79638874412130089552874778746627d-1 * self%s3 &
      + lg * (   -9724.47923442373252100878744386137 + 7.37365361944446462416635768022388d-2 * &
      self%s3 + lg * (  -8831.38352079528132776431448291987 + &
      2.65107511409190976792160654440522d-1 * self%s3 + lg * &
      ( -3102.77671484495165543648909078911 + 8.88007970603036866918955638539046d-2*self%s3 &
      + (-384.846273144581285663434755406342 + 1.88791814905350663345018347172299d-2 * &
      self%s3) * lg )  )   )    )     )

    end if

      CumThrustNS3loopFit = CumThrustNS3loopFit/8

  end function CumThrustNS3loopFit

!ccccccccccccccc

  double precision function CumHJMSing2loop(self, t)
    class(NonSingular), intent(in) :: self
    real (dp), intent(in) :: t
    real (dp) :: L

    if (t > 0) then

      L = log(t)

      CumHJMSing2loop = 8 * L**4/9 + L**3 * (19._dp/3 - 2 * self%nf/9._dp) + &
       L**2 * ( -65._dp/36 + 11 * self%nf/54._dp - 7 * Pi2/27)             + &
      L * ( - 53._dp/4 + 5 * self%nf/6._dp - 4 * Pi2/9 + 38 * Zeta3/9)      + &
      22.025732322734132 - 0.5864750817738028 * self%nf + 4 * self%s2rho

    else

      CumHJMSing2loop = 0._dp

    end if

  end

!ccccccccccccccc

end module NonSingularClass
