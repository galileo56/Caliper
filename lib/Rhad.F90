
module SigmaClass
  use AnomDimClass;  use RunningClass;  use ElectroWeakClass; use Chaplin
  use Constants, only: dp, Pi, Pi2, Zeta3;  implicit none;  private

  real (dp), dimension(2)               :: EWfact
  public                                :: setEWfact, Pi0, Pi1, Pi0Der, Pi1Der, &
  P2, P2Der

!ccccccccccccccc

  type, public :: Sigma
   private
   integer                        :: nf
   real (dp), dimension(0:3,4)    :: RhadCoef
   real (dp)                      :: m
   type (AnomDim)                 :: andim
   type (Running)                 :: run
   type (ElectroWeak)             :: EW
   character (len = 6)            :: scheme

   contains

    procedure, pass(self), public :: RhadMass, SetAlpha, SetMTop, SetMBottom, &
    Rhad, SetMCharm, RHadCoefs

  end type Sigma

!ccccccccccccccc

  interface Sigma
    module procedure InitSigma
  end interface Sigma

  contains

!ccccccccccccccc

   type (Sigma) function InitSigma(run, EW)
    type (Running), intent(in) :: run
    type (ElectroWeak)         :: EW
    integer                    :: nf

    InitSigma%run = run           ;   InitSigma%andim = run%adim()
    nf = InitSigma%andim%numFlav();   InitSigma%nf    = nf
    InitSigma%m = run%scales('mL');   if (nf < 5) InitSigma%m = 0
    InitSigma%EW = EW             ;   InitSigma%scheme = run%scheme()

    InitSigma%RhadCoef = 0;  InitSigma%RhadCoef(0,1) = 1
    InitSigma%RhadCoef(0,2) = 1.985707398577798_dp   - 0.11529539789360388_dp * nf
    InitSigma%RhadCoef(0,3) = - 6.636935585488629_dp - 1.2001340534564595_dp * nf &
    - 0.0051783553199245685_dp * nf**2

    InitSigma%RhadCoef(0,4) = - 156.61_dp + 18.77_dp * nf - 0.7974_dp * nf**2 &
    + 0.0215_dp * nf**3

    call InitSigma%andim%expandAlpha(InitSigma%RhadCoef)

   end function InitSigma

!ccccccccccccccc

  subroutine SetAlpha(self, alpha)
    class (Sigma), intent(inout) :: self
    real (dp)    , intent(in   ) :: alpha

    call self%run%setAlpha(alpha)

  end subroutine SetAlpha

!ccccccccccccccc

  subroutine SetMtop(self, mT, muT)
    class (Sigma), intent(inout) :: self
    real (dp)    , intent(in   ) :: mT, muT

    call self%run%SetMtop(mT, muT)

  end subroutine SetMtop

!ccccccccccccccc

  subroutine SetMBottom(self, mB, muB)
    class (Sigma), intent(inout) :: self
    real (dp)    , intent(in   ) :: mB, muB

    call self%run%SetMBottom(mB, muB)

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, mC, muC)
    class (Sigma), intent(inout) :: self
    real (dp)    , intent(in   ) :: mC, muC

    call self%run%SetMCharm(mC, muC)

  end subroutine SetMCharm

!ccccccccccccccc

  function RHadCoefs(self) result(coefs)
    class (Sigma), intent(in) :: self
    real (dp), dimension(4)   :: coefs

    coefs = self%RhadCoef(0,:)

  end function RHadCoefs

!ccccccccccccccc

  pure real (dp) function Rhad(self, order, mu, Q)
    class (Sigma)               , intent(in) :: self
    Integer                     , intent(in) :: order
    real (dp)                   , intent(in) :: mu, Q
    real (dp), dimension(0:min(order,4))     :: alphaList, rQ
    real (dp), dimension(0:min(order,4) - 1) :: logList
    integer                                  :: i, n

    rQ = 1;  logList(0) = 1; alphaList(0) = 1 ; n = min(order,4) - 1
    if (order > 1)   logList(1:) = PowList( log(mu/Q), n)
    if (order > 0) alphaList(1:) = PowList( self%run%alphaQCD(mu)/Pi, min(order,4) )

    alpha_loop: do i = 2, min(order,4)
      rQ(i) = dot_product( self%RhadCoef(:n,i), logList )
    end do alpha_loop

    Rhad = dot_product( alphaList, rQ )

  end function Rhad

!ccccccccccccccc

  real (dp) function RhadMass(self, current, order, mu, Q)
    class (Sigma)      , intent(in) :: self
    Integer            , intent(in) :: order
    real (dp)          , intent(in) :: mu, Q
    character (len = *), intent(in) :: current
    real (dp)                       :: v, mb, m, alphaRH
    real (dp), dimension(4)         :: delta
    real (dp), dimension(2)         :: EWfactors

    if ( self%scheme(:4) == 'pole'  ) mb = self%m
    if ( self%scheme(:5) == 'MSbar' ) mb = self%run%MSbarMass(mu)

    m = mb/Q;  RhadMass = 0

    if ( 2 * m < 1) then

      EWfactors = self%EW%EWfactors(self%nf, Q)
      call setEWfact( EWfactors/sum(EWfactors) )
      v = sqrt( 1 - m**2 );  alphaRH = self%run%alphaQCD(mu)/Pi

      if (order >= 0) RhadMass = Rtree(current, v)
      if (order >= 1) then

        RhadMass = RhadMass + alphaRH * R1loop(current, v)

        if ( self%scheme(:5) == 'MSbar' ) then

          delta = self%run%MSbarDeltaMu(mu)
          RhadMass = RhadMass + A0MS(current, m) * delta(1)/mb

        end if
      end if
    end if

  end function RhadMass

!ccccccccccccccc

  pure recursive function Rtree(current, v) result(tree)
    real (dp)          , intent(in) :: v
    character (len = *), intent(in) :: current

    real (dp) :: tree

    tree = 0

    if (current(:6) == 'vector') then
      tree = v * (3 - v**2)/2
    else if (current(:5) == 'axial') then
      tree = v**3
    else if (current(:3) == 'all') then
      tree = EWfact(1) * Rtree('vector', v) + EWfact(2) * Rtree('axial', v)
    end if

  end function Rtree

!ccccccccccccccc

  recursive function R1loop(current, v) result(tree)
    real (dp)          , intent(in) :: v
    character (len = *), intent(in) :: current

    real (dp) :: tree, x, lx

    tree = 0

    if ( current(:3) == 'all' ) then
      tree = EWfact(1) * R1loop('vector', v) + EWfact(2) * R1loop('axial', v)
    else

      x = 1 - v

      if (x < 1e-3_dp) then

        lx = log(x/2)

        if ( current(:6) == 'vector' ) then
          tree =  1 + 6 * x - x**2 * (0.5 + 6 * lx) - x**3 * (161._dp/54 + 4 * lx/9)
        else if ( current(:5) == 'axial' ) then
          tree = 1 - x * (3 + 6 * lx) - x**2 * (1 - 9 * lx) + x**3 * (353._dp/108  - 40 * lx/9)
        end if
      else

        if ( current(:6) == 'vector' .or. current(:5) == 'axial' ) tree = 4 * Rtree(current,v) * KV(current,v)/3

      end if
    end if

  end function R1loop

!ccccccccccccccc

  pure real (dp) function PV(current, v)
    real (dp)          , intent(in) :: v
    character (len = *), intent(in) :: current

    PV = 0

    if (current(:6) == 'vector') then
      PV = 33._dp/24 + 11 * v**2/12 - 7 * v**4/24
    else if (current(:5) == 'axial') then
      PV = 21._dp/32 + 59 * v**2/32 + 19 * v**4/32 - 3 * v**6/32
    end if

  end function PV

!ccccccccccccccc

  pure real (dp) function QV(current, v)
    real (dp)          , intent(in) :: v
    character (len = *), intent(in) :: current

    QV = 0

    if (current(:6) == 'vector') then
      QV = 5 * v/4 - 3 * v**3/4
    else if (current(:5) == 'axial') then
      QV = - 21 * v/16 + 15 * v**3/8 + 3 * v**5/16
    end if

  end function QV

!ccccccccccccccc

  real (dp)  function KV(current, v)
    real (dp)          , intent(in) :: v
    character (len = *), intent(in) :: current
    real (dp)                       :: b, c, d

    c = 1 + v;  b = (1 - v)/c
    if (current(:6) == 'vector') d = 1/(1 - v**2/3)
    if (current(:5) == 'axial')  d = 1/v**2

    KV = ( Amass(v) - PV(current,v) * d * Log(b) + d * QV(current,v) )/v

  end function KV

!ccccccccccccccc

  real (dp) function Amass(v)
    real (dp), intent(in) :: v

    real (dp) :: b, c, Dilog

    c = 1 + v;  b = (1 - v)/c

    Amass = ( DiLog(b**2) + 2 * DiLog(b) - log(b) * log(c**3/8/v**2) ) * &
    (1 + v**2) + 3 * v * log( (1 - v**2)/4/v ) - v * log(v)

  end function Amass

!ccccccccccccccc

  pure recursive function A0MS(current, m) result(A0)
    real (dp)          , intent(in) :: m
    character (len = *), intent(in) :: current
    real (dp)                       :: v, A0

    A0 = 0

    if ( current(:3) == 'all' ) then
      A0 = EWfact(1) * A0MS('vector', m) + EWfact(2) * A0MS('axial', m)
    else

      if (m > 1e-3_dp) then

        v = sqrt(1 - 4 * m**2)
        if (current(:6) == 'vector') A0 = - 24 * m**4 / v
        if (current(:5) == 'axial' ) A0 = - 12 * m**2 * v

      else
        if (current(:6) == 'vector') A0 = -  24 * m**3 -    4 * m**5  - 144 * m**7 &
                                          - 480 * m**9 - 1680 * m**11
        if (current(:5) == 'axial' ) A0 = - 12 * m + 24 * m**3 + 24 * m**5 + 48 * m**7 &
                                          + 120 * m**9 + 336 * m**11
      end if
    end if

  end function A0MS

!ccccccccccccccc

  subroutine setEWfact(EW)
    real (dp), dimension(2) :: EW
    EWfact = EW
  end subroutine setEWfact

!ccccccccccccccc

  complex (dp) function u(z)
    complex (dp), intent(in) :: z

    u = 1 + 2 * Sqrt(z - 1) *  Sqrt(z) - 2 * z

  end function u

!ccccccccccccccc

  complex (dp) function uDer(i, u)
    integer     , intent(in) :: i
    complex (dp), intent(in) :: u

    uDer = 0

    if ( i == 0 ) then
      uder = u
    else if (i == 1) then
      uDer = 4 * u**2/(1 - u**2)
    else if (i == 2) then
      uDer = 32 * u**3/(1 - u**2)**3
    else if (i == 3) then
      uDer = 384 * u**4*(1 + u**2)/(1 - u**2)**5
    end if

  end function uDer

!ccccccccccccccc

  complex (dp) function Gz(u)
    complex (dp), intent(in) :: u

    Gz = ( 2 * u * Log(u) )/(u**2 - 1)

  end function Gz

!ccccccccccccccc

  complex (dp) function GzDer(i, u)
    integer     , intent(in) :: i
    complex (dp), intent(in) :: u

    GzDer = 0

    if ( i == 0 ) then
      GzDer = Gz(u)
    else if (i == 1) then
      GzDer = - 2 * ( 1 - u**2 + (1 + u**2) * log(u) )/(1 - u**2)**2
    else if (i == 2) then
      GzDer = ( 2 + 4 * u**2 - 6 * u**4 + 4 * u**2 * (3 + u**2) * log(u) )&
      /u/(u**2 - 1)**3
    else if (i == 3) then
      GzDer = - 2 * ( 15 * u**2 - 1 - 3 * u**4 - 11 * u**6 + 6 * (u**2 + &
      6 * u**4 + u**6) * log(u) )/u**2/(u**2 - 1)**4
    end if

  end function GzDer

!ccccccccccccccc

  complex (dp) function Iz(u)
    complex (dp), intent(in) :: u

    Iz = 6 * ( Zeta3 + 4 * cli3(-u) + 2 * cli3(u) ) - 8 * ( 2 * cli2(-u) + &
    cli2(u) ) * Log(u) - 2 * ( 2 * Log(1 + u) + Log(1 - u) ) * Log(u)**2

  end function Iz

!ccccccccccccccc

  complex (dp) function IzDer(i, u)
    integer     , intent(in) :: i
    complex (dp), intent(in) :: u

    if (i == 0) then

      IzDer = Iz(u)

    else if (i == 1) then

      IzDer = 2 * (  Log(u) * ( 2 * (u**2 - 1) * Log(1 - u) + (1 - 3 * u) * u &
      * Log(u) + 4 * (u**2 - 1) * Log(1 + u) ) + 4 * (u**2 - 1) * cli2(-u) +  &
      2 * (u**2 - 1) * cli2(u)  )/u/(u**2 - 1)

    else if (i == 2) then

      IzDer = (  2 * Log(u) * ( - 2 * (u**2 - 1)**2 * Log(1 - u) + u**2 * &
      ( 3 + u * (3 * u - 2) ) * Log(u) - 4 * (u**2 - 1)**2 * Log(1 + u) ) + &
      4 * (u**2 - 1)**2 * cli2(u) - 4 * (u**2 - 1)**2 * cli2(u**2)  )&
      /u**2/(u**2 - 1)**2

    else if (i == 3) then

      IzDer = 4 * (  Log(u) * ( 2 * (u**2 - 1)**3 * Log(1 - u) + u * (1 - &
      u * (6 + (u - 6) * u**2) + u**2 * ( 1 - 3 * u * (3 + (u - 1) * u) ) * &
      Log(u)) + 4 * (u**2 - 1)**3 * Log(1 + u) ) - 2 * (u**2 - 1)**3 * cli2(u) &
      + 2 * (u**2 - 1)**3 * cli2(u**2)  )/u**3/(u**2 - 1)**3

    end if

  end function IzDer

!ccccccccccccccc

complex (dp) function DIz(z, u)
  complex (dp), intent(in) :: z, u
  complex (dp) :: lu, u21, sz, c2, c21, lu1, lu2

  lu = log(u); u21 = u**2 - 1; sz = Sqrt(z - 1); c2 = cli2(u); c21 = cli2(-u)
  lu1 = Log(1 - u); lu2 = Log(1 + u)

  DIz = 2 * (   u21 * sz * lu**2 * ( lu1 + 2 * lu2 ) + 4 * u21 * sz * lu * &
  ( 2 * c21 + c2 ) - (  Sqrt(z) * (lu * ( 4 * u21 * lu1 + (3 - 4 * u - 3 * u**2 &
  + 12 * u * z) * lu + 8 * u21 * lu2 ) + 8 * u21 * c21 + 4 * u21 * c2)  )/2 - &
  3 * u21 * sz * (4 * cli3(-u) + 2 * cli3(u) + Zeta3)   )/u21/sz/z

end function DIz

!ccccccccccccccc

complex (dp) function Pi0(z)
  complex (dp), intent(in) :: z

  Pi0 = 3 * (  20._dp/9 + 4/z/3 - 4 * (1 - z) * (1 + 2 * z)/3/z * Gz( u(z) )  )/16/Pi2

end function Pi0

!ccccccccccccccc

complex (dp) function Pi0Der(i, z)
  integer     , intent(in) :: i
  complex (dp), intent(in) :: z
  complex (dp) :: uu, uuDer, uuDer2, uuDer3

  Pi0Der = 0; if (i > 0) uu = u(z); if (i > 0) uuDer = uDer(1,uu)
  if (i > 1) uuDer2 = uDer(2,uu);  if (i > 2) uuDer3 = uDer(3,uu)

  if ( i == 0 ) then

    Pi0Der = Pi0(z)

  else if (i == 1) then

    Pi0Der = 4 * ( (1 + 2 * z**2) * Gz(uu) - 1 + z * (2 * z**2 - 1 - z) * &
    Gzder(1,uu) * uuDer )/3/z**2

  else if (i == 2) then

    Pi0Der = 4 * (  2 - 2 * Gz(uu) + (z - 1) * z**2 * (1 + 2 * z) * uuDer**2 * &
    Gzder(2,uu) + z * Gzder(1,uu) * ( 2 * (1 + 2 * z**2) * uuDer + (z - 1) * &
    z * (1 + 2 * z) * uuDer2 )  )/3/z**3

  else if (i == 3) then

    Pi0Der = 4 * (   6 * Gz(uu) - 6 + z * (   z * ( uuDer3 * (z - 1) * z * &
    (1 + 2 * z) + uuDer2 * (3 + 6 * z**2) ) - 6 * uuDer  ) * GZDer(1,uu) + &
    uuDer * z**2 * (  3 * ( uuDer2 * (z - 1) * z * (1 + 2 * z) + uuDer * &
    (1 + 2 * z**2) ) * GZDer(2, uu) + uuDer**2 * (z - 1) * z * (1 + 2 * z) *&
     GZDer(3, uu)  )   )/3/z**4

  end if

  if (i > 0) Pi0Der = 3 * Pi0Der/16/Pi2

end function Pi0Der

!ccccccccccccccc

complex (dp) function Pi1Der(i, z)
  integer     , intent(in) :: i
  complex (dp), intent(in) :: z
  complex (dp) :: uu, uuDer, uuDer2, uuDer3

  Pi1Der = 0; if (i > 0) uu = u(z); if (i > 0) uuDer = uDer(1,uu)
  if (i > 0) uuDer2 = uDer(2,uu);  if (i > 1) uuDer3 = uDer(3,uu)

  if ( i == 0 ) then

    Pi1Der = Pi1(z)

  else if (i == 1) then

    Pi1Der = (z*(-1 + 16*z**2)*Gz(uu)**2 + 2*z*Gz(uu)*(9 + 6*z**2 + &
    uuDer*(-1 + z)*z*(-1 + 16*z)*GZDer(1,uu)) - 2*Iz(U(z)) +  &
    z*(-13 + 2*uuDer2*(-1 + z)*z*(1 + 2*z)*IzDer(1,uu) + uuDer*(3*IzDer(1,uu) &
    + 2*(-1 + z)*z*((9 + 6*z)*GZDer(1,uu) + uuDer*(1 + 2*z)*IzDer(2,uu)))))/ &
    (6*z**3)

  else if (i == 2) then

    Pi1Der = (2*z*Gz(uu)**2 + 2*z*Gz(uu)*(-18 + z*(uuDer2*(-1 + z)*z*(-1 + 16*z)  &
    + uuDer*(-2 + 32*z**2))*GZDer(1,uu) + uuDer**2*(-1 + z)*z**2*(-1 + 16*z)* &
    GZDer(2,uu)) + 6*Iz(U(z)) + z*(26 + 6*z*(uuDer2*z*(-3 + z + 2*z**2) + &
    uuDer*(6 + 4*z**2))*GZDer(1,uu) + 2*uuDer**2*(-1 + z)*z**2*(-1 + 16*z)* &
    GZDer(1,uu)**2 + (-8*uuDer + z*(2*uuDer3*(-1 + z)*z*(1 + 2*z) + &
    uuDer2*(5 + 4*z**2)))*IzDer(1,uu) + uuDer*z*(6*uuDer2*(-1 + z)*z*(1 + 2*z) &
    *IzDer(2,uu) + uuDer*(6*(-1 + z)*z*(3 + 2*z)*GZDer(2,uu) + &
    (5 + 4*z**2)*IzDer(2,uu)) + 2*uuDer**2*(-1 + z)*z*(1 + 2*z)*IzDer(3,uu))))/(6*z**4)

  end if

  if (i > 0) Pi1Der = 3 * Pi1Der/16/Pi2

end function Pi1Der

!ccccccccccccccc

complex (dp) function Pi1(z)
  complex (dp), intent(in) :: z
  complex (dp)             :: uu, Gzuu, IIz

  uu = u(z); Gzuu = Gz(uu); IIz = Iz(uu)

  Pi1 = 3 * (  5._dp/6 + 13/z/6 - (1 - z) * (3 + 2 * z) * Gzuu/z + &
  (1 - z) * (1 - 16 * z) * Gzuu**2/z/6 - (1 + 2 * z)/z/6 * ( IIz/z + &
  2 * (1 - z) * ( IzDer(1,uu) * uDer(1,uu) - IIz/z )  )   )/16/Pi2

end function Pi1

!ccccccccccccccc

complex (dp) function P2Der(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: r32
  complex (dp)                :: r0, rz, r1, r2, r3, l1, l2, l3, PL1, PL2, PL3, &
  PL4, l4, d1, d2

  rz = sqrt(z); r1 = sqrt(1 - z); r2 = 1 + r1; r3 = 1 - r1; r0 = sqrt(z - 1)
  d2 = 1 + 2 * r0 * rz - 2 * z; d1 = d2**2 - 1; l1 = Log(d2); PL2 = cli2(-d2)
  l2 = Log(1 + r0 * rz - z); l3 = Log( 8 * (z - r0 * rz) ); PL1 = cli2(d2)
  PL3 = cli3(d2); PL4 = cli3(-d2); l4 = Log(1 + d2)

  r32 = PowList(r3/r2, 12)

  P2Der = 2/(r2 - z)**2 * (1 + 1/r1/2) * (0.688472275231591_dp + &
  0.1891443911186843_dp * r3/r2 - 0.030243753917029492_dp * r3**2/r2**2 + &
  0.0034204077615457664_dp * r3**3/r2**3 + 0.006118445792320977_dp * r3**4/r2**4 &
  + 0.0016455445946186606_dp * r3**5/r2**5 + 0.001252836593493741_dp * r3**6/r2**6 &
  + 0.00046831306017028774_dp * r3**7/r2**7 + 0.0018087001167161186_dp * r3**8/r2**8 &
  + 0.00048205327014239514_dp * r3**9/r2**9 - 0.0024098794332697245_dp * r3**10/r2**10 &
  - 0.0001825139597445611_dp * r3**11/r2**11 + 0.0013086431124221074_dp * r3**12/r2**12)&
  + 2/(r2 - z)/r1/r2 * ( 0.09457219555934215_dp + 0.06432844164231266_dp * r32(1) - &
  0.025113142274710842_dp * r32(2) + 0.017367503226960603_dp * r32(3) + &
  0.016350753071188606_dp * r32(4) + 0.007872371267027875_dp * r32(5) + &
  0.00539760549107723_dp * r32(6) + 0.008873896177460481_dp * r32(7) + &
  0.009404040182505253_dp * r32(8) - 0.009880157450707843_dp * r32(9) - &
  0.013053223944943708_dp * r32(10) + 0.006848031895937559_dp * r32(11) + &
  0.007851858674532645_dp * r32(12) )

end function P2Der

!ccccccccccccccc

complex (dp) function P2(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: r32
  complex (dp)                :: rz, r1, r2, r3, l1, l2, l3, PL1, PL2, PL3, PL4, &
  d1

  rz = sqrt(z); r1 = sqrt(1 - z); r2 = 1 + r1; r3 = 1 - r1; r1 = sqrt(z - 1)
  d1 = 1 + 2 * r1 * rz - 2 * z; l1 = Log(d1); l2 = Log(1 + r1 * rz - z)
  l3 = Log( 8 * (z - r1 * rz) ); PL1 = cli2(d1); PL2 = cli2(-d1); PL3 = cli3(d1)
  PL4 = cli3(-d1); r32 = PowList(r3/r2, 12)

  P2 = - 0.688472275231591_dp + 2 * ( 0.688472275231591_dp + &
  0.1891443911186843_dp * r32(1) - 0.030243753917029492_dp * r32(2) + &
  0.0034204077615457664_dp * r32(3) + 0.006118445792320977_dp * r32(4) + &
  0.0016455445946186606_dp * r32(5) + 0.001252836593493741_dp * r32(6) + &
  0.00046831306017028774_dp * r32(7) + 0.0018087001167161186_dp * r32(8) + &
  0.00048205327014239514_dp * r32(9) - 0.0024098794332697245_dp * r32(10) &
  - 0.0001825139597445611_dp * r32(11) + 0.0013086431124221074_dp * r32(12) )/(r2 - z) + &
  (0.0035539160594078058_dp + 0.012811625542535428_dp * z - &
  0.011957520634128492_dp * z**2 + 1.5887591794330882_dp * z**3 + &
  0.555408076785252_dp * z**4 + 0.17805743238532604_dp * z**5)/z**5 + &
  (1.74901919894622e-2_dp*PL1**2*r1**2*rz**2*(0.5_dp+z)**2)/z**5+  &
  (6.99607679578488e-2_dp*PL2**2*r1**2*rz**2*(0.5_dp+z)**2)/z**5+  &
  (PL3**2*(9.83823299407249e-3_dp-7.870586395257993e-2_dp*z**2+  &
  1.5741172790515985e-1_dp*z**4))/z**5 +  &
  (PL4**2*(3.9352931976289964e-2_dp-3.148234558103197e-1_dp*z**2+  &
  6.296469116206394e-1_dp*z**4))/z**5+(PL4*(2.3652231770834637e-2_dp  &
  +4.263234297431412e-2_dp*z-1.7282079917655628e-1_dp*z**2-  &
  5.329825490307057e-1_dp*z**3+3.1284748837287095e-1_dp*z**4+  &
  1.4498127085337966_dp*z**5))/z**5+  &
  PL1*((6.99607679578488e-2_dp*PL2*r1**2*rz**2*(0.5_dp+z)**2)/z**5  &
  +(1.0494115193677322e-1_dp*PL3*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5+(2.098823038735464e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(r1*rz*(-7.88407725694488e-3_dp-  &
  2.9978935505327797e-2_dp*z-2.3509379518035e-3_dp*z**2+  &
  1.7295897377329488e-1_dp*z**3+  &
  2.4163545142229945e-1_dp*z**4))/z**5)+  &
  PL2*((2.098823038735464e-1_dp*PL3*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5+(4.197646077470928e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(r1*rz*(-1.5768154513889758e-2_dp-  &
  5.9957871010655595e-2_dp*z-4.701875903607e-3_dp*z**2+  &
  3.4591794754658975e-1_dp*z**3+4.832709028445989e-1_dp*z**4))/z**5)  &
  +PL3*((PL4*(3.9352931976289964e-2_dp-3.148234558103197e-1_dp*z**2+  &
  6.296469116206394e-1_dp*z**4))/z**5+(1.1826115885417319e-2_dp+  &
  2.131617148715706e-2_dp*z-8.641039958827815e-2_dp*z**2-  &
  2.6649127451535284e-1_dp*z**3+1.5642374418643548e-1_dp*z**4+  &
  7.249063542668983e-1_dp*z**5)/z**5)+  &
  l1**3*((-6.996076795784882e-2_dp*l2**2*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5-(1.7490191989462205e-2_dp*l3**2*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5-  &
  (34.295003315177963_dp*(-6.275297446905026e-5_dp*r1*rz**2*(1-  &
  4._dp*(0.5_dp+r1*rz-1._dp*z)**2)**3-  &
  2.1515305532245805e-4_dp*r1**2*rz**3*(1-4._dp*(0.5_dp+r1*rz-  &
  1._dp*z)**2)**3-1.1355300142018618e-4_dp*r1*rz**2*(1-4._dp*(0.5_dp  &
  +r1*rz-1._dp*z)**2)**3*z-1.4343537021497204e-4_dp*r1**2*rz**3*(1  &
  -4._dp*(0.5_dp+r1*rz-1._dp*z)**2)**3*z-  &
  4.303061106449161e-4_dp*r1**2*rz**2*(1-4._dp*(0.5_dp+r1*rz-  &
  1._dp*z)**2)**3*z**1.5_dp+3.8249432057325876e-4_dp*r1*rz**2*(1-  &
  4._dp*(0.5_dp+r1*rz-1._dp*z)**2)**3*z**2-  &
  2.8687074042994407e-4_dp*r1**2*rz**2*(1-4._dp*(0.5_dp+r1*rz-  &
  1._dp*z)**2)**3*z**2.5_dp-1.1724466331119023e-3_dp*r1*(1-  &
  4._dp*(0.5_dp+r1*rz-1._dp*z)**2)**3*z**3-  &
  4.019817027812235e-3_dp*r1**2*rz*(1-4._dp*(0.5_dp+r1*rz-  &
  1._dp*z)**2)**3*z**3+2.8687074042994407e-4_dp*r1*rz**2*(1-  &
  4._dp*(0.5_dp+r1*rz-1._dp*z)**2)**3*z**3+1.25e-1_dp*rz*z**4+  &
  7.5e-1_dp*r1*rz**2*z**4+1.5_dp*r1**2*rz**3*z**4+r1**3*rz**4*z**4  &
  -1.3399390092707453e-3_dp*r1*(1-4._dp*(0.5_dp+r1*rz-  &
  1._dp*z)**2)**3*z**4-8.03963405562447e-3_dp*r1**2*(1-4._dp*(0.5_dp  &
  +r1*rz-1._dp*z)**2)**3*z**4.5_dp-1._dp*rz*z**5-  &
  4.5_dp*r1*rz**2*z**5-6._dp*r1**2*rz**3*z**5-  &
  2._dp*r1**3*rz**4*z**5+8.03963405562447e-3_dp*r1*(1-4._dp*(0.5_dp+  &
  r1*rz-1._dp*z)**2)**3*z**5+3.125_dp*rz*z**6+  &
  9.75_dp*r1*rz**2*z**6+7.499999999999998_dp*r1**2*rz**3*z**6+  &
  r1**3*rz**4*z**6-4.75_dp*rz*z**7-9._dp*r1*rz**2*z**7-  &
  3._dp*r1**2*rz**3*z**7+3.5_dp*rz*z**8+3._dp*r1*rz**2*z**8-  &
  1._dp*rz*z**9))/(rz*(1-4._dp*(0.5_dp+r1*rz-1._dp*z)**2)**3*z**5)  &
  +(PL1*(-1.9129897488474286e-3_dp-2.1862739986827755e-3_dp*z+  &
  2.0769602987486366e-2_dp*z**2+8.745095994731102e-3_dp*z**3-  &
  5.24705759683866e-2_dp*z**4+r1*(rz*(-6.558821996048325e-3_dp+  &
  2.62352879841933e-2_dp*z**2)+z**1.5_dp*(-1.311764399209665e-2_dp+  &
  5.24705759683866e-2_dp*z**2))))/z**5+  &
  (PL2*(-3.825979497694857e-3_dp-4.372547997365551e-3_dp*z+  &
  4.153920597497273e-2_dp*z**2+1.74901919894622e-2_dp*z**3-  &
  1.049411519367732e-1_dp*z**4+r1*(rz*(-1.311764399209665e-2_dp+  &
  5.24705759683866e-2_dp*z**2)+z**1.5_dp*(-2.62352879841933e-2_dp+  &
  1.049411519367732e-1_dp*z**2))))/z**5+  &
  l3*((PL1*(2.1862739986827755e-3_dp-1.74901919894622e-2_dp*z**2+  &
  3.49803839789244e-2_dp*z**4))/z**5+(PL2*(4.372547997365551e-3_dp-  &
  3.49803839789244e-2_dp*z**2+6.99607679578488e-2_dp*z**4))/z**5+  &
  (r1*(r1*rz**3*(-6.558821996048325e-3_dp-1.311764399209665e-2_dp*z)  &
  -4.595321276507649e-2_dp*z**3+1.8381285106030596e-1_dp*z**5+  &
  rz**2*(-4.372547997365551e-3_dp-7.651958995389714e-3_dp*z-  &
  1.311764399209665e-2_dp*r1*z**1.5_dp+1.8583328988803591e-2_dp*z**2-  &
  2.62352879841933e-2_dp*r1*z**2.5_dp+  &
  3.2794109980241632e-2_dp*z**3)))/(rz*z**5))+  &
  l2*((-6.996076795784882e-2_dp*l3*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5+(PL2*(8.745095994731102e-3_dp-  &
  6.99607679578488e-2_dp*z**2+1.399215359156976e-1_dp*z**4))/z**5+  &
  (PL1*(4.372547997365551e-3_dp-3.49803839789244e-2_dp*z**2+  &
  6.99607679578488e-2_dp*z**4))/z**5+  &
  (r1*(r1*rz**3*(-1.311764399209665e-2_dp-2.62352879841933e-2_dp*z)-  &
  9.190642553015298e-2_dp*z**3+3.676257021206119e-1_dp*z**5+  &
  rz**2*(-8.745095994731102e-3_dp-1.5303917990779428e-2_dp*z-  &
  2.62352879841933e-2_dp*r1*z**1.5_dp+3.7166657977607183e-2_dp*z**2-  &
  5.24705759683866e-2_dp*r1*z**2.5_dp+  &
  6.5588219960483265e-2_dp*z**3)))/(rz*z**5)))+  &
  l1*((-6.996076795784882e-2_dp*PL1**2*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5-(2.7984307183139516e-1_dp*PL2**2*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(PL4*r1*(5.514385531809179e-1_dp*z**3-  &
  2.2057542127236713_dp*z**5+rz**2*(2.951469898221747e-2_dp+  &
  1.9676465988144982e-2_dp*z-1.1805879592886988e-1_dp*z**2-  &
  7.870586395257993e-2_dp*z**3)))/(rz*z**5)+  &
  (PL3*r1*(2.7571927659045894e-1_dp*z**3-  &
  1.1028771063618357_dp*z**5+rz**2*(1.4757349491108736e-2_dp+  &
  9.83823299407249e-3_dp*z-5.902939796443494e-2_dp*z**2-  &
  3.9352931976289964e-2_dp*z**3)))/(rz*z**5)+  &
  l3*((3.49803839789244e-2_dp*PL1*r1**2*rz**2*(0.5_dp+z)**2)/z**5+  &
  (6.99607679578488e-2_dp*PL2*r1**2*rz**2*(0.5_dp+z)**2)/z**5+  &
  (1.0494115193677322e-1_dp*PL3*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5+(2.098823038735464e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(r1*rz*(-7.88407725694488e-3_dp-  &
  2.9978935505327797e-2_dp*z-2.3509379518035e-3_dp*z**2+  &
  1.7295897377329488e-1_dp*z**3+  &
  2.4163545142229945e-1_dp*z**4))/z**5)+  &
  l2*((6.99607679578488e-2_dp*PL1*r1**2*rz**2*(0.5_dp+z)**2)/z**5+  &
  (1.399215359156976e-1_dp*PL2*r1**2*rz**2*(0.5_dp+z)**2)/z**5+  &
  (2.098823038735464e-1_dp*PL3*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5+(4.197646077470928e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(r1*rz*(-1.5768154513889758e-2_dp-  &
  5.9957871010655595e-2_dp*z-4.701875903607e-3_dp*z**2+  &
  3.4591794754658975e-1_dp*z**3+4.832709028445989e-1_dp*z**4))/z**5)  &
  +(rz*z**3*(-7.148623368126624e-1_dp+1.3109525063667173_dp*z+  &
  3.9981229997518186e-1_dp*z**2+3.468306076973711e-1_dp*z**3-  &
  1.3427330772266073_dp*z**4)+  &
  r1**3*rz**2*(z**3*(1.6571512987986534e-1_dp+  &
  2.9869588297299714e-1_dp*z-5.4797748760677e-1_dp*z**2)+  &
  rz**2*(8.869586914062989e-3_dp+2.1900186558076458e-2_dp*z-  &
  1.8671366291378122e-2_dp*z**2-1.5547290944834784e-1_dp*z**3-  &
  9.06132942833623e-2_dp*z**4))+  &
  r1**2*rz*(z**3*(1.6571512987986534e-1_dp-3.273437678673352e-2_dp*z  &
  -1.1453692535527642_dp*z**2+1.09595497521354_dp*z**3)+  &
  rz**2*(8.869586914062989e-3_dp+4.16101272995048e-3_dp*z-  &
  6.247173940753104e-2_dp*z**2-1.1813017686559162e-1_dp*z**3+  &
  2.203325246133334e-1_dp*z**4+1.8122658856672458e-1_dp*z**5))+  &
  r1*z*(z**3*(-1.6571512987986534e-1_dp-1.3298075309313182e-1_dp*z+  &
  8.466733705797672e-1_dp*z**2-5.4797748760677e-1_dp*z**3)+  &
  rz**2*(-8.869586914062989e-3_dp-1.303059964401347e-2_dp*z-  &
  1.3891531207758698_dp*z**2-1.0074279136024535e-1_dp*z**3+  &
  2.596763157509472e-1_dp*z**4+  &
  1.2521197829432449_dp*z**5)))/(rz*z**5*(r1**2*rz**2+  &
  r1*rz*(1-2._dp*z)-1._dp*z+z**2))+  &
  PL2*((PL4*(-5.24705759683866e-2_dp+4.197646077470928e-1_dp*z**2-  &
  8.395292154941856e-1_dp*z**4))/z**5+(PL3*(-2.62352879841933e-2_dp+  &
  2.098823038735464e-1_dp*z**2-4.197646077470928e-1_dp*z**4))/z**5+  &
  (-1.5768154513889758e-2_dp-2.8421561982876082e-2_dp*z+  &
  1.1521386611770421e-1_dp*z**2+3.5532169935380375e-1_dp*z**3-  &
  2.0856499224858065e-1_dp*z**4-9.665418056891978e-1_dp*z**5+  &
  r1**2*((-3.676257021206119e-1_dp-7.352514042412238e-1_dp*z)*z**3+  &
  rz**2*(-1.9676465988144982e-2_dp-5.24705759683866e-2_dp*z-  &
  2.62352879841933e-2_dp*z**2)))/z**5)+  &
  PL1*((-2.7984307183139516e-1_dp*PL2*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5+(PL4*(-2.62352879841933e-2_dp+  &
  2.098823038735464e-1_dp*z**2-4.197646077470928e-1_dp*z**4))/z**5+  &
  (PL3*(-1.311764399209665e-2_dp+1.049411519367732e-1_dp*z**2-  &
  2.098823038735464e-1_dp*z**4))/z**5+(-7.88407725694488e-3_dp-  &
  1.4210780991438041e-2_dp*z+5.76069330588521e-2_dp*z**2+  &
  1.7766084967690188e-1_dp*z**3-1.0428249612429032e-1_dp*z**4-  &
  4.832709028445989e-1_dp*z**5+r1**2*((-1.8381285106030596e-1_dp-  &
  3.676257021206119e-1_dp*z)*z**3+rz**2*(-9.83823299407249e-3_dp-  &
  2.62352879841933e-2_dp*z-1.311764399209665e-2_dp*z**2)))/z**5))+  &
  l1**2*((6.99607679578488e-2_dp*l2**2*r1**2*rz**2*(0.5_dp+  &
  z)**2)/z**5+(1.74901919894622e-2_dp*l3**2*r1**2*rz**2*(0.5_dp+  &
  z)**2)/z**5+(PL1**2*(4.372547997365551e-3_dp-  &
  3.49803839789244e-2_dp*z**2+6.99607679578488e-2_dp*z**4))/z**5+  &
  (PL2**2*(1.74901919894622e-2_dp-1.399215359156976e-1_dp*z**2+  &
  2.798430718313952e-1_dp*z**4))/z**5+  &
  (PL2*r1*(r1*rz**3*(-1.311764399209665e-2_dp-  &
  2.62352879841933e-2_dp*z)-3.676257021206119e-1_dp*z**3+  &
  1.4705028084824476_dp*z**5+rz**2*(-2.3502445485839836e-2_dp-  &
  2.5142150984851916e-2_dp*z-2.62352879841933e-2_dp*r1*z**1.5_dp+  &
  9.619605594204211e-2_dp*z**2-5.24705759683866e-2_dp*r1*z**2.5_dp+  &
  1.049411519367732e-1_dp*z**3)))/(rz*z**5)+  &
  (PL4*(5.738969246542284e-3_dp+6.5588219960483265e-3_dp*z-  &
  6.2308808962459095e-2_dp*z**2-2.6235287984193305e-2_dp*z**3+  &
  1.5741172790515985e-1_dp*z**4+  &
  r1*(z**1.5_dp*(3.9352931976289964e-2_dp-  &
  1.5741172790515985e-1_dp*z**2)+rz*(1.9676465988144982e-2_dp-  &
  7.870586395257993e-2_dp*z**2))))/z**5+  &
  (PL3*(2.869484623271142e-3_dp+3.2794109980241632e-3_dp*z-  &
  3.1154404481229547e-2_dp*z**2-1.3117643992096653e-2_dp*z**3+  &
  7.870586395257993e-2_dp*z**4+  &
  r1*(z**1.5_dp*(1.9676465988144982e-2_dp-  &
  7.870586395257993e-2_dp*z**2)+rz*(9.83823299407249e-3_dp-  &
  3.9352931976289964e-2_dp*z**2))))/z**5+  &
  PL1*((PL2*(1.74901919894622e-2_dp-1.399215359156976e-1_dp*z**2+  &
  2.798430718313952e-1_dp*z**4))/z**5+  &
  (r1*(r1*rz**3*(-6.558821996048325e-3_dp-1.311764399209665e-2_dp*z)  &
  -1.8381285106030596e-1_dp*z**3+7.352514042412238e-1_dp*z**5+  &
  rz**2*(-1.1751222742919918e-2_dp-1.2571075492425958e-2_dp*z-  &
  1.311764399209665e-2_dp*r1*z**1.5_dp+4.809802797102106e-2_dp*z**2-  &
  2.62352879841933e-2_dp*r1*z**2.5_dp+  &
  5.24705759683866e-2_dp*z**3)))/(rz*z**5))+  &
  (r1**8*rz**6*((2.067894574428442e-1_dp+  &
  1.3785963829522947e-1_dp*z)*z**3+rz**2*(5.534006059165773e-3_dp+  &
  7.378674745554367e-3_dp*z+2.4595582485181224e-3_dp*z**2))+  &
  r1**7*rz**5*(z**3*(6.203683723285325e-1_dp-  &
  8.271578297713765e-1_dp*z-8.271578297713768e-1_dp*z**2)+  &
  rz*z**1.5_dp*(1.1826115885417319e-2_dp+2.131617148715706e-2_dp*z-  &
  3.910593604660887e-2_dp*z**2-1.8122658856672458e-1_dp*z**3)+  &
  rz**2*(2.2515076120205983e-2_dp-4.0992637475301805e-4_dp*z-  &
  5.644634175107628e-2_dp*z**2-1.05370643774471e-1_dp*z**3))+  &
  r1**6*rz**4*(z**3*(6.203683723285325e-1_dp-  &
  2.688262946756974_dp*z+1.0339472872142208_dp*z**2+  &
  2.067894574428442_dp*z**3)+rz*z**1.5_dp*(3.5478347656251956e-2_dp  &
  -7.008180851032737e-3_dp*z-2.4521483706276896e-1_dp*z**2-  &
  3.090441494205205e-1_dp*z**3+1.0873595314003475_dp*z**4)+  &
  rz**2*(3.606583390558e-2_dp-5.9298529420226584e-2_dp*z-  &
  1.5687514372601394e-1_dp*z**2-1.349980322489558e-1_dp*z**3+  &
  5.89474644046767e-1_dp*z**4+1.8122658856672458e-1_dp*z**5))+  &
  z**3*(-1.7246418999566924e-3_dp+9.429804375678598e-5_dp*z+  &
  1.3973695628629281e-1_dp*z**2-6.953544212480171e-1_dp*z**3+  &
  1.2636666694692562_dp*z**4-6.99884105612706e-1_dp*z**5-  &
  5.7854811751174235e-1_dp*z**6+8.535692914002645e-1_dp*z**7-  &
  2.815559289271477e-1_dp*z**8+  &
  rz*z**1.5_dp*(-1.5516080288480205e-1_dp+1.0298946905183168_dp*z-  &
  2.1024883245708406_dp*z**2+6.071112253981175e-1_dp*z**3+  &
  2.2576491274556427_dp*z**4-1.6370059159164343_dp*z**5)+  &
  rz**2*(1.5516080288480205e-1_dp-1.0298946905183168_dp*z+  &
  2.1024883245708406_dp*z**2-6.071112253981175e-1_dp*z**3-  &
  2.2576491274556427_dp*z**4+1.6370059159164343_dp*z**5))+  &
  r1**5*rz**3*(z**3*(2.067894574428442e-1_dp-  &
  2.3436138510189006_dp*z+4.5493680637425715_dp*z**2-  &
  2.7571927659045894_dp*z**4)+  &
  rz*z**1.5_dp*(3.5478347656251956e-2_dp-1.1344322381978862e-1_dp*z-  &
  2.596686421659227e-1_dp*z**2+3.626518473063152e-1_dp*z**3+  &
  2.1318097878017355_dp*z**4-2.7183988285008684_dp*z**5)+  &
  rz**2*(2.8447105587161827e-2_dp-1.1085997830572953e-1_dp*z-  &
  1.2230555287586704e-1_dp*z**2+1.8356144873766764e-1_dp*z**3+  &
  1.344398890267224_dp*z**4-9.181198412335519e-1_dp*z**5-  &
  1.0873595314003475_dp*z**6))+  &
  r1**4*rz**2*(z**4*(-6.203683723285325e-1_dp+  &
  3.308631319085507_dp*z-3.722210233971195_dp*z**2-  &
  1.033947287214221_dp*z**3+2.067894574428442_dp*z**4)+  &
  rz*z**1.5_dp*(1.1826115885417319e-2_dp-1.2059721913785078e-1_dp*z+  &
  5.988348267002598e-2_dp*z**2+6.910074708989473e-1_dp*z**3+  &
  5.752175516592876e-1_dp*z**4-4.654678936069561_dp*z**5+  &
  3.6245317713344916_dp*z**6)+rz**2*(1.1086983642578736e-2_dp-  &
  8.753137327743323e-2_dp*z+1.5163930887920032e-2_dp*z**2+  &
  4.279728688732304e-1_dp*z**3+3.617476994346818e-1_dp*z**4-  &
  2.231446955308871_dp*z**5-2.7293448232965707e-1_dp*z**6+  &
  2.255616311006996_dp*z**7))+  &
  r1**3*rz*(z**5*(6.203683723285325e-1_dp-2.0678945744284416_dp*z+  &
  1.4475262020999093_dp*z**2+8.271578297713768e-1_dp*z**3-  &
  8.271578297713768e-1_dp*z**4)+  &
  rz*z**2.5_dp*(-3.5478347656251956e-2_dp+1.4892157147604057e-1_dp*z+  &
  1.4622541834613412e-1_dp*z**2-6.223204894722381e-1_dp*z**3-  &
  1.7691579404954203_dp*z**4+4.850208616302604_dp*z**5-  &
  2.7183988285008684_dp*z**6)+rz**2*(1.7246418999566924e-3_dp-  &
  3.3355248971492997e-2_dp*z+6.787015939650393e-2_dp*z**2-  &
  5.013122264757835e-2_dp*z**3+3.434358805567287e-1_dp*z**4-  &
  6.605520133844134e-1_dp*z**5-7.694053954306895e-1_dp*z**6+  &
  2.745121611836966_dp*z**7-1.7734017013590027_dp*z**8))+  &
  r1**2*z*(rz**3*z**3.5_dp*(-6.206432115392082e-1_dp+  &
  1.6370059159164343_dp*z+6.206432115392082e-1_dp*z**2-  &
  1.6370059159164343_dp*z**3)+rz**4*z**2*(6.206432115392082e-1_dp-  &
  1.6370059159164343_dp*z-6.206432115392082e-1_dp*z**2+  &
  1.6370059159164343_dp*z**3)+z**5*(-2.067894574428442e-1_dp+  &
  4.82508734033303e-1_dp*z-2.067894574428441e-1_dp*z**2-  &
  2.067894574428442e-1_dp*z**3+1.3785963829522947e-1_dp*z**4)+  &
  rz*z**2.5_dp*(3.5478347656251956e-2_dp-7.796487616353664e-2_dp*z-  &
  1.9572012770445157e-1_dp*z**2+1.7437734385398471e-1_dp*z**3+  &
  1.4602329931786189_dp*z**4-2.4837632122212154_dp*z**5+  &
  1.0873595314003475_dp*z**6)+rz**2*(-5.173925699870077e-3_dp+  &
  3.3543845059006565e-2_dp*z-2.732149001661349e-1_dp*z**2+  &
  1.1773687739492913_dp*z**3-1.2674879771538905_dp*z**4-  &
  2.005013536505061_dp*z**5+4.876031382656149_dp*z**6-  &
  2.4777573856771244_dp*z**7-5.829627646236535e-2_dp*z**8))+  &
  r1*z**2*(rz**3*z*(6.206432115392082e-1_dp-2.878292338994851_dp*z  &
  +2.6533686202936604_dp*z**2+2.878292338994851_dp*z**3-  &
  3.2740118318328686_dp*z**4)+  &
  rz**2*z**2.5_dp*(-6.206432115392082e-1_dp+2.878292338994851_dp*z-  &
  2.6533686202936604_dp*z**2-2.878292338994851_dp*z**3+  &
  3.2740118318328686_dp*z**4)+z**2.5_dp*(-1.1826115885417319e-2_dp+  &
  1.4162176169094898e-2_dp*z+6.757610285182808e-2_dp*z**2+  &
  1.178638185084411e-2_dp*z**3-4.050457860731901e-1_dp*z**4+  &
  5.045738296535649e-1_dp*z**5-1.8122658856672458e-1_dp*z**6)+  &
  rz*(5.173925699870077e-3_dp-8.895027921625012e-2_dp*z+  &
  4.84478594547966e-1_dp*z**2-7.039239940781732e-1_dp*z**3-  &
  9.682849779437941e-1_dp*z**4+3.3811430503505404_dp*z**5-  &
  2.4569641054175904_dp*z**6-4.1644275251771035e-1_dp*z**7+  &
  7.637705385751417e-1_dp*z**8)))/(z**5*(r1**2*rz**2+r1*rz*(1-  &
  2._dp*z)+(-1._dp+z)*z)**3)+  &
  l2*((6.99607679578488e-2_dp*l3*r1**2*rz**2*(0.5_dp+z)**2)/z**5-  &
  (1.7490191989462203e-1_dp*PL1*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5-(3.498038397892441e-1_dp*PL2*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(PL4*(-1.311764399209665e-2_dp+  &
  1.049411519367732e-1_dp*z**2-2.098823038735464e-1_dp*z**4))/z**5+  &
  (PL3*(-6.558821996048325e-3_dp+5.24705759683866e-2_dp*z**2-  &
  1.049411519367732e-1_dp*z**4))/z**5+(-3.94203862847244e-3_dp-  &
  7.10539049571902e-3_dp*z+2.880346652942605e-2_dp*z**2+  &
  8.883042483845093e-2_dp*z**3-5.2141248062145165e-2_dp*z**4-  &
  2.4163545142229945e-1_dp*z**5+r1**2*((-3.676257021206119e-1_dp-  &
  7.352514042412238e-1_dp*z)*z**3+rz**2*(-1.9676465988144982e-2_dp-  &
  5.24705759683866e-2_dp*z-2.62352879841933e-2_dp*z**2)))/z**5)+  &
  l3*((-8.745095994731102e-2_dp*PL1*r1*rz*(-0.5_dp+z)*(0.5_dp+  &
  z)**2)/z**5-(1.7490191989462203e-1_dp*PL2*r1*rz*(-0.5_dp+  &
  z)*(0.5_dp+z)**2)/z**5+(PL4*(-6.558821996048325e-3_dp+  &
  5.24705759683866e-2_dp*z**2-1.049411519367732e-1_dp*z**4))/z**5+  &
  (PL3*(-3.2794109980241624e-3_dp+2.62352879841933e-2_dp*z**2-  &
  5.24705759683866e-2_dp*z**4))/z**5+(-1.97101931423622e-3_dp-  &
  3.55269524785951e-3_dp*z+1.4401733264713026e-2_dp*z**2+  &
  4.441521241922547e-2_dp*z**3-2.6070624031072582e-2_dp*z**4-  &
  1.2081772571114973e-1_dp*z**5+r1**2*((-1.8381285106030596e-1_dp-  &
  3.676257021206119e-1_dp*z)*z**3+rz**2*(-9.83823299407249e-3_dp-  &
  2.62352879841933e-2_dp*z-1.311764399209665e-2_dp*z**2)))/z**5))+  &
  l1**4*((l3**2*(2.732842498353469e-4_dp-  &
  2.1862739986827755e-3_dp*z**2+4.372547997365551e-3_dp*z**4))/z**5+  &
  (l2**2*(1.0931369993413877e-3_dp-8.745095994731102e-3_dp*z**2+  &
  1.74901919894622e-2_dp*z**4))/z**5+(2.092332537801875e-4_dp+  &
  4.782474372118571e-4_dp*z-2.5962003734357957e-3_dp*z**2-  &
  3.2794109980241632e-3_dp*z**3+9.83823299407249e-3_dp*z**4+  &
  r1**2*(2.4595582485181224e-3_dp*rz**2+  &
  9.83823299407249e-3_dp*rz*z**1.5_dp+9.83823299407249e-3_dp*z**3)+  &
  r1*(z**1.5_dp*(2.869484623271142e-3_dp+3.2794109980241632e-3_dp*z-  &
  1.9676465988144982e-2_dp*z**2)+rz*(1.434742311635571e-3_dp+  &
  1.6397054990120816e-3_dp*z-9.83823299407249e-3_dp*z**2)))/z**5+  &
  (l3*(-4.782474372118571e-4_dp-5.465684996706938e-4_dp*z+  &
  5.192400746871591e-3_dp*z**2+2.1862739986827755e-3_dp*z**3-  &
  1.311764399209665e-2_dp*z**4+r1*(rz*(-1.6397054990120812e-3_dp+  &
  6.558821996048325e-3_dp*z**2)+z**1.5_dp*(-3.2794109980241624e-3_dp+  &
  1.311764399209665e-2_dp*z**2))))/z**5+  &
  l2*((l3*(1.0931369993413877e-3_dp-8.745095994731102e-3_dp*z**2+  &
  1.74901919894622e-2_dp*z**4))/z**5+(-9.564948744237142e-4_dp-  &
  1.0931369993413877e-3_dp*z+1.0384801493743183e-2_dp*z**2+  &
  4.372547997365551e-3_dp*z**3-2.62352879841933e-2_dp*z**4+  &
  r1*(rz*(-3.2794109980241624e-3_dp+1.311764399209665e-2_dp*z**2)+  &
  z**1.5_dp*(-6.558821996048325e-3_dp+  &
  2.62352879841933e-2_dp*z**2)))/z**5))

end function P2

!ccccccccccccccc

recursive function fact(i) result(prod)
  integer      , intent(in) :: i
  real (dp)                 :: prod

  select case(i)
  case(:-1);     prod = 0
  case(0:1);     prod = 1
  case default;  prod = i * fact(i - 1)
  end select

end function fact

end module SigmaClass
