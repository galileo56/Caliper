
module SigmaClass
  use AnomDimClass;  use RunningClass;  use ElectroWeakClass; use Chaplin
  use Constants, only: dp, Pi, Pi2, Zeta3;  implicit none;  private

  real (dp), dimension(2)               :: EWfact
  public                                :: setEWfact, Pi0, Pi0Der, Pi1Der, &
  P2, P2Der, Pi3, Pi1

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
    Rhad, SetMCharm, RHadCoefs, RQCD

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
    if (order > 1)   logList(1:) = PowList( log(mu/Q), n )
    if (order > 0) alphaList(1:) = PowList( self%run%alphaQCD(mu)/Pi, n + 1 )

    alpha_loop: do i = 2, min(order,4)
      rQ(i) = dot_product( self%RhadCoef(:n,i), logList )
    end do alpha_loop

    Rhad = dot_product( alphaList, rQ )

  end function Rhad

!ccccccccccccccc

  real (dp) function RQCD(self, order, gt, h, Q)
    class (Sigma)      , intent(in)          :: self
    Integer            , intent(in)          :: order
    real (dp)          , intent(in)          :: h, Q, gt
    complex (dp)                             :: z
    real (dp)                                :: mu
    integer                                  :: i, n
    real (dp), dimension(0:min(order,3))     :: alphaList, rQ
    real (dp), dimension(0:min(order,3) - 1) :: logList
    real (dp), dimension(0:min(order,3) - 1, min(order,3)) :: Rcoef

    mu = h * self%m; z = ( Q + (0,1) * gt )**2/4/self%m**2 ; n = min(order,3)

    Rcoef= 0; rQ = PiCoef(z,n); Rcoef(0,:n) = rQ(1:); rQ(1:) = 0

    logList(0) = 1; alphaList(0) = 1
    if (order > 0) alphaList(1:) = PowList( self%run%alphaQCD(mu)/Pi, n )
    n = n + 1;  if (order > 1) logList(1:) = PowList( log(mu/Q), n )

    call self%andim%expandAlpha(Rcoef)

    do i = 1, n
      rQ(i) = dot_product( Rcoef(:n,i), logList )
    end do

    RQCD = rQ(0) + dot_product( alphaList, rQ )

  end function RQCD

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

  function PiCoef(z, n) result(res)
  complex (dp), intent(in)           :: z
  integer     , intent(in)           :: n
  real (dp), dimension( 0:min(n,3) ) :: res

  if (n >= 0) res(0) = ImagPart( Pi0(z) )
  if (n >= 1) res(2) = ImagPart( Pi1(z) )
  if (n >= 2) res(3) = ImagPart( P2(z)  )
  if (n >= 3) res(4) = ImagPart( Pi3(z) )

  end function PiCoef

!ccccccccccccccc

complex (dp) function Pi1Der(i, z)
  integer     , intent(in) :: i
  complex (dp), intent(in) :: z
  complex (dp)             :: uu, uuDer, uuDer2, uuDer3

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
  PL4, l4, d1, d2, l5, z05, z05m

  rz = sqrt(z); r1 = sqrt(1 - z); r2 = 1 + r1; r3 = 1 - r1; r0 = sqrt(z - 1)
  d2 = 1 + 2 * r0 * rz - 2 * z; d1 = d2**2 - 1; l1 = Log(d2); PL2 = cli2(-d2)
  l2 = Log(1 + r0 * rz - z); l3 = Log( 8 * (z - r0 * rz) ); PL1 = cli2(d2)
  PL3 = cli3(d2); PL4 = cli3(-d2); l4 = Log(1 + d2); z05 = 0.5_dp + z
  z05m = z - 0.5_dp

  r32 = PowList(r3/r2, 12); l5 = Log(2 * z - 2 * r0 * rz)

  P2Der = (2 + 1/r1)/(r2 - z)**2 * ( 0.688472275231591_dp + &
  0.1891443911186843_dp * r32(1) - 0.030243753917029492_dp * r32(2) + &
  0.0034204077615457664_dp * r32(3) + 0.006118445792320977_dp * r32(4) &
  + 0.0016455445946186606_dp * r32(5) + 0.001252836593493741_dp * r32(6) &
  + 0.00046831306017028774_dp * r32(7) + 0.0018087001167161186_dp * r32(8) &
  + 0.00048205327014239514_dp * r32(9) - 0.0024098794332697245_dp * r32(10) &
  - 0.0001825139597445611_dp * r32(11) + 0.0013086431124221074_dp * r32(12) )&
  + 2/(r2 - z)/r1/r2 * ( 0.09457219555934215_dp + 0.06432844164231266_dp * r32(1) - &
  0.025113142274710842_dp * r32(2) + 0.017367503226960603_dp * r32(3) + &
  0.016350753071188606_dp * r32(4) + 0.007872371267027875_dp * r32(5) + &
  0.00539760549107723_dp * r32(6) + 0.008873896177460481_dp * r32(7) + &
  0.009404040182505253_dp * r32(8) - 0.009880157450707843_dp * r32(9) - &
  0.013053223944943708_dp * r32(10) + 0.006848031895937559_dp * r32(11) + &
  0.007851858674532645_dp * r32(12) ) + &
(2.6854661544532146_dp*(r0**2-2*r0*rz+rz**2)*(-1._dp+  &
z)*(1.0647869616635923_dp+1.2416980235463153_dp*z+  &
z**2))/(d1*r0*rz*z**2)+(PL4**2*(-1.9676465988144982e-1_dp+  &
9.44470367430959e-1_dp*z**2-6.296469116206391e-1_dp*z**4))/z**6+  &
(-1.7769580297039034e-2_dp-5.124650217014172e-2_dp*z+  &
3.58725619023855e-2_dp*z**2-3.1775183588661764_dp*z**3-  &
5.55408076785252e-1_dp*z**4)/z**6+  &
(PL3**2*(-4.919116497036246e-2_dp+2.3611759185773975e-1_dp*z**2-  &
1.5741172790515978e-1_dp*z**4))/z**6+  &
PL4*((-1.1826115885417319e-1_dp-1.7052937189725645e-1_dp*z+  &
5.18462397529669e-1_dp*z**2+1.0659650980614113_dp*z**3-  &
3.1284748837287113e-1_dp*z**4)/z**6-  &
(7.870586395257993e-2_dp*(r0**2-2*r0*rz+rz**2)*(-2.5e-1_dp+  &
z**2)*(2.802528429206536e1_dp*z**3+rz**2*(1.5_dp+  &
z)))/(d2*rz**2*z**5))+PL3*((PL4*(-1.9676465988144982e-1_dp+  &
9.44470367430959e-1_dp*z**2-6.296469116206391e-1_dp*z**4))/z**6+  &
(-5.91305794270866e-2_dp-8.526468594862823e-2_dp*z+  &
2.592311987648345e-1_dp*z**2+5.329825490307057e-1_dp*z**3-  &
1.5642374418643556e-1_dp*z**4)/z**6-  &
(3.9352931976289964e-2_dp*(r0**2-2*r0*rz+rz**2)*(-2.5e-1_dp+  &
z**2)*(2.802528429206536e1_dp*z**3+rz**2*(1.5_dp+  &
z)))/(d2*rz**2*z**5))+PL2**2*((-2.798430718313952e-1_dp*(r0**2-  &
2*r0*rz+rz**2)*z05m*z05**2)/(d2*z**5)+((1+  &
2*z)*(6.99607679578488e-2_dp*r0**2*rz**2*z-  &
1.7490191989462203e-1_dp*r0**2*rz**2*z05+  &
3.49803839789244e-2_dp*r0**2*z*z05+  &
3.49803839789244e-2_dp*rz**2*z*z05+  &
(1.049411519367732e-1_dp*(r0**2-2*r0*rz+rz**2)*z*(-2.5e-1_dp+  &
z**2))/(r0*rz - z05m)))/z**6)+  &
(r0**2*(z**3*(1.6571512987986532e-1_dp+2.986958829729971e-1_dp*z-  &
5.4797748760677e-1_dp*z**2)+rz**2*(8.869586914062989e-3_dp+  &
2.1900186558076458e-2_dp*z-1.8671366291378115e-2_dp*z**2-  &
1.5547290944834788e-1_dp*z**3-9.06132942833623e-2_dp*z**4))+  &
rz**2*(z**3*(1.6571512987986532e-1_dp+2.986958829729971e-1_dp*z-  &
5.4797748760677e-1_dp*z**2)+rz**2*(8.869586914062989e-3_dp+  &
2.1900186558076458e-2_dp*z-1.8671366291378115e-2_dp*z**2-  &
1.5547290944834788e-1_dp*z**3-9.06132942833623e-2_dp*z**4))+  &
r0*rz*(z**3*(-3.3143025975973064e-1_dp-5.973917659459942e-1_dp*z+  &
1.09595497521354_dp*z**2)+rz**2*(-1.7739173828125978e-2_dp-  &
4.3800373116152915e-2_dp*z+3.734273258275623e-2_dp*z**2+  &
3.1094581889669577e-1_dp*z**3+  &
1.8122658856672458e-1_dp*z**4)))/(d2*rz**2*z**5)+  &
l3*((3.49803839789244e-2_dp*PL1*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)+  &
(6.99607679578488e-2_dp*PL2*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)+  &
(1.049411519367732e-1_dp*PL3*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/(d2*z**5)+  &
(2.098823038735464e-1_dp*PL4*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/(d2*z**5)+(r0*rz*(1.5768154513889756e-2_dp+  &
5.9957871010655595e-2_dp*z+4.701875903607003e-3_dp*z**2-  &
3.459179475465897e-1_dp*z**3-4.832709028445988e-1_dp*z**4)+  &
r0**2*(-7.884077256944878e-3_dp-2.9978935505327797e-2_dp*z-  &
2.3509379518035014e-3_dp*z**2+1.7295897377329486e-1_dp*z**3+  &
2.416354514222994e-1_dp*z**4)+rz**2*(-7.884077256944878e-3_dp-  &
2.9978935505327797e-2_dp*z-2.3509379518035014e-3_dp*z**2+  &
1.7295897377329486e-1_dp*z**3+  &
2.416354514222994e-1_dp*z**4))/(d2*z**5))+  &
l5*((-3.49803839789244e-2_dp*PL1*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)-  &
(6.99607679578488e-2_dp*PL2*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)-  &
(1.049411519367732e-1_dp*PL3*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/(d2*z**5)-  &
(2.098823038735464e-1_dp*PL4*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/(d2*z**5)+(r0**2*(7.884077256944878e-3_dp+  &
2.9978935505327797e-2_dp*z+2.3509379518035014e-3_dp*z**2-  &
1.7295897377329486e-1_dp*z**3-2.416354514222994e-1_dp*z**4)+  &
rz**2*(7.884077256944878e-3_dp+2.9978935505327797e-2_dp*z+  &
2.3509379518035014e-3_dp*z**2-1.7295897377329486e-1_dp*z**3-  &
2.416354514222994e-1_dp*z**4)+r0*rz*(-1.5768154513889756e-2_dp-  &
5.9957871010655595e-2_dp*z-4.701875903607003e-3_dp*z**2+  &
3.459179475465897e-1_dp*z**3+  &
4.832709028445988e-1_dp*z**4))/(d2*z**5))+  &
l4*((-3.49803839789244e-2_dp*PL1*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/((r0*rz - z05m)*z**5)-  &
(6.99607679578488e-2_dp*PL2*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/((r0*rz - z05m)*z**5)-  &
(1.049411519367732e-1_dp*PL3*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/((r0*rz - z05m)*z**5)-  &
(2.098823038735464e-1_dp*PL4*(r0**2-2*r0*rz+rz**2)*z05 &
*(-2.5e-1_dp+z**2))/((r0*rz - z05m)*z**5)+  &
(r0**2*(7.884077256944878e-3_dp+2.9978935505327797e-2_dp*z+  &
2.3509379518035014e-3_dp*z**2-1.7295897377329486e-1_dp*z**3-  &
2.416354514222994e-1_dp*z**4)+rz**2*(7.884077256944878e-3_dp+  &
2.9978935505327797e-2_dp*z+2.3509379518035014e-3_dp*z**2-  &
1.7295897377329486e-1_dp*z**3-2.416354514222994e-1_dp*z**4)+  &
r0*rz*(-1.5768154513889756e-2_dp-5.9957871010655595e-2_dp*z-  &
4.701875903607003e-3_dp*z**2+3.459179475465897e-1_dp*z**3+  &
4.832709028445988e-1_dp*z**4))/((r0*rz - z05m)*z**5))+  &
l2*((6.99607679578488e-2_dp*PL1*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)+  &
(1.399215359156976e-1_dp*PL2*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)+  &
(2.098823038735464e-1_dp*PL3*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/(d2*z**5)+  &
(4.197646077470928e-1_dp*PL4*(r0**2-2*r0*rz+rz**2)*z05  &
*(-2.5e-1_dp+z**2))/(d2*z**5)+(r0*rz*(3.153630902777951e-2_dp+  &
1.199157420213112e-1_dp*z+9.403751807214006e-3_dp*z**2-  &
6.918358950931794e-1_dp*z**3-9.665418056891976e-1_dp*z**4)+  &
r0**2*(-1.5768154513889756e-2_dp-5.9957871010655595e-2_dp*z-  &
4.701875903607003e-3_dp*z**2+3.459179475465897e-1_dp*z**3+  &
4.832709028445988e-1_dp*z**4)+rz**2*(-1.5768154513889756e-2_dp-  &
5.9957871010655595e-2_dp*z-4.701875903607003e-3_dp*z**2+  &
3.459179475465897e-1_dp*z**3+  &
4.832709028445988e-1_dp*z**4))/(d2*z**5))+  &
PL1**2*((-6.996076795784882e-2_dp*r0*rz*z05m*z05**2  &
+r0**2*(-4.372547997365551e-3_dp-8.745095994731102e-3_dp*z+  &
1.74901919894622e-2_dp*z**2+3.49803839789244e-2_dp*z**3)+  &
rz**2*(-4.372547997365551e-3_dp-8.745095994731102e-3_dp*z+  &
1.74901919894622e-2_dp*z**2+  &
3.49803839789244e-2_dp*z**3))/(d2*z**5)+  &
(1.74901919894622e-2_dp*rz**2*z*z05**2+  &
r0**2*(1.74901919894622e-2_dp*z*z05**2+  &
rz**2*(-2.1862739986827755e-2_dp-6.99607679578488e-2_dp*z-  &
5.247057596838661e-2_dp*z**2)))/z**6)+  &
PL2*((rz**2*(-1.5768154513889756e-2_dp-2.842156198287608e-2_dp*z+  &
1.152138661177042e-1_dp*z**2+3.5532169935380375e-1_dp*z**3-  &
2.085649922485806e-1_dp*z**4-9.665418056891976e-1_dp*z**5)+  &
r0*rz*(3.153630902777951e-2_dp+5.684312396575216e-2_dp*z-  &
2.304277322354084e-1_dp*z**2-7.106433987076075e-1_dp*z**3+  &
4.171299844971612e-1_dp*z**4+1.9330836113783954_dp*z**5)+  &
r0**4*((-3.6762570212061183e-1_dp-7.352514042412237e-1_dp*z)*z**3+  &
rz**2*(-1.9676465988144982e-2_dp-5.247057596838662e-2_dp*z-  &
2.62352879841933e-2_dp*z**2))+r0**2*(-1.5768154513889756e-2_dp-  &
2.842156198287608e-2_dp*z+1.152138661177042e-1_dp*z**2+  &
3.5532169935380375e-1_dp*z**3+rz**2*(-3.6762570212061183e-1_dp-  &
7.352514042412237e-1_dp*z)*z**3-2.085649922485806e-1_dp*z**4-  &
9.665418056891976e-1_dp*z**5+rz**4*(-1.9676465988144982e-2_dp-  &
5.247057596838662e-2_dp*z-2.62352879841933e-2_dp*z**2))+  &
r0**3*(rz*z**3*(7.352514042412237e-1_dp+1.4705028084824474_dp*z)  &
+rz**3*(3.9352931976289964e-2_dp+1.0494115193677325e-1_dp*z+  &
5.24705759683866e-2_dp*z**2)))/(d2*r0*rz*z**5)+  &
PL3*((-2.62352879841933e-2_dp*(-2._dp+r0/rz+rz/r0)*(1-2*z)*(1  &
+2*z)*(1-4*z**2))/(d2*z**5)+  &
(4.197646077470928e-1_dp*r0*rz*z**2*z05+  &
(3.9352931976289964e-2_dp*(2-r0/rz-rz/r0)*z*(1-  &
4*z**2)**2)/(-1._dp-2*r0*rz+2*z)+  &
2.098823038735464e-1_dp*r0*rz*z*(-2.5e-1_dp+z**2)-  &
1.0494115193677322_dp*r0*rz*z05*(-2.5e-1_dp+z**2)+  &
(1.049411519367732e-1_dp*r0*z*z05*(-2.5e-1_dp+z**2))/rz+  &
(1.049411519367732e-1_dp*rz*z*z05*(-2.5e-1_dp+  &
z**2))/r0)/z**6)+PL4*((-5.24705759683866e-2_dp*(-2._dp+r0/rz+  &
rz/r0)*(1-2*z)*(1+2*z)*(1-4*z**2))/(d2*z**5)+  &
(8.395292154941856e-1_dp*r0*rz*z**2*z05+  &
(7.870586395257993e-2_dp*(2-r0/rz-rz/r0)*z*(1-  &
4*z**2)**2)/(-1._dp-2*r0*rz+2*z)+  &
4.197646077470928e-1_dp*r0*rz*z*(-2.5e-1_dp+z**2)-  &
2.0988230387354645_dp*r0*rz*z05*(-2.5e-1_dp+z**2)+  &
(2.098823038735464e-1_dp*r0*z*z05*(-2.5e-1_dp+z**2))/rz+  &
(2.098823038735464e-1_dp*rz*z*z05*(-2.5e-1_dp+  &
z**2))/r0)/z**6)+(rz**2*z*(7.884077256944878e-3_dp+  &
1.4210780991438035e-2_dp*z-5.7606933058852094e-2_dp*z**2-  &
1.7766084967690188e-1_dp*z**3+1.042824961242903e-1_dp*z**4+  &
4.83270902844599e-1_dp*z**5)+r0*rz*z*(-2.3652231770834637e-2_dp-  &
4.263234297431412e-2_dp*z+1.7282079917655626e-1_dp*z**2+  &
5.329825490307057e-1_dp*z**3-3.128474883728709e-1_dp*z**4-  &
1.4498127085337966_dp*z**5+rz**2*(-7.884077256944878e-3_dp-  &
2.9978935505327797e-2_dp*z-2.3509379518035014e-3_dp*z**2+  &
1.7295897377329486e-1_dp*z**3+2.416354514222994e-1_dp*z**4))+  &
r0**3*rz*(rz**2*(7.884077256944879e-2_dp+2.398314840426224e-1_dp*z  &
+1.4105627710820967e-2_dp*z**2-6.918358950931797e-1_dp*z**3-  &
4.832709028445988e-1_dp*z**4)+z*(-7.884077256944878e-3_dp-  &
2.9978935505327797e-2_dp*z-2.3509379518035014e-3_dp*z**2+  &
1.7295897377329486e-1_dp*z**3+2.416354514222994e-1_dp*z**4))+  &
r0**2*(rz**2*(3.9420386284724396e-2_dp+4.10749694518624e-2_dp*z-  &
2.327786701872119e-1_dp*z**2-3.600235752574108e-1_dp*z**3+  &
4.502004436708803e-1_dp*z**4+4.832709028445988e-1_dp*z**5)+  &
z*(7.884077256944878e-3_dp+1.4210780991438035e-2_dp*z-  &
5.7606933058852094e-2_dp*z**2-1.7766084967690188e-1_dp*z**3+  &
1.042824961242903e-1_dp*z**4+  &
4.83270902844599e-1_dp*z**5)))/(r0*rz*(0.5_dp+r0*rz-  &
z)*z**6))+l1**4*((l2**2*(-5.465684996706939e-3_dp+  &
2.6235287984193305e-2_dp*z**2-1.749019198946221e-2_dp*z**4))/z**6+  &
(l3**2*(-1.3664212491767347e-3_dp+6.5588219960483265e-3_dp*z**2-  &
4.372547997365553e-3_dp*z**4))/z**6+(-1.0461662689009374e-3_dp-  &
2.391237186059285e-3_dp*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)-1.3664212491767347e-3_dp*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)**2+(9.56494874423714e-4_dp*(-1._dp+  &
r0/2/rz+rz/2/r0)*z*(-1._dp+4*z**2))/(1+r0*rz-  &
z)+(4.78247437211857e-4_dp*(8._dp-(4*r0)/rz-  &
(4*rz)/r0)*z*(-1._dp+4*z**2))/(-8*r0*rz+8*z)+  &
(1.0931369993413877e-3_dp*(-1._dp+r0/2/rz+  &
rz/2/r0)*z*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)*(-1._dp+4*z**2))/(1+r0*rz - z)+  &
(5.465684996706938e-4_dp*(8._dp-(4*r0)/rz-  &
(4*rz)/r0)*z*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)*(-1._dp+4*z**2))/(-8*r0*rz+8*z)+  &
4.78247437211857e-4_dp*z*(1+(1.5_dp*r0)/rz+9*r0*rz-1.2e1_dp*z+  &
(1.5_dp*rz+3*z**1.5_dp)/r0)+  &
5.465684996706938e-4_dp*z*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)*(1+(1.5_dp*r0)/rz+9*r0*rz-1.2e1_dp*z+(1.5_dp*rz+  &
3*z**1.5_dp)/r0))/z**6+(l3*(2.391237186059285e-3_dp-  &
5.738969246542284e-3_dp*z**2+  &
4.372547997365551e-3_dp*z**2*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)+(1.0931369993413877e-3_dp*(-1._dp+r0/2/rz+  &
rz/2/r0)*z*(1-4*z**2)**2)/(1+r0*rz - z)+  &
(5.465684996706938e-4_dp*(8._dp-(4*r0)/rz-(4*rz)/r0)*z*(1-  &
4*z**2)**2)/(-8*r0*rz+8*z)-  &
2.7328424983534694e-3_dp*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)*(-1._dp+4*z**2)+5.465684996706938e-4_dp*z*(-1._dp+  &
4*z**2)*(1+(1.5_dp*r0)/rz+9*r0*rz-1.2e1_dp*z+(1.5_dp*rz+  &
3*z**1.5_dp)/r0)))/z**6+l2*((l3*(-5.465684996706939e-3_dp+  &
2.6235287984193305e-2_dp*z**2-1.749019198946221e-2_dp*z**4))/z**6+  &
(4.78247437211857e-3_dp-1.147793849308457e-2_dp*z**2+  &
8.745095994731102e-3_dp*z**2*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)+(2.1862739986827755e-3_dp*(-1._dp+r0/2/rz+  &
rz/2/r0)*z*(1-4*z**2)**2)/(1+r0*rz - z)+  &
(1.0931369993413877e-3_dp*(8._dp-(4*r0)/rz-(4*rz)/r0)*z*(1  &
-4*z**2)**2)/(-8*r0*rz+8*z)-  &
5.465684996706939e-3_dp*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)*(-1._dp+4*z**2)+1.0931369993413877e-3_dp*z*(-1._dp+  &
4*z**2)*(1+(1.5_dp*r0)/rz+9*r0*rz-1.2e1_dp*z+(1.5_dp*rz+  &
3*z**1.5_dp)/r0))/z**6))+  &
PL1*((r0*rz*(-7.884077256944878e-3_dp-1.4210780991438035e-2_dp*z+  &
5.7606933058852094e-2_dp*z**2+1.7766084967690188e-1_dp*z**3-  &
1.042824961242903e-1_dp*z**4-4.832709028445988e-1_dp*z**5)+  &
rz**2*(3.942038628472439e-3_dp+7.105390495719018e-3_dp*z-  &
2.8803466529426047e-2_dp*z**2-8.883042483845093e-2_dp*z**3+  &
5.214124806214516e-2_dp*z**4+2.416354514222994e-1_dp*z**5)+  &
r0**4*((-1.8381285106030592e-1_dp-3.6762570212061183e-1_dp*z)*z**3  &
+rz**2*(-9.83823299407249e-3_dp-2.623528798419331e-2_dp*z-  &
1.311764399209665e-2_dp*z**2))+r0**2*(3.942038628472439e-3_dp+  &
7.105390495719018e-3_dp*z-2.8803466529426047e-2_dp*z**2-  &
8.883042483845093e-2_dp*z**3+rz**2*(-1.8381285106030592e-1_dp-  &
3.6762570212061183e-1_dp*z)*z**3+5.214124806214516e-2_dp*z**4+  &
2.416354514222994e-1_dp*z**5+rz**4*(-9.83823299407249e-3_dp-  &
2.623528798419331e-2_dp*z-1.311764399209665e-2_dp*z**2))+  &
r0**3*(rz*(3.6762570212061183e-1_dp+  &
7.352514042412237e-1_dp*z)*z**3+rz**3*(1.9676465988144982e-2_dp+  &
5.247057596838662e-2_dp*z+  &
2.62352879841933e-2_dp*z**2)))/(d2*r0*rz*z**5)+  &
PL2*(((-6.996076795784882e-2_dp*r0**2+  &
1.3992153591569765e-1_dp*r0*rz-  &
6.996076795784882e-2_dp*rz**2)*z05m*z05**2)/(d2*z**5)+((1+  &
2*z)*(6.99607679578488e-2_dp*r0**2*rz**2*z-  &
1.7490191989462203e-1_dp*r0**2*rz**2*z05+  &
3.49803839789244e-2_dp*r0**2*z*z05+  &
3.49803839789244e-2_dp*rz**2*z*z05+  &
(5.24705759683866e-2_dp*(r0**2-2*r0*rz+rz**2)*z*(-2.5e-1_dp+  &
z**2))/(r0*rz - z05m)))/z**6)+  &
(rz**2*z*(-3.942038628472439e-3_dp-1.4989467752663899e-2_dp*z-  &
1.1754689759017507e-3_dp*z**2+8.647948688664743e-2_dp*z**3+  &
1.208177257111497e-1_dp*z**4)+  &
r0**2*(rz**2*(3.9420386284724396e-2_dp+1.199157420213112e-1_dp*z+  &
7.052813855410483e-3_dp*z**2-3.4591794754658984e-1_dp*z**3-  &
2.416354514222994e-1_dp*z**4)+z*(-3.942038628472439e-3_dp-  &
1.4989467752663899e-2_dp*z-1.1754689759017507e-3_dp*z**2+  &
8.647948688664743e-2_dp*z**3+  &
1.208177257111497e-1_dp*z**4)))/(r0*rz*z**6)+  &
PL4*((r0*rz*(-2.623528798419331e-2_dp+2.098823038735465e-1_dp*z**2  &
-4.19764607747093e-1_dp*z**4)+r0**2*(1.3117643992096655e-2_dp-  &
1.0494115193677325e-1_dp*z**2+2.098823038735465e-1_dp*z**4)+  &
rz**2*(1.3117643992096655e-2_dp-1.0494115193677325e-1_dp*z**2+  &
2.098823038735465e-1_dp*z**4))/(d2*r0*rz*z**5)+  &
(1.0494115193677322e-1_dp*rz**2*z05m*z*z05**2+  &
r0**2*(1.0494115193677322e-1_dp*z05m*z*z05**2+  &
rz**2*(1.3117643992096653e-1_dp+2.0988230387354645e-1_dp*z-  &
3.148234558103197e-1_dp*z**2-  &
4.19764607747093e-1_dp*z**3)))/(r0*rz*z**6))+  &
PL3*((r0*rz*(-1.3117643992096655e-2_dp+  &
1.0494115193677325e-1_dp*z**2-2.098823038735465e-1_dp*z**4)+  &
r0**2*(6.558821996048328e-3_dp-5.247057596838662e-2_dp*z**2+  &
1.0494115193677325e-1_dp*z**4)+rz**2*(6.558821996048328e-3_dp-  &
5.247057596838662e-2_dp*z**2+  &
1.0494115193677325e-1_dp*z**4))/(d2*r0*rz*z**5)+  &
(5.24705759683866e-2_dp*rz**2*z05m*z*z05**2+  &
r0**2*(5.24705759683866e-2_dp*z05m*z*z05**2+  &
rz**2*(6.5588219960483265e-2_dp+1.0494115193677322e-1_dp*z-  &
1.5741172790515985e-1_dp*z**2-  &
2.098823038735465e-1_dp*z**3)))/(r0*rz*z**6)))+  &
l1*((1.399215359156976e-1_dp*l2**2*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)+  &
(3.49803839789244e-2_dp*l3**2*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)+(1.0476837861865178e2_dp*(r0**2-  &
2*r0*rz+rz**2)*(r0*rz - z05m)*(-1._dp+  &
z)*(-3.7913315126399993e-1_dp+z)*(1+z)*(rz+r0*z05m-  &
z**1.5_dp))/(d1**3*r0*z**2)+  &
((-3.702260139950978_dp*d2*(-2._dp+r0/rz+rz/r0)*(1-  &
z)**2*(1+z))/z-(5.370932308906429_dp*d2**2*(r0**2-  &
2*r0*rz+rz**2)*(-1._dp+z)*(1.0647869616635923_dp+  &
1.2416980235463153_dp*z+z**2))/(r0*rz*z**2))/d1**2+  &
((2.6854661544532146_dp*(r0**2-2*r0*rz+rz**2)*(-1._dp+  &
z)*(1.0647869616635923_dp+1.2416980235463153_dp*z+  &
z**2))/(r0*rz*z**2)+(d2*(5.7188986945012985_dp+  &
4.750886690344309e-1_dp*z+2.6854661544532137_dp*z**3))/z**3)/d1+  &
PL2**2*((5.596861436627904e-1_dp*(r0**2-2*r0*rz+  &
rz**2)*(2.5e-1_dp - z**2)**2)/(d2*r0*rz*z**5)+((1+  &
2*z)*(-2.798430718313952e-1_dp*r0*rz*z05m*z-  &
1.399215359156976e-1_dp*r0*rz*z*z05-  &
(5.24705759683866e-2_dp*(2-r0/rz-rz/r0)*(1-  &
2*z)*z*(1-4*z**2))/(-1._dp-2*r0*rz+2*z)+  &
r0*rz*(-1.7490191989462203e-1_dp+6.996076795784881e-1_dp*z**2)+  &
(1.74901919894622e-2_dp*r0*z-6.99607679578488e-2_dp*r0*z**3)/rz+  &
(1.74901919894622e-2_dp*rz*z-  &
6.99607679578488e-2_dp*rz*z**3)/r0))/z**6)+  &
(rz**2*(3.449283799913384e-3_dp+1.015925531222658e-2_dp*z-  &
2.795273928870987e-2_dp*z**2-1.085254099884784e-1_dp*z**3+  &
1.7803009237642873e-2_dp*z**4+3.624531771334491e-1_dp*z**5)+  &
r0**4*((4.135789148856884e-1_dp+2.7571927659045894e-1_dp*z)*z**3+  &
rz**2*(1.1068012118331547e-2_dp+1.4757349491108736e-2_dp*z+  &
4.919116497036245e-3_dp*z**2))+r0*rz*(-6.898567599826768e-3_dp-  &
2.031851062445316e-2_dp*z+5.590547857741974e-2_dp*z**2+  &
2.170508199769568e-1_dp*z**3-3.5606018475285746e-2_dp*z**4-  &
7.249063542668982e-1_dp*z**5+rz*z**1.5_dp*(2.3652231770834637e-2_dp  &
+4.263234297431412e-2_dp*z-7.8211872093217725e-2_dp*z**2-  &
3.624531771334491e-1_dp*z**3)+rz**2*(1.1826115885417319e-2_dp+  &
2.131617148715706e-2_dp*z-3.9105936046608862e-2_dp*z**2-  &
1.8122658856672456e-1_dp*z**3))+  &
r0**3*(rz**3*(-2.2136024236663094e-2_dp-2.951469898221747e-2_dp*z-  &
9.83823299407249e-3_dp*z**2)+z**1.5_dp*(2.3652231770834637e-2_dp+  &
4.263234297431412e-2_dp*z-7.8211872093217725e-2_dp*z**2-  &
3.624531771334491e-1_dp*z**3)+rz*(1.1826115885417319e-2_dp+  &
2.131617148715706e-2_dp*z-3.9105936046608862e-2_dp*z**2-  &
1.0083844183381014_dp*z**3-5.514385531809179e-1_dp*z**4))+  &
r0**2*(3.449283799913384e-3_dp+1.015925531222658e-2_dp*z-  &
2.795273928870987e-2_dp*z**2-1.085254099884784e-1_dp*z**3+  &
1.7803009237642873e-2_dp*z**4+3.624531771334491e-1_dp*z**5+  &
rz**4*(1.1068012118331547e-2_dp+1.4757349491108736e-2_dp*z+  &
4.919116497036245e-3_dp*z**2)+  &
rz*z**1.5_dp*(-4.7304463541669275e-2_dp-8.526468594862823e-2_dp*z+  &
1.5642374418643545e-1_dp*z**2+7.249063542668982e-1_dp*z**3)+  &
rz**2*(-2.3652231770834637e-2_dp-4.263234297431412e-2_dp*z+  &
7.8211872093217725e-2_dp*z**2+7.760320920191375e-1_dp*z**3+  &
2.7571927659045894e-1_dp*z**4)))/(d2*r0*rz*z**5)+  &
(r0*rz**3*z*(4.927548285590554e-4_dp+8.177471943766133e-2_dp*z-  &
3.8101238150423375e-3_dp*z**2-5.417438213948252e-1_dp*z**3+  &
5.1952301496464415e-1_dp*z**4-2.7183988285008684e-1_dp*z**5)+  &
rz**2*z**2*(-4.434793457031494e-3_dp-8.937286476193938e-2_dp*z-  &
4.620460012183861e-2_dp*z**2+4.917374568683684e-1_dp*z**3-  &
3.064185513858778e-1_dp*z**4-4.530664714168115e-2_dp*z**5)+  &
r0**4*rz**2*(rz**2*(-4.434793457031494e-2_dp-  &
8.760074623230583e-2_dp*z+5.6014098874134355e-2_dp*z**2+  &
3.109458188966958e-1_dp*z**3+9.06132942833623e-2_dp*z**4)+  &
z*(-7.3913224283858225e-3_dp-3.401830997895347e-2_dp*z-  &
3.4429234983312496e-1_dp*z**2-1.169938770372288e-1_dp*z**3+  &
3.171465299917682e-1_dp*z**4))+  &
r0**3*rz*(rz**2*(-4.434793457031494e-2_dp+  &
2.4747354679158713e-2_dp*z+3.2115239785472935e-1_dp*z**2+  &
2.0597043500383756e-1_dp*z**3-1.0501552648299135_dp*z**4-  &
9.06132942833623e-1_dp*z**5)+z*(4.927548285590554e-4_dp-  &
1.0828455022713415e-3_dp*z-3.188731951814062e-1_dp*z**2+  &
3.9510544192802084e-1_dp*z**3+5.689372933038687e-1_dp*z**4-  &
2.7183988285008676e-1_dp*z**5))+  &
r0**2*z*(rz**4*(-7.3913224283858225e-3_dp+  &
4.883925496097918e-2_dp*z+1.364858514131042e-1_dp*z**2-  &
9.228673786761659e-2_dp*z**3+3.17146529991768e-1_dp*z**4)+  &
z*(-4.434793457031494e-3_dp-6.5152998220067335e-3_dp*z+  &
3.517160361844579e-1_dp*z**2+3.5666394791751443e-2_dp*z**3-  &
3.311256905554902e-1_dp*z**4-4.530664714168098e-2_dp*z**5)+  &
rz**2*(5.223201182725983e-2_dp+4.957951539648404e-2_dp*z-  &
2.3120071367062014e-1_dp*z**2-4.349435076512668e-1_dp*z**3+  &
4.975739945109186e-1_dp*z**4+  &
8.155196485502605e-1_dp*z**5)))/(r0*rz*(r0*rz - z)*(1+r0*rz-  &
z)*z**6)+l5*((1.399215359156976e-1_dp*PL1*(r0**2-  &
2*r0*rz+rz**2)*z05m*z05**2)/(d2*z**5)+  &
(2.798430718313952e-1_dp*PL2*(r0**2-2*r0*rz+rz**2)*(-0.5_dp+  &
z)*z05**2)/(d2*z**5)+(1.311764399209665e-2_dp*PL3*(-2._dp+  &
r0/rz+rz/r0)*(1-2*z)*(1+2*z)*(1-4*z**2))/(d2*z**5)+  &
(2.62352879841933e-2_dp*PL4*(-2._dp+r0/rz+rz/r0)*(1-2*z)*(1+  &
2*z)*(1-4*z**2))/(d2*z**5)+  &
(r0*rz*(-1.5768154513889756e-2_dp-2.842156198287608e-2_dp*z+  &
1.152138661177042e-1_dp*z**2+3.5532169935380375e-1_dp*z**3-  &
2.085649922485806e-1_dp*z**4-9.665418056891976e-1_dp*z**5)+  &
rz**2*(7.884077256944878e-3_dp+1.421078099143804e-2_dp*z-  &
5.7606933058852094e-2_dp*z**2-1.7766084967690188e-1_dp*z**3+  &
1.042824961242903e-1_dp*z**4+4.832709028445988e-1_dp*z**5)+  &
r0**3*(rz*(-3.6762570212061183e-1_dp-  &
7.352514042412237e-1_dp*z)*z**3+rz**3*(-1.9676465988144982e-2_dp-  &
5.247057596838662e-2_dp*z-2.62352879841933e-2_dp*z**2))+  &
r0**4*((1.8381285106030592e-1_dp+3.6762570212061183e-1_dp*z)*z**3+  &
rz**2*(9.83823299407249e-3_dp+2.623528798419331e-2_dp*z+  &
1.311764399209665e-2_dp*z**2))+r0**2*(7.884077256944878e-3_dp+  &
1.421078099143804e-2_dp*z-5.7606933058852094e-2_dp*z**2-  &
1.7766084967690188e-1_dp*z**3+rz**2*(1.8381285106030592e-1_dp+  &
3.6762570212061183e-1_dp*z)*z**3+1.042824961242903e-1_dp*z**4+  &
4.832709028445988e-1_dp*z**5+rz**4*(9.83823299407249e-3_dp+  &
2.623528798419331e-2_dp*z+  &
1.311764399209665e-2_dp*z**2)))/(d2*r0*rz*z**5))+  &
l4*((1.399215359156976e-1_dp*PL1*(r0**2-2*r0*rz+  &
rz**2)*z05m*z05**2)/((r0*rz - z05m)*z**5)+  &
(2.798430718313952e-1_dp*PL2*(r0**2-2*r0*rz+rz**2)*(-0.5_dp+  &
z)*z05**2)/((r0*rz - z05m)*z**5)+  &
(2.62352879841933e-2_dp*PL3*(2-r0/rz-rz/r0)*(1-  &
2*z)*(1+2*z)*(1-4*z**2))/(z**5*(-1._dp-2*r0*rz+  &
2*z))+(5.24705759683866e-2_dp*PL4*(2-r0/rz-  &
rz/r0)*(1-2*z)*(1+2*z)*(1-  &
4*z**2))/(z**5*(-1._dp-2*r0*rz+2*z))+  &
(r0*rz*(-1.5768154513889756e-2_dp-2.842156198287608e-2_dp*z+  &
1.152138661177042e-1_dp*z**2+3.5532169935380375e-1_dp*z**3-  &
2.085649922485806e-1_dp*z**4-9.665418056891976e-1_dp*z**5)+  &
rz**2*(7.884077256944878e-3_dp+1.421078099143804e-2_dp*z-  &
5.7606933058852094e-2_dp*z**2-1.7766084967690188e-1_dp*z**3+  &
1.042824961242903e-1_dp*z**4+4.832709028445988e-1_dp*z**5)+  &
r0**3*(rz*(-3.6762570212061183e-1_dp-  &
7.352514042412237e-1_dp*z)*z**3+rz**3*(-1.9676465988144982e-2_dp-  &
5.247057596838662e-2_dp*z-2.62352879841933e-2_dp*z**2))+  &
r0**4*((1.8381285106030592e-1_dp+3.6762570212061183e-1_dp*z)*z**3+  &
rz**2*(9.83823299407249e-3_dp+2.623528798419331e-2_dp*z+  &
1.311764399209665e-2_dp*z**2))+r0**2*(7.884077256944878e-3_dp+  &
1.421078099143804e-2_dp*z-5.7606933058852094e-2_dp*z**2-  &
1.7766084967690188e-1_dp*z**3+rz**2*(1.8381285106030592e-1_dp+  &
3.6762570212061183e-1_dp*z)*z**3+1.042824961242903e-1_dp*z**4+  &
4.832709028445988e-1_dp*z**5+rz**4*(9.83823299407249e-3_dp+  &
2.623528798419331e-2_dp*z+  &
1.311764399209665e-2_dp*z**2)))/(r0*rz*(0.5_dp+r0*rz-  &
z)*z**5))+PL1**2*((r0**2*(-4.372547997365551e-3_dp+  &
3.49803839789244e-2_dp*z**2-6.99607679578488e-2_dp*z**4)+  &
rz**2*(-4.372547997365551e-3_dp+3.49803839789244e-2_dp*z**2-  &
6.99607679578488e-2_dp*z**4)+r0*rz*(8.745095994731102e-3_dp-  &
6.99607679578488e-2_dp*z**2+  &
1.399215359156976e-1_dp*z**4))/(d2*r0*rz*z**5)+  &
(rz**2*z*(4.372547997365551e-3_dp+8.745095994731102e-3_dp*z-  &
1.74901919894622e-2_dp*z**2-3.49803839789244e-2_dp*z**3)+  &
r0**2*(-3.498038397892442e-2_dp*z05m*z*z05**2+  &
rz**2*(-4.372547997365551e-2_dp-6.99607679578488e-2_dp*z+  &
1.0494115193677322e-1_dp*z**2+  &
1.3992153591569765e-1_dp*z**3)))/(r0*rz*z**6))+PL3*(((-2._dp+  &
r0/rz+rz/r0)*(1-4*z**2)*(5.738969246542284e-3_dp+  &
1.9676465988144975e-2_dp*r0*rz+6.558821996048325e-3_dp*z+  &
3.935293197628995e-2_dp*r0*z**1.5_dp-  &
3.935293197628995e-2_dp*z**2))/(d2*z**5)+  &
(rz**2*z**2*(-7.378674745554367e-3_dp-1.3540008004671134e-1_dp*z+  &
1.7229345377448315e-1_dp*z**2+5.416003201868453e-1_dp*z**3-  &
5.711150191690628e-1_dp*z**4)+r0*rz**3*z*(8.198527495060407e-4_dp+  &
1.3458022729720531e-1_dp*z-2.494839886062656e-1_dp*z**2-  &
5.3832090918882125e-1_dp*z**3+9.848183104329657e-1_dp*z**4)+  &
r0**4*rz**2*(rz**2*(-7.378674745554368e-2_dp-  &
3.9352931976289964e-2_dp*z+1.7708819389330483e-1_dp*z**2+  &
7.870586395257993e-2_dp*z**3)+z*(-1.229779124259061e-2_dp-  &
3.4433815479253713e-2_dp*z-5.022473882105555e-1_dp*z**2+  &
1.3773526191701486e-1_dp*z**3))+  &
r0**3*rz*(rz**2*(-7.378674745554368e-2_dp+  &
1.475734949110873e-1_dp*z+3.344999217984646e-1_dp*z**2-  &
4.328822517391895e-1_dp*z**3-4.722351837154795e-1_dp*z**4)+  &
z*(8.198527495060407e-4_dp-3.279410998024162e-3_dp*z-  &
5.252032651967245e-1_dp*z**2+1.1159947503539325_dp*z**3-  &
1.1805879592886988e-1_dp*z**4))+  &
r0**2*z*(rz**4*(-1.229779124259061e-2_dp+  &
1.0342582281597574e-1_dp*z+4.919116497036244e-2_dp*z**2-  &
4.13703291263903e-1_dp*z**3)+z*(-7.378674745554367e-3_dp+  &
2.459558248518121e-3_dp*z+5.858723686601715e-1_dp*z**2-  &
5.612767861749903e-1_dp*z**3-1.9676465988144982e-2_dp*z**4)+  &
rz**2*(8.690439144764033e-2_dp-4.755145947135036e-2_dp*z-  &
3.4761756579056122e-1_dp*z**2+1.508529059091115e-1_dp*z**3+  &
3.9352931976289964e-1_dp*z**4)))/(r0*rz*(r0*rz - z)*(1+r0*rz  &
 - z)*z**6))+PL4*(((-2._dp+r0/rz+rz/r0)*(1-  &
4*z**2)*(1.147793849308457e-2_dp+3.935293197628995e-2_dp*r0*rz+  &
1.311764399209665e-2_dp*z+7.87058639525799e-2_dp*r0*z**1.5_dp-  &
7.87058639525799e-2_dp*z**2))/(d2*z**5)+  &
(rz**2*z**2*(-1.4757349491108736e-2_dp-2.7080016009342267e-1_dp*z+  &
3.445869075489663e-1_dp*z**2+1.0832006403736907_dp*z**3-  &
1.1422300383381256_dp*z**4)+r0*rz**3*z*(1.6397054990120816e-3_dp  &
+2.6916045459441063e-1_dp*z-4.989679772125312e-1_dp*z**2-  &
1.0766418183776425_dp*z**3+1.9696366208659315_dp*z**4)+  &
r0**4*rz**2*(rz**2*(-1.4757349491108736e-1_dp-  &
7.870586395257993e-2_dp*z+3.5417638778660967e-1_dp*z**2+  &
1.5741172790515985e-1_dp*z**3)+z*(-2.459558248518122e-2_dp-  &
6.886763095850743e-2_dp*z-1.004494776421111_dp*z**2+  &
2.7547052383402972e-1_dp*z**3))+  &
r0**3*rz*(rz**2*(-1.4757349491108736e-1_dp+  &
2.951469898221746e-1_dp*z+6.689998435969292e-1_dp*z**2-  &
8.65764503478379e-1_dp*z**3-9.44470367430959e-1_dp*z**4)+  &
z*(1.6397054990120816e-3_dp-6.558821996048324e-3_dp*z-  &
1.0504065303934491_dp*z**2+2.231989500707865_dp*z**3-  &
2.3611759185773975e-1_dp*z**4))+  &
r0**2*z*(rz**4*(-2.459558248518122e-2_dp+2.068516456319515e-1_dp*z  &
+9.838232994072488e-2_dp*z**2-8.27406582527806e-1_dp*z**3)+  &
z*(-1.4757349491108736e-2_dp+4.919116497036242e-3_dp*z+  &
1.1717447373203431_dp*z**2-1.1225535723499807_dp*z**3-  &
3.9352931976289964e-2_dp*z**4)+rz**2*(1.7380878289528068e-1_dp-  &
9.510291894270072e-2_dp*z-6.9523513158112245e-1_dp*z**2+  &
3.01705811818223e-1_dp*z**3+  &
7.870586395257993e-1_dp*z**4)))/(r0*rz*(r0*rz - z)*(1+r0*rz-  &
z)*z**6))+PL1*((PL3*(6.5588219960483265e-2_dp-  &
3.148234558103197e-1_dp*z**2+2.098823038735465e-1_dp*z**4))/z**6+  &
(PL4*(1.3117643992096653e-1_dp-6.296469116206394e-1_dp*z**2+  &
4.19764607747093e-1_dp*z**4))/z**6+  &
(r0**3*rz**2*(rz*(-1.311764399209665e-2_dp-  &
2.62352879841933e-2_dp*z)+(-2.62352879841933e-2_dp-  &
5.24705759683866e-2_dp*z)*z**1.5_dp)+  &
r0*rz*(rz**4*(-1.311764399209665e-2_dp-2.62352879841933e-2_dp*z)+  &
rz**3*(-2.62352879841933e-2_dp-5.24705759683866e-2_dp*z)*z**1.5_dp+  &
1.8381285106030587e-1_dp*z**3-7.352514042412235e-1_dp*z**5+  &
rz**2*(1.749019198946221e-2_dp+3.060783598155885e-2_dp*z-  &
7.433331595521439e-2_dp*z**2-1.311764399209665e-1_dp*z**3))+  &
rz**2*(-9.190642553015294e-2_dp*z**3+3.6762570212061174e-1_dp*z**5  &
+rz**2*(-8.745095994731106e-3_dp-1.5303917990779425e-2_dp*z+  &
3.7166657977607196e-2_dp*z**2+6.558821996048325e-2_dp*z**3))+  &
r0**2*(rz**4*(2.62352879841933e-2_dp+5.24705759683866e-2_dp*z)+  &
rz**3*(5.24705759683866e-2_dp+1.049411519367732e-1_dp*z)*z**1.5_dp-  &
9.190642553015294e-2_dp*z**3+3.6762570212061174e-1_dp*z**5+  &
rz**2*(-8.745095994731106e-3_dp-1.5303917990779425e-2_dp*z+  &
3.7166657977607196e-2_dp*z**2+  &
6.558821996048325e-2_dp*z**3)))/(d2*rz**2*z**5)+PL2*(((1+  &
2*z)*(-2.798430718313952e-1_dp*r0*rz*z05m*z-  &
1.399215359156976e-1_dp*r0*rz*z*z05-  &
(2.62352879841933e-2_dp*(2-r0/rz-rz/r0)*(1-  &
2*z)*z*(1-4*z**2))/(-1._dp-2*r0*rz+2*z)+  &
r0*rz*(-1.7490191989462203e-1_dp+6.996076795784881e-1_dp*z**2)+  &
(1.74901919894622e-2_dp*r0*z-6.99607679578488e-2_dp*r0*z**3)/rz+  &
(1.74901919894622e-2_dp*rz*z-  &
6.99607679578488e-2_dp*rz*z**3)/r0))/z**6+  &
(r0*rz*(-1.74901919894622e-2_dp+1.399215359156976e-1_dp*z**2-  &
2.798430718313952e-1_dp*z**4)+r0**2*(8.745095994731102e-3_dp-  &
6.99607679578488e-2_dp*z**2+1.399215359156976e-1_dp*z**4)+  &
rz**2*(8.745095994731102e-3_dp-6.99607679578488e-2_dp*z**2+  &
1.399215359156976e-1_dp*z**4))/(d2*r0*rz*z**5))+  &
(rz**2*z*(-3.9420386284724396e-2_dp+(-1.7422737681027756e-2_dp+  &
9.83823299407249e-3_dp*rz**2)*z+(2.296639231423085e-1_dp+  &
1.0830348052027379e-1_dp*rz**2)*z**2+(1.8250090017724732e-1_dp+  &
7.87887815380563e-2_dp*rz**2)*z**3+(-3.676977699479411e-1_dp-  &
1.9693049505240257e-1_dp*rz**2)*z**4+1.9618892165444357e-1_dp*z**5  &
-1.838128510603061e-1_dp*z**6)+r0*rz**3*(3.9420386284724396e-2_dp+  &
(-2.1997648603696645e-2_dp-5.46568499670694e-3_dp*rz**2)*z+  &
(-2.8650704710806063e-1_dp-9.409269952883573e-2_dp*rz**2)*z**2+  &
(-9.680101000691055e-3_dp+4.372547997365555e-3_dp*rz**2)*z**3+  &
(7.230194693017449e-1_dp+3.4139041413641857e-1_dp*rz**2)*z**4-  &
2.0856499224858105e-1_dp*z**5+3.676257021206122e-1_dp*z**6)+  &
r0**4*rz**2*((-9.190642553015296e-2_dp-  &
1.8381285106030592e-1_dp*z)*z**4+rz**4*(4.919116497036246e-2_dp+  &
1.0494115193677322e-1_dp*z+3.935293197628995e-2_dp*z**2)+  &
rz**2*z*(3.2794109980241597e-3_dp+2.623528798419329e-2_dp*z+  &
4.069786340969018e-1_dp*z**2+3.6762570212061183e-1_dp*z**3))+  &
r0**3*rz*(-9.190642553015296e-2_dp*z**4+  &
3.6762570212061183e-1_dp*z**6+rz**4*(4.919116497036246e-2_dp-  &
1.9676465988144975e-2_dp*z-2.7547052383402963e-1_dp*z**2-  &
1.836470158893531e-1_dp*z**3)+rz**2*z*(-5.46568499670694e-3_dp-  &
2.1862739986827755e-3_dp*z+3.7199825011797745e-1_dp*z**2-  &
3.9386099010480526e-1_dp*z**3-7.352514042412237e-1_dp*z**4))+  &
r0**2*(z**5*(9.190642553015296e-2_dp+9.190642553015296e-2_dp*z-  &
1.8381285106030592e-1_dp*z**2)+rz**6*z*(3.2794109980241597e-3_dp-  &
6.567113754595965e-2_dp*z-1.4445991908401594e-1_dp*z**2)+  &
rz**2*z**2*(9.83823299407249e-3_dp+1.6397054990120816e-2_dp*z-  &
3.807433461127085e-1_dp*z**2-1.3117643992096653e-2_dp*z**3+  &
3.6762570212061183e-1_dp*z**4)+rz**4*(3.9420386284724396e-2_dp-  &
1.0931369993413889e-3_dp*z-2.373158821376982e-1_dp*z**2-  &
2.1977271143547163e-1_dp*z**3+1.566701545072005e-1_dp*z**4-  &
1.838128510603061e-1_dp*z**5)))/(rz**2*(r0*rz - z)*(1+r0*rz-  &
z)*z**6))+PL2*((PL3*(1.3117643992096653e-1_dp-  &
6.296469116206394e-1_dp*z**2+4.19764607747093e-1_dp*z**4))/z**6+  &
(PL4*(2.6235287984193305e-1_dp-1.2592938232412787_dp*z**2+  &
8.39529215494186e-1_dp*z**4))/z**6-(5.24705759683866e-2_dp*(r0**2-  &
2*r0*rz+rz**2)*z05*(r0*rz**3+(2.802528429206536e1_dp-  &
5.605056858413072e1_dp*z)*z**3+rz**2*(1.791666666666667_dp-  &
1.666666666666667_dp*z+2*r0*z**1.5_dp-  &
4*z**2)))/(d2*rz**2*z**5)+  &
(rz**2*z*(-3.9420386284724396e-2_dp+(6.141803488842104e-2_dp-  &
4.919116497036243e-3_dp*rz**2)*z+(2.64509398504364e-1_dp+  &
9.354613102916504e-2_dp*rz**2)*z**2+(-2.768269461073695e-1_dp-  &
6.895054854398383e-2_dp*rz**2)*z**3+(-1.0084188468928947_dp-  &
3.741845241166602e-1_dp*rz**2)*z**4+(1.207303738140785_dp+  &
3.5450805812851516e-1_dp*rz**2)*z**5+5.266864119926429e-1_dp*z**6-  &
7.352514042412235e-1_dp*z**7)+  &
r0**5*rz**3*((-1.8381285106030592e-1_dp-  &
3.6762570212061183e-1_dp*z)*z**4+rz**4*(9.838232994072492e-2_dp+  &
2.0988230387354645e-1_dp*z+7.87058639525799e-2_dp*z**2)+  &
rz**2*z*(6.558821996048319e-3_dp+5.247057596838658e-2_dp*z+  &
8.139572681938036e-1_dp*z**2+7.352514042412237e-1_dp*z**3))+  &
r0**4*rz**2*(rz**4*(1.4757349491108736e-1_dp-  &
3.2794109980241624e-2_dp*z-7.214704195653159e-1_dp*z**2-  &
4.459998957312862e-1_dp*z**3)+rz**2*z*(7.105390495719011e-3_dp+  &
2.5142150984851916e-2_dp*z+1.0394751604000352_dp*z**2-  &
1.273406478259092_dp*z**3-2.2057542127236713_dp*z**4))+  &
r0**2*(z**5*(-1.8381285106030596e-1_dp+1.8381285106030596e-1_dp*z+  &
7.352514042412238e-1_dp*z**2-7.352514042412238e-1_dp*z**3)+  &
rz**6*z*(7.105390495719011e-3_dp-2.5057712560560703e-1_dp*z-  &
6.340194596180043e-2_dp*z**2+9.323477344645791e-1_dp*z**3)+  &
rz**2*z**2*(-4.919116497036243e-3_dp+1.6397054990120783e-3_dp*z-  &
8.961083783153605e-1_dp*z**2+1.8315696886070105_dp*z**3+  &
2.560262270852186_dp*z**4-5.146759829688566_dp*z**5)+  &
rz**4*(1.1826115885417319e-1_dp-1.5344390575840092e-1_dp*z-  &
7.687907703788467e-1_dp*z**2+3.929105787437024e-1_dp*z**3+  &
2.160562064666933_dp*z**4-1.0716948724770288_dp*z**5))+  &
r0**3*rz*(rz**6*z*(6.558821996048319e-3_dp-  &
1.3134227509191931e-1_dp*z-2.889198381680319e-1_dp*z**2)+  &
z**4*(1.8381285106030596e-1_dp-1.8381285106030596e-1_dp*z-  &
5.51438553180918e-1_dp*z**2+1.102877106361836_dp*z**3)+  &
rz**4*(1.2803193753981126e-1_dp-1.497597689097701e-1_dp*z-  &
7.3042582212128115e-1_dp*z**2+4.5807404836632815e-2_dp*z**3+  &
7.593402047456873e-1_dp*z**4-3.676257021206122e-1_dp*z**5)+  &
rz**2*z*(9.291664494401793e-3_dp+8.745095994731118e-3_dp*z+  &
3.304590441430047e-1_dp*z**2-2.372076871794515_dp*z**3+  &
1.0494115193677314e-1_dp*z**4+4.411508425447343_dp*z**5))+  &
r0*rz*(z**5*(5.514385531809179e-1_dp-5.514385531809179e-1_dp*z-  &
2.2057542127236713_dp*z**2+2.2057542127236713_dp*z**3)+  &
rz**4*z*(9.291664494401793e-3_dp-8.316132953542185e-2_dp*z+  &
3.304590441430047e-1_dp*z**2+5.6892874517038e-1_dp*z**3-  &
9.979359544250627e-1_dp*z**4)+rz**2*(3.9420386284724396e-2_dp-  &
1.7967919374259422e-1_dp*z-2.4784252628050547e-1_dp*z**2+  &
1.0128236065059748_dp*z**3+1.2453654863310664_dp*z**4-  &
2.902085091952581_dp*z**5+1.5296228751740504e-1_dp*z**6+  &
1.1028771063618361_dp*z**7)))/(rz**2*(r0*rz - z)*(0.5_dp+  &
r0*rz - z)*(1+r0*rz - z)*z**6))+  &
l2*((1.399215359156976e-1_dp*l3*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)-  &
(6.99607679578488e-2_dp*l5*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)-  &
(6.99607679578488e-2_dp*l4*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/((r0*rz - z05m)*z**5)+  &
PL2*(((-6.996076795784881e-1_dp*r0**2+  &
1.3992153591569758_dp*r0*rz-  &
6.996076795784881e-1_dp*rz**2)*z05m*z05**2)/(d2*z**5)+((1+  &
2*z)*(1.399215359156976e-1_dp*r0**2*rz**2*z-  &
3.4980383978924405e-1_dp*r0**2*rz**2*z05+  &
6.99607679578488e-2_dp*r0**2*z*z05+  &
6.99607679578488e-2_dp*rz**2*z*z05+  &
(1.049411519367732e-1_dp*(r0**2-2*r0*rz+rz**2)*z*(-2.5e-1_dp+  &
z**2))/(r0*rz - z05m)))/z**6)-  &
(5.24705759683866e-2_dp*(r0**2-2*r0*rz+  &
rz**2)*(1.5025711289494927e-1_dp+2.708333333333333e-1_dp*z-  &
1.0978902364929275_dp*z**2-3.3859138459608698_dp*z**3+  &
1.9874471396525217_dp*z**4+9.210322050510145_dp*z**5+  &
r0**2*(z**3*(1.401264214603268e1_dp+2.802528429206536e1_dp*z)+  &
rz**2*(7.5e-1_dp+2 * z+  &
z**2))))/(d2*r0*rz*z**5)+(rz**2*z*(-7.884077256944878e-3_dp-  &
2.9978935505327797e-2_dp*z-2.3509379518035014e-3_dp*z**2+  &
1.7295897377329486e-1_dp*z**3+2.416354514222994e-1_dp*z**4)+  &
r0**2*(rz**2*(7.884077256944879e-2_dp+2.398314840426224e-1_dp*z+  &
1.4105627710820967e-2_dp*z**2-6.918358950931797e-1_dp*z**3-  &
4.832709028445988e-1_dp*z**4)+z*(-7.884077256944878e-3_dp-  &
2.9978935505327797e-2_dp*z-2.3509379518035014e-3_dp*z**2+  &
1.7295897377329486e-1_dp*z**3+  &
2.416354514222994e-1_dp*z**4)))/(r0*rz*z**6)+  &
PL1*(((-1.3992153591569765e-1_dp*r0**2+  &
2.7984307183139516e-1_dp*r0*rz- 1.3992153591569765e-1_dp*rz**2)*z05m*z05**2)&
/(d2*z**5)+(6.99607679578488e-2_dp*rz**2*z*z05**2+  &
r0**2*(6.99607679578488e-2_dp*z*z05**2+  &
rz**2*(-8.745095994731102e-2_dp-2.798430718313952e-1_dp*z-  &
2.0988230387354645e-1_dp*z**2)))/z**6)+  &
PL4*((-4.197646077470928e-1_dp*(r0**2-2*r0*rz+rz**2)*(2.5e-1_dp  &
 - z**2)**2)/(d2*r0*rz*z**5)+  &
(2.098823038735464e-1_dp*rz**2*z05m*z*z05**2+  &
r0**2*(2.098823038735464e-1_dp*z05m*z*z05**2+  &
rz**2*(2.6235287984193305e-1_dp+4.197646077470929e-1_dp*z-  &
6.296469116206394e-1_dp*z**2-  &
8.39529215494186e-1_dp*z**3)))/(r0*rz*z**6))+  &
PL3*((-2.098823038735464e-1_dp*(r0**2-2*r0*rz+rz**2)*(2.5e-1_dp  &
 - z**2)**2)/(d2*r0*rz*z**5)+  &
(1.0494115193677322e-1_dp*rz**2*z05m*z*z05**2+  &
r0**2*(1.0494115193677322e-1_dp*z05m*z*z05**2+  &
rz**2*(1.3117643992096653e-1_dp+2.0988230387354645e-1_dp*z-  &
3.148234558103197e-1_dp*z**2-  &
4.19764607747093e-1_dp*z**3)))/(r0*rz*z**6)))+  &
l3*((-3.49803839789244e-2_dp*l5*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/(d2*z**5)-  &
(3.49803839789244e-2_dp*l4*r0*rz*(r0**2-2*r0*rz+  &
rz**2)*z05**2)/((r0*rz - z05m)*z**5)-  &
(2.62352879841933e-2_dp*(r0**2-2*r0*rz+  &
rz**2)*(1.5025711289494927e-1_dp+2.708333333333333e-1_dp*z-  &
1.0978902364929275_dp*z**2-3.3859138459608698_dp*z**3+  &
1.9874471396525217_dp*z**4+9.210322050510145_dp*z**5+  &
r0**2*(z**3*(1.401264214603268e1_dp+2.802528429206536e1_dp*z)+  &
rz**2*(7.5e-1_dp+2 * z+  &
z**2))))/(d2*r0*rz*z**5)+PL2*(((1+  &
2*z)*(6.99607679578488e-2_dp*r0**2*rz**2*z-  &
1.7490191989462203e-1_dp*r0**2*rz**2*z05+  &
3.49803839789244e-2_dp*r0**2*z*z05+  &
3.49803839789244e-2_dp*rz**2*z*z05+  &
(5.24705759683866e-2_dp*(r0**2-2*r0*rz+rz**2)*z*(-2.5e-1_dp+  &
z**2))/(r0*rz - z05m)))/z**6+  &
(6.996076795784881e-1_dp*r0*rz*z05m*z05**2+  &
r0**2*(4.37254799736555e-2_dp+8.7450959947311e-2_dp*z-  &
1.74901919894622e-1_dp*z**2-3.49803839789244e-1_dp*z**3)+  &
rz**2*(4.37254799736555e-2_dp+8.7450959947311e-2_dp*z-  &
1.74901919894622e-1_dp*z**2-  &
3.49803839789244e-1_dp*z**3))/(d2*z**5))+  &
(rz**2*z*(-3.942038628472439e-3_dp-1.4989467752663899e-2_dp*z-  &
1.1754689759017507e-3_dp*z**2+8.647948688664743e-2_dp*z**3+  &
1.208177257111497e-1_dp*z**4)+  &
r0**2*(rz**2*(3.9420386284724396e-2_dp+1.199157420213112e-1_dp*z+  &
7.052813855410483e-3_dp*z**2-3.4591794754658984e-1_dp*z**3-  &
2.416354514222994e-1_dp*z**4)+z*(-3.942038628472439e-3_dp-  &
1.4989467752663899e-2_dp*z-1.1754689759017507e-3_dp*z**2+  &
8.647948688664743e-2_dp*z**3+  &
1.208177257111497e-1_dp*z**4)))/(r0*rz*z**6)+  &
PL1*(((-6.996076795784882e-2_dp*r0**2+  &
1.3992153591569765e-1_dp*r0*rz-  6.996076795784882e-2_dp*rz**2)*z05m*&
z05**2)/(d2*z**5)+(3.49803839789244e-2_dp*rz**2*z*z05**2+  &
r0**2*(3.49803839789244e-2_dp*z*z05**2+  &
rz**2*(-4.372547997365551e-2_dp-1.399215359156976e-1_dp*z-  &
1.0494115193677322e-1_dp*z**2)))/z**6)+  &
PL4*((-2.098823038735464e-1_dp*(r0**2-2*r0*rz+rz**2)*(2.5e-1_dp  &
 - z**2)**2)/(d2*r0*rz*z**5)+  &
(1.0494115193677322e-1_dp*rz**2*z05m*z*z05**2+  &
r0**2*(1.0494115193677322e-1_dp*z05m*z*z05**2+  &
rz**2*(1.3117643992096653e-1_dp+2.0988230387354645e-1_dp*z-  &
3.148234558103197e-1_dp*z**2-  &
4.19764607747093e-1_dp*z**3)))/(r0*rz*z**6))+  &
PL3*((-1.049411519367732e-1_dp*(r0**2-2*r0*rz+rz**2)*(2.5e-1_dp  &
 - z**2)**2)/(d2*r0*rz*z**5)+  &
(5.24705759683866e-2_dp*rz**2*z05m*z*z05**2+  &
r0**2*(5.24705759683866e-2_dp*z05m*z*z05**2+  &
rz**2*(6.5588219960483265e-2_dp+1.0494115193677322e-1_dp*z-  &
1.5741172790515985e-1_dp*z**2-  &
2.098823038735465e-1_dp*z**3)))/(r0*rz*z**6))))+  &
l1**3*((d2**3*(4.286875414397245_dp-4.286875414397245_dp/z**2)  &
+(1.2860626243191735e1_dp*d2**2*(-2._dp+r0/rz+rz/r0)*(1-  &
z)**2)/z)/d1**3-(2.572125248638347e1_dp*d2**4*(-2._dp+r0/rz  &
+rz/r0)*(1 - z)**2)/(d1**4*z)+(l4*(2-r0/rz-  &
rz/r0)*(1-2*z)*(1+2*z)*(3.8259794976948562e-3_dp+  &
1.311764399209665e-2_dp*r0*rz+4.372547997365551e-3_dp*z+  &
2.62352879841933e-2_dp*r0*z**1.5_dp-  &
2.62352879841933e-2_dp*z**2))/(z**5*(-1._dp-2*r0*rz+2*z))+  &
(l5*(-2._dp+r0/rz+rz/r0)*(1-2*z)*(1+  &
2*z)*(1.9129897488474281e-3_dp+6.558821996048325e-3_dp*r0*rz+  &
2.1862739986827755e-3_dp*z+1.311764399209665e-2_dp*r0*z**1.5_dp-  &
1.311764399209665e-2_dp*z**2))/(d2*z**5)+((-2._dp+r0/rz+  &
rz/r0)*(8.3693301512075e-4_dp+1.9129897488474281e-3_dp*(3*r0*rz  &
+z+6*r0*z**1.5_dp-6*z**2)+  &
1.0931369993413877e-3_dp*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)**2))/(d2*z**5)+(PL2*(1.9129897488474281e-2_dp+  &
1.5303917990779425e-2_dp*z05m*z-  &
7.6519589953897125e-2_dp*z**2+1.5303917990779425e-2_dp*z*z05  &
-8.745095994731102e-3_dp*(1-2*z)*z*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)+2.1862739986827755e-2_dp*(1-  &
2*z)*(1+2*z)*(3*r0*rz+z+6*r0*z**1.5_dp-6*z**2)  &
+8.745095994731102e-3_dp*z*(1+2*z)*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)-(8.745095994731102e-3_dp*(-1._dp+  &
r0/2/rz+rz/2/r0)*(1-2*z)*z*(1+2*z)*(-1._dp+  &
4*z**2))/(1+r0*rz - z)-(4.372547997365551e-3_dp*(8._dp-  &
(4*r0)/rz-(4*rz)/r0)*(1-2*z)*z*(1+2*z)*(-1._dp+  &
4*z**2))/(-8*r0*rz+8*z)-4.372547997365551e-3_dp*(1-  &
2*z)*z*(1+2*z)*(1+(1.5_dp*r0)/rz+9*r0*rz-1.2e1_dp*z+  &
(1.5_dp*rz+3*z**1.5_dp)/r0)))/z**6+  &
(PL1*(9.56494874423714e-3_dp+7.6519589953897125e-3_dp*z05m*z  &
-3.8259794976948562e-2_dp*z**2+7.6519589953897125e-3_dp*z*z05 &
-4.372547997365551e-3_dp*(1-2*z)*z*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)+1.0931369993413877e-2_dp*(1-  &
2*z)*(1+2*z)*(3*r0*rz+z+6*r0*z**1.5_dp-6*z**2)  &
+4.372547997365551e-3_dp*z*(1+2*z)*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)-(4.372547997365551e-3_dp*(-1._dp+  &
r0/2/rz+rz/2/r0)*(1-2*z)*z*(1+2*z)*(-1._dp+  &
4*z**2))/(1+r0*rz - z)-(2.1862739986827755e-3_dp*(8._dp-  &
(4*r0)/rz-(4*rz)/r0)*(1-2*z)*z*(1+2*z)*(-1._dp+  &
4*z**2))/(-8*r0*rz+8*z)-2.1862739986827755e-3_dp*(1-  &
2*z)*z*(1+2*z)*(1+(1.5_dp*r0)/rz+9*r0*rz-1.2e1_dp*z+  &
(1.5_dp*rz+3*z**1.5_dp)/r0)))/z**6+  &
(rz**3*z**2*(-1.0760567337266782e-3_dp-2.097562413107113e-2_dp*z+  &
5.633896478307894e-3_dp*z**2+1.5919653917875565e-1_dp*z**3-  &
1.427787547922657e-1_dp*z**4)+  &
r0**5*rz**2*(rz**3*(-7.378674745554368e-2_dp-  &
3.935293197628995e-2_dp*z)*z**1.5_dp+6.892981914761471e-2_dp*z**4+  &
rz*z**2.5_dp*(-1.229779124259061e-2_dp-3.4433815479253713e-2_dp*z-  &
5.514385531809179e-1_dp*z**2)+rz**4*(-3.689337372777184e-2_dp+  &
2.4595582485181264e-3_dp*z+1.4757349491108736e-2_dp*z**2)+  &
rz**2*z*(-2.4595582485181207e-3_dp-1.4757349491108736e-2_dp*z-  &
2.7571927659045894e-1_dp*z**2+4.135789148856883e-1_dp*z**3))+  &
r0**4*rz*(rz**3*z**1.5_dp*(-7.378674745554368e-2_dp+  &
1.4757349491108736e-1_dp*z+1.574117279051598e-1_dp*z**2)+  &
z**4*(-1.3785963829522943e-1_dp*z+2.7571927659045894e-1_dp*z**2)+  &
rz**4*(-4.765394106503862e-2_dp+8.0345569451592e-2_dp*z+  &
8.85440969466524e-2_dp*z**2-9.838232994072486e-3_dp*z**3)+  &
rz*z**2.5_dp*(8.198527495060425e-4_dp-3.279410998024165e-3_dp*z-  &
5.2192385419870035e-1_dp*z**2+1.1028771063618357_dp*z**3)+  &
rz**2*z*(-1.3835015147914427e-3_dp-1.608961020905605e-2_dp*z-  &
3.2498299444811316e-1_dp*z**2+9.633363037748557e-1_dp*z**3-  &
8.271578297713766e-1_dp*z**4))+r0*rz*z*((-1.3785963829522943e-1_dp  &
+1.3785963829522943e-1_dp*z)*z**5.5_dp+  &
rz**2*z**2.5_dp*(-1.4757349491108738e-2_dp-  &
1.329405217981932e-1_dp*z+1.4769787128930196e-1_dp*z**2)+  &
rz*z**3*(-2.297660638253824e-2_dp+1.608362446777677e-1_dp*z**2-  &
2.7571927659045894e-1_dp*z**3)+rz**3*(-1.1102172649560968e-3_dp+  &
1.5253735150143553e-2_dp*z-7.338623525471595e-2_dp*z**2-  &
1.1324332641367914e-1_dp*z**3+2.3636634461416897e-1_dp*z**4))+  &
r0**3*(rz**6*(-2.4595582485181207e-3_dp+5.4172469656506e-2_dp*z)*z  &
+rz**5*(-4.919116497036241e-3_dp+  &
1.08344939313012e-1_dp*z)*z**2.5_dp+  &
rz*z**3.5_dp*(-7.378674745554368e-3_dp+2.4595582485181233e-3_dp*z+  &
5.563576696779541e-1_dp*z**2-5.514385531809179e-1_dp*z**3)+  &
z**4*(-2.297660638253824e-2_dp+1.608362446777677e-1_dp*z**2-  &
2.7571927659045894e-1_dp*z**3)+  &
rz**3*z**2.5_dp*(8.690439144764033e-2_dp-4.755145947135037e-2_dp*z-  &
1.1805879592886988e-1_dp*z**2+1.3785963829522943e-1_dp*z**3)+  &
rz**4*(-1.0760567337266784e-2_dp+6.251377214983561e-2_dp*z+  &
4.755145947135034e-2_dp*z**2-1.6643010814972627e-1_dp*z**3-  &
5.404809327829138e-2_dp*z**4)+rz**2*z*(-1.1102172649560968e-3_dp-  &
4.8507954345774085e-3_dp*z-6.764208365908138e-2_dp*z**2+  &
5.301016522973917e-1_dp*z**3-6.367446979222839e-1_dp*z**4-  &
1.3785963829522951e-1_dp*z**5))+  &
r0**2*rz*z*(rz*(1.3785963829522943e-1_dp-  &
2.7571927659045885e-1_dp*z)*z**4.5_dp+  &
rz**3*z**1.5_dp*(8.19852749506041e-3_dp+1.247419943031328e-1_dp*z-  &
2.5604281060231395e-1_dp*z**2)+rz**4*(-1.3835015147914427e-3_dp+  &
7.294473952327962e-2_dp*z-8.372862743146168e-2_dp*z**2-  &
9.358758982190327e-2_dp*z**3)+rz**2*(1.5133115334632334e-2_dp-  &
6.490500933589472e-4_dp*z-9.646934019187745e-2_dp*z**2+  &
4.263234297431411e-2_dp*z**3-1.9800842366359572e-2_dp*z**4+  &
2.7571927659045894e-1_dp*z**5)+z*(-1.0760567337266782e-3_dp-  &
8.710935463501686e-4_dp*z+1.348773073800855e-1_dp*z**2-  &
1.7396425336804888e-1_dp*z**3-2.3468518032241885e-1_dp*z**4+  &
5.514385531809179e-1_dp*z**5)))/(r0*rz**2*(r0*rz - z)*(1+  &
r0*rz - z)*z**6)+l3**2*((1.0931369993413877e-3_dp*(-2._dp+  &
r0/rz+rz/r0)*(1-4*z**2)**2)/(d2*z**5)+  &
(-8.745095994731102e-3_dp*rz**2*z05m*z*z05**2+  &
r0**2*(-8.745095994731102e-3_dp*z05m*z*z05**2+  &
rz**2*(-1.0931369993413877e-2_dp-1.74901919894622e-2_dp*z+  &
2.6235287984193305e-2_dp*z**2+  &
3.498038397892441e-2_dp*z**3)))/(r0*rz*z**6))+  &
l2**2*((4.372547997365551e-3_dp*(-2._dp+r0/rz+rz/r0)*(1-  &
4*z**2)**2)/(d2*z**5)+(rz**2*z*(4.372547997365551e-3_dp+  &
8.745095994731102e-3_dp*z-1.74901919894622e-2_dp*z**2-  &
3.49803839789244e-2_dp*z**3)+  &
r0**2*(-3.498038397892442e-2_dp*z05m*z*z05**2+  &
rz**2*(-4.372547997365551e-2_dp-6.99607679578488e-2_dp*z+  &
1.0494115193677322e-1_dp*z**2+  &
1.3992153591569765e-1_dp*z**3)))/(r0*rz*z**6))+  &
l3*((2.1862739986827755e-3_dp*l5*(-2._dp+r0/rz+rz/r0)*(1-  &
2*z)*(1+2*z)*(-1._dp+4*z**2))/(d2*z**5)+  &
(4.372547997365551e-3_dp*l4*(2-r0/rz-rz/r0)*(1-  &
2*z)*(1+2*z)*(-1._dp+4*z**2))/(z**5*(-1._dp-  &
2*r0*rz+2*z))+((-2._dp+r0/rz+  &
rz/r0)*(1.9129897488474281e-3_dp+6.558821996048325e-3_dp*r0*rz+  &
2.1862739986827755e-3_dp*z+1.311764399209665e-2_dp*r0*z**1.5_dp-  &
1.311764399209665e-2_dp*z**2)*(-1._dp+4*z**2))/(d2*z**5)+  &
(PL2*(-2.1862739986827755e-2_dp+1.049411519367732e-1_dp*z**2-  &
6.996076795784884e-2_dp*z**4))/z**6+  &
(PL1*(-1.0931369993413877e-2_dp+5.24705759683866e-2_dp*z**2-  &
3.498038397892442e-2_dp*z**4))/z**6+  &
(rz**2*z**2*(2.1862739986827755e-3_dp+2.461631188155032e-2_dp*z-  &
3.6094250374634895e-2_dp*z**2-9.901181602587197e-2_dp*z**3+  &
1.0830348052027379e-1_dp*z**4)+  &
r0**5*rz**2*(rz*(-6.558821996048325e-3_dp-  &
1.311764399209665e-2_dp*z)*z+rz**2*(6.5588219960483265e-2_dp+  &
1.0494115193677322e-1_dp*z)*z**1.5_dp+(-6.5588219960483265e-3_dp-  &
1.3117643992096653e-2_dp*z)*z**2.5_dp+  &
rz**3*(3.2794109980241632e-2_dp+3.2794109980241637e-2_dp*z-  &
3.935293197628995e-2_dp*z**2))+  &
r0*rz**2*z**2*(z**1.5_dp*(1.311764399209665e-2_dp+  &
1.311764399209665e-2_dp*z-2.62352879841933e-2_dp*z**2)+  &
rz*(-1.8057489885501994e-2_dp+4.759291826408856e-2_dp*z+  &
8.534760353410464e-2_dp*z**2-1.6413638507216097e-1_dp*z**3))+  &
r0**4*rz*(rz**2*z**1.5_dp*(6.5588219960483265e-2_dp-  &
2.6235287984193305e-2_dp*z-2.0988230387354645e-1_dp*z**2)+  &
z**2.5_dp*(-6.5588219960483265e-3_dp+2.6235287984193305e-2_dp*z**2)  &
+rz*z*(-2.1862739986827746e-3_dp+9.291664494401795e-3_dp*z+  &
1.011980900245548e-1_dp*z**2-3.6073520978265794e-2_dp*z**3)+  &
rz**3*(5.465684996706939e-2_dp-2.186273998682771e-3_dp*z-  &
1.6069113890318398e-1_dp*z**2+1.3117643992096626e-2_dp*z**3))+  &
r0**2*z*(rz**3*z**1.5_dp*(-1.311764399209665e-2_dp+  &
5.24705759683866e-2_dp*z**2)+rz**4*(-2.1862739986827746e-3_dp-  &
1.3684941888136448e-2_dp*z+9.291664494401795e-3_dp*z**2+  &
5.583290455188718e-2_dp*z**3)+rz**2*(-2.62352879841933e-2_dp-  &
4.372547997365553e-3_dp*z+1.3008330292162515e-1_dp*z**2-  &
7.651958995389708e-3_dp*z**3-1.7052937189725648e-1_dp*z**4)+  &
z*(2.1862739986827755e-3_dp+1.6397054990120818e-3_dp*z-  &
1.0502406952224963e-1_dp*z**2+8.480103503443393e-2_dp*z**3+  &
1.6397054990120812e-2_dp*z**4))+  &
r0**3*(rz**5*(-6.558821996048325e-3_dp-  &
1.311764399209665e-2_dp*z)*z+rz**4*(-1.311764399209665e-2_dp-  &
2.62352879841933e-2_dp*z)*z**2.5_dp+  &
z**3.5_dp*(6.5588219960483265e-3_dp+6.5588219960483265e-3_dp*z-  &
1.3117643992096653e-2_dp*z**2)+  &
rz**2*z**2.5_dp*(-6.5588219960483265e-2_dp-  &
3.9352931976289964e-2_dp*z+1.0494115193677322e-1_dp*z**2)+  &
rz*z**2*(4.919116497036243e-3_dp+9.354613102916504e-2_dp*z-  &
1.9037167305635425e-1_dp*z**2+1.9676465988144975e-2_dp*z**3)+  &
rz**3*(2.1862739986827755e-2_dp-5.902939796443494e-2_dp*z-  &
1.4320094691372178e-1_dp*z**2+1.7052937189725648e-1_dp*z**3+  &
1.9676465988144982e-1_dp*z**4)))/(r0*rz*(r0*rz - z)*(1+r0*rz  &
 - z)*z**6))+l2*((4.372547997365551e-3_dp*l5*(-2._dp+r0/rz+  &
rz/r0)*(1-2*z)*(1+2*z)*(-1._dp+4*z**2))/(d2*z**5)+  &
(8.745095994731102e-3_dp*l4*(2-r0/rz-rz/r0)*(1-  &
2*z)*(1+2*z)*(-1._dp+4*z**2))/(z**5*(-1._dp-  &
2*r0*rz+2*z))+((-2._dp+r0/rz+  &
rz/r0)*(3.8259794976948562e-3_dp+1.311764399209665e-2_dp*r0*rz+  &
4.372547997365551e-3_dp*z+2.62352879841933e-2_dp*r0*z**1.5_dp-  &
2.62352879841933e-2_dp*z**2)*(-1._dp+4*z**2))/(d2*z**5)+  &
(PL2*(-4.372547997365551e-2_dp+2.098823038735464e-1_dp*z**2-  &
1.3992153591569767e-1_dp*z**4))/z**6+  &
(PL1*(-2.1862739986827755e-2_dp+1.049411519367732e-1_dp*z**2-  &
6.996076795784884e-2_dp*z**4))/z**6+  &
(rz**2*z**2*(4.372547997365551e-3_dp+4.923262376310064e-2_dp*z-  &
7.218850074926979e-2_dp*z**2-1.9802363205174394e-1_dp*z**3+  &
2.1660696104054757e-1_dp*z**4)+  &
r0**5*rz**2*(rz*(-1.311764399209665e-2_dp-  &
2.62352879841933e-2_dp*z)*z+rz**2*(1.3117643992096653e-1_dp+  &
2.0988230387354645e-1_dp*z)*z**1.5_dp+(-1.3117643992096653e-2_dp-  &
2.6235287984193305e-2_dp*z)*z**2.5_dp+  &
rz**3*(6.5588219960483265e-2_dp+6.558821996048327e-2_dp*z-  &
7.87058639525799e-2_dp*z**2))+  &
r0*rz**2*z**2*(z**1.5_dp*(2.62352879841933e-2_dp+  &
2.62352879841933e-2_dp*z-5.24705759683866e-2_dp*z**2)+  &
rz*(-3.6114979771003988e-2_dp+9.518583652817712e-2_dp*z+  &
1.7069520706820929e-1_dp*z**2-3.2827277014432195e-1_dp*z**3))+  &
r0**4*rz*(rz**2*z**1.5_dp*(1.3117643992096653e-1_dp-  &
5.247057596838661e-2_dp*z-4.197646077470929e-1_dp*z**2)+  &
z**2.5_dp*(-1.3117643992096653e-2_dp+5.247057596838661e-2_dp*z**2)+  &
rz*z*(-4.372547997365549e-3_dp+1.8583328988803591e-2_dp*z+  &
2.023961800491096e-1_dp*z**2-7.214704195653159e-2_dp*z**3)+  &
rz**3*(1.0931369993413877e-1_dp-4.372547997365542e-3_dp*z-  &
3.2138227780636797e-1_dp*z**2+2.623528798419325e-2_dp*z**3))+  &
r0**2*z*(rz**3*z**1.5_dp*(-2.62352879841933e-2_dp+  &
1.049411519367732e-1_dp*z**2)+rz**4*(-4.372547997365549e-3_dp-  &
2.7369883776272896e-2_dp*z+1.8583328988803591e-2_dp*z**2+  &
1.1166580910377435e-1_dp*z**3)+rz**2*(-5.24705759683866e-2_dp-  &
8.745095994731106e-3_dp*z+2.601666058432503e-1_dp*z**2-  &
1.5303917990779417e-2_dp*z**3-3.4105874379451295e-1_dp*z**4)+  &
z*(4.372547997365551e-3_dp+3.2794109980241637e-3_dp*z-  &
2.1004813904449926e-1_dp*z**2+1.6960207006886787e-1_dp*z**3+  &
3.2794109980241624e-2_dp*z**4))+  &
r0**3*(rz**5*(-1.311764399209665e-2_dp-2.62352879841933e-2_dp*z)*z  &
+rz**4*(-2.62352879841933e-2_dp-5.24705759683866e-2_dp*z)*z**2.5_dp  &
+z**3.5_dp*(1.3117643992096653e-2_dp+1.3117643992096653e-2_dp*z-  &
2.6235287984193305e-2_dp*z**2)+  &
rz**2*z**2.5_dp*(-1.3117643992096653e-1_dp-  &
7.870586395257993e-2_dp*z+2.0988230387354645e-1_dp*z**2)+  &
rz*z**2*(9.838232994072486e-3_dp+1.8709226205833007e-1_dp*z-  &
3.807433461127085e-1_dp*z**2+3.935293197628995e-2_dp*z**3)+  &
rz**3*(4.372547997365551e-2_dp-1.1805879592886988e-1_dp*z-  &
2.8640189382744357e-1_dp*z**2+3.4105874379451295e-1_dp*z**3+  &
3.9352931976289964e-1_dp*z**4)))/(r0*rz*(r0*rz - z)*(1+r0*rz  &
 - z)*z**6)+l3*((4.372547997365551e-3_dp*(-2._dp+r0/rz+  &
rz/r0)*(1-4*z**2)**2)/(d2*z**5)+  &
(rz**2*z*(4.372547997365551e-3_dp+8.745095994731102e-3_dp*z-  &
1.74901919894622e-2_dp*z**2-3.49803839789244e-2_dp*z**3)+  &
r0**2*(-3.498038397892442e-2_dp*z05m*z*z05**2+  &
rz**2*(-4.372547997365551e-2_dp-6.99607679578488e-2_dp*z+  &
1.0494115193677322e-1_dp*z**2+  &
1.3992153591569765e-1_dp*z**3)))/(r0*rz*z**6))))+  &
l1**2*((-3.143051358559554e2_dp*d2**2*(r0**2-2*r0*rz+  &
rz**2)*(r0*rz - z05m)*(-1._dp+  &
z)*(-3.7913315126399993e-1_dp+z)*(1+z)*(rz+r0*z05m-  &
z**1.5_dp))/(d1**4*r0*z**2)+  &
(PL2**2*(-8.745095994731102e-2_dp+4.197646077470929e-1_dp*z**2-  &
2.7984307183139534e-1_dp*z**4))/z**6+  &
(PL1**2*(-2.1862739986827755e-2_dp+1.0494115193677322e-1_dp*z**2-  &
6.996076795784884e-2_dp*z**4))/z**6+(d2**2*(1.851130069975489_dp  &
+1.851130069975489_dp/z**2-3.702260139950978_dp*z)-  &
(3.702260139950978_dp*d2*(-2._dp+r0/rz+rz/r0)*(1 - z)**2*(1  &
+z))/z)/d1**2+(1.4757349491108736e-2_dp*(r0**2-2*r0*rz+  &
rz**2)*(2.802528429206536e1_dp*r0*rz*z**3+  &
r0*rz**3*(1.5_dp+z)+z**3*(8.174041251852396_dp+  &
9.341761430688454_dp*z+5.605056858413072e1_dp*r0*z**1.5_dp-  &
5.605056858413072e1_dp*z**2)+rz**2*(4.374999999999999e-1_dp+  &
7.916666666666665e-1_dp*z+2.999999999999999_dp*r0*z**1.5_dp-  &
2.666666666666666_dp*z**2+2*r0*z**2.5_dp-  &
2*z**3)))/(d2*rz**2*z**5)+(PL3*(-1.434742311635571e-2_dp+  &
3.4433815479253713e-2_dp*z**2-  &
2.62352879841933e-2_dp*z**2*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)-1.6397054990120816e-2_dp*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)*(1-4*z**2)-  &
(2.62352879841933e-2_dp*(r0**2-2*r0*rz+rz**2)*z*(2.5e-1_dp-  &
z**2)**2)/(r0*rz*(r0*rz - z))-  &
(5.24705759683866e-2_dp*(r0**2-2*r0*rz+rz**2)*z*(2.5e-1_dp-  &
z**2)**2)/(r0*rz*(1+r0*rz - z))+  &
3.2794109980241624e-3_dp*z*(1-4*z**2)*(1+(1.5_dp*r0)/rz+  &
9*r0*rz-1.2e1_dp*z+(1.5_dp*rz+3*z**1.5_dp)/r0)))/z**6+  &
(PL4*(-2.869484623271142e-2_dp+6.886763095850743e-2_dp*z**2-  &
5.24705759683866e-2_dp*z**2*(3*r0*rz+z+6*r0*z**1.5_dp-  &
6*z**2)-3.2794109980241632e-2_dp*(3*r0*rz+z+  &
6*r0*z**1.5_dp-6*z**2)*(1-4*z**2)-  &
(5.24705759683866e-2_dp*(r0**2-2*r0*rz+rz**2)*z*(2.5e-1_dp-  &
z**2)**2)/(r0*rz*(r0*rz - z))-  &
(1.049411519367732e-1_dp*(r0**2-2*r0*rz+rz**2)*z*(2.5e-1_dp-  &
z**2)**2)/(r0*rz*(1+r0*rz - z))+  &
6.558821996048325e-3_dp*z*(1-4*z**2)*(1+(1.5_dp*r0)/rz+  &
9*r0*rz-1.2e1_dp*z+(1.5_dp*rz+3*z**1.5_dp)/r0)))/z**6+  &
l5*((-1.399215359156976e-1_dp*PL1*(r0**2-2*r0*rz+  &
rz**2)*(2.5e-1_dp - z**2)**2)/(d2*r0*rz*z**5)-  &
(2.798430718313952e-1_dp*PL2*(r0**2-2*r0*rz+rz**2)*(2.5e-1_dp-  &
z**2)**2)/(d2*r0*rz*z**5)+(1.311764399209665e-2_dp*(r0**2-  &
2*r0*rz+rz**2)*z05*(r0*rz**3+(2.802528429206536e1_dp-  &
5.605056858413072e1_dp*z)*z**3+rz**2*(1.791666666666667_dp-  &
1.666666666666667_dp*z+2*r0*z**1.5_dp-  &
4*z**2)))/(d2*rz**2*z**5))+  &
(r0**5*rz**2*((-1.7232454786903677e-1_dp-  &
4.825087340333031e-1_dp*z)*z**4+rz**4*(-2.767003029582887e-2_dp-  &
2.951469898221747e-2_dp*z-7.378674745554367e-3_dp*z**2)+  &
rz**2*z*(-9.223343431942958e-3_dp-3.197425723073559e-2_dp*z-  &
4.307958226253152e-1_dp*z**2-1.3785963829522947e-1_dp*z**3))+  &
r0**4*rz*(z**4*(1.1488303191269131e-2_dp-4.59532127650765e-2_dp*z+  &
4.135789148856883e-1_dp*z**2)+  &
rz**3*z**1.5_dp*(-5.913057942708661e-2_dp-8.526468594862823e-2_dp*z  &
+1.1731780813982664e-1_dp*z**2+3.6245317713344916e-1_dp*z**3)+  &
rz**4*(-5.723532000937217e-2_dp+3.04468914454696e-2_dp*z+  &
2.209897484721094e-1_dp*z**2+1.7667796596420995e-1_dp*z**3-  &
2.7183988285008684e-1_dp*z**4)+rz**2*z*(3.5714185334838575e-3_dp+  &
3.2794109980241624e-3_dp*z-4.028590801596895e-1_dp*z**2+  &
1.2101874470064926_dp*z**3+1.3785963829522947_dp*z**4))+  &
r0*rz**2*z*(1.0594228814019682e-2_dp+9.334938429820554e-3_dp*z-  &
8.730743856580952e-2_dp*z**2-6.7806313680500025e-2_dp*z**3+  &
1.733484471746268e-1_dp*z**4+6.816925857291789e-2_dp*z**5-  &
2.9352335798583438e-1_dp*z**6+rz*z**1.5_dp*(5.913057942708658e-3_dp  &
-1.1680301418387882e-3_dp*z-4.0869139510461485e-2_dp*z**2-  &
5.150735823675343e-2_dp*z**3+1.8122658856672458e-1_dp*z**4)+  &
rz**2*(-1.1447064001874434e-2_dp*z-1.040713672659719e-1_dp*z**2+  &
8.025316558130509e-2_dp*z**3+1.2289673563288628e-1_dp*z**4-  &
1.8122658856672458e-1_dp*z**5))+  &
r0**3*(rz**6*z*(-9.223343431942958e-3_dp+7.142047149068651e-2_dp*z  &
+5.171291140798787e-2_dp*z**2)+z**5*(-1.033947287214221e-1_dp+  &
3.4464909573807354e-2_dp*z+6.8929819147614735e-2_dp*z**2)+  &
rz**3*z**1.5_dp*(-5.913057942708661e-2_dp+3.299647290554497e-2_dp*z  &
+2.878471800370831e-1_dp*z**2+1.278175608537958e-1_dp*z**3-  &
7.249063542668983e-1_dp*z**4)+rz**2*z*(-1.1447064001874434e-2_dp*z  &
-6.766385445498013e-4_dp*z**2+6.43180021953492e-1_dp*z**3-  &
4.055652111654932e-1_dp*z**4-1.4219633332237898_dp*z**5)+  &
rz**4*(-3.8188499213326765e-2_dp+5.142716294634666e-2_dp*z+  &
1.8091386621445444e-1_dp*z**2-3.789747231315446e-2_dp*z**3-  &
7.449329110534121e-1_dp*z**4+6.117103166687077e-2_dp*z**5))+  &
rz**2*z*(z**2.5_dp*(-5.913057942708658e-3_dp-  &
4.74502780086987e-3_dp*z+3.0211053766882956e-2_dp*z**2+  &
7.106032626005786e-2_dp*z**3-9.06132942833623e-2_dp*z**4)+  &
rz*(-9.855096571181098e-4_dp-1.7763476239297544e-3_dp*z+  &
1.0157395603710841e-2_dp*z**2+1.5710533195984677e-2_dp*z**3-  &
4.412796751434556e-2_dp*z**4-6.660957395064714e-2_dp*z**5+  &
1.8122658856672458e-1_dp*z**6))+  &
r0**2*rz*(rz**4*z*(3.5714185334838575e-3_dp+  &
1.0667413971944626e-1_dp*z-1.2713980356923065e-1_dp*z**2-  &
1.6840893594580186e-1_dp*z**3)+  &
rz**3*z**2.5_dp*(5.913057942708658e-3_dp+1.065808574357853e-2_dp*z-  &
1.955296802330443e-2_dp*z**2-9.06132942833623e-2_dp*z**3)+  &
rz*z**2.5_dp*(5.913057942708661e-2_dp+2.6134106521541627e-2_dp*z-  &
2.025824940884549e-1_dp*z**2-2.4513536899362247e-1_dp*z**3+  &
3.6245317713344916e-1_dp*z**4)+z*(-9.855096571181098e-4_dp-  &
1.7763476239297544e-3_dp*z+1.0157395603710841e-2_dp*z**2+  &
1.5710533195984677e-2_dp*z**3-4.412796751434556e-2_dp*z**4-  &
6.660957395064714e-2_dp*z**5+1.8122658856672458e-1_dp*z**6)+  &
rz**2*(-8.62320949978346e-3_dp+3.2406256031365723e-2_dp*z+  &
8.855209535819451e-2_dp*z**2-1.3406433811862728e-1_dp*z**3-  &
3.796441818583479e-1_dp*z**4+4.444692359126824e-1_dp*z**5+  &
5.041922091690504e-1_dp*z**6)))/(r0*rz**2*(r0*rz - z)*(1+  &
r0*rz - z)*z**6)+l2**2*((-2.098823038735464e-1_dp*(r0**2-  &
2*r0*rz+rz**2)*z05*(-2.5e-1_dp+z**2))/(d2*z**5)+  &
(6.99607679578488e-2_dp*rz**2*z*z05**2+  &
r0**2*(6.99607679578488e-2_dp*z*z05**2+  &
rz**2*(-8.745095994731102e-2_dp-2.798430718313952e-1_dp*z-  &
2.0988230387354645e-1_dp*z**2)))/z**6)+  &
l3**2*((-5.24705759683866e-2_dp*(r0**2-2*r0*rz+rz**2)*(0.5_dp  &
+z)*(-2.5e-1_dp+z**2))/(d2*z**5)+  &
(1.74901919894622e-2_dp*rz**2*z*z05**2+  &
r0**2*(1.74901919894622e-2_dp*z*z05**2+  &
rz**2*(-2.1862739986827755e-2_dp-6.99607679578488e-2_dp*z-  &
5.247057596838661e-2_dp*z**2)))/z**6)+  &
l4*((-1.74901919894622e-2_dp*PL1*(2-r0/rz-  &
rz/r0)*(1-4*z**2)**2)/(z**5*(-1._dp-2*r0*rz+  &
2*z))-(3.49803839789244e-2_dp*PL2*(2-r0/rz-  &
rz/r0)*(1-4*z**2)**2)/(z**5*(-1._dp-2*r0*rz+  &
2*z))+(r0**3*rz**2*(rz*(6.558821996048325e-3_dp+  &
1.311764399209665e-2_dp*z)+(1.311764399209665e-2_dp+  &
2.62352879841933e-2_dp*z)*z**1.5_dp)+  &
rz**2*(1.8381285106030592e-1_dp*z**3-7.352514042412237e-1_dp*z**5+  &
rz**2*(1.1751222742919918e-2_dp+1.2571075492425958e-2_dp*z-  &
4.8098027971021065e-2_dp*z**2-5.24705759683866e-2_dp*z**3))+  &
r0**2*(rz**4*(-1.311764399209665e-2_dp-2.62352879841933e-2_dp*z)+  &
rz**3*(-2.62352879841933e-2_dp-5.24705759683866e-2_dp*z)*z**1.5_dp+  &
1.8381285106030592e-1_dp*z**3-7.352514042412237e-1_dp*z**5+  &
rz**2*(1.1751222742919918e-2_dp+1.2571075492425958e-2_dp*z-  &
4.8098027971021065e-2_dp*z**2-5.24705759683866e-2_dp*z**3))+  &
r0*rz*(rz**4*(6.558821996048325e-3_dp+1.311764399209665e-2_dp*z)+  &
rz**3*(1.311764399209665e-2_dp+2.62352879841933e-2_dp*z)*z**1.5_dp-  &
3.6762570212061183e-1_dp*z**3+1.4705028084824474_dp*z**5+  &
rz**2*(-2.3502445485839836e-2_dp-2.5142150984851916e-2_dp*z+  &
9.619605594204213e-2_dp*z**2+  &
1.049411519367732e-1_dp*z**3)))/(rz**2*(0.5_dp+r0*rz-  &
z)*z**5))+PL1*((PL2*(-8.745095994731102e-2_dp+  &
4.197646077470929e-1_dp*z**2-2.7984307183139534e-1_dp*z**4))/z**6+  &
(3.935293197628995e-2_dp*(r0**2-2*r0*rz+  &
rz**2)*(-7.291666666666668e-2_dp-8.333333333333336e-2_dp*z+  &
7.916666666666666e-1_dp*z**2+3.333333333333334e-1_dp*z**3-  &
2*z**4+r0*(rz*(-2.5e-1_dp+z**2)+z**1.5_dp*(-0.5_dp+  &
2*z**2))))/(d2*r0*rz*z**5)+  &
(rz**2*z**2*(5.875611371459959e-3_dp+9.231635190490595e-2_dp*z-  &
1.2224097726187646e-1_dp*z**2-3.698119761192946e-1_dp*z**3+  &
3.9386099010480513e-1_dp*z**4)+  &
r0**5*rz**2*(rz*(-6.558821996048325e-3_dp-  &
1.311764399209665e-2_dp*z)*z+rz**2*(6.5588219960483265e-2_dp+  &
1.0494115193677322e-1_dp*z)*z**1.5_dp+(-6.5588219960483265e-3_dp-  &
1.3117643992096653e-2_dp*z)*z**2.5_dp+  &
rz**3*(3.2794109980241632e-2_dp+3.2794109980241637e-2_dp*z-  &
3.935293197628995e-2_dp*z**2))+  &
r0**4*rz*(rz**2*z**1.5_dp*(6.5588219960483265e-2_dp-  &
2.623528798419331e-2_dp*z-2.0988230387354645e-1_dp*z**2)+  &
z**2.5_dp*(-6.5588219960483265e-3_dp+2.6235287984193305e-2_dp*z**2)  &
+rz*z*(3.9626216226125277e-3_dp+2.650857223402865e-2_dp*z+  &
3.5232178412983246e-1_dp*z**2-1.049411519367732e-1_dp*z**3)+  &
rz**3*(9.155022369484122e-2_dp+1.7490191989462205e-2_dp*z-  &
2.4923523584983642e-1_dp*z**2-2.623528798419331e-2_dp*z**3))+  &
r0*rz**2*z*(z**2.5_dp*(1.311764399209665e-2_dp+  &
1.311764399209665e-2_dp*z-2.62352879841933e-2_dp*z**2)+  &
rz*(-4.099263747530217e-4_dp-8.534760353410464e-2_dp*z+  &
1.7233491256722133e-1_dp*z**2+3.5450805812851516e-1_dp*z**3-  &
6.565455402886437e-1_dp*z**4))+  &
r0**2*z*(rz**3*z**1.5_dp*(-1.311764399209665e-2_dp+  &
5.24705759683866e-2_dp*z**2)+rz**4*(3.9626216226125277e-3_dp-  &
6.53978532961243e-2_dp*z-1.5303917990779432e-2_dp*z**2+  &
2.6268455018383863e-1_dp*z**3)+rz**2*(-6.968748370801345e-2_dp+  &
1.9403181738309612e-2_dp*z+3.0389208581690577e-1_dp*z**2-  &
8.307841194994543e-2_dp*z**3-3.672940317787062e-1_dp*z**4)+  &
z*(5.875611371459959e-3_dp+4.099263747530189e-4_dp*z-  &
3.979602538523353e-1_dp*z**2+3.6543942812192904e-1_dp*z**3+  &
2.62352879841933e-2_dp*z**4))+  &
r0**3*(rz**5*(-6.558821996048325e-3_dp-  &
1.311764399209665e-2_dp*z)*z+rz**4*(-1.311764399209665e-2_dp-  &
2.62352879841933e-2_dp*z)*z**2.5_dp+  &
z**3.5_dp*(6.5588219960483265e-3_dp+6.5588219960483265e-3_dp*z-  &
1.3117643992096653e-2_dp*z**2)+  &
rz**2*z**2.5_dp*(-6.5588219960483265e-2_dp-  &
3.935293197628995e-2_dp*z+1.0494115193677322e-1_dp*z**2)+  &
rz*z*(-4.099263747530217e-4_dp+6.5588219960483265e-3_dp*z+  &
3.5614776362752725e-1_dp*z**2-7.483690482333203e-1_dp*z**3+  &
7.87058639525799e-2_dp*z**4)+rz**3*(5.875611371459958e-2_dp-  &
1.3281614541997857e-1_dp*z-3.104509078129541e-1_dp*z**2+  &
3.869704977668512e-1_dp*z**3+  &
4.328822517391895e-1_dp*z**4)))/(r0*rz*(r0*rz - z)*(1+r0*rz-  &
z)*z**6))+PL2*(((-2._dp+r0/rz+rz/r0)*(1-2*z)*(1+  &
2*z)*(-1.147793849308457e-2_dp-3.935293197628995e-2_dp*r0*rz-  &
1.311764399209665e-2_dp*z-7.87058639525799e-2_dp*r0*z**1.5_dp+  &
7.87058639525799e-2_dp*z**2))/(d2*z**5)+  &
(rz**2*z**2*(3.006126748188817e-3_dp+8.015520278723303e-2_dp*z-  &
2.7243986559243467e-1_dp*z**2-1.433667820846745e-1_dp*z**3+  &
1.0416614343987176_dp*z**4-7.090161162570303e-1_dp*z**5)+  &
r0**6*rz**3*(rz*(-1.311764399209665e-2_dp-  &
2.62352879841933e-2_dp*z)*z+rz**2*(1.3117643992096653e-1_dp+  &
2.0988230387354645e-1_dp*z)*z**1.5_dp+(-1.3117643992096653e-2_dp-  &
2.6235287984193305e-2_dp*z)*z**2.5_dp+  &
rz**3*(6.5588219960483265e-2_dp+6.558821996048327e-2_dp*z-  &
7.87058639525799e-2_dp*z**2))+  &
r0**5*rz**2*(rz**2*z**1.5_dp*(1.9676465988144982e-1_dp-  &
7.870586395257995e-2_dp*z-6.296469116206392e-1_dp*z**2)+  &
rz*z*(1.1204654243249221e-2_dp+5.30171444680573e-2_dp*z+  &
6.915259242675683e-1_dp*z**2-2.098823038735464e-1_dp*z**3)+  &
rz**3*(2.1589455736992407e-1_dp+2.1862739986827933e-3_dp*z-  &
6.0341162363644605e-1_dp*z**2+2.623528798419328e-2_dp*z**3))+  &
r0**4*rz*(rz**5*(-1.311764399209665e-2_dp-  &
2.62352879841933e-2_dp*z)*z+rz**4*(-2.62352879841933e-2_dp-  &
5.24705759683866e-2_dp*z)*z**2.5_dp+  &
z**2.5_dp*(1.311764399209665e-2_dp-1.3117643992096648e-2_dp*z-  &
3.9352931976289947e-2_dp*z**2+7.870586395257989e-2_dp*z**3)+  &
rz**2*z**1.5_dp*(6.5588219960483265e-2_dp-  &
3.2794109980241637e-1_dp*z-2.3611759185773975e-1_dp*z**2+  &
7.87058639525799e-1_dp*z**3)+rz*z*(1.5850486490450115e-2_dp+  &
1.5303917990779425e-2_dp*z+9.410928304593103e-1_dp*z**2-  &
2.2407345967025956_dp*z**3+4.459998957312862e-1_dp*z**4)+  &
rz**3*(2.090624511240404e-1_dp-4.509190122283225e-1_dp*z-  &
9.05117435454669e-1_dp*z**2+1.3248820432017618_dp*z**3+  &
9.182350794467656e-1_dp*z**4))+  &
r0**3*(rz**4*z**2.5_dp*(-1.9676465988144973e-2_dp+  &
7.870586395257989e-2_dp*z**2)+  &
rz**2*z**2.5_dp*(-1.0494115193677322e-1_dp+  &
1.7052937189725648e-1_dp*z+3.4105874379451295e-1_dp*z**2-  &
5.24705759683866e-1_dp*z**3)+z**3.5_dp*(-1.311764399209665e-2_dp+  &
1.311764399209665e-2_dp*z+5.24705759683866e-2_dp*z**2-  &
5.24705759683866e-2_dp*z**3)+rz**5*z*(1.1204654243249221e-2_dp-  &
1.307957065922486e-1_dp*z-4.372547997365551e-2_dp*z**2+  &
5.2536910036767726e-1_dp*z**3)+rz**3*(5.875611371459958e-2_dp-  &
4.15118775499892e-1_dp*z+2.6781856483864e-2_dp*z**2+  &
1.7566711579416099_dp*z**3-6.383920076153702e-1_dp*z**4-  &
1.7577642949409513_dp*z**5)+rz*z*(2.45955824851812e-3_dp+  &
6.832106245883679e-3_dp*z+3.1597497890173134e-1_dp*z**2-  &
2.1680409862463934_dp*z**3+2.371911036623562_dp*z**4-  &
2.62352879841933e-1_dp*z**5))+  &
r0**2*z*(rz**3*(6.5588219960483265e-3_dp+  &
1.3117643992096648e-2_dp*z)*z**1.5_dp+  &
rz*z**2.5_dp*(3.935293197628995e-2_dp-3.935293197628995e-2_dp*z-  &
1.574117279051598e-1_dp*z**2+1.574117279051598e-1_dp*z**3)+  &
rz**4*(1.5850486490450115e-2_dp-2.604153585996794e-1_dp*z+  &
3.8965427727839232e-1_dp*z**2+1.067896722382911_dp*z**3-  &
1.7597543169923848_dp*z**4)+z*(3.006126748188817e-3_dp-  &
1.1751222742919918e-2_dp*z-3.643462911225877e-1_dp*z**2+  &
1.1433231753374669_dp*z**3-7.964670762043414e-1_dp*z**4+  &
2.6235287984193305e-2_dp*z**5)+rz**2*(-7.542645295455575e-2_dp+  &
1.8337373163951776e-1_dp*z+3.208357093066973e-1_dp*z**2-  &
8.679507774770618e-1_dp*z**3-3.323136477997819e-1_dp*z**4+  &
1.049411519367732_dp*z**5))+  &
r0*rz*z*(rz*z**2.5_dp*(-6.5588219960483265e-3_dp+  &
6.5588219960483265e-3_dp*z+2.6235287984193305e-2_dp*z**2-  &
2.6235287984193305e-2_dp*z**3)+z*(5.738969246542284e-3_dp+  &
8.1985274950604e-4_dp*z-6.886763095850741e-2_dp*z**2+  &
3.6073520978265785e-2_dp*z**3+1.836470158893531e-1_dp*z**4-  &
1.574117279051598e-1_dp*z**5)+rz**2*(2.45955824851812e-3_dp-  &
8.507431928426927e-2_dp*z+4.997878299620371e-1_dp*z**2-  &
1.4609962458302805e-1_dp*z**3-2.03959738882378_dp*z**4+  &
1.9434013328817377_dp*z**5)))/(r0*rz*(r0*rz - z)*(0.5_dp+  &
r0*rz - z)*(1+r0*rz - z)*z**6))+  &
((1.2860626243191735e1_dp*d2**2*(-2._dp+r0/rz+rz/r0)*(1-  &
z)**2)/z+(7.404520279901956_dp*d2**3*(-2._dp+r0/rz+  &
rz/r0)*(1 - z)**2*(1+z))/z+(5.238418930932589e1_dp*(r0**2-  &
2*r0*rz+rz**2)*(r0*rz - z05m)*(-1._dp+  &
z)*(-3.7913315126399993e-1_dp+z)*(1+z)*(rz+r0*z05m-  &
z**1.5_dp))/(r0*z**2)+  &
(d2*(rz**2*z*(-2.4825728461568324_dp+1.6478315048293066e1_dp*z-  &
3.3639813193133445e1_dp*z**2+9.713779606369878_dp*z**3+  &
3.6122386039290273e1_dp*z**4-2.619209465466294e1_dp*z**5+  &
rz*z**1.5_dp*(-9.930291384627331_dp+2.6192094654662945e1_dp*z+  &
9.930291384627331_dp*z**2-2.6192094654662945e1_dp*z**3)+  &
rz**2*(9.930291384627331_dp-2.6192094654662945e1_dp*z-  &
9.930291384627331_dp*z**2+2.6192094654662945e1_dp*z**3))+  &
r0**3*rz*(z*(-9.930291384627331_dp+4.605267742391761e1_dp*z-  &
4.245389792469856e1_dp*z**2-4.605267742391761e1_dp*z**3+  &
5.238418930932589e1_dp*z**4)+rz**2*(1.9860582769254662e1_dp-  &
4.605267742391762e1_dp*z-4.605267742391762e1_dp*z**3+  &
1.0476837861865176e2_dp*z**4))+  &
r0**2*(rz*z**2.5_dp*(-1.9860582769254662e1_dp+  &
5.238418930932589e1_dp*z+1.9860582769254662e1_dp*z**2-  &
5.238418930932589e1_dp*z**3)+  &
rz**3*z**1.5_dp*(3.9721165538509324e1_dp-5.238418930932589e1_dp*z-  &
5.238418930932588e1_dp*z**3)+rz**4*(-3.9721165538509324e1_dp+  &
2.2593315155443907e1_dp*z+7.857628396398883e1_dp*z**2+  &
8.217506346320786e1_dp*z**3-7.857628396398883e1_dp*z**4)+  &
rz**2*(9.930291384627331_dp-3.165755942704145_dp*z-  &
7.857628396398884e1_dp*z**2-1.0363314941142248e1_dp*z**3+  &
2.2306582812115e2_dp*z**4-1.5715256792797763e2_dp*z**5)+  &
z*(-2.482572846156833_dp+1.6478315048293066e1_dp*z-  &
3.3639813193133445e1_dp*z**2+9.713779606369878_dp*z**3+  &
3.6122386039290277e1_dp*z**4-2.6192094654662945e1_dp*z**5))+  &
r0*(rz*z*(9.930291384627331_dp-4.60526774239176e1_dp*z+  &
4.245389792469856e1_dp*z**2+4.60526774239176e1_dp*z**3-  &
5.238418930932589e1_dp*z**4)+z**2.5_dp*(-4.965145692313666_dp+  &
2.3026338711958805e1_dp*z-2.122694896234928e1_dp*z**2-  &
2.3026338711958805e1_dp*z**3+2.6192094654662945e1_dp*z**4)+  &
rz**2*z**1.5_dp*(1.9860582769254662e1_dp-4.605267742391761e1_dp*z-  &
4.6052677423917645e1_dp*z**3+1.0476837861865178e2_dp*z**4)+  &
rz**3*(-1.9860582769254662e1_dp+2.122694896234928e1_dp*z+  &
1.15131693559794e2_dp*z**2-6.008206738782878e1_dp*z**3-  &
2.1990007217844583e2_dp*z**4+  &
1.309604732733147e2_dp*z**5))))/(r0*rz*z**3))/d1**3+  &
l2*((PL3*(3.2794109980241632e-2_dp-1.5741172790515985e-1_dp*z**2+  &
1.049411519367733e-1_dp*z**4))/z**6+(PL4*(6.5588219960483265e-2_dp  &
-3.148234558103197e-1_dp*z**2+2.098823038735466e-1_dp*z**4))/z**6-  &
(7.87058639525799e-2_dp*(r0**2-2*r0*rz+  &
rz**2)*(3.5031605365081706_dp*z**3-1.4012642146032683e1_dp*z**5+  &
r0*rz**3*z05+rz**2*(3.333333333333334e-1_dp+  &
5.833333333333334e-1_dp*z+r0*z**1.5_dp-1.4166666666666665_dp*z**2  &
+2*r0*z**2.5_dp-2.5_dp*z**3)))/(d2*rz**2*z**5)+  &
(l5*(1.7490191989462198e-1_dp*r0**2*z05m*z05**2+  &
1.7490191989462198e-1_dp*rz**2*z05m*z05**2+  &
r0*rz*(4.37254799736555e-2_dp+8.7450959947311e-2_dp*z-  &
1.74901919894622e-1_dp*z**2-  &
3.49803839789244e-1_dp*z**3)))/(d2*z**5)+  &
(l4*(1.7490191989462198e-1_dp*r0**2*z05m*z05**2+  &
1.7490191989462198e-1_dp*rz**2*z05m*z05**2+  &
r0*rz*(4.37254799736555e-2_dp+8.7450959947311e-2_dp*z-  &
1.74901919894622e-1_dp*z**2-3.49803839789244e-1_dp*z**3)))/((0.5_dp  &
+r0*rz - z)*z**5)+(rz**2*z*(-1.9710193142362198e-2_dp+  &
(-8.711368840513877e-3_dp+1.9676465988144982e-2_dp*rz**2)*z+  &
(1.1483196157115425e-1_dp+2.1660696104054757e-1_dp*rz**2)*z**2+  &
(9.125045008862365e-2_dp+1.575775630761126e-1_dp*rz**2)*z**3+  &
(-4.598924667874115e-2_dp-3.9386099010480513e-1_dp*rz**2)*z**4+  &
2.359540991224512e-1_dp*z**5-3.6762570212061196e-1_dp*z**6)+  &
r0*rz**3*(1.9710193142362198e-2_dp+(-1.0998824301848322e-2_dp-  &
1.0931369993413877e-2_dp*rz**2)*z+(-1.4325352355403032e-1_dp-  &
1.8818539905767144e-1_dp*rz**2)*z**2+(-4.840050500345527e-3_dp+  &
8.74509599473111e-3_dp*rz**2)*z**3+(2.236500963556431e-1_dp+  &
6.8278082827283715e-1_dp*rz**2)*z**4-1.0428249612429052e-1_dp*z**5  &
+7.352514042412239e-1_dp*z**6)+  &
r0**4*rz**2*((-1.8381285106030592e-1_dp-  &
3.6762570212061183e-1_dp*z)*z**4+rz**4*(9.838232994072492e-2_dp+  &
2.0988230387354645e-1_dp*z+7.87058639525799e-2_dp*z**2)+  &
rz**2*z*(6.558821996048319e-3_dp+5.247057596838658e-2_dp*z+  &
8.139572681938036e-1_dp*z**2+7.352514042412237e-1_dp*z**3))+  &
r0**3*rz*(-1.8381285106030592e-1_dp*z**4+  &
7.352514042412237e-1_dp*z**6+rz**4*(9.838232994072492e-2_dp-  &
3.935293197628995e-2_dp*z-5.509410476680593e-1_dp*z**2-  &
3.672940317787062e-1_dp*z**3)+rz**2*z*(-1.0931369993413877e-2_dp-  &
4.372547997365551e-3_dp*z+7.439965002359549e-1_dp*z**2-  &
7.877219802096105e-1_dp*z**3-1.4705028084824474_dp*z**4))+  &
r0**2*(z**5*(1.8381285106030592e-1_dp+1.8381285106030592e-1_dp*z-  &
3.6762570212061183e-1_dp*z**2)+rz**6*z*(6.558821996048319e-3_dp-  &
1.3134227509191931e-1_dp*z-2.889198381680319e-1_dp*z**2)+  &
rz**2*z**2*(1.9676465988144982e-2_dp+3.2794109980241632e-2_dp*z-  &
7.61486692225417e-1_dp*z**2-2.623528798419328e-2_dp*z**3+  &
7.352514042412237e-1_dp*z**4)+rz**4*(1.9710193142362198e-2_dp-  &
8.745095994731102e-2_dp*z-2.1540056551056193e-1_dp*z**2+  &
9.343712615976228e-2_dp*z**3+1.569165648279657e-1_dp*z**4-  &
3.6762570212061196e-1_dp*z**5)))/(rz**2*(r0*rz - z)*(1+r0*rz  &
 - z)*z**6)+l3*((-2.098823038735464e-1_dp*(r0**2-2*r0*rz+  &
rz**2)*z05*(-2.5e-1_dp+z**2))/(d2*z**5)+  &
(6.99607679578488e-2_dp*rz**2*z*z05**2+  &
r0**2*(6.99607679578488e-2_dp*z*z05**2+  &
rz**2*(-8.745095994731102e-2_dp-2.798430718313952e-1_dp*z-  &
2.0988230387354645e-1_dp*z**2)))/z**6)+  &
PL1*((r0*rz*(-1.311764399209665e-2_dp+1.049411519367732e-1_dp*z**2  &
-2.098823038735464e-1_dp*z**4)+r0**2*(6.558821996048325e-3_dp-  &
5.24705759683866e-2_dp*z**2+1.049411519367732e-1_dp*z**4)+  &
rz**2*(6.558821996048325e-3_dp-5.24705759683866e-2_dp*z**2+  &
1.049411519367732e-1_dp*z**4))/(d2*r0*rz*z**5)+  &
(-8.745095994731102e-2_dp*rz**2*z05m*z*z05**2+  &
r0**2*(-8.745095994731102e-2_dp*z05m*z*z05**2+  &
rz**2*(-1.0931369993413877e-1_dp-1.7490191989462203e-1_dp*z+  &
2.6235287984193305e-1_dp*z**2+  &
3.498038397892441e-1_dp*z**3)))/(r0*rz*z**6))+  &
PL2*((-2.62352879841933e-2_dp*(-2._dp+r0/rz+rz/r0)*(1-2*z)*(1  &
+2*z)*(-1._dp+4*z**2))/(d2*z**5)+  &
(rz**2*z*(4.372547997365549e-3_dp-3.498038397892439e-2_dp*z**2+  &
6.996076795784878e-2_dp*z**4)+r0*rz*z*(1.311764399209665e-2_dp-  &
1.049411519367732e-1_dp*z**2+2.098823038735464e-1_dp*z**4-  &
1.7490191989462198e-1_dp*rz**2*z05m*z05**2)+  &
r0**3*rz*(-1.7490191989462198e-1_dp*z05m*z*z05**2+  &
rz**2*(-2.1862739986827755e-1_dp-3.4980383978924405e-1_dp*z+  &
5.247057596838661e-1_dp*z**2+6.996076795784882e-1_dp*z**3))+  &
r0**2*(4.372547997365549e-3_dp*z-3.498038397892439e-2_dp*z**3+  &
6.996076795784878e-2_dp*z**5+rz**2*(-1.0931369993413877e-1_dp+  &
4.37254799736555e-2_dp*z+6.121567196311771e-1_dp*z**2-  &
1.74901919894622e-1_dp*z**3-  &
6.996076795784882e-1_dp*z**4)))/(r0*rz*(0.5_dp+r0*rz-  &
z)*z**6)))+l3*((l5*(8.745095994731102e-2_dp*r0**2- 1.7490191989462198e-1_dp*&
r0*rz+ 8.745095994731102e-2_dp*rz**2)*z05m*z05**2)&
/(d2*z**5)+(l4*(8.745095994731102e-2_dp*r0**2- 1.7490191989462198e-1_dp*r0*rz+  &
8.745095994731102e-2_dp*rz**2)*z05m*z05**2)/((0.5_dp  &
+r0*rz - z)*z**5)+(PL3*(1.6397054990120816e-2_dp-  &
7.870586395257993e-2_dp*z**2+5.247057596838665e-2_dp*z**4))/z**6+  &
(PL4*(3.2794109980241632e-2_dp-1.5741172790515985e-1_dp*z**2+  &
1.049411519367733e-1_dp*z**4))/z**6-  &
(3.935293197628995e-2_dp*(r0**2-2*r0*rz+  &
rz**2)*(3.5031605365081706_dp*z**3-1.4012642146032683e1_dp*z**5+  &
r0*rz**3*z05+rz**2*(3.333333333333334e-1_dp+  &
5.833333333333334e-1_dp*z+r0*z**1.5_dp-1.4166666666666665_dp*z**2  &
+2*r0*z**2.5_dp-2.5_dp*z**3)))/(d2*rz**2*z**5)+  &
(rz**2*z*(-9.8550965711811e-3_dp+(-4.355684420256939e-3_dp+  &
9.83823299407249e-3_dp*rz**2)*z+(5.741598078557713e-2_dp+  &
1.0830348052027379e-1_dp*rz**2)*z**2+(4.562522504431183e-2_dp+  &
7.87887815380563e-2_dp*rz**2)*z**3+(-2.2994623339370577e-2_dp-  &
1.9693049505240257e-1_dp*rz**2)*z**4+1.179770495612256e-1_dp*z**5-  &
1.8381285106030598e-1_dp*z**6)+r0*rz**3*(9.8550965711811e-3_dp+  &
(-5.499412150924162e-3_dp-5.46568499670694e-3_dp*rz**2)*z+  &
(-7.162676177701516e-2_dp-9.409269952883573e-2_dp*rz**2)*z**2+  &
(-2.4200252501727637e-3_dp+4.372547997365555e-3_dp*rz**2)*z**3+  &
(1.1182504817782155e-1_dp+3.4139041413641857e-1_dp*rz**2)*z**4-  &
5.214124806214526e-2_dp*z**5+3.6762570212061196e-1_dp*z**6)+  &
r0**4*rz**2*((-9.190642553015296e-2_dp-  &
1.8381285106030592e-1_dp*z)*z**4+rz**4*(4.919116497036246e-2_dp+  &
1.0494115193677322e-1_dp*z+3.935293197628995e-2_dp*z**2)+  &
rz**2*z*(3.2794109980241597e-3_dp+2.623528798419329e-2_dp*z+  &
4.069786340969018e-1_dp*z**2+3.6762570212061183e-1_dp*z**3))+  &
r0**3*rz*(-9.190642553015296e-2_dp*z**4+  &
3.6762570212061183e-1_dp*z**6+rz**4*(4.919116497036246e-2_dp-  &
1.9676465988144975e-2_dp*z-2.7547052383402963e-1_dp*z**2-  &
1.836470158893531e-1_dp*z**3)+rz**2*z*(-5.46568499670694e-3_dp-  &
2.1862739986827755e-3_dp*z+3.7199825011797745e-1_dp*z**2-  &
3.9386099010480526e-1_dp*z**3-7.352514042412237e-1_dp*z**4))+  &
r0**2*(z**5*(9.190642553015296e-2_dp+9.190642553015296e-2_dp*z-  &
1.8381285106030592e-1_dp*z**2)+rz**6*z*(3.2794109980241597e-3_dp-  &
6.567113754595965e-2_dp*z-1.4445991908401594e-1_dp*z**2)+  &
rz**2*z**2*(9.83823299407249e-3_dp+1.6397054990120816e-2_dp*z-  &
3.807433461127085e-1_dp*z**2-1.311764399209664e-2_dp*z**3+  &
3.6762570212061183e-1_dp*z**4)+rz**4*(9.8550965711811e-3_dp-  &
4.372547997365551e-2_dp*z-1.0770028275528096e-1_dp*z**2+  &
4.671856307988114e-2_dp*z**3+7.845828241398285e-2_dp*z**4-  &
1.8381285106030598e-1_dp*z**5)))/(rz**2*(r0*rz - z)*(1+r0*rz  &
 - z)*z**6)+PL1*((r0*rz*(-6.558821996048325e-3_dp+  &
5.24705759683866e-2_dp*z**2-1.049411519367732e-1_dp*z**4)+  &
r0**2*(3.2794109980241624e-3_dp-2.62352879841933e-2_dp*z**2+  &
5.24705759683866e-2_dp*z**4)+rz**2*(3.2794109980241624e-3_dp-  &
2.62352879841933e-2_dp*z**2+  &
5.24705759683866e-2_dp*z**4))/(d2*r0*rz*z**5)+  &
(rz**2*z*(5.465684996706938e-3_dp+1.0931369993413875e-2_dp*z-  &
2.186273998682775e-2_dp*z**2-4.37254799736555e-2_dp*z**3)+  &
r0**2*(-4.372547997365551e-2_dp*z05m*z*z05**2+  &
rz**2*(-5.465684996706939e-2_dp-8.745095994731102e-2_dp*z+  &
1.3117643992096653e-1_dp*z**2+  &
1.7490191989462205e-1_dp*z**3)))/(r0*rz*z**6))+  &
PL2*((-1.311764399209665e-2_dp*(-2._dp+r0/rz+rz/r0)*(1-  &
2*z)*(1+2*z)*(-1._dp+4*z**2))/(d2*z**5)+  &
(rz**2*z*(2.1862739986827746e-3_dp-1.7490191989462196e-2_dp*z**2+  &
3.498038397892439e-2_dp*z**4)+r0*rz*z*(6.558821996048325e-3_dp-  &
5.24705759683866e-2_dp*z**2+1.049411519367732e-1_dp*z**4-  &
8.745095994731102e-2_dp*rz**2*z05m*z05**2)+  &
r0**3*rz*(-8.745095994731102e-2_dp*z05m*z*z05**2+  &
rz**2*(-1.0931369993413877e-1_dp-1.7490191989462203e-1_dp*z+  &
2.6235287984193305e-1_dp*z**2+3.498038397892441e-1_dp*z**3))+  &
r0**2*(2.1862739986827746e-3_dp*z-1.7490191989462196e-2_dp*z**3+  &
3.498038397892439e-2_dp*z**5+rz**2*(-5.465684996706939e-2_dp+  &
2.186273998682775e-2_dp*z+3.0607835981558855e-1_dp*z**2-  &
8.7450959947311e-2_dp*z**3-  &
3.498038397892441e-1_dp*z**4)))/(r0*rz*(0.5_dp+r0*rz-  &
z)*z**6))))

end function P2Der

!ccccccccccccccc

complex (dp) function P2(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: r32
  complex (dp)                :: rz, r1, r2, r3, l1, l2, l3, PL1, PL2, PL3, PL4, &
  d1, z05, z05m

  rz = sqrt(z); r1 = sqrt(1 - z); r2 = 1 + r1; r3 = 1 - r1; r1 = sqrt(z - 1)
  d1 = 1 + 2 * r1 * rz - 2 * z; l1 = Log(d1); l2 = Log(1 + r1 * rz - z)
  l3 = Log( 8 * (z - r1 * rz) ); PL1 = cli2(d1); PL2 = cli2(-d1); PL3 = cli3(d1)
  PL4 = cli3(-d1); r32 = PowList(r3/r2, 12); z05 = 0.5_dp + z; z05m = z - 0.5_dp

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
  (1.74901919894622e-2_dp*PL1**2*r1**2*rz**2*z05**2)/z**5+  &
  (6.99607679578488e-2_dp*PL2**2*r1**2*rz**2*z05**2)/z**5+  &
  (PL3**2*(9.83823299407249e-3_dp-7.870586395257993e-2_dp*z**2+  &
  1.5741172790515985e-1_dp*z**4))/z**5 +  &
  (PL4**2*(3.9352931976289964e-2_dp-3.148234558103197e-1_dp*z**2+  &
  6.296469116206394e-1_dp*z**4))/z**5+(PL4*(2.3652231770834637e-2_dp  &
  +4.263234297431412e-2_dp*z-1.7282079917655628e-1_dp*z**2-  &
  5.329825490307057e-1_dp*z**3+3.1284748837287095e-1_dp*z**4+  &
  1.4498127085337966_dp*z**5))/z**5+  &
  PL1*((6.99607679578488e-2_dp*PL2*r1**2*rz**2*z05**2)/z**5  &
  +(1.0494115193677322e-1_dp*PL3*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5+(2.098823038735464e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(r1*rz*(-7.88407725694488e-3_dp-  &
  2.9978935505327797e-2_dp*z-2.3509379518035e-3_dp*z**2+  &
  1.7295897377329488e-1_dp*z**3+  &
  2.4163545142229945e-1_dp*z**4))/z**5)+  &
  PL2*((2.098823038735464e-1_dp*PL3*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5+(4.197646077470928e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(r1*rz*(-1.5768154513889758e-2_dp-  &
  5.9957871010655595e-2_dp*z-4.701875903607e-3_dp*z**2+  &
  3.4591794754658975e-1_dp*z**3+4.832709028445989e-1_dp*z**4))/z**5)  &
  +PL3*((PL4*(3.9352931976289964e-2_dp-3.148234558103197e-1_dp*z**2+  &
  6.296469116206394e-1_dp*z**4))/z**5+(1.1826115885417319e-2_dp+  &
  2.131617148715706e-2_dp*z-8.641039958827815e-2_dp*z**2-  &
  2.6649127451535284e-1_dp*z**3+1.5642374418643548e-1_dp*z**4+  &
  7.249063542668983e-1_dp*z**5)/z**5)+  &
  l1**3*((-6.996076795784882e-2_dp*l2**2*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5-(1.7490191989462205e-2_dp*l3**2*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5-  &
  (34.295003315177963_dp*(-6.275297446905026e-5_dp*r1*rz**2*(1-  &
  4*(0.5_dp+r1*rz - z)**2)**3-  &
  2.1515305532245805e-4_dp*r1**2*rz**3*(1-4*(0.5_dp+r1*rz-  &
  z)**2)**3-1.1355300142018618e-4_dp*r1*rz**2*(1-4*(0.5_dp  &
  +r1*rz - z)**2)**3*z-1.4343537021497204e-4_dp*r1**2*rz**3*(1  &
  -4*(0.5_dp+r1*rz - z)**2)**3*z-  &
  4.303061106449161e-4_dp*r1**2*rz**2*(1-4*(0.5_dp+r1*rz-  &
  z)**2)**3*z**1.5_dp+3.8249432057325876e-4_dp*r1*rz**2*(1-  &
  4*(0.5_dp+r1*rz - z)**2)**3*z**2-  &
  2.8687074042994407e-4_dp*r1**2*rz**2*(1-4*(0.5_dp+r1*rz-  &
  z)**2)**3*z**2.5_dp-1.1724466331119023e-3_dp*r1*(1-  &
  4*(0.5_dp+r1*rz - z)**2)**3*z**3-  &
  4.019817027812235e-3_dp*r1**2*rz*(1-4*(0.5_dp+r1*rz-  &
  z)**2)**3*z**3+2.8687074042994407e-4_dp*r1*rz**2*(1-  &
  4*(0.5_dp+r1*rz - z)**2)**3*z**3+1.25e-1_dp*rz*z**4+  &
  7.5e-1_dp*r1*rz**2*z**4+1.5_dp*r1**2*rz**3*z**4+r1**3*rz**4*z**4  &
  -1.3399390092707453e-3_dp*r1*(1-4*(0.5_dp+r1*rz-  &
  z)**2)**3*z**4-8.03963405562447e-3_dp*r1**2*(1-4*(0.5_dp  &
  +r1*rz - z)**2)**3*z**4.5_dp - rz*z**5-  &
  4.5_dp*r1*rz**2*z**5-6*r1**2*rz**3*z**5-  &
  2*r1**3*rz**4*z**5+8.03963405562447e-3_dp*r1*(1-4*(0.5_dp+  &
  r1*rz - z)**2)**3*z**5+3.125_dp*rz*z**6+  &
  9.75_dp*r1*rz**2*z**6+7.499999999999998_dp*r1**2*rz**3*z**6+  &
  r1**3*rz**4*z**6-4.75_dp*rz*z**7-9*r1*rz**2*z**7-  &
  3*r1**2*rz**3*z**7+3.5_dp*rz*z**8+3*r1*rz**2*z**8-  &
  rz*z**9))/(rz*(1-4*(0.5_dp+r1*rz - z)**2)**3*z**5)  &
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
  l2*((-6.996076795784882e-2_dp*l3*r1*rz*z05m*(0.5_dp+  &
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
  l1*((-6.996076795784882e-2_dp*PL1**2*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5-(2.7984307183139516e-1_dp*PL2**2*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(PL4*r1*(5.514385531809179e-1_dp*z**3-  &
  2.2057542127236713_dp*z**5+rz**2*(2.951469898221747e-2_dp+  &
  1.9676465988144982e-2_dp*z-1.1805879592886988e-1_dp*z**2-  &
  7.870586395257993e-2_dp*z**3)))/(rz*z**5)+  &
  (PL3*r1*(2.7571927659045894e-1_dp*z**3-  &
  1.1028771063618357_dp*z**5+rz**2*(1.4757349491108736e-2_dp+  &
  9.83823299407249e-3_dp*z-5.902939796443494e-2_dp*z**2-  &
  3.9352931976289964e-2_dp*z**3)))/(rz*z**5)+  &
  l3*((3.49803839789244e-2_dp*PL1*r1**2*rz**2*z05**2)/z**5+  &
  (6.99607679578488e-2_dp*PL2*r1**2*rz**2*z05**2)/z**5+  &
  (1.0494115193677322e-1_dp*PL3*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5+(2.098823038735464e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(r1*rz*(-7.88407725694488e-3_dp-  &
  2.9978935505327797e-2_dp*z-2.3509379518035e-3_dp*z**2+  &
  1.7295897377329488e-1_dp*z**3+  &
  2.4163545142229945e-1_dp*z**4))/z**5)+  &
  l2*((6.99607679578488e-2_dp*PL1*r1**2*rz**2*z05**2)/z**5+  &
  (1.399215359156976e-1_dp*PL2*r1**2*rz**2*z05**2)/z**5+  &
  (2.098823038735464e-1_dp*PL3*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5+(4.197646077470928e-1_dp*PL4*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(r1*rz*(-1.5768154513889758e-2_dp-  &
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
  r1*rz*(1-2*z) - z+z**2))+  &
  PL2*((PL4*(-5.24705759683866e-2_dp+4.197646077470928e-1_dp*z**2-  &
  8.395292154941856e-1_dp*z**4))/z**5+(PL3*(-2.62352879841933e-2_dp+  &
  2.098823038735464e-1_dp*z**2-4.197646077470928e-1_dp*z**4))/z**5+  &
  (-1.5768154513889758e-2_dp-2.8421561982876082e-2_dp*z+  &
  1.1521386611770421e-1_dp*z**2+3.5532169935380375e-1_dp*z**3-  &
  2.0856499224858065e-1_dp*z**4-9.665418056891978e-1_dp*z**5+  &
  r1**2*((-3.676257021206119e-1_dp-7.352514042412238e-1_dp*z)*z**3+  &
  rz**2*(-1.9676465988144982e-2_dp-5.24705759683866e-2_dp*z-  &
  2.62352879841933e-2_dp*z**2)))/z**5)+  &
  PL1*((-2.7984307183139516e-1_dp*PL2*r1*rz*z05m*(0.5_dp+  &
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
  2*z)+(-1._dp+z)*z)**3)+  &
  l2*((6.99607679578488e-2_dp*l3*r1**2*rz**2*z05**2)/z**5-  &
  (1.7490191989462203e-1_dp*PL1*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5-(3.498038397892441e-1_dp*PL2*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(PL4*(-1.311764399209665e-2_dp+  &
  1.049411519367732e-1_dp*z**2-2.098823038735464e-1_dp*z**4))/z**5+  &
  (PL3*(-6.558821996048325e-3_dp+5.24705759683866e-2_dp*z**2-  &
  1.049411519367732e-1_dp*z**4))/z**5+(-3.94203862847244e-3_dp-  &
  7.10539049571902e-3_dp*z+2.880346652942605e-2_dp*z**2+  &
  8.883042483845093e-2_dp*z**3-5.2141248062145165e-2_dp*z**4-  &
  2.4163545142229945e-1_dp*z**5+r1**2*((-3.676257021206119e-1_dp-  &
  7.352514042412238e-1_dp*z)*z**3+rz**2*(-1.9676465988144982e-2_dp-  &
  5.24705759683866e-2_dp*z-2.62352879841933e-2_dp*z**2)))/z**5)+  &
  l3*((-8.745095994731102e-2_dp*PL1*r1*rz*z05m*(0.5_dp+  &
  z)**2)/z**5-(1.7490191989462203e-1_dp*PL2*r1*rz*(-0.5_dp+  &
  z)*z05**2)/z**5+(PL4*(-6.558821996048325e-3_dp+  &
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

complex (dp) function Pi3(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: rat
  complex (dp)                :: sz, r1, r2, l1, l2, l3, PL1, PL2, PL3, &
  PL4, d1, d2, z1, z2, lz, d3, r
  real (dp), parameter :: f = 1.2020569031595942_dp

  sz = sqrt(z); r1 = sqrt(1 - z); r2 = sqrt(z - 1); d1 = 1 + 2 * r2 * sz - 2 * z
  d2 = d1**2 - 1; l1 = Log(d1); PL2 = cli2(-d1); l2 = Log(d1 + 1); lz = log(z)
  l3 = Log(1 - d1); PL1 = cli2(d1); PL3 = cli3(d1); PL4 = cli3(-d1)
  z1 = 1 - z; z2 = 1 + z; d3 = r2/sz + sz/r2 - 2; r = (1 - r1)/(1 + r1)

  rat = powList( r, 7 )

  Pi3 = 3.344517218879539_dp + (1 + r)**4/4/(1 - r)**2/r * ( &
  - 1.0550389626825085_dp + 1.5804670902096216_dp * r + 2.75436533583785_dp * &
  rat(2) - 1.2823166814463995_dp * rat(3) + 0.5368852704121565_dp * rat(4) +  &
  0.22668405128072777_dp * rat(5) - 0.7432166595042755_dp * rat(6) + &
  0.6185763562892763_dp * rat(7) ) + &
  (-1.7477904178303265e1_dp*d1**4*l1**4*z1**2)/d2**4+  &
  (4*d1**2*l1**2*(-6.044677596146622_dp-  &
  2.962967887415238_dp/z)*z1)/d2**2+  &
  (7.436968010168568e-1_dp*d1**2*l1**2*z1**2)/d2**2+  &
  (2*d1*l1*(3.6179344180753232_dp+6.0958158898685975e-1_dp/z**2  &
  +7.215572574811175_dp/z)*z1)/d2+(2.6375974067062713e-1_dp*(1+  &
  r)**2)/r+2.500611735451309_dp/z**2-2.922987438453376_dp/z+  &
  (8*d1**3*l1**3*(2.190634021761004_dp-  &
  3.1101933244381694_dp/z)*z1**2*z2)/(d2**3*z)-  &
  (1.6549126661581834_dp*d1**3*l1**3*z1**3*z2)/(d2**3*z)+  &
  (3.1842906213744557e-4_dp*((5+13/z)/6+  &
  (2*d1**2*l1**2*(1 - 16 * z)*z1)/(3*d2**2*z)-  &
  (2*d1*l1*(3._dp+2*z)*z1)/(d2*z)-  &
  ((1 + 2*z)/6*((-2*l1**2*(2*l2+  &
  l3)-8*l1*(PL1+2*PL2)+6*(f+  &
  2*PL3+4*PL4))/z+2*((-(-2*l1**2*(2*l2+l3)  &
  -8*l1*(PL1+2*PL2)+6*(f+2*PL3+  &
  4*PL4)))/z**2+(-8*l1*((-2*d3*l2)/d1-  &
  (d3*l3)/d1)+6*((2*d3*PL1)/d1+(4*d3*PL2)/d1)-  &
  (4*d3*l1*(2*l2+l3))/d1-(8*d3*(PL1+2*PL2))/d1-  &
  2*((-d3)/(1-d1)+(2*d3)/(1 +  &
  d1))*l1**2)/z)*z*z1))/z)**3)/z+  &
  1.959018998158646_dp*((5+13/z)/6+  &
  (2*d1**2*l1**2*(1 - 16 * z)*z1)/(3*d2**2*z)-  &
  (2*d1*l1*(3._dp+2*z)*z1)/(d2*z)-((1 +  &
  2*z)*((-2*l1**2*(2*l2+l3)-8*l1*(PL1+2*PL2)+  &
  6*(f+2*PL3+4*PL4))/z+  &
  2*((-(-2*l1**2*(2*l2+l3)-8*l1*(PL1+  &
  2*PL2)+6*(f+2*PL3+  &
  4*PL4)))/z**2+(-8*l1*((-2*d3*l2)/d1-  &
  (d3*l3)/d1)+6*((2*d3*PL1)/d1+(4*d3*PL2)/d1)-  &
  (4*d3*l1*(2*l2+l3))/d1-(8*d3*(PL1+2*PL2))/d1-  &
  2*((-d3)/(1-d1)+(2*d3)/(1+  &
  d1))*l1**2)/z)*z*z1))/(6*z))+  &
  (23*d1*l1*((5+13/z)/6+  &
  (2*d1**2*l1**2*(1 - 16 * z)*z1)/(3*d2**2*z)-  &
  (2*d1*l1*(3._dp+2*z)*z1)/(d2*z)-((1+  &
  2*z)*((-2*l1**2*(2*l2+l3)-8*l1*(PL1+2*PL2)+  &
  6*(f+2*PL3+4*PL4))/z+  &
  2*((-(-2*l1**2*(2*l2+l3)-8*l1*(PL1+  &
  2*PL2)+6*(f+2*PL3+  &
  4*PL4)))/z**2+(-8*l1*((-2*d3*l2)/d1-  &
  (d3*l3)/d1)+6*((2*d3*PL1)/d1+(4*d3*PL2)/d1)-  &
  (4*d3*l1*(2*l2+l3))/d1-(8*d3*(PL1+2*PL2))/d1-  &
  2*((-d3)/(1-d1)+(2*d3)/(1+  &
  d1))*l1**2)/z)*z*z1))/(6*z)))/(2.7e1_dp*d2*z)+  &
  3.14573896293457e-2_dp*((5+13/z)/6+  &
  (2*d1**2*l1**2*(1 - 16 * z)*z1)/(3*d2**2*z)-  &
  (2*d1*l1*(3._dp+2*z)*z1)/(d2*z)-((1+  &
  2*z)*((-2*l1**2*(2*l2+l3)-8*l1*(PL1+2*PL2)+  &
  6*(f+2*PL3+4*PL4))/z+  &
  2*((-(-2*l1**2*(2*l2+l3)-8*l1*(PL1+  &
  2*PL2)+6*(f+2*PL3+  &
  4*PL4)))/z**2+(-8*l1*((-2*d3*l2)/d1-  &
  (d3*l3)/d1)+6*((2*d3*PL1)/d1+(4*d3*PL2)/d1)-  &
  (4*d3*l1*(2*l2+l3))/d1-(8*d3*(PL1+2*PL2))/d1-  &
  2*((-d3)/(1-d1)+(2*d3)/(1+  &
  d1))*l1**2)/z)*z*z1))/(6*z))**2

end function Pi3

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
