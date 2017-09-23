
module SigmaClass
  use AnomDimClass;  use RunningClass;  use ElectroWeakClass; use Chaplin
  use Constants, only: dp, Pi, Pi2, Zeta3;  implicit none;  private

  real (dp), dimension(2)               :: EWfact
  public                                :: setEWfact, VacPol0, VacPol0Der, &
  VacPol2, VacPol2Der, VacPol3, VacPol1, PiCoef, VacPol1Der, PiCoefDer1,   &
  PiCoefDer2, PiCoefDer3

!ccccccccccccccc

  type, public :: Sigma
   private
   integer                        :: nf
   real (dp), dimension(0:3,4)    :: RhadCoef
   real (dp)                      :: m, mH
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
    InitSigma%mH = run%scales('mH')

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
    real (dp), dimension( 0:min(order,4) )   :: rQ
    integer                                  :: n

    rQ = 1 ; n = min(order,4) - 1
    if (order > 1) rQ(2:) = matmul( PowList0( log(mu/Q), n ), self%RhadCoef(:n,2:n) )

    Rhad = dot_product( PowList0( self%run%alphaQCD(mu)/Pi, n + 1 ), rQ )

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

complex (dp) function VacPol0(z)
  complex (dp), intent(in) :: z

  VacPol0 = 3 * (  20._dp/9 + 4/z/3 - 4 * (1 - z) * (1 + 2 * z)/3/z * Gz( u(z) )  )/16/Pi2

end function VacPol0

!ccccccccccccccc

complex (dp) function VacPol0Der(i, z)
  integer     , intent(in) :: i
  complex (dp), intent(in) :: z
  complex (dp) :: uu, uuDer, uuDer2, uuDer3

  VacPol0Der = 0; if (i > 0) uu = u(z); if (i > 0) uuDer = uDer(1,uu)
  if (i > 1) uuDer2 = uDer(2,uu);  if (i > 2) uuDer3 = uDer(3,uu)

  if ( i == 0 ) then

    VacPol0Der = VacPol0(z)

  else if (i == 1) then

    VacPol0Der = 4 * ( (1 + 2 * z**2) * Gz(uu) - 1 + z * (2 * z**2 - 1 - z) * &
    Gzder(1,uu) * uuDer )/3/z**2

  else if (i == 2) then

    VacPol0Der = 4 * (  2 - 2 * Gz(uu) + (z - 1) * z**2 * (1 + 2 * z) * uuDer**2 * &
    Gzder(2,uu) + z * Gzder(1,uu) * ( 2 * (1 + 2 * z**2) * uuDer + (z - 1) * &
    z * (1 + 2 * z) * uuDer2 )  )/3/z**3

  else if (i == 3) then

    VacPol0Der = 4 * (   6 * Gz(uu) - 6 + z * (   z * ( uuDer3 * (z - 1) * z * &
    (1 + 2 * z) + uuDer2 * (3 + 6 * z**2) ) - 6 * uuDer  ) * GZDer(1,uu) + &
    uuDer * z**2 * (  3 * ( uuDer2 * (z - 1) * z * (1 + 2 * z) + uuDer * &
    (1 + 2 * z**2) ) * GZDer(2, uu) + uuDer**2 * (z - 1) * z * (1 + 2 * z) *&
     GZDer(3, uu)  )   )/3/z**4

  end if

  if (i > 0) VacPol0Der = 3 * VacPol0Der/16/Pi2

end function VacPol0Der

!ccccccccccccccc

  function PiCoef(z, n) result(res)
    complex (dp), intent(in)           :: z
    integer     , intent(in)           :: n
    real (dp), dimension( 0:min(n,3) ) :: res

    if (n >= 0) res(0) = ImagPart( z * VacPol0(z) )
    if (n >= 1) res(1) = 4 * ImagPart( z * VacPol1(z) )/3
    if (n >= 2) res(2) = ImagPart( z * VacPol2(z)  )
    if (n >= 3) res(3) = ImagPart( z * VacPol3(z) )

  end function PiCoef

!ccccccccccccccc

  function PiCoefDer1(z, n) result(res)
    complex (dp), intent(in)           :: z
    integer     , intent(in)           :: n
    real (dp), dimension( 0:min(n,3) - 1 ) :: res

    if (n >= 1) res(0) = ImagPart( z**2 * VacPol0Der(1,z) )
    if (n >= 2) res(1) = 4 * ImagPart( z**2 * VacPol1Der(1,z) )/3
    if (n >= 3) res(2) = ImagPart( z**2 * VacPol2Der(z)  )

 end function PiCoefDer1

!ccccccccccccccc

 function PiCoefDer2(z, n) result(res)
   complex (dp), intent(in)           :: z
   integer     , intent(in)           :: n
   real (dp), dimension( 0:min(n,3) - 2 ) :: res

   if (n >= 2) res(0) = ImagPart( z**3 * VacPol0Der(2,z) )
   if (n >= 3) res(1) = 4 * ImagPart( z**3 * VacPol1Der(2,z) )/3

end function PiCoefDer2

!ccccccccccccccc

real (dp) function PiCoefDer3(z)
  complex (dp), intent(in) :: z

  PiCoefDer3 = ImagPart( z**4 * VacPol0Der(3,z) )

end function PiCoefDer3

!ccccccccccccccc

complex (dp) function VacPol1Der(i, z)
  integer     , intent(in) :: i
  complex (dp), intent(in) :: z
  complex (dp)             :: uu, uuDer, uuDer2, uuDer3

  VacPol1Der = 0; if (i > 0) uu = u(z); if (i > 0) uuDer = uDer(1,uu)
  if (i > 0) uuDer2 = uDer(2,uu);  if (i > 1) uuDer3 = uDer(3,uu)

  if ( i == 0 ) then

    VacPol1Der = VacPol1(z)

  else if (i == 1) then

    VacPol1Der = (z*(-1 + 16*z**2)*Gz(uu)**2 + 2*z*Gz(uu)*(9 + 6*z**2 + &
    uuDer*(-1 + z)*z*(-1 + 16*z)*GZDer(1,uu)) - 2*Iz(U(z)) +  &
    z*(-13 + 2*uuDer2*(-1 + z)*z*(1 + 2*z)*IzDer(1,uu) + uuDer*(3*IzDer(1,uu) &
    + 2*(-1 + z)*z*((9 + 6*z)*GZDer(1,uu) + uuDer*(1 + 2*z)*IzDer(2,uu)))))/ &
    (6*z**3)

  else if (i == 2) then

    VacPol1Der = (2*z*Gz(uu)**2 + 2*z*Gz(uu)*(-18 + z*(uuDer2*(-1 + z)*z*(-1 + 16*z)  &
    + uuDer*(-2 + 32*z**2))*GZDer(1,uu) + uuDer**2*(-1 + z)*z**2*(-1 + 16*z)* &
    GZDer(2,uu)) + 6*Iz(U(z)) + z*(26 + 6*z*(uuDer2*z*(-3 + z + 2*z**2) + &
    uuDer*(6 + 4*z**2))*GZDer(1,uu) + 2*uuDer**2*(-1 + z)*z**2*(-1 + 16*z)* &
    GZDer(1,uu)**2 + (-8*uuDer + z*(2*uuDer3*(-1 + z)*z*(1 + 2*z) + &
    uuDer2*(5 + 4*z**2)))*IzDer(1,uu) + uuDer*z*(6*uuDer2*(-1 + z)*z*(1 + 2*z) &
    *IzDer(2,uu) + uuDer*(6*(-1 + z)*z*(3 + 2*z)*GZDer(2,uu) + &
    (5 + 4*z**2)*IzDer(2,uu)) + 2*uuDer**2*(-1 + z)*z*(1 + 2*z)*IzDer(3,uu))))/(6*z**4)

  end if

  if (i > 0) VacPol1Der = 3 * VacPol1Der/16/Pi2

end function VacPol1Der

!ccccccccccccccc

complex (dp) function VacPol1(z)
  complex (dp), intent(in) :: z
  complex (dp)             :: uu, Gzuu, IIz

  uu = u(z); Gzuu = Gz(uu); IIz = Iz(uu)

  VacPol1 = 3 * (  5._dp/6 + 13/z/6 - (1 - z) * (3 + 2 * z) * Gzuu/z + &
  (1 - z) * (1 - 16 * z) * Gzuu**2/z/6 - (1 + 2 * z)/z/6 * ( IIz/z + &
  2 * (1 - z) * ( IzDer(1,uu) * uDer(1,uu) - IIz/z )  )   )/16/Pi2

end function VacPol1

!ccccccccccccccc

complex (dp) function VacPol2Der(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: r32
  complex (dp)                :: rz, r1, r2, r3, Pi1z, GGzDer, uu, GGz, Pi1zDer

  rz = sqrt(z); r1 = sqrt(1 - z); r2 = 1 + r1; r3 = 1 - r1; uu = u(z)
  Pi1z = VacPol1(z); GGz = Gz(uu); GGzDer = GzDer(1,uu) * uDer(1,uu)
  Pi1zDer = VacPol1Der(1,z);  r32 = PowList(r3/r2, 12)

  VacPol2Der = (2 + 1/r1)/(r2 - z)**2 * ( 0.688472275231591_dp + &
  0.1891443911186843_dp * r32(1) - 0.030243753917029492_dp * r32(2) + &
  0.0034204077615457664_dp * r32(3) + 0.006118445792320977_dp * r32(4) &
  + 0.0016455445946186606_dp * r32(5) + 0.001252836593493741_dp * r32(6) &
  + 0.00046831306017028774_dp * r32(7) + 0.0018087001167161186_dp * r32(8) &
  + 0.00048205327014239514_dp * r32(9) - 0.0024098794332697245_dp * r32(10) &
  - 0.0001825139597445611_dp * r32(11) + 0.0013086431124221074_dp * r32(12) ) &
  + 2/(r2 - z)/r1/r2 * ( 0.09457219555934215_dp + 0.06432844164231266_dp * r32(1) - &
  0.025113142274710842_dp * r32(2) + 0.017367503226960603_dp * r32(3) + &
  0.016350753071188606_dp * r32(4) + 0.007872371267027875_dp * r32(5) + &
  0.00539760549107723_dp * r32(6) + 0.008873896177460481_dp * r32(7) + &
  0.009404040182505253_dp * r32(8) - 0.009880157450707843_dp * r32(9) - &
  0.013053223944943708_dp * r32(10) + 0.006848031895937559_dp * r32(11) + &
  0.007851858674532645_dp * r32(12) ) + (-3.4800925587898566_dp - &
  0.7128763906962368_dp * z + GGz * (  2.859449347250649_dp + 0.2375443345172152_dp &
  * z + 1.342733077226607_dp * z**3 + GGz * ( 0.6206432115392082_dp - &
  0.3557204404643448_dp * z - 0.3557204404643448_dp * z**3 - &
  0.9255650349877445_dp * z**4 + z * (-0.5358594267996556_dp + &
  0.5358594267996556_dp * z**2) * GGz )  ) + 14.513280990412627_dp * z**3 * &
  GGz * Pi1z - 6.814814814814815_dp * z * Pi1z**2 + z * (-1.4297246736253244_dp &
  + z * (-0.23754433451721538_dp + z * (0.3245359309159326_dp + &
  1.342733077226607_dp * z)) + GGz * (-0.6206432115392082_dp + z * &
  (0.7114408809286896_dp + z * (1.5462082465269527_dp + &
  (-0.7114408809286896_dp - 0.9255650349877445_dp * z) * z)) + &
  1.607578280398967_dp * (1 - z)**2 * z * GGz) + z**2 * &
  (-14.513280990412621_dp + 14.513280990412621_dp * z) * Pi1z) * GGzDer + &
   z**2 * (-4.76969262963348_dp * z + z * (-14.513280990412627_dp + &
   14.513280990412627_dp * z) * GGz + 13.62962962962963_dp * Pi1z) * Pi1zDer)/z**3

end function VacPol2Der

!ccccccccccccccc

complex (dp) function VacPol2(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: r32
  complex (dp)                :: z1, r1, r2, r3, GGz, Pi1z

  z1 = 1 - z; r1 = sqrt(z1); r2 = 1 + r1; r3 = 1 - r1; Pi1z =  VacPol1(z)
  r32 = PowList(r3/r2, 12); GGz = Gz( u(z) )

  VacPol2 = - 0.8705931079221861_dp - 4.76969262963348_dp * Pi1z - &
  14.513280990412627_dp * GGz * Pi1z * z1 - GGz * (1.342733077226607_dp &
  + 1.4297246736253244_dp/z**2 + 1.6672690081425396_dp/z) * z1 + &
  2 * ( 0.688472275231591_dp + &
  0.1891443911186843_dp * r32(1) - 0.030243753917029492_dp * r32(2) + &
  0.0034204077615457664_dp * r32(3) + 0.006118445792320977_dp * r32(4) + &
  0.0016455445946186606_dp * r32(5) + 0.001252836593493741_dp * r32(6) + &
  0.00046831306017028774_dp * r32(7) + 0.0018087001167161186_dp * r32(8) + &
  0.00048205327014239514_dp * r32(9) - 0.0024098794332697245_dp * r32(10) - &
  0.0001825139597445611_dp * r32(11) + 0.0013086431124221074_dp * r32(12) ) &
  /(r2 - z) + 1.7400462793949283_dp/z**2 + 0.7128763906962368_dp/z + &
  (6.814814814814815_dp * VacPol1(z)**2)/z + 0.5358594267996556_dp * GGz**3 * &
  z1**2/z + GGz**2 * (0.818502957958217_dp - 0.3103216057696041_dp/z) * &
  (1 - z**2)/z - 0.4627825174938722_dp * GGz**2 * z1**2 * (1 + z)/z

end function VacPol2

!ccccccccccccccc

complex (dp) function VacPol3(z)
  complex (dp), intent(in)    :: z
  complex (dp), dimension(12) :: rat
  complex (dp)                :: z1, z2, r1, r2, r, GGz, Pi1z

  z1 = 1 - z; r1 = sqrt(z1); r2 = sqrt(z - 1);  r = (1 - r1)/(1 + r1)
  Pi1z =  VacPol1(z); GGz = Gz( u(z) ); z2 = 1 + z

  rat = powList( r, 7 )

  VacPol3 = 1.0550389626825085_dp/z + (1 + r)**4/4/(1 - r)**2/r * ( &
  - 1.0550389626825085_dp + 1.5804670902096216_dp * r + 2.75436533583785_dp * &
  rat(2) - 1.2823166814463995_dp * rat(3) + 0.5368852704121565_dp * rat(4) +  &
  0.22668405128072777_dp * rat(5) - 0.7432166595042755_dp * rat(6) + &
  0.6185763562892763_dp * rat(7) ) + 3.3445172188795387_dp + &
  103.11862680556928_dp * Pi1z + 87.16048298942633_dp * Pi1z**2 - &
  GGz**2 * (6.044677596146621_dp + 2.962967887415238_dp/z) * z1 + &
  GGz * (3.6179344180753232_dp + 0.6095815889868598_dp/z**2 + &
  7.215572574811175_dp/z) * z1 + 0.1859242002542142_dp * GGz**2 * z1**2 - &
  1.0923690111439541_dp * GGz**4 * z1**2 + 2.500611735451309_dp/z**2 - &
  2.922987438453376_dp/z + 22.41984209630175_dp * GGz * Pi1z/z + &
  46.44170096021948_dp * Pi1z**3/z + GGz**3 * (2.190634021761004_dp - &
  3.11019332443817_dp/z) * z1**2 * z2/z - &
  0.20686408326977293_dp * GGz**3 * z1**3 * z2/z

end function VacPol3

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
