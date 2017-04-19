
module SigmaClass
  use AnomDimClass;  use RunningClass;  use ElectroWeakClass; use Constants, only: dp, Pi
  implicit none
  private
  real (dp), dimension(2)               :: EWfact
  public                                :: setEWfact

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

end module SigmaClass
