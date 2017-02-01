
module AlphaClass
  use AnomDimClass;  use Constants, only: dp, Pi, d1mach; implicit none
  private

  interface alphaGeneric
    module procedure :: alphaGenericReal, alphaGenericComplex
  end interface alphaGeneric

!ccccccccccccccc

  type, public :: Alpha
    private
    character (len = 5)             :: str
    integer                         :: order, run, n
    type (AnomDim) , dimension(3:6) :: andim
    real (dp)      , dimension(3:6) :: alphaRef, muRef
    real (dp)                       :: mT, mB, mC
    logical                         :: QmT, QmB, QmC, QmuT, QmuC, QmuB

   contains

   generic                          :: alphaQCD => alphaQCDReal, alphaQCDComplex
   procedure, pass(self)            :: alphaQCDReal, alphaQCDComplex
   procedure, pass(self), public    :: scales, orders, adim, scheme, SetAlpha, SetMTop, &
                                       SetMBottom, SetMCharm
  end type Alpha

!ccccccccccccccc

  interface Alpha
    module procedure  InitAlpha
  end interface Alpha
  
  contains

!ccccccccccccccc

  type (Alpha) function InitAlpha(str, order, run, G4, mZ, amZ, mT, muT, mB, muB, mC, muC)
    integer  , intent(in)                 :: order
    integer  , intent(in)                 :: run
    real (dp), intent(in), optional       :: amZ
    real (dp), intent(in), optional       :: mZ, mT, muT, mB, muB, mC, muC
    real (dp), intent(in), dimension(3:6) :: G4
    character (len = *), intent(in):: str

    InitAlpha%muRef = 0; InitAlpha%run = run; InitAlpha%n = order - 1
    InitAlpha%muRef(5) = mZ

! initialising all Anomalous Dimension Objects

    InitAlpha%andim(5) = AnomDim( str, 5, G4(5) )

    if ( present(mB) ) then; InitAlpha%mB = mB; InitAlpha%QmB = .true.; end if
    if ( present(mC) ) then; InitAlpha%mC = mC; InitAlpha%QmC = .true.; end if
    if ( present(mT) ) then; InitAlpha%mT = mT; InitAlpha%QmT = .true.; end if

    if ( present(muB) ) then
      InitAlpha%muRef(4) = muB; InitAlpha%andim(4) = AnomDim( str, 4, G4(4) )
      InitAlpha%QmuB = .true.
    end if

    if ( present(muT) ) then
      InitAlpha%muRef(6) = muT; InitAlpha%andim(6) = AnomDim( str, 6, G4(6) )
      InitAlpha%QmuT = .true.
    end if
    
    if ( present(muC) ) then
      InitAlpha%muRef(3) = muC; InitAlpha%andim(3) = AnomDim( str, 3, G4(3) )
      InitAlpha%QmuC = .true.
    end if

    InitAlpha%str = str  ;  InitAlpha%order = order ;  InitAlpha%alphaRef = 0

    if ( .not. present(mZ) .or. mZ <= d1mach(1) ) return

    call InitAlpha%setAlpha(amZ)

  end function InitAlpha

!ccccccccccccccc

  subroutine SetAlpha(self, amZ)
    class (Alpha), intent(inout) :: self
    real (dp)    , intent(in)    :: amZ

    self%alphaRef(5) = amZ

    if (self%QmuT .and. self%QmT) call self%SetMTop(self%mT, self%muRef(6) )

!   Running from mZ to muB, with nf = 5 flavors, and matching

    if (self%QmuB .and. self%QmB) then
      call self%SetMBottom( self%mB, self%muRef(4) )

    end if
  
  end subroutine SetAlpha

!ccccccccccccccc

  subroutine SetMTop(self, mT, muT)
    class (Alpha), intent(inout) :: self
    real (dp)    , intent(in   ) :: mT, muT

    real (dp)                             :: aS
    real (dp)        , dimension(0:3,0:3) :: tab
    real (dp)      , dimension(0:self%n) :: lgList, alphaList

    self%mT = mT; self%muRef(6) = muT; alphaList(0) = 1; lgList(0) = 1

    tab = self%andim(6)%alphaMatchingInverse(6)
    aS = alphaGeneric(self%run, self%andim(5)%betaQCD('beta'), self%muRef(5), self%alphaRef(5), muT)
    lgList(1:) = PowList(2 * log(muT/mT), self%n); alphaList(1:) = PowList(aS/Pi, self%n)
    self%alphaRef(6) = aS * dot_product( alphaList, matmul(tab(:self%n,:self%n), lgList) )

  end subroutine SetMTop

!ccccccccccccccc

  subroutine SetMBottom(self, mB, muB)
    class (Alpha), intent(inout) :: self
    real (dp)    , intent(in   ) :: mB, muB

    real (dp)                             :: aS
    real (dp)        , dimension(0:3,0:3) :: tab
    real (dp)      , dimension(0:self%n) :: lgList, alphaList

    self%mB = mB; self%muRef(4) = muB; alphaList(0) = 1; lgList(0) = 1

    aS  = alphaGeneric(self%run, self%andim(5)%betaQCD('beta'), self%muRef(5), self%alphaRef(5), muB )
    tab = self%andim(5)%alphaMatching(5)
    lgList(1:) = PowList(2 * log(muB/mB), self%n); alphaList(1:) = PowList(aS/Pi, self%n)
    self%alphaRef(4) = aS * dot_product( alphaList, matmul(tab(:self%n,:self%n), lgList) )

!   Running from mB to muC, with nf = 4 flavors, and matching
 
    if (self%QmuC .and. self%QmC) call self%SetMCharm( self%mC, self%muRef(3) )

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, mC, muC)
    class (Alpha), intent(inout) :: self
    real (dp)    , intent(in   ) :: mC, muC

    real (dp)                             :: aS
    real (dp)        , dimension(0:3,0:3) :: tab
    real (dp)       , dimension(0:self%n) :: lgList, alphaList

    self%mC = mC; self%muRef(3) = muC; alphaList(0) = 1; lgList(0) = 1

    tab = self%andim(4)%alphaMatching(4)
    aS  = alphaGeneric(self%run, self%andim(4)%betaQCD('beta'), self%muRef(4), self%alphaRef(4), muC)
    lgList(1:) = PowList(2 * log(muC/mC), self%n); alphaList(1:) = PowList(aS/Pi, self%n)
    self%alphaRef(3) = aS * dot_product(  alphaList, matmul( tab(:self%n,:self%n), lgList )  )

  end subroutine SetMCharm

!ccccccccccccccc

   pure character (len = 5) function str(self)
    class (Alpha), intent(in) :: self
    str = self%str
   end function str

!ccccccccccccccc

   pure real (dp) function alphaQCDReal(self, nf, mu)
    class (Alpha)  , intent(in) :: self
    real (dp)      , intent(in) :: mu
    integer        , intent(in) :: nf

    integer                     :: n

    n = nf; if (nf < 4) n = 3

    if ( self%muRef(5) <= d1mach(1) ) then
      alphaQCDReal = Pi
    else
      alphaQCDReal = alphaGeneric( self%run, self%andim(n)%betaQCD('beta'), self%muRef(n), &
                                   self%alphaRef(n), mu )
    end if

   end function alphaQCDReal

!ccccccccccccccc

   pure complex (dp) function alphaQCDComplex(self, nf, mu)
    class   (Alpha), intent(in) :: self
    complex (dp)   , intent(in) :: mu
    integer        , intent(in) :: nf

    integer                        :: n

    n = nf; if (nf < 4) n = 3

    if ( self%muRef(5) <= d1mach(1) ) then
      alphaQCDComplex = Pi
    else
      alphaQCDComplex = alphaGeneric( self%run, self%andim(n)%betaQCD('beta'), &
                                             self%muRef(n), self%alphaRef(n), mu )
    end if

   end function alphaQCDComplex

!ccccccccccccccc

  pure real (dp) function scales(self, str)
    class (Alpha)      , intent(in) :: self
    character (len = *), intent(in) :: str
    
    scales = 0

    if ( str(:3) == 'muC' ) scales = self%muRef(3)
    if ( str(:3) == 'muB' ) scales = self%muRef(4)
    if ( str(:2) == 'mZ'  ) scales = self%muRef(5)
    if ( str(:3) == 'muT' ) scales = self%muRef(6)
    if ( str(:3) == 'amZ' ) scales = self%alphaRef(5)
    if ( str(:2) == 'mC'  ) scales = self%mC
    if ( str(:2) == 'mB'  ) scales = self%mB
    if ( str(:2) == 'mT'  ) scales = self%mT

  end function scales

!ccccccccccccccc

  pure integer function orders(self, str)
    class (Alpha)      , intent(in) :: self
    character (len = *), intent(in) :: str

    orders = 0

    if ( str(:5) == 'order' ) orders = self%order
    if ( str(:3) == 'run'   ) orders = self%run

  end function orders

!ccccccccccccccc

 pure real (dp) function alphaGenericReal(order, beta, mZ, amZ, mu)
    integer                  , intent(in) :: order
    real (dp)                , intent(in) :: mZ, amZ, mu
    real (dp), dimension(0:3), intent(in) :: beta
    real (dp)                             :: L, arg, LG

    if ( max( amZ, mZ, mu) <= d1mach(1) ) then
      alphaGenericReal = 0; return
    end if

    L = log(mu/mZ);  arg = amZ * L * beta(0)/2/Pi;  LG  = log(1 + arg)
    alphaGenericReal = 1

    if (order > 0) alphaGenericReal = alphaGenericReal + arg
    if (order > 1) alphaGenericReal = alphaGenericReal + beta(1)/beta(0) * amZ * LG/4/Pi

    if (order > 2) alphaGenericReal = alphaGenericReal + 1/( 16 * Pi**2 + &
     8 * Pi * amZ * L * beta(0) ) * (  amZ**2 * LG * beta(1)**2/beta(0)**2 - &
     ( amZ**3 * L * beta(1)**2/2/Pi )/beta(0) + amZ**3 * L * beta(2)/2/Pi  )

    if (order > 3) alphaGenericReal = alphaGenericReal + ( -amZ**3 * LG**2 * beta(1)**3/beta(0)**3 + &
     amZ**5 * L**2 * beta(1)**3/beta(0)/4/Pi**2 - amZ**5 * L**2 * beta(1) * beta(2)/2/Pi**2 + &
     2 * amZ**3 * LG * beta(1) * beta(2)/beta(0)**2 - amZ**4 * L * beta(1) * beta(2)/beta(0)/Pi + &
     amZ**4 * L * beta(3)/Pi + amZ**5 * L**2 * beta(0) * beta(3)/4/Pi**2 )/32/Pi/( 2 * Pi + amZ * L * beta(0) )**2

    alphaGenericReal = amZ/alphaGenericReal

  end function alphaGenericReal

!ccccccccccccccc

  pure complex (dp) function alphaGenericComplex(order, beta, mZ, amZ, mu)
    integer                  , intent(in) :: order
    real    (dp)             , intent(in) :: mZ, amZ
    complex (dp)             , intent(in) :: mu
    real (dp), dimension(0:3), intent(in) :: beta
    complex (dp)                          :: L, arg, LG

    if ( max( amZ, mZ) <= d1mach(1) ) then
      alphaGenericComplex = 0; return
    end if

    L = log(mu/mZ);  arg = amZ * L * beta(0)/2/Pi;  LG  = log(1 + arg)
    alphaGenericComplex = 1

    if (order > 0) alphaGenericComplex = alphaGenericComplex + arg
    if (order > 1) alphaGenericComplex = alphaGenericComplex + beta(1)/beta(0) * amZ * LG/4/Pi

    if (order > 2) alphaGenericComplex = alphaGenericComplex + 1/( 16 * Pi**2 + &
     8 * Pi * amZ * L * beta(0) ) * (  amZ**2 * LG * beta(1)**2/beta(0)**2 - &
     ( amZ**3 * L * beta(1)**2/2/Pi )/beta(0) + amZ**3 * L * beta(2)/2/Pi  )

    if (order > 3) alphaGenericComplex = alphaGenericComplex + ( -amZ**3 * LG**2 * beta(1)**3/beta(0)**3 + &
     amZ**5 * L**2 * beta(1)**3/beta(0)/4/Pi**2 - amZ**5 * L**2 * beta(1) * beta(2)/2/Pi**2 + &
     2 * amZ**3 * LG * beta(1) * beta(2)/beta(0)**2 - amZ**4 * L * beta(1) * beta(2)/beta(0)/Pi + &
     amZ**4 * L * beta(3)/Pi + amZ**5 * L**2 * beta(0) * beta(3)/4/Pi**2 )/32/Pi/( 2 * Pi + amZ * L * beta(0) )**2

    alphaGenericComplex = amZ/alphaGenericComplex

  end function alphaGenericComplex

!ccccccccccccccc

  pure type (AnomDim) function adim(self, nf)
    class (Alpha), intent(in) :: self
    integer      , intent(in) :: nf

    integer                   :: n

    n = nf; if (nf > 6) n = 6; if (nf < 4) n = 4;  adim = self%andim(n)

  end function adim

!ccccccccccccccc

   pure character (len = 5) function scheme(self)
    class (Alpha), intent(in) :: self
    scheme = self%str
   end function scheme

!ccccccccccccccc

end module AlphaClass
