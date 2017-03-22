
module AlphaClass
  use AnomDimClass;  use Constants, only: dp, Pi, Pi2, d1mach; implicit none
  private

  real (dp), parameter :: fourPi = 4*Pi

  interface PiBeta
    module procedure :: PiBetaReal, PiBetaComplex
  end interface PiBeta

  interface Iter
    module procedure :: IterReal, IterComplex
  end interface Iter

!ccccccccccccccc

  type, public :: Alpha
    private
    character (len = 9)             :: str, method
    integer                         :: order, run, n
    type (AnomDim) , dimension(3:6) :: andim
    real (dp)      , dimension(3:6) :: alphaRef, muRef
    real (dp)                       :: mT, mB, mC
    logical                         :: QmT, QmB, QmC, QmuT, QmuC, QmuB

   contains

   generic                          :: alphaQCD     => alphaQCDReal, alphaQCDComplex
   generic, public                  :: alphaGeneric => alphaGenericReal, alphaGenericComplex

   procedure, pass(self)            :: alphaQCDReal, alphaQCDComplex, alphaGenericReal, alphaGenericComplex
   procedure, pass(self), public    :: scales, orders, scheme, SetAlpha, SetMTop, &
   SetMBottom, SetMCharm, adim


  end type Alpha

!ccccccccccccccc

  interface Alpha
    module procedure  InitAlpha
  end interface Alpha

  contains

!ccccccccccccccc

  type (Alpha) function InitAlpha(str, order, run, G4, mZ, amZ, mT, muT,&
    mB, muB, mC, muC, method)
    character (len = *), optional, intent(in) :: method
    character (len = *)      , intent(in) :: str
    integer                  , intent(in) :: order
    integer                  , intent(in) :: run
    real (dp), optional      , intent(in) :: amZ
    real (dp), optional      , intent(in) :: mZ, mT, muT, mB, muB, mC, muC
    real (dp), dimension(3:6), intent(in) :: G4

    InitAlpha%muRef = 0; InitAlpha%run = run; InitAlpha%n = order - 1
    InitAlpha%muRef(5) = mZ

    if ( present(method) ) then
      InitAlpha%method = method
    else
      InitAlpha%method = 'analytic'
    end if

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
    aS  = self%alphaGeneric(self%andim(5), self%muRef(5), self%alphaRef(5), muT)
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

    aS  = self%alphaGeneric(self%andim(5), self%muRef(5), self%alphaRef(5), muB )
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
    aS  = self%alphaGeneric(self%andim(4), self%muRef(4), self%alphaRef(4), muC)
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
      alphaQCDReal = self%alphaGeneric(self%andim(n), self%muRef(n), self%alphaRef(n), mu )
    end if

   end function alphaQCDReal

!ccccccccccccccc

   pure complex (dp) function alphaQCDComplex(self, nf, mu)
    class   (Alpha), intent(in) :: self
    complex (dp)   , intent(in) :: mu
    integer        , intent(in) :: nf
    integer                     :: n

    n = nf; if (nf < 4) n = 3

    if ( self%muRef(5) <= d1mach(1) ) then
      alphaQCDComplex = Pi
    else
      alphaQCDComplex = self%alphaGeneric(self%andim(n), self%muRef(n), self%alphaRef(n), mu )
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

  pure real (dp) function PiBetaReal(beta, alpha)
    real (dp), dimension(:), intent(in) :: beta
    real (dp)              , intent(in) :: alpha
    real (dp)                           :: a

    PiBetaReal = 0; if ( size(beta) < 1 ) return;  a = alpha/4/Pi

    PiBetaReal = - 2 * alpha * dot_product(  beta, powList( a, size(beta) )  )

  end function PiBetaReal

!ccccccccccccccc

  pure complex (dp) function PiBetaComplex(beta, alpha)
    real (dp), dimension(:), intent(in) :: beta
    complex (dp)           , intent(in) :: alpha
    complex (dp)                        :: a

    PiBetaComplex = 0; if ( size(beta) < 1 ) return;  a = alpha/4/Pi

    PiBetaComplex = - 2 * dcmplx(0,1) * alpha * &
    dot_product(  beta, powList( a, size(beta) )  )

  end function PiBetaComplex

!ccccccccccccccc

 pure real (dp) function alphaGenericReal(self, adim, mZ, amZ, mu)
   class (Alpha)            , intent(in) :: self
   type (AnomDim)           , intent(in) :: adim
   real (dp)                , intent(in) :: mZ, amZ, mu
   real (dp), dimension(0:4)             :: beta
   real (dp)                             :: L, arg, LG, h, k1, k2, k3, k4, aLLInv, corr
   integer                               :: n, i, ord
   real (dp), dimension(0:4)             :: bCoef, cCoef
   real (dp), dimension(0:20)            :: ExpandCoeff

    if ( max( amZ, mZ, mu ) <= d1mach(1) ) then
      alphaGenericReal = 0; return
    end if

    beta  = adim%betaQCD('beta')
    bCoef = adim%betaQCD('bCoef')

    if ( self%method(:8) == 'analytic' ) then

      L = log(mu/mZ);  arg = amZ * L * beta(0)/2/Pi;  LG  = log(1 + arg)
      alphaGenericReal = 1

      if (self%run > 0) alphaGenericReal = alphaGenericReal + arg
      if (self%run > 1) alphaGenericReal = alphaGenericReal + beta(1)/beta(0) * amZ * LG/4/Pi

      if (self%run > 2) alphaGenericReal = alphaGenericReal + 1/( 16 * Pi2 + &
       8 * Pi * amZ * L * beta(0) ) * (  amZ**2 * LG * beta(1)**2/beta(0)**2 - &
       ( amZ**3 * L * beta(1)**2/2/Pi )/beta(0) + amZ**3 * L * beta(2)/2/Pi  )

      if (self%run > 3) alphaGenericReal = alphaGenericReal + ( -amZ**3 * LG**2 * beta(1)**3/beta(0)**3 + &
       amZ**5 * L**2 * beta(1)**3/beta(0)/4/Pi2 - amZ**5 * L**2 * beta(1) * beta(2)/2/Pi2 + &
       2 * amZ**3 * LG * beta(1) * beta(2)/beta(0)**2 - amZ**4 * L * beta(1) * beta(2)/beta(0)/Pi + &
       amZ**4 * L * beta(3)/Pi + amZ**5 * L**2 * beta(0) * beta(3)/4/Pi2 )/32/Pi/( 2 * Pi + amZ * L * beta(0) )**2

      if (self%run > 4) alphaGenericReal = alphaGenericReal + 4 * amZ**4 * ( - arg * (   6 * beta(0)**2 * &
      ( beta(1) * beta(3) - beta(0) * beta(4) ) + 2 * arg**2 * (  beta(1)**4 - 3 * beta(0) * beta(1)**2 * beta(2) &
      + 2 * beta(0)**2 * beta(1) * beta(3) + beta(0)**2 * ( beta(2)**2 - beta(0) * beta(4) )  ) + &
      3 * arg * (  beta(1)**4 - 4 * beta(0) * beta(1)**2 * beta(2) + 3 * beta(0)**2 * beta(1) * beta(3) + &
      2 * beta(0)**2 * ( beta(2)**2 - beta(0) * beta(4) )  )   ) + beta(1) * LG * (6 * arg * ( beta(1)**3 &
      - beta(0) * beta(1) * beta(2)) + 6 * beta(0)**2 * beta(3) - 6 * beta(0) * beta(1) * beta(2) * LG +&
      beta(1)**3 * LG * (2 * LG - 3)))/(6144 * (1 + arg)**3 * beta(0)**4 * Pi2**2)

      alphaGenericReal = amZ/alphaGenericReal

    else if ( self%method(:7) == 'numeric' ) then

      L = Log(Mz/mu);  h = 0.04_dp; n = max(  1, Abs(  Nint(L/h)  )  )

      h = - L/n;  alphaGenericReal = amZ; ord = min(self%run,5)

      do i = 1, n

       k1 = h * PiBeta( beta(:ord-1), alphaGenericReal        )
       k2 = h * PiBeta( beta(:ord-1), alphaGenericReal + k1/2 )
       k3 = h * PiBeta( beta(:ord-1), alphaGenericReal + k2/2 )
       k4 = h * PiBeta( beta(:ord-1), alphaGenericReal + k3   )

       alphaGenericReal = alphaGenericReal + (k1 + k4)/6 + (k2 + k3)/3

      end do

    else if ( self%method(:9) == 'iterative' ) then

      if (self%run <= 0) then
        alphaGenericReal = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      cCoef = adim%betaQCD('cCoef')
      cCoef(1:self%run) = cCoef(1:self%run)/PowList(4 * Pi,self%run)

      alphaGenericReal = 1/aLLInv

      if ( self%run <= 1) return

      do i = 1, 100
        corr = 1/iter(cCoef(1:self%run), alphaGenericReal, amZ, aLLInv)
        if ( abs(corr - alphaGenericReal) < 1e-10_dp ) return
        alphaGenericReal = corr
      end do

    else if ( self%method(:4) == 'root' ) then

      if (self%run <= 0) then
        alphaGenericReal = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      alphaGenericReal = 1/aLLInv

      if ( self%run <= 1) return

      do i = 1, 100
        corr = 1/adim%root(self%run - 1, alphaGenericReal, amZ, aLLInv)
        if ( abs(corr - alphaGenericReal) < 1e-10_dp ) return
        alphaGenericReal = corr
      end do

    else if ( self%method(:6) == 'expand' ) then

      if (self%run <= 0) then
        alphaGenericReal = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      alphaGenericReal = 1/aLLInv

      if ( self%run <= 1) return

      ExpandCoeff = adim%cCoeff(self%run - 1, 20)

      do i = 1, 100
        corr = 1/expand(amZ, alphaGenericReal)
        if ( abs(corr - alphaGenericReal) < 1e-10_dp ) return
        alphaGenericReal = corr
      end do

    end if

  contains

    pure real (dp) function expand(a0, aMu)
      real (dp), intent(in) :: a0, aMu
      real (dp)             :: corr
      integer               :: i

      expand = aLLInv + ExpandCoeff(1) * log(aMu/a0)

      do i = 2, 20
        corr = ExpandCoeff(i) * ( aMu**(i-1) - a0**(i-1) )/(i - 1)
        if ( abs(corr) <= 1e-10_dp ) return
        expand = expand + corr
      end do

    end function

  end function alphaGenericReal

!ccccccccccccccc

  pure real (dp) function root(c, aMu, a0, aLLInv)
    real (dp)              , intent(in) :: aMu, a0, aLLInv
    real (dp), dimension(:), intent(in) :: c
    real (dp)                           :: radical, raiz, aSuPi, aSuPi0

    aSuPi = aMu/fourPi; aSuPi0 = a0/fourPi

    root = aLLInv - c(1)/fourPi * log(aMu/a0)

    if ( size(c) == 1 ) then

      root = root + c(1)/fourPi * log( (1 + c(1) * aSuPi)/(1 + c(1) * aSuPi0)  )

    else if ( size(c) == 2 ) then

      radical = c(1)**2 - 4 * c(2); raiz = sqrt( abs(radical) )

      root = root + c(1)/fourPi/2 * log( (1 + c(1) * aSuPi + c(2) * aSuPi**2)/&
      (1 + c(1) * aSuPi0 + c(2) * aSuPi0**2)  )

      if (radical < 0) then

        root = root + ( c(1)**2 - 2 * c(2) )/fourPi/raiz * &
        (  atan( (c(1) + 2 * c(2) * aSuPi)/raiz ) - atan( (c(1) + 2 * c(2) * aSuPi0)/raiz )  )

      else

        root = root - ( c(1)**2 - 2 * c(2) )/fourPi/raiz * &
        (  atanh( (c(1) + 2 * c(2) * aSuPi)/raiz ) - atanh( (c(1) + 2 * c(2) * aSuPi0)/raiz )  )

      end if

    end if

  end function root

!ccccccccccccccc

  pure real (dp) function iterReal(cCoef, aMu, a0, aLLInv)
    real (dp)              , intent(in) :: aMu, a0, aLLInv
    real (dp), dimension(:), intent(in) :: cCoef
    integer                             :: i

    iterReal = aLLInv

    if ( size(cCoef) > 0 ) iterReal = iterReal + cCoef(1) * Log(aMu/a0)

    if ( size(cCoef) < 2 ) return

    do i = 1, size(cCoef) - 1
      iterReal = iterReal + cCoef(i+1) * (aMu**i - a0**i)/i
    end do

  end function iterReal

!ccccccccccccccc

  pure complex (dp) function iterComplex(cCoef, aMu, a0, aLLInv)
    complex (dp)           , intent(in) :: aMu, aLLInv
    real (dp)              , intent(in) :: a0
    real (dp), dimension(:), intent(in) :: cCoef
    integer                             :: i

    iterComplex = aLLInv

    if ( size(cCoef) > 0 ) iterComplex = iterComplex + cCoef(1) * Log(aMu/a0)

    if ( size(cCoef) < 2 ) return

    do i = 1, size(cCoef) - 1
      iterComplex = iterComplex + cCoef(i+1) * (aMu**i - a0**i)/i
    end do

  end function iterComplex

!ccccccccccccccc

  pure complex (dp) function alphaGenericComplex(self, adim, mZ, amZ, mu)
    class (Alpha)            , intent(in) :: self
    type (AnomDim)           , intent(in) :: adim
    real    (dp)             , intent(in) :: mZ, amZ
    complex (dp)             , intent(in) :: mu
    real (dp), dimension(0:4)             :: beta
    real (dp)                             :: h, theta, mod
    integer                               :: n, i, ord
    real (dp), dimension(0:4)             :: bCoef, cCoef
    complex (dp)                          :: L, arg, LG, k1, k2, k3, k4, aLLinv, corr

    if ( max( amZ, mZ ) <= d1mach(1) ) then
      alphaGenericComplex = 0; return
    end if

    beta  = adim%betaQCD('beta')
    bCoef = adim%betaQCD('bCoef')

    if ( self%method(:8) == 'analytic' ) then

    L = log(mu/mZ);  arg = amZ * L * beta(0)/2/Pi;  LG  = log(1 + arg)
    alphaGenericComplex = 1

    if (self%run > 0) alphaGenericComplex = alphaGenericComplex + arg
    if (self%run > 1) alphaGenericComplex = alphaGenericComplex + beta(1)/beta(0) * amZ * LG/4/Pi

    if (self%run > 2) alphaGenericComplex = alphaGenericComplex + 1/( 16 * Pi2 + &
     8 * Pi * amZ * L * beta(0) ) * (  amZ**2 * LG * beta(1)**2/beta(0)**2 - &
     ( amZ**3 * L * beta(1)**2/2/Pi )/beta(0) + amZ**3 * L * beta(2)/2/Pi  )

    if (self%run > 3) alphaGenericComplex = alphaGenericComplex + ( -amZ**3 * LG**2 * beta(1)**3/beta(0)**3 + &
     amZ**5 * L**2 * beta(1)**3/beta(0)/4/Pi2 - amZ**5 * L**2 * beta(1) * beta(2)/2/Pi2 + &
     2 * amZ**3 * LG * beta(1) * beta(2)/beta(0)**2 - amZ**4 * L * beta(1) * beta(2)/beta(0)/Pi + &
     amZ**4 * L * beta(3)/Pi + amZ**5 * L**2 * beta(0) * beta(3)/4/Pi2 )/32/Pi/( 2 * Pi + amZ * L * beta(0) )**2

    if (self%run > 4) alphaGenericComplex = alphaGenericComplex + 4 * amZ**4 * ( - arg * (   6 * beta(0)**2 * &
    ( beta(1) * beta(3) - beta(0) * beta(4) ) + 2 * arg**2 * (  beta(1)**4 - 3 * beta(0) * beta(1)**2 * beta(2) &
    + 2 * beta(0)**2 * beta(1) * beta(3) + beta(0)**2 * ( beta(2)**2 - beta(0) * beta(4) )  ) + &
    3 * arg * (  beta(1)**4 - 4 * beta(0) * beta(1)**2 * beta(2) + 3 * beta(0)**2 * beta(1) * beta(3) + &
    2 * beta(0)**2 * ( beta(2)**2 - beta(0) * beta(4) )  )   ) + beta(1) * LG * (6 * arg * ( beta(1)**3 &
    - beta(0) * beta(1) * beta(2)) + 6 * beta(0)**2 * beta(3) - 6 * beta(0) * beta(1) * beta(2) * LG +&
    beta(1)**3 * LG * (2 * LG - 3)))/(6144 * (1 + arg)**3 * beta(0)**4 * Pi2**2)

    alphaGenericComplex = amZ/alphaGenericComplex

    else if ( self%method(:7) == 'numeric' ) then

      mod = abs(mu); theta = IMAGPART( log(mu) - log(mod) ); ord = min(self%run,5)

      alphaGenericComplex = self%alphaGenericReal(adim, mZ, amZ, mod)

      h = 0.04_dp; n = max(  1, Abs( Nint(theta/h) )  )

      h = theta/n

      do i = 1, n

       k1 = h * PiBeta( beta(:ord-1), alphaGenericComplex        )
       k2 = h * PiBeta( beta(:ord-1), alphaGenericComplex + k1/2 )
       k3 = h * PiBeta( beta(:ord-1), alphaGenericComplex + k2/2 )
       k4 = h * PiBeta( beta(:ord-1), alphaGenericComplex + k3   )

       alphaGenericComplex = alphaGenericComplex + (k1 + k4)/6 + (k2 + k3)/3

      end do

  else if ( self%method(:9) == 'iterative' ) then

    if (self%run <= 0) then
      alphaGenericComplex = aMz; return
    end if

    aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi
    alphaGenericComplex = 1/aLLInv

    cCoef = adim%betaQCD('cCoef')
    cCoef(1:self%run) = cCoef(1:self%run)/PowList(4 * Pi,self%run)

    if ( self%run <= 1) return

    do i = 1, 100

      corr = 1/iter( cCoef(1:self%run), alphaGenericComplex, amZ, aLLInv)
      if ( abs(corr - alphaGenericComplex) < 1e-10_dp ) return
      alphaGenericComplex = corr

    end do

    end if

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
