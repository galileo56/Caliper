
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
    character (len = 9)             :: method, str
    integer                         :: order, run, n
    type (AnomDim) , dimension(3:6) :: andim
    real (dp)      , dimension(3:6) :: alphaRef, muRef
    real (dp)                       :: mT, mB, mC
    logical                         :: QmT, QmB, QmC, QmuT, QmuC, QmuB

   contains

   generic                          :: alphaQCD     => alphaQCDReal, alphaQCDComplex
   generic, public                  :: alphaGeneric => alphaGenericReal, alphaGenericComplex
   generic, public                  :: alphaGenericFlavor => alphaGenericRealFlavor, alphaGenericComplexFlavor

   procedure, pass(self)            :: alphaQCDReal, alphaQCDComplex, alphaGenericReal, &
   alphaGenericComplex, thresholdMatching, alphaGenericRealFlavor, alphaGenericComplexFlavor

   procedure, pass(self), public    :: scales, orders, scheme, SetAlpha, SetMTop, &
   SetMBottom, SetMCharm, adim

  end type Alpha

!ccccccccccccccc

  interface Alpha
    module procedure  InitAlpha
  end interface Alpha

  contains

!ccccccccccccccc

  type (Alpha) function InitAlpha(andimList, order, run, mZ, amZ, mT, muT,&
    mB, muB, mC, muC, method)
    character (len = *), optional, intent(in) :: method
    real (dp)         , optional , intent(in) :: amZ
    real (dp)         , optional , intent(in) :: mZ, mT, muT, mB, muB, mC, muC
    type (AnomDim) , dimension(4), intent(in) :: andimlist
    integer                      , intent(in) :: order
    integer                      , intent(in) :: run

    InitAlpha%andim = andimlist

    InitAlpha%muRef = 0    ; InitAlpha%run = min(5,run)
    InitAlpha%muRef(5) = mZ; InitAlpha%n   = min(4,order - 1)

    if ( present(method) ) then
      InitAlpha%method = method
    else
      InitAlpha%method = 'analytic'
    end if

! initialising all Anomalous Dimension Objects

    if ( present(mB) ) then; InitAlpha%mB = mB; InitAlpha%QmB = .true.; end if
    if ( present(mC) ) then; InitAlpha%mC = mC; InitAlpha%QmC = .true.; end if
    if ( present(mT) ) then; InitAlpha%mT = mT; InitAlpha%QmT = .true.; end if

    if ( present(muB) ) then
      InitAlpha%muRef(4) = muB; InitAlpha%QmuB = .true.
    end if

    if ( present(muT) ) then
      InitAlpha%muRef(6) = muT; InitAlpha%QmuT = .true.
    end if

    if ( present(muC) ) then
      InitAlpha%muRef(3) = muC; InitAlpha%QmuC = .true.
    end if

    InitAlpha%str   = andimlist(3)%scheme()
    InitAlpha%order = order ;  InitAlpha%alphaRef = 0

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
    class (Alpha)  , intent(inout) :: self
    real (dp)      , intent(in   ) :: mT, muT
    real (dp)                      :: aS
    real (dp), dimension(0:self%n) :: alphaList

    self%mT = mT; self%muRef(6) = muT; alphaList(0) = 1

    aS  = self%alphaGeneric(self%method, self%andim(5), self%muRef(5), self%alphaRef(5), muT)

    self%alphaRef(6) = pi * sum(   PowList(aS/Pi, self%n+1) * &
    self%thresholdMatching( 'up', 6, log(muT/mT) )  )

  end subroutine SetMTop

!ccccccccccccccc

  subroutine SetMBottom(self, mB, muB)
    class (Alpha)  , intent(inout) :: self
    real (dp)      , intent(in   ) :: mB, muB
    real (dp)                      :: aS
    real (dp), dimension(0:self%n) :: alphaList

    self%mB = mB; self%muRef(4) = muB; alphaList(0) = 1

    aS  = self%alphaGeneric(self%method, self%andim(5), self%muRef(5), self%alphaRef(5), muB )

    self%alphaRef(4) = pi * sum(  self%thresholdMatching( 'down', 5, log(muB/mB) ) &
    * PowList(aS/Pi, self%n+1)  )

!   Running from mB to muC, with nf = 4 flavors, and matching

    if (self%QmuC .and. self%QmC) call self%SetMCharm( self%mC, self%muRef(3) )

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, mC, muC)
    class (Alpha)  , intent(inout) :: self
    real (dp)      , intent(in   ) :: mC, muC
    real (dp)                      :: aS
    real (dp), dimension(0:self%n) :: alphaList

    self%mC = mC; self%muRef(3) = muC; alphaList(0) = 1

    aS  = self%alphaGeneric(self%method, self%andim(4), self%muRef(4), self%alphaRef(4), muC)

    self%alphaRef(3) = pi * sum(  self%thresholdMatching( 'down', 4, log(muC/mC) ) &
    * PowList(aS/Pi, self%n+1)  )

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
      alphaQCDReal = self%alphaGeneric(self%method, self%andim(n), self%muRef(n), self%alphaRef(n), mu )
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
      alphaQCDComplex = self%alphaGeneric(self%method, self%andim(n), self%muRef(n), self%alphaRef(n), mu )
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

pure real (dp) function alphaGenericRealFlavor(self, method, nf, mZ, amZ, mu)
  class (Alpha)      , intent(in) :: self
  character (len = *), intent(in) :: method
  real (dp)          , intent(in) :: mZ, amZ, mu
  integer            , intent(in) :: nf

  alphaGenericRealFlavor = self%alphaGeneric( method, self%adim(nf), mZ, amZ, mu )

end function alphaGenericRealFlavor

!ccccccccccccccc

complex (dp) function alphaGenericComplexFlavor(self, method, nf, mZ, amZ, mu)
  class (Alpha)      , intent(in) :: self
  character (len = *), intent(in) :: method
  real (dp)          , intent(in) :: mZ, amZ
  complex (dp)       , intent(in) :: mu
  integer            , intent(in) :: nf

  alphaGenericComplexFlavor = self%alphaGeneric( method, self%adim(nf), mZ, amZ, mu )

end function alphaGenericComplexFlavor

!ccccccccccccccc

 pure real (dp) function alphaGenericReal(self, method, adim, mZ, amZ, mu)
   class (Alpha)            , intent(in) :: self
   character (len = *)      , intent(in) :: method
   type (AnomDim)           , intent(in) :: adim
   real (dp)                , intent(in) :: mZ, amZ, mu
   integer                               :: n, i, ord
   real (dp), dimension(0:4)             :: beta
   real (dp), dimension(0:20)            :: a0List, aLList
   real (dp), dimension(-20:2,0:20)      :: b ! (2,:) corresponds to log expansion
   real (dp), dimension(0:20)            :: cCoef
   real (dp)                             :: L, h, k1, k2, k3, k4, aLL, aLLInv, corr, a0

    alphaGenericReal = 0;  if ( max( amZ, mZ, mu ) <= d1mach(1) ) return

    if (self%run == 0) then; alphaGenericReal = amZ; return; end if

    beta  = adim%betaQCD('beta') ; cCoef(:4) = adim%betaQCD('cCoef')

    if ( method(:6) == 'series' ) then
      a0 = amZ/fourPi;  b = 0; b(:1,0) = 1
      cCoef = adim%cCoeff(self%run - 1, 20)

      if (self%run > 0) then
         aLL = fourPi/( 1/a0 + 2 * log(mu/mZ) * beta(0) ); b(1,0) = 1
         alphaGenericReal = 1; if ( self%run == 1 ) go to 20
       end if

      if (self%run > 1) then
        b(1,1) = aLL * cCoef(1) * log(aLL/amZ)
        alphaGenericReal = 1 + b(1,1)
        a0List(0) = 1; a0List(1:) = powList(amZ , 20)
        aLList(0) = 1; aLList(1:) = powList(aLL, 20)
      end if

      do n = 2, 20

        b(2, n - 1) = b(1, n - 1) - sum(  b(1,n-2:1:-1) * b(2,1:n-2) * &
        [ (i, i = 1, n - 2) ]  )/(n - 1)

        do i = 1, n
          b(-i, n - 1) = b(1 - i, n - 1) - Sum( b(1, n - 1:1:-1) * b(-i, :n - 2) )
        end do

        do i = 1, n - 1
          b(1 - n, i) = b(2 - n, i) - Sum( b(1, i:1:-1) * b(1 - n, :i - 1) )
        end do

        do i = 2, n
          b(1, n) = b(1, n) + cCoef(i) * aLList(i) * b(1 - i, n - i)/(i - 1)
        end do

        b(1, n) = b(1, n) - cCoef(n) * aLL * a0List(n - 1)/(n - 1) &
        - cCoef(1) * aLL * b(2, n - 1)

        alphaGenericReal = alphaGenericReal + b(1,n)

        if (  abs( b(1,n) ) <= 1e-10_dp  ) exit

      end do

  20    alphaGenericReal = aLL/alphaGenericReal

    end if

    if ( method(:8) == 'analytic' .or. method(:7) == 'inverse' ) then

      a0 = amZ/fourPi;  b = 0; b(:1,0) = 1

      if (self%run > 0) then
         aLL = 1/( 1/a0 + 2 * log(mu/mZ) * beta(0) ); b(1,0) = 1
       end if

      if (self%run > 1) b(1,1) = aLL * cCoef(1) * log(aLL/a0)

      if (self%run > 2) then
        a0List(0) = 1; a0List(1:self%run - 1) = powList(a0 , self%run - 1)
        aLList(0) = 1; aLList(1:self%run - 1) = powList(aLL, self%run - 1)
      end if

      do n = 2, self%run - 1

        b(2, n - 1) = b(1, n - 1) - sum(  b(1,n-2:1:-1) * b(2,1:n-2) * &
        [ (i, i = 1, n - 2) ]  )/(n - 1)

        do i = 1, self%run - 2
          b(-i, n - 1) = b(1 - i, n - 1) - Sum( b(1, n - 1:1:-1) * b(-i, :n-2) )
        end do

        do i = 2, n
          b(1, n) = b(1, n) + cCoef(i) * aLList(i) * b(1 - i, n - i)/(i - 1)
        end do

        b(1, n) = b(1, n) - cCoef(n) * aLL * a0List(n - 1)/(n - 1) &
        - cCoef(1) * aLL * b(2, n - 1)

      end do

      if ( method(:8) == 'analytic') then
        alphaGenericReal = fourpi * aLL/Sum( b(1,:self%run-1) )
      else

        if (self%run > 1) b(-1, self%run - 1) = - sum( b(-1,:self%run - 2) * &
        b(1,self%run - 1:1:-1) )

        alphaGenericReal = fourpi * aLL * Sum( b(-1,:self%run-1) )

      end if

    else if ( method(:7) == 'numeric' ) then

      L = Log(Mz/mu);  h = 0.04_dp; n = max(  1, Abs(  Nint(L/h)  )  )

      h = - L/n;  alphaGenericReal = amZ; ord = min(self%run,5)

      do i = 1, n

       k1 = h * PiBeta( beta(:ord-1), alphaGenericReal        )
       k2 = h * PiBeta( beta(:ord-1), alphaGenericReal + k1/2 )
       k3 = h * PiBeta( beta(:ord-1), alphaGenericReal + k2/2 )
       k4 = h * PiBeta( beta(:ord-1), alphaGenericReal + k3   )

       alphaGenericReal = alphaGenericReal + (k1 + k4)/6 + (k2 + k3)/3

      end do

    else if ( method(:9) == 'iterative' ) then

      if (self%run <= 0) then
        alphaGenericReal = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      cCoef(1:self%run) = cCoef(1:self%run)/PowList(fourPi,self%run)

      alphaGenericReal = 1/aLLInv

      if ( self%run <= 1) return

      do i = 1, 100
        corr = 1/iter(cCoef(1:self%run), alphaGenericReal, amZ, aLLInv)
        if ( abs(corr - alphaGenericReal) < 1e-10_dp ) return
        alphaGenericReal = corr
      end do

    else if ( method(:4) == 'root' ) then

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

    else if ( method(:6) == 'expand' ) then

      if (self%run <= 0) then
        alphaGenericReal = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      alphaGenericReal = 1/aLLInv

      if ( self%run <= 1) return

      cCoef = adim%cCoeff(self%run - 1, 20)

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

      expand = aLLInv + cCoef(1) * log(aMu/a0)

      do i = 2, 20
        corr = cCoef(i) * ( aMu**(i-1) - a0**(i-1) )/(i - 1)
        if ( abs(corr) <= 1e-10_dp ) return
        expand = expand + corr
      end do

    end function

  end function alphaGenericReal

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

  pure complex (dp) function alphaGenericComplex(self, method, adim, mZ, amZ, mu)
    class (Alpha)            , intent(in) :: self
    character (len = *)      , intent(in) :: method
    type (AnomDim)           , intent(in) :: adim
    real    (dp)             , intent(in) :: mZ, amZ
    complex (dp)             , intent(in) :: mu
    real (dp)                             :: h, theta, mod
    integer                               :: n, i, ord
    complex (dp), dimension(0:20)         :: a0List, aLList
    real (dp)   , dimension(0:4)          :: beta
    real (dp)   , dimension(0:20)         :: cCoef
    complex (dp), dimension(-20:2,0:20)   :: b ! (:,2) corresponds to log expansion
    complex (dp)                          :: a0, aLL, k1, k2, k3, k4, aLLinv, corr

    alphaGenericComplex = 0; if ( max( amZ, mZ ) <= d1mach(1) ) return

    beta = adim%betaQCD('beta'); cCoef(:4) = adim%betaQCD('cCoef')

    if ( method(:6) == 'series' ) then
      a0 = amZ/fourPi;  b = 0; b(:1,0) = 1
      cCoef = adim%cCoeff(self%run - 1, 20)

      if (self%run > 0) then
         aLL = fourPi/( 1/a0 + 2 * log(mu/mZ) * beta(0) ); b(1,0) = 1
         alphaGenericComplex = 1; if ( self%run == 1 ) go to 30
       end if

      if (self%run > 1) then
        b(1,1) = aLL * cCoef(1) * log(aLL/amZ)
        alphaGenericComplex = 1 + b(1,1)
        a0List(0) = 1; a0List(1:) = powList(amZ , 20)
        aLList(0) = 1; aLList(1:) = powList(aLL, 20)
      end if

      do n = 2, 20

        b(2, n - 1) = b(1, n - 1) - sum(  b(1,n-2:1:-1) * b(2,1:n-2) * &
        [ (i, i = 1, n - 2) ]  )/(n - 1)

        do i = 1, n
          b(-i, n - 1) = b(1 - i, n - 1) - Sum( b(1, n - 1:1:-1) * b(-i, :n - 2) )
        end do

        do i = 1, n - 1
          b(1 - n, i) = b(2 - n, i) - Sum( b(1, i:1:-1) * b(1 - n, :i - 1) )
        end do

        do i = 2, n
          b(1, n) = b(1, n) + cCoef(i) * aLList(i) * b(1 - i, n - i)/(i - 1)
        end do

        b(1, n) = b(1, n) - cCoef(n) * aLL * a0List(n - 1)/(n - 1) &
        - cCoef(1) * aLL * b(2, n - 1)

        alphaGenericComplex = alphaGenericComplex + b(1,n)

        if (  abs( b(1,n) ) <= 1e-10_dp  ) exit

      end do

  30    alphaGenericComplex = aLL/alphaGenericComplex

    end if

    if ( method(:8) == 'analytic' .or. method(:7) == 'inverse' ) then

      a0 = amZ/fourPi;  b = 0; b(:1,0) = 1

      if (self%run > 0) then
         aLL = 1/( 1/a0 + 2 * log(mu/mZ) * beta(0) ); b(1,0) = 1
       end if

      if (self%run > 1) b(1,1) = aLL * cCoef(1) * log(aLL/a0)

      if (self%run > 2) then
        a0List(0) = 1; a0List(1:self%run - 1) = powList(a0 , self%run - 1)
        aLList(0) = 1; aLList(1:self%run - 1) = powList(aLL, self%run - 1)
      end if

      do n = 2, self%run - 1

        b(2, n - 1) = b(1, n - 1) - sum(  b(1,n-2:1:-1) * b(2,1:n-2) * &
        [ (i, i = 1, n - 2) ]  )/(n - 1)

        do i = 1, self%run - 2
          b(-i, n - 1) = b(1 - i, n - 1) - Sum( b(1, n - 1:1:-1) * b(-i, :n-2) )
        end do

        do i = 2, n
          b(1, n) = b(1, n) + cCoef(i) * aLList(i) * b(1 - i, n - i)/(i - 1)
        end do

        b(1, n) = b(1, n) - cCoef(n) * aLL * a0List(n - 1)/(n - 1) &
        - cCoef(1) * aLL * b(2, n - 1)

      end do

      if ( method(:8) == 'analytic') then
        alphaGenericComplex = fourpi * aLL/Sum( b(1,:self%run-1) )
      else

        if (self%run > 1) b(-1, self%run - 1) = - sum( b(-1,:self%run - 2) * &
        b(1,self%run - 1:1:-1) )

        alphaGenericComplex = fourpi * aLL * Sum( b(-1,:self%run-1) )

      end if

    else if ( method(:7) == 'numeric' ) then

      mod = abs(mu); theta = IMAGPART( log(mu) - log(mod) ); ord = min(self%run,5)

      alphaGenericComplex = self%alphaGenericReal(method, adim, mZ, amZ, mod)

      h = 0.04_dp; n = max(  1, Abs( Nint(theta/h) )  )

      h = theta/n

      do i = 1, n

       k1 = h * PiBeta( beta(:ord-1), alphaGenericComplex        )
       k2 = h * PiBeta( beta(:ord-1), alphaGenericComplex + k1/2 )
       k3 = h * PiBeta( beta(:ord-1), alphaGenericComplex + k2/2 )
       k4 = h * PiBeta( beta(:ord-1), alphaGenericComplex + k3   )

       alphaGenericComplex = alphaGenericComplex + (k1 + k4)/6 + (k2 + k3)/3

      end do

    else if ( method(:9) == 'iterative' ) then

      if (self%run <= 0) then
        alphaGenericComplex = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi
      alphaGenericComplex = 1/aLLInv

      cCoef(1:self%run) = cCoef(1:self%run)/PowList(fourPi,self%run)

      if ( self%run <= 1) return

      do i = 1, 100

        corr = 1/iter( cCoef(1:self%run), alphaGenericComplex, amZ, aLLInv)
        if ( abs(corr - alphaGenericComplex) < 1e-10_dp ) return
        alphaGenericComplex = corr

      end do

    else if ( method(:4) == 'root' ) then

      if (self%run <= 0) then
        alphaGenericComplex = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      alphaGenericComplex = 1/aLLInv

      if ( self%run <= 1) return

      do i = 1, 100
        corr = 1/adim%root(self%run - 1, alphaGenericComplex, amZ, aLLInv)
        if ( abs(corr - alphaGenericComplex) < 1e-10_dp ) return
        alphaGenericComplex = corr
      end do

    else if ( method(:6) == 'expand' ) then

      if (self%run <= 0) then
        alphaGenericComplex = aMz; return
      end if

      aLLInv = 1/amZ + log(mu/mZ) * beta(0)/2/Pi

      alphaGenericComplex = 1/aLLInv

      if ( self%run <= 1) return

      cCoef = adim%cCoeff(self%run - 1, 20)

      do i = 1, 100
        corr = 1/expand(amZ, alphaGenericComplex)
        if ( abs(corr - alphaGenericComplex) < 1e-10_dp ) return
        alphaGenericComplex = corr
      end do

    end if

  contains

    pure complex (dp) function expand(a0, aMu)
      complex (dp), intent(in) :: aMu
      real (dp)   , intent(in) :: a0
      complex (dp)             :: corr
      integer                  :: i

      expand = aLLInv + cCoef(1) * log(aMu/a0)

      do i = 2, 20
        corr = cCoef(i) * ( aMu**(i-1) - a0**(i-1) )/(i - 1)
        if ( abs(corr) <= 1e-10_dp ) return
        expand = expand + corr
      end do

    end function

  end function alphaGenericComplex

!ccccccccccccccc

  function thresholdMatching(self, direction, nf, lg) result(e)
    class (Alpha)                  , intent(in) :: self
    character (len = *)            , intent(in) :: direction
    integer                        , intent(in) :: nf
    real (dp)                      , intent(in) :: lg
    real (dp), dimension(self%n + 1)            :: e
    real (dp), dimension(0:self%n  )            :: lgList
    real (dp), dimension(0:4, 5    )            :: b

    lgList(0) = 1; lgList(1:) = powList(lg, self%n)
    b = self%andim(nf)%MatchingAlphaLog(direction)

    e = matmul( lgList, b(:self%n, :self%n + 1) )

  end function thresholdMatching

!ccccccccccccccc

  ! function thresholdMatching2(self, nf, lg) result(e) ! smarter but slower
  !   class (Alpha)                  , intent(in) :: self
  !   integer                        , intent(in) :: nf
  !   real (dp)                      , intent(in) :: lg
  !   real (dp), dimension(self%n + 1)            :: e
  !   real (dp), dimension(0:self%n  )            :: lgList
  !   real (dp), dimension(0:4       , 0:5      ) :: b, c
  !   real (dp), dimension(self%n + 1,self%n + 1) :: ePow
  !   integer                                     :: n, i
  !
  !   lgList(0) = 1; lgList(1:) = powList(lg, self%n); e = 0; e(1) = 1; ePow = 0
  !
  !   b = 0; c = 0; c(0,1) = 1 ;  b(0,1:) = self%andim(nf)%MatchingAlpha()
  !   call self%andim(nf    )%expandAlpha( b(:,1:) )
  !   call self%andim(nf - 1)%expandAlpha( c(:,1:) )
  !
  !   b(0,:self%n + 1) = matmul( lgList, b(:self%n,:self%n + 1) )
  !   c(0,:self%n + 1) = matmul( lgList, c(:self%n,:self%n + 1) )
  !
  !   do n = 2, self%n + 1
  !
  !     ePow(1,n - 1) = e(n - 1)
  !
  !     do i = 2, n
  !       ePow(i,n) = sum( e(:n + 1 - i) * ePow(i - 1,n - 1:i - 1:-1) )
  !     end do
  !
  !     e(n) = b(0,n) - sum( c(0,2:n) * ePow(2:n,n) )
  !
  !   end do
  !
  ! end function thresholdMatching2

!ccccccccccccccc

  pure type (AnomDim) function adim(self, nf)
    class (Alpha), intent(in) :: self
    integer      , intent(in) :: nf
    integer                   :: n

    n = nf; if (nf > 6) n = 6; if (nf < 3) n = 3;  adim = self%andim(n)

  end function adim

!ccccccccccccccc

   pure character (len = 5) function scheme(self)
    class (Alpha), intent(in) :: self
    scheme = self%str
   end function scheme

!ccccccccccccccc

end module AlphaClass
