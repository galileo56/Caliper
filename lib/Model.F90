
module ModelClass
  use QuadPack, only: qagi; use Constants, only: dp, Pi, ExpEuler, d1mach, prec
  use AnomDimClass, only: powList; use MatrixElementsClass, only: factList
  use MCtopClass; use adapt; use QuadPack; implicit none
  private

  public    :: Model2D, BreitWigner, BreitModel2D

!ccccccccccccccc

  type, public                             :: Model
    private
    integer                                :: Nmax, n, m, clen
    real (dp)                              :: lambda
    real (dp),              dimension(5)   :: Omega, OmegaCum
    real (dp), allocatable, dimension(:,:) :: EmodCoef
    real (dp), allocatable, dimension(:  ) :: c, fourList, twoList, Zetalist, DModList, ListFact
    character (len = 5)                    :: str

  contains

    final                                 :: delete_object
    procedure, pass(self), public         :: BreitModel, SoftFourier, SetLambda, &
    CumMomentModel, ShapeFun, Taylor, ModelUnstable

    procedure, pass(self)                 :: CumSoft, MomentModel, EModCoefSum, &
    DModCoefSum, PlusSoft, ModCoef, DModCoef, BinomialInteger, BinomialReal,    &
    Binom, fact, CModCoef, Zeta

    generic, private  :: Binomial => BinomialInteger, BinomialReal

  end type Model

!ccccccccccccccc

  interface Model
    module procedure InModel
  end interface Model

  contains

!ccccccccccccccc

  subroutine delete_object(this)
    type(Model) :: this
    if ( allocated(this%EmodCoef) ) deallocate(this%EmodCoef)
    if ( allocated(this%c)        ) deallocate(this%c       )
    if ( allocated(this%FourList) ) deallocate(this%FourList)
    if ( allocated(this%TwoList ) ) deallocate(this%TwoList)
    if ( allocated(this%ZetaList) ) deallocate(this%ZetaList)
    if ( allocated(this%DmodList) ) deallocate(this%DmodList)
    if ( allocated(this%ListFact) ) deallocate(this%ListFact)
  end subroutine delete_object

!ccccccccccccccc

   type (Model) function InModel(lambda, c, piece, str)
     real (dp)          , intent(in), dimension(:) :: c
     integer            , intent(in), dimension(2) :: piece
     real (dp)          , intent(in)               :: lambda
     character (len = *), intent(in)               :: str
     integer                                       :: i, l

     InModel%n = piece(1);  InModel%m = piece(2); InModel%clen = size(c)

     allocate( InModel%c(InModel%clen) )

     InModel%Nmax = 0 ;  InModel%c = 0; InModel%lambda = lambda; InModel%str = str

     if ( str(:3) == "sum"   ) InModel%Nmax = 2 * (InModel%clen - 1)
     if ( str(:5) == "piece" ) InModel%Nmax = InModel%n + InModel%m
     if ( str(:3) == "sum"   ) InModel%c    = c

     allocate( InModel%EmodCoef(0:InModel%Nmax,0:3*InModel%Nmax) )
     allocate( InModel%FourList(0:3*InModel%Nmax), InModel%TwoList(0:InModel%Nmax) )
     allocate( InModel%ListFact(0:3*InModel%Nmax + 8) )

     InModel%ListFact = FactList(3*InModel%Nmax + 8)

     InModel%FourList(0) = 1; InModel%FourList(1:) = powList(4._dp/3, 3*InModel%Nmax)
     InModel%TwoList(0)  = 1; InModel%TwoList(1:)  = powList(    - 2,   InModel%Nmax)

     if ( str(:3) == "sum" ) then

       allocate( InModel%Zetalist(0:InModel%clen - 1) )
       allocate( InModel%DmodList(0:2 * (InModel%clen - 1)) )

       do i = 0, InModel%clen - 1
         InModel%Zetalist(i) = InModel%Zeta(i)
       end do

       do i = 0, 2 * (InModel%clen - 1)
         InModel%DmodList(i) = InModel%DModCoefSum(i)
       end do

     else

       allocate( InModel%DmodList(0:InModel%n/2 + InModel%m/2) )

       do i = 0, InModel%n/2 + InModel%m/2
         InModel%DmodList(i) = InModel%DModCoef(i)
       end do

     end if

     InModel%EmodCoef = 0

     do i = 0, InModel%Nmax
       do l = 0, 3 * i

         if ( str(:3) == "sum" ) InModel%EmodCoef(i,l) = 128 * InModel%EModCoefSum(i, l)/3
         if ( str(:5) == "piece" ) InModel%EmodCoef(i,l) = 2**(7 + InModel%Nmax) * &
         sqrt( (2._dp * InModel%n + 1) * (2 * InModel%m + 1) ) * InModel%ModCoef(i, l)/3

       end do
     end do

     do i = 1, 5
       InModel%Omega(i) = InModel%MomentModel(i)
     end do

     InModel%OmegaCum = InModel%CumMomentModel(InModel%Omega)

   end function InModel

!ccccccccccccccc

  subroutine SetLambda(self, lam)
    class (Model), intent(inout) :: self
    real (dp)    , intent(in)    :: lam
    self%lambda = lam
  end subroutine SetLambda

!ccccccccccccccc

  real (dp) function ModelUnstable(self, MC, k, p, p2)
    class (Model)     , intent(in) :: self
    type (MCtop)      , intent(in) :: MC
    real(dp)          , intent(in) :: p
    real(dp), optional, intent(in) :: p2
    integer           , intent(in) :: k
    real (dp)                      :: plim, err
    ! integer                        :: neval, ier

    ModelUnstable = 0; plim = p

    if ( present(p2) ) then
      if (p2 <= p) return; plim = p2
    end if

    if ( MC%shape(:6) == 'Cparam' ) then
      ModelUnstable = dGauss( InteUns, 0._dp, min( MC%maxES(), plim/MC%Qval() ), prec )
    else
      call DAdapt(InteUns, 0._dp, min( MC%maxES(), plim/MC%Qval() ), 1, prec, &
      prec, ModelUnstable, ERR)
    end if

    ! call qags( InteUns, 0._dp, min( MC%maxES(), plim/MC%Qval() ), prec, &
    ! prec, ModelUnstable, err, neval, ier )

  contains

    real (dp) function InteUns(x)
      real(dp), intent(in) :: x

      if ( .not. present(p2) ) then
        InteUns = self%ShapeFun( k, p - x * MC%Qval() ) * MC%Distribution(x)
      else
        if ( p - x * MC%Qval() <= 0 ) then
          InteUns = self%ShapeFun( k, p2 - x * MC%Qval() ) * MC%Distribution(x)
        else
          InteUns = self%ShapeFun( k, p - x * MC%Qval(), p2 - x * MC%Qval() ) * &
          MC%Distribution(x)
        end if
      end if

    end function InteUns

  end function ModelUnstable

!ccccccccccccccc

  real (dp) function Taylor(self, k)
    class (Model), intent(in) :: self
    integer      , intent(in) :: k
    integer                   :: i, l
    real (dp)                 :: SoftPieceMat2

    Taylor = 0

    do i = 0, self%Nmax

      SoftPieceMat2 = 0

      do l = 0, min(3 * i, k - 3)
        SoftPieceMat2  = SoftPieceMat2 + self%EmodCoef(i,l) &
        * ( - 4 * (i + 1) )**(k - l - 3)/self%fact(k - l - 3)
      enddo

      Taylor = Taylor + self%TwoList(i) * SoftPieceMat2

    enddo

    Taylor = Taylor/self%lambda**(k + 1)

  end function Taylor

!ccccccccccccccc

  real (dp) function ShapeFun(self, k, p, p2)
    class (Model)      , intent(in) :: self
    integer            , intent(in) :: k
    real (dp)          , intent(in) :: p
    real (dp), optional, intent(in) :: p2
    integer                         :: i, l, s
    real (dp)                       :: x, x2, SoftPieceMat2, SoftPieceMat3, a, &
    SoftPieceMat3B, SoftPieceMat2B, b

    ShapeFun = 0

    if (k < -3) return
    if (  ( .not. present(p2) ) .and. p < 0 ) return
    if (  present(p2) .and. (p2 <= p .or. p2 < 0) ) return

    x = p/self%lambda; if ( present(p2) ) x2 = p2/self%lambda

    if (k == - 1) then
      if ( .not. present(p2) ) ShapeFun = self%CumSoft(x)
      if (       present(p2) ) ShapeFun = self%CumSoft(x,x2); return
    else if (k < - 1) then
      if (p > 0)         ShapeFun = self%PlusSoft(k + 2, p)
      if ( present(p2) ) ShapeFun = self%PlusSoft(k + 2, p, p2); return
    end if

    do i = 0, self%Nmax

      SoftPieceMat2 = 0; if ( present(p2) ) SoftPieceMat2B = 0

      do l = 0, 3 * i

        SoftPieceMat3 = 0; if ( present(p2) ) SoftPieceMat3B = 0

        do s = 0, min(k, l + 3)
          a = ( - 4 * (i + 1) )**(k - s)/self%fact(s)/self%fact(k - s)/self%fact(l + 3 - s)
          if (x > 0 .or. l + 3 > s) then
            if (p >= 0)        SoftPieceMat3  = SoftPieceMat3  +  x**(l + 3 - s) * a
            if ( present(p2) ) SoftPieceMat3B = SoftPieceMat3B + x2**(l + 3 - s) * a
          else
            if (p >= 0)        SoftPieceMat3  = SoftPieceMat3  + a
            if ( present(p2) ) SoftPieceMat3B = SoftPieceMat3B + a
          end if
        end do

        b = self%Listfact(l + 3) * self%EmodCoef(i,l)
        if (p >= 0)        SoftPieceMat2  = SoftPieceMat2  + b * SoftPieceMat3
        if ( present(p2) ) SoftPieceMat2B = SoftPieceMat2B + b * SoftPieceMat3B

      enddo

      if (p >= 0)        ShapeFun = ShapeFun + self%TwoList(i) * Exp( - 4 * (i + 1) * x  ) * SoftPieceMat2
      if ( present(p2) ) ShapeFun = ShapeFun - self%TwoList(i) * Exp( - 4 * (i + 1) * x2 ) * SoftPieceMat2B

    enddo

    ShapeFun = self%fact(k) * ShapeFun/self%lambda**(k + 1)
    if ( present(p2) ) ShapeFun = - ShapeFun

  end function ShapeFun

!ccccccccccccccc

  real (dp) function CumSoft(self, x, x2)
    class (Model)      , intent(in) :: self
    real (dp)          , intent(in) :: x
    real (dp), optional, intent(in) :: x2
    integer                         :: i, j, l
    real (dp)                       :: CumSoft2, CumSoft2B, CumSoft3, CumSoft3B, &
    a, b, c

    CumSoft = 0
    if (  ( .not. present(x2) ) .and. x <= 0 ) return
    if (  present(x2) .and. (x2 <= x .or. x2 < 0) ) return

    do i = 0, self%Nmax

      CumSoft2 = 0;  if ( present(x2) ) CumSoft2B = 0

      do l = 0, 3 * i

        CumSoft3 = 0;  if ( present(x2) ) CumSoft3B = 0

        do j = 0, l + 3
          a = ( 4._dp * (i + 1) )**(- 1 - j)/self%Listfact(l + 3 - j)
          if (x > 0)         CumSoft3  = CumSoft3  +  x**(l + 3 - j) * a
          if ( present(x2) ) CumSoft3B = CumSoft3B + x2**(l + 3 - j) * a
        end do

        b = self%EmodCoef(i,l) * self%Listfact(l + 3); c = ( 4._dp * (i + 1) )**(- 4 - l)
        if (x > 0) CumSoft2 = CumSoft2 + b * (  c - exp( - 4 * (i + 1) * x ) * CumSoft3 )
        if ( present(x2) ) CumSoft2B = CumSoft2B + b * (  c - exp( - 4 * (i + 1) * x2 ) * CumSoft3B )

      enddo

      if (x > 0)         CumSoft = CumSoft + self%TwoList(i) * CumSoft2
      if ( present(x2) ) CumSoft = CumSoft - self%TwoList(i) * CumSoft2B

    enddo

    if ( present(x2) ) CumSoft = - CumSoft

  end function CumSoft

!ccccccccccccccc

  real (dp) function PlusSoft(self, cum, x, x2)
    class (Model)      , intent(in) :: self
    real (dp)          , intent(in) :: x
    real (dp), optional, intent(in) :: x2
    integer            , intent(in) :: cum
    real (dp)                       :: abserr, xmid, fac, fac2, corr
    integer                         :: neval, ier, i

    xmid = x; PlusSoft = 0

    if ( present(x2) ) then
      xmid = (x + x2)/2; if ( x2 < x .or. x2 < 0) return
    end if

    if (cum == -1) then
      fac  = 1
      if ( present(x2) ) then
        PlusSoft = self%MomentModel(0) * log(x2/x); fac2 = 1
      else
        PlusSoft = self%MomentModel(0) * log(x)
      end if
    else
      fac = 1/x
      if ( present(x2) ) then
        fac2 = 1/x2; PlusSoft = self%MomentModel(0) * (fac2 - fac)
      else
        PlusSoft = self%MomentModel(0) * fac
      end if
    end if

    do i = 1, 80

      fac = fac/x

      if ( present(x2) ) then
        fac2 = fac2/x2; corr = (fac2 - fac) * self%MomentModel(i)
      else
        corr = fac * self%MomentModel(i)
      end if

      if ( cum < 0 ) corr = - corr/i

      if ( abs(corr) > abs(PlusSoft) ) exit
      if ( abs(corr/PlusSoft) <= prec) return;  PlusSoft = PlusSoft + corr

    end do

    call qagi( IntePlusSoft, 0._dp, 1, prec, prec, PlusSoft, abserr, neval, ier )

    PlusSoft = - PlusSoft/pi/xmid

    contains

!ccccccccccccccc

    real (dp) function IntePlusSoft(z)
      real (dp), intent(in) :: z
      complex (dp)          :: y, C11, t, lgc, funct, CI, h

      CI = (0, 1);  C11 = (1, 1);  h = ( C11 * z - CI/10 );  y = h/xmid;  t = CI * h
      h  = CI * y;  lgc = log(h * ExpEuler)

      if ( .not. present(x2) ) then
        funct = C11 * lgc * exp(t) * self%SoftFourier(y)
      else
        funct = C11 * lgc * ( exp(h * x2) - exp(h * x) ) * self%SoftFourier(y)
      end if

      if (cum == -1) funct = funct/h;   IntePlusSoft = real(funct)

    end function IntePlusSoft

  end function PlusSoft

!ccccccccccccccc

  complex (dp) function SoftFourier(self, y)
    class (Model), intent(in) :: self
    complex (dp) , intent(in) :: y
    integer                   :: i, l
    complex (dp)              :: SoftPieceMat2, x

    x = y * self%lambda;  SoftFourier = 0

    do i = 0, self%Nmax

      SoftPieceMat2 = 0

      do l = 0, 3 * i
        SoftPieceMat2 = SoftPieceMat2 + self%EmodCoef(i,l) * self%Listfact(l + 3) / &
        ( 4 * (i + 1) + (0, 1) * x )**(l + 4)
      enddo

      SoftFourier = SoftFourier + self%TwoList(i) * SoftPieceMat2

    enddo

   end function SoftFourier

!ccccccccccccccc

  real (dp) function BreitModel(self, width, i, l, l2)
    class (Model)      , intent(in) :: self
    real (dp)          , intent(in) :: l
    real (dp), optional, intent(in) :: l2, width
    integer            , intent(in) :: i
    real (dp)                       :: abserr, corr, A, A2, Lg, Lg2
    real (dp), dimension(0:32)      :: derBreit, derBreit2, derPlus, derPlus2
    integer                         :: neval, ier, j, k, ii

    BreitModel = 0; if ( present(l2) .and. l2 < l ) return
    if ( i < -3 ) return; ii = i

    if (  ( .not. present(width) ) .or. width <= d1mach(1)  ) then
      if ( .not. present(l2) ) then
        BreitModel = self%ShapeFun(i, l)
      else
        BreitModel = self%ShapeFun(i, l, l2)
      end if
      return
    else if ( self%lambda <= d1mach(1) ) then
      BreitModel = BreitWigner(width, i, l)
      if ( present(l2) ) BreitModel = BreitWigner(width, i, l, l2)
      BreitModel = BreitModel * self%MomentModel(0)
      return
    end if

    if ( abs(l) <= d1mach(1)                    ) go to 20
    if ( present(l2) .and. abs(l2) <= d1mach(1) ) go to 20

    derBreit(0) = BreitWigner(width, i, l)

    if ( .not. present(l2) ) then
      BreitModel   = derBreit(0) * self%MomentModel(0)
    else
      derBreit2(0) = BreitWigner(width, i, l2)
      BreitModel   = ( derBreit2(0) - derBreit(0) ) * self%MomentModel(0)
      if (i == 0) derBreit2(0) = Pi * derBreit2(0)/width
    end if

    if (i == 0) derBreit(0) = Pi * derBreit(0)/width

    if (i == -1 .or. i == -3) then

      derBreit(0) = BreitWigner(width, i + 1, l)

      if ( .not. present(l2) ) then
        BreitModel  = BreitModel - derBreit(0) * self%MomentModel(1)
      else
        derBreit2(0) = BreitWigner(width, i + 1, l2)
        BreitModel   = BreitModel - ( derBreit2(0) - derBreit(0) ) * self%MomentModel(1)
        if (i == -1) derBreit2(0) = Pi * derBreit2(0)/width
      end if

      if (i == -1) derBreit(0) = Pi * derBreit(0)/width

    else if (i > 0) then
      derBreit(0) = Pi * BreitWigner(width, 0, l)/width
      if ( present(l2) ) derBreit2(0) = Pi * BreitWigner(width, 0, l2)/width
    end if

    derPlus = 0; derPlus2 = 0

    if (i < - 1) then
      derBreit(0) = Pi * BreitWigner(width, 0, l)/width
      A = ATan(l/width); Lg = log(l**2 + width**2)/2; ii = i + 2

      if ( present(l2) ) then
        derBreit2(0) = Pi * BreitWigner(width, 0, l2)/width
        A2 = ATan(l2/width); Lg2 = log(l2**2 + width**2)/2
      end if
    end if

    do j = 1, 32

      derBreit(j) = - 2 * l * derBreit(0) * derBreit(j - 1)
      if ( present(l2) ) derBreit2(j) = - 2 * l2 * derBreit2(0) * derBreit2(j - 1)

      do  k = 0, j - 2
        derBreit(j) = derBreit(j) - 2 * ( l * derBreit(j - 1 - k) + &
                      derBreit(j - 2 - k)  ) * derBreit(k)

        if ( present(l2) ) derBreit2(j) = derBreit2(j) - 2 * ( l2 * derBreit2(j - 1 - k) &
                     + derBreit2(j - 2 - k)  ) * derBreit2(k)
      end do

      derBreit(j) = derBreit(j)/j; if ( present(l2) ) derBreit2(j) = derBreit2(j)/j

      if (i < - 1) then

        do k = 2, j
          derPlus(j) = derPlus(j) + ( l * derBreit(k - 1) + derBreit(k - 2) ) * &
                                      derBreit(j - k)/k/(j + 1 - k)
          if ( present(l2) ) derPlus2(j) = derPlus2(j) + ( l2 * derBreit2(k - 1) +  &
                               derBreit2(k - 2) ) * derBreit2(j - k)/k/(j + 1 - k)
        end do

        derPlus(j) = (0.5_dp + A/Pi) * (l * derBreit(j) + derBreit(j - 1) ) + &
        width/Pi * ( (j + 1) * derPlus(j) + derBreit(j) * Lg + l * (1 + j)  * &
        derBreit(0) * derBreit(j - 1)/j )

        if ( present(l2) ) derPlus2(j) = (0.5_dp + A2/Pi) * (l2 * derBreit2(j) + derBreit2(j - 1) ) + &
        width/Pi * ( (j + 1) * derPlus2(j) + derBreit2(j) * Lg2 + l2 * (1 + j)  * &
        derBreit2(0) * derBreit2(j - 1)/j )

      else
        derPlus(j) = width/Pi * derBreit(j)
        if ( present(l2) ) derPlus2(j) = width/Pi * derBreit2(j)
      end if

      if (j <= i) cycle

      if ( .not. present(l2) ) then
        corr = self%MomentModel(j - ii) * derPlus(j) * (-1)**(j + i)
      else
        corr = self%MomentModel(j - ii) * ( derPlus2(j) - derPlus(j) ) * (-1)**(j + i)
      end if

      if (i == -1 .or. i == -3) then
        corr = corr/(j + 1)
      else if (i == 1) then
        corr = corr * j
      else if (i == 2) then
        corr = corr * j * (j - 1)
      end if

      if ( abs(corr) <= prec ) return
      if ( ( j > 5 + i .and. abs(corr) >= abs(BreitModel) ) .or. isNaN(corr) ) exit

      BreitModel = BreitModel + corr

    end do

    20 continue

    call qagi( inteBreitModel, 0._dp, 1, prec, prec, BreitModel, abserr, neval, ier )

    contains

!ccccccccccccccc

    real (dp) function inteBreitModel(p)
      real (dp), intent(in) :: p

      if (.not. present(l2) ) then
        inteBreitModel = BreitWigner(width, i, l - p) * self%ShapeFun(0, p)
      else
        inteBreitModel = BreitWigner(width, i, l - p, l2 - p) * self%ShapeFun(0, p)
      end if

    end function inteBreitModel

  end

!ccccccccccccccc

  function CumMomentModel(self, OmegaMat) result(cum)
    class (Model)          , intent(in)  :: self
    real (dp), dimension(:), intent(in)  :: OmegaMat

    real (dp), dimension( size(OmegaMat) ) :: cum
    integer                                :: k, j, n

    n = size(OmegaMat); cum = omegaMat

    do k = 2, n
      do j = 1, k - 1
        cum = cum - self%Binom(k - 1, j - 1)
      end do
    end do

  end function CumMomentModel

!ccccccccccccccc

  real (dp) function MomentModel(self, k)
    class (Model), intent(in) :: self
    integer      , intent(in) :: k
    integer                   :: i, l
    real (dp)                 :: MomentPiece2

    MomentModel = 0

    do i = 0, self%Nmax

      MomentPiece2 = 0

      do l = 0, 3 * i

        MomentPiece2 = MomentPiece2 + self%EmodCoef(i,l) * self%fact(k + l + 3) * &
                       ( 4._dp * (i + 1) )**(- 1 - k - l - 3)
      enddo

      MomentModel = MomentModel + self%TwoList(i) * MomentPiece2

    enddo

    MomentModel = MomentModel * self%lambda**k

   end function MomentModel

!ccccccccccccccc

   real (dp) function ModCoef(self, i, l)
     class (Model), intent(in) :: self
     integer      , intent(in) :: i, l
     integer                   :: h, kn, km

     ModCoef = 0; if ( i < 0 .or. i > self%Nmax .or. l < 0 .or. l > 3 * i ) return

     kn = mod(self%n, 2);  km = mod(self%m, 2)

     do h = max(0, ceiling(i/2.) - kn - km ), self%n/2 + self%m/2
       ModCoef = ModCoef + self%DmodList(h) * self%CModCoef(2 * h + kn + km, i, l)
     end do

     ModCoef = self%fourList(l) * ModCoef

   end function ModCoef

!ccccccccccccccc

   real (dp) function EModCoefSum(self, i, l)
     class (Model), intent(in) :: self
     integer      , intent(in) :: i, l
     integer                   :: h

     EModCoefSum = 0; if ( i < 0 .or. i > 2 * self%clen - 2 .or. l < 0 .or. l > 3 * i ) return

     do h = i, 2 * (self%clen - 1)
       EModCoefSum = EModCoefSum + self%DmodList(h) * self%CModCoef(h, i, l)
     end do

     EModCoefSum = self%fourList(l) * EModCoefSum

   end function EModCoefSum

!ccccccccccccccc

   real (dp) function DModCoef(self, h)
     class (Model), intent(in) :: self
     integer      , intent(in) :: h
     integer                   :: i, kn, km

     kn = self%n/2; km = self%m/2; DModCoef = 0; if ( (h > km + kn) .or. (h < 0) ) return

     do i = max( 0, h - km ), min( h, kn )

       DModCoef = DModCoef + self%Binomial( self%n, 2 * i + mod(self%n, 2) ) * &
       self%Binomial( self%m, 2 * (h - i) + mod(self%m, 2) ) * &
       self%Binomial( ceiling(self%n/2.) - 0.5_dp + i, self%n )      * &
       self%Binomial( ceiling(self%m/2.) - 0.5_dp + h - i, self%m)

     end do
   end function DModCoef

!ccccccccccccccc

   pure real (dp) function DModCoefSum(self, h)
     class (Model), intent(in) :: self
     integer      , intent(in) :: h
     integer                   :: j

     DModCoefSum = 0

     do j = max(0, h - self%clen), min(h, self%clen - 1)
       if ( h - j < 0 .or. h - j > self%clen - 1 ) cycle
       DModCoefSum = DModCoefSum + self%ZetaList(j) * self%ZetaList(h - j)
     end do

   end function DModCoefSum

!ccccccccccccccc

   real (dp) function Zeta(self, j)
     class (Model), intent(in) :: self
     integer      , intent(in) :: j
     integer                   :: k, mj

     Zeta = 0; if (j > self%clen - 1) return; mj = mod(j,2)

     do k = j/2, ( self%clen - 1 - mj )/2
       Zeta = Zeta + 4**k * sqrt(4 * k + 2 * mj + 1._dp) * self%c(2 * k + mj + 1) * &
       self%Binomial(2 * k + mj, j) * self%Binomial( ceiling(j/2.) - 0.5_dp + k, 2 * k + mj )
     end do

     Zeta = 2**mod(j,2) * Zeta

   end function Zeta

!ccccccccccccccc

   real (dp) function BinomialReal(self, i,j)
     class (Model), intent(in) :: self
     real (dp)    , intent(in) :: i
     integer      , intent(in) :: j
     integer                   :: k

    select case(j)
    case(0)  ;  BinomialReal = 1
    case(:-1);  BinomialReal = 0
    case default

      BinomialReal = i

      do k = 1, j - 1
        BinomialReal = BinomialReal * (i - k)
      end do

      BinomialReal = BinomialReal/self%Listfact(j)

    end select

   end function BinomialReal

!ccccccccccccccc

   real (dp) function BinomialInteger(self, i, j)
     class (Model), intent(in) :: self
     integer      , intent(in) :: i, j
     integer                   :: k

    select case(j)
    case(0)  ;  BinomialInteger = 1
    case(:-1);  BinomialInteger = 0
    case default

      BinomialInteger = i

      do k = 1, j - 1
        BinomialInteger = BinomialInteger * (i - k)
      end do

      BinomialInteger = BinomialInteger/self%Listfact(j)

    end select

   end function BinomialInteger

!ccccccccccccccc

   real (dp) function Binom(self, i, j)
     class (Model), intent(in) :: self
     integer      , intent(in) :: i, j
     integer                   :: k

    select case(j)
    case(0)  ;  Binom = 1
    case(:-1);  Binom = 0
    case default

      Binom = i

      do k = 1, j - 1
        Binom = Binom * (i - k)
      end do

      Binom = Binom/self%fact(j)

    end select

   end function Binom

!ccccccccccccccc

   real (dp) function CModCoef(self, n, i, m)
     class (Model), intent(in) :: self
     integer, intent(in) :: n, i, m
     integer             :: j, k
     real (dp)           :: CModCoef2

     CModCoef = 0; if ( i < 0 .or. i > n .or. m < 0 .or. m > 3 * i ) return

     do j = ceiling(m/3.), min(i,m)

       CModCoef2 = 0

       do k = ceiling( (m - j)/2. ), min(m - j, j)
         CModCoef2 = CModCoef2 + 1.5_dp**k/( self%Listfact(j - k) * &
                     self%Listfact(m - j - k) * self%Listfact(2 * k - m + j) )
       end do

       CModCoef = CModCoef + 3**j * CModCoef2/self%Listfact(i - j)

     end do

     CModCoef = CModCoef * self%Listfact(n)/self%Listfact(n - i)

   end function CModCoef

!ccccccccccccccc

  recursive function fact(self, i) result(prod)
    class (Model), intent(in) :: self
    integer      , intent(in) :: i
    real (dp)                 :: prod

    if ( i >= 0 .and. i <= 3 * self%Nmax + 8) then
      prod = self%Listfact(i); return
    end if

    select case(i)
    case(:-1);     prod = 0
    case(0:1);     prod = 1
    case default;  prod = i * self%fact(i - 1)
    end select

  end function fact

!ccccccccccccccc

  real (dp) function Model2D(c, modList, k1, k2, l1, l2)
    real (dp)   , dimension(:,:), intent(in) :: c
    type (Model), dimension(:,:), intent(in) :: modList
    integer                     , intent(in) :: k1, k2
    real (dp)                   , intent(in) :: l1, l2
    integer                                  :: i, j, k, l, clen

    Model2D = 0; clen = size(c,1)

    do i = 1, clen
      do j = 1, clen
        do k = 1, clen
          do l = 1, clen
            Model2D = Model2D + c(i,j) * c(k,l) * modList(i,k)%ShapeFun(k1, l1) * &
                                                  modList(j,l)%ShapeFun(k2, l2)
          end do
        end do
      end do
    end do

  end function Model2D

!ccccccccccccccc

  real (dp) function BreitModel2D(width, c, modList, k1, k2, l1, l2)
    real (dp)   , dimension(:,:), intent(in) :: c
    type (Model), dimension(:,:), intent(in) :: modList
    real (dp)                   , intent(in) :: l1, l2, width
    integer                     , intent(in) :: k1, k2
    integer                                  :: clen, i, j, k, l

    BreitModel2D = 0; clen = size(c,1)

    do i = 1, clen
      do j = 1, clen
        do k = 1, clen
          do l = 1, clen
            BreitModel2D = BreitModel2D + modList(i,k)%BreitModel(width, k1, l1) * &
                        c(i,j) * c(k,l) * modList(j,l)%BreitModel(width, k2, l2)
          end do
        end do
      end do
    end do

  end function BreitModel2D

!ccccccccccccccc

  real (dp) function SemiBreitModel2D(width, c, clen, modList, k1, k2, l1, l2)
    real (dp)   ,    dimension(clen,clen), intent(in) :: c
    type (Model),    dimension(clen,clen), intent(in) :: modList
    real (dp)                            , intent(in) :: l1, l2, width
    integer                              , intent(in) :: clen, k1, k2
    integer                                           :: i, j, k, l

    SemiBreitModel2D = 0

    do i = 1, clen
      do j = 1, clen
        do k = 1, clen
          do l = 1, clen
            SemiBreitModel2D = SemiBreitModel2D + modList(i,k)%ShapeFun(k1, l1) * &
                               c(i,j) * c(k,l) * modList(j,l)%BreitModel(width, k2, l2)
          end do
        end do
      end do
    end do

  end  function SemiBreitModel2D

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

end module ModelClass

!ccccccccccccccccc
