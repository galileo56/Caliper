
module KernelsClass
  use DeriGamma; use Constants, only: dp, d1mach;  implicit none
  private

  public :: KernelsCombine, KernelMatrixList, LogMatrixSum, KerList, DerAdd1, &
            NGLSum, KernelsSum, AddBinomial, KernelWidth, SumKernelWidth,     &
            DerSubtract1

!ccccccccccccccc

  type, public                           ::  Kernels
    private
    integer                              :: n
    real (dp), allocatable, dimension(:) :: res, derinvlist
    real (dp)                            :: w

   contains

   procedure, pass(self), public         :: KernelList, dimList, DerList
   final                                 :: delete_object

  end type Kernels

!ccccccccccccccc

  interface Kernels
    module procedure InitKernels
  end interface Kernels

  contains

!ccccccccccccccc

 subroutine delete_object(this)
   type (Kernels) :: this
     if ( allocated(this%res       ) ) deallocate(this%res       )
     if ( allocated(this%DerInvList) ) deallocate(this%DerInvList)
  end subroutine delete_object

!ccccccccccccccc

   type (Kernels) function InitKernels(n, width, w, mu, p)
     integer  , intent(in)           :: n
     real (dp), intent(in)           :: w
     real (dp), intent(in), optional :: width, mu, p

     real (dp)                       :: lg, x, y
     real (dp), dimension(3,0:n)     :: DerInvList

     integer                         :: i, j

     allocate( InitKernels%res(0:n) )        ; InitKernels%res        = 0
     allocate( InitKernels%DerInvList(0:n) ) ; InitKernels%DerInvList = 0

     InitKernels%n = n; InitKernels%w = w

     if ( .not. present(width) .or. width <= d1mach(1) ) then

       InitKernels%DerInvList = GammaDerComputer(-1, n, w)

       if ( present(p) .and. present(mu) ) x = p

     else

       DerInvList(1,:) = GammaDerComputer(1, n, 2 + w)
       if ( .not. present(p) ) InitKernels%DerInvList = DerInvList(1,:)

       if ( present(p) ) then

       y = 0.5_dp + Atan(p/width)/Pi

       DerInvList(2:,:)       = KernelWidth(n, w, y)
       InitKernels%DerInvList = SumKernelWidth(DerInvList, y)

       if ( present(mu) ) x = Sqrt(width**2 + p**2)

       end if
     end if

     if ( present(mu) ) then

       y = mu * ExpEuler/x; lg = log(y)

       do i = 0, n
         do j = 0, i
           InitKernels%res(i) = InitKernels%res(i) + binomial(i,j) * lg**j * &
                                                     InitKernels%DerInvList(i - j)
         end do
       end do

       InitKernels%res = y**w * InitKernels%res/x

     end if

   end function InitKernels

!ccccccccccccccc

  pure function SumKernelWidth(DerInvList, y) result(InvList)
    real (dp)                , intent(in) :: y
    real (dp), dimension(:,:), intent(in) :: DerInvList

    integer                                    :: n, i, j
    real (dp), dimension( size(DerInvList,2) ) :: SumList
    real (dp), dimension( size(DerInvList,2) ) :: InvList

    n = size(DerInvList,2) - 1; InvList = 0; SumList = 0

    do i = 0, n
      do j = 0, i
        SumList(i + 1) = SumList(i + 1) + (-1)**j * binomial(i,j) *  &
                                        DerInvList(2,i - j + 1) * DerInvList(3,j + 1)
      end do
    end do

       do i = 0, n
         do j = 0, i
           InvList(i + 1) = InvList(i + 1) + y**(j + 1) * binomial(i,j) * &
                                  DerInvList(1,i - j + 1) * SumList(j + 1)
         end do
       end do

   end function SumKernelWidth

!ccccccccccccccc

  function KernelWidth(n, w, y) result(DerInvList)
    real (dp)             , intent(in) :: y, w
    integer               , intent(in) :: n
    real (dp), dimension(2:3,0:n)      :: DerInvList
    real (dp)                          :: t

    t = (1 + w) * y

    DerInvList(2,:) = GammaDerComputer(2, n, 1 + t)
    DerInvList(3,:) = GammaDerComputer(2, n, 1 - t)

   end function KernelWidth

!ccccccccccccccc

  function KernelList(self) result(beta)
    class (Kernels), intent(in) :: self
    real (dp)                   :: beta(0:self%n)

    beta = self%res

   end function KernelList

!ccccccccccccccc

  function DerList(self) result(beta)
    class (Kernels), intent(in) :: self
    real (dp)                   :: beta(0:self%n)

    beta = self%derinvlist

   end function DerList

!ccccccccccccccc

  pure function KernelMatrixList(DerInvList, Alist) result(beta)
    real (dp), intent(in), dimension(:)                 :: DerInvList, Alist
    real (dp)            , dimension(0:size(Alist) - 1) :: beta
    integer                                             :: i, j, n

    n = size(DerInvList) - 1; beta = 0

    do j = 0, n
      do i = j, n
        beta(j) = beta(j) + binomial(i,j) * AList(i + 1) * DerInvList(i - j + 1)
      end do
    end do

   end function KernelMatrixList

!ccccccccccccccc

  pure function KerList(DerInvList, Alist) result(beta)
    real (dp), intent(in), dimension(:)                 :: DerInvList, Alist
    real (dp)            , dimension(0:size(Alist) - 1) :: beta
    integer                                             :: i, j, n

    n = size(Alist) - 1; beta = 0

    do i = 0, n
      do j = 0, i
        beta(i) = beta(i) + binomial(i,j) * AList(j + 1) * DerInvList(i - j + 1)
      end do
    end do

   end function KerList

!ccccccccccccccc

  pure function AddBinomial(DerInvList, Alist) result(beta)
    real (dp), intent(in), dimension(:)   :: Alist
    real (dp), intent(in), dimension(:,:) :: DerInvList
    integer                               :: i, j, n
    real (dp), dimension( 0:size(Alist,1) - 1,0:size(Alist) - 1 ) :: beta

    n = size(Alist) - 1; beta = 0

    do i = 0, n
      do j = i, n
        beta(:,i) = beta(:,i) + binomial(j,i) * AList(j + 1) * DerInvList(:,j - i + 1)
      end do
    end do

   end function AddBinomial

!ccccccccccccccc

  pure function NGLSum(SoftNGLList, NGL2loop) result(SoftNGL)
    real (dp), intent(in), dimension(:)          :: NGL2loop
    real (dp), intent(in), dimension(:)          :: SoftNGLList
    real (dp), dimension( 0:2 * size(NGL2loop) ) :: SoftNGL
    integer                                      :: i, j, isoft

    isoft = size(NGL2loop); SoftNGL = 0

    do i = 0, 2 * isoft
      do j = max( 1, ceiling(i/2.) ), isoft

        SoftNGL(i) = SoftNGL(i) + binomial(2 * j, i) * SoftNGLList(2*j - i + 1) * NGL2loop(j)

      end do
    end do

   end function NGLSum

!ccccccccccccccc

  pure real (dp) function LogMatrixSum(DerInvList, w, mu, p, width)
    real (dp), intent(in), dimension(:) :: DerInvList
    real (dp), intent(in)               :: w, mu, p
    real (dp), intent(in), optional     :: width
    real (dp)                           :: lg, y, x
    integer                             :: j, n

    x = p; if ( present(width) ) x = sqrt(p**2 + width**2)

    n = size(DerInvList) - 1;  y = mu/x * ExpEuler; lg = log(y)
    LogMatrixSum = DerInvList(1)

    do j = 1, n
      LogMatrixSum = LogMatrixSum + lg**j * DerInvList(j + 1)
    end do

    LogMatrixSum = LogMatrixSum * y**w/x

  end function LogMatrixSum

!ccccccccccccccc

   integer function DimList(self)
     class(Kernels), intent(in) :: self

     DimList = self%n

   end function DimList

!ccccccccccccccc

   pure function KernelsCombine(DerInvList1, DerInvList2, n) result(res)
     real (dp), intent(in), dimension(:) :: DerInvList1
     real (dp), intent(in), dimension(:) :: DerInvList2
     integer  , intent(in)               :: n
     real (dp), dimension(n)             :: res
     integer                             :: i, j, nn, n1, n2

     nn = 2 * n; n1 = size(DerInvList1) - 1 - nn; n2 = size(DerInvList2) - 1 - nn; res = 0

     do i = 2, nn, 2
       do j = 0, i
         res(i/2) = res(i/2) + binomial(i,j) * (-1)**j * DerInvList1(i - j + n1 + 1) * &
                                                         DerInvList2(j + n2 + 1)
       end do
     end do

   end function KernelsCombine

!ccccccccccccccc

   pure function DerAdd1(pow, list, w) result(res)
     real (dp), dimension(:)     , intent(in) :: list
     real (dp)                   , intent(in) :: w
     integer                     , intent(in) :: pow
     real (dp), dimension( 0:size(list) - 1 ) :: res
     real (dp)                                :: w2
     integer                                  :: i, n

     n = size(list) - 1; w2 = w; if (pow < 0) w2 = w + 1

     res(0) = w2 * list(1)

     do i = 1, n
       res(i) = i * list(i) + w2 * list(i + 1)
     end do

     if (pow < 0) res = - res

   end function DerAdd1

!ccccccccccccccc

   pure function DerSubtract1(pow, list, w) result(res)
     real (dp), dimension(:)     , intent(in) :: list
     real (dp)                   , intent(in) :: w
     integer                     , intent(in) :: pow
     real (dp), dimension( 0:size(list) - 1 ) :: res
     real (dp)                                :: w2, fact
     integer                                  :: i, j, n

     n = size(list) - 1 ; w2 = w - 1 ; if (pow < 0) w2 = w ; res = list/w2

     do i = 1, n
       fact = 1/w2
       do j = 1, i
         fact = - fact * (i - j + 1)/w2;  res(i) = res(i) + fact * list(i - j + 1)
       end do
     end do

     if (pow < 0) res = - res

   end function DerSubtract1

!ccccccccccccccc

   pure function KernelsSum(DerInvList1, DerInvList2) result(res)
     real (dp), intent(in), dimension(:) :: DerInvList1
     real (dp), intent(in), dimension(:) :: DerInvList2

     integer                             :: i, j, n, k

     real (dp), dimension( 0:size(DerInvList2) - 1 ):: res

     n = size(DerInvList2) - 1; k = size(DerInvList1) - 1 - n; res = 0

     do i = 0, n
       do j = 0, i
         res(i) = res(i) + binomial(i,j) * (-1)**j * DerInvList1(i + k - j + 1) * &
                                                     DerInvList2(j + 1)
       end do
     end do

   end function KernelsSum

!ccccccccccccccc

  function GammaDerComputer(pow, n, w) result(res)
   integer  , intent(in)     :: n, pow
   real (dp), intent(in)     :: w
   real (dp), dimension(0:n) :: res, PolyGammaList, DerInvList
   real (dp)                 :: w2
   integer                   :: i, j, power, NZ, IERR

   res = 0; PolyGammaList = 0; DerInvList = 0;

   w2 = w;  power = pow;  if (pow <  0) w2 = 4 - w;  if (pow == 2) power = - 1

   call DPSIFN(w2, 0, 1, n + 1, PolyGammaList, NZ, IERR)

   DerInvList(0) = Gamma(w2)**power

   do i = 1, n
     do j = 0, i - 1
       DerInvList(i) = DerInvList(i) + power * binomial(i - 1,j) * &
       (-1)**(i - j) * fact(i - 1 - j) * PolyGammaList(i - 1 - j) * DerInvList(j)
     end do
   end do

   if (pow < 0) then

     do i = 0, n
       do j = max(0,i - 4), i
         res(i) = res(i) + binomial(i,j) * Take4( i - j, -w ) * DerInvList(j)
       end do

       res(i) = (-1)**i * res(i)

     end do

   else

     res = DerInvList

   end if

   end function GammaDerComputer

!ccccccccccccccc

  pure recursive function fact(i) result(prod)
    integer, intent(in) :: i
    real (dp)           :: prod

    select case(i)
    case(:-1)
      prod = 0
    case(0:1)
      prod = 1
    case default
      prod = i * fact(i - 1)
    end select

  end function fact

!ccccccccccccccc

  pure real (dp) function Take4(n, w)
    real (dp), intent(in) :: w
    integer  , intent(in) :: n

    select case(n)
    case(0)
      Take4 = (3 + w) * (2 + w) * (1 + w) * w
    case(1)
      Take4 = 2 * (3 + 2 * w) * (1 + 3 * w + w**2)
    case(2)
      Take4 = 2 * (11 + 18 * w + 6 * w**2)
    case(3)
      Take4 = 12 * (3 + 2 * w)
    case(4)
      Take4 = 24
    case default
      Take4 = 0
    end select

  end function Take4

!ccccccccccccccc

  pure real (dp) function Binomial(i,j)
    integer, intent(in) :: i, j
    integer             :: k

    if (j < 0) then
      Binomial = 0
    elseif (j == 0) then
      Binomial = 1
    else
      Binomial = i
     do k = 1, j - 1
        Binomial = Binomial * (i - k)
      end do
      Binomial = Binomial/fact(j)
    end if

  end function Binomial

!ccccccccccccccc

end module KernelsClass
