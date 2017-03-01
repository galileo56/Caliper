
module CumulantClass
  use ProfilesClass; use SingularClass; use ModelClass; use MassiveNSClass
  use constants, only: dp, prec, d1mach; use adapt; use QuadPack, only: qags; implicit none;  private

!ccccccccccccccc

  type, public                            :: CumulantMassless
    private

    class (SingularMassless), allocatable :: Sing
    class (Profiles)        , allocatable :: Prof
    real (dp)                             :: width, Q, muH
    real (dp), dimension(2)               :: tScen

    contains

    final                                 :: delete_object
    procedure, pass(self), private        :: SetScales
    procedure, pass(self), public         :: ListDistPiece, ListDist, ListBin, BinPiece, &
                                             Bin, ListBinPiece, FindOrigin
  end type CumulantMassless

!ccccccccccccccc

  type, extends (CumulantMassless), public :: CumulantMass
    private
    real (dp)                              :: xiJ, xiB
    real (dp)                              :: muM
    type (MassiveScales)                   :: MassNS
  end type CumulantMass

!ccccccccccccccc

  interface CumulantMassless
    module procedure InCum
  end interface CumulantMassless

!ccccccccccccccc

  interface CumulantMass
    module procedure InCumMass
  end interface CumulantMass

  contains

!ccccccccccccccc

 subroutine delete_object(this)
   type (CumulantMassless) :: this
     if ( allocated(this%Sing ) ) deallocate(this%Sing)
     if ( allocated(this%Prof ) ) deallocate(this%Prof)
  end subroutine delete_object

!ccccccccccccccc

  type (CumulantMassless) function InCum(Prof, Sing)
    type (ProfilesMassless), intent(in) :: Prof
    type (SingularScales)  , intent(in) :: Sing

    InCum%Q = Prof%Energy();  InCum%muH = Prof%HardScale(); InCum%width = 0

    allocate( ProfilesMassless :: InCum%Prof ); allocate( SingularScales :: InCum%Sing )

    Select Type (selector => InCum%Prof)
    type is (ProfilesMassless); selector = Prof
    end select

    Select Type (selector => InCum%Sing)
    type is (SingularScales); selector = Sing
    end select

    call InCum%Sing%SetHard( InCum%Q, InCum%muH )

    InCum%tScen = [-2._dp, -2._dp]

  end function InCum

!ccccccccccccccc

  type (CumulantMass) function InCumMass(Prof, Sing, xiJ, xiB, width)
    type (ProfilesPythia)    , intent(in) :: Prof
    type (SingularMassScales), intent(in) :: Sing
    real (dp)                , intent(in) :: xiJ, xiB
    real (dp), optional      , intent(in) :: width

    InCumMass%width = 0   ; if ( present(width) ) InCumMass%width = width
    InCumMass%xiJ = xiJ   ; InCumMass%xiB = xiB;  InCumMass%Q = Prof%Energy()
    InCumMass%muM = Prof%MassScale(); InCumMass%muH = Prof%HardScale()

    allocate( ProfilesPythia      :: InCumMass%Prof )
    allocate( SingularMassScales  :: InCumMass%Sing )

    Select Type (selector => InCumMass%Prof)
    type is (ProfilesPythia); selector = Prof
    end select

    InCumMass%tScen = Prof%Sectors()

    Select Type (selector => InCumMass%Sing)
    type is (SingularMassScales)
      selector = Sing
      call selector%SetHard( InCumMass%Q, InCumMass%muH )
      call selector%SetHardMass( InCumMass%muM )
    end select

  end function InCumMass

!ccccccccccccccc

  subroutine SetScales(self, muNS, muJ, muS, R, Rmass)
    class (CumulantMassless), intent(inout) :: self
    real (dp)               , intent(in   ) :: muNS, muJ, muS, R, Rmass

    select type (self)
    type is (CumulantMassless)
      call self%Sing%SetMat(muJ, muS); call self%Sing%SetRunning(muJ, muS, R, muJ)
    type is (CumulantMass)
      select type (sing => self%Sing)
      class is (SingularMassScales)
        call sing%SetEverything(muNS, self%Q, self%muH, muJ, muS, R, Rmass, self%muM, self%width)
        self%MassNS = Sing%NSScales()
      end select
      call self%Sing%SetRunning(muJ, muS, R, muJ)
      call self%Sing%SetMat(muJ, muS, self%xiJ, self%xiB)
    end select

  end subroutine SetScales

!ccccccccccccccc

  function ListDist(self, terms, cum, Mod, setup, gap, space, order, R0, mu0, &
                         delta0, h, tauList) result(crossList)
    class (CumulantMassless), intent(inout) :: self
    real (dp), dimension(:) , intent(in   ) :: tauList
    type (Model)            , intent(in   ) :: Mod
    character (len = *)     , intent(in   ) :: setup, space, gap, cum, terms
    real (dp)               , intent(in   ) :: R0, mu0, delta0, h
    integer                 , intent(in   ) :: order
    real (dp), dimension( size(tauList) )   :: crossList
    integer                                 :: i

    do i = 1, size(tauList)
      crossList(i) = self%Bin( terms, cum, Mod, setup, gap, space, order, R0, &
      mu0, delta0, h, 0, tauList(i) )
    end do

  end function ListDist

!ccccccccccccccc

  function ListBin(self, terms, cum, Mod, setup, gap, space, order, R0, mu0, &
                         delta0, h, tauList) result(crossList)
    class (CumulantMassless) , intent(inout) :: self
    real (dp), dimension(:,:), intent(in   ) :: tauList
    type (Model)             , intent(in   ) :: Mod
    character (len = *)      , intent(in   ) :: setup, space, gap, cum, terms
    real (dp)                , intent(in   ) :: R0, mu0, delta0, h
    integer                  , intent(in   ) :: order
    real (dp), dimension( size(tauList,2) )  :: crossList
    integer                                  :: i

    do i = 1, size(tauList,2)
      crossList(i) = self%Bin( terms, cum, Mod, setup, gap, space, order, R0, &
      mu0, delta0, h, 0, tauList(1,i), tauList(2,i) )
    end do

  end function ListBin

!ccccccccccccccc

  function ListDistPiece(self, terms, cum, ModList, setup, gap, space, order, R0, mu0, &
                         delta0, h, tauList) result(crossList)
    class (CumulantMassless)  , intent(inout) :: self
    real (dp), dimension(:)   , intent(in   ) :: tauList
    type (Model), dimension(:), intent(in   ) :: ModList
    character (len = *)       , intent(in   ) :: space, gap, cum, terms, setup
    real (dp)                 , intent(in   ) :: R0, mu0, delta0, h
    integer                   , intent(in   ) :: order

    real (dp), dimension( size(ModList), size(tauList) ) :: crossList
    integer                                              :: i

    do i = 1, size(tauList)
      crossList(:,i) = self%BinPiece( terms, cum, ModList, setup, gap, space, order, R0, &
      mu0, delta0, h, 0, tauList(i) )
    end do

  end function ListDistPiece

!ccccccccccccccc

  function ListBinPiece(self, terms, cum, ModList, setup, gap, space, order, R0, mu0, &
                         delta0, h, tauList) result(crossList)
    class (CumulantMassless)  , intent(inout) :: self
    real (dp), dimension(:,:) , intent(in   ) :: tauList
    type (Model), dimension(:), intent(in   ) :: ModList
    character (len = *)       , intent(in   ) :: space, gap, cum, terms, setup
    real (dp)                 , intent(in   ) :: R0, mu0, delta0, h
    integer                   , intent(in   ) :: order
    integer                                   :: i

    real (dp), dimension( size(ModList), size(tauList,2) ) :: crossList

    do i = 1, size(tauList,2)
      crossList(:,i) = self%BinPiece( terms, cum, ModList, setup, gap, space, order, R0, &
      mu0, delta0, h, 0, tauList(1,i), tauList(2,i) )
    end do

  end function ListBinPiece

!ccccccccccccccc

  recursive function BinPiece(self, terms, cum, ModList, setup, gap, space, order, &
                              R0, mu0, delta0, h, pow, t0, t1) result(resList)
    class (CumulantMassless)  , intent(inout) :: self
    real (dp)                 , intent(in   ) :: t0
    real (dp)       , optional, intent(in   ) :: t1
    type (Model), dimension(:), intent(in   ) :: ModList
    character (len = *)       , intent(in   ) :: space, gap, cum, terms, setup
    real (dp)                 , intent(in   ) :: R0, mu0, delta0, h
    integer                   , intent(in   ) :: order
    integer, optional         , intent(in   ) :: pow
    real (dp)   , dimension( size(ModList) )  :: resList
    real (dp)                                 :: muNS, muH, muJ, muS, R, tmin
    integer                                   :: i, power, run
    character (len = 4)                       :: cumul

    run = self%sing%runOrd();  tmin = self%FindOrigin(gap, order, R0, mu0, delta0, h)

    if ( .not. present(pow) ) then;  power = 0;  else; power = pow;  end if

    if ( .not. present(t1) ) then
      power = 0; cumul = cum; resList = CrossSing(t0); return
    end if

    if ( self%width < d1mach(1) .and. t0 < tmin  .and. cum(:3) /= 'cum') then

      resList = self%BinPiece(terms, cum, ModList, setup, gap, space, order, &
                          R0, mu0, delta0, h, pow, tmin, t1);  return

    else if ( t0 < self%tScen(1) .and. t1 > self%tScen(1)) then

      resList = self%BinPiece( terms, cum, ModList, setup, gap, space, order, R0, mu0, delta0, &
      h, pow, t0, self%tScen(1) - prec ) + self%BinPiece( terms, cum, ModList, setup, gap, &
      space, order, R0, mu0, delta0, h, pow, self%tScen(1) + prec , t1)

      return

    else if ( t0 < self%tScen(2) .and. t1 > self%tScen(2) ) then

      resList = self%BinPiece( terms, cum, ModList, setup, gap, space, order, R0, mu0, delta0, &
      h, pow, t0, self%tScen(2) - prec ) + self%BinPiece( terms, cum, ModList, setup, gap, &
      space, order, R0, mu0, delta0, h, pow, self%tScen(2) + prec , t1)

      return

    end if

    if ( cum(:8) == 'MidPoint' ) then

      cumul = 'cum';  call self%Prof%Scales( (t1 + t0)/2, muNS, muH, muJ, muS, R )
      call self%SetScales(muNS, muJ, muS, R, muJ); power = 0; resList = CrossProf(t0, t1)

    else if ( cum(:3) == 'cum' ) then
      power = 0; cumul = 'cum';  resList = CrossSing(t1) - CrossSing(t0)
    else if ( cum(:9) == 'Integrate' ) then

      do i = 1, size(ModList)
        resList(i) = self%Bin(terms, 'Integrate', ModList(i), 'Model', gap, space, &
        order, R0, mu0, delta0, h, pow, t0, t1)
      end do
    end if

    contains

!ccccccccccccccc

    function CrossSing(tau) result(list1D)
      real (dp)                 , intent(in) :: tau
      real (dp), dimension( size(ModList) )  :: list1D

      call self%Prof%Scales(tau, muNS, muH, muJ, muS, R)
      call self%SetScales(muNS, muJ, muS, R, muJ)

      list1D = CrossProf(tau) * tau**power

    end function CrossSing

!ccccccccccccccc

    function CrossProf(tau, tau2) result(list)
      real (dp)                , intent(in) :: tau
      real (dp),       optional, intent(in) :: tau2
      real (dp), dimension( size(ModList) ) :: list

      list = 0

      if ( terms(:3) /= 'Non' ) then

        if (self%width > 0) then
          select type ( sing => self%Sing )
          type is (SingularMassScales)

            if ( present(tau2) ) then
              list = Sing%SingleSingWidth(ModList, setup, gap, space, 'cum', order, &
              R0, mu0, delta0, h, tau, tau2)
            else
              list = Sing%SingleSingWidth(ModList, setup, gap, space, cumul, order, &
              R0, mu0, delta0, h, tau)
            end if
          end select

        else

          if ( present(tau2) ) then
            list = self%Sing%SingleSing(ModList, gap, space, 'cum', order, &
            R0, mu0, delta0, h, tau, tau2)
          else
            list = self%Sing%SingleSing(ModList, gap, space, cumul, order, &
            R0, mu0, delta0, h, tau)
          end if
        end if
      end if

      select type ( sing => self%Sing )
      type is (SingularMassScales)

        if ( terms(:4) /= 'Sing' .and. terms(:7) /= 'NonSCET' ) then

          if ( present(tau2) ) then
            list = list + Sing%NonDist(ModList, space, gap, 'cum', order, R0, mu0, delta0, h, tau, tau2)
          else
            list = list + Sing%NonDist(ModList, space, gap, cumul, order, R0, mu0, delta0, h, tau)
          end if
        end if

        if ( terms(:3) == 'All' .or. terms(:7) == 'NonSCET'  .or. terms(:7) == 'NonSing' ) then

          select type (self)
          type is (CumulantMass)

             if ( present(tau2) ) then
              list = list + self%MassNS%NSMass(ModList, setup, gap, 'cum', &
              order, run, R0, mu0, delta0, h, tau, tau2)
            else
              list = list + self%MassNS%NSMass(ModList, setup, gap, cumul, &
              order, run,  R0, mu0, delta0, h, tau)
            end if

          end select

        end if
      end select

    end function CrossProf

  end function BinPiece

!ccccccccccccccc

  recursive real (dp) function Bin(self, terms, cum, Mod, setupIn, gap, space, &
                             order, R0, mu0, delta0, h, pow, t0, t1) result(res)
    class (CumulantMassless), intent(inout) :: self
    real (dp)               , intent(in) :: t0
    real (dp), optional     , intent(in) :: t1
    type (Model)            , intent(in) :: Mod
    character (len = *)     , intent(in) :: setupIn, space, gap, cum, terms
    real (dp)               , intent(in) :: R0, mu0, delta0, h
    integer                 , intent(in) :: order
    integer, optional       , intent(in) :: pow
    real (dp), dimension(2)              :: pScen
    integer                              :: power, run, shapeCum, neval, ier
    character (len = 4)                  :: cumul
    character (len = 15)                 :: setup
    logical                              :: dobsing
    real (dp)                            :: muNS, muH, muJ, muS, R, abserr, p, &
                                            p2, tmin, res0, res1, pmin

    run = self%sing%runOrd();  setup = setupIn
    tmin = self%FindOrigin(gap, order, R0, mu0, delta0, h)

    if ( .not. present(pow) ) then;  power = 0;  else; power = pow;  end if

    if ( .not. present(t1) .and. setupIn(:8) /= 'Profiles' ) then
      power = 0; cumul = cum; res = CrossSing(t0); return
    end if

    if ( present(t1) .and. t1 <= t0 ) return
    if ( .not. present(t1) .and. t0 < tmin ) return

    if ( setupIn(:8) == 'Profiles' ) then

      power = 0; setup = 'NoModel'; cumul = 'diff'; p = self%Q * t0; p2 = p
      pScen = self%Q * self%tScen
      pmin = self%Q * tmin; shapeCum = 0; if ( cum(:3) == 'cum') shapeCum = -1
      if ( present(t1) ) p2 = self%Q * t1
      if ( present(t1) .and. t1 <= tmin ) p = p2
      dobsing = .not. present(t1) .or. ( present(t1) .and. abs(p2 - p) <= d1mach(1) )

      if ( p2 < pScen(1) .or. pmin > pScen(1) ) then
        call qags( ProfCrossSing, pmin, p2, prec, prec, res, abserr, neval, ier )
      else if ( p2 > pScen(1) .and. p2 < pScen(2) ) then
        call qags( ProfCrossSing, pmin, pScen(1), prec, prec, res, abserr, neval, ier )
        call qags( ProfCrossSing, pScen(1), p2, prec, prec, res0, abserr, neval, ier )
        res = res + res0
      else if ( p2 > pScen(2) ) then
        call qags( ProfCrossSing, pmin, pScen(1), prec, prec, res, abserr, neval, ier )
        call qags( ProfCrossSing, pScen(1), pScen(2), prec, prec, res0, abserr, neval, ier )
        call qags( ProfCrossSing, pScen(2), p2, prec, prec, res1, abserr, neval, ier )
        res = res + res0 + res1
      end if

      res = res * self%Q**shapeCum; return

    end if

    if ( self%width < d1mach(1) .and. t0 < tmin  .and. cum(:3) /= 'cum') then

      res = self%Bin(terms, cum, Mod, setup, gap, space, order, R0, mu0, &
      delta0, h, pow, tmin, t1);  return

    else if ( t0 < self%tScen(1) .and. t1 > self%tScen(1)) then

      res = self%Bin( terms, cum, Mod, setup, gap, space, order, R0, mu0,delta0, &
      h, pow, t0, self%tScen(1) - prec ) + self%Bin( terms, cum, Mod, setup, &
      gap, space, order, R0, mu0, delta0, h, pow, self%tScen(1) + prec , t1)

      return

    else if ( t0 < self%tScen(2) .and. t1 > self%tScen(2) ) then

      res = self%Bin( terms, cum, Mod, setup, gap, space, order, R0, mu0,delta0, &
      h, pow, t0, self%tScen(2) - prec ) + self%Bin( terms, cum, Mod, setup, &
      gap, space, order, R0, mu0, delta0, h, pow, self%tScen(2) + prec , t1)

      return

    end if

    if ( cum(:8) == 'MidPoint' ) then

      cumul = 'cum';  call self%Prof%Scales( (t1 + t0)/2, muNS, muH, muJ, muS, R )
      call self%SetScales(muNS, muJ, muS, R, muJ); power = 0;  res = CrossProf(t0, t1)

    else if ( cum(:3) == 'cum' ) then
      power = 0; cumul = 'cum';  res = CrossSing(t1) - CrossSing(t0)
    else if ( cum(:9) == 'Integrate' ) then
      cumul = 'diff'
      call DAdapt(CrossSing, t0, t1, 1, prec, prec, res, abserr)
    end if

    contains

!ccccccccccccccc

    real (dp) function ProfCrossSing(l)
      real (dp), intent(in) :: l

      if ( dobsing ) then
        ProfCrossSing = CrossSing(l/self%Q) * Mod%ShapeFun(shapeCum, p2 - l)
      else
        ProfCrossSing = CrossSing(l/self%Q) * Mod%ShapeFun(shapeCum, p - l, p2 - l)
      end if

    end function ProfCrossSing

!ccccccccccccccc

    real (dp) function CrossSing(tau)
      real (dp), intent(in) :: tau

      call self%Prof%Scales(tau, muNS, muH, muJ, muS, R)
      call self%SetScales(muNS, muJ, muS, R, muJ)

      CrossSing = CrossProf(tau) * tau**power

    end function CrossSing

!ccccccccccccccc

    real (dp) function CrossProf(tau, tau2)
      real (dp)          , intent(in) :: tau
      real (dp), optional, intent(in) :: tau2

      CrossProf = 0

      if ( terms(:3) /= 'Non' ) then

        if (self%width > 0) then
          select type ( sing => self%Sing )
          type is (SingularMassScales)

            if ( present(tau2) ) then
              CrossProf = Sing%SingleSingWidth(Mod, setup, gap, space, 'cum', order, &
              R0, mu0, delta0, h, tau, tau2)
            else
              CrossProf = Sing%SingleSingWidth(Mod, setup, gap, space, cumul, order, &
              R0, mu0, delta0, h, tau)
            end if
          end select

        else

          if ( present(tau2) ) then
            CrossProf = self%Sing%SingleSing(Mod, setup, gap, space, 'cum', order, &
                                            R0, mu0, delta0, h, tau, tau2)
          else
            CrossProf = self%Sing%SingleSing(Mod, setup, gap, space, cumul, order, &
                                             R0, mu0, delta0, h, tau)
          end if
        end if
      end if

      select type ( sing => self%Sing )
      type is (SingularMassScales)

        if ( terms(:4) /= 'Sing' .and. terms(:7) /= 'NonSCET' ) then

          if ( present(tau2) ) then
            CrossProf = CrossProf + Sing%NonDist(Mod, setup, gap, space, 'cum', order,&
                                                 R0, mu0, delta0, h, tau, tau2)
          else
            CrossProf = CrossProf + Sing%NonDist(Mod, setup, gap, space, cumul, order,&
                                                 R0, mu0, delta0, h, tau)
          end if
        end if

        if ( terms(:3) == 'All' .or. terms(:7) == 'NonSCET'  .or. terms(:7) == 'NonSing' ) then

          select type (self)
          type is (CumulantMass)

             if ( present(tau2) ) then
              CrossProf = CrossProf + self%MassNS%NSMass(Mod, setup, gap, 'cum', &
                                      order, run, R0, mu0, delta0, h, tau, tau2)
            else
              CrossProf = CrossProf + self%MassNS%NSMass(Mod, setup, gap, cumul, &
                                      order, run, R0, mu0, delta0, h, tau)
            end if

          end select
        end if
      end select

    end function CrossProf

  end function Bin

!ccccccccccccccc

  real (dp) function FindOrigin(self, gap, order, R0, mu0, delta0, h)
    class (CumulantMassless), intent(inout) :: self
    character (len = *)     , intent(in   ) :: gap
    real (dp)               , intent(in   ) :: R0, mu0, delta0, h
    integer                 , intent(in   ) :: order
    real (dp)                               :: muNS, muH, muJ, muS, R, a

    select type (self)
    type is (CumulantMassless)
      if ( gap(:6) /= 'thrust' .and. gap(:6) /= 'Cparam' ) then
        FindOrigin = 0; return
      end if
    end select

      FindOrigin = 0

    do

      a = Origin(FindOrigin); if ( abs(a - FindOrigin) <= 1.e-10_dp ) exit

      FindOrigin = a

    end do

!ccccccccccccccc

    contains

    real (dp) function Origin(tau)
      real (dp), intent(in) :: tau

      call self%Prof%Scales(tau, muNS, muH, muJ, muS, R)
      call self%SetScales(muNS, muJ, muS, R, muJ)

      Origin = self%Sing%ESmin(order, gap, delta0, R0, mu0, h)

    end function Origin


  end function FindOrigin

!ccccccccccccccc

end module CumulantClass
