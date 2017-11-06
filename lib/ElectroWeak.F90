
module ElectroWeakClass
  use Constants, only: dp;  implicit none
  private
  public                          :: charge

!ccccccccccccccc

  type, public                    :: ElectroWeak
   private
   real (dp)                      :: ve, ae, vUp, aUp, vDown, mZ, gammaZ

   contains

   procedure, pass(self), public :: EWfactors

  end type ElectroWeak

!ccccccccccccccc

  interface ElectroWeak
    module procedure InitElectroWeak
  end interface ElectroWeak

  contains

!ccccccccccccccc

   type (ElectroWeak) function InitElectroWeak(mZ, gammaZ, sin2ThetaW)
    real (dp), intent(in)           :: gammaZ, sin2ThetaW
    real (dp), intent(in), optional :: mZ
    real (dp)                       :: sinWein, CosWein, sCWein

    if ( .not. present(mZ) ) return

    InitElectroWeak%mZ = mZ; InitElectroWeak%gammaZ = gammaZ

    sinWein = sqrt(sin2ThetaW) ; CosWein = sqrt( 1 - sin2ThetaW )
    sCWein = 2 * sinWein * CosWein

    InitElectroWeak%ve = (2 * sin2ThetaW - 0.5_dp)/sCWein ;
    InitElectroWeak%ae = - 1/sCWein/2 ;

    InitElectroWeak%vUp = ( - 4 * sin2ThetaW/3 + 0.5_dp)/sCWein ;
    InitElectroWeak%aUp = 1/sCWein/2 ;

    InitElectroWeak%vDown = (2 * sin2ThetaW/3 - 0.5_dp)/sCWein ;

   end function InitElectroWeak

!ccccccccccccccc

  function EWfactors(self, nf, Q) result(EW)
    class (ElectroWeak), intent(in) :: self
    Integer            , intent(in) :: nf
    real (dp)          , intent(in) :: Q
    real (dp), dimension(2)         :: eW
    real (dp)                       :: Q2, Q4, prop, prop2, ratio, Qquark, vt, electron

    Q2 = Q**2; Q4 = Q2**2; prop = Q2 - self%mZ**2; prop2 = prop**2
    ratio = prop2 + Q4 * (self%gammaZ/self%mZ)**2

    select case(nf)
      case(2,4,6) ;  Qquark = 2._dp/3;   vt = self%vUp
      case(1,3,5) ;  Qquark = - 1._dp/3; vt = self%vDown
      case default;  Qquark = 0; vt = 0;
    end select

    electron = Q4 * ( self%ve**2 + self%ae**2 )

    EW(1) = Qquark**2 + ( electron * vt**2 - 2 * Q2 * self%ve * vt * Qquark * prop )/ratio
    EW(2) = electron * self%aUp**2/ratio

  end function EWfactors

!ccccccccccccccc

  real (dp) function charge(nf)
    integer, intent(in) :: nf
     select case(nf)
      case(2,4,6)
       charge = 2._dp/3
    case(1,3,5)
       charge = 1
    case default
       charge = - 1._dp/3
    end select
  end function charge

!ccccccccccccccc

end module ElectroWeakClass
