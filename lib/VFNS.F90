module VFNSMSRClass
  use AnomDimClass;  use AlphaClass;  use RunningClass
  use Constants, only: dp, Pi, ExpEuler, d1mach, prec
  use QuadPack , only: qags         ; implicit none   ;  private

  public :: VFNSMSR

!ccccccccccccccc

  type, public                   :: VFNSMSR
    private

    integer                      :: nf, nl
    type (Running), dimension(2) :: Run
    real (dp)                    :: mQ, mL

    contains

    procedure, pass(self), public :: MSRMass!setMtop, setMbottom, setMcharm

  end type VFNSMSR

!ccccccccccccccc

  interface VFNSMSR
    module procedure  InitMSR
  end interface VFNSMSR

contains

!ccccccccccccccc

   type (VFNSMSR) function InitMSR(Run)
     type (Running), dimension(2), intent(in) :: run

     InitMSR%run = run
     InitMSR%nf  = Run(2)%numFlav()   ; InitMSR%nl = Run(1)%numFlav()
     InitMSR%mQ  = Run(2)%scales('mH'); InitMSR%mL = Run(1)%scales('mH')

   end function InitMSR

!ccccccccccccccc

  recursive real (dp) function MSRMass(self, type, order, R, lambda, method) result(res)
    class (VFNSMSR)                    , intent(in) :: self
    real (dp)                          , intent(in) :: R
    real (dp), optional                , intent(in) :: lambda
    character (len = *)                , intent(in) :: type
    character (len = *), optional      , intent(in) :: method
    integer                            , intent(in) :: order

    if (R >= self%mL) then
    end if

    ! res = self%

  end function MSRMass

end module VFNSMSRClass
