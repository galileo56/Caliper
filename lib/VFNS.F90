module VFNSMSRClass
  use AnomDimClass;  use AlphaClass;  use RunningClass
  use Constants, only: dp, Pi, ExpEuler, d1mach, prec, pi2
  use QuadPack , only: qags         ; implicit none   ;  private

  public :: VFNSMSR

!ccccccccccccccc

  type, public                   :: VFNSMSR
    private

    integer                      :: nf, nl
    type (Running), dimension(2) :: Run
    type (AnomDim)               :: AnDim
    real (dp)                    :: mQ, mL, b1, beta0

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
     type (AnomDim)                           :: AnDim
     real (dp), dimension(0:4)                :: bHat

     InitMSR%run = run
     InitMSR%nf  = Run(2)%numFlav()     ; InitMSR%nl = Run(1)%numFlav()
     InitMSR%mQ  = Run(2)%scales('mH')  ; InitMSR%mL = Run(1)%scales('mH')
     AnDim       = Run(2)%adim()        ; bHat       = AnDim%betaQCD('bHat')
     InitMSR%b1  = bHat(1)              ; InitMSR%AnDim = AnDim
     bHat        = AnDim%betaQCD('beta'); InitMSR%beta0 = bHat(0)

   end function InitMSR

!ccccccccccccccc

  recursive real (dp) function MSRMass(self, type, order, R, lambda, method) result(res)
    class (VFNSMSR)                    , intent(in) :: self
    real (dp)                          , intent(in) :: R
    real (dp), optional                , intent(in) :: lambda
    character (len = *)                , intent(in) :: type
    character (len = *), optional      , intent(in) :: method
    integer                            , intent(in) :: order
    integer                                         :: neval, ier
    real (dp)                                       :: corr, abserr, lambdaQCD, &
    t0, t1

    if (R < self%mL) then
      res = self%MSRMass(type, order, self%mL, lambda, method) + &
      self%run(1)%MSREvol(type, order, R, lambda, method)
      return
    end if

    res = self%run(2)%MSREvol(type, order, R, lambda, method)

    if ( type(:7) == 'numeric' ) then

      call qags( InteNum, R/lambda, self%mQ/lambda, prec, prec, corr, &
      abserr, neval, ier )

      res = res + corr

    else if ( type == 'analytic' ) then

      t0 = - 2 * pi/self%beta0;  t1 = t0/self%run(2)%alphaQCD(self%mQ)
      t0 = t0/self%run(2)%alphaQCD(R)
      lambdaQCD = self%run(2)%lambdaQCD( self%run(2)%orders('runMass') )

      call qags( InteAn, t1, t0, prec, prec, corr, abserr, neval, ier )

      res = res + corr/4/self%beta0**2

    end if

  contains

!ccccccccccccccc

    real (dp) function InteNum(mu)
      real (dp), intent(in) :: mu

      InteNum = gammaRcharm2(mu) * self%run(2)%alphaQCD(self%mQ/mu)**2/16/Pi2

    end function InteNum

!ccccccccccccccc

    real (dp) function InteAn(t)
      real (dp), intent(in) :: t

      InteAn = Exp(-t) * (-t)**(- 2 - self%b1) * &
      gammaRcharm2(self%mQ/lambdaQCD * self%Andim%Gfun(self%run(2)%orders('runMass'),t) )

    end function InteAn

  end function MSRMass

end module VFNSMSRClass
