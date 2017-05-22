module VFNSMSRClass
  use AnomDimClass;  use AlphaClass;  use RunningClass
  use Constants, only: dp, Pi, ExpEuler, d1mach, prec, pi2
  use QuadPack , only: qags         ; implicit none   ;  private

  public :: VFNSMSR

!ccccccccccccccc

  type, public                   :: VFNSMSR
    private

    integer                      :: nf, nl, run
    type (Running), dimension(2) :: AlphaMass
    type (AnomDim)               :: AnDim
    real (dp)                    :: mH, mL, b1, beta0

    contains

    procedure, pass(self), public :: MSRMass!setMtop, setMbottom, setMcharm

  end type VFNSMSR

!ccccccccccccccc

  interface VFNSMSR
    module procedure  InitMSR
  end interface VFNSMSR

contains

!ccccccccccccccc

   type (VFNSMSR) function InitMSR(AlphaMass)
     type (Running), dimension(2), intent(in) :: AlphaMass
     type (AnomDim)                           :: AnDim
     real (dp), dimension(0:4)                :: bHat

     InitMSR%AlphaMass = AlphaMass      ; InitMSR%run = AlphaMass(2)%orders('runMass')
     InitMSR%nf  = AlphaMass(2)%numFlav()     ; InitMSR%nl = AlphaMass(1)%numFlav()
     InitMSR%mH  = AlphaMass(2)%scales('mH')  ; InitMSR%mL = AlphaMass(1)%scales('mH')
     AnDim       = AlphaMass(2)%adim()        ; bHat       = AnDim%betaQCD('bHat')
     InitMSR%b1  = bHat(1)              ; InitMSR%AnDim = AnDim
     bHat        = AnDim%betaQCD('beta'); InitMSR%beta0 = bHat(0)

   end function InitMSR

!ccccccccccccccc

  recursive real (dp) function MSRMass(self, type, order, R, lambda, method) result(res)
    class (VFNSMSR)              , intent(in) :: self
    real (dp)                    , intent(in) :: R
    real (dp), optional          , intent(in) :: lambda
    character (len = *)          , intent(in) :: type
    character (len = *), optional, intent(in) :: method
    real (dp)                  , dimension(4) :: a
    integer                      , intent(in) :: order
    integer                                   :: neval, ier, i
    real (dp)                                 :: corr, abserr, lambdaQCD, t0, t1, &
    delta, mu, Rstep, delta1, delta2, alpha2, matching, alphaM

    res = 0

    if (R < self%mL) then

      if ( type(:4) == 'MSRp' ) then
        alphaM = self%AlphaMass(2)%alphaQCD(self%mL)/Pi
        a = self%AlphaMass(1)%MSRMatching('charm'); i = min(order, 4)
        matching = dot_product( a(:i), PowList(alphaM, i) ) - 1
      else
        matching = - 1
      end if

      res = self%MSRMass(type, order, self%mL, lambda, method)  + &
      self%AlphaMass(1)%MSRMass(type, order, R, lambda, method) + &
      self%mL * matching
      return
    end if

    res = self%AlphaMass(2)%MSRMass(type, order, R, lambda, method)

    if ( method(:7) == 'numeric' ) then

      call qags( InteNum, R/lambda, self%mH/lambda, prec, prec, corr, &
      abserr, neval, ier )

      res = res + corr

    else if ( method(:8) == 'analytic' ) then

      t0 = - 2 * pi/self%beta0/lambda;  t1 = t0/self%AlphaMass(2)%alphaQCD(self%mH)
      t0 = t0/self%AlphaMass(2)%alphaQCD(R)
      lambdaQCD = self%AlphaMass(2)%lambdaQCD( self%run )

      call qags( InteAn, t1, t0, prec, prec, corr, abserr, neval, ier )

      res = res + lambdaQCD * corr/4/self%beta0**2

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( ( self%mH - R )/delta )  )
      delta = (self%mH - R)/neval; corr = 0; Rstep = self%mH + delta

      do ier = 1, neval

        Rstep = Rstep - delta; mu = (Rstep - delta/2)/lambda
        alpha2 = ( self%AlphaMass(2)%alphaQCD(mu)/Pi )**2

        delta1 = Rstep * alpha2 * deltaCharm2(self%mL/Rstep)

        delta2 = (Rstep - delta) * alpha2 * deltaCharm2( self%mL/(Rstep - delta )  )


        corr = corr + delta1 - delta2

      end do

      res = res + corr

    end if

  contains

!ccccccccccccccc

    real (dp) function InteNum(mu)
      real (dp), intent(in) :: mu

      InteNum = gammaRcharm2(self%mL/mu) * self%AlphaMass(2)%alphaQCD(mu)**2/16/Pi2

    end function InteNum

!ccccccccccccccc

    real (dp) function InteAn(t)
      real (dp), intent(in) :: t

      InteAn = Exp(-t) * (-t)**(- 2 - self%b1) * &
      gammaRcharm2(self%mL/lambdaQCD * self%Andim%Gfun(self%run,t) )

    end function InteAn

  end function MSRMass

end module VFNSMSRClass
