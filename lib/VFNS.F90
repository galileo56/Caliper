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
    integer                      , intent(in) :: order
    integer                                   :: neval, ier
    ! real (dp), dimension(self%run,self%run)   :: gammaLog
    ! real (dp), dimension(4)                   :: gammaR
    ! real (dp), dimension(self%run)            :: alphaList, lglist
    real (dp)                                 :: corr, abserr, lambdaQCD, t0, t1, &
    delta, mu, Rstep, delta1, delta2

    res = 0

    if (R < self%mL) then
      res = self%MSRMass(type, order, self%mL, lambda, method) + &
      self%AlphaMass(1)%MSREvol(type, order, R, lambda, method)
      return
    end if

    ! if ( type(:4) == 'diff' )
    res = self%AlphaMass(2)%MSREvol(type, order, R, lambda, method)

    if ( type(:7) == 'numeric' ) then

      call qags( InteNum, R/lambda, self%mH/lambda, prec, prec, corr, &
      abserr, neval, ier )

      res = res + corr

    else if ( type == 'analytic' ) then

      t0 = - 2 * pi/self%beta0/lambda;  t1 = t0/self%AlphaMass(2)%alphaQCD(self%mH)
      t0 = t0/self%AlphaMass(2)%alphaQCD(R)
      lambdaQCD = self%AlphaMass(2)%lambdaQCD( self%run )

      call qags( InteAn, t1, t0, prec, prec, corr, abserr, neval, ier )

      res = res + corr/4/self%beta0**2

    else if ( method(:4) == 'diff' ) then

      ! gammaR = self%AlphaMass(2)%gammaR(type)
      delta = 0.1_dp;   neval = abs(  Nint( ( self%mH - R )/delta )  )
      delta = (self%mH - R)/neval; corr = 0; Rstep = self%mH + delta
      ! gammaLog(1,:) = gammaR(:self%run)
      ! call self%andim%expandAlpha( gammaLog )

      do ier = 1, neval

        Rstep = Rstep - delta; mu = (Rstep - delta/2)/lambda
        ! alphaList = powList( self%AlphaMass(2)%alphaQCD(mu)/Pi, self%run )

        ! lgList = powList( log(mu/Rstep), self%run )
        delta1 = Rstep * ( self%AlphaMass(2)%alphaQCD(mu)/Pi )**2 * deltaCharm2(self%mL/Rstep)
        ! delta1 = Rstep * (  alphaList(2) * deltaCharm2(self%mL/Rstep) + &
        ! dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )  )

        ! lgList = powList( log(mu/(Rstep - delta) ), self%run )
        delta2 = (Rstep - delta) * ( self%AlphaMass(2)%alphaQCD(mu)/Pi )**2 * deltaCharm2( self%mL/(Rstep - delta )  )

        ! delta2 = (Rstep - delta) * (  alphaList(2) * deltaCharm2( self%mL/(Rstep - delta )  ) &
        ! + dot_product( DeltaComputer(gammaLog, lgList, 0), alphaList )   )

        corr = corr + delta1 - delta2

      end do

      res = res + corr

    end if

  contains

!ccccccccccccccc

    real (dp) function InteNum(mu)
      real (dp), intent(in) :: mu

      InteNum = gammaRcharm2(self%mL/mu) * self%AlphaMass(2)%alphaQCD(self%mH/mu)**2/16/Pi2

    end function InteNum

!ccccccccccccccc

    real (dp) function InteAn(t)
      real (dp), intent(in) :: t

      InteAn = Exp(-t) * (-t)**(- 2 - self%b1) * &
      gammaRcharm2(self%mL/lambdaQCD * self%Andim%Gfun(self%run,t) )

    end function InteAn

  end function MSRMass

end module VFNSMSRClass
