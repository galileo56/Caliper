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
    type (Alpha)                 :: AlphaOb
    real (dp)                    :: mH, mL, beta0, rat
    real (dp), dimension(0:4)    :: bHat

    contains

    procedure, pass(self), public :: MSRMass, setMass, RunArray, numFlav, alphaAll

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
     AnDim       = AlphaMass(2)%adim()        ; InitMSR%bHat       = AnDim%betaQCD('bHat')
     InitMSR%AnDim = AnDim; bHat = AnDim%betaQCD('beta'); InitMSR%beta0 = bHat(0)
     InitMSR%rat = InitMSR%mL/InitMSR%mH
     InitMSR%AlphaOb = AlphaMass(1)%AlphaAll()

   end function InitMSR

!ccccccccccccccc

  function RunArray(self,i) result(res)
    class (VFNSMSR) , intent(in) :: self
    integer         , intent(in) :: i
    type (Running)               :: res

    res = self%AlphaMass(i)

  end function RunArray

!ccccccccccccccc

  function AlphaAll(self) result(res)
    class (VFNSMSR) , intent(in) :: self
    type (Alpha)                 :: res

    res = self%AlphaOb

  end function AlphaAll

!ccccccccccccccc

  integer function numFlav(self)
    class (VFNSMSR) , intent(in) :: self

    numFlav = self%nf

  end function numFlav

!ccccccccccccccc

  recursive real (dp) function MSRMass(self, up, type, order, R, lambda, method) result(res)
    class (VFNSMSR)              , intent(in) :: self
    real (dp)                    , intent(in) :: R
    real (dp), optional          , intent(in) :: lambda
    character (len = *)          , intent(in) :: type, up
    character (len = *), optional, intent(in) :: method
    real (dp)                  , dimension(4) :: a
    real (dp)                  , dimension(2) :: bet2
    integer                      , intent(in) :: order
    integer                                   :: neval, ier, i
    real (dp)                                 :: corr, abserr, lambdaQCD, t0, &
    delta, mu, Rstep, delta1, delta2, alpha, alpha2, alpha3, matching, alphaM,&
    rat1, rat2, t1, b1, b2, mol

    res = 0; b1 = self%bHat(1); b2 = self%bHat(2)

    if ( up(:4) == 'down' ) then

      if ( type(:4) == 'MSRp' ) then
        alphaM = self%AlphaMass(2)%alphaQCD(self%mL)/Pi
        a = self%AlphaMass(2)%MSRMatching('charm'); i = min(order, 4)
        matching = dot_product( a(:i), PowList(alphaM, i) ) - 1
      else
        matching = - 1
      end if

      res = self%MSRMass('up', type, order, self%mL, lambda, method)  + &
      self%AlphaMass(1)%MSRMass(type, order, R, lambda, method) + &
      self%mL * matching
      return
    end if

    res = self%AlphaMass(2)%MSRMass(type, order, R, lambda, method)

    if ( self%mL < tiny(1._dp) ) return

    if ( type(:4) == 'MSRn' .and. order >= 3 ) then

      res = res + self%mH * deltaCharmNh(self%rat) * &
      ( self%AlphaOb%alphaQCD(self%nf + 1, self%mH)/Pi )**3

    end if

    if (self%run < 2) return

    if ( method(:7) == 'numeric' ) then

      call qags( InteNum, R/lambda, self%mH/lambda, prec, prec, corr, &
      abserr, neval, ier )

      res = res + corr

    else if ( method(:8) == 'analytic' ) then

      t0 = - 2 * pi/self%beta0/lambda;  t1 = t0/self%AlphaMass(2)%alphaQCD(self%mH)
      t0 = t0/self%AlphaMass(2)%alphaQCD(R);  bet2(2) = 2 * self%beta0
      lambdaQCD = self%AlphaMass(2)%lambdaQCD( self%run )
      mol = self%mL/lambdaQCD; bet2(1) = bet2(2)**2; bet2(2) = Product(bet2)

      call qags( InteAn, t1, t0, prec, prec, corr, abserr, neval, ier )

      res = res + lambdaQCD * corr

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( ( self%mH - R )/delta )  )
      delta = (self%mH - R)/neval; corr = 0; Rstep = self%mH + delta

      do ier = 1, neval

        Rstep  = Rstep - delta; mu = (Rstep - delta/2)/lambda
        alpha  = self%AlphaMass(2)%alphaQCD(mu)/Pi;  alpha2 = alpha**2
        rat1 = self%mL/Rstep; rat2 = self%mL/(Rstep - delta )

        delta1 = alpha2 * deltaCharm2(rat1)
        delta2 = alpha2 * deltaCharm2(rat2)

        if (self%run >= 3) then

          alpha3 = alpha * alpha2

          if ( type(:4) == 'MSRp' ) then
            delta1 = delta1 + deltaCharm3(self%nf, 1, rat1) * alpha3
            delta2 = delta2 + deltaCharm3(self%nf, 1, rat2) * alpha3
          else if ( type(:4) == 'MSRn' ) then
            delta1 = delta1 + deltaCharm3(self%nf, 0, rat1) * alpha3
            delta2 = delta2 + deltaCharm3(self%nf, 0, rat2) * alpha3
          end if

        end if

        corr = corr + Rstep * delta1 - (Rstep - delta) * delta2

      end do

      res = res + corr

    end if

  contains

!ccccccccccccccc

    real (dp) function InteNum(mu)
      real (dp), intent(in) :: mu
      real (dp)             :: aPi, api2, api3

      aPi = self%AlphaMass(2)%alphaQCD(mu)/4/Pi; aPi2 = aPi**2

      InteNum = gammaRcharm2(self%mL/mu) * aPi2

      if (self%run >= 3) then

        aPi3 = aPi * aPi2

        if ( type(:4) == 'MSRp' ) then
          InteNum = InteNum + gammaRcharm3(self%nf, 1, self%mL/mu) * aPi3
        else if ( type(:4) == 'MSRn' ) then
          InteNum = InteNum + gammaRcharm3(self%nf, 0, self%mL/mu) * aPi3
        end if

      end if

    end function InteNum

!ccccccccccccccc

    real (dp) function InteAn(t)
      real (dp), intent(in) :: t
      real (dp)             :: gammaR, gammaR2

      gammaR = gammaRcharm2(mol * self%Andim%Gfun(self%run,t) )/bet2(1)

      InteAn = (-t)**(- 2 - b1) * gammaR

      if (self%run >= 3) then

        if ( type(:4) == 'MSRp' ) then
          gammaR2 = gammaRcharm3(self%nf, 1, mol * self%Andim%Gfun(self%run,t) )
        else if ( type(:4) == 'MSRn' ) then
          gammaR2 = gammaRcharm3(self%nf, 0, mol * self%Andim%Gfun(self%run,t) )
        end if

        gammaR = gammaR2/bet2(2) - (b1 + b2) * gammaR
        InteAn = InteAn + (-t)**(- 3 - b1) * gammaR

      end if

      InteAn = Exp(-t) * InteAn

    end function InteAn

  end function MSRMass

!ccccccccccccccc

  subroutine setMass(self, m, mu)
    class (VFNSMSR), intent(inout) :: self
    real (dp)      , intent(in)    :: m, mu

    self%mH = m

    if (self%nl == 5) then
      call self%alphaMass(1)%SetMTop(m, mu)
      call self%alphaMass(2)%SetMTop(m, mu)
    else if (self%nl == 4) then
      call self%alphaMass(1)%SetMBottom(m, mu)
      call self%alphaMass(2)%SetMBottom(m, mu)
    end if

  end subroutine setMass

!ccccccccccccccc

end module VFNSMSRClass
