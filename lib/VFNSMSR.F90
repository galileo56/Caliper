module VFNSMSRClass
  use AnomDimClass;  use AlphaClass;  use RunningClass
  use Constants, only: dp, Pi, ExpEuler, d1mach, prec, pi2
  use QuadPack , only: qags         ; implicit none   ;  private

  public :: VFNSMSR

!ccccccccccccccc

  type, public                   :: VFNSMSR
    private

    integer                      :: nf, nl, run
    type (Running), dimension(3) :: AlphaMass
    type (AnomDim), dimension(3) :: AnDim
    type (Alpha)                 :: AlphaOb
    real (dp)                    :: mH, mL, beta0, rat
    real (dp), dimension(0:4)    :: bHat

    contains

    procedure, pass(self), public :: MSRMass, setMass, RunArray, numFlav, &
    MSRDelta, alphaAll, deltaM, mass, DeltaM2, setCharm, SetAlpha, OptimalR

    procedure, pass(self), private :: LowMass, LowRevol

  end type VFNSMSR

!ccccccccccccccc

  type, extends(VFNSMSR), public :: topMSR
    private

    integer                      :: nT
    real (dp), dimension(0:4)    :: bHTop
    real (dp)                    :: mC, mB, mT, betaTop, ratBottom, ratCharm

    contains

    procedure, pass(self), private :: MSRTop

  end type topMSR

!ccccccccccccccc

  interface VFNSMSR
    module procedure  InitMSR
  end interface VFNSMSR

!ccccccccccccccc

  interface topMSR
    module procedure  InitMSRtop
  end interface topMSR

contains

!ccccccccccccccc

  type (VFNSMSR) function InitMSR(AlphaMass)
    type (Running), dimension(2), intent(in) :: AlphaMass
    type (AnomDim), dimension(2)             :: AnDim
    real (dp), dimension(0:4)                :: bHat

    InitMSR%AlphaMass(:2) = AlphaMass; InitMSR%nf = AlphaMass(2)%numFlav()
    InitMSR%run = AlphaMass(2)%orders('runMass')
    InitMSR%nl  = AlphaMass(1)%numFlav(); InitMSR%mL = AlphaMass(1)%scales('mH')
    InitMSR%mH  = AlphaMass(2)%scales('mH')
    AnDim       = [ AlphaMass(1)%adim(), AlphaMass(2)%adim() ]
    InitMSR%bHat = AnDim(2)%betaQCD('bHat'); InitMSR%AnDim(:2) = AnDim
    bHat = AnDim(2)%betaQCD('beta'); InitMSR%beta0 = bHat(0)
    InitMSR%rat = InitMSR%mL/InitMSR%mH
    InitMSR%AlphaOb = AlphaMass(1)%AlphaAll()

   end function InitMSR

!ccccccccccccccc

   type (topMSR) function InitMSRtop(AlphaMass)
    type (Running), dimension(3), intent(in) :: AlphaMass
    type (AnomDim), dimension(3)             :: AnDim
    real (dp), dimension(0:4)                :: bHat, bHTop

    InitMSRtop%AlphaMass = AlphaMass
    InitMSRtop%run = AlphaMass(2)%orders('runMass')
    InitMSRtop%nf  = AlphaMass(2)%numFlav()
    InitMSRtop%nl  = AlphaMass(1)%numFlav()
    InitMSRtop%nT  = AlphaMass(3)%numFlav()
    InitMSRtop%mT  = AlphaMass(3)%scales('mH')
    InitMSRtop%mH  = AlphaMass(2)%scales('mH') ; InitMSRtop%mB = InitMSRtop%mH
    InitMSRtop%mL  = AlphaMass(1)%scales('mH') ; InitMSRtop%mC = InitMSRtop%mL
    AnDim = [ AlphaMass(1)%adim(), AlphaMass(2)%adim(), AlphaMass(3)%adim() ]
    InitMSRtop%bHat = AnDim(2)%betaQCD('bHat'); InitMSRtop%AnDim = AnDim
    bHat = AnDim(2)%betaQCD('beta'); InitMSRtop%beta0 = bHat(0)
    InitMSRtop%rat = InitMSRtop%mL/InitMSRtop%mH
    InitMSRtop%ratCharm = InitMSRtop%mC/InitMSRtop%mT
    InitMSRtop%ratBottom = InitMSRtop%mB/InitMSRtop%mT
    InitMSRtop%AlphaOb = AlphaMass(1)%AlphaAll()
    bHTop = AnDim(3)%betaQCD('beta'); InitMSRtop%betaTop = bHTop(0)
    InitMSRtop%bHTop = AnDim(3)%betaQCD('bHat')

  end function InitMSRtop

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

  real (dp) function DeltaM(self, up, R)
    class (VFNSMSR)    , intent(in) :: self
    character (len = *), intent(in) :: up
    real (dp)          , intent(in) :: R

    DeltaM = 0

    select type (self)
    type is (VFNSMSR)
      if ( self%mL <= tiny(1._dp) ) return
      if ( up(:2) == 'up' ) DeltaM = deltaCharm2(self%mL/R)
    type is (topMSR)
      if ( self%mB <= tiny(1._dp) .and. self%mC <= tiny(1._dp) ) return
      if ( up(:2) == 'up'   ) DeltaM = deltaCharm2(self%mB/R) + deltaCharm2(self%mC/R)
      if ( up(:4) == 'down' ) DeltaM = deltaCharm2(self%mC/R)
    end select

  end function DeltaM

!ccccccccccccccc

  real (dp) function DeltaM2(self, type, up, R)
    class (VFNSMSR)    , intent(in) :: self
    character (len = *), intent(in) :: up, type
    real (dp)          , intent(in) :: R

    DeltaM2 = 0

    select type (self)
    type is (VFNSMSR)

      if ( self%mL <= tiny(1._dp) ) return

      if ( up(:2) == 'up' ) then

        if ( type(:4) == 'MSRn' ) then
          DeltaM2 = deltaCharm3(self%nf, 0, self%mL/R)
        else if ( type(:4) == 'MSRp' .or. type(:5) == 'MSbar' ) then
          DeltaM2 = deltaCharm3(self%nf, 1, self%mL/R)
        end if
      end if

    type is (topMSR)

      if ( up(:2) == 'up' ) then

        if ( type(:4) == 'MSRn' ) then

          DeltaM2 = deltaCharm3(self%nT, 0, self%mB/R) + &
          deltaCharm3(self%nf, 0, self%mC/R)

        else if ( type(:4) == 'MSRp' .or. type(:5) == 'MSbar' ) then

          DeltaM2 = deltaCharm3(self%nT, 1, self%mB/R) + &
          deltaCharm3(self%nT, 1, self%mC/R)

        end if

        DeltaM2 = DeltaM2 + deltaBottomCharm(self%mB/R, self%mC/R)

      else if ( up(:4) == 'down' ) then

        if ( type(:4) == 'MSRn' ) then
          DeltaM2 = deltaCharm3(self%nf, 0, self%mC/R)
        else if ( type(:4) == 'MSRp' .or. type(:5) == 'MSbar' ) then
          DeltaM2 = deltaCharm3(self%nf, 1, self%mC/R)
        end if

      end if

    end select

  end function DeltaM2

!ccccccccccccccc

  real (dp) function mass(self)
    class (VFNSMSR), intent(in) :: self

    select type (self)
    type is (VFNSMSR); mass = self%mH
    type is (topMSR) ; mass = self%mT
    end select

  end function mass

!ccccccccccccccc

  function MSRDelta(self, up, type, R) result(res)
    class (VFNSMSR)    , intent(in) :: self
    character (len = *), intent(in) :: type, up
    real (dp)          , intent(in) :: R
    real (dp)        , dimension(4) :: res
    real (dp)                       :: rat, ratC

    select type (self)
    type is (VFNSMSR)

      if ( up(:2) == 'up' ) then

        res = self%AnDim(2)%MSRDelta(type)        ; rat = self%mL/R
        res(2) = res(2) + deltaCharm2(rat)

        if ( type(:4) == 'MSRn' ) then
          res(3) = res(3) + deltaCharm3(self%nf, 0, rat)
        else if ( type(:4) == 'MSRp' ) then
          res(3) = res(3) + deltaCharm3(self%nf, 1, rat)
        end if

      else
        res = self%AnDim(1)%MSRDelta(type)
      end if

    type is (topMSR)

      if ( up(:2) == 'up' ) then
        res = self%AnDim(3)%MSRDelta(type) ; rat = self%mB/R ; ratC = self%mC/R
        res(2) = res(2) + deltaCharm2(rat) + deltaCharm2(ratC)

        if ( type(:4) == 'MSRn' ) then
          res(3) = res(3) + deltaCharm3(self%nT, 0, rat) + deltaCharm3(self%nT, 0, ratC)
        else if ( type(:4) == 'MSRp' ) then
          res(3) = res(3) + deltaCharm3(self%nT, 1, rat) + deltaCharm3(self%nT, 1, ratC)
        end if

        res(3) = res(3) + deltaBottomCharm(rat, ratC)

      else if ( up(:4) == 'down' ) then

        res = self%AnDim(2)%MSRDelta(type)        ; rat = self%mC/R
        res(2) = res(2) + deltaCharm2(rat)

        if ( type(:4) == 'MSRn' ) then
          res(3) = res(3) + deltaCharm3(self%nf, 0, rat)
        else if ( type(:4) == 'MSRp' ) then
          res(3) = res(3) + deltaCharm3(self%nf, 1, rat)
        end if

      end if
    end select

  end function MSRDelta

!ccccccccccccccc

  real (dp) function MSRMass(self, up, type, order, R, lambda, method)
    class (VFNSMSR)    , intent(in) :: self
    real (dp)          , intent(in) :: R
    real (dp), optional, intent(in) :: lambda
    character (len = *), intent(in) :: type, up, method
    integer            , intent(in) :: order

    select type (self)
    type is (VFNSMSR)
      MSRMass = self%LowMass(up, type, order, R, lambda, method)
    type is (topMSR)
      MSRMass = self%MSRTop(up, type, order, R, lambda, method)
    end select

  end function MSRMass

!ccccccccccccccc

  recursive real (dp) function MSRTop(self, up, type, order, R, lambda, method) result(res)
    class (topMSR)       , intent(in) :: self
    real (dp)            , intent(in) :: R
    real (dp), optional  , intent(in) :: lambda
    character (len = *)  , intent(in) :: type, up, method
    integer              , intent(in) :: order
    real (dp)          , dimension(4) :: a
    real (dp)          , dimension(2) :: bet2
    character (len = 5)               :: upLow
    integer                           :: neval, ier, i
    real (dp)                         :: corr, abserr, lambdaQCD, t0, ratC, mu,&
    delta, Rstep, delta1, delta2, alpha, alpha2, alpha3, matching, alphaM, rat,&
    rat1, rat2, t1, b1, b2, mol, molC, rat2C, ratB, rat1C, alphaQ

    res = 0; b1 = self%bHTop(1); b2 = self%bHTop(2)

    if ( up(:2) /= 'up' ) then

      if ( up(:4) == 'down' ) then
        upLow = 'up'
      else
        upLow = 'down'
      end if

      if ( type(:4) == 'MSRp' ) then

        alphaM = self%AlphaMass(3)%alphaQCD(self%mB)/Pi
        a = self%AlphaMass(3)%MSRMatching('charm'); i = min(order, 4)
        matching = dot_product( a(:i), PowList(alphaM, i) ) - 1
        if (i > 2) matching = matching + deltaCharmNh(self%ratCharm) * alphaM**3

      else if ( type(:4) == 'MSRh' ) then
        if ( type(5:5) == 'V' ) then
          res = self%MSRTop('up', 'MSRnV', order, self%mB, lambda, method) + &
          self%LowMass(upLow, 'MSRpV', order, R, lambda, method) - self%mB
        else
          res = self%MSRTop('up', 'MSRn', order, self%mB, lambda, method) + &
          self%LowMass(upLow, 'MSRp', order, R, lambda, method) - self%mB
        end if

      else
        matching = - 1
      end if

      if ( type(5:5) == 'V' .and. order > 1 ) then

        rat = self%mB/self%mT;  ratC = self%mC/self%mT;  ratB = self%mC/self%mB
        mu = 2 * (  (2 * lambda - 1) * self%mT/2 + (2 - lambda) * self%mB  )/3
        alphaQ = self%AlphaMass(3)%alphaQCD(mu)/Pi

        res = res + (  self%mT * ( deltaCharm2(rat) + deltaCharm2(ratC) ) &
        - self%mB * ( deltaCharm2(1._dp) + deltaCharm2(ratB) )  ) * alphaQ**2

        if (order > 2) then
          if ( type(:4) == 'MSRp' .or. type(:4) == 'MSRh' ) then

            res = res + (   self%mT * (  deltaCharm3(self%nT, 1, rat) + &
            deltaCharm3(self%nT, 1, ratC) + self%betaTop * log(mu/self%mT) &
            * ( deltaCharm2(rat) + deltaCharm2(ratC) )  ) &
            - self%mB * (  deltaCharm3(self%nT, 1, 1._dp) + &
            deltaCharm3(self%nT, 1, ratB) + self%betaTop * log(mu/self%mB) &
            * ( deltaCharm2(1._dp) + deltaCharm2(ratB) )  )   ) * alphaQ**3

          else if ( type(:4) == 'MSRn' ) then

            res = res + (   self%mT * ( deltaCharm3(self%nT, 0, rat) + &
            deltaCharm3(self%nT, 0, ratC) + self%betaTop * log(mu/self%mT) &
            * ( deltaCharm2(rat) + deltaCharm2(ratC) )  ) &
            - self%mB * (  deltaCharm3(self%nT, 0, 1._dp) + &
            deltaCharm3(self%nT, 0, ratB) + self%betaTop * log(mu/self%mB) &
            * ( deltaCharm2(1._dp) + deltaCharm2(ratB) )  )   ) * alphaQ**3

          end if

          res = res + ( self%mT * deltaBottomCharm(rat, ratC) - &
          self%mB * deltaBottomCharm(1._dp, ratB) ) * alphaQ**3

        end if
      end if

      if ( type(:4) == 'MSRh' ) return

      res = res + self%MSRTop('up', type, order, self%mB, lambda, method) + &
      self%LowMass(upLow, type, order, R, lambda, method) + &
      self%mB * matching

      return

    end if

    res = self%AlphaMass(3)%MSRMass(type, order, R, lambda, method)

    if ( self%mB < tiny(1._dp) .and. self%mC < tiny(1._dp) ) return

    if ( type(:4) == 'MSRn' .and. order > 2 ) then

      res = res + self%mT * ( deltaCharmNh(self%ratBottom) + deltaCharmNh(self%ratCharm) ) * &
      ( self%AlphaOb%alphaQCD(self%nT + 1, self%mT)/Pi )**3

    end if

    if ( self%run < 2 .or. type(5:5) == 'V' ) return

    if ( method(:7) == 'numeric' ) then

      call qags( InteNum, R/lambda, self%mT/lambda, prec, prec, corr, &
      abserr, neval, ier )

      res = res + lambda * corr

    else if ( method(:8) == 'analytic' ) then

      t0 = - 2 * pi/self%betaTop;  t1 = t0/self%AlphaMass(3)%alphaQCD(self%mT/lambda)
      t0 = t0/self%AlphaMass(3)%alphaQCD(R/lambda);  bet2(2) = 2 * self%betaTop
      lambdaQCD = lambda * self%AlphaMass(3)%lambdaQCD(self%run)
      mol = self%mB/lambdaQCD; bet2(1) = bet2(2)**2; bet2(2) = Product(bet2)
      molC = self%mC/lambdaQCD

      call qags( InteAn, t1, t0, prec, prec, corr, abserr, neval, ier )

      res = res + lambda * lambdaQCD * corr

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( ( self%mT - R )/delta )  )
      delta = (self%mT - R)/neval; corr = 0; Rstep = self%mT + delta

      do ier = 1, neval

        Rstep  = Rstep - delta; mu = (Rstep - delta/2)/lambda
        alpha  = self%AlphaMass(3)%alphaQCD(mu)/Pi;  alpha2 = alpha**2
        rat1 = self%mB/Rstep; rat2 = self%mB/(Rstep - delta )
        rat1C = self%mC/Rstep; rat2C = self%mC/(Rstep - delta )

        delta1 = alpha2 * ( deltaCharm2(rat1) + deltaCharm2(rat1C) )
        delta2 = alpha2 * ( deltaCharm2(rat2) + deltaCharm2(rat2C) )

        if (self%run >= 3) then

          alpha3 = alpha * alpha2

          if ( type(:4) == 'MSRp' ) then

            delta1 = delta1 + alpha3 * (  deltaCharm3(self%nT, 1, rat1) + &
            deltaCharm3(self%nT, 1, rat1C) )

            delta2 = delta2 + alpha3 * ( deltaCharm3(self%nT, 1, rat2) + &
            deltaCharm3(self%nT, 1, rat2C) )

          else if ( type(:4) == 'MSRn' ) then

            delta1 = delta1 + alpha3 * ( deltaCharm3(self%nT, 0, rat1) + &
            deltaCharm3(self%nT, 0, rat1C) )

            delta2 = delta2 + alpha3 * ( deltaCharm3(self%nT, 0, rat2) + &
            deltaCharm3(self%nT, 0, rat2C) )

          end if

          delta1 = delta1 + alpha3 * deltaBottomCharm(rat1, rat1C)
          delta2 = delta2 + alpha3 * deltaBottomCharm(rat2, rat2C)

        end if

        corr = corr + Rstep * delta1 - (Rstep - delta) * delta2

      end do

      res = res + corr

    end if

  contains

!ccccccccccccccc

    real (dp) function InteNum(R)
      real (dp), intent(in) :: R
      real (dp)             :: aPi, api2, api3

      aPi = self%AlphaMass(3)%alphaQCD(R)/4/Pi; aPi2 = aPi**2

      InteNum = ( gammaRcharm2(self%mB/R/lambda) + gammaRcharm2(self%mC/R/lambda) ) * aPi2

      if (self%run >= 3) then

        if ( abs(lambda - 1) >= tiny(1._dp) ) &
        InteNum = InteNum - 4 * log(lambda) * InteNum * aPi

        aPi3 = aPi * aPi2

        if ( type(:4) == 'MSRp' ) then
          InteNum = InteNum + ( gammaRcharm3(self%nT, 1, self%mB/R/lambda) + &
          gammaRcharm3(self%nT, 1, self%mC/R/lambda) ) * aPi3
        else if ( type(:4) == 'MSRn' ) then
          InteNum = InteNum + ( gammaRcharm3(self%nT, 0, self%mB/R/lambda) + &
          gammaRcharm3(self%nT, 0, self%mC/R/lambda) ) * aPi3
        end if

        InteNum = InteNum + aPi3 * GammaRBottomCharm(self%mB/R/lambda, self%mC/R/lambda)

      end if

    end function InteNum

!ccccccccccccccc

    real (dp) function InteAn(t)
      real (dp), intent(in) :: t
      real (dp)             :: gammaR, gammaR2

      gammaR = gammaRcharm2(mol * self%Andim(3)%Gfun(self%run,t) ) + &
      gammaRcharm2(molC * self%Andim(3)%Gfun(self%run,t) )

      InteAn = (-t)**(- 2 - b1) * gammaR/bet2(1)

      if (self%run >= 3) then

        if ( type(:4) == 'MSRp' ) then
          gammaR2 = gammaRcharm3(self%nT, 1, mol * self%Andim(3)%Gfun(self%run,t) ) &
          + gammaRcharm3(self%nT, 1, molC * self%Andim(3)%Gfun(self%run,t) )
        else if ( type(:4) == 'MSRn' ) then
          gammaR2 = gammaRcharm3(self%nT, 0, mol * self%Andim(3)%Gfun(self%run,t) ) &
          + gammaRcharm3(self%nT, 0, molC * self%Andim(3)%Gfun(self%run,t) )
        end if

        gammaR2 = gammaR2 + GammaRBottomCharm(mol * self%Andim(3)%Gfun(self%run,t), &
        molC * self%Andim(3)%Gfun(self%run,t) )

        if ( abs(lambda - 1) >= tiny(1._dp) ) &
        gammaR2 = gammaR2 - 4 * log(lambda) * gammaR

        gammaR = gammaR2/bet2(2) - (b1 + b2) * gammaR/bet2(1)
        InteAn = InteAn + (-t)**(- 3 - b1) * gammaR

      end if

      InteAn = Exp(-t) * InteAn

    end function InteAn

  end function MSRTop

!ccccccccccccccc

  recursive real (dp) function LowMass(self, up, type, order, R, lambda, method) result(res)
    class (VFNSMSR)     , intent(in) :: self
    real (dp)           , intent(in) :: R
    real (dp), optional , intent(in) :: lambda
    character (len = *) , intent(in) :: type, up, method
    integer             , intent(in) :: order
    real (dp)         , dimension(4) :: a
    integer                          :: i
    real (dp)                        :: matching, alphaM, rat, alphaQ, mu

    res = 0

    if ( up(:4) == 'down' ) then

      if ( type(:4) == 'MSRp' ) then
        alphaM = self%AlphaMass(2)%alphaQCD(self%mL)/Pi
        a = self%AlphaMass(2)%MSRMatching('charm'); i = min(order, 4)
        matching = dot_product( a(:i), PowList(alphaM, i) ) - 1
      else if ( type(:4) == 'MSRh' ) then

        if ( type(5:5) == 'V' ) then
          res = self%LowMass('up', 'MSRnV', order, self%mL, lambda, method)
        else
          res = self%LowMass('up', 'MSRn', order, self%mL, lambda, method)
        end if

        res = res + self%AlphaMass(1)%MSREvol('MSRp', R, lambda, method)

      else
        matching = - 1
      end if

      if ( type(5:5) == 'V' .and. order > 1 ) then

        rat = self%mL/self%mH
        mu = 2 * (  (2 * lambda - 1)/2 * self%mH + (2 - lambda) * self%mL  )/3
        alphaQ = self%AlphaMass(2)%alphaQCD(mu)/Pi

        res = res + ( self%mH * deltaCharm2(rat) - &
        self%mL * deltaCharm2(1._dp) ) * alphaQ**2

        if (order > 2) then
          if ( type(:4) == 'MSRp' .or. type(:4) == 'MSRh' ) then

            res = res + (  self%mH * ( deltaCharm3(self%nf, 1, rat) + &
            self%beta0 * deltaCharm2(rat) * log(mu/self%mH) )  - &
            self%mL * ( deltaCharm3(self%nf, 1, 1._dp) + &
            self%beta0 * deltaCharm2(1._dp) * log(mu/self%mL) )  ) * alphaQ**3

          else if ( type(:4) == 'MSRn' ) then

            res = res + (  self%mH * ( deltaCharm3(self%nf, 0, rat) + &
            self%beta0 * deltaCharm2(rat) * log(mu/self%mH) )  - &
            self%mL * ( deltaCharm3(self%nf, 0, 1._dp) + &
            self%beta0 * deltaCharm2(1._dp) * log(mu/self%mL) )  ) * alphaQ**3

          end if
        end if

      end if

      if ( type(:4) == 'MSRh' ) return

      res = res + self%LowMass('up', type, order, self%mL, lambda, method) + &
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

    if (self%run < 2 .or. type(5:5) == 'V' ) return

    res = res + self%LowRevol(type, self%mH, R, lambda, method)

  end function LowMass

!ccccccccccccccc

  real (dp) function LowRevol(self, type, R0, R, lambda, method)
    class (VFNSMSR)     , intent(in) :: self
    real (dp)           , intent(in) :: R0, R
    real (dp), optional , intent(in) :: lambda
    character (len = *) , intent(in) :: type, method
    real (dp)         , dimension(2) :: bet2
    integer                          :: neval, ier
    real (dp)                        :: abserr, t0, t1, lambdaQCD, b1, delta, &
    mu, Rstep, delta1, delta2, alpha, alpha2, alpha3, rat1, rat2, mol, b2

    LowRevol = 0; b1 = self%bHat(1); b2 = self%bHat(2)

    if ( method(:7) == 'numeric' ) then

      call qags( InteNum, R/lambda, R0/lambda, prec, prec, LowRevol, &
      abserr, neval, ier )

      LowRevol = lambda * LowRevol

    else if ( method(:8) == 'analytic' ) then

      t0 = - 2 * pi/self%beta0;  t1 = t0/self%AlphaMass(2)%alphaQCD(R0/lambda)
      t0 = t0/self%AlphaMass(2)%alphaQCD(R/lambda);  bet2(2) = 2 * self%beta0
      lambdaQCD = lambda * self%AlphaMass(2)%lambdaQCD(self%run)
      mol = self%mL/lambdaQCD; bet2(1) = bet2(2)**2; bet2(2) = Product(bet2)

      call qags( InteAn, t1, t0, prec, prec, LowRevol, abserr, neval, ier )

      LowRevol = lambdaQCD * LowRevol

    else if ( method(:4) == 'diff' ) then

      delta = 0.1_dp;   neval = abs(  Nint( (R0 - R)/delta )  )
      delta = (R0 - R)/neval; LowRevol = 0; Rstep = R0 + delta

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

        LowRevol = LowRevol + Rstep * delta1 - (Rstep - delta) * delta2

      end do
    end if

  contains

!ccccccccccccccc

    real (dp) function InteNum(R)
      real (dp), intent(in) :: R
      real (dp)             :: aPi, api2, api3

      aPi = self%AlphaMass(2)%alphaQCD(R)/4/Pi; aPi2 = aPi**2

      InteNum = gammaRcharm2(self%mL/R/lambda) * aPi2

      if (self%run >= 3) then

        if ( abs(lambda - 1) >= tiny(1._dp) ) &
        InteNum = InteNum - 4 * log(lambda) * InteNum * aPi

        aPi3 = aPi * aPi2

        if ( type(:4) == 'MSRp' ) then
          InteNum = InteNum + gammaRcharm3(self%nf, 1, self%mL/R/lambda) * aPi3
        else if ( type(:4) == 'MSRn' ) then
          InteNum = InteNum + gammaRcharm3(self%nf, 0, self%mL/R/lambda) * aPi3
        end if

      end if

    end function InteNum

!ccccccccccccccc

    real (dp) function InteAn(t)
      real (dp), intent(in) :: t
      real (dp)             :: gammaR, gammaR2

      gammaR = gammaRcharm2(mol * self%Andim(2)%Gfun(self%run,t) )

      InteAn = (-t)**(- 2 - b1) * gammaR/bet2(1)

      if (self%run >= 3) then

        if ( type(:4) == 'MSRp' ) then
          gammaR2 = gammaRcharm3(self%nf, 1, mol * self%Andim(2)%Gfun(self%run,t) )
        else if ( type(:4) == 'MSRn' ) then
          gammaR2 = gammaRcharm3(self%nf, 0, mol * self%Andim(2)%Gfun(self%run,t) )
        end if

        if ( abs(lambda - 1) >= tiny(1._dp) ) &
        gammaR2 = gammaR2 - 4 * log(lambda) * gammaR

        gammaR = gammaR2/bet2(2) - (b1 + b2) * gammaR/bet2(1)
        InteAn = InteAn + (-t)**(- 3 - b1) * gammaR

      end if

      InteAn = Exp(-t) * InteAn

    end function InteAn

  end function LowRevol

!ccccccccccccccc

  real (dp) function OptimalR(self, up, type, n, order, lambda, method)
    class (VFNSMSR)                    , intent(in) :: self
    character (len = *)                , intent(in) :: type, up
    real (dp)          , optional      , intent(in) :: lambda
    character (len = *), optional      , intent(in) :: method
    integer                            , intent(in) :: order
    real (dp)                          , intent(in) :: n
    integer                                         :: IFLAG
    real (dp)                                       :: a, b, c, alpha

    a = 0.25_dp; b = self%mH

    call DFZERO(root, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG); OptimalR = a

  contains

    real (dp) function root(x)
      real (dp), intent(in) :: x

      if ( up(:2) == 'up' ) then
        alpha = self%AlphaMass(2)%alphaQCD(x)
      else if ( up(:4) == 'down' ) then
        alpha = self%AlphaMass(1)%alphaQCD(x)
      end if

      root = n * x - 4 * self%MSRMass(up, type, order, x, lambda, method) * alpha/3

    end function root

  end function OptimalR

!ccccccccccccccc

  subroutine setMass(self, m, mu)
    class (VFNSMSR), intent(inout) :: self
    real (dp)      , intent(in)    :: m, mu

    self%mH = m ; self%rat = self%mL/m

    if (self%nf == 5) then

      call self%alphaMass(1)%SetMTop(m, mu)
      call self%alphaMass(2)%SetMTop(m, mu)

      select type (self)
      type is (topMSR)
        call self%alphaMass(3)%SetMTop(m, mu)
      end select

      call self%alphaOb%SetMTop(m, mu)

    else if (self%nf == 4) then
      call self%alphaMass(1)%SetMBottom(m, mu)
      call self%alphaMass(2)%SetMBottom(m, mu)
      call self%alphaOb%SetMBottom(m, mu)

    else if (self%nf == 3) then
      call self%alphaMass(1)%SetMCharm(m, mu)
      call self%alphaMass(2)%SetMCharm(m, mu)
      call self%alphaOb%SetMBottom(m, mu)
    end if

  end subroutine setMass

!ccccccccccccccc

  subroutine SetAlpha(self, alpha)
    class (VFNSMSR), intent(inout) :: self
    real (dp)      , intent(in   ) :: alpha

    call self%AlphaOb%SetAlpha(alpha, 0._dp)
    call self%alphaMass(1)%SetAlpha(alpha)
    call self%alphaMass(2)%SetAlpha(alpha)

  end subroutine SetAlpha

!ccccccccccccccc

  subroutine setCharm(self, m, mu)
    class (VFNSMSR), intent(inout) :: self
    real (dp)      , intent(in)    :: m, mu

    self%mL = m

    if (self%nf == 5) then

      call self%alphaMass(1)%SetMBottom(m, mu)
      call self%alphaMass(2)%SetMBottom(m, mu)

      select type (self)
      type is (topMSR)
        call self%alphaMass(3)%SetMBottom(m, mu)
        ! self%ratCharm = self%mC/self%mT
      end select

      call self%alphaOb%SetMBottom(m, mu)

    else if (self%nf == 4) then
      call self%alphaMass(1)%SetMCharm(m, mu)
      call self%alphaMass(2)%SetMCharm(m, mu)
      call self%alphaOb%SetMCharm(m, mu)
    end if

    self%rat = self%mL/self%mH

  end subroutine setCharm

!ccccccccccccccc

end module VFNSMSRClass
