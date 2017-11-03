
module ProfilesClass
  use constants, only: dp; implicit none
  private

!ccccccccccccccc

  type, abstract, public :: Profiles
    private
    integer              :: ns
    real (dp)            :: muH, mu0, R0, t0, t1, t2, cnt, ts, tmid1, tmid2, &
    c11, c12, c21, c22, d11, d12, d21, d22, e11, e12, e21, e31, e22, c31, c32,&
    d31, d32, Q, n0, n1, const, tau2, eH, slope, e32, eS, eJ, b

    contains

    procedure, pass(self), public   :: Scales, HardScale, Energy, setQ

  end type Profiles

!ccccccccccccccc

  type, extends (Profiles), public :: ProfilesMassless
    private
    real (dp)                      :: eR, tR
  end type ProfilesMassless

!ccccccccccccccc

  type, extends (Profiles), public :: ProfilesPythia
    private
    character (len = 6)             :: EShape
    character (len = 1)             :: def
    real (dp)                       :: x32, muM, a31, b31, b32, x31, a32, mass, beta, &
                                       delta0, delta1, eJ0, ts0, deltaLambda
    real (dp)            , public   :: JetIntersect, SoftIntersect

    contains

    procedure, pass(self)           :: MassIntersection
    procedure, pass(self), public   :: MassScale, setMassMuM, Sectors

  end type ProfilesPythia

!ccccccccccccccc

  interface ProfilesPythia
    module procedure InProfPy
  end interface ProfilesPythia

!ccccccccccccccc

  interface ProfilesMassless
    module procedure InProf
  end interface ProfilesMassless

  contains

!ccccccccccccccc

   type (ProfilesMassless) function InProf(Q, mu0, R0, n0, n1, tau2, tR, ts, slope, cnt, &
                                           eH, eS, eJ, eR, ns)
    integer            , intent(in) :: ns
    real (dp)          , intent(in) :: Q, mu0, R0, n0, n1, tau2, tR, ts, slope, cnt, eH, &
                                       eS, eJ, eR

    InProf%n0 = n0    ;  InProf%n1 = n1 ; InProf%const = cnt   ;  InProf%eR = eR
    InProf%tau2 = tau2;  InProf%tR = tR ; InProf%slope = slope ;  InProf%eH = eH
    InProf%eS  = eS   ; InProf%eJ  = eJ; InProf%ns  = ns       ;  InProf%R0  = R0
    InProf%mu0 = mu0  ; InProf%tS  = tS

    call InProf%setQ(Q)

   end function InProf

!ccccccccccccccc

   type (ProfilesPythia) function InProfPy(Q, beta, mu0, deltaLambda, R0, n0, delta0, &
            n1, delta1, tau2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, def, EShape)
    character (len = 6), optional, intent(in) :: Eshape
    character (len = 1), optional, intent(in) :: def
    integer                      , intent(in) :: ns
    real (dp)          , optional, intent(in) :: mass, muM, eS
    real (dp)                    , intent(in) :: Q, mu0, R0, n0, n1, tau2, slope, delta1, &
                                                 cnt, eH, eJ, deltaLambda, delta0, beta, ts

    InProfPy%mass   = 0    ; if ( present(mass)   ) InProfPy%mass   = mass
    InProfPy%def    = 'N'  ; if ( present(def)    ) InProfPy%def    = def(:1)
    InProfPy%EShape = 'No' ; if ( present(EShape) ) InProfPy%EShape = EShape(:6)
    InProfPy%muM    = 1    ; if ( present(muM)    ) InProfPy%muM    = muM
    InProfPy%eS     = 0    ; if ( present(eS )    ) InProfPy%eS     = eS

    InProfPy%n0 = n0  ; InProfPy%n1  = n1 ; InProfPy%const = cnt; InProfPy%eH     = eH
    InProfPy%R0  = R0 ; InProfPy%ts0 = tS ; InProfPy%tau2 = tau2; InProfPy%delta0 = delta0
    InProfPy%eJ0 = eJ ; InProfPy%ns  = ns ; InProfPy%slope = slope; InProfPy%beta = beta
    InProfPy%mu0 = mu0; InProfPy%delta1  = delta1 ; InProfPy%deltaLambda = deltaLambda

    call InProfPy%setQ(Q)

   end function InProfPy

!ccccccccccccccc

  subroutine setQ(self,Q)
    class (Profiles), intent(inout) :: self
    real (dp)       , intent(in   ) :: Q
    real (dp)                       :: tdiff, tSdiff2, tS2diff, t1S, t1R, mu1, mu2, moQ, &
                                       taumin, tmax, b, t00, t0, t1, t2, t2s2, t01, t2s, &
                                       t012s
    self%Q = Q;  self%muH = self%eH * Q

    select type (self)
    type is (ProfilesMassless)

    self%t0 = self%n0/Q;  self%t1 = abs(self%n1)/Q ;  self%cnt  = self%const

    if ( abs(self%tau2) <= self%t1 ) then
      self%t1 = abs(self%tau2);  self%t2 = abs(self%tau2)
    else
      self%t2 = abs(self%tau2)
    end if

    self%b = self%slope * self%muH
    self%tmid1 = (self%t0 + self%t1)/2;  self%tmid2 = self%t2 + self%ts

    tdiff = self%t1 - self%t0    ; tSdiff2 = self%ts - self%t2
    tS2diff = self%tmid2 * tSdiff2 ; tSdiff2 = tSdiff2**2; self%tmid2 = self%tmid2/2

    if ( self%n1 < 0 .and. self%tau2 >= 0) then

      self%cnt = self%mu0 -  self%b * self%tmid1

      self%c11 = ( 2 * self%mu0 * tdiff - self%b * self%t0**2 )/tdiff/2

      self%d11 = - self%b * self%t0/tdiff; self%e11 = self%b/tdiff/2;  self%d12 = self%d11

      self%c12 = self%c11 ; self%e12 = self%e11

      self%c21 = - (  2 * self%cnt * (self%t2**2 + 2 * self%t2 * self%ts - self%ts**2) &
      + self%t2**2 * ( - 4 * self%muH + self%b * (self%t2 + 3 * self%ts) )  )/tSdiff2/2

      self%c22 = (  2 * self%muH * (self%t2**2 - 2 * self%t2 * self%ts - self%ts**2) + &
        self%ts**2 * ( 4 * self%cnt + self%b * (3 * self%t2 + self%ts) )  )/tSdiff2/2

      self%d21 = ( 4 * self%cnt * self%t2 - 4 * self%muH * self%t2 + &
      self%b * (2 * self%t2**2 + self%t2 * self%ts + self%ts**2) )/tSdiff2

      self%d22 = - self%ts * ( 4 * self%cnt - 4 * self%muH + &
      self%b * (3 * self%t2 + self%ts) )/tSdiff2

       self%e21 = - ( 4 * self%cnt - 4 * self%muH + self%b * &
                          (self%t2 + 3 * self%ts) )/tSdiff2/2

       self%e22 =   ( 4 * self%cnt - 4 * self%muH + self%b * &
                          (3 * self%t2 + self%ts) )/tSdiff2/2

       else if ( self%tau2 < 0 .and. self%n1 >= 0) then

       self%b   = (self%muH - self%cnt)/self%tmid2
       self%c21 = (self%const * self%ts**2 - self%muH * self%t2**2)/tS2diff
       self%d21 = self%ts * (self%muH - self%const)/tS2diff
       self%e21 = (self%const - self%muH)/tS2diff

       self%e22 = self%e21; self%c22 = self%c21; self%d22 = self%d21

       self%c11 = (  - 2 * self%mu0 * (self%t0**2 + 2 * self%t0 * self%t1 - self%t1**2) &
       + self%t0**2 * ( 4 * self%const + self%b * (self%t0 + 3 * self%t1) )  )/tdiff**2/2

       self%c12 = ( 2 * self%const * (self%t0**2 - 2 * self%t0 * self%t1 - self%t1**2) - &
       self%t1**2 * ( - 4 * self%mu0 + self%b * (3 * self%t0 + self%t1) )  )/tdiff**2/2

       self%d11 = - self%t0 * ( 4 * (self%const - self%mu0) + &
       self%b * (self%t0 + 3 * self%t1) )/tdiff**2

       self%d12 = ( 4 * (self%const - self%mu0) * self%t1 + self%b * &
       (self%t0**2 + self%t0 * self%t1 + 2 * self%t1**2) )/tdiff**2

       self%e11 = ( 4 * (self%const - self%mu0) + self%b * &
                           (self%t0 + 3 * self%t1) )/tdiff**2/2

       self%e12 = - ( 4 * (self%const - self%mu0) + self%b * &
                          (3 * self%t0 + self%t1) )/tdiff**2/2

       else if (self%n1 < 0 .and. self%tau2 < 0) then

       self%b   = (self%mu0 - self%muH)/(self%tmid1 - self%tmid2)

       self%cnt = ( self%muH * self%tmid1 - self%mu0 * self%tmid2 )/&
                        (self%tmid1 - self%tmid2)

       self%c11 = ( 2 * self%mu0 * tdiff + self%b * self%t0**2 )/tdiff/2
       self%d11 = - self%b * self%t0/tdiff;  self%e11 = self%b/tdiff/2

       self%c12 = self%c11; self%e12 = self%e11

       self%c21 = (self%cnt * self%ts**2 - self%muH * self%t2**2)/tS2diff
       self%d21 = 2 * self%ts * (self%muH - self%cnt)/tS2diff
       self%e21 = (self%cnt - self%muH)/tS2diff

       self%c22 = self%c21; self%e22 = self%e21; self%d22 = self%d21; self%d12 = self%d11

       else

       self%c11 = (  - 2 * self%mu0 * (self%t0**2 + 2 * self%t0 * self%t1 - self%t1**2) &
       + self%t0**2 * ( 4 * self%const + self%b * (self%t0 + 3 * self%t1) )  )/tdiff**2/2

       self%c12 = ( 2 * self%const * (self%t0**2 - 2 * self%t0 * self%t1 - self%t1**2) - &
       self%t1**2 * ( - 4 * self%mu0 + self%b * (3 * self%t0 + self%t1) )  )/tdiff**2/2

       self%d11 = - self%t0 * ( 4 * (self%const - self%mu0) + self%b * &
       (self%t0 + 3 * self%t1) )/tdiff**2

       self%d12 = ( 4 * (self%const - self%mu0) * self%t1 + self%b * (self%t0**2&
        + self%t0 * self%t1 + 2 * self%t1**2) )/tdiff**2

       self%e11 =   ( 4 * (self%const - self%mu0) + self%b * (self%t0 + &
                              3 * self%t1) )/tdiff**2/2
       self%e12 = - ( 4 * (self%const - self%mu0) + self%b * (3 * self%t0 + &
                             self%t1) )/tdiff**2/2

       self%c21 = - (  2 * self%const * (self%t2**2 + 2 * self%t2 * self%ts - self%ts**2) &
       + self%t2**2 * ( - 4 * self%muH + self%b * (self%t2 + 3 * self%ts) )  )/tSdiff2/2

       self%c22 =  (  2 * self%muH * (self%t2**2 - 2 * self%t2 * self%ts - self%ts**2) &
        + self%ts**2 * ( 4 * self%const + self%b * (3 * self%t2 + self%ts) )  )/tSdiff2/2

       self%d21 = ( 4 * (self%const - self%muH) * self%t2 + &
       self%b * (2 * self%t2**2 + self%t2 * self%ts + self%ts**2) )/tSdiff2

       self%d22 = - self%ts * ( 4 * (self%const - self%muH) + self%b * &
                               (3 * self%t2 + self%ts) )/tSdiff2

       self%e21 = - ( 4 * (self%const - self%muH) + self%b * &
                          (self%t2 + 3 * self%ts) )/tSdiff2/2

       self%e22 =   ( 4 * (self%const - self%muH) + self%b * &
                           (3 * self%t2 + self%ts) )/tSdiff2/2

       end if

     ! Setting constants for R scale

     t1S = self%ts - self%t1; t1R = self%t1 - self%tR

     mu1 = (self%cnt + self%b * self%t1) * ( 1 + self%eR * t1R**2 ) * &
           ( 1 + self%eS * tdiff**2 * t1S**2 )

     mu2 = self%b * ( 1 + self%eR * t1R**2 ) * ( 1 + self%eS * tdiff**2 * t1S**2 ) + &
           2 * self%eR * (self%cnt + self%b * self%t1) * t1R * &
           ( 1 + self%eS * tdiff**2 * t1S**2 ) + 2 * self%eS * tdiff * t1S * (self%cnt + &
           self%b * self%t1) * ( 1 + self%eR * t1R**2 ) * (self%t0 - 2 * self%t1 + self%ts)

     call JunctionSet(self%R0, mu1 - mu2 * self%t1, 0._dp, mu2, self%t0, self%t1, &
             self%c31, self%c32, self%d31, self%d32, self%e31, self%e32)

    type is (ProfilesPythia)

    moQ = self%mass/Q

    if ( self%def(:1) /= 'N' .and. self%Eshape(:2) /= 'No' .and. self%mass > 0 ) then
      if ( self%def(:1) == 'Q' ) then
        if ( self%EShape(:6) == 'thrust' ) then
          taumin = ( 1 - sqrt(1 - 4 * moQ**2) );  tmax = moQ
        else if ( self%EShape(:6) == 'Cparam' ) then
          taumin = 12 * moQ**2 * (1 - moQ**2)  ;  tmax = moQ * (2 + moQ)
        else if ( self%EShape(:3) == 'HJM'    ) then
          taumin = moQ**2                      ;  tmax = 1.0459267339672078_dp * moQ
        end if
      end if
    else
      tmax = 0; taumin = 0
    end if

    taumin = taumin + self%deltaLambda/Q;  t00 = self%n0/Q + self%delta0/sqrt(Q)
    t0 = t00 + taumin;  t1 = abs(self%n1)/Q**self%beta + self%delta1/sqrt(Q) + taumin
    self%ts = self%ts0 + tmax

    self%eJ = self%eJ0 * ( (self%ts0 - t00)/(self%ts - t0) )**2

    b = self%slope * self%muH; b = b * ( 1 + self%eS * self%ts0/(self%ts - taumin) )

    self%cnt = self%const - taumin * b;  if ( self%EShape(:6) == 'Cparam' ) b = b/6

    t2 = abs(self%tau2) + tmax; if ( abs(self%tau2) + moQ <= t1 ) t1 = abs(self%tau2) + tmax

    self%t1 = t1; self%t2 = t2; self%b = b ; self%t0  = t0
    self%tmid1 = (t0 + t1)/2;  self%tmid2 = (t2 + self%ts)/2

    t01 = t0 - t1; t2s = t2 - self%ts; t2s2 = t2**2 - self%ts**2

    if ( self%n1 < 0 .and. self%tau2 >= 0) then

      self%cnt = self%mu0 - b * (t0 + t1)/2
      self%c11 = ( 2 * self%mu0 * t01 - b * t0**2 )/t01/2
      self%d11 = b * t0/t01;  self%e11 = - b/t01/2

      self%c21 = - (  2 * self%cnt * (t2**2 + 2 * t2 * self%ts - self%ts**2) + &
      t2**2 * ( - 4 * self%muH + b * (t2 + 3 * self%ts) )  )/t2s**2/2

      self%c22 = (  2 * self%muH * (self%t2**2 - 2 * t2 * self%ts - self%ts**2) + &
      self%ts**2 * ( 4 * self%cnt + b * (3 * t2 + self%ts) )  )/t2s**2/2

      self%d21 = ( (4 * self%cnt - self%muH) * t2 + &
      b * (2 * t2**2 + t2 * self%ts + self%ts**2) )/t2s**2

      self%e21 = - ( (4 * self%cnt - self%muH) + b * (t2 + 3 * self%ts) )/t2s**2/2
      self%e22 =   ( (4 * self%cnt - self%muH) + b * (3 * t2 + self%ts) )/t2s**2/2

      self%d12 = self%d11; self%c12 = self%c11
      self%e12 = self%e11; self%d22 = - 2 * self%ts * self%e22

    else if ( self%tau2 < 0 .and. self%n1 >= 0 ) then

      self%b = 2 * (self%muH - self%cnt)/(t2 + self%ts)
      self%c21 = (self%muH * t2**2 - self%cnt * self%ts**2)/t2s2
      self%d21 = 2 * self%ts * (self%cnt - self%muH)/t2s2
      self%e21 = (self%muH - self%cnt)/t2s2

      self%c11 = (  t0**2 * ( 4 * self%cnt + b * (t0 + 3 * t1) ) - &
      2 * self%mu0 * (t0**2 + 2 * t0 * t1 - t1**2) )/t01**2/2

      self%c12 = ( 2 * self%cnt * (t0**2 - 2 * t0 * t1 - t1**2) - &
      t1**2 * ( b * (3 * t0 + t1) - 4 * self%mu0 )  )/t01**2/2

      self%d12 = ( 4 * (self%cnt - self%mu0) * t1 + b * (t0**2 + t0 * t1 + 2 * t1**2) )/t01**2
      self%e11 =   ( 4 * (self%cnt - self%mu0) + b * (t0 + 3 * t1) )/t01**2/2
      self%e12 = - ( 4 * (self%cnt - self%mu0) + b * (3 * t0 + t1) )/t01**2/2

      self%d22 = self%d21; self%c22 = self%c21
      self%d11 = - 2 * t0 * self%e11; self%e22 = self%e21

    else if (self%n1 < 0 .and. self%tau2 < 0) then

      t012s = t0 + t1 - t2 - self%ts

      self%b   = 2 * (self%mu0 - self%muH)/t012s
      self%cnt = ( self%muH * (t0 + t1) - self%mu0 * (t2 + self%ts) )/t012s
      self%c11 = ( 2 * self%mu0 * t01 - b * t0**2 )/t01/2
      self%d11 = b * t0/t01
      self%e11 = - b/t01/2
      self%c21 = (self%muH * t2**2 - self%cnt * self%ts**2)/t2s2
      self%d21 = 2 * self%ts * (self%cnt - self%muH)/t2s2
      self%e21 = (self%muH - self%cnt)/t2s2

      self%e22 = self%e21; self%c12 = self%c11; self%d12 = self%d11; self%e12 = self%e11
      self%c22 = self%c21; self%d22 = self%d21

    else

      self%c11 = (  t0**2 * ( 4 * self%cnt + b * (t0 + 3 * t1) ) - &
      2 * self%mu0 * (t0**2 + 2 * t0 * t1 - t1**2) )/t01**2/2

      self%c12 = ( 2 * self%cnt * (t0**2 - 2 * t0 * t1 - t1**2) - &
      t1**2 * (  b * (3 * t0 + t1) - 4 * self%mu0 )  )/t01**2/2

      self%d12 = ( 4 * (self%cnt - self%mu0) * t1 + b * (t0**2 + t0 * t1 + 2 * t1**2) )/t01**2
      self%e11 =   ( 4 * (self%cnt - self%mu0) + b * (t0 + 3 * t1) )/t01**2/2
      self%e12 = - ( 4 * (self%cnt - self%mu0) + b * (3 * t0 + t1) )/t01**2/2
      self%d11 = - 2 * t0 * self%e11

      self%c21 = - (  2 * self%cnt * (t2**2 + 2 * t2 * self%ts - self%ts**2) + &
      t2**2 * ( b * (t2 + 3 * self%ts) - 4 * self%muH )  )/t2s**2/2

      self%c22 =  (  2 * self%muH * (t2**2 - 2 * t2 * self%ts - self%ts**2) + &
      self%ts**2 * ( 4 * self%cnt + b * (3 * t2 + self%ts) )  )/t2s**2/2

      self%d21 = ( 4 * (self%cnt - self%muH) * t2 + b * (2 * t2**2 + &
                         t2 * self%ts + self%ts**2) )/t2s**2

      self%e21 = - ( 4 * (self%cnt - self%muH) + b * (t2 + 3 * self%ts) )/t2s**2/2
      self%e22 =   ( 4 * (self%cnt - self%muH) + b * (3 * t2 + self%ts) )/t2s**2/2
      self%d22 = - 2 * self%ts * self%e22

    end if

    call JunctionSet(self%R0, self%cnt, 0._dp, b, t0, t1, self%c31, self%c32, self%d31, &
                     self%d32, self%e31, self%e32)

    call JunctionSet(1 + self%eJ * (self%ts - t0)**2, 1 + self%eJ * (self%ts**2 - t1**2) ,&
    0._dp, - 2 * self%eJ * (self%ts - t1), t0, t1, self%a31, self%a32, self%b31, self%b32,&
    self%x31, self%x32)

    call self%MassIntersection(self%SoftIntersect, self%JetIntersect)

    end select

  end subroutine setQ

!ccccccccccccccc

  real (dp) function HardScale(self)
    class (Profiles), intent(in) :: self
    HardScale = self%muH
  end function HardScale

!ccccccccccccccc

  real (dp) function Energy(self)
    class (Profiles), intent(in) :: self
    Energy = self%Q
  end function Energy

!ccccccccccccccc

  real (dp) function MassScale(self)
    class (ProfilesPythia), intent(in) :: self
    MassScale = self%muM
  end function MassScale

!ccccccccccccccc

  subroutine setMassMuM(self, mass, muM)
    class (ProfilesPythia), intent(inout) :: self
    real (dp)             , intent(in   ) :: mass, muM
      self%mass = mass; self%muM = muM;
  end subroutine setMassMuM

!ccccccccccccccc

  pure subroutine Scales(self, tau, muNS, muH, muJ, muS, R)
    class (Profiles), intent(in ) :: self
    real (dp)       , intent(in ) :: tau
    real (dp)       , intent(out) :: muNS, muH, muJ, muS, R
    real (dp)                     :: muJ0, taus2, tau02, taut0

    muH = self%muH;   tauS2 = (tau - self%ts)**2;   tauT0 = (tau - self%t0)**2

    select type (self)
    type is (ProfilesMassless)

!cccccc soft scale before trumpet

      if (tau < 0) then
        muS = 0
      else if (tau <= self%t0) then
        muS = self%mu0
      else if (tau <= self%tmid1) then
        muS = self%c11 + self%d11 * tau + self%e11 * tau**2
      else if (tau <= self%t1) then
        muS = self%c12 + self%d12 * tau + self%e12 * tau**2
      else if (tau <= self%t2) then
        muS = self%b * tau + self%cnt
      else if (tau <= self%tmid2) then
        muS = self%c21 + self%d21 * tau + self%e21 * tau**2
      else if (tau <= self%ts) then
        muS = self%c22 + self%d22 * tau + self%e22 * tau**2
      else
        muS = muH
      end if

!cccccc implementing trumpet scale

      if (tau <= self%ts) then
         muJ = (1 + self%eJ * tauS2) * sqrt(self%muH * muS)
        if (tau > self%t0) then
          muS = (1 + self%eS * tauS2 * tauT0) * muS
        end if
      else
         muS = muH;  muJ = muH
      end if

    type is (ProfilesPythia)

      tau02 = (self%t0 - self%ts)**2

!cccccc soft scale before trumpet

    if (tau <= self%t0) then
      muS = self%mu0;    muJ0 = sqrt(self%muH * self%mu0)
      muJ = ( 1 + self%eJ * tau02 ) * min( muJ0**2/self%muM, muJ0 )
    else if (tau <= self%tmid1) then
      muS  = self%c11 + self%d11 * tau + self%e11 * tau**2
      muJ0 = sqrt(self%muH * muS)
      muJ  = (self%a31 + self%b31 * tau + self%x31 * tau**2) * min(muJ0**2/self%muM, muJ0)
    else if (tau <= self%t1) then
      muS  = self%c12 + self%d12 * tau + self%e12 * tau**2
      muJ0 = sqrt(self%muH * muS)
      muJ  = (self%a32 + self%b32 * tau + self%x32 * tau**2) * min(muJ0**2/self%muM, muJ0)
    else if (tau <= self%t2) then
      muS = self%b * tau + self%cnt
      muJ0 = sqrt(self%muH * muS)
      muJ = ( 1 + self%eJ * tauS2 ) * min(muJ0**2/self%muM, muJ0)
    else if (tau <= self%tmid2) then
      muS = self%c21 + self%d21 * tau + self%e21 * tau**2
      muJ0 = sqrt(self%muH * muS)
      muJ = ( 1 + self%eJ * tauS2 ) * min(muJ0**2/self%muM, muJ0)
    else if (tau <= self%ts) then
      muS = self%c22 + self%d22 * tau + self%e22 * tau**2
      muJ0 = sqrt(self%muH * muS)
      muJ = ( 1 + self%eJ * tauS2 ) * min(muJ0**2/self%muM, muJ0)
    else
      muS = self%muH; muJ = min(self%muH**2/self%muM, self%muH)
    end if

    end select

!cccccc nonsingular scale

    if (self%ns == -1) then
      muNS = (muH + muJ)/2
    else if (self%ns == 1) then
      muNS = (3 * muH - muJ)/2
    else if (self%ns == 0) then
      muNS = muH
    else
      muNS = 0
    end if

!cccccc r scale

    if (tau < 0) then
      r = 0
    else if (tau < self%t0) then
      r = self%R0
    else if (tau < self%tmid1) then
      r = self%R0 + self%e31 * tauT0
    else if (tau < self%t1) then
      r = self%c32 + self%d32 * tau + self%e32 * tau**2
    else
      r = muS
    end if

    select type (self)
    type is (ProfilesMassless)
     if ( tau > self%t1 .and. tau < self%tR) r = ( 1 + self%eR * (tau - self%tR)**2 ) * muS
    end select

  end subroutine Scales

!ccccccccccccccc

  function Sectors(self) result(res)
    class(ProfilesPythia), intent(in ) :: self
    real (dp), dimension(2)            :: res

    res = [ self%SoftIntersect, self%JetIntersect ]

  end

!ccccccccccccccc

  subroutine MassIntersection(self, resJet, resSoft)
    class(ProfilesPythia), intent(in ) :: self
    real (dp)            , intent(out) :: resJet, resSoft
    integer                            :: IFLAG
    real (dp)                          :: a, b, c

    a = 0; b = 1; resJet = 0; resSoft = 0

    if ( solveJet(0._dp) * solveJet(1._dp) < 0 ) then
      call DFZERO(solveJet, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG);  resJet = b
    end if

      a = 0; b = 1

    if ( solveSoft(0._dp) * solveSoft(1._dp) < 0 ) then
      call DFZERO(solveSoft, a, b, c, 1e-9_dp, 1e-9_dp, IFLAG);  resSoft = b
    end if

    contains

!ccccccccccccccc

    real (dp) function solveJet(tau)
      real (dp), intent(in) :: tau
      real (dp)             :: muNS, muH, muJ, muS, R

      call self%Scales(tau, muNS, muH, muJ, muS, R);  solveJet = muJ - self%muM

    end function solveJet

!ccccccccccccccc

     real (dp) function solveSoft(tau)
      real (dp), intent(in) :: tau
      real (dp)             :: muNS, muH, muJ, muS, R

      call self%Scales(tau, muNS, muH, muJ, muS, R);  solveSoft = muS - self%muM

    end function solveSoft
  end subroutine MassIntersection

!ccccccccccccccc

   pure Subroutine JunctionSet(mu0, cnt, b, b2, t1, t2, c1, c2, d1, d2, e1, e2)
     real (dp), intent(in ) :: mu0, cnt, b, b2, t1, t2
     real (dp), intent(out) :: c1, c2, d1, d2, e1, e2
     real (dp)              :: mu02, cnt1, dift, difa, difb

     mu02 = t1**2 + 2 * t1 * t2 - t2**2; difa = mu0 - cnt
     cnt1 = t2**2 + 2 * t1 * t2 - t1**2; dift = (t1 - t2)**2; difb = b - b2

     c1 = (   t1**2 * ( 4 * cnt - difb * (t1 + 3 * t2) ) - 2 * mu0 * mu02  )/dift/2
     c2 = (   t2**2 * ( 4 * mu0 + difb * (t2 + 3 * t1) ) - 2 * cnt * cnt1 )/dift/2
     d1 = (  t1 * (  2 * b  * t1 - b2 * t1 + 4 * difa ) + (b  - 3 * b2) * t1 * t2 + b * t2**2  )/dift
     d2 = (  t2 * (  2 * b2 * t2 - b * t2  - 4 * difa ) + (b2 - 3 * b ) * t1 * t2 + b2 * t1**2 )/dift
     e1 = - ( 4 * difa + difb * (t1 + 3 * t2) )/dift/2
     e2 =   ( 4 * difa + difb * (t2 + 3 * t1) )/dift/2

   end Subroutine JunctionSet

!ccccccccccccccc

end module ProfilesClass

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine f90MassIntersection(Q, beta, mu0, deltaLambda, n0, delta0, n1, delta1, t2, &
                               ts, slope, cnt, eH, eS, eJ, mass, mum, def, EShape, res)
  use constants; use ProfilesClass;  implicit none
  character (len = *)    , intent(in ) :: EShape, def
  real (dp)              , intent(in ) :: Q, mu0, n0, delta0, n1, t2, cnt, ts, slope, eH, &
                                         eJ, mass, beta, delta1, deltaLambda, muM, eS
  real (dp), dimension(2), intent(out) :: res
  type (ProfilesPythia)                :: Prof

  Prof = ProfilesPythia(Q, beta, mu0, deltaLambda, mu0, n0, delta0, n1, delta1, t2, ts, &
                        slope,  cnt, eH, eS, eJ, mass, muM, 0, def, EShape)
  res = [ Prof%JetIntersect, Prof%SoftIntersect ]

end subroutine f90MassIntersection

!ccccccccccccccc

subroutine f90RunningScales(Q, beta, mu0, deltaLambda, R0, n0, delta0, n1, delta1, t2, &
                            ts, slope, cnt, eH, eJ, ns, tau, res)
  use constants; use ProfilesClass; implicit none
  integer                 , intent(in) :: ns
  real (dp)               , intent(in) :: Q, mu0, R0, n0, delta0, n1, delta1, t2, cnt, &
                                          eH, eJ, beta, deltaLambda, tau, ts, slope
  real (dp), dimension(5), intent(out) :: res
  real (dp)                            :: muS, muJ, muH, muNS, R
  type (ProfilesPythia)                :: Prof

  Prof = ProfilesPythia(Q = Q, beta = beta, mu0 = mu0, deltaLambda = deltaLambda, R0 = R0,&
  n0 = n0, delta0 = delta0, n1 = n1, delta1 = delta1, tau2 = t2, ts = ts, slope = slope,  &
  cnt = cnt, eH = eH, eJ = eJ, ns = ns)

  call Prof%Scales(tau, muNS, muH, muJ, muS, R)

  res = [ muH, muJ, muS, muNS, R ]

end subroutine f90RunningScales

!ccccccccccccccc
