
module MassiveNSClass
  use MatrixElementsClass; use RunningClass; use ElectroWeakClass; use AnomDimClass
  use ModelClass         ; use carlson_elliptic_module           ; use GapMassClass
  use Constants, only: dp, Pi, Pi2, Pio2, d1mach, sr2, l2, prec  ; use MCtopClass
  use QuadPack, only: qags; implicit none                        ; private

  real (dp), parameter :: mth = 0.39307568887871164_dp, mcrit = 0.29913473867076185_dp

!ccccccccccccccc

  type, public                :: MassiveNS
    private
    class (MatricesElementsMass), allocatable :: matEl
    type (ElectroWeak)        :: EWeak
    type (Running)            :: alphaMass, alphaMassNl
    type (MCtop)              :: MC
    real (dp), dimension(2)   :: EW, Dirac
    real (dp), dimension(3)   :: deltaM, deltaM2
    character (len = 10)      :: shape, ES, current, singular, scheme
    integer                   :: nf, massOrd
    real (dp)                 :: m4, m6, m8, m10, alphaJ, tmin, tmax, mPole,   &
    v, lm, m, thrustMin, cp, mm, A1Singular, tint, B1Singular, alphaMu, B1, Q, &
    alphaM, muM, GlueMax, GlueInt, width, AMS, alphaH, m2, mass, ESfac

   contains

   final                         :: delete_object
   procedure, pass(self), public :: JetMassExp, FOMass, HJMNSMass, SetMCharm, &
   MassDelta, MassDeltaShift, numFlav, matElementNS, SetEverything, SetAlpha, &
   SetMBottom, EShape, MassVar, SetMTop, SetMasses, SetAll

   procedure, pass(self)         :: Quark, Glue, EWAdd, f1Q, f2Q, f3Q, f4Q, B1NS, &
   MassSing1loop, CumMassSing1loop, A1loop, FunCp, A0MS, CparamFOMass, Hcorr, A0, &
   CoefsCparam, NSMassMod, NSMassList, A1NSDiff

   generic   , public            :: NSMass => NSMassMod, NSMassList

  end type MassiveNS

!ccccccccccccccc

  type, extends (MassiveNS), public  :: MassiveScales
    private

    contains

    procedure, pass(self)  , public  :: HJMNSMassScales, orderMass, matElementScales

  end type MassiveScales

!ccccccccccccccc

  interface MassiveNS
    module procedure InMassNS
  end interface MassiveNS

!ccccccccccccccc

  interface MassiveScales
    module procedure InScales
  end interface MassiveScales

  contains

!ccccccccccccccc

 subroutine delete_object(this)
   type (MassiveNS) :: this
     if ( allocated(this%MatEl ) ) deallocate(this%MatEl)
  end subroutine delete_object

!ccccccccccccccc

   type (MassiveScales) function InScales(shape, ES, scheme, singular, current, &
   massOrd, matEl, EW)
     character (len = *)        , intent(in) :: shape, ES, current, singular, scheme
     type (MatricesElementsMass), intent(in) :: matEl
     type (ElectroWeak)         , intent(in) :: EW
     integer                    , intent(in) :: massOrd
     type (AnomDim)                          :: andim

     InScales%alphaMass = matEl%run()    ; andim = InScales%alphaMass%adim()
     InScales%shape     = shape(:6)
     InScales%ES        = ES             ; InScales%singular    = singular
     InScales%current   = current        ; InScales%alphaMassNl = matEl%runNl()
     InScales%nf        = andim%numFlav()

     InScales%massOrd = massOrd;  InScales%scheme = scheme; InScales%EWeak = EW

     allocate( MatricesElementsMass :: InScales%MatEl )

     select type (selector => InScales%MatEl)
     type is (MatricesElementsMass)
       selector = matEl
     end select

     call InScales%SetMasses()

     InScales%ESfac = 1; if ( shape(:6) == 'Cparam' ) InScales%ESfac = 6

   end function InScales

!ccccccccccccccc

   type (MassiveNS) function InMassNS(shape, ES, scheme, singular, current, massOrd,  &
   matEl, EW, width, mu)
     character (len = *)      , intent(in) :: shape, scheme, ES, current, singular
     real (dp), optional      , intent(in) :: mu, width
     type (MatrixElementsMass), intent(in) :: matEl
     type (ElectroWeak)       , intent(in) :: EW
     integer                  , intent(in) :: massOrd
     type (AnomDim)                        :: andim
     type (Running)                        :: alphaMass

     alphaMass          = matEl%run()   ; InMassNS%alphaMassNl = matEl%runNl()
     InMassNS%alphaMass = alphaMass     ; andim = alphaMass%adim()
     InMassNS%alphaJ    = matEl%alphaScale('jet') ;  InMassNS%shape   = shape(:6)
     InMassNS%alphaH    = matEl%alphaScale('hard');  InMassNS%alphaMu = 0
     InMassNS%Q         = matEl%scales('Q')   ; InMassNS%scheme  = scheme
     InMassNS%singular  = singular       ; InMassNS%current = current;  InMassNS%ES = ES
     InMassNS%nf        = andim%numFlav(); InMassNS%muM = matEl%scales('muM')
     InMassNS%alphaM    = matEl%alphaScale('massNf'); InMassNS%massOrd  = massOrd

     if ( present(mu) ) then
       InMassNS%alphaMu = alphaMass%alphaQCD(mu)/Pi
       call matEl%MassDeltaScheme(scheme, massOrd, InMassNS%mass, InMassNS%deltaM)
       if ( present(width) ) then
         call InMassNS%SetAll(width)
       else
         call InMassNS%SetAll(0._dp)
       end if
     end if

     InMassNS%MC = MCtop(shape, InMassNS%mass, InMassNS%Q)

     InMassNS%EW = EW%EWfactors(InMassNS%nf, InMassNS%Q)
     InMassNS%EW = InMassNS%EW/Sum(InMassNS%EW)

     allocate( MatrixElementsMass :: InMassNS%MatEl )

     select type (selector => InMassNS%MatEl)
     type is (MatrixElementsMass)
       selector = matEl
     end select

     call InMassNS%SetMasses()

     InMassNS%ESfac = 1; if ( shape(:6) == 'Cparam' ) InMassNS%ESfac = 6

   end function InMassNS

!ccccccccccccccc

  integer function orderMass(self)
    class (MassiveScales), intent(in) :: self
    orderMass = self%massOrd
  end function orderMass

!ccccccccccccccc

  subroutine SetAlpha(self, alpha)
    class (MassiveNS), intent(inout) :: self
    real (dp)        , intent(in)    :: alpha

    call self%alphaMass%SetAlpha(alpha); call self%alphaMassNl%SetAlpha(alpha)
    call self%MatEl%SetAlpha(alpha)

  end subroutine SetAlpha

!ccccccccccccccc

  subroutine SetMTop(self, m, mu)
    class (MassiveNS), intent(inout) :: self
    real (dp)        , intent(in)    :: m, mu

    call self%alphaMass%SetMTop(m, mu); call self%alphaMassNl%SetMTop(m, mu)
    call self%MatEl%SetMTop(m, mu)

  end subroutine SetMTop

!ccccccccccccccc

  subroutine SetMBottom(self, m, mu)
    class (MassiveNS), intent(inout) :: self
    real (dp)        , intent(in)    :: m, mu

    call self%alphaMass%SetMBottom(m, mu); call self%alphaMassNl%SetMBottom(m, mu)
    call self%MatEl%SetMBottom(m, mu)

  end subroutine SetMBottom

!ccccccccccccccc

  subroutine SetMCharm(self, m, mu)
    class (MassiveNS), intent(inout) :: self
    real (dp)        , intent(in)    :: m, mu

    call self%alphaMass%SetMCharm(m, mu); call self%alphaMassNl%SetMCharm(m, mu)
    call self%MatEl%SetMCharm(m, mu)

  end subroutine SetMCharm

!ccccccccccccccc

  subroutine SetMasses(self)
    class (MassiveNS), intent(inout) :: self

    self%mm = self%matEl%scales('mm')

    if ( self%scheme(:4) == 'pole' ) then
      self%mPole = self%alphaMass%scales('mL')
    else
      self%mPole = self%alphaMass%PoleMass(self%massOrd, self%mm)
    end if

  end subroutine SetMasses

!ccccccccccccccc

  subroutine SetAll(self, width)
    class (MassiveNS), intent(inout) :: self
    real (dp)        , intent(in)    :: width
    real (dp)                        :: m, m2, v2, m8, m6, m4, root

     m = self%mass; self%deltaM = self%deltaM/m; self%cp = 0; m = m/self%Q
     m2 = m**2; self%GlueInt = m/(1 - m); m4 = m2**2; self%m = m
     self%m2 = m2; v2 = 1 - 4 * m2;  self%v = sqrt(v2); self%lm = log(m)
     self%m4 = m4; m6 = m2**3; self%m10 = m2**5;  self%width = 2 * width * m
     self%m6 = m6; m8 = m4**2; self%m8  = m8; self%thrustMin = 1 - self%v

     if ( self%shape(:6) /= 'Cparam' ) then
       root = Sqrt(1 - 3 * m2); self%GlueMax = (5 - 4 * root)/3
     end if

     if ( self%ES(:1) /= 'Q' ) then
       self%tmin = 0; self%deltaM2 = 0
     else

       if ( self%shape(:6) == 'thrust' ) then

         self%tmin = self%thrustMin;  self%tint = self%GlueInt;  self%tmax = self%GlueMax

         self%deltaM2 = 2 * m**2 * [  2 * self%deltaM(1), ( self%deltaM(1)**2 &
         + 2 * self%v**2 * self%deltaM(2) )/v2, ( self%deltaM(1) * &
         self%deltaM(2) + 2 * m2 * ( self%deltaM(1)**3 - 2 * self%deltaM(1) * &
         self%deltaM(2) - 4 * self%deltaM(3) ) + self%deltaM(3) + 16 * m**4 * &
         self%deltaM(3) )/v2**2  ]/self%v

       else if ( self%shape(:3) == 'HJM' ) then

         self%tmin = m2;  self%tint = m * (1 - m - m2)/(1 - m)
         self%tmax = (2 * root - 1)/3 + m2;  self%width = self%width/2

         self%deltaM2 = m**2 * [ 2 * self%deltaM(1), 2 * self%deltaM(2) &
         + self%deltaM(1)**2, self%deltaM(1) * self%deltaM(2) +  self%deltaM(3)]

       else if ( self%shape(:3) == 'SJM' ) then

         self%tmin = 2 * m2;  self%tint = m * (1 - 2 * m2)/(1 - m)
         self%tmax = (2 * root - 1)/3 + 2 * m2

         self%deltaM2 = 2 * m**2 * [ 2 * self%deltaM(1), 2 * self%deltaM(2) &
         + self%deltaM(1)**2, self%deltaM(1) * self%deltaM(2) + self%deltaM(3) ]

       else if ( self%shape(:6) == 'Cparam' ) then

         self%tmin = 12 * m2 * (1 - m2); self%tint = 4 * m2 * (1 + 2 * m2)/(1 + 4 * m2)**2
         self%cp   = ( 2 * m2 - 2 * m4 + Sqrt(m2 - 4 * m6 + 4 * m8) )/2

         self%deltaM2 = 2 * m**2 * [ 2 * (1 - 2 * m2) * self%deltaM(1), &
         self%deltaM(1)**2 * (1 - 6 * m2) + 2 * (1 - 2 * m2) * self%deltaM(2),  &
         2 * self%deltaM(1) *  self%deltaM(2) - 4 * m2 * self%deltaM(1)**3  &
         - 12 * m2 * self%deltaM(1) * self%deltaM(2) + 2 * self%deltaM(3)   &
         -  4 * m2 * self%deltaM(3) ]

         if (m >= mth) then
           self%tmax = self%tint
         else
           self%tmax = (1 + 16 * m2 + 32 * m4)/(1 + 2 * m2)**2/8
         end if
       end if
     end if

     self%AMS = self%A0MS() ; self%B1 = self%B1NS()
     self%A1Singular = (4 * self%lm + 16 * self%lm**2)/3

     if ( self%shape(:6) /= 'Cparam' ) self%A1Singular = self%A1Singular + 6.579736267392905_dp
     if ( self%shape(:6) == 'Cparam' ) self%A1Singular = self%A1Singular + 8.772981689857206_dp

     self%B1Singular = - 8 * ( 1 + 2 * self%lm )/3
     self%Dirac = [ self%A0(), self%A1loop() ]

  end subroutine SetAll

!ccccccccccccccc

  subroutine SetEverything(self, mu, Q, muH, muJ, muS, R, Rmass, muM, width)
    class (MassiveNS)  , intent(inout) :: self
    real (dp)          , intent(in   ) :: muJ, muS, R, Rmass, width, muM, Q, muH
    real (dp), optional, intent(in   ) :: mu
    real (dp), dimension(3)            :: alphaList
    real (dp)                          :: alphaJ

    self%EW     = self%EWeak%EWfactors(self%nf, Q); self%EW  = self%EW/Sum(self%EW)
    self%alphaM = self%alphaMass%alphaQCD(muM)/Pi ; self%muM = muM;
    self%alphaH = self%alphaMass%alphaQCD(muH)/Pi ; self%Q   = Q
    call self%matEl%SetScales('Q'    , Q  )  ; call self%matEl%SetScales('muH', muH)
    call self%matEl%SetScales('muM'  , muM)  ; call self%matEl%SetScales('muS', muS)
    call self%matEl%SetScales('muJ'  , muJ)  ; call self%matEl%SetScales('R'  , R  )
    call self%matEl%SetScales('Rmass', Rmass)

    if ( present(mu) ) self%alphaMu = self%alphaMass%alphaQCD(mu)/Pi

    if ( muS >= muM ) then
      alphaList = powList( self%alphaMass%alphaQCD(muS)/Pi, 3 )
    else
      alphaList = powList( self%alphaMassNl%alphaQCD(muS)/Pi, 3 )
    end if

    if ( muJ >= muM ) then
      self%alphaJ = self%alphaMass%alphaQCD(muJ)/Pi
    else
      self%alphaJ = self%alphaMassNl%alphaQCD(muJ)/Pi
    end if

    call self%matEl%MassDeltaSet(self%scheme, self%massOrd, muJ, Rmass, &
    self%mass, self%deltaM, alphaJ)

    call self%SetAll(width);  call self%matEl%SetDelta(muS, R, alphaList)

    self%MC = MCtop(self%shape, self%mass, Q)

  end subroutine SetEverything

!ccccccccccccccc

  real (dp) function NSMassMod(self, Mod, setup, gap, cum, order, run, R0, mu0, &
  delta0, h, t, t2)
    class (MassiveNS)    , intent(in) :: self
    character (len = *)  , intent(in) :: setup, gap, cum
    integer              , intent(in) :: order, run
    real (dp)            , intent(in) :: R0, mu0, delta0, h, t
    real (dp), optional  , intent(in) :: t2
    type (Model)         , intent(in) :: Mod
    type (GapMass)                    :: GapMassive
    logical                           :: dobsing
    integer                           :: neval, ier, cumul
    real (dp), dimension(order,order) :: delt
    real (dp), dimension(order)       :: delta
    real (dp), dimension(1)           :: resul
    real (dp)                         :: abserr, res, tau, Modelo, ModPlus, p2, &
    shift, ModDer, tshift, tshift2, p3, tau2, p3shift, p2shift, p

    GapMassive = GapMass(run, gap, self%matEl, R0, mu0); NSMassMod = 0

    if ( present(t2) ) then; if (t2 <= t) return;  end if

    if ( setup(:5) == 'Model' ) then
      if ( present(t2) ) then
        resul = NSMassList(self, [Mod], setup, gap, cum, order, run, R0, mu0, delta0, h, t, t2)
      else
        resul = NSMassList(self, [Mod], setup, gap, cum, order, run, R0, mu0, delta0, h, t)
      end if
      NSMassMod = resul(1); return
    end if

    call GapMassive%setGammaShift( order, self%shape, delta0, h, &
    self%matEl%scales('R'), self%matEl%scales('muS'), shift, delt )

    delta = delt(1,:); Modelo = 0; ModPlus = 0; ModDer = 0

    if ( self%shape(:3) == 'HJM' ) then; delta = delta/2; shift = shift/2;end if

    tshift = t - shift/self%Q;  tau = tshift - self%tmin; p = self%Q * tau
    tau2   = tau; p3 = p; p2 = p; tshift2 = tshift; p3shift = self%Q * tshift

    if ( present(t2) ) then
      tshift2 = t2 - shift/self%Q;  tau2    = tshift2 - self%tmin
      p2 = self%Q * tau2         ;  p2shift = self%Q * tshift2
    end if

    if ( present(t2) .and. p < 0 ) then
      tau = tau2; p = p2; tshift = tshift2
    end if

    dobsing = ( .not. present(t2) .or. ( present(t2) .and. abs(p2 - p) <= d1mach(1) ) )

    cumul = 0; if ( cum(:3) == 'cum' ) cumul = - 1

    if (tau2 <= 0 .and. self%width <= d1mach(1) ) return

    setup_if: if ( setup(:2) == 'FO' .and. tau2 > 0) then

      cum_if: if ( cum(:4) == 'diff' ) then

        if (order > 0) then

          if ( dobsing ) then
            NSMassMod = self%alphaMu * self%FOMass(t)
          else
            NSMassMod = self%alphaMu * self%FOMass(t2, t)
          end if

        end if

      else if ( cum(:3) == 'cum' ) then

        if (  order >= 0 .and. dobsing ) NSMassMod = 1 + self%Dirac(1)

        order_if: if (order > 0) then

          if ( dobsing ) then

            call qags( NS1lloop, 0._dp, tau, prec, prec, res, abserr, neval, ier )

            NSMassMod = NSMassMod + self%alphaMu * (  res + self%A1Singular + self%Dirac(2) + &
            ( self%B1Singular + self%B1 ) * log(tau/self%ESfac)  ) + self%deltaM(1) * self%AMS

         else

           call qags( NS1lloop, tau, tau2, prec, prec, res, abserr, neval, ier )

            NSMassMod = self%alphaMu * ( res + (self%B1Singular + self%B1) * log(tau2/tau) )

         end if

        end if order_if
      end if   cum_if

    else if ( setup(:7) == 'NoModel') then

      cum_if_2: if ( cum(:4) == 'diff' .and. tau2 > 0 ) then

        if (order > 0) then

          NSMassMod = self%FOMass(tshift) - self%MassSing1loop(tshift)
          if ( .not. dobsing ) NSMassMod = self%FOMass(tshift2) - self%MassSing1loop(tshift2) - NSMassMod

          if ( self%singular(:3) == 'abs' .or. self%width >= d1mach(1) ) then

            if ( dobsing ) then
              NSMassMod = NSMassMod - self%B1/tau
            else
              NSMassMod = NSMassMod + self%B1 * (1/tau - 1/tau2)
            end if

          end if

          NSMassMod = NSMassMod * self%alphaMu

        end if

      else if ( cum(:3) == 'cum'  .and. tau2 > 0 ) then

        if (  order >= 0 .and. self%singular(:3) /= 'abs' .and. self%width <= d1mach(1) &
        .and. dobsing  ) NSMassMod = self%Dirac(1)

        order_if_2: if (order > 0) then

          if ( dobsing ) then

            call qags( NS1lloop, 0._dp, tau, prec, prec, res, abserr, neval, ier )

            NSMassMod = NSMassMod + self%alphaMu * ( res - self%CumMassSing1loop(tshift) +  &
            self%A1Singular + self%B1Singular * log(tau/self%ESfac) )

            if ( self%singular(:3) /= 'abs' .and. self%width <= d1mach(1) )  &
            NSMassMod = NSMassMod + self%AMS * self%deltaM(1) + &
            self%alphaMu * ( self%Dirac(2) + self%B1 * log(tau/self%ESfac) )

          else

            call qags( NS1lloop, tau, tau2, prec, prec, res, abserr, neval, ier )

            NSMassMod = NSMassMod + self%alphaMu * ( res + self%CumMassSing1loop(tshift) -  &
            self%CumMassSing1loop(tshift2) +  self%B1Singular * log(tau2/tau) )

            if ( self%singular(:3) /= 'abs' .and. self%width <= d1mach(1) )  &
            NSMassMod = NSMassMod + self%alphaMu * self%B1 * log(tau2/tau)

          end if

        end if order_if_2
      end if   cum_if_2

      if ( self%singular(:3) /= 'abs' .and. self%width > tiny(1._dp) ) then

        if ( setup(8:15) /= 'Unstable' .or. self%shape(:6) == 'thrust' ) then

          if ( .not. present(t2) ) then

            Modelo = self%Q**(1 + cumul) * BreitWigner(self%width, cumul, p3)

            if (order > 0) then
              ModPlus = self%Q**(1 + cumul) * BreitWigner(self%width, cumul - 2, p3)
              ModDer  = self%Q**(2 + cumul) * BreitWigner(self%width, cumul + 1, p3)
            end if

          else

            Modelo = self%Q**(1 + cumul) * BreitWigner(self%width, cumul, p3, p2)

            if (order > 0) then
              ModPlus = self%Q**(1 + cumul) * BreitWigner(self%width, cumul - 2, p3, p2)
              ModDer  = self%Q**(2 + cumul) * BreitWigner(self%width, cumul + 1, p3, p2)
            end if
          end if

        else
          Modelo = 0; ModPlus = 0; ModDer = 0
        end if

        if ( setup(8:15) == 'Unstable' ) then
          if ( .not. present(t2) ) then

            Modelo = Modelo * self%MC%Delta() + self%Q**(1 + cumul) * &
            self%MC%BreitUnstable(self%width, cumul, p3shift)

            if (order > 0) then

              ModPlus = ModPlus * self%MC%Delta() + self%Q**(1 + cumul) * &
              self%MC%BreitUnstable(self%width, cumul - 2, p3shift)

              ModDer  = ModDer  * self%MC%Delta() + self%Q**(2 + cumul) * &
              self%MC%BreitUnstable(self%width, cumul + 1, p3shift)

            else

              ModPlus = ModPlus * self%MC%Delta() + self%Q**(1 + cumul) * &
              self%MC%BreitUnstable(self%width, cumul - 2, p3shift, p2shift)

              ModDer  = ModDer  * self%MC%Delta() + self%Q**(2 + cumul) * &
              self%MC%BreitUnstable(self%width, cumul + 1, p3shift, p2shift)

            end if
          end if
        end if
      end if
    end if   setup_if

    if (  self%singular(:3) /= 'abs' .and. ( setup(:5) == 'Model' .or. &
         ( setup(:7) == 'NoModel' .and. self%width > tiny(1._dp) ) )   ) then

      if (order >= 0) NSMassMod = NSMassMod + self%Dirac(1) * Modelo

      if (order > 0) then

        ModPlus = ModPlus - Modelo * log(self%Q)

        NSMassMod = NSMassMod + self%alphaMu * self%B1 * ModPlus + &
        ( self%AMS * self%deltaM(1) + self%alphaMu * self%Dirac(2) ) * &
        Modelo - ModDer * self%Dirac(1) * self%deltaM2(1)

        if ( gap(:5) /= 'NoGap' ) NSMassMod = NSMassMod - delta(1)/self%Q * self%Dirac(1) * ModDer

      end if
    end if

    contains

!ccccccccccccccc

    real (dp) function NS1lloop(tau)
      real (dp) , intent(in) :: tau
      NS1lloop = self%FOMass(tau + self%tmin) - ( self%B1 + self%B1Singular )/tau
    end function NS1lloop

  end function NSMassMod

!ccccccccccccccc

  function NSMassList(self, ModList, setup, gap, cum, order, run, R0, mu0, &
  delta0, h, t, t2) result(resList)
    class (MassiveNS)    , intent(in)      :: self
    character (len = *)  , intent(in)      :: setup, gap, cum
    integer              , intent(in)      :: order, run
    real (dp)            , intent(in)      :: R0, mu0, delta0, h, t
    real (dp), optional  , intent(in)      :: t2
    type (Model), dimension(:), intent(in) :: ModList
    type (Model)                           :: Mod
    type (GapMass)                         :: GapMassive
    logical                                :: dobsing
    integer                                :: neval, ier, cumul, i
    real (dp), dimension(order,order)      :: delt
    real (dp), dimension(order)            :: delta
    real (dp), dimension( size(ModList) )  :: resList, Modelo, ModPlus, ModDer
    real (dp)                              :: abserr, tau, p, tau2, p2, p3, &
    shift, tshift, tshift2, p3shift, p2shift

    GapMassive = GapMass(run, gap, self%matEl, R0, mu0); resList = 0

    if ( present(t2) ) then; if (t2 <= t) return;  end if

    call GapMassive%setGammaShift( order, self%shape, delta0, h, &
    self%matEl%scales('R'), self%matEl%scales('muS'), shift, delt )

    delta = delt(1,:); Modelo = 0; ModPlus = 0; ModDer = 0

    if ( self%shape(:3) == 'HJM' ) then; delta = delta/2; shift = shift/2;end if

    tshift = t - shift/self%Q;  tau = tshift - self%tmin; p = self%Q * tau
    tau2   = tau; p3 = p; p2 = p; tshift2 = tshift; p3shift = self%Q * tshift

    if ( present(t2) ) then
      tshift2 = t2 - shift/self%Q;  tau2    = tshift2 - self%tmin
      p2 = self%Q * tau2         ;  p2shift = self%Q * tshift2
    end if

    if ( present(t2) .and. p < 0 ) then
      tau = tau2; p = p2; tshift = tshift2
    end if

    dobsing = ( .not. present(t2) .or. ( present(t2) .and. abs(p2 - p) <= d1mach(1) ) )

    cumul = 0; if ( cum(:3) == 'cum' ) cumul = - 1

    if (tau2 <= 0 .and. self%width <= d1mach(1) ) return

    do i = 1, size(modList)

      Mod = ModList(i)

      abs_if: if ( self%singular(:3) /= 'abs' .and. order >= 0) then

        if ( self%width <= d1mach(1) ) then

          if ( dobsing ) then

            Modelo(i)  = self%Q**(1 + cumul) * Mod%ShapeFun(cumul, p)

            if (order > 0) then
              ModPlus(i) = self%Q**(1 + cumul) * Mod%ShapeFun(cumul - 2, p)
              ModDer(i)  = self%Q**(2 + cumul) * Mod%ShapeFun(cumul + 1, p)
            end if

          else

            Modelo(i)  = self%Q**(1 + cumul) * Mod%ShapeFun(cumul, p, p2)

            if (order > 0) then
              ModPlus(i) = self%Q**(1 + cumul) * Mod%ShapeFun(cumul - 2, p, p2)
              ModDer(i)  = self%Q**(2 + cumul) * Mod%ShapeFun(cumul + 1, p, p2)
            end if

          end if

        else

          if ( setup(6:13) /= 'Unstable' .or. self%shape(:6) == 'thrust' ) then

            if ( .not. present(t2) ) then

              Modelo(i)  = self%Q**(1 + cumul) * Mod%BreitModel(self%width, cumul, p)

              if (order > 0) then
                ModPlus(i) = self%Q**(1 + cumul) * Mod%BreitModel(self%width, cumul - 2, p)
                ModDer(i)  = self%Q**(2 + cumul) * Mod%BreitModel(self%width, cumul + 1, p)
              end if

            else

              Modelo(i) = self%Q**(1 + cumul) * Mod%BreitModel(self%width, cumul, p3, p2)

              if (order > 0) then
                ModPlus(i) = self%Q**(1 + cumul) * Mod%BreitModel(self%width, cumul - 2, p3, p2)
                ModDer(i)  = self%Q**(2 + cumul) * Mod%BreitModel(self%width, cumul + 1, p3, p2)
              end if

            end if
          else
            Modelo(i) = 0; ModPlus(i) = 0; ModDer(i) = 0
          end if

          if ( setup(6:13) == 'Unstable' ) then

            if ( .not. present(t2) ) then

              Modelo(i) = self%MC%Delta() * Modelo(i) + self%Q**(1 + cumul) * &
              Mod%BreitModUns(self%MC, self%width, cumul, p3shift)

              if (order > 0) then

                ModPlus(i) = self%MC%Delta() * ModPlus(i) + self%Q**(1 + cumul) * &
                Mod%BreitModUns(self%MC, self%width, cumul - 2, p3shift)

                ModDer(i)  = self%MC%Delta() * ModDer(i) + self%Q**(2 + cumul) * &
                Mod%BreitModUns(self%MC, self%width, cumul + 1, p3shift)

              end if

            else

              Modelo(i) = self%MC%Delta() * Modelo(i) + self%Q**(1 + cumul) * &
              Mod%BreitModUns(self%MC, self%width, cumul    , p3shift, p2shift)

              if (order > 0) then

                ModPlus(i) = self%MC%Delta() * ModPlus(i) + self%Q**(1 + cumul) * &
                Mod%BreitModUns(self%MC, self%width, cumul - 2, p3shift, p2shift)

                ModDer(i)  = self%MC%Delta() * ModDer(i) + self%Q**(2 + cumul) * &
                Mod%BreitModUns(self%MC, self%width, cumul + 1, p3shift, p2shift)

              end if
            end if
          end if
        end if

      end if abs_if

      order_if_3: if (order > 0 .and. p2 > 0) then

        call qags( inteNSMod, 0._dp, p2, prec, prec, resList(i), abserr, neval, ier )

      end if order_if_3
    end do

    resList = self%alphaMu * resList * self%Q**cumul

    if ( self%singular(:3) /= 'abs' ) then

      if (order >= 0) resList = resList + self%Dirac(1) * Modelo

      if (order > 0) then

        ModPlus = ModPlus - Modelo * log(self%Q)

        resList = resList + self%alphaMu * self%B1 * ModPlus + &
        ( self%AMS * self%deltaM(1) + self%alphaMu * self%Dirac(2) ) * &
        Modelo - ModDer * self%Dirac(1) * self%deltaM2(1)

        if ( gap(:5) /= 'NoGap' ) resList = resList - delta(1)/self%Q * self%Dirac(1) * ModDer

      end if
    end if

    contains

!ccccccccccccccc

    real (dp) function inteNSMod(l)
      real (dp), intent(in) :: l
      real (dp)             :: t2, t3

      t3 = l/self%Q; t2 = t3 + self%tmin

      inteNSMod = self%FOMass(t2) - ( self%B1 + self%B1Singular )/t3 - &
      2 * ( t3/(t3 + self%m2)**2 - 4 * log( t3/self%m2 + 1 )/t3 )/3

      if ( dobsing ) then
        inteNSMod = inteNSMod * Mod%ShapeFun(cumul, p - l)
      else
        inteNSMod = inteNSMod * Mod%ShapeFun(cumul, p - l, p2 - l)
      end if

    end function inteNSMod

  end function NSMassList

!ccccccccccccccc

  real (dp) function HJMNSMassScales(self, c, modList, setup, gap, cum, order, &
  run, mu, Q, muJ, muS, R, Rmass, muM, width, R0, mu0, delta0, h, t, t2)
    class (MassiveScales)     , intent(inout) :: self
    character (len = *)          , intent(in) :: setup, gap, cum
    integer                      , intent(in) :: order, run
    real (dp)        , optional  , intent(in) :: t2
    real (dp)    , dimension(:,:), intent(in) :: c
    type (Model) , dimension(:,:), intent(in) :: modList
    real (dp)                    , intent(in) :: R0, mu0, delta0, h, t, mu, muJ, muS, Q, &
    Rmass, width, R, muM

    call self%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

    if ( present(t2) ) then

      HJMNSMassScales = self%HJMNSMass(c, modList, setup, gap, cum, order, run, R0, &
      mu0, delta0, h, t, t2)

    else

      HJMNSMassScales = self%HJMNSMass(c, modList, setup, gap, cum, order, run, R0, &
      mu0, delta0, h, t)

    end if

  end function HJMNSMassScales

!ccccccccccccccc

  real (dp) function HJMNSMass(self, c, modList, setup, gap, cum, order, run, R0, &
                                mu0, delta0, h, t, t2)
    class (MassiveNS)           , intent(in) :: self
    character (len = *)         , intent(in) :: setup, gap, cum
    integer                     , intent(in) :: order, run
    real (dp)                   , intent(in) :: R0, mu0, delta0, h, t
    real (dp)   , optional      , intent(in) :: t2
    real (dp)   , dimension(:,:), intent(in) :: c
    type (Model), dimension(:,:), intent(in) :: modList
    integer                                  :: neval, ier, cumul
    logical                                  :: dobsing
    real (dp), dimension(order,order)        :: delt
    real (dp), dimension(order)              :: delta
    real (dp)                                :: abserr, res, tau, Modelo, ModPlus, tau2, &
                                                shift, ModDer, tshift, tshift2, p, p2

    HJMNSMass = 0; if ( self%shape(:3) /= 'HJM' ) return

    if ( present(t2) ) then
      if (t2 <= t) return
    end if

    if ( setup(:2) == 'FO' .or. setup(:7) == 'NoModel' ) then

      if ( .not. present(t2) ) HJMNSMass = self%NSMass( ModList(1,1), setup, gap, cum, &
      order, run, R0, mu0, delta0, h, t )

      if ( present(t2) ) HJMNSMass = self%NSMass( ModList(1,1), setup, gap, cum, order, &
      run, R0, mu0, delta0, h, t, t2 )

      return
    end if

    call self%matEl%setGammaShift( order, run, gap, self%shape, delta0, R0, mu0, h, shift, delt )

    delta = delt(1,:)/2; tshift = t - shift/self%Q/2;  tau = tshift - self%tmin
    p = self%Q * tau;  tshift2 = tshift; tau2 = tau;  p2 = p

    if ( present(t2) ) then
      tshift2 = t2 - shift/self%Q/2; tau2 = tshift2 - self%tmin; p2 = self%Q * tau
    end if

    if ( present(t2) .and. p < 0 ) then
      tau = tau2; p = p2; tshift = tshift2
    end if

    dobsing = ( .not. present(t2) .or. ( present(t2) .and. abs(p2 - p) <= d1mach(1) ) )

    if ( self%width <= d1mach(1) .and. tau2 <= 0 ) return

    cumul = 0; if ( cum(:3) == 'cum' ) cumul = - 1

    Modelo = 0; ModPlus = 0; ModDer = 0

    model_if: if ( self%singular(:3) /= 'abs' .and. order >= 0 ) then

      if ( self%width <= d1mach(1) ) then

        Modelo  = (2*self%Q)**(1 + cumul) * Model2D(c, modList, cumul, - 1, p, p)
        if ( .not. dobsing ) Modelo = (2*self%Q)**(1 + cumul) * &
                                       Model2D(c, modList, cumul, - 1, p2, p2) - Modelo

        if (order > 0) then

           ModPlus = self%Q**(1 + cumul)/2**(-cumul) * (        &
           Model2D(c, modList, - 2 + cumul, - 1  , p, p) +      &
           Model2D(c, modList, - 3        , cumul, p, p) )

           ModDer  = self%Q**(2 + cumul) * 2**(1 + cumul) * (  &
           Model2D(c, modList, cumul,         0, p, p) +       &
           Model2D(c, modList, -1   , 1 + cumul, p, p) )

           if ( .not. dobsing ) then

             ModPlus = self%Q**(1 + cumul)/2**(-cumul) * (                 &
             Model2D(c, modList, - 2 + cumul, - 1  , p2, p2) +             &
             Model2D(c, modList, - 3        , cumul, p2, p2) ) - ModPlus

             ModDer  = self%Q**(2 + cumul) * 2**(1 + cumul) * (            &
             Model2D(c, modList, cumul,         0, p2, p2) +               &
             Model2D(c, modList, -1   , 1 + cumul, p2, p2) ) - ModDer

           end if
        end if

        else

        Modelo = (2*self%Q)**(1 + cumul) * BreitModel2D(self%width, c, modList, &
        cumul, - 1, p, p)

        if ( .not. dobsing ) Modelo = (2*self%Q)**(1 + cumul) * BreitModel2D(self%width,  &
        c, modList, cumul, - 1, p2, p2) - Modelo

          if (order > 0) then

             ModPlus = self%Q**(1 + cumul)/2**(-cumul) * (                      &
             BreitModel2D(self%width, c, modList, - 2 + cumul, - 1  , p, p)  +  &
             BreitModel2D(self%width, c, modList, - 3        , cumul, p, p) )

             ModDer  = self%Q**(2 + cumul) * 2**(1 + cumul) * (                 &
             BreitModel2D(self%width, c, modList, cumul,         0, p, p)  +    &
             BreitModel2D(self%width, c, modList, -1   , 1 + cumul, p, p) )

             if ( .not. dobsing ) then

             ModPlus = self%Q**(1 + cumul)/2**(-cumul) * (                        &
             BreitModel2D(self%width, c, modList, - 2 + cumul, - 1  , p2, p2)  +  &
             BreitModel2D(self%width, c, modList, - 3        , cumul, p2, p2) ) - ModPlus

             ModDer  = self%Q**(2 + cumul) * 2**(1 + cumul) * (                   &
             BreitModel2D(self%width, c, modList, cumul,         0, p2, p2)  +    &
             BreitModel2D(self%width, c, modList, -1   , 1 + cumul, p2, p2) ) - ModDer

             end if
          end if
        end if

      HJMNSMass = self%Dirac(1) * Modelo;   ModPlus = ModPlus - Modelo * log(self%Q)

    end if model_if

    order_if: if (order > 0) then

      call qags( integrand, 0._dp, p2, prec, prec, res, abserr, neval, ier )

      HJMNSMass = HJMNSMass + self%alphaMu * self%Q**cumul * res

      if ( self%singular(:3) /= 'abs' ) HJMNSMass = HJMNSMass + self%alphaMu * self%B1  * &
      ModPlus + ( self%AMS * self%deltaM(1) + self%alphaMu * self%Dirac(2) ) * &
      Modelo - 2 * ModDer * self%m2/(1 - self%tmin) * self%Dirac(1) * self%deltaM(1)

      if ( gap(:5) /= 'NoGap' ) HJMNSMass = HJMNSMass - delta(1)/self%Q * self%Dirac(1) * ModDer

    end if order_if

    contains

    real (dp) function integrand(l)
      real (dp), intent(in) :: l

      real (dp)             :: tt3, t3, p3, p4

      t3 = l/self%Q; p3 = self%Q * tshift; tt3 = t3 + self%tmin

      if (dobsing) then

        integrand = 4 * self%Glue(tt3)/2**(1 - cumul)                * &
                  ( Model2D(c, modList, cumul,    -1, p - l, p3)     + &
                    Model2D(c, modList, -1   , cumul, p - l, p3) )/3

        integrand = integrand + (  4 * self%Quark(tt3)/3 - ( self%B1    + &
        self%B1Singular )/t3 - ( t3/(t3 + self%m2)**2 - 4/t3 * log(t3/self%m2 + 1) )/3  )&
         * 2 * ( Model2D(c, modList, cumul, -1, p - l, p)                        + &
                 Model2D(c, modList, -1, cumul, p - l, p) )/2**(-cumul)
      else

        p4 = self%Q * tshift2

        integrand = 4 * self%Glue(tt3)/2**(1 - cumul)                 * &
        ( Model2D(c, modList, cumul,    -1, p2 - l, p4)     - &
        Model2D(c, modList, cumul,    -1, p  - l, p3)     + &
        Model2D(c, modList, -1   , cumul, p2 - l, p4)     - &
        Model2D(c, modList, -1   , cumul, p  - l, p3) )/3

        integrand = integrand + (  4 * self%Quark(tt3)/3 - ( self%B1 + self%B1Singular &
         )/t3 - 2 * ( t3/(t3 + self%m2)**2 - 4/t3 * log(t3/self%m2 + 1) )/3  )&
        /2**(-cumul) * ( Model2D(c, modList, cumul, -1, p2 - l, p2) - &
        Model2D(c, modList, cumul, -1, p  - l, p ) + &
        Model2D(c, modList, -1, cumul, p2 - l, p2) - Model2D(c, modList, -1, cumul, p - l, p) )

      end if

    end function integrand

  end function HJMNSMass

!ccccccccccccccc

  recursive real (dp) function FOMass(self, t, t2) result(res)
    class (MassiveNS)  , intent(in) :: self
    real (dp)          , intent(in) :: t
    real (dp), optional, intent(in) :: t2

    res = 0

    if ( present(t2) ) then
      if ( t2 > t ) res = self%FOMass(t2) - self%FOMass(t);  return
    end if

    if ( self%m < d1mach(1) ) then

      if ( self%shape(:6) == 'Cparam' ) then
        res = CParamFO1Loop(t)
      else
        res = ThrustFO1loop(t)
      end if

    else

      if ( self%shape(:6) == 'Cparam' ) then
        res = self%CparamFOMass(t)
      else
        res = 2 * ( 2 * self%Quark(t) + self%Glue(t) )/3
      end if

    end if

  end function FOMass

!ccccccccccccccc

   real (dp) function MassSing1loop(self, tau)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: tau
    real (dp)                     :: t

    t = (tau - self%tmin)/self%ESfac; MassSing1loop = 0; if (t <= 0) return

    MassSing1loop = ( - 8/t + 2 * t/(t + self%m2)**2 - 8 * log(t + self%m2)/t )/3/self%ESfac

   end function MassSing1loop

!ccccccccccccccc

   real (dp) function CumMassSing1loop(self, tau)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: tau
    real (dp)                     :: t, DiLog, dilogpiece

    t = (tau - self%tmin)/self%ESfac; CumMassSing1loop = 0;  if (t <= 0) return

    if (self%m > 1e-4_dp) then
       dilogpiece = 8 * DiLog(-t/self%m2) + 4 * Log(t/self%m2)**2
    else
       dilogpiece = - 13.15947253478581 - 2 * self%m6**2/t**6/9 + 0.32 * self%m10/t**5 &
              - self%m8/t**4/2 + 8 * self%m6/t**3/9 - 2 * self%m4/t**2 + 8 * self%m2/t
    end if

    CumMassSing1loop = ( 2 * log(t + self%m2) + dilogpiece - 2 * t/(t + self%m2) &
                             - 4 * log(t)**2 - 8 * log(t) )/3

    if ( self%shape(:6) /= 'Cparam' ) CumMassSing1loop = CumMassSing1loop + 19.739208802178716_dp/3
    if ( self%shape(:6) == 'Cparam' ) CumMassSing1loop = CumMassSing1loop + 26.318945069571620_dp/3

   end function CumMassSing1loop

!ccccccccccccccc

  real (dp) function Quark(self, t)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: t
    real (dp)                     :: tr

    tr = t; if ( self%shape(:3) == 'SJM' ) tr = t - self%m2

    Quark = 0

    if (tr > self%tmin .and. tr <= self%tint) then
      Quark = self%f2Q(tr)
    else if (tr > self%tint .and. tr <= self%tmax) then
      Quark = self%f3Q(tr)
    end if

  end function Quark

!ccccccccccccccc

  real (dp) function Glue(self, t)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: t

    Glue = 0

    if (t > 4 * self%m2 .and. t <= self%GlueInt) then
      Glue = self%f1Q(t)
    else if (t > self%GlueInt .and. t <= self%GlueMax) then
      Glue = self%f4Q(t)
    end if

  end function Glue

!ccccccccccccccc

  real (dp) function CparamFOMass(self, c2)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: c2
    real (dp), dimension(4)       :: Coef
    real (dp)                     :: c, ypt, ymt, yp, ym, denom1, root1, n2, &
    denom2, root2, Thetap, Thetam, Xi, Xi2, h2, k2

    c = c2/6;  denom1 = 2 * ( c + (c - 1) * self%m2 ); denom2 = 2 * (1 + c)
    root2 = Sqrt( 1 + 16 * self%m2 + 32 * self%m4 - 8 * c * (1 + 2 * self%m2)**2 )
    root1 = self%m * Sqrt( (1 + 8 * c) * self%m2 - 4 * c**2 + 4 * self%m6 &
               - 4 * (1 + 2 * c) * self%m4 )

    ypt = 2 * c - 3 * self%m2 + 2 * self%m4
    ymt = (ypt - root1)/denom1;  ypt = (ypt + root1)/denom1

    yp = 1 + 4 * c - 8 * self%m2
    ym = (yp - root2)/denom2;  yp = (yp + root2)/denom2

    Thetap = ASin(   Sqrt(  ( yp * (ypt - ym) )/( ypt * (yp - ym) )  )   )
    Thetam = ASin(   Sqrt(  ( yp * (ymt - ym) )/( ymt * (yp - ym) )  )   )

    CparamFOMass = 0

    if (self%m < 0 .or. 2 * self%m >= 1 .or. c2 <= self%tmin .or. c >= self%tmax) return

    call self%CoefsCparam(c, Coef, n2, h2, k2, Xi, Xi2)

    if (c < self%tint) then

      CparamFOMass = self%FunCp(Coef, n2, h2, k2, Xi, Xi2, c, Thetam)

    else if (c < self%cp .and. self%m < mcrit) then

      CparamFOMass = self%FunCp(Coef, n2, h2, k2, Xi, Xi2, c, Thetam) - &
                     self%FunCp(Coef, n2, h2, k2, Xi, Xi2, c, Thetap) + &
                     self%FunCp(Coef, n2, h2, k2, Xi, Xi2, c, Pio2  )

    else if ( ( (mcrit < self%m .and. self%m < mth) .and. c >= self%tint ) .or.  &
    (self%m <= mcrit .and. c >= self%cp)  ) then

      CparamFOMass = self%FunCp(Coef, n2, h2, k2, Xi, Xi2, c, Pio2)

    end if

    CparamFOMass = CparamFOMass/9

  end function CparamFOMass

!ccccccccccccccc

  subroutine CoefsCparam(self, c, Coef, n2, h2, k2, Xi, Xi2)
    class (MassiveNS)      , intent(in)  :: self
    real (dp)              , intent(in)  :: c
    real (dp)              , intent(out) :: n2, h2, k2, Xi, Xi2
    real (dp), dimension(4), intent(out) :: Coef
    real (dp)                            :: CoefV1, CoefV2, CoefV3, CoefV4, c1, c2, &
    CoefA1, CoefA2, CoefA3, CoefA4

    Xi        = Sqrt( 1 + 16 * self%m2 + 32 * self%m4 - 8 * c * (1 + 2 * self%m2)**2 )
    c2        = c - 2 * self%m2 ; c1 = 1 + c
    Xi2       = Sqrt( (1 + Xi)**2 - 16 * c2**2 )
    n2        = 2 * Xi/(1 + 4 * c - 8 * self%m2 + Xi);    CoefV1 = 0
    h2        = 2 * Xi/(1 - 4 * c + 8 * self%m2 + Xi);    k2 = n2 * h2/Xi

    CoefV2 = 0; CoefV3 = 0; CoefV4 = 0; CoefA1 = 0; CoefA2 = 0; CoefA3 = 0; CoefA4 = 0

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      CoefV1 = - 4 * (3 * c**2 * (1 + 2 * c) + 2 * c * ( c * (8 * c - 3) - 4 ) * self%m2 &
            + 2 * (  2 + c * ( 2 * c * (1 + 7 * c) - 7 )  ) * self%m4 + 4 * self%m6    * &
            (5 - 11 * c + 6 * c**3) + 8 * (8 + c * (12 + c * (9 + 2 * c))) * self%m8   + &
            32 * (c**2 - 2) * self%m10) * (1 - 4 * c + 8 * self%m2 + Xi) / c1    / &
            c2**2 / (1 + 4 * c - 8 * self%m2 - Xi) / Xi2

      CoefV2 = 8 * (    c**2 - c**4 - 2 * c * (2 - 2 * c + 3 * c**3) * self%m2    +      &
            2 * (   2 + c * (  c * ( 7 + (2 - 7 * c) * c ) - 8  )   ) * self%m4  -       &
            4 * (   c * (  8 + c * ( c * (3 * c - 4) - 5 )  ) - 4  ) * self%m6   -       &
            8 * (   c * (  1 + c * ( 5 + c * (4 + c) )  ) - 3   ) * self%m8      -       &
            8 * (  2 + c * ( c * (3 * c - 7) - 15 )  ) * self%m10 - 8 * (  13 + c *      &
            ( 12 + c * (5 + 2 * c) )  ) * self%m6**2 - 32 * (c**2 - 2) * self%m2**7    ) &
            /c1 / c2**2 / ( c + 2 * self%m2 * (self%m2 - 1) ) / Xi2

      CoefV3 = - 4 * c1 * self%m4 * (1 - 4 * c + 8 * self%m2 - Xi) *  (c**2 + 10 * c**2 * &
            self%m**2 + 4 * ( c * (3 * c - 5) - 1) * self%m4 + 4 * (2 - 3 * c) * self%m6  &
             + 4 * (1 + 2 * c) * self%m8 + 16 * self%m10) / c2**3  / Xi2 &
            / ( c + 2 * self%m2 * (self%m2 - 1) )

      CoefV4 = (   2 * (  2 * c * (  4 + c * ( 2 * c * (4 + c) - 3 )  ) * self%m2 + 2 * ( c &
      * (  21 + c * ( c * (13 + 15 * c) - 23 )  ) - 4   ) * self%m4 + 4 * (   c * ( 34 +    &
      (c - 2) * c * (8 + 9 * c) ) - 11  ) * self%m6 + 4 * (  c * ( 49 + (c - 3) * c    *    &
      (2 * c - 7) ) - 27  ) * self%m8 + 16 * (  c * ( 3 + c * (7 + c) ) - 7  ) * self%m10 - &
      c**2 * (2 + c + 2 * c**2)  ) * (1 + 4 * c - 8 * self%m2 - Xi)   ) / c1**2 / c2**3 / Xi2

    end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      CoefA1 = -4 * (  2 * c * self%m2 * (4 - 33 * self%m2 + 20 * self%m4 + 172 * self%m**6 &
             - 80 * self%m8) + 2 * c**3 * (8 * self%m2 - 3 + 34 * self%m4 + 16 * self%m6 + &
            16 * self%m8) + 4 * self%m4 * (9 * self%m2 - 1 - 10 * self%m4 - 40 * self%m6 + &
            48 * self%m8) + c**2 * (32 * self%m2 - 3 - 48 * self%m4 - 208 * self%m6      + &
            152 * self%m8 + 80 * self%m10)  ) * (1 - 4 * c + 8 * self%m2 + Xi) / c1      / &
            c2**2 / (Xi - 1 - 4 * c + 8 * self%m**2)/ Xi2

      CoefA2 = 8 * ( 4 * c**3 * self%m4 * (24 * self%m4 - 7 - 38 * self%m**2 + 10 * self%m6  &
      + 8 * self%m8) + c**4 * (4 * self%m2 - 1 + 34 * self%m4 + 16 * self%m6 + 16 * self%m8) &
      + 4 * self%m4 * (1 - 2 * self%m2 - 6 * self%m4 + 16 * self%m6 + 30 * self%m8 - 84   *  &
      self%m10 + 48 * self%m6**2) - 4 * c * (self%m2 - 2 * self%m4 - 4 * self%m6 + 26     *  &
      self%m8 + 78 * self%m10 - 148 * self%m6**2 + 40 * self%m2**7) + c**2 * (1 - 2       *  &
      self%m2 + 2 * self%m4 + 72 * self%m6 + 288 * self%m8 - 408 * self%m10 + 72 * self%m6**2&
      + 80 * self%m2**7) ) / c2**2 / (c - 2 * self%m2 - 2 * c * self%m2       &
      + 2 * self%m4 + 2 * c * self%m4 + c**2) / Xi2

      CoefA3 = - 4 * c1 * self%m4 * ( 4 * c * self%m4 * (4 * self%m4 - 3 - 7 * self%m2) +    &
      c**2 * (4 * self%m2 + 24 * self%m4 - 1) + 4 * (self%m4 - 5 * self%m8 + 10 * self%m10) &
      ) * (4 * c - 8 * self%m2 + Xi - 1) / c2**3 / (c + 2 * self%m2 * (self%m2 - 1) ) / Xi2

      CoefA4 = 2 * ( 2 * c**4 * (1 + 17 * self%m4 + 36 * self%m6 + 8 * self%m**8)    + &
          4 * self%m4 * (2 - self%m2 - 7 * self%m4 + 26 * self%m6 + 80 * self%m**8) - &
          2 * c * self%m2 * (4 - 3 * self%m2 + 30 * self%m4 + 138 * self%m6 + 244   * &
          self%m8 + 32 * self%m10) + c**3 * (1 - 30 * self%m2 - 42 * self%m4 - 132  * &
          self%m6 - 172 * self%m8 + 40 * self%m10) + c**2 * (2 - 6 * self%m2 + 82   * &
          self%m4 + 220 * self%m6 + 444 * self%m8 + 344 * self%m10) )               * &
          (Xi - 1 - 4 * c + 8 * self%m2) / c1**2 / c2**3 / Xi2

    end if

    Coef(2) = self%EWAdd(CoefV1, CoefA1); Coef(1) = self%EWAdd(CoefV2, CoefA2)
    Coef(3) = self%EWAdd(CoefV3, CoefA3); Coef(4) = self%EWAdd(CoefV4, CoefA4)

  end subroutine CoefsCparam

!ccccccccccccccc

  real (dp) function FunCp(self, Coef, n2, h2, k2, Xi, Xi2, c, theta)
    class (MassiveNS)      , intent(in) :: self
    real (dp)              , intent(in) :: c, theta, n2, h2, k2, Xi, Xi2
    real (dp), dimension(4), intent(in) :: Coef
    real (dp), dimension(4)             :: res
    real (dp)                           :: IndTerm, Stheta, Ctheta, Xi2Ctheta, c1, &
    c2, IndTermV, IndTermA

    c2        = c - 2 * self%m2 ; c1 = 1 + c  ; Ctheta = Cos(2 * theta)
    Xi2Ctheta = Xi2**2 - 2 * (1 - Ctheta) * Xi; Stheta = Sin(2 * theta)

    IndTermV = 0 ; IndTermA = 0

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      IndTermV = Stheta * Xi * (    4 + (   2 * (  8 * self%m6 - 4 * self%m2 - 2 * self%m4  &
        + c * (3 - 2 * self%m2 - 12 * self%m4) + c**2 * (5 + 14 * self%m2 + 14 * self%m4 +  &
        8 * self%m6)  ) * Xi2Ctheta   ) / c2**2 / (1 + 4 * c - 8 * self%m2 + Ctheta * Xi) - &
        (  (1 + 2 * self%m2) * ( 3 + 32 * c**2 - 52 * self%m**2 + 112 * self%m4 + 2 * c   * &
        (13 - 56 * self%m2 + 8 * self%m4) + 3 * Ctheta * (1 + 2 * c - 4 * self%m2) * Xi ) * &
        Xi2Ctheta  ) / c2 / (1 + 4 * c - 8 * self%m2 + Ctheta * Xi)**2 - ( 128 * c1**4    * &
        self%m6 * (1 + 2 * self%m2) * (1 - 4 * c + 8 * self%m2 + Xi) * Xi2Ctheta ) / c2**2  &
        / (8 * self%m2 - 1 - 4 * c + Xi) / (1 - 4 * c + 8 * self%m2 + Ctheta * Xi) / Xi2**2 &
            ) / 4 / c1**2 / Sqrt(Xi2Ctheta)

    end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      IndTermA = STheta * Xi * (  4 * Xi2**2 - (  2 * (1 + 2 * self%m2) * ( 4 * self%m2   - &
      14 * self%m4 + c * (18 * self%m2 + 4 * self%m4 - 3) + c**2 * (6 * self%m2 - 5       + &
      10 * self%m4) ) * (1 + 4 * c - 8 * self%m2 + Xi) * (1 - 4 * c + 8 * self%m2 + Xi)   * &
      (2 * Xi * (CTheta - 1) + Xi2**2)  ) / ( c2**2 * (1 + 4 * c - 8 * self%m2 + CTheta   * &
      Xi) ) + (128 * c1**4 * self%m6 * (4 * self%m2 - 1) * (1 - 4 * c + 8 * self%m2 + Xi) * &
      Xi2Ctheta) / (  (Xi - 1 - 4 * c + 8 * self%m2) * (1 - 4 * c + 8 * self%m2 + CTheta  * &
      Xi) * c2**2 ) + (  (1 + 2 * self%m2)**2 * (8 * self%m2 - 1 - 4 * c - Xi) * (1 - 4 * c &
      + 8 * self%m2 + Xi) * (3 + 32 * c**2 + 112 * self%m4 + 3 * CTheta * Xi - 4 * self%m2  &
      * (13 + 3 * CTheta * Xi) + 2 * c * (13 - 56 * self%m2 + 8 * self%m4 + 3 * CTheta * Xi)&
      ) * Xi2Ctheta  )/ ( c2 * (1 + 4 * c - 8 * self%m2 + CTheta * Xi)**2 )  )/ 4 / c1**2 / &
      Xi2**2 / Sqrt(Xi2Ctheta)

    end if

    IndTerm = self%EWAdd(IndTermV, IndTermA)

    call f90Elliptic4(h2, n2, theta, k2, res)

    FunCp = 2 * ( IndTerm + dot_product(coef, res) )/c1

  end function FunCp

!ccccccccccccccc

  real (dp) function f1Q(self, t)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: t
    real (dp)                     :: t2, t1, f1VQ, f1AQ, lg

    t2 = Sqrt( 1 - 4 * self%m2/t ); t1 = 1 - t; lg = log( (1 - t2)/(1 + t2) )
    f1VQ = 0; f1AQ = 0

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      f1VQ = - 2/t1 * ( t2 * (1 + 4 * self%m2 * t + t**2 ) + ( 1 - 8 * self%m4 - &
              4 * self%m2 * t1 + t**2) * lg  )

    end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      f1AQ = - 2/t1 * ( t2 * (1 - 8 * self%m2 * t  + t**2) + ( 1 + 16 * self%m4 + t**2 - &
               2 * self%m2 * (1 + 6 * t - t**2) ) * lg  )

    end if

    f1Q = self%EWAdd(f1VQ, f1AQ)

  end function f1Q

!ccccccccccccccc

  real (dp) function f2Q(self, t)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: t
    real (dp)                     :: t1, t2, t11, t3, A, f2VQ, f2AQ, lg

    f2VQ = 0; f2AQ = 0;

    shape_if: if ( self%shape(:6) == 'thrust' ) then

      t2 = 1 - t;  t1 = Sqrt( t2**2 + 4 * self%m2 );  t11 = 1 - t1
      t3 = t2 - 2 * self%m2 + t1; A = t2/t1/t11; lg = log( t3/(t11 + self%m2)/2 )

      if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        f2VQ = A * ( 4 * t11 * (3 + 8 * self%m2 + t1) + t11**2 * t3**2/(t11 + self%m2)**2 + &
        32 * self%m2 * (1 + 2 * self%m2) * (t11 + self%m2)/t3 - 8 * (1 + 2 * self%m2) * t3  &
        + 8 * ( 1 - 8 * self%m4 - 4 * self%m2 * t11 + t1**2 ) * lg )/8

      end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        f2AQ = A * (  4 * t11 * (3 + t1 - 2 * self%m2 * (13 - 5 * t1) ) + (1 + 2 * self%m2)  * &
        t11**2 * t3**2/(t11 + self%m2)**2 + 32 * self%m2 * (1 - 4 * self%m2) * (t11 + self%m2) &
        /t3 - 8 * t3 * (t11 - 4 * self%m4 - self%m2 * (5 - 2 * (4 - t1) * t1) )/(t11 + self%m2)&
        + 8 * (  1 + 16 * self%m4 + t1**2 - 2 * self%m2 * ( 1 + (6 - t1) * t1 )  ) * lg )/8

      end if

    else if ( self%shape(2:3) == 'JM' ) then

      t1 = Sqrt( (1 - t + self%m2)**2 - 4 * self%m2 );  t2 = 1 + t1
      t11 = t - self%m2; lg = log( (t2 - self%m2 - t)/t/2 )

      if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        f2VQ = (  self%m8/t**2 - 2 * self%m6 * t2/t**2 + (t2 - 3 * t) * (t + t1 - 7) -  &
        self%m4 * ( 22 + 64 * t/(self%m2 + t - t2) - 2 * t2/t - t2**2/t**2 ) + self%m2 * &
        ( 8 * t * ( 7 - 4/(self%m2 + t - t2) ) - 2 * t2**2/t - 2 * (11 + 7 * t1) )  )/t11/8 &
        - ( 2 + 5 * self%m2 - t - (2 - 8 * self%m4)/t11 ) * lg

       end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        f2AQ = (t2 - self%m2 - 3 * t)/t**2/t11/(t2 - t - self%m2) * (   2 * self%m10 -      &
        t**2 * (t - t2) * (t + t1 - 7) - self%m8 * (3 + 20 * t + 4 * t1) + 2 * self%m6 *   &
        ( t * (11 + 8 * t) + t1 + 12 * t * t1 + t1**2 ) + self%m4 * ( 20 * t**3 - 4 * t * &
        t1 * t2 + t2**2 - 4 * t**2 * (1 + 9 * t1) ) - 2 * self%m2 * t * (  9 * t**3 + t**2 &
        * (7 - 8 * t1) + t2**2 - t * ( 12 + t1 * (17 + t1) )  )   )/8 + ( 2 * self%m2 * t    &
        - 2 + 7 * self%m2 - 2 * self%m4 + t + 2 * (1 - 6 * self%m2 + 8 * self%m4)/t11 ) * lg

      end if

    end if shape_if

    f2Q = self%EWAdd(f2VQ, f2AQ)

  end function f2Q

!ccccccccccccccc

  real (dp) function f3Q(self, t)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: t
    real (dp)                     :: t1, t2, t3, f3VQ, f3AQ, lg, A, t11

    f3VQ = 0; f3AQ = 0

    shape_if: if ( self%shape(:6) == 'thrust' ) then

      t2 = 1 - t; t1 = Sqrt( t2**2 + 4 * self%m2 ); t11 = 1 - t1; t3 = t - t1
      lg = t2/t1/t11 * log( -t11/t3 ); A = t2 * (1 + t - 2 * t1)/t1/t3/t11**2/2

      if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        f3VQ = A * ( t11 * t3 * (3 + t) - 4 * self%m2 * (1 + 2 * t1 * (t11 + t) - 3 * t) - &
        8 * self%m4 * t2 ) - (1 - 8 * self%m4 - 4 * self%m2 * t11 + t1**2) * lg

      end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        f3AQ = A * (  16 * self%m4 * t2 + t11 * t3 * (3 + t) - 2 * self%m2 * ( 2        - &
        4 * t1**3 + (11 - t) * t + t1**2 * (17 + 3 * t) - t1 * (13 + (16 - t) * t) )  ) - &
        (  1 + 16 * self%m4 + t1**2 - 2 * self%m2 * ( 1 + (6 - t1) * t1 )  ) * lg

      end if

    else if ( self%shape(2:3) == 'JM' ) then

      t1 = Sqrt( (1 - t + self%m2)**2 - 4 * self%m2 );  t11 = t - self%m2;  lg = log(t1/t11 - 1)

      if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        f3VQ = (  (t1 - 2 * t11) * ( 4 * (1 + 2 * self%m2) * t11**2 - ( self%m4 +  &
        6 * self%m2 * t + t * (4 + t) ) * t1 + t11 * t1**2 )  )/ ( 2 * t11**2 * &
        (t1 - t11) ) + ( t - 2 - 5 * self%m2 + (2 - 8 * self%m4)/t11 ) * lg

      end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        f3AQ = (t1 - 2 * t11)/t11**2/(t1 - t11)/2 * (  4 * t11**2 * ( 1 + 2 * self%m4 -  &
        2 * self%m2 * (2 + t) ) + ( 6 * self%m6 + 6 * self%m2 * t * (3 + t) - t *    &
        (4 + t) - self%m4 * (1 + 12 * t) ) * t1 + (1 + 2 * self%m2) * t11 * t1**2  ) + &
        (2 * self%m2 * t - 2 + 7 * self%m2 - 2 * self%m4 + t + ( 2 * (1 - 6 * self%m2  &
        + 8 * self%m4) )/t11) * lg

      end if

    end if shape_if

    f3Q = self%EWAdd(f3VQ, f3AQ)

  end function f3Q

!ccccccccccccccc

  real (dp) function f4Q(self, t)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: t
    real (dp)                     :: t1, t11, t2, t3, f4AQ, f4VQ, lg, A

    t3 = 1 - t  ; t1 = Sqrt( t3**2 + 4 * self%m2 );  t2 = t1 - t; f4VQ = 0; f4AQ = 0
    t11 = 1 - t1; lg = log(t11/t2)/t3; A = 2 * (1 - 2 * t1 + t)/t11/t2

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      f4VQ = A * (2 * self%m2 + 4 * self%m4 + t1 - t1**2 - t + t1 * t) &
          - 2 * (1 - 4 * self%m2 - 8 * self%m4 + 4 * self%m2 * t + t**2) * lg

    end if; if( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      f4AQ = A * ( t11 * t2 + 2 * self%m2 - 8 * self%m4 ) - 2 * ( 1 + 16 * self%m4 + &
             t**2 - 2 * self%m2 * (1 + 6 * t - t**2) ) * lg

    end if

    f4Q = self%EWAdd(f4VQ, f4AQ)

  end function f4Q

!ccccccccccccccc

  real (dp) function A0(self)
    class (MassiveNS), intent(in) :: self
    real (dp)                     :: AV, AA

    A0 = 0; AV = 0; AA = 0

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      AV = - 2 * self%m2 * (self%thrustMin + 4 * self%m2)/(1 + self%v)

    end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      AA = - 4 * self%m2 * (2 - 4 * self%m2 + self%v)/(1 + self%v)

    end if

    A0 = self%EWAdd(AV, AA)

  end function A0

!ccccccccccccccc

  real (dp) function B1NS(self)
    class (MassiveNS), intent(in) :: self
    real (dp)                     :: t2, lg, ln, B1A, B1V

    t2 = 1 + self%v; ln = log(t2/2); lg = ln - self%lm; B1A = 0; B1V = 0

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      B1V = 4 * self%m2/t2 - 2 * self%m2 * self%v + 2 * ln - 8 * self%m4 * lg

    end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      B1A = 4 * self%m2/t2 + 4 * self%m2 * self%v - 4 * self%m2 * (3 - 4 * self%m2) * lg &
            + 2 * ln

    end if

    B1NS = 8 * self%EWAdd(B1V, B1A)/3

  end

!ccccccccccccccc

  real (dp) function A0MS(self)
    class (MassiveNS), intent(in) :: self
    real (dp)                     :: A, V

    V = 0; A = 0

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

      if (self%m >  1e-3_dp) V = - 24 * self%m4 / self%v
      if (self%m <= 1e-3_dp) V = ( - 24 * self%m2 - 48 * self%m4 - 144 * self%m6 - &
                                    480 * self%m8 - 1680 * self%m10 ) * self%m

    end if; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

      if (self%m >  1e-3_dp) A = - 12 * self%m2 * self%v
      if (self%m <= 1e-3_dp) A = ( - 12 + 24 * self%m2 + 24 * self%m4 + 48 * self%m6 &
                                  + 120 * self%m8 + 336 * self%m10 )  * self%m

    end if

    A0MS = self%EWAdd(V, A)

  end function A0MS

!ccccccccccccccc

  real (dp) function A1loop(self)
    class (MassiveNS), intent(in) :: self
    real (dp)                     :: h, lh, A, V, dilog

    A = 0; V = 0

    massIf: if ( self%m > 0.499_dp ) then

      h = sqrt(0.5_dp - self%m); lh = 0; if (h > 1.d-6) lh = log(h)

      current_if_1: if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        shape_if_1: if ( self%shape(:6) == 'thrust' .or. self%shape(:6) == 'Cparam' ) then

          V =  42.595415652878366_dp * h**2 - 241.2319251819811_dp * h**4 + 4.954944900637917_dp   + &
          118.01211476041078_dp * h**6 - 289.658731751325_dp * h**8 - 311.89178067391487_dp * h**10 - &
          553.6278884247788_dp * h**12 + h**9 * (269.863816214773_dp + 42.74285714285714_dp * lh)  + &
          h**11 * (437.7395234630508_dp + 124.06435786435786 * lh) + h**3 * (236.83518400168563_dp &
          + 128 * lh) - h**5 * (129.3500602687456_dp + 473.6_dp/3 * lh)  - 48 * h    + &
          h**7 * (216.09698685182715_dp + 216.22857142857143_dp * lh)

        else if ( self%shape(2:3) == 'JM' ) then

          V = 4.954944900637917_dp + 42.595415652878366_dp * h**2 - 553.6278884247788_dp * h**12 + &
          148.11234489001265_dp * h**3 - 241.2319251819811_dp * h**4 + 44.07477463565101_dp * h**5 &
          + 118.01211476041078_dp * h**6 + 19.2854288762748_dp * h**7 - 289.658731751325_dp * h**8 &
          + 330.21767767807745_dp * h**9 - 311.89178067391487_dp * h**10 - 48 * h                  &
          + 416.86212391885135_dp * h**11

        end if shape_if_1

      end if current_if_1; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        shape_if_2: if ( self%shape(:6) == 'thrust' .or. self%shape(:6) == 'Cparam' ) then

          A = - 24.65386830263016_dp + 42.595415652878366_dp * h**2 + 36.88227950096423_dp * h**8  &
         + 136.5090860703081_dp * h**4 - 552.1365744108343 * h**6 - 64 * h**3                &
         + h**11 * (416.3450777973517_dp + 124.13968253968254_dp * lh)                          &
         - 567.4437614406518_dp * h**12 + h**7 * (129.9422250626211_dp - 3072 * lh/10)            &
         + h**9 * (148.53285087521925_dp + 132.8761904761905_dp * lh) +                         &
         h**5 * (364.89382400449506_dp + 1024 * lh/3) - 321.15654257867675_dp * h**10

        else if ( self%shape(2:3) == 'JM' ) then

          A = -24.65386830263016_dp + 42.595415652878366_dp * h**2 - 64 * h**3                  + &
          136.5090860703081_dp * h**4 + 128.29958637336705_dp * h**5 - 552.1365744108342_dp * h**6 + &
          513.543705597303_dp * h**7 + 36.88227950096423_dp * h**8 - 11.836572583565095_dp * h**9  - &
          321.15654257867686_dp * h**10 + 376.8249909763481_dp * h**11 - 567.4437614406518_dp * h**12

        end if shape_if_2
      end if

    else if ( self%m > 1e-3_dp ) then

      current_if_2: if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        shape_if_3: if ( self%shape(:6) == 'thrust' .or. self%shape(:6) == 'Cparam' ) then

          V = - 4 * (5 - 8 * self%m2 - 24 * self%m4) * log( self%thrustMin/2/self%m )   + &
          2 * (1 - 4 * self%m4) * (  - 2 * self%v * self%thrustMin/(1 - 2 * self%m2)    + &
          Pi2 + 2 * l2 * log(8 * self%m4) - 4 * log(8 * self%m2) * log(self%thrustMin)  + &
          (   4 * self%v  * log(self%m * self%thrustMin/2/self%v)  )/(1 - 2 * self%m2)  + &
          6 * log(self%thrustMin)**2 + 8 * log( self%thrustMin/2/self%m )               * &
          log( self%v * self%thrustMin/2/self%m**2 ) - 4 * Dilog( self%thrustMin/2 )    + &
          8 * Dilog( self%thrustMin/(1 + self%v) )  )

        else if ( self%shape(2:3) == 'JM' ) then

          V = - 4 * (5 - 8 * self%m2 - 24 * self%m4) * log( self%thrustMin/self%m/2 ) + 2 *  &
          (1 - 4 * self%m4) * (  - 2 * self%v * self%thrustMin/(1 - 2 * self%m2) + Pi2 + 6 * &
          log(self%thrustMin)**2 + (  4 * self%v * log(self%m * self%thrustMin/self%v/2)  )/ &
          (1 - 2 * self%m2) - 4 * log(self%thrustMin) * (3 * l2 + 2 * self%lm) + 2 * l2    * &
          (3 * l2 + 4 * self%lm) +  8 * log( self%v * self%thrustMin/self%m2/2 )           * &
          log( self%thrustMin/2/self%m ) - 4 * Dilog( self%thrustMin/2 )                   + &
          8 * Dilog( self%thrustMin/(1 + self%v) )  ) + 8 * (1 + 2 * self%m2) * log(self%v)  &
          * ( 2 * (1 - 2 * self%m2) * log(self%thrustMin/2/self%m) + self%v)

        end if shape_if_3

      V = V - 2 * Pi2 - 4 * self%lm - 16 * self%lm**2

      end if current_if_2; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        shape_if_4: if ( self%shape(:6) == 'thrust' .or. self%shape(:6) == 'Cparam' ) then

          A = 4 * self%v**2 * (   self%v * (  self%v - 2 * log( self%v*(1 + self%v)/2     / &
          self%m**3 ) - 1  ) - (5 - 8 * self%m2) * log( self%thrustMin/2/self%m )         + &
          2 * (1 - 2 * self%m2) * (2 * Dilog( self%thrustMin/(1 + self%v) ) + 7           * &
          log( self%thrustMin/2 )**2/2 + 2 * log(self%v) * log( self%thrustMin/2/self%m ) - &
          8 * self%lm * log( self%thrustMin/2/Sqrt(self%m) ) + Pi2/4                      - &
          Dilog( self%thrustMin/2 )  )   )

        else if ( self%shape(2:3) == 'JM' ) then

          A = 4 * self%v**2 * (   self%v * (  self%v - 2 * log( self%v * (1 + self%v)/2     / &
          self%m**3 ) - 1  ) - (5 - 8 * self%m2) * log( self%thrustMin/2/self%m )           + &
          2 * (1 - 2 * self%m2) * ( - Dilog(self%thrustMin/2)  + 2 * Dilog( self%thrustMin  / &
          (1 + self%v) ) + 3.5 * log( self%thrustMin/2 )**2 + Pi2/4 + 2 * log(self%v)       * &
          log( self%thrustMin/2/self%m ) - 8 * self%lm * log( self%thrustMin/2/Sqrt(self%m) ) &
            )   ) + 8 * self%v**2 * log(self%v) * (  self%v + 2 * (1 - 2 * self%m2)         * &
          log( self%thrustMin/2/self%m )  )

        end if shape_if_4

        A = A - 2 * Pi2 - 4 * self%lm - 16 * self%lm**2

      end if

    else

      current_if_3: if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then

        shape_if_5: if ( self%shape(:6) == 'thrust' .or. self%shape(:6) == 'Cparam' ) then

          V = (4 + 48 * self%lm) * self%m2 + (848/9._dp - 800 * self%lm/3) * self%m6 &
          + (89._dp/2 - 556 * self%lm) * self%m8 + self%m4 * (60 - 40 * self%lm      &
          - 64 * self%lm**2 - 8 * Pi2)

        else if ( self%shape(2:3) == 'JM' ) then

          V = 4 * (4 * self%lm - 3) * self%m2 - 32 * (2 + 87 * self%lm)*self%m6/9     - &
          (821/6._dp + 812 * self%lm) * self%m8 - 4 * self%m4 * (1 + 2 * self%lm * (13 + &
          8 * self%lm) + 2 * Pi2)

        end if shape_if_5

      end if current_if_3; if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then

        shape_if_6: if ( self%shape(:6) == 'thrust' .or. self%shape(:6) == 'Cparam' ) then

          A = (- 628/9._dp + 496 * self%lm/3) * self%m6 + (180 * self%lm - 37/6._dp)  * &
          self%m8 + self%m2    * (4 - 16 * self%lm - 96 * self%lm**2 - 12 * Pi2)     + &
          self%m4 * ( 128 * self%lm**2 + 16 * Pi2 - 4 - 72 * self%lm )

        else if ( self%shape(2:3) == 'JM' ) then

          A = 4 * (47 + 276 * self%lm) * self%m6/9 + (36.5_dp + 180 * self%lm) * self%m8  - &
          12 * self%m2 * (1 + 4 * self%lm + 8 * self%lm**2 + Pi2) + 4 * self%m4 * (7 + &
          2 * self%lm * (7 + 16 * self%lm) + 4 * Pi2)

        end if shape_if_6
      end if
    end if massif

    A1loop = self%EWAdd(V, A)/3 + self%A1NSDiff()

  end

!ccccccccccccccc

  real (dp) function A1NSDiff(self)
    class (MassiveNS), intent(in) :: self
    real (dp)                     :: h2, h3, factV, factA, dilog

    A1NSDiff = 0; if ( self%shape(:6) /= 'Cparam' ) return

    if ( self%current(:6) == 'vector' .or. self%current(:3) == 'all' ) then
      factV = 1 + 2 * self%m2
    else if ( self%current(:5) == 'axial' .or. self%current(:3) == 'all' ) then
      factA = 1 - 4 * self%m2
    end if

    if ( self%m > 1e-3_dp ) then

      h2 = Sqrt(1 - 8 * self%m4);  h3 = 1 - 2 * self%m2

      A1NSdiff = 1 - 4 * self%m2 - 2 * (1 - 2 * self%m4)/self%m2 * &
      log( (1 - self%v)/(2 * self%m) ) + self%v * (  1 + 2 * &
      log( (1 - self%v)/(2 * self%m2 * h3 * self%v) )  ) +   &
      2 * h2/self%m2 * Log(  (h2 - self%v)/( 2 * self%m * Sqrt(h3) )  ) + &
      h3 * ( -Pi2/6 + 13 * l2**2 + 36 * l2 * self%lm + 24 * self%lm**2  + &
      4 * Log(2 * self%m) * Log(self%v) + 2 * ( 5 * l2 + 6 * self%lm ) * Log(h3) - &
      2 * Log(2 * self%m4 * self%v**2 * h3**2) * Log(1 - self%v) + Log(1 - self%v)**2 -&
      4 * Log(4 * self%m**2 * h3) * Log(1 - h2) + (- 12 * l2 + 8 * &
      Log( (1 - h2)/self%m2) ) * Log(h2 -self%v) + 2 * Dilog( (1 + self%v)/2 ) - &
      2 * Dilog( (self%v - h2)/(1 - h2) ) + 2 * Dilog( (h2 - self%v)/(1 + h2) ) + &
      2 * Dilog( (self%v + h2)/(h2 - 1) ) - 2 * Dilog( (self%v + h2)/(1 + h2) ) )

    else

      A1NSdiff = (119._dp/6 + 16 * self%lm) * self%m4 - (29._dp/9 + 112 * self%lm/3) &
      * self%m6 + (10471._dp/360 + 272 * self%lm/3) * self%m8 - (2569._dp/300 + &
      32 * self%lm/5) * self%m10 - Pi2/6 + self%m2 * (4 + 4 * self%lm + Pi2/3)

    end if

    A1NSdiff = - 4 * (self%EWAdd(factV, factA) * A1NSdiff + Pi2/6)/3

  end function

!ccccccccccccccc

  real (dp) function Hcorr(self)
    class (MassiveNS), intent(in) :: self

    Hcorr = self%Dirac(2) - self%Dirac(1) * self%A1Singular &
      + 2 * ( self%B1 - self%Dirac(1) * self%B1Singular ) * self%lm

  end function Hcorr

!ccccccccccccccc

  function JetMassExp(self, hard, xiJ, xiB) result(res)
    class (MassiveNS)   , intent(in) :: self
    character (len = *) , intent(in) :: hard
    real  (dp)          , intent(in) :: xiJ
    real  (dp), optional, intent(in) :: xiB
    real (dp),      dimension(0:4,3) :: res
    real (dp)                        :: A0, B0, lgJ
    integer                          :: hardFac

    res = 0; if ( self%singular(:3) /= 'abs' .and. present(xiB) ) return
    A0 = self%MassVar('A0'); B0 = 0; hardFac = 0; if ( hard(:6) == 'expand' ) hardFac = 1

    if ( self%singular(:3) == 'abs' ) B0 = self%B1

    if ( .not. present(xiB) ) then

      if ( self%ES(:1) == 'Q' ) res(0,1) = A0 * self%alphaJ * (1 + Pi2)/3
      if ( self%ES(:1) == 'P' ) res(0,1) = A0 * self%alphaJ * (5 * Pi2/3 - 1)

      if ( self%singular(:3) == 'abs' ) res(0,1) = res(0,1) + self%deltaM(1) * self%AMS/2 &
       +  ( self%alphaJ * xiJ + hardFac * self%alphaH * (1 - xiJ) )/2 * self%Hcorr()

      res(1:2,1) = self%alphaJ * [ ( A0/3 - B0/2 + (A0 - 1) * self%B1Singular/2 ), 2 * A0/3 ]

    else

      res(:1,1) = [ (  xiJ * ( self%alphaJ * xiB + hardFac * self%alphaM * (1 - xiB) ) + &
      hardFac * (1 - xiJ) * self%alphaH  ) * self%Hcorr() + self%deltaM(1) * self%AMS, &
      self%alphaJ * ( (A0 - 1) * self%B1Singular - B0 ) ]/2

    end if

    lgJ = log( self%mass**2/self%Q/self%matEl%scales('muS') )

    call self%MatEl%AddLogs(res(:2,:1), lgJ)

  end function JetMassExp

!ccccccccccccccc

  recursive real (dp) function MassVar(self, str) result(res)
    class (MassiveNS)   , intent(in) :: self
    character (len = *) , intent(in) :: str

    if      ( str(:5) == 'mPole' ) then
      res = self%mPole
    else if ( str(:4) == 'mass' ) then
      res = self%mass
    else if ( str(:2) == 'mm' ) then
      res = self%mm
    else if ( str(:2) == 'm2'   ) then
      res = self%m2
    else if ( str(:2) == 'lm'   ) then
      res = self%lm
    else if ( str(:1) == 'm'    ) then
      res = self%m
    else if ( str(:4) == 'tmin' ) then
      res = self%tmin
    else if ( str(:5) == 'width' ) then
      res = self%width
    else if ( str(:6) == 'alphaH' ) then
      res = self%alphaH
    else if ( str(:6) == 'alphaM' ) then
      res = self%alphaM
    else if ( str(:5) == 'Hcorr' ) then
      res = 0
      if (self%singular(:3) == 'abs' ) res = self%Hcorr()/self%MassVar('A0')
    else if ( str(:2) == 'A0' ) then
      res = 1
      if ( self%singular(:3) == 'abs' ) res = 1 + self%Dirac(1)
    end if

  end function MassVar

!ccccccccccccccc

  integer function numFlav(self)
    class (MassiveNS), intent(in) :: self
    numFlav = self%nf
  end function numFlav

!ccccccccccccccc

  character (len = 6) function EShape(self, str)
    class (MassiveNS)  , intent(in) :: self
    character (len = *), intent(in) :: str

    if ( str(:5) == 'shape' ) then
      EShape = self%shape
    else if ( str(:6) == 'EShape' ) then
      EShape = self%ES
    end if

  end function EShape

!ccccccccccccccc

  function MassDelta(self) result(delta)
    class (MassiveNS), intent(in) :: self
    real (dp)      , dimension(3) :: delta
    delta = self%deltaM
  end function MassDelta

!ccccccccccccccc

  function MassDeltaShift(self) result(delta)
    class (MassiveNS), intent(in) :: self
    real (dp)      , dimension(3) :: delta
    delta = self%deltaM2
  end function MassDeltaShift

!ccccccccccccccc

  type (MatricesElementsMass) function matElementScales(self)
    class (MassiveScales), intent(in) :: self

     select type (or => self%MatEl)
     type is (MatricesElementsMass)
       matElementScales = or
     end select

  end function matElementScales

!ccccccccccccccc

  type (MatrixElementsMass) function matElementNS(self)
    class (MassiveNS), intent(in) :: self

     select type (or => self%MatEl)
     type is (MatrixElementsMass)
       matElementNS = or
     end select

  end function matElementNS

!ccccccccccccccc

  real (dp) function EWAdd(self, V, A)
    class (MassiveNS), intent(in) :: self
    real (dp)        , intent(in) :: V, A

    EWAdd = 0

    if ( self%current(:3) == 'all'         ) then
      EWAdd = dot_product( self%EW, [V, A] )
    else if ( self%current(:6) == 'vector' ) then
      EWAdd = V
    else if ( self%current(:5) == 'axial'  ) then
      EWAdd = A
    end if

  end function EWAdd

!ccccccccccccccc

  real (dp) function ThrustFO1loop(t)
    real (dp), intent(in) :: t
    real (dp)             :: t2

    t2 = t**2; ThrustFO1loop = 0;  if (t < 0 .or. 3 * t > 1) return

    ThrustFO1loop = (  4/t/3 * (3 - 9 * t + 9 * t**3 + log(1/t - 2) * &
                    ( 6 * t - 4 - 6 * t2) - 3 * t2 )  )/(t - 1)/2

  end function ThrustFO1loop

!ccccccccccccccc

  real (dp) function CParamFO1Loop(c2)
    real (dp), intent(in)   :: c2
    real (dp)               :: c, c1, root1, root2, k, n
    real (dp), dimension(3) :: res, coef

    c = c2/6;  c1 = 1 + c;  CparamFO1loop = 0;  if (c < 0 .or. c > 0.125 ) return

    root1 = sqrt(1 - 8 * c);  root2 = Sqrt( 1 + root1 - 4 * c * (1 + 2 * c) )

    coef = [  8 * sr2 * (1 - c) / c / c1 / root2                            ,  &
            - 3 * (1 + 2 * c) * root2 / sr2 / c / c1**3                     ,  &
            - 4 * sr2 * (2 + c + 2 * c**2) * (root1 - 1 + 4 * c + 8 * c**2) /  &
             c / c1**3 / (root1 - 1 + 4 * c) / root2 ]

    k = 2 * root1/(1 + root1 - 4 * c - 8 * c**2); n = 2 * root1/(1 + root1 + 4 * c)

    call f90Elliptic(n, Pio2, k, res);  CParamFO1Loop = dot_product(coef, res)/9

  end function CParamFO1Loop

!ccccccccccccccc

  subroutine f90Elliptic(c, phi, k, res)
    real (dp), intent(in )               :: phi, k, c
    real (dp), intent(out), dimension(3) :: res
    integer                              :: ier
    real (dp)                            :: Sphi, Sphi2, Sphi3, Cphi2

    Sphi = Sin(phi); Sphi2 = Sphi**2; Sphi3 = Sphi**3; Cphi2 = Cos(phi)**2

    res(1)  = Sphi * drf(Cphi2, 1 - k * Sphi2, 1._dp, ier)

    res(2:) = res(1) - Sphi3/3 * [k * drj(Cphi2, 1 - k * Sphi2, 1._dp, 1._dp, ier),      &
                         c * drj( Cphi2, 1._dp, 1 - k * Sphi2, 1 - c * Sphi2, ier) ]

  end subroutine f90Elliptic

!ccccccccccccccc

  subroutine f90Elliptic4(c, c2, phi, k, res)
    real (dp), intent(in )               :: phi, k, c, c2
    real (dp), intent(out), dimension(4) :: res
    integer                              :: ier
    real (dp)                            :: Sphi, Sphi2, Sphi3, Cphi2

    Sphi = Sin(phi); Sphi2 = Sphi**2; Sphi3 = Sphi**3; Cphi2 = Cos(phi)**2

    res(1)  = Sphi * drf(Cphi2, 1 - k * Sphi2, 1._dp, ier)

    res(2:) = res(1) - Sphi3/3 * &
    [  k * drj(Cphi2,        1 - k * Sphi2, 1._dp, 1._dp  , ier) , &
    -  c * drj(Cphi2, 1._dp, 1 - k * Sphi2, 1 - c  * Sphi2, ier) , &
    - c2 * drj(Cphi2, 1._dp, 1 - k * Sphi2, 1 - c2 * Sphi2, ier) ]

end subroutine f90Elliptic4

!ccccccccccccccc

end module MassiveNSClass
