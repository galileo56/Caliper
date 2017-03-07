
module GapMassClass
  use MatrixElementsClass ; use AnomDimClass; use Constants, only: dp, Pi
  use QuadPack, only: qags; use RunningClass; implicit none; private

!ccccccccccccccc

  type, public              ::  GapMass
    private
    real (dp)               :: mm, DeltaRef, rRef, muRef, lambda, muS, muM
    real (dp), dimension(3) :: sCoef
    character (len = 12)    :: gap
    type (Running)          :: run
    integer                 :: order
    class (MatricesElementsMass), allocatable :: MatEl

    contains

   procedure, pass(self), public :: DeltaGapMass, setGammaShift

  end type GapMass

!ccccccccccccccc

  interface GapMass
    module procedure  InitGapMass
  end interface GapMass

  contains

!ccccccccccccccc

   type (GapMass) function InitGapMass(order, gap, MatEl, R0, mu0)
     character (len = *)         , intent(in) :: gap
     integer                     , intent(in) :: order
     class (MatricesElementsMass), intent(in) :: MatEl
     real (dp)                   , intent(in) :: R0, mu0
     real (dp)                                :: mm
     real (dp), dimension(3  )                :: sCoef
     real (dp), dimension(0:4)                :: beta
     type (Running)                           :: run
     type (AnomDim)                           :: andim

     mm = MatEl%scales('mm')  ;  InitGapMass%mm   = mm  ;  InitGapMass%DeltaRef = 0
     InitGapMass%gap = gap    ;  InitGapMass%rRef = R0  ;  InitGapMass%muRef    = mu0
     run = MatEl%run()        ;  InitGapMass%run  = run ;  InitGapMass%order    = order
     InitGapMass%muS  = MatEl%scales('muS'); InitGapMass%muM = MatEl%scales('muM')

     select type (MatEl)
     type is (MatrixElementsMass)

       allocate( MatrixElementsMass :: InitGapMass%MatEl  )
       select type (selector => InitGapMass%MatEl)
       type is (MatrixElementsMass)
       selector = MatEl
       end select

     type is (MatricesElementsMass)

       allocate( MatricesElementsMass :: InitGapMass%MatEl  )
       select type (selector => InitGapMass%MatEl)
       type is (MatricesElementsMass)
       selector = MatEl
       end select

     end select

     if ( InitGapMass%muM <= InitGapMass%muS ) then

       if (R0 < mm) then

         InitGapMass%rRef     = mm ;  InitGapMass%muRef = mm

         InitGapMass%DeltaRef = MatEl%DiffDeltaGapNl(gap, order, R0, mm, mu0, mm) &
                                - run%DeltaGapMatching()

       end if

       sCoef = MatEl%sCoefNl(gap); andim = MatEl%adim('nl'); beta = andim%betaQCD('beta')

       InitGapMass%sCoef(:2) = [ 0._dp, sCoef(2) * beta(0)**2 ]

       sCoef = MatEl%sCoef(gap)  ; andim = MatEl%adim('nf'); beta = andim%betaQCD('beta')

       InitGapMass%sCoef(2:) = [ InitGapMass%sCoef(2)/beta(0)**2, sCoef(3) ]

       InitGapMass%lambda    = run%lambdaQCD(order)

     end if

   end function InitGapMass

!ccccccccccccccc

   real (dp) function DeltaGapMass(self, R, mu)
    class (GapMass), intent(in) :: self
    real (dp)      , intent(in) :: R, mu

    real (dp)                   :: mass

    if (self%muM <= self%muS) then

      mass = self%run%MSbarMass(mu)

      DeltaGapMass = self%DeltaRef + self%lambda * self%run%DiffRMass(self%order, mass, self%rRef, R) &
       + self%run%DiffDelta(self%sCoef, cuspConst(self%gap), self%order, self%rRef, R, self%muRef, mu)

    else
      DeltaGapMass = self%MatEl%DiffDeltaGapNl(self%gap, self%order, self%rRef, R, self%muRef, mu)
    end if

   end function DeltaGapMass

!ccccccccccccccc

  subroutine setGammaShift(self, order, shape, delta0, h, R, mu, shift, delta)
    class (GapMass)    , intent(in)                    :: self
    character (len = *), intent(in)                    :: shape
    integer            , intent(in)                    :: order
    real (dp)          , intent(in)                    :: h, delta0, R, mu
    real (dp)          , intent(out), dimension(order) :: delta
    real (dp)          , intent(out)                   :: shift
    real (dp)                       , dimension(3)     :: del

    del = 0;  delta = 0;  shift = 0

    if ( self%gap(:6) == 'thrust' .or. self%gap(:6) == 'Cparam' .or. self%gap(2:3) == 'JM' ) then

      del   = self%MatEl%delta(self%gap);  delta = del(:order)
      shift = self%DeltaGapMass(R, mu) + gapCons(self%gap) * delta0 + (1 - h) * Sum(delta)
      delta = h * delta

      if (  ( shape(:6) == 'thrust' .or. shape(2:3) == 'JM' ) .and. self%gap(:6) == 'Cparam'  ) then
        shift = 4 * shift/Pi;  delta = 4 * delta/Pi
      else if (  ( self%gap(:6) == 'thrust' .or. self%gap(2:3) == 'JM' ) .and. shape(:6) == 'Cparam'  ) then
        shift = Pi * shift/4;  delta = Pi * delta/4
      end if
    end if

  end subroutine setGammaShift

!ccccccccccccccc

end module GapMassClass
