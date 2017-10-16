

! cccccccccccccccccccc

subroutine f90Delta1S(nl, orderAlp, runAlp, orderM, runM, muLam, xlam, method, &
mZ, aMz, mt, R, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none
  integer  , intent(in)                  :: nl, orderAlp, runAlp, runM, orderM
  real (dp), intent(in)                  :: muLam, xlam, mZ, aMz, mt, R
  character (len = *)       , intent(in) :: method
  real (dp), intent(out), dimension(0:4) :: res
  type (ElectroWeak)                     :: EW
  integer                                :: i
  type (RNRQCD)                          :: NRQCD
  type (Alpha)                           :: alphaAll
  type (Running), dimension(2)           :: alphaMass
  type (VFNSMSR)                         :: MSR
  type (AnomDim), dimension(3:6)         :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, aMz, mt, mt, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, runM, alphaAll, muLam), &
  Running(nl, runM, alphaAll, muLam) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'MSRn', method(:8), 0._dp, orderM, 4, R, xlam)
  res = NRQCD%Delta1S( )

end subroutine f90Delta1S

! cccccccccccccccccccc

subroutine f90Qswitch(nl, orderAlp, runAlp, orderM, runM, ord1S, muLam, xlam, &
method, mZ, aMz, mt, gt, R, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none
  integer            , intent(in) :: nl, orderAlp, runAlp, runM, orderM, ord1S
  real (dp)          , intent(in) :: muLam, xlam, mZ, aMz, mt, R, gt
  character (len = *), intent(in) :: method
  real (dp), intent(out)          :: res
  integer                         :: i
  type (ElectroWeak)              :: EW
  type (RNRQCD)                   :: NRQCD
  type (Alpha)                    :: alphaAll
  type (Running), dimension(2)    :: alphaMass
  type (VFNSMSR)                  :: MSR
  type (AnomDim), dimension(3:6)  :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, aMz, mt, mt, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, runM, alphaAll, muLam), &
  Running(nl, runM, alphaAll, muLam) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'MSRn', method(:8), gt, orderM, ord1S, R, xlam)
  res =  NRQCD%Switch()

end subroutine f90Qswitch

! cccccccccccccccccccc

subroutine f90RNRQCD(nl, order, scheme, method, orderAlp, runAlp, orderMass, &
  runMass, ord1S, R1S, muLam, xlam, mZ, aMz, Q, mtpole, gt, h, nu, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none
  real (dp)          , intent(in)  :: Q, mtpole, gt, h, nu, mZ, aMz, xlam, muLam, R1S
  character (len = *), intent(in)  :: scheme, method
  integer            , intent(in)  :: nl, orderAlp, runAlp, runMass, orderMass, order, ord1S
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, aMz, mtpole, mtpole, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, runMass, alphaAll, muLam), &
  Running(nl, runMass, alphaAll, muLam) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, scheme(:4), method(:8), gt, orderMass, ord1S, R1S, xlam )

  res = NRQCD%Xsec( order, q, h, nu )

end subroutine f90RNRQCD

! cccccccccccccccccccc

subroutine f90A1Pole(nl, order, En, mtpole, gamtop, asoft, VcsNNLL, musoft, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: En, mtpole, gamtop, asoft, VcsNNLL, musoft
  integer            , intent(in)  :: nl, order
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, mtpole, mtpole, 0._dp, 0._dp, &
  0._dp, 0._dp, 'analytic', 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', gamtop, 0, 0, 0._dp, 0._dp)

  res = NRQCD%A1pole(order, En, mtpole, asoft, VcsNNLL, musoft)

end subroutine f90A1Pole

! cccccccccccccccccccc

subroutine f90CDigamma(zr, zi, res)
  use Constants; use DeriGamma; implicit none
  real (dp)              , intent(in)  :: zr, zi
  real (dp), dimension(2), intent(out) :: res
  complex (dp)                         :: z
  z = DIGAM( cmplx(zr, zi, kind = dp) ); res = [ realpart(z), imagpart(z) ]
end subroutine f90CDigamma

! cccccccccccccccccccc

subroutine f90CTrigamma(zr, zi, res)
  use Constants; use DeriGamma; implicit none
  real (dp)              , intent(in)  :: zr, zi
  real (dp), dimension(2), intent(out) :: res
  complex (dp)                         :: z
  z = TrIGAM( cmplx(zr, zi, kind = dp) ); res = [ realpart(z), imagpart(z) ]
end subroutine f90CTrigamma

! cccccccccccccccccccc

subroutine f90qFromV(v, m, gt, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp), intent(in)  :: v, m, gt
  real (dp), intent(out) :: res
  res = qFromV(v, m, gt)
end subroutine f90qFromV

! cccccccccccccccccccc

subroutine f90vC(q, m, gt, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp)              , intent(in)  :: q, m, gt
  real (dp), dimension(2), intent(out) :: res
  complex (dp)                         :: v
  v = vC(q, m, m, gt); res = [ realpart(v), imagpart(v) ]
end subroutine f90vC

! cccccccccccccccccccc

subroutine f90vStar(q, m, gt, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp), intent(in)  :: q, m, gt
  real (dp), intent(out) :: res
  res = vStar(q, m, gt)
end subroutine f90vStar

! cccccccccccccccccccc

subroutine f90vRootStar(q, m, gt, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp), intent(in)  :: q, m, gt
  real (dp), intent(out) :: res
  res = vRootStar(q, m, gt)
end subroutine f90vRootStar

! cccccccccccccccccccc

subroutine f90SwitchOff(q, m, gt, v0, v1, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp), intent(in)  :: q, m, gt, v0, v1
  real (dp), intent(out) :: res
  res = SwitchOff(q, m, gt, v0, v1)
end subroutine f90SwitchOff

! cccccccccccccccccccc

subroutine f90VssLL(nl, ah, as, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass;  use VFNSMSRClass; use ElectroWeakClass; implicit none
  real (dp)          , intent(in)  :: ah, as
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%VssLL(ah, as)

end subroutine f90VssLL

! cccccccccccccccccccc

subroutine f90Vk1sLL(nl, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass;  use VFNSMSRClass; use ElectroWeakClass; implicit none
  real (dp)          , intent(in)  :: au, as
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (RNRQCD)                    :: NRQCD
  type (ElectroWeak)               :: EW
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%Vk1sLL(as, au)

end subroutine f90Vk1sLL

! cccccccccccccccccccc

subroutine f90Vk2sLL(nl, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: au, as
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (RNRQCD)                    :: NRQCD
  type (ElectroWeak)               :: EW
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%Vk2sLL(as, au)

end subroutine f90Vk2sLL

! cccccccccccccccccccc

subroutine f90VkeffsLL(nl, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: au, as
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (RNRQCD)                    :: NRQCD
  type (ElectroWeak)               :: EW
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%VkeffsLL(as, au)

end subroutine f90VkeffsLL

! cccccccccccccccccccc

subroutine f90VcsLL(as, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp), intent(in)  :: as
  real (dp), intent(out) :: res

  res = VcsLL(as)

end subroutine f90VcsLL

! cccccccccccccccccccc

subroutine f90VrsLL(nl, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: as, au
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%VrsLL(as, au)

end subroutine f90VrsLL

! cccccccccccccccccccc

subroutine f90V2sLL(nl, ah, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: as, au, ah
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%V2sLL(ah, as, au)

end subroutine f90V2sLL

! cccccccccccccccccccc

subroutine f90VceffsNNLL(nl, asNNLL, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: as, au, asNNLL
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%VceffsNNLL(asNNLL, as, au)

end subroutine f90VceffsNNLL

! cccccccccccccccccccc

subroutine f90XiNLL(nl, ah, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, as, au
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%XiNLL(ah, as, au)

end subroutine f90XiNLL

! cccccccccccccccccccc

subroutine f90XiNNLLmixUsoft(nl, ah, as, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, as
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%XiNNLLmixUsoft(ah, as)

end subroutine f90XiNNLLmixUsoft

! cccccccccccccccccccc

subroutine f90MLLc2(nl, ah, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, au
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%MLLc2(ah, au)

end subroutine f90MLLc2

! cccccccccccccccccccc

subroutine f90MNLLc1(nl, ah, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, as, au
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%MNLLc1(ah, as, au)

end subroutine f90MNLLc1

! cccccccccccccccccccc

subroutine f90MNLLplusNNLLnonmixc1(nl, ah, as, au, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, as, au
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%MNLLplusNNLLnonmixc1(ah, as, au)

end subroutine f90MNLLplusNNLLnonmixc1

! cccccccccccccccccccc

subroutine f90MNNLLAllc1InclSoftMixLog(nl, ah, as, au, nu, hh, ss, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, as, au, nu, hh, ss
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = Running(5, 0, alphaAll, 0._dp)

  MSR = VFNSMSR(alphaMass)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%MNNLLAllc1InclSoftMixLog(ah, as, au, nu, hh, ss)

end subroutine f90MNNLLAllc1InclSoftMixLog

! cccccccccccccccccccc

subroutine f90XiNNLLSoftMixLogc1(ah, nu, hh, res)
  use Constants; use RNRQCDClass; implicit none
  real (dp), intent(in)  :: ah, nu, hh
  real (dp), intent(out) :: res

  res = XiNNLLSoftMixLogc1(ah, nu, hh)

end subroutine f90XiNNLLSoftMixLogc1

! cccccccccccccccccccc

subroutine f90XiNNLLnonmix(nl, ah, as, au, hh, ss, res)
  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass;  implicit none
  real (dp)          , intent(in)  :: ah, as, au, hh, ss
  integer            , intent(in)  :: nl
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = [ Running(5, 0, alphaAll, 0._dp), Running(5, 0, alphaAll, 0._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(91.187_dp, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, 'pole', 'analytic', 0._dp, 0, 0, 0._dp, 0._dp)

  res = NRQCD%XiNNLLnonmix(ah, as, au, hh, ss)

end subroutine f90XiNNLLnonmix

! cccccccccccccccccccc

subroutine f90Toppik(energy, tm, tg, alphas, scale, cutn, cutv, c0, c1, c2, &
cdeltc, cdeltl, cfullc, cfulll, crm2, kincm, kinca, jknflg, jgcflg, xkincv, &
jvflg, res)
  use constants, only: dp; use ToppikClass; implicit none
  integer  , intent(in)                :: jknflg, jgcflg, jvflg
  real (dp), intent(out), dimension(2) :: res
  real (dp), intent(in)                :: energy, tm, tg, alphas, scale, cutn, &
  c0, c1, c2, cdeltc, cdeltl, cfullc, cfulll, crm2, cutv

  real (dp)                            :: kincm, kinca, xkincv
  type (Toppik)                        :: ttbar

  ttbar = Toppik(5, energy, tm, tg, alphas, scale, cutn, cutv, c0, c1, &
  c2, cdeltc, cdeltl, cfullc, cfulll, crm2, kincm, kinca, jknflg,   &
  jgcflg, xkincv, jvflg)

  res = ttbar%CrossSec()

end subroutine f90Toppik

! cccccccccccccccccccc

subroutine f90ToppikList(energy, tm, tg, alphas, scale, cutn, cutv, c0,  &
c1, c2, cdeltc, cdeltl, cfullc, cfulll, crm2, kincm, kinca, jknflg, &
jgcflg, xkincv, jvflg, res, list)
  use constants, only: dp; use ToppikClass; implicit none
  integer  , intent(in)  :: jknflg, jgcflg, jvflg
  real (dp), intent(out), dimension(2) :: res
  real (dp), intent(out), dimension(5,nmax) :: list
  real (dp), intent(in)  :: energy, tm, tg, alphas, scale, cutn, cutv, &
  c0, c1, c2, cdeltc, cdeltl, cfullc, cfulll, crm2

  real (dp)                     :: kincm, kinca, xkincv
  type (Toppik)                 :: ttbar

  ttbar = Toppik(5, energy, tm, tg, alphas, scale, cutn, cutv, c0, c1, &
  c2, cdeltc, cdeltl, cfullc, cfulll, crm2, kincm, kinca, jknflg,   &
  jgcflg, xkincv, jvflg)

  res = ttbar%CrossSec()

  list(1,:) = ttbar%pGrid(); list(2,:) = ttbar%Weight(); list(3,:) = ttbar%Sigma()
  list(4,:) = realpart( ttbar%Kernel() ); list(5,:) = imagpart( ttbar%Kernel() )

end subroutine f90ToppikList

!ccccccccccccccc

subroutine f90Pi0(zr, zi, res)
  use SigmaClass, only: VacPol0; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol0( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi0

!ccccccccccccccc

subroutine f90Pi0Der(i, zr, zi, res)
  use SigmaClass, only: VacPol0Der; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  integer  , intent(in) :: i
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol0Der( i, cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi0Der

!ccccccccccccccc

subroutine f90Pi1Der(i, zr, zi, res)
  use SigmaClass, only: VacPol1Der; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  integer  , intent(in) :: i
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol1Der( i, cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi1Der

!ccccccccccccccc

subroutine f90Pi1(zr, zi, res)
  use SigmaClass, only: VacPol1; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol1( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi1

!ccccccccccccccc

subroutine f90Pi3(zr, zi, res)
  use SigmaClass, only: VacPol3; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol3( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi3

!ccccccccccccccc

subroutine f90Pi2(zr, zi, res)
  use SigmaClass, only: VacPol2; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol2( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi2

!ccccccccccccccc

subroutine f90Pi2Der(zr, zi, res)
  use SigmaClass, only: VacPol2Der; use constants, only: dp; implicit none
  real (dp), intent(in) :: zr, zi
  real (dp), intent(out), dimension(2) :: res
  complex (dp) :: resul

  resul = VacPol2Der( cmplx(zr, zi, kind = dp) )

  res = [ RealPart(resul), ImagPart(resul) ]

end subroutine f90Pi2Der

!ccccccccccccccc

subroutine f90DeltaBottomCharm(r1, r2, res)
  use constants, only: dp; use AnomDimClass, only: DeltaBottomCharm; implicit none
  real (dp), intent(in ) :: r1, r2
  real (dp), intent(out) :: res

  res = DeltaBottomCharm(r1,r2)

end subroutine f90DeltaBottomCharm

!ccccccccccccccc

subroutine f90gammaRBottomCharm(r1, r2, res)
  use constants, only: dp; use AnomDimClass, only: GammaRBottomCharm; implicit none
  real (dp), intent(in ) :: r1, r2
  real (dp), intent(out) :: res

  res = GammaRBottomCharm(r1,r2)

end subroutine f90gammaRBottomCharm

!ccccccccccccccc

subroutine f90DeltaCharm2(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm2; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharm2(r)

end subroutine f90DeltaCharm2

!ccccccccccccccc

subroutine f90P2(r, res)
  use constants, only: dp; use AnomDimClass, only: p; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = p(r)

end subroutine f90P2

!ccccccccccccccc

subroutine f90P2Int(r, res)
  use constants, only: dp; use AnomDimClass, only: pInt; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = pInt(r)

end subroutine f90P2Int

!ccccccccccccccc

subroutine f90P2Double(r1, r2, res)
  use constants, only: dp; use AnomDimClass, only: P2Int; implicit none
  real (dp), intent(in ) :: r1, r2
  real (dp), intent(out) :: res

  res = P2Int(r1, r2)

end subroutine f90P2Double

!ccccccccccccccc

subroutine f90DeltaCharmNh(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNh; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNh(r)

end subroutine f90DeltaCharmNh

!ccccccccccccccc

subroutine f90DeltaCharmGlue(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmGlue; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmGlue(r)

end subroutine f90DeltaCharmGlue

!ccccccccccccccc

subroutine f90DeltaCharmGlueDer(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmGlueDer; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmGlueDer(r)

end subroutine f90DeltaCharmGlueDer

!ccccccccccccccc

subroutine f90DeltaCharmNhDer(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNhDer; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNhDer(r)

end subroutine f90DeltaCharmNhDer

!ccccccccccccccc

subroutine f90DeltaCharmNL(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNL; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNL(r)

end subroutine f90DeltaCharmNL

!ccccccccccccccc

subroutine f90DeltaCharmNLDer(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharmNLDer; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharmNLDer(r)

end subroutine f90DeltaCharmNLDer

!ccccccccccccccc

subroutine f90DeltaCharm3(nl, nh, r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm3; implicit none
  real (dp), intent(in ) :: r
  integer  , intent(in)  :: nl, nh
  real (dp), intent(out) :: res

  res = deltaCharm3(nl, nh, r)

end subroutine f90DeltaCharm3

!ccccccccccccccc

subroutine f90DeltaCharm3Der(nl, nh, r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm3Der; implicit none
  real (dp), intent(in ) :: r
  integer  , intent(in)  :: nl, nh
  real (dp), intent(out) :: res

  res = deltaCharm3Der(nl, nh, r)

end subroutine f90DeltaCharm3Der

!ccccccccccccccc

subroutine f90GammaRCharm3(nl, nh, r, res)
  use constants, only: dp; use AnomDimClass, only: GammaRCharm3; implicit none
  real (dp), intent(in ) :: r
  integer  , intent(in)  :: nl, nh
  real (dp), intent(out) :: res

  res = GammaRCharm3(nl, nh, r)

end subroutine f90GammaRCharm3

!ccccccccccccccc

subroutine f90GammaRCharm2(r, res)
  use constants, only: dp; use AnomDimClass, only: gammaRcharm2; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = gammaRCharm2(r)

end subroutine f90GammaRCharm2

!ccccccccccccccc

subroutine f90DeltaCharm2Der(r, res)
  use constants, only: dp; use AnomDimClass, only: deltaCharm2Der; implicit none
  real (dp), intent(in ) :: r
  real (dp), intent(out) :: res

  res = deltaCharm2Der(r)

end subroutine f90DeltaCharm2Der

!ccccccccccccccc

subroutine f90LegendreList(n, k, x, res)
  use constants, only: dp; use Legendre; implicit none
  integer                  , intent(in)  :: n, k
  real (dp)                , intent(in)  :: x
  real (dp), dimension(0:n), intent(out) :: res

  res = LegendreList(n,k,x)

end subroutine f90LegendreList

!ccccccccccccccc

subroutine f90QLegendreList(n, x, res)
  use constants, only: dp; use Legendre; implicit none
  integer                  , intent(in)  :: n
  real (dp)                , intent(in)  :: x
  real (dp), dimension(0:n), intent(out) :: res

  res = QLegendreList(n,x)

end subroutine f90QLegendreList

!ccccccccccccccc

subroutine f90MCtop(shape, mt, Q, n, k, x, res)
  use constants, only: dp; use hyper; use MCtopClass; implicit none
  character (len = *), intent(in)  :: shape
  real (dp)          , intent(in)  :: mt, Q, x
  integer            , intent(in)  :: n, k
  real (dp)          , intent(out) :: res
  type (MCScales)                  :: MC

  MC = MCScales( shape(:6), n );     call MC%setMass(mt, Q)
  res = MC%Distribution(k,x)

end subroutine f90MCtop

!ccccccccccccccc

subroutine f90BreitUnstable(shape, mt, Q, gamma, n, k, x, res)
  use constants, only: dp; use hyper; use MCtopClass; implicit none
  character (len = *), intent(in)  :: shape
  real (dp)          , intent(in)  :: mt, Q, x, gamma
  integer            , intent(in)  :: n, k
  real (dp)          , intent(out) :: res
  type (MCtop)                     :: MC

  MC = MCtop( shape(:6), mt, Q, n );  res = MC%BreitUnstable(gamma, k,x)

end subroutine f90BreitUnstable

!ccccccccccccccc

subroutine f90ModelUnstable(shape, mt, Q, c, clen, lambda, n, k, p, res)
  use constants, only: dp; use hyper; use MCtopClass; use ModelClass; implicit none
  character (len = *)       , intent(in)  :: shape
  real (dp)                 , intent(in)  :: mt, Q, p, lambda
  real (dp), dimension(clen), intent(in)  :: c
  integer                   , intent(in)  :: n, clen, k
  real (dp)                 , intent(out) :: res
  type (MCtop)                            :: MC
  type (Model)                            :: Mod

  MC = MCtop( shape(:6), mt, Q, n ); Mod = Model(lambda, c, [0, 0], 'sum')
  res = Mod%ModelUnstable(MC, k, p)

end subroutine f90ModelUnstable

!ccccccccccccccc

subroutine f90BreitModelUnstable(shape, mt, Q, gamma, c, clen, lambda, n, k, p, res)
  use constants, only: dp; use hyper; use MCtopClass; use ModelClass; implicit none
  character (len = *)       , intent(in)  :: shape
  real (dp)                 , intent(in)  :: mt, Q, p, lambda, gamma
  real (dp), dimension(clen), intent(in)  :: c
  integer                   , intent(in)  :: n, clen, k
  real (dp)                 , intent(out) :: res
  type (MCtop)                            :: MC
  type (Model)                            :: Mod

  MC = MCtop( shape(:6), mt, Q, n ); Mod = Model(lambda, c, [0, 0], 'sum')
  res = Mod%BreitModUns(MC, gamma, k, p)

end subroutine f90BreitModelUnstable

!ccccccccccccccc

subroutine f90ModelUnstablediff(shape, mt, Q, c, clen, lambda, n, k, p, VacPol2, res)
  use constants, only: dp; use hyper; use MCtopClass; use ModelClass; implicit none
  character (len = *)       , intent(in)  :: shape
  real (dp)                 , intent(in)  :: mt, Q, p, VacPol2, lambda
  real (dp), dimension(clen), intent(in)  :: c
  integer                   , intent(in)  :: n, clen, k
  real (dp)                 , intent(out) :: res
  type (MCtop)                            :: MC
  type (Model)                            :: Mod

  MC = MCtop( shape(:6), mt, Q, n ); Mod = Model(lambda, c, [0, 0], 'sum')
  res = Mod%ModelUnstable(MC, k, p, VacPol2)

end subroutine f90ModelUnstablediff

!ccccccccccccccc

subroutine f90DeltaMCtop(shape, mt, Q, res)
  use constants, only: dp; use hyper; use MCtopClass; implicit none
  character (len = *), intent(in)  :: shape
  real (dp)          , intent(in)  :: mt, Q
  real (dp)          , intent(out) :: res
  type (MCtop)                     :: MC

  MC = MCtop( shape(:6), mt, Q );  res = MC%delta()

end subroutine f90DeltaMCtop

!ccccccccccccccc

subroutine f90pfq(A, IP, B, IQ, Z, res)
  use constants, only: dp; use hyper; implicit none
  real (dp), dimension(IP), intent(in)  :: A
  real (dp), dimension(IQ), intent(in)  :: B
  integer                 , intent(in)  :: IP, IQ
  real (dp)               , intent(in)  :: Z
  real (dp)               , intent(out) :: res

  res = pFq(a, b, z)

end subroutine f90pfq

!ccccccccccccccc

subroutine f90PolyGamma(n, w, res)
  use constants, only: dp; use DeriGamma; implicit none
  integer      , intent(in) :: n
  real (dp)    , intent(in) :: w
  real (dp), dimension(0:n) :: res
  integer                   ::  NZ, IERR

  call DPSIFN(w, 0, 1, n + 1, res, NZ, IERR)

end subroutine f90PolyGamma

!ccccccccccccccc

subroutine f90MassiveProfList(terms, hard, shape, EShape, setup, gap, space, cum, scheme, &
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia
  use AnomDimClass; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, &
  cum, abs, current, Eshape, terms
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, ns, &
  runMass, orderMass, n, clen
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, &
  muC, Q, G3, lambda, R0, mu0, delta0, h, s3, muR0, j3, Rat0, n0, n1, t2, cnt, eS,&
  slope, eH, eJ, width, ts, xi, xiB, muLambda2, beta, deltaLambda, del0, delta1, &
  muM, mass, gammaZ, sin2ThetaW

  real (dp), dimension(n), intent(in)   :: tauList
  real (dp), dimension(n), intent(out)  :: res
  type (AnomDim), dimension(3:6)        :: AnDim
  type (ElectroWeak)                    :: EW
  type (Alpha)                          :: alphaAll
  type (MatricesElementsMass)           :: MatEl
  type (SingularMassScales)             :: Sing
  type (Model)                          :: Mod
  type (MassiveScales)                  :: nonSing
  type (ProfilesPythia)                 :: Prof
  type (CumulantMass)                   :: Cumul
  character (len = 5)                   :: alphaScheme
  integer                               :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, &
  delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), &
  current(:6), orderMass, matEl, EW )

  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListDist(terms(:7), cum(:4), Mod, setup(:15), gap(:12), space(:6), &
  order, R0, muR0, del0, h, tauList)

end subroutine f90MassiveProfList

!ccccccccccccccc

subroutine f90MassiveProfPieceList(terms, hard, shape, EShape, setup, gap, space, cum, scheme,  &
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3,&
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,&
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, &
  cum, abs, current, setup, Eshape
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, &
  runMass, orderMass, n, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC,&
  G3, lambda, R0, mu0, delta0, s3, muR0, j3, muC, Q, Rat0, n0, n1, t2, cnt, eS,  &
   slope, eH, eJ, width, ts, xi, xiB, muLambda2, beta, deltaLambda, del0, delta1,&
   muM, mass, gammaZ, sin2ThetaW, h
  real (dp), dimension(n), intent(in)   :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out)  :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListDistPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, tauList)

end subroutine f90MassiveProfPieceList

!ccccccccccccccc

subroutine f90MassivePieceBin(terms, hard, shape, EShape, setup, gap, space, cum, scheme,  &
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3,&
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,&
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, Eshape, &
                                        cum, abs, current, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, &
                                        runMass, orderMass, n
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q,&
                                        G3, lambda, R0, mu0, delta0, h, s3, muR0, j3, &
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW
  real (dp), dimension(2,n), intent(in) :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out) :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListBinPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, tauList)

end subroutine f90MassivePieceBin

!ccccccccccccccc

subroutine f90MassiveBinList(terms, hard, shape, EShape, setup, gap, space, cum, scheme, &
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *)        , intent(in) :: shape, scheme, setup, space, gap, &
  cum, abs, current, Eshape, terms, hard
  integer                    , intent(in) :: orderAlp, order, runAlp, run, &
  runMass, orderMass, n, clen, nf, ns
  real (dp)                  , intent(in) :: mZ, amZ, muLambda1, mT, muT, mB,   &
  muC, Q, G3, lambda, R0, mu0, delta0, h, s3, muR0, j3, Rat0, n0, n1, t2, cnt,  &
  slope, eH, eJ, width, ts, xiB, muLambda2, beta, deltaLambda, del0, delta1,&
  muM, mass, gammaZ, sin2ThetaW, muB, mC, eS, xi
  real (dp), dimension(2,n), intent(in) :: tauList
  real (dp), dimension(n) , intent(out) :: res
  integer                               :: i
  type (AnomDim), dimension(3:6)        :: AnDim
  type (ElectroWeak)                    :: EW
  type (Alpha)                          :: alphaAll
  type (MatricesElementsMass)           :: MatEl
  type (SingularMassScales)             :: Sing
  type (Model)                          :: Mod
  type (MassiveScales)                  :: nonSing
  type (ProfilesPythia)                 :: Prof
  type (CumulantMass)                   :: Cumul
  character (len = 5)                   :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
  ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
  orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%ListBin(terms(:7), cum(:9), Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, del0, h, tauList)

end subroutine f90MassiveBinList

!ccccccccccccccc

subroutine f90MassiveProf(terms, hard, shape, EShape, setup, gap, space, cum, scheme, &
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, terms, &
                                        cum, abs, current, Eshape
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (Model)                       :: Mod
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%Bin(terms(:7), cum(:4), Mod, setup(:15), gap(:12), space(:6), order, R0, &
  muR0, del0, h, 0, tau)

end subroutine f90MassiveProf

!ccccccccccccccc

subroutine f90MassOrigin(shape, EShape, gap, scheme, orderAlp, runAlp, orderMass,  &
 runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, &
 beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, &
 mass, muM, R0, muR0, del0, h, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ElectroWeakClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia; use AnomDimClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, gap, Eshape
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        R0, mu0, delta0, h, muR0, Rat0, n0, n1, t2, cnt,  &
                                        eS, slope, eH, ts, muLambda2, beta, deltaLambda,  &
                                        del0, delta1, muM, mass, eJ
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (MatricesElementsMass)        :: MatEl
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, 0._dp, 0._dp)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, 0, EShape(:1), shape(:6) )
  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), 'abs', 'all', &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, 'expand' )
  Cumul     = CumulantMass(Prof, Sing, 0._dp, 0._dp, 0._dp)

  res = Cumul%FindOrigin( gap(:12), order, R0, muR0, del0, h )

end subroutine f90MassOrigin

!ccccccccccccccc

subroutine f90MassiveProfPiece(terms, hard, shape, EShape, setup, gap, space, cum, scheme, abs, &
 current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3, &
 mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,    &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, &
                                        cum, abs, current, Eshape, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%BinPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, 0, tau)

end subroutine f90MassiveProfPiece

!ccccccccccccccc

subroutine f90MassiveProfDiffPiece(terms, hard, shape, EShape, setup, gap, space, cum, scheme,  &
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3,&
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda,&
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width,     &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, scheme, space, gap, hard, terms, &
                                        cum, abs, current, Eshape, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW, tau2
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%BinPiece(terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, del0, h, 0, tau, tau2)

end subroutine f90MassiveProfDiffPiece

!ccccccccccccccc

subroutine f90MassiveProfDiff(terms, hard, shape, EShape, setup, gap, space, cum, scheme,&
 abs, current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, terms, &
                                        cum, abs, current, Eshape
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, &
                                        runMass, orderMass
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW, tau2
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (Model)                       :: Mod
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%Bin(terms(:7), cum(:9), Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, del0, h, 0, tau, tau2)

end subroutine f90MassiveProfDiff

!ccccccccccccccc

subroutine f90MassiveMoment(terms, hard, shape, EShape, setup, gap, space, scheme, abs,  &
 current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3,     &
 G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, deltaLambda, &
 Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, c,   &
 clen, lambda, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau, tau2, pow, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: ProfilesPythia  ; use ElectroWeakClass
  use CumulantClass, only: CumulantMass; use MassiveNSClass; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, scheme, setup, space, gap, hard, terms, &
                                        abs, current, Eshape
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, &
                                        runMass, orderMass, pow
  real (dp)          , intent(in)    :: mZ, amZ, muLambda1, mT, muT, mB, muB, mC, muC, Q, &
                                        G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, j3,&
                                        Rat0, n0, n1, t2, cnt, eS, slope, eH, eJ, width, &
                                        ts, xi, xiB, muLambda2, beta, deltaLambda, del0, &
                                        delta1, muM, mass, gammaZ, sin2ThetaW, tau2
  real (dp)          , intent(out)   :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (Model)                       :: Mod
  type (MassiveScales)               :: nonSing
  type (ProfilesPythia)              :: Prof
  type (CumulantMass)                :: Cumul
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  Prof      = ProfilesPythia( Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, &
                       ts, slope, cnt, eH, eS, eJ, mass, muM, ns, EShape(:1), shape(:6) )
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), EShape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')
  Cumul     = CumulantMass(Prof, Sing, xi, xiB, width)

  res = Cumul%Bin(terms(:7), 'Integrate', Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, del0, h, pow, tau, tau2)

end subroutine f90MassiveMoment

!ccccccccccccccc

subroutine f90MasslessProfList(terms, hard, shape, setup, gap, space, cum, orderAlp, &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda,   &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda,    &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ
  real (dp), dimension(n), intent(in ) :: tauList
  real (dp), dimension(n), intent(out) :: res
  type (Alpha)                        :: alphaAll
  type (MatricesElements)             :: MatEl
  type (SingularScales)               :: Sing
  type (Model)                        :: Mod
  type (ProfilesMassless)             :: Prof
  type (CumulantMassless)             :: Cumul
  type (AnomDim), dimension(3:6)      :: AnDim
  integer                             :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, &
  eS, eJ, eR, ns)

  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListDist( terms(:7), cum(:4), Mod, setup(:15), gap(:12), space(:6), &
  order, R0, muR0, delta0, h, tauList )

end subroutine f90MasslessProfList

!ccccccccccccccc

subroutine f90MasslessProfPieceList(terms, hard, shape, setup, gap, space, cum, orderAlp, &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda,     &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ
  real (dp), dimension(n), intent(in ) :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out) :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListDistPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, delta0, h, tauList )

end subroutine f90MasslessProfPieceList

!ccccccccccccccc

subroutine f90MasslessPieceBin(terms, hard, shape, setup, gap, space, cum, orderAlp, &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda,   &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda,    &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS, &
  eR, tR, slope, eH, eJ, muC, j3
  real (dp), dimension(2,n), intent(in ) :: tauList
  real (dp), dimension((clen + 2) * (clen + 1)/2,n), intent(out) :: res

  type (AnomDim), dimension(3:6)     :: AnDim
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListBinPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), order, R0, &
                  muR0, delta0, h, tauList )

end subroutine f90MasslessPieceBin

!ccccccccccccccc

subroutine f90MasslessBinList(terms, hard, shape, setup, gap, space, cum, orderAlp, &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda,  &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda,   &
 R0, muR0, delta0, h, tauList, n, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, n, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC,&
  Q, G3, lambda, ts, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS, &
  eR, tR, slope, eH, eJ, muC, j3
  real (dp), dimension(2,n), intent(in ) :: tauList
  real (dp), dimension(n)  , intent(out) :: res
  type (Alpha)                           :: alphaAll
  type (MatricesElements)                :: MatEl
  type (SingularScales)                  :: Sing
  type (Model)                           :: Mod
  type (ProfilesMassless)                :: Prof
  type (CumulantMassless)                :: Cumul
  type (AnomDim), dimension(3:6)         :: AnDim
  integer                                :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, &
  eS, eJ, eR, ns)

  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%ListBin( terms(:7), cum(:9), Mod, setup(:8), gap(:12), space(:6), order, &
                       R0, muR0, delta0, h, tauList )

end subroutine f90MasslessBinList

!ccccccccccccccc

subroutine f90MasslessProf(terms, hard, shape, setup, gap, space, cum, orderAlp,   &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda,  &
 R0, muR0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS, &
  eR, tR, slope, eH, eJ, ts, muC, j3
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, &
  eS, eJ, eR, ns)

  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%Bin( terms(:7), cum(:4), Mod, setup(:8), gap(:12), space(:6), order, R0, &
  muR0, delta0, h, 0, tau )

end subroutine f90MasslessProf

!ccccccccccccccc

subroutine f90MasslessProfPiece(terms, hard, shape, setup, gap, space, cum, orderAlp,   &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda,  &
 R0, muR0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ, &
                                        ts
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  integer                            :: i, j, k
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%BinPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), &
  order, R0, muR0, delta0, h, 0, tau )

end subroutine f90MasslessProfPiece

!ccccccccccccccc

subroutine f90MasslessProfDiffPiece(terms, hard, shape, setup, gap, space, cum, orderAlp, &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, &
 Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambda, R0, &
 muR0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, cum, space, gap, hard, terms, setup
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ, &
                                        ts, tau2
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out)   :: res

  integer                            :: i, j, k
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%BinPiece( terms(:7), cum(:4), Mod, setup, gap(:12), space(:6), &
  order, R0, muR0, delta0, h, 0, tau, tau2 )

end subroutine f90MasslessProfDiffPiece

!ccccccccccccccc

subroutine f90MasslessProfDiff(terms, hard, shape, setup, gap, space, cum, orderAlp,  &
 runAlp, order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
 mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda, R0,    &
 muR0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, Rat0, n0, n1, t2, cnt, eS,   &
  eR, tR, slope, eH, eJ, ts, tau2, muC, j3
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%Bin(terms(:7), cum(:9), Mod, setup(:15), gap(:12), space(:6), &
  order, R0, muR0, delta0, h, 0, tau, tau2)

end subroutine f90MasslessProfDiff

!ccccccccccccccc

subroutine f90MasslessMoment(terms, hard, shape, setup, gap, space, orderAlp, runAlp,&
 order, run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, &
 n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, clen, lambda, R0, muR0, delta0,  &
 h, tau, tau2, pow, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, setup, space, gap, hard, terms
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen, ns, pow
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, j3, &
                                        Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muR0, &
                                        Rat0, n0, n1, t2, cnt, eS, eR, tR, slope, eH, eJ, &
                                        ts, tau2
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl    = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod      = Model(lambda, c, [0,0], 'sum')
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%Bin(terms(:7), 'Integrate', Mod, setup(:15), gap(:12), space(:6), order, R0, &
                  muR0, delta0, h, pow, tau, tau2)

end subroutine f90MasslessMoment

!ccccccccccccccc

subroutine f90FindOrigin(shape, gap, orderAlp, runAlp, order, run, nf, mZ, amZ, mT, &
 muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, &
 eR, R0, muR0, delta0, h, res)

  use AlphaClass;  use MatrixElementsClass;  use SingularClass
  use constants, only: dp; use ProfilesClass, only: profilesmassless
  use CumulantClass, only: CumulantMassless; use AnomDimClass; implicit none

  character (len = *), intent(in)    :: shape, gap
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, muC, &
                                        Q, R0, mu0, delta0, h, muR0, ts, eH, slope, n0, &
                                        Rat0, n1, t2, cnt, eS, eR, tR
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (ProfilesMassless)            :: Prof
  type (CumulantMassless)            :: Cumul
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  do i = 3, 6
    AnDim = AnomDim('MSbar', i, 0._dp)
  end do

  Prof     = ProfilesMassless(Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, 0._dp, eR, 0)
  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl    = MatricesElements(alphaAll, nf, 0._dp, 0._dp, 0._dp, muLambda)
  Sing     = SingularScales( MatEl, run, shape(:6), 'expand' )
  Cumul    = CumulantMassless(Prof, Sing)

  res = Cumul%FindOrigin( gap(:12), order, R0, muR0, delta0, h)

end subroutine f90FindOrigin

!ccccccccccccccc

subroutine f90Profiles(Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, &
                       tau, res)
  use constants, only: dp; use ProfilesClass, only: profilesmassless; implicit none
  integer                 , intent(in) :: ns
  real (dp)               , intent(in) :: mu0, R0, n0, n1, t2, cnt, eS, eR, tR, &
  Q, eH, eJ, tau, ts, slope
  real (dp), dimension(5), intent(out) :: res

  real (dp)                            :: muS, muJ, muH, muNS, R
  type (ProfilesMassless)              :: Prof

  Prof = ProfilesMassless(Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns)

  call Prof%Scales(tau, muNS, muH, muJ, muS, R)

  res = [ muH, muJ, muS, muNS, R ]

end subroutine f90Profiles

!ccccccccccccccc

subroutine f90ProfilesMass(Q, beta, mu0, deltaLambda, R0, n0, delta0, n1, delta1, &
  t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, def, EShape, tau, res)

  use constants, only: dp; use ProfilesClass, only: profilesPythia; implicit none

  character (len = *)    , intent(in ) :: EShape, def
  integer                , intent(in ) :: ns
  real (dp)              , intent(in ) :: Q, mu0, R0, n0, delta0, n1, t2, cnt, ts, slope, &
                                          eJ, mass, beta, delta1, deltaLambda, muM,   &
                                          eH, eS, tau
  real (dp), dimension(6), intent(out) :: res
  real (dp)                            :: muS, muJ, muH, muNS, R
  type (ProfilesPythia)                :: Prof

  Prof = ProfilesPythia(Q, beta, mu0, deltaLambda, R0, n0, delta0, n1, delta1, t2, ts, &
                        slope, cnt, eH, eS, eJ, mass, muM, ns, def, EShape)

  call Prof%Scales(tau, muNS, muH, muJ, muS, R)

  res = [ muH, muJ, muS, mum, muNS, R ]

end subroutine f90ProfilesMass

!ccccccccccccccc

subroutine f90SingularMass(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,      &
  current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3, &
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,   &
  muM, mu, width, c, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass     ; use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass    ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, space, gap, Eshape, abs, hard,&
                                        current
  integer            , intent(in   ) :: orderAlp, order, runAlp, run, nf, clen,      &
                                        runMass, orderMass
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, xiB, gammaZ, sin2ThetaW, width, mC, muC, j3, &
                                        Rmass
  real (dp)          , intent(inout) :: muH, muJ, muS, R, muM
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC )

  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), &
  current(:6), orderMass, matEl, EW )

  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model( lambda, c, [0,0], 'sum' )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                               mu0, delta0, h, tau )
  else

    res = Sing%SingleSingWidth( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, setup(:15), gap(:12), space(:6), cum(:4), &
                                           order, R0, mu0, delta0, h, tau)

end subroutine f90SingularMass

!ccccccccccccccc

subroutine f90SingularMassDiff(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, &
  current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3,&
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,  &
  muM, mu, width, c, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, space, gap, Eshape, abs, hard,&
                                        current
  integer            , intent(in   ) :: orderAlp, order, runAlp, run, nf, clen,      &
                                        runMass, orderMass
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, xiB, gammaZ, sin2ThetaW, width, mC, muC, j3, &
                                        tau2, Rmass
  real (dp)          , intent(inout) :: muH, muJ, muS, R, muM
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  character (len = 5)                :: alphaScheme
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model( lambda, c, [0,0], 'sum' )

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                               mu0, delta0, h, tau, tau2 )
  else

    res = Sing%SingleSingWidth( Mod, setup(:15), gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau, tau2 )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, setup(:15), gap(:12), space(:6), cum(:4), order,&
                                               R0, mu0, delta0, h, tau, tau2)
end subroutine f90SingularMassDiff

!ccccccccccccccc

subroutine f90SingularMassPiece(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,        &
  current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3, &
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,   &
  muM, mu, width, piece, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, abs, current,  &
                                        hard, scheme, setup
  integer              , intent(in ) :: orderAlp, order, runAlp, run, nf, runMass, &
                                        orderMass
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, gammaZ, sin2ThetaW, width, mC, muC, j3, muH, &
                                        muJ, muS, R, muM, xiB, Rmass
  real (dp)            , intent(out) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
  muT, mB, muB, mC, muC )

  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
  muLambda1, muLambda2 )

  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
  orderMass, matEl, EW )

  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], piece, 'piece')

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0,   &
                           delta0, h, tau )
  else

    res = Sing%SingleSingWidth( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0, delta0, &
                                h, tau )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, 'Model', gap(:12), space(:6), cum(:4), order,&
  R0, mu0, delta0, h, tau)

end subroutine f90SingularMassPiece

!ccccccccccccccc

subroutine f90SingularMassList(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,&
  xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,   &
  width, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, abs, current,  &
                                        hard, scheme, setup
  integer              , intent(in ) :: orderAlp, order, runAlp, run, nf, runMass, &
                                        orderMass, clen
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, gammaZ, sin2ThetaW, width, mC, muC, j3, muH, &
                                        muJ, muS, R, muM, xiB, Rmass

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  integer                                              :: i, j, k
  type (ElectroWeak)                                   :: EW
  type (Alpha)                                         :: alphaAll
  type (MatricesElementsMass)                          :: MatEl
  type (SingularMassScales)                            :: Sing
  type (MassiveScales)                                 :: nonSing
  type (AnomDim), dimension(3:6)                       :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)                                  :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, gap(:12), space(:6), cum(:4), order, R0,   &
    mu0, delta0, h, tau )
  else

    res = Sing%SingleSingWidth( Mod, setup, gap(:12), space(:6), cum(:4), order, R0,   &
    mu0, delta0, h, tau )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, gap(:12), space(:6), cum(:4), order, &
  R0, mu0, delta0, h, tau)

end subroutine f90SingularMassList

!ccccccccccccccc

subroutine f90SingularMassDiffList(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,&
  current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3, &
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,   &
  muM, mu, width, clen, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, current,  &
  hard, scheme, setup, abs
  integer              , intent(in ) :: orderAlp, order, runAlp, run, nf, runMass, &
                                        orderMass, clen
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, mB, muB, &
                                        mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, &
                                        xi, gammaZ, sin2ThetaW, width, mC, muC, j3, muH, &
                                        muJ, muS, R, muM, tau2, xiB, Rmass

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  integer                                              :: i, j, k
  type (ElectroWeak)                                   :: EW
  type (Alpha)                                         :: alphaAll
  type (MatricesElementsMass)                          :: MatEl
  type (SingularMassScales)                            :: Sing
  type (MassiveScales)                                 :: nonSing
  type (AnomDim), dimension(3:6)                       :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)                                  :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau, tau2 )
  else

    res = Sing%SingleSingWidth( Mod, setup, gap(:12), space(:6), cum(:4), order, R0,   &
                                mu0, delta0, h, tau, tau2 )
  end if
  if (muJ >= muM) res = res + Sing%NonDist(Mod, gap(:12), space(:6), cum(:4), order, &
                                               R0, mu0, delta0, h, tau, tau2)
end subroutine f90SingularMassDiffList

!ccccccccccccccc

subroutine f90SingularMassDiffPiece(hard, shape, Eshape, setup, gap, space, cum, scheme, abs,   &
  current, xi, xiB, orderAlp, runAlp, orderMass, runMass, order, run, nf, j3, s3, G3,&
  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,  &
  muM, mu, width, piece, lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass     ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, cum, space, gap, Eshape, current,  &
  hard, scheme, setup, abs
  integer              , intent(in ) :: orderAlp, order, runAlp, nf, runMass, &
  orderMass, run
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, muT, &
  mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, mB, muB, xi, gammaZ, &
  sin2ThetaW, width, mC, muC, j3, muH, muJ, muS, R, muM, tau2, xiB, Rmass
  real (dp)            , intent(out) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,    &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, s3, s3, j3, j3,     &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), abs(:3), current(:6), &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], piece, 'piece')

    call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, width)
    call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
    call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, xi, xiB)

  if ( width <= d1mach(1) ) then

    res = Sing%SingleSing( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0, &
                           delta0, h, tau, tau2 )
  else

    res = Sing%SingleSingWidth( Mod, 'Model'//setup, gap(:12), space(:6), cum(:4), order, R0, mu0, delta0, &
                                h, tau, tau2 )
  end if

  if (muJ >= muM) res = res + Sing%NonDist(Mod, 'Model', gap(:12), space(:6), cum(:4), order, &
                                               R0, mu0, delta0, h, tau)
end subroutine f90SingularMassDiffPiece

!ccccccccccccccc

subroutine f90MassNonDist(hard, shape, Eshape, setup, gap, space, cum, scheme,  &
  orderAlp, runAlp, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c, &
  clen, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass     ; use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ; use ElectroWeakClass    ;  use constants, only: dp, d1mach
  use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, gap, hard, Eshape, space
  integer            , intent(in   ) :: orderAlp, order, runAlp, run, clen, runMass, &
                                        orderMass, nf
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, &
                                        mC, muC, muT
  real (dp)          , intent(inout) :: muH, muS, R, muM, muJ, Rmass
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW       = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl    = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
                                    muLambda1, muLambda2 )
  nonSing  = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                              orderMass, matEl, EW )
  Sing     = SingularMassScales( nonSing, run, hard(:6) )
  Mod      = Model( lambda, c, [0,0], 'sum' )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res = Sing%NonDist(Mod, setup(:15), gap(:12), space(:3), cum(:4), order, R0, &
                     mu0, delta0, h, tau)

end subroutine f90MassNonDist

!ccccccccccccccc

subroutine f90MassNonDistDiff(hard, shape, Eshape, setup, gap, space, cum, scheme, &
  orderAlp, runAlp, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,  &
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c,  &
  clen, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in   ) :: shape, cum, setup, gap, hard, Eshape, space
  integer            , intent(in   ) :: orderAlp, order, runAlp, run, clen, runMass, &
                                        orderMass, nf
  real (dp)          , intent(in   ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, h, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, tau2,&
                                        mC, muC, muT
  real (dp)          , intent(inout) :: muH, muS, R, muM, muJ, Rmass
  character (len = *), intent(inout) :: scheme
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu; muM = mu
    if ( scheme(:4) /= 'pole' ) scheme = 'MSbar'
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
                                    muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                              orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model( lambda, c, [0,0], 'sum' )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res = Sing%NonDist(Mod, setup(:15), gap(:12), space(:3), cum(:4), order, R0, &
                     mu0, delta0, h, tau, tau2)

end subroutine f90MassNonDistDiff

!ccccccccccccccc

subroutine f90MassNonDistPiece(hard, shape, Eshape, gap, space, cum, scheme, &
  orderAlp, runAlp, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, &
  mu, c, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer, dimension(2), intent(in ) :: c
  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlp, runAlp, run, nf, runMass, &
                                        orderMass, order
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, &
                                        mC, muC, muT, muH, muS, R, muM, muJ, Rmass
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,     &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], c, 'piece')

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, 'Model', gap(:12), space(:3), cum(:4), order, &
  R0, mu0, delta0, h, tau )

end subroutine f90MassNonDistPiece

!ccccccccccccccc

subroutine f90MassNonDistList(hard, shape, Eshape, gap, space, cum, scheme, &
  orderAlp, runAlp, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, &
  mu, clen, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlp, runAlp, run, nf, runMass, &
                                        orderMass, clen, order
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muB, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, &
                                        mC, muC, muT, muH, muS, R, muM, muJ, Rmass
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out  ) :: res

  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (AnomDim), dimension(3:6)     :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT, &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, &
                            delta0, h, tau )

end subroutine f90MassNonDistList

!ccccccccccccccc

subroutine f90MassNonDistDiffList(hard, shape, Eshape, gap, space, cum, scheme,  &
  orderAlp, runAlp, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,&
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,   &
  clen, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlp, runAlp, run, nf, runMass, &
                                        orderMass, clen, order
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muJ, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, tau2,   &
                                        mC, muC, muT, muH, muS, R, muM, muB, Rmass, h
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out  ) :: res

  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i, j, k
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,     &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, delta0, h, tau, tau2 )

end subroutine f90MassNonDistDiffList

!ccccccccccccccc

subroutine f90MassNonDistDiffPiece(hard, shape, Eshape, gap, space, cum, scheme, &
  orderAlp, runAlp, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,&
  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c,&
  lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass     ;  use MatrixElementsClass ;  use SingularClass ;  use ModelClass
  use MassiveNSClass ;  use ElectroWeakClass    ;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer, dimension(2), intent(in ) :: c
  character (len = *)  , intent(in ) :: hard, shape, cum, gap, Eshape, scheme, space
  integer              , intent(in ) :: orderAlp, order, runAlp, run, nf, runMass, &
                                        orderMass
  real (dp)            , intent(in ) :: mZ, amZ, muLambda1, muLambda2, mT, mB, muJ, &
                                        mu, Q, G3, lambda, tau, R0, mu0, delta0, h, tau2,&
                                        mC, muC, muT, muH, muS, R, muM, muB, Rmass
  real (dp)          , intent(out  ) :: res
  type (ElectroWeak)                 :: EW
  type (Alpha)                       :: alphaAll
  type (MatricesElementsMass)        :: MatEl
  type (SingularMassScales)          :: Sing
  type (MassiveScales)               :: nonSing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  character (len = 5)                :: alphaScheme
  integer                            :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, G3)
  end do

  EW        = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,     &
                     muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass( alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp,  &
               muLambda1, muLambda2 )
  nonSing   = MassiveScales( shape(:6), Eshape(:1), scheme(:10), 'noAbs', 'vector', &
                             orderMass, matEl, EW )
  Sing      = SingularMassScales( nonSing, run, hard(:6) )
  Mod       = Model(lambda, [1._dp], c, 'piece')

  call Sing%SetEverything(muH, Q, muH, muJ, muS, R, Rmass, muM, 0._dp)
  call Sing%setHard(Q, muH);  call Sing%setHardMass(muM)
  call Sing%SetRunning(muJ, muS, R, mu); call Sing%SetMat(muJ, muS, 0._dp, 0._dp)

  res       = Sing%NonDist( Mod, 'Model', gap(:12), space(:3), cum(:4), order, R0, mu0, delta0, h, &
                            tau, tau2 )

end subroutine f90MassNonDistDiffPiece

!ccccccccccccccc

subroutine f90NSMass(shape, setup, gap, cum, scheme, abs, current, orderAlp, runAlp, &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h, &
  gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass;  use ModelClass
  use constants, only: dp; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in ) :: shape, current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nf, clen, runMass, &
                                      orderMass
  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, mu,&
                                      muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muJ,  &
                                      delta0, h, R, width, mB, muB, muM, muLambda2
  real (dp)          , intent(out) :: res
  type (MassiveScales)             :: nonSing
  type (ElectroWeak)               :: EW
  type (MatricesElementsMass)      :: MatEl
  type (Alpha)                     :: alphaAll
  type (Model)                     :: Mod
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, c, [0,0], 'sum')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:10), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass( Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, delta0, &
                        h, t )

end subroutine f90NSMass

!ccccccccccccccc

subroutine f90NSMassDiff(shape, setup, gap, cum, scheme, abs, current, orderAlp, runAlp, &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h, &
  gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use ModelClass; use AnomDimClass
  use constants, only: dp; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in ) :: shape, current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nf, clen, runMass, &
                                      orderMass
  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, t2, &
                                      muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muB, h,&
                                      delta0, R, width, mB, muJ, muM, muLambda2, mu
  real (dp)          , intent(out) :: res
  type (MassiveScales)             :: nonSing
  type (ElectroWeak)               :: EW
  type (MatricesElementsMass)      :: MatEl
  type (Alpha)                     :: alphaAll
  type (Model)                     :: Mod
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, c, [0,0], 'sum')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, &
  muT, mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
  muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), &
  orderMass, matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass( Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, &
   delta0, h, t, t2 )

end subroutine f90NSMassDiff

!ccccccccccccccc

subroutine f90NSMassPiece(shape, gap, cum, scheme, abs, current, orderAlp, runAlp,   &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, piece, width, lambda, R0, mu0, delta0, h,   &
  gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme
  integer              , intent(in ) :: orderAlp, runAlp, order, run, nf, runMass,   &
                                        orderMass
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda ,&
                                        muJ, delta0, R, h, mT, mB, width, muM,      &
                                        muLambda2
  real (dp)           , intent(out) :: res
  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  character (len = 5)               :: alphaScheme
  integer                           :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, [1._dp], piece, 'piece')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, 'Model', gap(:12), cum(:4), order, &
                       run, R0, mu0, delta0, h, t)

end subroutine f90NSMassPiece

!ccccccccccccccc

subroutine f90NSMassList(shape, setup, gap, cum, scheme, abs, current, orderAlp, &
runAlp, orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC,  &
muLambda1, muLambda2, Q, mu, muM, muJ, muS, R, Rmass, clen, width, lambda, R0, mu0,&
delta0, h, gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme, setup
  integer              , intent(in ) :: orderAlp, runAlp, order, run, nf, runMass,   &
                                        orderMass, clen
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda , &
                                        muJ, delta0, R, h, mT, mB, width, muM,       &
                                        muLambda2

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (AnomDim), dimension(3:6)    :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)               :: alphaScheme
  integer                           :: i, j, k

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, delta0, h, t)

end subroutine f90NSMassList

!ccccccccccccccc

subroutine f90NSMassDiffList(shape, setup, gap, cum, scheme, abs, current, orderAlp, &
runAlp, orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, &
muLambda1, muLambda2, Q, mu, muM, muJ, muS, R, Rmass, clen, width, lambda, R0, mu0,&
delta0, h, gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme, setup
  integer              , intent(in ) :: orderAlp, runAlp, order, run, nf, &
                                        orderMass, runMass, clen
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda, &
                                        muJ, delta0, R, h, mT, mB, width, muM, t2, &
                                        muLambda2

  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (AnomDim), dimension(3:6)    :: AnDim
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod
  character (len = 5)               :: alphaScheme
  integer                           :: i, j, k

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, setup(:15), gap(:12), cum(:4), order, run, R0, mu0, delta0, h, t, t2)

end subroutine f90NSMassDiffList

!ccccccccccccccc

subroutine f90NSMassDiffPiece(shape, gap, cum, scheme, abs, current, orderAlp, runAlp,   &
  orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,&
  muLambda2, Q, mu, muM, muJ, muS, R, Rmass, piece, width, lambda, R0, mu0, delta0, h,   &
  gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  integer, dimension(2), intent(in ) :: piece
  character (len = *)  , intent(in ) :: shape, current, cum, gap, abs, scheme
  integer              , intent(in ) :: orderAlp, runAlp, order, run, nf, runMass,   &
                                        orderMass
  real (dp)            , intent(in ) :: t, muT, mC, muC, Q, mZ, gammaZ, sin2ThetaW, muB, &
                                        muS, Rmass, R0, mu0, muLambda1, mu, amZ, lambda, &
                                        muJ, delta0, R, h, mT, mB, width, muM, t2,  &
                                        muLambda2
  real (dp)           , intent(out) :: res
  type (MassiveScales)              :: nonSing
  type (ElectroWeak)                :: EW
  type (MatricesElementsMass)       :: MatEl
  type (Alpha)                      :: alphaAll
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  character (len = 5)               :: alphaScheme
  integer                           :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  Mod       = Model(lambda, [1._dp], piece, 'piece')
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha( AnDim, orderAlp, runAlp, mZ, amZ, mT,  &
                    muT, mB, muB, mC, muC )
  MatEl     = MatricesElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, 0._dp, &
                                    muLambda1, muLambda2)
  nonSing   = MassiveScales( shape(:6), 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
                             matEl, EW)

  call nonSing%SetEverything(mu, Q, mu, muJ, muS, R, Rmass, muM, width)

  res = nonSing%NSMass(Mod, 'Model', gap(:12), cum(:4), order, &
                       run, R0, mu0, delta0, h, t, t2)

end subroutine f90NSMassDiffPiece

!ccccccccccccccc

subroutine f90HJMNSMass(setup, gap, cum, scheme, abs, current, orderAlp, runAlp,     &
 orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, &
 muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h,  &
 gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass; use MatrixElementsClass
  use ModelClass         ; use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in ) :: current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nf, clen, &
  orderMass, runMass

  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, sin2ThetaW, &
  mu, muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muB, delta0, R, width, mB, &
  muJ, muM, muLambda2, h, gammaZ

  real (dp)             , intent(out) :: res
  integer                             :: i, j
  type (AnomDim), dimension(3:6)      :: AnDim
  type (MassiveNS)                    :: nonSing
  type (ElectroWeak)                  :: EW
  type (Alpha)                        :: alphaAll
  type (MatrixElementsMass)           :: MatEl
  type (Model), dimension(clen, clen) :: Mod
  character (len = 5)                 :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp,  &
  0._dp, Q, mu, muJ, muM, muS, R, Rmass, muLambda1, muLambda2)

  nonSing = MassiveNS('HJM', 'Q', scheme(:7), abs(:3), current(:6), orderMass, &
  matEl, EW, width, mu)

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = nonSing%HJMNSMass(c, Mod, setup(:15), gap(:12), cum(:4), &
  order, run, R0, mu0, delta0, h, t)

end subroutine f90HJMNSMass

!ccccccccccccccc

subroutine f90HJMNSMassDiff(setup, gap, cum, scheme, abs, current, orderAlp, runAlp, &
 orderMass, runMass, order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1,     &
 muLambda2, Q, mu, muM, muJ, muS, R, Rmass, width, c, clen, lambda, R0, mu0, delta0, h, &
 gammaZ, sin2ThetaW, t, t2, res)

  use MassiveNSClass     ; use ElectroWeakClass; use AlphaClass
  use MatrixElementsClass; use AnomDimClass    ; use ModelClass
  use constants, only: dp; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in ) :: current, cum, setup, gap, abs, scheme
  integer            , intent(in ) :: orderAlp, runAlp, order, run, clen, &
  orderMass, runMass, nf

  real (dp)          , intent(in ) :: t, mT, muT, mC, muC, Q, mZ, gammaZ, mu, &
  muS, Rmass, R0, mu0, muLambda1, amZ, lambda, muB, h, delta0, R, width, mB, &
  muJ, muM, muLambda2, t2, sin2ThetaW

  real (dp)             , intent(out) :: res
  type (AnomDim), dimension(3:6)      :: AnDim
  integer                             :: i, j
  type (MassiveNS)                    :: nonSing
  type (ElectroWeak)                  :: EW
  type (Alpha)                        :: alphaAll
  type (MatrixElementsMass)           :: MatEl
  type (Model), dimension(clen, clen) :: Mod
  character (len = 5)                 :: alphaScheme

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, &
  0._dp, Q, mu, muJ, muM, muS, R, Rmass, muLambda1, muLambda2)

  nonSing = MassiveNS('HJM', 'Q', scheme(:7), abs(:3), current(:6), orderMass,  &
  matEl, EW,width, mu)

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = nonSing%HJMNSMass(c, Mod, setup(:15), gap(:12), cum(:4), &
  order, run, R0, mu0, delta0, h, t, t2)

end subroutine f90HJMNSMassDiff

!ccccccccccccccc

subroutine f90FOMass(shape, current, m, Q, mZ, gammaZ, sin2ThetaW, t, res)

  use MassiveNSClass     ;  use ElectroWeakClass; use MatrixElementsClass
  use AlphaClass         ;  use constants, only: dp; use AnomDimClass
  implicit none

  real (dp)          , intent(in ) :: t, m, Q, mZ, gammaZ, sin2ThetaW
  character (len = *), intent(in ) :: shape, current
  real (dp)          , intent(out) :: res
  type (MassiveNS)                 :: nonSing
  type (ElectroWeak)               :: EW
  type (MatrixElementsMass)        :: MatEl
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, m, m, m, m, m, m, 'analytic', 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, Q, &
  Q, Q, tiny(1._dp), Q, Q, Q, 0._dp, 0._dp)

  nonSing   = MassiveNS( shape(:6), 'Q', 'pole', 'noAbs', current(:6), 0, matEl, &
  EW, 0._dp, Q)

  res = nonSing%FOMass(t)

end subroutine f90FOMass

!ccccccccccccccc

subroutine f90DiffDeltaGapMass(gap, order, R0, R, mu0, mu, muM, muLambda1, muLambda2, &
  orderAlp, runAlp, runMass, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, res)

  use MatrixElementsClass; use AlphaClass  ; use GapMassClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: gap
  integer            , intent(in ) :: order, nf, orderAlp, runAlp, runMass
  real (dp)          , intent(in ) :: R0, R, mu0, mu, muLambda1, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, muLambda2, muM
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaAll
  type (MatrixElementsMass)        :: MatEl
  type (AnomDim), dimension(3:6)   :: AnDim
  type (GapMass)                   :: MassGap
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl    = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp,  &
  0._dp, mu, mu, mu, muM, mu, mu, mu, muLambda1, muLambda2)

  MassGap = GapMass(order, gap(:12), MatEl, R0, mu0)

  res = MassGap%DeltaGapMass(R, mu)

end subroutine f90DiffDeltaGapMass

!ccccccccccccccc

subroutine f90ThrustNS1loop(t, res)
  use NonSingularClass; use constants, only: dp; implicit none
  real (dp), intent(in )               :: t
  real (dp), intent(out), dimension(3) :: res

  res = ThrustNS1loop(t)

end subroutine f90ThrustNS1loop

!ccccccccccccccc

subroutine f90ThrustNS2loop(er, t, res)
  use NonSingularClass; use AlphaClass;  use RunningClass; use AnomDimClass
  use constants, only: dp; implicit none

  real (dp), intent(in )               :: t, er
  real (dp), intent(out), dimension(2) :: res
  type (Alpha)                         :: alphaAll
  type (Running)                       :: alphaMass
  type (NonSingular)                   :: nonSing
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = Running(5, 0, alphaAll, 0._dp)
  nonSing   = NonSingular('thrust', alphaMass, 0._dp, 0._dp, 0)
  res       = nonSing%NS2loop(er, t)

end subroutine f90ThrustNS2loop

!ccccccccccccccc

subroutine f90EWfactors(nf, Q, mZ, gammaZ, sin2ThetaW, res)
  use ElectroWeakClass;  use constants, only: dp; implicit none
  integer                , intent(in ) :: nf
  real (dp)              , intent(in ) :: Q, mZ, gammaZ, sin2ThetaW
  real (dp), dimension(2), intent(out) :: res
  type (ElectroWeak)                   :: EW

  EW = ElectroWeak(mZ, gammaZ, sin2ThetaW); res = EW%EWfactors(nf, Q)

end subroutine f90EWfactors

!ccccccccccccccc

subroutine f90PolyLog(n, z, res)
  use Poly;  use constants, only: dp; implicit none
  integer  , intent(in ) :: n
  real (dp), intent(in ) :: z
  real (dp), intent(out) :: res

  res = polylog(n,z)

end subroutine f90PolyLog

!ccccccccccccccc

subroutine f90DiLog(z, res)
  use constants, only: dp; implicit none
  real (dp), intent(in ) :: z
  real (dp), intent(out) :: res
  real (dp)              :: dilog

  res = dilog(z)

end subroutine f90DiLog

!ccccccccccccccc

subroutine f90ComplexPolyLog(n, zre, zim, res)
  use Poly; use constants, only: dp; implicit none
  integer  , intent(in ) :: n
  real (dp), intent(in ) :: zre, zim
  real (dp), intent(out) :: res(2)
  complex (dp)           :: pl

  pl  = ComplexPolyLog( n, cmplx(zre, zim, dp) )
  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90ComplexPolyLog

!ccccccccccccccc

subroutine f90NGLSoft(nf, zre, zim, res)
  use MatrixElementsClass;  use constants, only: dp; implicit none
  integer  , intent(in )     :: nf
  real (dp), intent(in )     :: zre, zim
  real (dp), intent(out)     :: res(2)
  complex (dp)               :: pl

  pl = NGLSoft( nf, cmplx(zre, zim, dp) );  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90NGLSoft

!ccccccccccccccc

subroutine f90CLi2(zre, zim, res)
  use Chaplin;  use constants, only: dp; implicit none
  real (dp), intent(in ) :: zre, zim
  real (dp), intent(out) :: res(2)
  complex (dp)           :: pl

  pl = Cli2( cmplx(zre, zim, dp) );  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90CLi2

!ccccccccccccccc

subroutine f90CLi3(zre, zim, res)
  use Chaplin;  use constants, only: dp; implicit none
  real (dp), intent(in ) :: zre, zim
  real (dp), intent(out) :: res(2)
  complex (dp)           :: pl

  pl = Cli3( cmplx(zre, zim, dp) );  res = [ RealPart(pl), ImagPart(pl) ]

end subroutine f90CLi3

!ccccccccccccccc

subroutine f90Model(c, clen, lambda, k, l, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, l, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%ShapeFun(k, l)

end subroutine f90Model

!ccccccccccccccc

subroutine f90ModelDiff(c, clen, lambda, k, l, l2, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer                   , intent(in)  :: k, clen
  real (dp)                 , intent(in)  :: lambda, l, l2
  real (dp), dimension(clen), intent(in)  :: c
  real (dp)                 , intent(out) :: res
  type (Model)                            :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%ShapeFun(k, l, l2)

end subroutine f90ModelDiff

!ccccccccccccccc

subroutine f90BreitModel(c, clen, lambda, width, k, l, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, l, width, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%BreitModel(width, k, l)

end subroutine f90BreitModel

!ccccccccccccccc

subroutine f90BreitModelDiff(c, clen, lambda, width, k, l, l2, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, l, l2, width, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%BreitModel(width, k, l, l2)

end subroutine f90BreitModelDiff

!ccccccccccccccc

subroutine f90Taylor(c, clen, lambda, k, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%Taylor(k)

end subroutine f90Taylor

!ccccccccccccccc

subroutine f90TaylorPiece(c, lambda, k, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in) :: k, c(2)
  real (dp), intent(in) :: lambda
  real(dp), intent(out) :: res
  type (Model)          :: Mod

  Mod = Model(lambda, [0._dp, 0._dp], c, 'piece');  res = Mod%Taylor(k)

end subroutine f90TaylorPiece

!ccccccccccccccc

subroutine f90MomentModel(c, clen, lambda, k, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in)  :: k, clen
  real (dp), intent(in)  :: lambda, c(clen)
  real (dp), intent(out) :: res
  type (Model)           :: Mod

  Mod = Model(lambda, c, [0, 0], 'sum');  res = Mod%MomentModel(k)

end subroutine f90MomentModel

!ccccccccccccccc

subroutine f90ModelPiece(c, lambda, k, l, res)
  use ModelClass;  use constants, only: dp; implicit none
  integer  , intent(in) :: k, c(2)
  real (dp), intent(in) :: lambda, l
  real(dp), intent(out) :: res
  type (Model)          :: Mod

  Mod = Model(lambda, [0._dp, 0._dp], c, 'piece');  res = Mod%ShapeFun(k, l)

end subroutine f90ModelPiece

!ccccccccccccccc

subroutine f90MSbarDeltaPiece(nl, nh, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                    , intent(in )  :: nl, nh
  real (dp), dimension(0:4,4), intent(out)  :: beta

  beta = MSbarDeltaPiece(nl, nh)

end subroutine f90MSbarDeltaPiece

!ccccccccccccccc

subroutine f90AlphaMatchingLog(str, direction, nf, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)          , intent(in )  :: str, direction
  integer                      , intent(in )  :: nf
  real (dp), dimension(0:4,5), intent(out)    :: beta
  type (AnomDim)                              :: run

  run = AnomDim(str, nf, 0._dp)
  beta = run%MatchingAlphaLog(direction)

end subroutine f90AlphaMatchingLog

!ccccccccccccccc

subroutine f90AlphaMatchingInverse(str, nf, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)      , intent(in )  :: str
  integer                  , intent(in )  :: nf
  real (dp), dimension(5), intent(out)    :: beta
  type (AnomDim)                          :: run

  run = AnomDim(str, nf, 0._dp)
  beta = getInverse( run%Matchingalpha() )

end subroutine f90AlphaMatchingInverse

!ccccccccccccccc

subroutine f90AnomDim(str, nf, G4, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in )  :: str
  integer                , intent(in )  :: nf
  real (dp)              , intent(in )  :: G4
  real (dp), dimension(5), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, G4);  beta = run%betaQCD(str)

end subroutine f90AnomDim

!ccccccccccccccc

subroutine f90AfromS(str, nf, lambda, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in )  :: str
  real (dp)              , intent(in )  :: lambda
  integer                , intent(in )  :: nf
  real (dp), dimension(4), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%aFromS(str, lambda)

end subroutine f90AfromS

!ccccccccccccccc

subroutine f90aNasy(str, nf, order, lambda, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in )  :: str
  real (dp)              , intent(in )  :: lambda
  integer                , intent(in )  :: nf, order
  real (dp), dimension(4), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%aNasy(str, order, lambda)

end subroutine f90aNasy

!ccccccccccccccc

subroutine f90Ql(str, nf, order, lambda, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in )  :: str
  real (dp)              , intent(in )  :: lambda
  integer                , intent(in )  :: nf, order
  real (dp), dimension(4), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%Ql(str, order, lambda)

end subroutine f90Ql

!ccccccccccccccc

subroutine f90anLambda(str, nf, lambda, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in )  :: str
  real (dp)              , intent(in )  :: lambda
  integer                , intent(in )  :: nf
  real (dp), dimension(4), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%anLambda(str, lambda)

end subroutine f90anLambda

!ccccccccccccccc

subroutine f90cCoef(nf, order, n, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                  , intent(in ) :: nf, order, n
  real (dp), dimension(0:n), intent(out) :: beta
  type (AnomDim)                         :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%cCoeff(order, n)

end subroutine f90cCoef

!ccccccccccccccc

subroutine f90PScoef(nf, lg, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                , intent(in )  :: nf
  real (dp)              , intent(in )  :: lg
  real (dp), dimension(4), intent(out)  :: beta
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%PScoef(lg)

end subroutine f90PScoef

!ccccccccccccccc

subroutine f90N12(str, order, nf, lambda, err, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: nf, order
  real (dp)          , intent(in ) :: lambda, err
  real (dp)          , intent(out) :: res
  type (AnomDim)                   :: run

  run = AnomDim('MSbar', nf, 0._dp, err);  res = run%N12SumRule(order, str, lambda)

end subroutine f90N12

!ccccccccccccccc

subroutine f90N12Ratio(str, order, nf, lambda, err, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: nf, order
  real (dp)          , intent(in ) :: lambda, err
  real (dp)          , intent(out) :: res
  type (AnomDim)                   :: run

  run = AnomDim('MSbar', nf, 0._dp, err);  res = run%N12Ratio(order, str, lambda)

end subroutine f90N12Ratio

!ccccccccccccccc

subroutine f90N12Residue(str, order, nf, lambda, err, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: nf, order
  real (dp)          , intent(in ) :: lambda, err
  real (dp)          , intent(out) :: res
  type (AnomDim)                   :: run

  run = AnomDim('MSbar', nf, 0._dp, err);  res = run%N12Residue(order, str, lambda)

end subroutine f90N12Residue

!ccccccccccccccc

subroutine f90N12RS(str, order, n, nf, lambda, err, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: nf, n, order
  real (dp)          , intent(in ) :: lambda, err
  real (dp)          , intent(out) :: res
  type (AnomDim)                   :: run

  run = AnomDim('MSbar', nf, 0._dp, err);  res = run%N12RS(order, n, str, lambda)

end subroutine f90N12RS

!ccccccccccccccc

subroutine f90N12Generic(a, order, nf, lambda, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  real (dp), dimension(4), intent(in ) :: a
  integer                , intent(in ) :: nf, order
  real (dp)              , intent(in ) :: lambda
  real (dp)              , intent(out) :: res
  type (AnomDim)                       :: run

  run = AnomDim('MSbar', nf, 0._dp, 0._dp);  res = run%N12Generic(order, a, lambda)

end subroutine f90N12Generic

!ccccccccccccccc

subroutine f90Scoef(str, nf, beta)
  use AlphaClass;  use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass; implicit none

  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp), dimension(3), intent(out) :: beta
  real (dp)                            :: mu = 20._dp
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, &
  0._dp, 0._dp, 0._dp, 0._dp, 0._dp)

  MatEl     = MatrixElementsMass(alphaAll, nf, 0, 0._dp, 0._dp, 0._dp, mu, &
  0._dp, mu, mu, tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)

  beta = MatEl%Scoef(str)

end subroutine f90Scoef

!ccccccccccccccc

subroutine f90sCoefLambda(str, nf, lambda, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp)              , intent(in ) :: lambda
  real (dp), dimension(4), intent(out) :: beta
  type (AnomDim)                       :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%sCoefLambda( str, lambda )

end subroutine f90sCoefLambda

!ccccccccccccccc

subroutine f90sCoefGamma(gama, n, nf, beta)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                , intent(in ) :: nf, n
  real (dp), dimension(n), intent(in ) :: gama
  real (dp), dimension(n), intent(out) :: beta
  type (AnomDim)                       :: run

  run = AnomDim('MSbar', nf, 0._dp);  beta = run%sCoef( gama )

end subroutine f90sCoefGamma

!ccccccccccccccc

subroutine f90wTilde(order, nf, gamma, a0, a1, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer  , intent(in )                 :: order, nf
  real (dp), intent(in )                 :: a0, a1
  real (dp), intent(in ), dimension(0:3) :: gamma
  real (dp), intent(out)                 :: res
  type (AnomDim)                         :: run

  run = AnomDim('MSbar', nf, 0._dp);  res = run%wTilde(order, gamma, a0, a1)

end subroutine f90wTilde

!ccccccccccccccc

subroutine f90kTilde(order, nf, gamma, a0, a1, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  integer                 , intent(in) :: order, nf
  real (dp)               , intent(in ) :: a0, a1
  real (dp)               , intent(out) :: res
  real (dp), intent(in), dimension(0:3) :: gamma
  type (AnomDim)                        :: run

  run = AnomDim('MSbar', nf, 0._dp);  res = run%kTilde(order, gamma, a0, a1)

end subroutine f90kTilde

!ccccccccccccccc

subroutine f90InteCorre(b, x0, x1, res)
  use AnomDimClass;  use constants, only: dp; implicit none
  real (dp), intent(in ) :: b, x0, x1
  real (dp), intent(out) :: res

  res = InteCorre(b, x0, x1)

end subroutine f90InteCorre

!ccccccccccccccc

subroutine f90alphaQCD(str, method, order, run, nf, mZ, amZ, mT, muT, mB, muB, &
mC, muC, mu, res)
  use AlphaClass;  use constants, only: dp; use AnomDimClass; implicit none
  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaOb
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaOb = alpha( Andim, order, run, mZ, amZ, mT, muT, mB, &
  muB, mC, muC, method(:9), 0._dp )

  res     = alphaOb%alphaQCD(nf, mu)

end subroutine f90alphaQCD

!ccccccccccccccc

subroutine f90alphaQED(nf, mZ, amZQED, mT, muT, mB, muB, mC, muC, mu, res)
  use AlphaClass;  use constants, only: dp; use AnomDimClass; implicit none
  integer            , intent(in ) :: nf
  real (dp)          , intent(in ) :: mZ, amZQED, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaOb
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim( 'MSbar', i, 0._dp)
  end do

  alphaOb = alpha( Andim, 0, 0, mZ, amZQED, mT, muT, mB, muB, mC, muC, &
  'analytic', amZQED )

  res = alphaOb%alphaQED(nf, mu)

end subroutine f90alphaQED

!ccccccccccccccc

subroutine f90alphaComplex(str, method, order, run, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, muR, muI, res)
  use AlphaClass;  use constants, only: dp; use AnomDimClass; implicit none
  character (len = *)    , intent(in ) :: str, method
  integer                , intent(in ) :: order, run, nf
  real (dp)              , intent(in ) :: mZ, amZ, muR, muI, mT, muT, mB, muB, mC, muC
  real (dp), dimension(2), intent(out) :: res
  type (Alpha)                         :: alphaOb
  type (AnomDim), dimension(3:6)       :: AnDim
  complex (dp)                         :: alp
  integer                              :: i

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaOb = alpha( AnDim, order, run, mZ, amZ, mT, muT, mB, &
  muB, mC, muC, method(:9) )

  alp = alphaOb%alphaQCD( nf, cmplx(muR, muI, dp) )

  res = [ RealPart(alp), ImagPart(alp) ]

end subroutine f90alphaComplex

!ccccccccccccccc

subroutine f90MSbarMass(orderAlp, runAlp, run, nf, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass; implicit none

  integer          , intent(in ) :: orderAlp, runAlp, run, nf
  real (dp)        , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)        , intent(out) :: res
  type (Running)                 :: alphaMass
  type (Alpha)                   :: alphaAll
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, muC);  res = alphaMass%MSbarMass(mu)

end subroutine f90MSbarMass

!ccccccccccccccc

subroutine f90PoleMass(orderAlp, runAlp, order, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer            , intent(in ) :: orderAlp, runAlp, run, order, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, muC);  res = alphaMass%PoleMass(order, mu)

end subroutine f90PoleMass

!ccccccccccccccc

subroutine f90MSbarMassLow(orderAlp, runAlp, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer          , intent(in ) :: orderAlp, runAlp, run, nf
  real (dp)        , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)        , intent(out) :: res
  type (Running)                 :: alphaMass
  type (Alpha)                   :: alphaAll
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, muC); res = alphaMass%MSbarMassLow(mu)

end subroutine f90MSbarMassLow

!ccccccccccccccc

subroutine f90OptimalR(type, n, method, orderAlp, runAlp, order, run, nf, &
mZ, amZ, mT, muT, mB, muB, mC, muC, lambda, mu, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in ) :: method, type
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, lambda, n
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, mu)
  res = alphaMass%OptimalR( type, n, order, lambda, method(:8) )

end subroutine f90OptimalR


!ccccccccccccccc

subroutine f90OptimalR2(n, orderAlp, runAlp, order, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mass, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer            , intent(in ) :: orderAlp, runAlp, order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, muC, n, mass
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, 100._dp)
  res = alphaMass%OptimalR2( n, mass )

end subroutine f90OptimalR2

!ccccccccccccccc

subroutine f90UpsilonDeltaCharm(n, l, alp, mb, mc, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none
  real (dp)        , intent(in ) :: alp, mb, mc
  integer          , intent(in ) :: n, l
  real (dp)        , intent(out) :: res
  type (NRQCD)                   :: Upsilon
  type (VFNSMSR)                 :: MSR
  type (Alpha)                   :: alphaAll
  type (Running), dimension(2)   :: alphaMass
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, &
  0._dp, 0._dp, mc, mc, mc, mc, 'analytic', 0._dp)

  alphaMass = Running(4, 0, alphaAll, 100._dp)
  MSR       = VFNSMSR(alphaMass)
  Upsilon   = NRQCD( 'up', 'MSbar', 'no', MSR, n, l, 0, 0 )

  res = Upsilon%DeltaCharm(alp, mb)

end subroutine f90UpsilonDeltaCharm

!ccccccccccccccc

subroutine f90UpsilonDeltaCharmBin(n, l, alp, mb, mc, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none
  real (dp)        , intent(in ) :: alp, mb, mc
  integer          , intent(in ) :: n, l
  real (dp)        , intent(out) :: res
  type (NRQCD)                   :: Upsilon
  type (VFNSMSR)                 :: MSR
  type (Alpha)                   :: alphaAll
  type (Running), dimension(2)   :: alphaMass
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, &
  0._dp, 0._dp, mc, mc, mc, mc, 'analytic', 0._dp)

  alphaMass = Running(4, 0, alphaAll, 100._dp)
  MSR       = VFNSMSR(alphaMass)
  Upsilon   = NRQCD( 'up', 'MSbar', 'no', MSR, n, l, 0, 0 )

  res = Upsilon%DeltaCharmBin(alp, mb)

end subroutine f90UpsilonDeltaCharmBin

!ccccccccccccccc

subroutine f90DeltaCharmExact(charm, type, scheme, average, n, l, j, s, nl, &
  mH, mL, mu, alp, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: type, scheme, average, charm
  integer            , intent(in ) :: nl, n, l, j, s
  real (dp)          , intent(in ) :: mH, mL, alp, mu
  real (dp)          , intent(out) :: res
  real (dp)                        :: mB, mC, mT
  type (NRQCD)                     :: Upsilon
  type (Alpha)                     :: alphaAll
  type (VFNSMSR)                   :: MSR
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  if (nl == 5) then
    mT = mH; mB = mL
  else if (nl == 4) then
    mB = mH; mC = mL
  else
    res = 0; return
  end if

  do i = 3, 6
    AnDim(i) = AnomDim(scheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 4, 4, 91.187_dp, 0.118_dp, &
  mT, mT, mB, mB, mC, mC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, 0, alphaAll, mL), &
  Running(nl, 0, alphaAll, mH) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm, scheme(:4), average(:3), MSR, n, l, j, s )

  res = Upsilon%DeltaCharmExact(type(:5), mu, alp, mH)

end subroutine f90DeltaCharmExact

!ccccccccccccccc

subroutine f90NRQCD(n, l, j, s, charm, scheme, average, method, counting, &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
lambda1, lambda2, lam, mu, R, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, charm, average, counting
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nl, n,&
  l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, &
  lambda1, lambda2, lam, R
  real (dp)          , intent(out) :: res(5)
  character (len = 5)              :: alphaScheme
  type (NRQCD)                     :: Upsilon
  type (Alpha)                     :: alphaAll
  type (VFNSMSR)                   :: MSR
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%En( order, mu, R, lam, method(:8), counting(:5) )

end subroutine f90NRQCD

!ccccccccccccccc

subroutine f90NRQCDDerCharm(n, l, j, s, charm, scheme, average, method, counting, &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
lambda1, lambda2, lam, mu, R, eps, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, charm, average, counting
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nl, n,&
  l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, &
  lambda1, lambda2, lam, R, eps
  real (dp)          , intent(out) :: res(5)
  character (len = 5)              :: alphaScheme
  type (NRQCD)                     :: Upsilon
  type (Alpha)                     :: alphaAll
  type (VFNSMSR)                   :: MSR
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%EnDerCharm( order, mu, R, lam, method(:8), counting(:5), eps )

end subroutine f90NRQCDDerCharm

!ccccccccccccccc

subroutine f90NRQCDDerAlpha(n, l, j, s, charm, scheme, average, method, counting, &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
lambda1, lambda2, lam, mu, R, eps, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, charm, average, counting
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nl, n,&
  l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, &
  lambda1, lambda2, lam, R, eps
  real (dp)          , intent(out) :: res(5)
  character (len = 5)              :: alphaScheme
  type (NRQCD)                     :: Upsilon
  type (Alpha)                     :: alphaAll
  type (VFNSMSR)                   :: MSR
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%EnDerAlpha( order, mu, R, lam, method(:8), counting(:5), eps )

end subroutine f90NRQCDDerAlpha

!ccccccccccccccc

subroutine f90FindMass(ord, n, l, j, s, iter, charm, scheme, average, method, &
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, &
muC, mass, lambda1, lambda2, lam, mu, R, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, charm, iter, average, counting
  integer            , intent(in ) :: ord, orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
  lambda1, lam, R, mass, lambda2
  real (dp)          , intent(out) :: res
  character (len = 5)              :: alphaScheme
  type (NRQCD)                     :: Upsilon
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%MassFitter( iter(:10), ord, order, mu, R, mass, &
  lam, method(:8), counting(:5) )

end subroutine f90FindMass

!ccccccccccccccc

subroutine f90FindEnergy(ord, n, l, j, s, iter, charm, scheme, average, method, &
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, &
muC, lambda1, lambda2, lam, mu, R, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, charm, iter, average, counting
  integer            , intent(in ) :: ord, orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
  lambda1, lam, R, lambda2
  real (dp)          , intent(out) :: res
  character (len = 5)              :: alphaScheme
  type (NRQCD)                     :: Upsilon
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%EnFitter( iter(:10), ord, order, mu, R, lam, method(:8), counting(:5) )

end subroutine f90FindEnergy

!ccccccccccccccc

subroutine f90MassError(ord, n, l, j, s, iter, charm, scheme, average, method,&
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC,&
mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, x, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *)    , intent(in ) :: method, scheme, charm, average, iter, counting
  integer                , intent(in ) :: ord, orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, x, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2, muC
  real (dp), dimension(2), intent(out) :: res
  character (len = 5)                  :: alphaScheme
  type (NRQCD)                         :: Upsilon
  type (VFNSMSR)                       :: MSR
  type (Alpha)                         :: alphaAll
  type (Running), dimension(2)         :: alphaMass
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%MassError( iter(:10), ord, order, mu0, mu1, &
  deltaMu, R0, R1, deltaR, x, mass, lam, method(:8), counting(:5) )

end subroutine f90MassError

!ccccccccccccccc

subroutine f90MassList(ord, n, l, j, s, iter, charm, scheme, average, method,  &
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *)    , intent(in ) :: method, scheme, charm, average, iter, counting
  integer                , intent(in ) :: ord, orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2, muC

  real (dp), dimension( 3, 0:Floor( (mu1 - mu0)/deltaMu ), &
  0:Floor( (R1 - R0)/deltaR ) ), intent(out) :: res

  character (len = 5)                  :: alphaScheme
  type (NRQCD)                         :: Upsilon
  type (VFNSMSR)                       :: MSR
  type (Alpha)                         :: alphaAll
  type (Running), dimension(2)         :: alphaMass
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%MassList( iter(:10), ord, order, mu0, mu1, &
  deltaMu, R0, R1, deltaR, mass, lam, method(:8), counting(:5) )

end subroutine f90MassList

!ccccccccccccccc

subroutine f90NRQCDList(n, l, j, s, iter, charm, scheme, average, method, &
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *)    , intent(in ) :: method, scheme, charm, average, iter, counting
  integer                , intent(in ) :: orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2, muC

  real (dp), dimension( 7, 0:Floor( (mu1 - mu0)/deltaMu ), &
  0:Floor( (R1 - R0)/deltaR ) ), intent(out) :: res

  character (len = 5)                  :: alphaScheme
  type (NRQCD)                         :: Upsilon
  type (VFNSMSR)                       :: MSR
  type (Alpha)                         :: alphaAll
  type (Running), dimension(2)         :: alphaMass
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%NRQCDList( iter(:10), order, mu0, mu1, &
  deltaMu, R0, R1, deltaR, mass, lam, method(:8), counting(:5) )

end subroutine f90NRQCDList

!ccccccccccccccc

subroutine f90UpsilonList(n, l, j, s, charm, scheme, average, method, &
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm, &
res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *)    , intent(in ) :: method, scheme, charm, average, counting
  integer                , intent(in ) :: orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)              , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, lambda2, epsCharm,   &
  epsAlpha

  real (dp), dimension( 5, 3, 0:Floor( (mu1 - mu0)/deltaMu ), &
  0:Floor( (R1 - R0)/deltaR ) ), intent(out) :: res

  character (len = 5)                  :: alphaScheme
  type (NRQCD)                         :: Upsilon
  type (VFNSMSR)                       :: MSR
  type (Alpha)                         :: alphaAll
  type (Running), dimension(2)         :: alphaMass
  type (AnomDim), dimension(3:6)       :: AnDim
  integer                              :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%UpsilonList( order, mu0, mu1, deltaMu, R0, R1, deltaR, &
  lam, method(:8), counting(:5), epsAlpha, epsCharm )

end subroutine f90UpsilonList

!ccccccccccccccc

subroutine f90Chi2NRQCD(qnList, datalist, m, iter, charm, scheme, average, method,  &
counting, orderAlp, runAlp, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC, &
muC, lambda1, lambda2, lam, muList, RList, ndim, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  real (dp)               , intent(out) :: res
  character (len = *)      , intent(in) :: method, scheme, charm, average, counting, iter
  integer , dimension(4,m) , intent(in) :: qnList
  real(dp), dimension(2,m) , intent(in) :: dataList
  real(dp), dimension(ndim), intent(in) :: Rlist, muList
  integer                  , intent(in) :: orderAlp, runAlp, order, run, m, nl, n, ndim
  real (dp)                , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, lambda2

  character (len = 5)                   :: alphaScheme
  type (NRQCD), dimension(m)            :: Upsilon
  type (VFNSMSR)                        :: MSR
  type (Alpha)                          :: alphaAll
  type (Running), dimension(2)          :: alphaMass
  type (AnomDim), dimension(3:6)        :: AnDim
  integer                               :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  do i = 1, m
    Upsilon(i) = NRQCD( charm(:4), scheme(:5), average(:3), MSR, qnList(1,i), &
    qnList(2,i), qnList(3,i), qnList(4,i) )
  end do

  res = Chi2NRQCD(Upsilon, datalist, m, iter, order, n, muList, RList, ndim, &
  lam, method(:8), counting(:5))

end subroutine f90Chi2NRQCD

!ccccccccccccccc

subroutine f90Chi2MinNRQCD(qnList, datalist, m, iter, charm, scheme, average, method,&
counting, orderAlp, runAlp, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC,  &
muC, lambda1, lambda2, lam, muList, RList, ndim, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  real (dp), dimension(3) , intent(out) :: res
  character (len = *)      , intent(in) :: method, scheme, charm, average, counting, iter
  integer , dimension(4,m) , intent(in) :: qnList
  real(dp), dimension(2,m) , intent(in) :: dataList
  integer                  , intent(in) :: orderAlp, runAlp, order, run, m, nl, n, ndim
  real(dp), dimension(ndim), intent(in) :: Rlist, muList
  real (dp)                , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, lambda2

  character (len = 5)                   :: alphaScheme
  type (NRQCD), dimension(m)            :: Upsilon
  type (VFNSMSR)                        :: MSR
  type (Alpha)                          :: alphaAll
  type (Running), dimension(2)          :: alphaMass
  type (AnomDim), dimension(3:6)        :: AnDim
  integer                               :: i
  real (dp)                             :: mH

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  do i = 1, m
    Upsilon(i) = NRQCD( charm(:4), scheme(:5), average(:3), MSR, qnList(1,i), &
    qnList(2,i), qnList(3,i), qnList(4,i) )
  end do

  if (nl == 5) mH = mT;  if (nl == 4) mH = mB;  if (nl == 3) mH = mC

  res = Chi2MinNRQCD(Upsilon, datalist, m, iter, order, n, mH, muList, RList, ndim, &
  lam, method(:8), counting(:5))

end subroutine f90Chi2MinNRQCD

!ccccccccccccccc

subroutine f90Chi2MinAlphaNRQCD(qnList, datalist, m, iter, charm, scheme, average, method,&
counting, orderAlp, runAlp, order, run, n, nl, mZ, amZ, mT, muT, mB, muB, mC,  &
muC, lambda1, lambda2, lam, muList, RList, ndim, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  real (dp), dimension(3) , intent(out) :: res
  character (len = *)      , intent(in) :: method, scheme, charm, average, counting, iter
  integer , dimension(4,m) , intent(in) :: qnList
  real(dp), dimension(2,m) , intent(in) :: dataList
  integer                  , intent(in) :: orderAlp, runAlp, order, run, m, nl, n, ndim
  real(dp), dimension(ndim), intent(in) :: Rlist, muList
  real (dp)                , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, lambda2

  character (len = 5)                   :: alphaScheme
  type (NRQCD), dimension(m)            :: Upsilon
  type (VFNSMSR)                        :: MSR
  type (Alpha)                          :: alphaAll
  type (Running), dimension(2)          :: alphaMass
  type (AnomDim), dimension(3:6)        :: AnDim
  integer                               :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  do i = 1, m
    Upsilon(i) = NRQCD( charm(:4), scheme(:5), average(:3), MSR, qnList(1,i), &
    qnList(2,i), qnList(3,i), qnList(4,i) )
  end do

  res = Chi2MinAlphaNRQCD( Upsilon, datalist, m, iter, order, n, amZ, muList, RList, &
  ndim, lam, method(:8), counting(:5) )

end subroutine f90Chi2MinAlphaNRQCD

!ccccccccccccccc

subroutine f90Chi2MinAlphaMbNRQCD(qnList, datalist, m, iter, charm, scheme, average, &
method, counting, orderAlp, runAlp, order, run, n, nl, mZ, amZ, mT, muT, mB, &
muB, mC, muC, lambda1, lambda2, lam, muList, RList, ndim, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  real (dp), dimension(6) , intent(out) :: res
  character (len = *)      , intent(in) :: method, scheme, charm, average, counting, iter
  integer , dimension(4,m) , intent(in) :: qnList
  real(dp), dimension(2,m) , intent(in) :: dataList
  integer                  , intent(in) :: orderAlp, runAlp, order, run, m, nl, n, ndim
  real(dp), dimension(ndim), intent(in) :: Rlist, muList
  real (dp)                , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, lambda2

  character (len = 5)                   :: alphaScheme
  type (NRQCD), dimension(m)            :: Upsilon
  type (VFNSMSR)                        :: MSR
  type (Alpha)                          :: alphaAll
  type (Running), dimension(2)          :: alphaMass
  type (AnomDim), dimension(3:6)        :: AnDim
  integer                               :: i
  real (dp)                             :: mH

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  do i = 1, m
    Upsilon(i) = NRQCD( charm(:4), scheme(:5), average(:3), MSR, qnList(1,i), &
    qnList(2,i), qnList(3,i), qnList(4,i) )
  end do

  if (nl == 5) mH = mT;  if (nl == 4) mH = mB;  if (nl == 3) mH = mC

  res = Chi2MinAlphaMbNRQCD( Upsilon, datalist, m, iter, order, n, mH, amZ, muList, &
  RList, ndim, lam, method(:8), counting(:5) )

end subroutine f90Chi2MinAlphaMbNRQCD

!ccccccccccccccc

subroutine f90CorrMat(qnList, m, charm, scheme, average, method, counting,     &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, &
lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm, massList, &
corMat)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *)     , intent(in) :: method, scheme, charm, average, counting
  integer, dimension(4,m) , intent(in) :: qnList
  integer                 , intent(in) :: orderAlp, runAlp, order, run, m, nl
  real (dp)               , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, lambda2, epsCharm,   &
  epsAlpha

  real (dp), dimension( m, 5, 5 ), intent(out) :: massList
  real (dp), dimension( m, m, 5 ), intent(out) :: corMat
  character (len = 5)                          :: alphaScheme
  type (NRQCD), dimension(m)                   :: Upsilon
  type (VFNSMSR)                               :: MSR
  type (Alpha)                                 :: alphaAll
  type (Running), dimension(2)                 :: alphaMass
  type (AnomDim), dimension(3:6)               :: AnDim
  integer                                      :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  do i = 1, m
    Upsilon(i) = NRQCD( charm(:4), scheme(:5), average(:3), MSR, qnList(1,i), &
    qnList(2,i), qnList(3,i), qnList(4,i) )
  end do

  call CorrMattNRQCD( Upsilon, m, order, mu0, mu1, deltaMu, R0, R1, deltaR, &
  lam, method(:8), counting(:5), epsAlpha, epsCharm, massList, CorMat )

end subroutine f90CorrMat

!ccccccccccccccc

subroutine f90ErrMat(qnList, m, charm, scheme, average, method, counting,     &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, &
lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm, massList, &
corMat)

  use constants, only: dp; use NRQCDClass;  implicit none

  character (len = *)     , intent(in) :: method, scheme, charm, average, counting
  integer, dimension(4,m) , intent(in) :: qnList
  integer                 , intent(in) :: orderAlp, runAlp, order, run, m, nl
  real (dp)               , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, lambda2, epsCharm,   &
  epsAlpha

  real (dp), dimension( m, 4, 5 ), intent(out) :: massList
  real (dp), dimension( m, m, 5 ), intent(out) :: corMat
  real (dp), dimension( m, 5, 5 )              :: massList2
  integer                                      :: i, j

  call f90CorrMat(qnList, m, charm, scheme, average, method, counting,     &
  orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, &
  lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm, massList2,&
  CorMat)

  massList(:,:2,:) = massList2(:,:2,:)
  massList(: ,3,:) = epsAlpha * massList2(:,4,:)
  massList(: ,4,:) = epsCharm * massList2(:,5,:)

  do i = 1, m
    do j = 1, m
      corMat(i,j,:) = massList2(i,3,:) * massList2(j,3,:) * corMat(i,j,:)
    end do
  end do

end subroutine f90ErrMat

!ccccccccccccccc

subroutine f90ErrMatrices(qnList, m, charm, scheme, average, method, counting, &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, &
lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm, massList, &
corMat)

  use constants, only: dp; use NRQCDClass;  implicit none

  character (len = *)     , intent(in) :: method, scheme, charm, average, counting
  integer, dimension(4,m) , intent(in) :: qnList
  integer                 , intent(in) :: orderAlp, runAlp, order, run, m, nl
  real (dp)               , intent(in) :: mZ, amZ, mT, muT, mB, muB, mC, muC, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, lambda2, epsCharm,   &
  epsAlpha

  real (dp), dimension( m, 2, 5    ), intent(out) :: massList
  real (dp), dimension( m, m, 3, 5 ), intent(out) :: corMat
  real (dp), dimension( m, 4, 5 )                 :: massList2
  integer                                         :: i, j

  call f90ErrMat( qnList, m, charm, scheme, average, method, counting, orderAlp,&
  runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, &
  lam, mu0, mu1, deltaMu, R0, R1, deltaR, epsAlpha, epsCharm, massList2,        &
  CorMat(:,:,1,:) )

  massList(:,:2,:) = massList2(:,:2,:)

  do i = 1, m
    do j = 1, m
      corMat(i,j,2:3,:) = massList2(i,3:4,:) * massList2(j,3:4,:)
    end do
  end do

end subroutine f90ErrMatrices

!ccccccccccccccc

subroutine f90NRQCDError(n, l, j, s, iter, charm, scheme, average, method, &
counting, orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, &
mass, lambda1, lambda2, lam, mu0, mu1, deltaMu, R0, R1, deltaR, x, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, charm, average, iter, counting
  integer            , intent(in ) :: orderAlp, runAlp, order, run, &
  nl, n, l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, muC, x, &
  lambda1, lam, mu0, mu1, deltaMu, R0, R1, deltaR, mass, lambda2
  real (dp), dimension(2,0:4), intent(out) :: res
  character (len = 5)              :: alphaScheme
  type (NRQCD)                     :: Upsilon
  type (VFNSMSR)                   :: MSR
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%EnError( iter(:10), order, mu0, mu1, &
  deltaMu, R0, R1, deltaR, x, mass, lam, method(:8), counting(:5) )

end subroutine f90NRQCDError

!ccccccccccccccc

subroutine f90MassIter(n, l, j, s, charm, scheme, average, method, orderAlp, &
counting, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass,  &
lambda1, lambda2, lam, mu, R, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, average, charm, counting
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nl, &
  n, l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
  lambda1, lambda2, lam, R, mass
  real (dp)          , intent(out) :: res(5)
  character (len = 5)              :: alphaScheme
  type (VFNSMSR)                   :: MSR
  type (NRQCD)                     :: Upsilon
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%MassIter( order, mu, R, mass, lam, method(:8), counting(:5) )

end subroutine f90MassIter

!ccccccccccccccc

subroutine f90MassExpand(n, l, j, s, charm, scheme, average, method, counting, &
orderAlp, runAlp, order, run, nl, mZ, amZ, mT, muT, mB, muB, mC, muC, mass,    &
lambda1, lambda2, lam, mu, R, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = *), intent(in ) :: method, scheme, average, charm, counting
  integer            , intent(in ) :: orderAlp, runAlp, order, run, nl, &
  n, l, j, s
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC,&
  lambda1, lambda2, lam, R, mass
  real (dp)          , intent(out) :: res(5)
  character (len = 5)              :: alphaScheme
  type (VFNSMSR)                   :: MSR
  type (NRQCD)                     :: Upsilon
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR     = VFNSMSR(alphaMass)
  Upsilon = NRQCD( charm(:4), scheme(:5), average(:3), MSR, n, l, j, s )

  res = Upsilon%EnExpand( order, mu, R, mass, lam, method(:8), counting(:5) )

end subroutine f90MassExpand

!ccccccccccccccc

subroutine f90mmfromMSR(type, orderAlp, runAlp, order, run, nf, mZ, amZ, &
mT, muT, mB, muB, mC, muC, mu, R, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  integer            , intent(in ) :: orderAlp, runAlp, order, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, R
  character (len = *), intent(in ) :: type
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (AnomDim), dimension(3:6)   :: AnDim
  type (Alpha)                     :: alphaAll
  real (dp)                        :: mass
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  if (nf == 5) mass = mT; if (nf == 4) mass = mB; if (nf == 3) mass = mC

  alphaMass = Running(nf, run, alphaAll, mu)
  res = alphaMass%mmFromMSR(type, mass, order, R)

end subroutine f90mmfromMSR

!ccccccccccccccc

subroutine f90MSRMass(type, method, orderAlp, runAlp, order, run, nf, mZ, &
amZ, mT, muT, mB, muB, mC, muC, lambda, mu, R, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in) :: method, type
  integer           , intent(in ) :: orderAlp, runAlp, order, run, nf
  real (dp)         , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, R, lambda
  real (dp)         , intent(out) :: res
  type (Running)                  :: alphaMass
  type (Alpha)                    :: alphaAll
  type (AnomDim), dimension(3:6)  :: AnDim
  integer                         :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, mu)
  res       = alphaMass%MSRMass( type, order, R, lambda, method(:8) )

end subroutine f90MSRMass

!ccccccccccccccc

subroutine f90RSMass(type, method, orderAlp, runAlp, order, run, nf, mZ, &
amZ, mT, muT, mB, muB, mC, muC, lambda, R, res)
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in) :: type, method
  integer           , intent(in ) :: orderAlp, runAlp, order, run, nf
  real (dp)         , intent(in ) :: mZ, amZ, mT, muT, mB, muB, mC, muC, R, lambda
  real (dp)         , intent(out) :: res
  type (Running)                  :: alphaMass
  type (Alpha)                    :: alphaAll
  type (AnomDim), dimension(3:6)  :: AnDim
  integer                         :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, run, alphaAll, 1._dp)
  res       = alphaMass%RSMass( type, order, R, lambda, method )

end subroutine f90RSMass

!ccccccccccccccc

subroutine f90MSRVFNS(up, type, method, orderAlp, runAlp, order, run, nf, mZ, &
amZ, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, R, res)
  use RunningClass;  use AlphaClass   ;  use constants, only: dp
  use AnomDimClass;  use VFNSMSRClass ;  implicit none

  character (len = *), intent(in) :: method, type, up
  integer           , intent(in ) :: orderAlp, runAlp, order, run, nf
  real (dp)         , intent(in ) :: mZ, amZ, mu1, mu2, mT, muT, mB, muB, mC, &
  muC, R, lambda
  real (dp)         , intent(out) :: res
  type (Running), dimension(2)    :: alphaMass
  type (AnomDim), dimension(3:6)  :: AnDim
  type (Alpha)                    :: alphaAll
  type (VFNSMSR)                  :: MSRVFNS
  integer                         :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass(1) = Running(nf - 1, run, alphaAll, mu2)
  alphaMass(2) = Running(nf    , run, alphaAll, mu1)
  MSRVFNS      = VFNSMSR( alphaMass )
  res          = MSRVFNS%MSRMass( up(:4), type(:5), order, R, lambda, method(:8) )

end subroutine f90MSRVFNS

!ccccccccccccccc

subroutine f90MSRTop(up, type, method, orderAlp, runAlp, order, run, nf, mZ, &
amZ, mT, muT, mB, muB, mC, muC, lambda, mu1, mu2, mu3, R, res)
  use RunningClass;  use AlphaClass   ;  use constants, only: dp
  use AnomDimClass;  use VFNSMSRClass ;  implicit none

  character (len = *), intent(in) :: method, type, up
  integer            , intent(in) :: orderAlp, runAlp, order, run, nf
  real (dp)          , intent(in) :: mZ, amZ, mu1, mu2, mT, muT, mB, muB, mC, &
  muC, R, lambda, mu3
  real (dp)         , intent(out) :: res
  type (Running), dimension(3)    :: alphaMass
  type (AnomDim), dimension(3:6)  :: AnDim
  type (Alpha)                    :: alphaAll
  type (topMSR)                   :: MSRVFNS
  integer                         :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass(1) = Running(nf - 2, run, alphaAll, mu3)
  alphaMass(2) = Running(nf - 1, run, alphaAll, mu2)
  alphaMass(3) = Running(nf    , run, alphaAll, mu1)
  MSRVFNS      = topMSR( alphaMass )
  res          = MSRVFNS%MSRMass( up(:5), type(:5), order, R, lambda, method(:8) )

end subroutine f90MSRTop

!ccccccccccccccc

subroutine f90JetMass(orderAlp, runAlp, order, run, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, muLambda, R, mu, res)
 use AlphaClass;    use MatrixElementsClass;  use constants, only: dp
 use AnomDimClass;  implicit none

  integer  , intent(in )      :: orderAlp, order, runAlp, run, nf
  real (dp), intent(in )      :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, R
  real (dp), intent(out)      :: res
  type (Alpha)                :: alphaAll
  type (MatricesElementsMass) :: MatEl
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                     :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl     = MatricesElementsMass(alphaAll, nf, run, 0._dp, 0._dp, 0._dp, 0._dp,  &
  muLambda, muLambda)

  res = MatEl%JetMass(order, run, R, mu)

end subroutine f90JetMass

!ccccccccccccccc

subroutine f90mmFromJetMass(orderAlp, runAlp, order, run, nf, mZ, amZ, mT, &
muT, mB, muB, mC, muC, muLambda, R, mu, res)
 use AlphaClass;   use MatrixElementsClass;  use constants, only: dp
 use AnomDimClass; implicit none

  integer  , intent(in )         :: orderAlp, order, runAlp, run, nf
  real (dp), intent(in )         :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, R
  real (dp), intent(out)         :: res
  type (Alpha)                   :: alphaAll
  type (MatricesElementsMass)    :: MatEl
  real (dp)                      :: mass
  type (AnomDim), dimension(3:6) :: AnDim
  integer                        :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl     = MatricesElementsMass(alphaAll, nf, run, 0._dp, 0._dp, 0._dp, 0._dp,  &
  muLambda, muLambda)

  if (nf == 6) mass = mT; if (nf == 5) mass = mB; if (nf == 4) mass = mC

  res = MatEl%mmFromJetMass(order, run, R, mu, mass)

end subroutine f90mmFromJetMass

!ccccccccccccccc

subroutine f90Singular(hard, shape, setup, gap, space, cum, orderAlp, runAlp, order, &
 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, &
 mu, c, clen, lambda, R0, mu0, delta0, h, tau, res)
  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, gap, hard
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC, &
  mu,  Q, G3, lambda, tau, R0, mu0, delta0, h, s3, muC, j3
  real (dp)          , intent(inout) :: muH, muJ, muS, R
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  integer                            :: i
  type (AnomDim), dimension(3:6)     :: AnDim

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')

    call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
    call Sing%SetRunning(muJ, muS, R, mu)

  res       = Sing%SingleSing(Mod, setup(:15), gap(:12), space(:6), cum(:4), &
  order, R0, mu0, delta0, h, tau)

end subroutine f90Singular

!ccccccccccccccc

subroutine f90SingularDiff(hard, shape, setup, gap, space, cum, orderAlp, runAlp, order, run, &
 nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, &
 c, clen, lambda, R0, mu0, delta0, h, tau1, tau2, res)
  use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension (clen), intent(in) :: c
  character (len = *), intent(in)    :: shape, cum, setup, space, hard, gap
  integer            , intent(in)    :: orderAlp, order, runAlp, run, nf, clen
  real (dp)          , intent(in)    :: mZ, amZ, muLambda, mT, muT, mB, muB, mC,&
  mu,  Q, G3, lambda, tau1, tau2, R0, mu0, delta0, h, s3, j3, muC
  real (dp)          , intent(inout) :: muH, muJ, muS, R
  real (dp)          , intent(out)   :: res
  type (Alpha)                       :: alphaAll
  type (MatricesElements)            :: MatEl
  type (SingularScales)              :: Sing
  type (Model)                       :: Mod
  type (AnomDim), dimension(3:6)     :: AnDim
  integer                            :: i

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 3, 6
    AnDim  = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')

  call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
  call Sing%SetRunning(muJ, muS, R, mu)

  res       = Sing%SingleSing(Mod, setup(:15), gap(:12), space(:6), cum(:4), &
  order, R0, mu0, delta0, h, tau1, tau2)

end subroutine f90SingularDiff

!ccccccccccccccc

subroutine f90SingularHJM(hard, setup, gap, space, cum, orderAlp, runAlp, order, run, &
  isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
  muH, muJ, muS, R, mu, c, clen, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass; use MatrixElementsClass; use ModelClass; use SingularClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in) :: cum, setup, space, hard, gap
  integer            , intent(in) :: orderAlp, order, runAlp, run, nf, clen, isoft
  real (dp)      ,  intent(inout) :: muH, muJ, muS, R
  real (dp)          , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC,&
  muC, j3, Q, G3, lambda, tau, R0, mu0, delta0, h, s31, s32, s3

  real (dp)          ,   intent(out)  :: res
  type (Alpha)                        :: alphaAll
  type (MatrixElements)               :: MatEl
  type (SingularMassless)             :: Sing
  type (Model), dimension(clen, clen) :: Mod
  type (AnomDim), dimension(3:6)      :: AnDim
  integer                             :: i, j

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = Sing%HJMSing(c, Mod, 'sum', setup(:15), gap(:12), space(:3), cum(:4), order, &
                     isoft, s31, s32, R0, mu0, delta0, h, tau)

end subroutine f90SingularHJM

!ccccccccccccccc

subroutine f90SingularHJM1D(hard, gap, cum, orderAlp, runAlp, order, run, isoft, nf, &
 j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, &
 muS, R, mu, c, clen, lambda, R0, mu0, delta0, h, tau, res)
  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), intent(in), dimension (clen) :: c
  character (len = *), intent(in ) :: cum, hard, gap
  integer            , intent(in ) :: orderAlp, order, runAlp, run, nf, clen, isoft
  real (dp)          , intent(in ) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC , &
                                     s3, muH, Q, muJ, muS, R, G3, lambda, tau, R0, mu0 , &
                                     delta0, h, s31, s32, j3
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaAll
  type (MatrixElements)            :: MatEl
  type (SingularMassless)          :: Sing
  type (Model)                     :: Mod
  type (AnomDim), dimension(3:6)   :: AnDim
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6) )
  Mod       = Model(lambda, c, [0,0], 'sum')

  res       = Sing%HJMSing1D(Mod, gap(:12), cum(:4), order, isoft, s31, s32, R0, mu0, &
                             delta0, h, tau)

end subroutine f90SingularHJM1D

!ccccccccccccccc

subroutine f90SingularHJM1DPiece(hard, gap, cum, orderAlp, runAlp, order, run, isoft, nf, &
 j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS , &
 R, mu, c, lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in ) :: cum, hard, gap
  integer              , intent(in ) :: orderAlp, order, runAlp, run, nf, isoft
  integer, dimension(2), intent(in ) :: c
  real (dp)            , intent(in ) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC , &
                                       s3, muH, Q, muJ, muS, R, G3, lambda, tau, mu0, &
                                       delta0, h, s31, s32, j3, muC, R0
  real (dp)           , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatrixElements)             :: MatEl
  type (SingularMassless)           :: Sing
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  integer                           :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))
  Mod       = Model(lambda, [1._dp], c, 'piece')

  res       = Sing%HJMSing1D(Mod, gap(:12), cum(:4), order, isoft, s31, s32, R0, mu0, &
                             delta0, h, tau)

end subroutine f90SingularHJM1DPiece

!ccccccccccccccc

subroutine f90SingularDouble(hard, setup, gap, space, cum1, cum2, orderAlp, runAlp, order, &
run, isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,   &
 muH, muJ, muS, R, mu, c, clen, lambda, R0, mu0, delta0, h, rho1, rho2, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  real (dp), dimension(clen, clen), intent(in) :: c
  character (len = *), intent(in) :: hard, cum1, cum2, setup, space, gap
  integer            , intent(in) :: orderAlp, order, runAlp, run, nf, clen, isoft
  real (dp)          , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC , &
                                     j3, s3, h, Q, G3, lambda, rho1, rho2, R0, mu0, s31,&
                                     s32, delta0
  real (dp)       , intent(inout) :: muH, muJ, muS, R
  real (dp)       , intent(out)   :: res
  type (Alpha)                        :: alphaAll
  type (MatrixElements)               :: MatEl
  type (SingularMassless)             :: Sing
  type (AnomDim), dimension(3:6)      :: AnDim
  type (Model), dimension(clen, clen) :: Mod
  integer                             :: i, j

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  if ( setup(:2) == 'FO' ) then
    muS = mu; muJ = mu; muH = mu; R = mu
  end if

  do i = 1, clen
    do j = 1, i
      Mod(i,j) = Model(lambda, [1._dp], [i - 1,j - 1], 'piece')
    enddo
  end do

  do i = 1, clen
    do j = i + 1, clen
      Mod(i,j) = Mod(j,i)
    enddo
  end do

  res = Sing%DoubleSing(c, Mod, 'sum', setup(:15), gap(:12), space(:3), cum1(:4), &
                        cum2(:4), order, isoft, s31, s32, R0, mu0, delta0, h, rho1, rho2)

end subroutine f90SingularDouble

!ccccccccccccccc

subroutine f90SingularHJMPiece(hard, gap, space, cum, orderAlp, runAlp, order, run , &
 isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
 muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, tau, res)

 use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
 use constants, only: dp; use AnomDimClass; implicit none

 character (len = *)  , intent(in) :: cum, space, hard, gap
 integer              , intent(in) :: orderAlp, order, runAlp, run, nf, isoft
 integer, dimension(4), intent(in) :: c
 real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                     j3, s3, muH, Q, muJ, muS, R, G3, lambda, tau, R0 , &
                                     mu0, delta0, h, s31, s32

 real (dp)     , intent(out) :: res

  type (Alpha)                     :: alphaAll
  type (MatrixElements)            :: MatEl
  type (SingularMassless)          :: Sing
  type (AnomDim), dimension(3:6)   :: AnDim
  type (Model)   , dimension(2,2)  :: ModList
  real (dp)      , dimension(2,2)  :: clist = 1
  integer                          :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  ModList(1,1) = Model(lambda, [1._dp], c(:2), 'piece'); ModList(2,2) = ModList(1,1)
  ModList(1,2) = Model(lambda, [1._dp], c(3:), 'piece'); ModList(2,1) = ModList(1,2)

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  res = Sing%HJMSing(clist, ModList, 'piece', 'Model', &
  gap(:12), space(:3), cum(:4), order, isoft, s31, s32, R0, mu0, delta0, h, tau)

end subroutine f90SingularHJMPiece

!ccccccccccccccc

subroutine f90SingularDoublePiece(hard, gap, space, cum1, cum2, orderAlp, runAlp, order , &
 run, isoft, nf, j3, s3, s31, s32, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, &
 muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, rho1, rho2, res)

 use AlphaClass;  use MatrixElementsClass;  use SingularClass; use ModelClass
 use constants, only: dp; use AnomDimClass; implicit none

 integer, dimension(4), intent(in ) :: c
 character (len = *)  , intent(in ) :: cum1, cum2, space, hard, gap
 integer              , intent(in ) :: orderAlp, order, runAlp, run, nf, isoft
 real (dp)            , intent(in ) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                       j3, s3, muH, muJ, muS, R, G3, lambda, rho1, rho2, &
                                       R0, mu0, Q, delta0, h, s31, s32
 real (dp)            , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatrixElements)             :: MatEl
  type (SingularMassless)           :: Sing
  type (AnomDim) , dimension(3:6)   :: AnDim
  type (Model)   , dimension(2,2)   :: ModList
  real (dp)      , dimension(2,2)   :: clist = 1
  integer                           :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  ModList(1,1) = Model(lambda, [1._dp], c(:2), 'piece'); ModList(2,2) = ModList(1,1)
  ModList(1,2) = Model(lambda, [1._dp], c(3:), 'piece'); ModList(2,1) = ModList(1,2)

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless(MatEl, run, mu, 'thrust', hard(:6))

  res = Sing%DoubleSing(clist, modList, 'piece'  , &
                        'Model', gap(:12), space(:3), cum1(:4), cum2(:4), order, isoft, &
                         s31, s32, R0, mu0, delta0, h, rho1, rho2)

end subroutine f90SingularDoublePiece

!ccccccccccccccc

subroutine f90SingularPiece(hard, shape, gap, space, cum, orderAlp, runAlp, order, run, nf, &
j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, piece, &
lambda, R0, mu0, delta0, h, tau, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in) :: shape, cum, space, hard, gap
  integer              , intent(in) :: orderAlp, order, runAlp, run, nf
  integer, dimension(2), intent(in) :: piece
  real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, &
  Q, muJ, muS, R, G3, lambda, tau, R0, mu0, delta0, j3, s3, muH, h, mC, muC
  real (dp)           , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatricesElements)           :: MatEl
  type (SingularScales)             :: Sing
  type (Model)                      :: Mod
  type (AnomDim),  dimension(3:6)   :: AnDim
  integer                           :: i

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  Mod       = Model(lambda, [1._dp], piece, 'piece')
  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
 mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )

  call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
  call Sing%SetRunning(muJ, muS, R, mu)

  res = Sing%SingleSing(Mod, 'Model', gap(:12), space(:3), cum(:4), order, R0, &
  mu0, delta0, h, tau)

end subroutine f90SingularPiece

!ccccccccccccccc

subroutine f90SingularList(hard, shape, gap, space, cum, orderAlp, runAlp, order  , &
run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, &
mu, clen, lambda, R0, mu0, delta0, h, tau, res)

 use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
 use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in) :: shape, cum, space, hard, gap
  integer              , intent(in) :: orderAlp, order, runAlp, run, nf, clen
  real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                      Q, muJ, muS, R, G3, lambda, tau, R0, mu0, delta0,  &
                                      j3, s3, muH, h
  real (dp), dimension( (clen + 2) * (clen + 1)/2 ), intent(out) :: res

  integer                          :: i, j, k
  type (AnomDim), dimension(3:6)   :: AnDim
  type (Alpha)                     :: alphaAll
  type (MatrixElements)            :: MatEl
  type (SingularMassless)          :: Sing
  type (Model), dimension( (clen + 2) * (clen + 1)/2 ) :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i
      Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1
    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElements(alphaAll, nf, s3, s3, j3, Q, muH, muJ, muS, R, muLambda)
  Sing      = SingularMassless( MatEl, run, mu, shape(:6), hard(:6) )

  res = Sing%SingleSing(Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, delta0, h, tau)

end subroutine f90SingularList

!ccccccccccccccc

subroutine f90SingularDiffPiece(hard, shape, gap, space, cum, orderAlp, runAlp, order, run,&
nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu,   &
piece, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *)  , intent(in) :: shape, cum, space, hard, gap
  integer              , intent(in) :: orderAlp, order, runAlp, run, nf
  integer, dimension(2), intent(in) :: piece
  real (dp)            , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, muC, &
                                      Q, muJ, muS, R, G3, lambda, tau, R0, mu0, delta0, &
                                      j3, s3, muH, h, tau2
  real (dp)           , intent(out) :: res
  type (Alpha)                      :: alphaAll
  type (MatricesElements)           :: MatEl
  type (SingularScales)             :: Sing
  type (Model)                      :: Mod
  type (AnomDim), dimension(3:6)    :: AnDim
  integer                           :: i

 do i = 3, 6
   AnDim(i) = AnomDim('MSbar', i, G3)
 end do

 Mod       = Model(lambda, [1._dp], piece, 'piece')
 alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                   mB, muB, mC, muC, 'analytic', 0._dp)
 MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
 Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )

 call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
 call Sing%SetRunning(muJ, muS, R, mu)

 res = Sing%SingleSing(Mod, 'Model', gap(:12), space(:3), cum(:4), order, R0, mu0, &
                       delta0, h, tau, tau2)

end subroutine f90SingularDiffPiece

!ccccccccccccccc

subroutine f90SingularDiffList(hard, shape, gap, space, cum, orderAlp, runAlp, order,&
 run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, &
 mu, clen, lambda, R0, mu0, delta0, h, tau, tau2, res)

  use AlphaClass;  use MatrixElementsClass; use SingularClass; use ModelClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in) :: shape, cum, space, hard, gap
  integer            , intent(in) :: orderAlp, order, runAlp, run, nf, clen
  real (dp)          , intent(in) :: mZ, amZ, mu, muLambda, mT, muT, mB, muB, mC, &
                                     Q, muJ, muS, R, G3, lambda, tau, R0, delta0, &
                                     j3, s3, muH, h, tau2, mu0, muC
  integer                         :: i, j, k
  type (Alpha)                    :: alphaAll
  type (MatricesElements)         :: MatEl
  type (SingularScales)           :: Sing
  type (AnomDim), dimension(3:6)  :: AnDim
  real (dp)   , dimension( (clen + 2) * (clen + 1)/2 ) , intent(out) :: res
  type (Model), dimension( (clen + 2) * (clen + 1)/2 )               :: Mod

  k = 1

  do i = 0, clen
    do j = 0, i

    Mod(k) = Model(lambda, [1._dp], [i,j], 'piece'); k = k + 1

    end do
  end do

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, G3)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
                    mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatricesElements(alphaAll, nf, s3, s3, j3, muLambda)
  Sing      = SingularScales( MatEl, run, shape(:6), hard(:6) )

   call Sing%setHard(Q, muH); call sing%SetMat(muJ, muS)
   call Sing%SetRunning(muJ, muS, R, mu)

  res = Sing%SingleSing(Mod, gap(:12), space(:3), cum(:4), order, R0, mu0, &
  delta0, h, tau, tau2)

end subroutine f90SingularDiffList

!ccccccccccccccc

subroutine f90Rhad(str, orderAlp, runAlp, order, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, mu, Q, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: order, runAlp, orderAlp, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, Q
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)
  alphaMass = Running(nf, runAlp, alphaAll, muC)
  EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)

  MatEl = Sigma(alphaMass, EW);  res = MatEl%Rhad(order, mu, Q)

end subroutine f90Rhad

!ccccccccccccccc

subroutine f90SigmaHad(str, curr, orderAlp, runAlp, order, nf, mZ, gammaZ, &
sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, mu, Q, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, &
  Q, gammaZ, sin2ThetaW, amZQED
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runAlp, alphaAll, muC)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaHad( curr(:6), order, mu, Q )

end subroutine f90SigmaHad

!ccccccccccccccc

subroutine f90SigmaRad(str, curr, orderAlp, runAlp, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x, theta, &
res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  Q, gammaZ, sin2ThetaW, amZQED, x, theta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runAlp, alphaAll, muC)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaRad( curr(:6), order, eH, Q, x, theta )

end subroutine f90SigmaRad

!ccccccccccccccc

subroutine f90SigmaRadCum(str, curr, orderAlp, runAlp, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x0, x1, &
theta, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  Q, gammaZ, sin2ThetaW, amZQED, x0, x1, theta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runAlp, alphaAll, muC)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaRadCum( curr(:6), order, eH, Q, x0, x1, theta )

end subroutine f90SigmaRadCum

!ccccccccccccccc

subroutine f90SigmaRadCone(str, curr, orderAlp, runAlp, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x, theta, &
deltaTheta, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  Q, gammaZ, sin2ThetaW, amZQED, x, theta, deltaTheta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runAlp, alphaAll, muC)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaRadCone( curr(:6), order, eH, Q, x, theta,&
  deltaTheta )

end subroutine f90SigmaRadCone

!ccccccccccccccc

subroutine f90SigmaRadConeCum(str, curr, orderAlp, runAlp, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x0, x1, &
theta, deltaTheta, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  Q, gammaZ, sin2ThetaW, amZQED, x0, x1, theta, deltaTheta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runAlp, alphaAll, muC)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaRadConeCum( curr(:6), order, eH, Q, x0, x1, &
  theta, deltaTheta )

end subroutine f90SigmaRadConeCum

!ccccccccccccccc

subroutine f90RhadCoefs(nf, res)

  use RunningClass; use AlphaClass; use SigmaClass;  use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  integer                , intent(in ) :: nf
  real (dp), dimension(4), intent(out) :: res
  type (Running)                       :: alphaMass
  type (Alpha)                         :: alphaAll
  type (Sigma)                         :: MatEl
  type (ElectroWeak)                   :: EW
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, &
  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp)
  alphaMass = Running(nf, 0, alphaAll, 0._dp)
  EW = ElectroWeak(0._dp, 2.4952_dp, 0.23119_dp)

  MatEl = Sigma(alphaMass, EW);  res = MatEl%RhadCoefs()

end subroutine f90RhadCoefs

!ccccccccccccccc

subroutine f90SigmaMass(str, curr, orderAlp, runAlp, runMass, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, mu, Q, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, &
  gammaZ, sin2ThetaW, Q, amZQED
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaMass(curr(:6), order, mu, Q)

end subroutine f90SigmaMass

!ccccccccccccccc

subroutine f90SigmaMassRad(str, curr, orderAlp, runAlp, runMass, order, &
nf, mZ, gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x, &
theta, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  gammaZ, sin2ThetaW, Q, amZQED, x, theta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaMassRad(curr(:6), order, eH, Q, x, theta)

end subroutine f90SigmaMassRad

!ccccccccccccccc

subroutine f90SigmaMassRadCum(str, curr, orderAlp, runAlp, runMass, order, &
nf, mZ, gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x0, &
x1, theta, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  gammaZ, sin2ThetaW, Q, amZQED, x0, x1, theta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaMassRadCum(curr(:6), order, eH, Q, x0, x1, theta)

end subroutine f90SigmaMassRadCum

!ccccccccccccccc

subroutine f90SigmaMassRadCone(str, curr, orderAlp, runAlp, runMass, order, &
nf, mZ, gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, Q, x, &
theta, deltaTheta, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  gammaZ, sin2ThetaW, Q, amZQED, x, theta, deltaTheta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaMassRadCone(curr(:6), order, eH, Q, x, theta, deltaTheta)

end subroutine f90SigmaMassRadCone

!ccccccccccccccc

subroutine f90SigmaMassRadConeCum(str, curr, orderAlp, runAlp, runMass,  &
order, nf, mZ, gammaZ, sin2ThetaW, amZ, amZQED, mT, muT, mB, muB, mC, muC, eH, &
Q, x0, x1, theta, deltaTheta, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, eH, mT, muT, mB, muB, mC, muC, &
  gammaZ, sin2ThetaW, Q, amZQED, x0, x1, theta, deltaTheta
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, &
  muC, 'analytic', amZQED)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%SigmaMassRadConeCum(curr(:6), order, eH, Q, x0, x1, &
  theta, deltaTheta)

end subroutine f90SigmaMassRadConeCum

!ccccccccccccccc

subroutine f90RhadMass(str, curr, orderAlp, runAlp, runMass, order, nf, mZ, &
gammaZ, sin2ThetaW, amZ, mT, muT, mB, muB, mC, muC, mu, Q, res)

  use RunningClass; use AlphaClass; use SigmaClass; use ElectroWeakClass
  use constants, only: dp; use AnomDimClass; implicit none

  character (len = *), intent(in ) :: str, curr
  integer            , intent(in ) :: order, runAlp, orderAlp, nf, runMass
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC, gammaZ, &
  sin2ThetaW, Q
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  type (Sigma)                     :: MatEl
  type (ElectroWeak)               :: EW
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(str(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  alphaMass = Running(nf, runMass, alphaAll, 0._dp)
  EW        = ElectroWeak(mZ, gammaZ, sin2ThetaW)
  MatEl     = Sigma(alphaMass, EW)
  res       = MatEl%RhadMass(curr(:6), order, mu, Q)

end subroutine f90RhadMass

!ccccccccccccccc

subroutine f90RQCD(str, runAlp, runMass, ordMass, order, ord1S, R1S, method, &
lambda, gt, mZ, amZ, mT, h, Q, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, h, mT, Q, gt, lambda, R1S
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%RQCD(order, h, Q)

end subroutine f90RQCD

!ccccccccccccccc

subroutine f90Rexp(str, runAlp, runMass, ordMass, order, ord1S, R1S, method, &
lambda, gt, mZ, amZ, mT, mu, nu, Q, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%RExp(order, mu, nu, Q)

end subroutine f90Rexp

!ccccccccccccccc

subroutine f90Rmatched(str, runAlp, runMass, ordMass, order, ord1S, R1S, method, &
lambda, gt, mZ, amZ, mT, mu, nu, v1, v2, Q, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S, v1, v2
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%Rmatched(order, mu, nu, v1, v2, Q)

end subroutine f90Rmatched

!ccccccccccccccc

subroutine f90SigmaMatched(str, runAlp, runMass, ordMass, order, ord1S, R1S, method, &
lambda, gt, mZ, gammaZ, sinW, amZ, amZQED, mT, mu, nu, v1, v2, Q, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S, &
  v1, v2, gammaZ, sinW, amZQED
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', amZQED)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, gammaZ, sinW)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%SigmaMatched(order, mu, nu, v1, v2, Q)

end subroutine f90SigmaMatched

!ccccccccccccccc

subroutine f90SigmaMatchedRad(str, runAlp, runMass, ordMass, order, &
ord1S, R1S, method, lambda, gt, mZ, gammaZ, sinW, amZ, amZQED, mT, mu, nu, &
v1, v2, Q, x, theta, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S, &
  v1, v2, gammaZ, sinW, amZQED, x, theta
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', amZQED)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, gammaZ, sinW)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%SigmaMatchedRad(order, mu, nu, v1, v2, Q, x, theta)

end subroutine f90SigmaMatchedRad

!ccccccccccccccc

subroutine f90SigmaMatchedRadCum(str, runAlp, runMass, ordMass, order, &
ord1S, R1S, method, lambda, gt, mZ, gammaZ, sinW, amZ, amZQED, mT, mu, nu, &
v1, v2, Q, x0, x1, theta, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S, &
  v1, v2, gammaZ, sinW, amZQED, x0, x1, theta
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', amZQED)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, gammaZ, sinW)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%SigmaMatchedRadCum(order, mu, nu, v1, v2, Q, x0, x1, theta)

end subroutine f90SigmaMatchedRadCum

!ccccccccccccccc

subroutine f90SigmaMatchedRadCone(str, runAlp, runMass, ordMass, order, &
ord1S, R1S, method, lambda, gt, mZ, gammaZ, sinW, amZ, amZQED, mT, mu, nu, v1,&
v2, Q, x, theta, deltaTheta, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S, &
  v1, v2, gammaZ, sinW, amZQED, x, theta, deltaTheta
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', amZQED)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, gammaZ, sinW)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%SigmaMatchedRadCone(order, mu, nu, v1, v2, Q, x, theta, &
  deltaTheta)

end subroutine f90SigmaMatchedRadCone

!ccccccccccccccc

subroutine f90SigmaMatchedRadConeCum(str, runAlp, runMass, ordMass, order, &
ord1S, R1S, method, lambda, gt, mZ, gammaZ, sinW, amZ, amZQED, mT, mu, nu, v1,&
v2, Q, x0, x1, theta, deltaTheta, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, mu, nu, mT, Q, gt, lambda, R1S, &
  v1, v2, gammaZ, sinW, amZQED, x0, x1, theta, deltaTheta
  real (dp)          , intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', amZQED)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, gammaZ, sinW)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%SigmaMatchedRadConeCum(order, mu, nu, v1, v2, Q, x0, x1, &
  theta, deltaTheta)

end subroutine f90SigmaMatchedRadConeCum

!ccccccccccccccc

subroutine f90RmatchedList(str, runAlp, runMass, ordMass, order, ord1S, R1S, &
method, lambda, gt, mZ, amZ, mT, h, hnu, v1, v2, Q0, Q1, deltaQ, res)

  use Constants; use RNRQCDClass; use AnomDimClass; use AlphaClass
  use RunningClass; use VFNSMSRClass; use ElectroWeakClass; implicit none

  character (len = *), intent(in ) :: str, method
  integer            , intent(in ) :: order, runAlp, runMass, ordMass, ord1S
  real (dp)          , intent(in ) :: mZ, amZ, h, hnu, mT, Q0, Q1, deltaQ, gt, &
  lambda, R1S, v1, v2
  real (dp), dimension( 2,0:Floor( (Q1 - Q0)/deltaQ ) ), intent(out) :: res
  integer                          :: i
  type (ElectroWeak)               :: EW
  type (RNRQCD)                    :: NRQCD
  type (Alpha)                     :: alphaAll
  type (Running), dimension(2)     :: alphaMass
  type (VFNSMSR)                   :: MSR
  type (AnomDim), dimension(3:6)   :: AnDim
  character (len = 5)              :: alphaScheme

  alphaScheme = 'pole'; if ( str(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, runAlp, runAlp, mZ, aMz, mT, mT, &
  0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  alphaMass = [ Running(4, runMass, alphaAll, 1._dp), &
  Running(5, runMass, alphaAll, 1._dp) ]

  MSR = VFNSMSR(alphaMass);   EW = ElectroWeak(mZ, 2.4952_dp, 0.23119_dp)
  NRQCD = RNRQCD(MSR, EW, str, method, gt, ordMass, ord1S, R1S, lambda)

  res = NRQCD%RmatchedList(order, h, hnu, v1, v2, Q0, Q1, deltaQ)

end subroutine f90RmatchedList

!ccccccccccccccc

subroutine f90lambdaQCD(str, order, runAlp, run, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, mu, res)

  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none

  character (len = *), intent(in ) :: str
  integer            , intent(in ) :: order, runAlp, run, nf
  real (dp)          , intent(in ) :: mZ, amZ, mu, mT, muT, mB, muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Running)                   :: alphaMass
  type (Alpha)                     :: alphaAll
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim( str(:5), i, 0._dp )
  end do

  alphaAll  = Alpha(AnDim, order, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)
  alphaMass = Running(nf, run, alphaAll, mu)
  res       = alphaMass%lambdaQCD(run)

end subroutine f90lambdaQCD

!ccccccccccccccc

subroutine f90delta(str, nf, mu, R, res)
  use AlphaClass; use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass; implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp)              , intent(in ) :: mu, R
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, &
  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 'analytic', 0._dp)

  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, &
  mu, mu, mu, tiny(1._dp), mu, R, R, 0._dp, 0._dp)

  res       = MatEl%delta( str(:12) )

end subroutine f90delta

!ccccccccccccccc

subroutine f90PSdelta(orderAlp, runAlp, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, R, lg, res)
  use AlphaClass;    use RunningClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer                , intent(in ) :: nf, orderAlp, runAlp
  real (dp)              , intent(in ) :: mu, R, mZ, amZ, mT, muT, mB, muB, mC, muC, lg
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (Running)                       :: alphaMass
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = Running(nf, 0, alphaAll, 0._dp)

  res       = alphaMass%PSdelta( R, mu, lg )

end subroutine f90PSdelta

!ccccccccccccccc

subroutine f90deltaGap(str, orderAlp, runAlp, runMass, nf, mZ, amZ, mT, muT, &
mB, muB, mC, muC, mu, R, res)
  use AlphaClass;    use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf, orderAlp, runAlp, runMass
  real (dp)              , intent(in ) :: mu, R, mZ, amZ, mT, muT, mB, muB, mC, muC
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  real (dp)                            :: muM
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  muM = tiny(1._dp)

  if ( str(:3) == 'MSR' .or. str(:2) == 'RS' .or. str(:8) == 'MSbarLow' &
  .or. str(:7) == 'JetMass' ) muM = 2 * mu

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl     = MatrixElementsMass(alphaAll, nf, runMass, 0._dp, 0._dp, 0._dp, &
  0._dp, mu, mu, mu, muM, mu, R, R, 0._dp, 0._dp)

  res       = MatEl%delta( str(:12) )

end subroutine f90deltaGap

!ccccccccccccccc

subroutine f90deltaMSbar(orderAlp, runAlp, run, nf, mZ, amZ, mT, muT, mB, &
muB, mC, muC, mu, res)
  use RunningClass; use AlphaClass;  use MatrixElementsClass
  use constants, only: dp; use AnomDimClass; implicit none
  integer                , intent(in ) :: orderAlp, runAlp, run, nf
  real (dp)              , intent(in ) :: mZ, amZ, mu, mT, muT, muB, mB, mC, muC
  real (dp), dimension(4), intent(out) :: res
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, mT, muT, &
  mB, muB, mC, muC, 'analytic', 0._dp)

  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, mu, mu, &
  tiny(1._dp), mu, mu, mu, muC, muC)

  res       = MatEl%delta('MSbar')

end subroutine f90deltaMSbar

!ccccccccccccccc

subroutine f90CoefMat(str, nf, s3, res)
  use AlphaClass; use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  character (len = *)        , intent(in ) :: str
  integer                    , intent(in ) :: nf
  real (dp)                  , intent(in ) :: s3
  real (dp), dimension(0:4,3), intent(out) :: res
  real (dp)                                :: mu = 20._dp
  type (Alpha)                             :: alphaAll
  type (MatrixElementsMass)                :: MatEl
  integer                                  :: i
  type (AnomDim), dimension(3:6)           :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                     0._dp, 0._dp, 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, s3, s3, s3, s3, mu, mu, mu, tiny(1._dp),&
                             mu, mu, mu, 0._dp, 0._dp)
  res = MatEl%CoefMat(str)

end subroutine f90CoefMat

!ccccccccccccccc

subroutine f90NGLintegral(nf, pow, w1, w2, res)
  use AlphaClass; use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer     , intent(in ) :: nf, pow
  real (dp)   , intent(in ) :: w1, w2
  real (dp)   , intent(out) :: res
  real (dp)                 :: mu = 20._dp
  type (MatrixElementsMass) :: MatEl
  type (Alpha)              :: alphaAll
  integer                   :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                     0._dp, 0._dp, 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, mu, mu, &
                         tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)

 res = MatEl%NGLIntegral(pow, w1, w2)

end subroutine f90NGLintegral

!ccccccccccccccc

subroutine f90NGLDoubleintegral(nf, pow, w1, w2, ratio, res)
  use MatrixElementsClass;  use constants, only: dp; implicit none
  integer        , intent(in ) :: nf, pow
  real (dp)      , intent(in ) :: w1, w2, ratio
  real (dp)      , intent(out) :: res

  res = NGLDoubleIntegral(nf, pow, w1, w2, ratio)

end subroutine f90NGLDoubleintegral

!ccccccccccccccc

subroutine f90NGLfunction(nf, x, res)
  use AlphaClass;  use MatrixElementsClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  integer     , intent(in ) :: nf
  real (dp)   , intent(in ) :: x
  real (dp)   , intent(out) :: res
  real (dp)                 :: mu = 20._dp
  type (Alpha)              :: alphaAll
  type (MatrixElementsMass) :: MatEl
  integer                   :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                     0._dp, 0._dp, 0._dp)
  MatEl    = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, mu, mu, &
                         tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)
  res      = MatEl%NGLfunction(x)

end subroutine f90NGLfunction

!ccccccccccccccc

subroutine f90DiffDeltaGap(str, scheme, order, R0, R, mu0, mu, muLambda, orderAlp, &
                           runAlp, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, res)
  use MatrixElementsClass; use AlphaClass;  use constants, only: dp
  use AnomDimClass;  implicit none
  character (len = *), intent(in ) :: str, scheme
  integer            , intent(in ) :: order, nf, orderAlp, runAlp
  real (dp)          , intent(in ) :: R0, R, mu0, mu, muLambda, mZ, amZ, mT, muT, mB, &
                                      muB, mC, muC
  real (dp)          , intent(out) :: res
  type (Alpha)                     :: alphaAll
  type (MatrixElementsMass)        :: MatEl
  integer                          :: i
  type (AnomDim), dimension(3:6)   :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim(scheme(:5), i, 0._dp)
  end do

  alphaAll  = Alpha(Andim, orderAlp, runAlp, mZ, amZ, mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, 5, 0, 0._dp, 0._dp, 0._dp, 0._dp, mu, &
  mu, mu, tiny(1._dp), mu, mu, mu, muLambda, muLambda)

  res       = MatEl%DiffDeltaGap( str(:12), order, R0, R, mu0, mu )

end subroutine f90DiffDeltaGap

!ccccccccccccccc

subroutine f90kernels(n, width, w, mu, p, res)
  use KernelsClass;  use constants, only: dp; implicit none
  integer  , intent(in ) :: n
  real (dp), intent(in ) :: width, w, mu, p
  real (dp), intent(out) :: res(0:n)
  type (Kernels)         :: ker

  ker = Kernels(n, width, w, mu, p);  res = ker%KernelList()

end subroutine f90kernels

!ccccccccccccccc

subroutine f90GammaDerList(n, w, res)
  use KernelsClass;  use constants, only: dp; implicit none
  integer   , intent(in) :: n
  real (dp), intent(in ) :: w
  real (dp), intent(out) :: res(0:n)
  type (Kernels)         :: ker

 ker = Kernels(n, 0._dp, w, 3._dp, 3._dp)
 res = ker%DerList()

end subroutine f90GammaDerList

!ccccccccccccccc

subroutine f90NGLkernels(n, n1, n2, width, w, mu, p, res)
  use KernelsClass;  use constants, only: dp; implicit none
  integer                , intent(in ) :: n, n1, n2
  real (dp), dimension(2), intent(in ) :: w, p
  real (dp)              , intent(in ) :: mu, width
  real (dp)              , intent(out) :: res(n)
  type (Kernels)                       :: ker1, ker2
  real (dp), dimension(0:2*n+n1)       :: kerList1
  real (dp), dimension(0:2*n+n2)       :: kerList2

  ker1     = Kernels( 2*n + n1, width, w(1), mu, p(1) )
  ker2     = Kernels( 2*n + n2, width, w(2), mu, p(2) )

  kerList1 = ker1%KernelList(); kerList2 = ker2%KernelList()

  res      = KernelsCombine(kerList1, kerList2, n)

end subroutine f90NGLkernels

!ccccccccccccccc

subroutine f90hyper2f1 ( a, b, c, x, res )
  use hyper; use constants, only: dp; implicit none
  real (dp), intent(in)  :: a, b, c, x
  real (dp), intent(out) :: res

  res = hyper2f1(a, b, c, x)

end subroutine f90hyper2f1

!ccccccccccccccc

subroutine f90hyperf32exact(w2, x2, res)
  use hyper; use constants, only: dp; implicit none
  real (dp), intent(in ) :: w2, x2
  real (dp), intent(out) :: res

  res = H3F2Exact(w2, x2)

end subroutine f90hyperf32exact

!ccccccccccccccc

subroutine f90GammaR(str, nf, beta)
  use MatrixElementsClass;  use AlphaClass; use constants
  use AnomDimClass;  implicit none
  character (len = *)    , intent(in ) :: str
  integer                , intent(in ) :: nf
  real (dp), dimension(3), intent(out) :: beta
  real (dp)                            :: mu = 20._dp
  type (Alpha)                         :: alphaAll
  type (MatrixElementsMass)            :: MatEl
  integer                              :: i
  type (AnomDim), dimension(3:6)       :: AnDim

  do i = 3, 6
    AnDim(i) = AnomDim('MSbar', i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, 0, 0, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                    0._dp, 0._dp, 0._dp)
  MatEl     = MatrixElementsMass(alphaAll, nf, 0, 0._dp, 0._dp, 0._dp, mu, 0._dp, mu, mu, &
                            tiny(1._dp), mu, mu, mu, 0._dp, 0._dp)

  beta = MatEl%GammaR( str(:7) )

end subroutine f90GammaR

!ccccccccccccccc
