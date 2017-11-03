Program Chi2MbNRQCD
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = 200)                   :: dataFile
  character (len = 15)                    :: method, scheme, charm, average, &
  counting, iter, alphaScheme, Fit
  real (dp), dimension(3)                 :: res
  real (dp), dimension(6)                 :: res2
  real(dp), dimension(:,:), allocatable   :: dataList
  real(dp), dimension(:)  , allocatable   :: RList, muList
  integer , dimension(:,:), allocatable   :: qnList
  type (NRQCD), dimension(:), allocatable :: Upsilon
  type (VFNSMSR)                          :: MSR
  type (Alpha)                            :: alphaAll
  type (Running), dimension(2)            :: alphaMass
  type (AnomDim), dimension(3:6)          :: AnDim
  integer                                 :: orderAlp, runAlp, order, run, m, &
  nl, n, ndim, i, j, imax, jmax
  real (dp)                               :: mZ, amZ, mT, muT, mB, muB, mC, &
  muC, lambda1, lam, lambda2, mH, R0, R1, deltaR, mu0, mu1, deltaMu, R, mu

  read*, dataFile
  read*, fit
  read*, iter, charm, scheme, average, method, counting
  read*, orderAlp, runAlp, order, run, nl, n
  read*, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam
  read*, R0, R1, deltaR, mu0, mu1, deltaMu

  print*, dataFile
  print*, iter, charm, scheme, average, method, counting
  print*, orderAlp, runAlp, order, run, nl, n
  write(*,'(11F10.4)') mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam

  imax = Floor( (mu1 - mu0)/deltaMu ) + 1;  jmax = Floor( (R1 - R0)/deltaR ) + 1

  call readData(); ndim = maxval( qnList(1,:) ); res = 0; res2 = 0

  allocate( RList(ndim), muList(ndim) )

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  allocate( Upsilon(m) )

  do i = 1, m
    Upsilon(i) = NRQCD( charm(:7), scheme(:5), average(:3), MSR, qnList(1,i), &
    qnList(2,i), qnList(3,i), qnList(4,i) )
  end do

  if (nl == 5) mH = mT;  if (nl == 4) mH = mB;  if (nl == 3) mH = mC

  print*,
  if ( fit(:4) == 'mass' ) then
    print*, '  R     mu    m(m)   DeltaM       chi2'
    res(1) = mH
  else if ( fit(:9) == 'alphaMass' ) then
    print*, '  R     mu   m(m)      DeltaM     alpha    DeltaAlpha   rho            chi2'
    res2(1:3:2) = [mH, amZ]
  else if ( fit(:5) == 'alpha' ) then
    print*, '  R     mu     alpha   DeltaAlpha     chi2'
    res(1) = amZ
  end if
  print*,

  do i = 1, imax
      mu = mu0 + (i - 1) * deltaMu
      muList = 1.5_dp + 2.5 * (mu - 1)/3
      if (ndim  > 1) muList(2) = mu
    do j = 1, jmax
      R = R0 + (j - 1) * deltaR
      RList = 1.5_dp + 2.5 * (R - 1)/3
      if (ndim  > 1) RList(2) = R

    if ( fit(:4) == 'mass' ) then

      res = Chi2MinNRQCD( Upsilon, dataList, m, iter(:10), order, n, res(1),  &
      muList, RList, ndim, lam, method(:8), counting(:5) )
      write( *, '(2F6.2,2F9.5,F25.13)' ) mu, R, res

    else if ( fit(:9) == 'alphaMass' ) then

      res2 = Chi2MinAlphaMbNRQCD( Upsilon, dataList, m, iter(:10), order, n, &
      res2(1), res2(3), muList, RList, ndim, lam, method(:8), counting(:5) )
      write( *, '(2F6.2,5F11.7,F25.13)' ) mu, R, res2

    else if ( fit(:5) == 'alpha' ) then

      res = Chi2MinAlphaNRQCD( Upsilon, dataList, m, iter(:10), order, n, res(1), &
      muList, RList, ndim, lam, method(:8), counting(:5) )
      write( *, '(2F6.2,2F11.7,F25.13)' ) mu, R, res

    end if

    end do
  end do

  deallocate( dataList, Upsilon, muList, Rlist )

contains

!ccccccccccccccc

subroutine readData()
  use constants; implicit none
  integer                           :: ierror, j
  character (len = 100)             :: readString

  OPEN (UNIT=1, FILE = dataFile, ACCESS = 'SEQUENTIAL', STATUS = 'OLD')

  read(1,*) m; read(1,*) ; allocate( dataList(2,m), qnList(4,m) )

  do j = 1, m
    read(1,'(A)', iostat = ierror) readString
    if (ierror == -1) exit
    read(readString,*) qnList(:,j), dataList(:,j)
  end do

  close(1)

end subroutine readData

end program Chi2MbNRQCD
