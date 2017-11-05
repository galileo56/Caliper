Program MbFinderNRQCD
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = 15)                    :: method, scheme, charm, average, &
  counting, iter, alphaScheme
  type (NRQCD)                            :: Upsilon
  type (VFNSMSR)                          :: MSR
  type (Alpha)                            :: alphaAll
  type (Running), dimension(2)            :: alphaMass
  type (AnomDim), dimension(3:6)          :: AnDim
  integer                                 :: orderAlp, runAlp, order, run, nl, &
  ord, i, k, imax, kmax, n, l, j, s
  real (dp)                               :: mZ, amZ, mT, muT, mB, muB, mC, R, &
  muC, lambda1, lam, lambda2, R0, R1, deltaR, mu0, mu1, deltaMu, mu, res, mass

  read*, iter, charm, scheme, average, method, counting
  read*, n, l, j, s
  read*, orderAlp, runAlp, order, run, nl, ord
  read*, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam, mass
  read*, R0, R1, deltaR, mu0, mu1, deltaMu

  print*, n, l, j, s
  print*, iter, charm, scheme, average, method, counting
  print*, orderAlp, runAlp, order, run, nl, ord
  write(*,'(11F10.4)') mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam
  write(*,'(F12.8)'  ) mass

  imax = Floor( (mu1 - mu0)/deltaMu ) + 1;  kmax = Floor( (R1 - R0)/deltaR ) + 1

  alphaScheme = 'pole'; if ( scheme(:4) /= 'pole' ) alphaScheme = 'MSbar'

  do i = 3, 6
    AnDim(i) = AnomDim(alphaScheme, i, 0._dp)
  end do

  alphaAll  = Alpha(AnDim, orderAlp, runAlp, mZ, amZ, &
  mT, muT, mB, muB, mC, muC, 'analytic', 0._dp)

  alphaMass = [ Running(nl - 1, run, alphaAll, lambda2), &
  Running(nl, run, alphaAll, lambda1) ]

  MSR = VFNSMSR(alphaMass)

  Upsilon = NRQCD( charm(:7), scheme(:5), average(:3), MSR, n, l, j, s )

  print*,; print*, '  R     mu    m(m)'; print*,

  do i = 0, imax
    mu = mu0 + i * deltaMu
    do k = 0, kmax
      R = R0 + k * deltaR

      res = Upsilon%MassFitter( iter(:10), ord, order, mu, R, mass, &
      lam, method(:8), counting(:5) )

      write( *, '(2F6.2,F11.7)' ) mu, R, res

    end do
  end do

end program MbFinderNRQCD
