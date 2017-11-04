Program NRQCDGrid
  use RunningClass;  use AlphaClass;  use constants, only: dp
  use AnomDimClass;  use NRQCDClass;  use VFNSMSRClass;  implicit none

  character (len = 15)                    :: method, scheme, charm, average, &
  counting, alphaScheme
  real (dp), dimension(0:4)               :: res
  type (NRQCD)                            :: Upsilon
  type (VFNSMSR)                          :: MSR
  type (Alpha)                            :: alphaAll
  type (Running), dimension(2)            :: alphaMass
  type (AnomDim), dimension(3:6)          :: AnDim
  integer                                 :: orderAlp, runAlp, order, run, nl, &
  i, k, m, imax, kmax, mMax, n, l, j, s
  real (dp)                               :: mZ, amZ, mT, muT, mB, muB, mC, R, &
  muC, lambda1, lam, lambda2, R0, R1, deltaR, mu0, mu1, deltaMu, mu, m0, m1,   &
  deltaM, mass

  read*, charm, scheme, average, method, counting
  read*, n, l, j, s
  read*, orderAlp, runAlp, order, run, nl
  read*, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam
  read*, R0, R1, deltaR, mu0, mu1, deltaMu, m0, m1, deltaM

  print*, n, l, j, s
  print*, charm, scheme, average, method, counting
  print*, orderAlp, runAlp, order, run, nl
  write(*,'(11F10.4)') mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam
  write(*,'(F12.8)'  ) mass

  imax = Floor( (mu1 - mu0)/deltaMu ) + 1;  kmax = Floor( (R1 - R0)/deltaR ) + 1
  mMax = Floor( (m1 - m0)/deltaM ) + 1

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

  print*,; print*, '  R     mu    m(m)        Tree ', &
  '         O(e0)         O(e1)         O(e2)         O(e3)'; print*,

  do i = 1, imax
      mu = mu0 + (i - 1) * deltaMu
    do k = 1, kmax
      R = R0 + (k - 1) * deltaR

      do m = 1, mMax

        mass = m0 + (m - 1) * deltaM;  call Upsilon%setMass( mass )
        res = Upsilon%En( order, mu, R, lam, method(:8), counting(:5) )

        write( *, '(2F6.2, F10.6, 5F14.8)' ) mu, R, mass, res

      end do

    end do
  end do

end program NRQCDGrid
