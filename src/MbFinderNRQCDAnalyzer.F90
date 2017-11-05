Program MbFinderNRQCDAnalyzer

  use constants, only: dp; implicit none
  character (len = 15)                    :: method, scheme, charm, average, &
  counting, iter, arg
  character (len = 200)                   :: line
  integer                                 :: orderAlp, runAlp, order, run, &
  nl, n, l, j, s, ord, ierror, i
  real (dp)                               :: mZ, amZ, mT, muT, mB, muB, mC, &
  muC, lambda1, lam, lambda2, mass, mean, maxim, minim, Rmin, R, mu

  call getArg(1, arg); read(arg, *) Rmin

  read*, n, l, j, s
  read*, iter, charm, scheme, average, method, counting
  read*, orderAlp, runAlp, order, run, nl, ord
  read*, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam
  read*, mass

  print*, n, l, j, s
  print*, iter, charm, scheme, average, method, counting
  print*, orderAlp, runAlp, order, run, nl, ord
  write(*,'(11F10.4)') mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam
  write(*,'(F12.8)'  ) mass

  read*, ; read*, ; read*,; print*,
  print*, '  m(m)      simga'; print*,

  i = 0; maxim = 0; minim = 1000; mean = 0

  do

    read(*,'(A)', iostat = ierror) line
    if (ierror == -1) exit
    i = i + 1
    read(line,*) R, mu, mass; if (R < Rmin .or. mu < Rmin) cycle
    mean = mean + mass; maxim = max(mass,maxim); minim = min(minim, mass)
  end do

  write(*,'(2F11.7,F25.13)') mean/i, (maxim - minim)/2

end program MbFinderNRQCDAnalyzer
