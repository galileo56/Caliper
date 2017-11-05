Program Chi2NRQCDAnalyzer

  use constants, only: dp; implicit none
  character (len = 15)                    :: method, scheme, charm, average, &
  counting, iter
  character (len = 200)                   :: dataFile, line
  integer                                 :: orderAlp, runAlp, order, run, &
  nl, n, dof, ierror, i
  real (dp)                               :: mZ, amZ, mT, muT, mB, muB, mC, &
  muC, lambda1, lam, lambda2, mass, err, chi2, mean, maxim, minim, chi, sigma

  read*, dataFile
  read*, iter, charm, scheme, average, method, counting
  read*, orderAlp, runAlp, order, run, nl, n
  read*, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam

  print*, trim( dataFile )
  print*, iter, charm, scheme, average, method, counting
  print*, orderAlp, runAlp, order, run, nl, n
  write(*,'(11F10.4)') mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam

  read*, ; read*, ; read*,; print*,
  print*, '  m(m)      simgaExp   sigmaPert      chi2/dof'; print*,

  OPEN (UNIT = 1, FILE = dataFile, ACCESS = 'SEQUENTIAL', STATUS = 'OLD')
  read(1,*) dof; dof = dof - 1; close(2)

  i = 0; maxim = 0; minim = 1000; mean = 0; chi = 0; sigma = 0

  do

    read(*,'(A)', iostat = ierror) line
    if (ierror == -1) exit
    i = i + 1
    read(line,*) mass, mass, mass, err, chi2
    mean = mean + mass; chi = chi + chi2; sigma = sigma + sqrt(chi2/dof) * err
    maxim = max(mass,maxim); minim = min(minim, mass)
  end do

  write(*,'(3F11.7,F25.13)') mean/i, sigma/i, (maxim - minim)/2, chi2/dof

end program Chi2NRQCDAnalyzer
