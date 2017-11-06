Program Chi2NRQCDAnalyzer

  use constants, only: dp; implicit none
  character (len = 15)                    :: method, scheme, charm, average, &
  counting, iter, fit
  character (len = 200)                   :: dataFile, line
  integer                                 :: orderAlp, runAlp, order, run, &
  nl, n, dof, ierror, i
  real (dp)                               :: mZ, amZ, mT, muT, mB, muB, mC,    &
  muC, lambda1, lam, lambda2, mass, err, chi2, mean, maxim, minim, chi, sigma, &
  alpha, rho, errAlp, fact, sigmaAlp, maxAlp, minAlp, meanAlp, meanRho, corr,  &
  meanMass2, meanAlp2, meanMassAlp, R, mu

  read*, dataFile
  read*, fit
  read*, iter, charm, scheme, average, method, counting
  read*, orderAlp, runAlp, order, run, nl, n
  read*, mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam

  print*, trim( dataFile )
  print*, trim( fit )
  print*, iter, charm, scheme, average, method, counting
  print*, orderAlp, runAlp, order, run, nl, n
  write(*,'(11F10.4)') mZ, amZ, mT, muT, mB, muB, mC, muC, lambda1, lambda2, lam

  read*, ; read*, ; read*, ; print*,

  if ( fit(:4) == 'mass'      ) then
    print*, '  m(m)      sigmaExp   sigmaPert      chi2/dof'
  else if ( fit(:9) == 'alphaMass' ) then
    print*, '  alpha     sigmaExp   sigmaPert      chi2/dof'
  else if ( fit(:5) == 'alpha'     ) then
    print*, '  alpha     sigmaExp   sigmaPert      chi2/dof'
  end if

  print*,

  OPEN (UNIT = 1, FILE = dataFile, ACCESS = 'SEQUENTIAL', STATUS = 'OLD')
  read(1,*) dof; dof = dof - 1; close(2)

  i = 0; maxim = 0; minim = 1000; mean = 0; chi = 0; sigma = 0; meanAlp = 0
  maxAlp = 0; minAlp = 1000; sigmaAlp = 0; meanRho = 0; meanMass2 = 0
  meanAlp2 = 0; meanMassAlp = 0

  do

    read(*,'(A)', iostat = ierror) line;  if (ierror == -1) exit

    if ( fit(:9) /= 'alphaMass' ) then
      read(line,*) R, mu, mass, err, chi2
    else
      read(line,*) R, mu, mass, err, alpha, errAlp, rho, chi2
    end if

    if (R > 4 .or. mu > 4) cycle; i = i + 1

    fact = sqrt(chi2/dof)
    mean = mean + mass; chi = chi + chi2; sigma = sigma + fact * err
    maxim = max(mass,maxim)   ; minim = min(minim, mass)

    if ( fit(:9) == 'alphaMass' ) then
      meanAlp = meanAlp + alpha; sigmaAlp = sigmaAlp + fact * errAlp
      maxAlp = max(alpha,maxAlp); minAlp = min(minAlp, alpha)
      meanRho = meanRho + rho; meanMass2 = meanMass2 + mass**2
      meanAlp2 = meanAlp2 + alpha**2; meanMassAlp = meanMassAlp + mass * alpha
    end if

  end do

  if ( fit(:9) /= 'alphaMass' ) then
    write(*,'(3F11.7,F25.13)') mean/i, sigma/i, (maxim - minim)/2, chi2/dof
  else
    corr = (meanMassAlp - mean * meanAlp/i)/sqrt(meanMass2 - mean**2/i)/&
    sqrt(meanAlp2 - meanAlp**2/i)

    write(*,'(8F14.7,F25.13)') mean/i, sigma/i, (maxim - minim)/2, &
    meanAlp/i, sigmaAlp/i, (maxAlp - minAlp)/2, meanRho/i, corr, chi2/dof

  end if

end program Chi2NRQCDAnalyzer
