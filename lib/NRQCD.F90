module NRQCDClass
  use Constants, only: dp, pi2, d1mach, Pi ;  use RunningClass; use AlphaClass
  use AnomDimClass; use VFNSMSRClass; implicit none ; private

!ccccccccccccccc

  type, public                    :: NRQCD
    private
    real (dp), dimension(0:4,0:3) :: c
    real (dp), dimension(2)       :: cnl, cnf
    real (dp), dimension(0:4)     :: beta
    character (len = 5)           :: scheme
    character (len = 4)           :: up
    character (len = 3)           :: average
    real (dp)                     :: mH, harm, h2, rat, mC
    type (Running)                :: alphaMass
    type (Alpha)                  :: alphaOb
    type (AnomDim)                :: Andim
    type (VFNSMSR)                :: MSR
    integer                       :: n, l, j, s, nl, nf
    integer, dimension(0:3)       :: listFact

  contains

    procedure, pass(self), private :: Binomial, EnInv
    procedure, pass(self), public  :: En, MassFitter, setMass, DeltaCharm, &
    ZeroBin, DeltaCharmBin, MassIter, EnExpand, DeltaCharmBin3, DeltaCharmDer, &
    DeltaCharmExact, DeltaCharmDerBin, MassError, EnError, MassList

  end type NRQCD

!ccccccccccccccc

  interface NRQCD
    module procedure InNRQCD
  end interface NRQCD

  contains

!ccccccccccccccc

  type (NRQCD) function InNRQCD(up, scheme, average, MSR, n, l, j, s)
    integer            , intent(in) :: n
    integer, optional  , intent(in) :: l, j, s
    character (len = *), intent(in) :: scheme, up, average
    type (VFNSMSR)     , intent(in) :: MSR
    real (dp)      , dimension(0:4) :: beta
    real (dp)                       :: h3, c2h, c2val, harm
    character (len = 5)             :: alphaScheme
    integer                         :: nf, nl, i, jj, k, fac

    if ( up(:2) == 'up' ) then
      InNRQCD%alphaMass = MSR%RunArray(2)
    else if ( up(:4) == 'down' ) then
      InNRQCD%alphaMass = MSR%RunArray(1)
    end if

    InNRQCD%n = n; nf = MSR%numFlav() ; alphaScheme = InNRQCD%alphaMass%scheme()
    InNRQCD%Andim = InNRQCD%alphaMass%adim(); InNRQCD%up = up; InNRQCD%l = l
    InNRQCD%mH = MSR%mass(); InNRQCD%nf = nf; InNRQCD%c = 0; InNRQCD%j = j
    InNRQCD%listFact = factList(3); InNRQCD%cnl = 0; InNRQCD%s = s; h3 = 0
    InNRQCD%alphaOb = MSR%AlphaAll(); beta = InNRQCD%Andim%betaQCD('beta')
    InNRQCD%beta = beta;  InNRQCD%average = average; InNRQCD%h2 = 0; c2h = 0
    InNRQCD%cnf = 0

    if ( up(:2) == 'up' ) then
      nl = nf
    else if ( up(:4) == 'down' ) then
      nl = nf - 1
    end if

    InNRQCD%cnl(1) = 103._dp/18 - 5._dp * nf/9
    InNRQCD%cnf(1) = 31._dp/6   - 5._dp * nf/9

    InNRQCD%nl = nl; InNRQCD%MSR = MSR

    if (nf == 5) then
      InNRQCD%rat = InNRQCD%alphaOb%scales('muT')/InNRQCD%mH
    else if (nf == 4) then
      InNRQCD%rat = InNRQCD%alphaOb%scales('muB')/InNRQCD%mH
    end if

    if (  alphaScheme(:4) == 'pole' ) then
      InNRQCD%scheme = 'pole'
    else
      InNRQCD%scheme = scheme
    end if

    InNRQCD%c(:1,0) = [ 1._dp, 31._dp/6 - 5._dp * nl/9 ]

    if ( average(:2) == 'no' ) then

      InNRQCD%c(2:,0) = [ c2(nl, n, l, j, s), c3(nl, n, l, j, s), &
      c3log(nl, n, l, j, s) ]

      InNRQCD%cnl(2) = c2(nf - 1, n, l, j, s)
      InNRQCD%cnf(2) = c2(nf    , n, l, j, s)

      InNRQCD%harm = Harmonic(n + l)

    else if ( average(:1) == 'j' ) then

      InNRQCD%harm = Harmonic(n + l)

      do jj = abs(l - s), l + s

        fac = 2 * jj + 1

        InNRQCD%c(2:,0) = InNRQCD%c(2:,0) + fac * [ c2(nl, n, l, jj, s), &
        c3(nl, n, l, jj, s), c3log(nl, n, l, jj, s) ]

        InNRQCD%cnl(2) = InNRQCD%cnl(2) + fac * c2(nf - 1, n, l, jj, s)
        InNRQCD%cnf(2) = InNRQCD%cnf(2) + fac * c2(nf    , n, l, jj, s)

      end do

      fac = (2 * l + 1) * (2 * s + 1)

      InNRQCD%c(2:,0) = InNRQCD%c(2:,0)/fac
      InNRQCD%cnf(2)  = InNRQCD%cnf(2)/fac
      InNRQCD%cnl(2)  = InNRQCD%cnl(2)/fac

    else if ( average(:3) == 'yes' ) then

      InNRQCD%harm = 0.5_dp - 1._dp/n/2 + Harmonic(n)

      do i = 0, n - 1
        harm = Harmonic(n + i)
        InNRQCD%h2 = InNRQCD%h2 + (2*i + 1) * harm**2
        h3         = h3         + (2*i + 1) * harm**3
      end do

      InNRQCD%h2 = InNRQCD%h2/n**2 - InNRQCD%harm**2
      h3         = h3/n**2         - InNRQCD%harm**3 - 3 * InNRQCD%harm * InNRQCD%h2

      do i = 0, n - 1

        harm = Harmonic(n + i)

        do k = 0, 1
          do jj = abs(i - k), i + k

            fac = 2 * jj + 1; c2val = c2(nl, n, i, jj, k)

            InNRQCD%c(2:,0) = InNRQCD%c(2:,0) + fac * [ c2val, &
            c3(nl, n, i, jj, k), c3log(nl, n, i, jj, k) ]

            InNRQCD%cnl(2) = InNRQCD%cnl(2) + fac * c2(nf - 1, n, i, jj, k)
            InNRQCD%cnf(2) = InNRQCD%cnf(2) + fac * c2(nf    , n, i, jj, k)

            c2h = c2h + fac * harm * c2val

          end do
        end do
      end do

      fac = 4 * n**2

      InNRQCD%c(2:,0) = InNRQCD%c(2:,0)/fac
      InNRQCD%cnf(2)  = InNRQCD%cnf(2)/fac
      InNRQCD%cnl(2)  = InNRQCD%cnl(2)/fac

      c2h = c2h/fac - InNRQCD%harm * InNRQCD%c(2,0)

      ! InNRQCD%c(3,0) = InNRQCD%c(3,0) + 2 * beta(0) * c2h
      ! InNRQCD%c(2,0) = InNRQCD%c(2,0) + 3 * beta(0)**2 * InNRQCD%h2/4

    end if

    do k = 1, 3
      do jj = 0, k - 1

        do i = jj + 1, k - 1

          InNRQCD%c(k,jj + 1) = InNRQCD%c(k,jj + 1) + beta(k - 1 - i) * &
          ( (i + 2) * InNRQCD%c(i, jj) - (jj + 1) * InNRQCD%c(i, jj + 1) )/4**(k - i)

        end do

        InNRQCD%c(k,jj + 1) = 2 * (  InNRQCD%c(k,jj + 1) + &
        (jj + 2) * InNRQCD%c(jj,jj) * beta(k - 1 - jj)/4**(k - jj)  )/(jj + 1)

      end do
    end do

    InNRQCD%c(3,0) = InNRQCD%c(3,0) + InNRQCD%c(3,2) * InNRQCD%h2 + &
    InNRQCD%c(3,3) * h3 + 2 * beta(0) * c2h

    InNRQCD%c(3,1) = InNRQCD%c(3,1) + 3 * InNRQCD%c(3,3) * InNRQCD%h2
    InNRQCD%c(2,0) = InNRQCD%c(2,0) +     InNRQCD%c(2,2) * InNRQCD%h2

    if (nf == 5) InNRQCD%mc = InNRQCD%alphaOb%scales('mB')
    if (nf == 4) InNRQCD%mc = InNRQCD%alphaOb%scales('mC')

  end function InNRQCD

!ccccccccccccccc

  subroutine setMass(self, m, mu)
    class (NRQCD), intent(inout) :: self
    real (dp)    , intent(in)    :: m, mu

    self%mH = m;    call self%MSR%setMass(m, mu)

    if (self%nf == 5) then
      call self%alphaMass%SetMTop(m, mu)
    else if (self%nf == 4) then
      call self%alphaMass%SetMBottom(m, mu)
    end if

  end subroutine setMass

!ccccccccccccccc

  function MassList(self, iter, charm, n, order, mu0, mu1, deltaMu, R0, R1, &
  deltaR, mUpsilon, lambda, method) result(list)
    class (NRQCD)      , intent(inout) :: self
    character (len = *), intent(in)    :: method, charm, iter
    integer            , intent(in)    :: order, n
    real (dp)          , intent(in)    :: lambda, mUpsilon, mu0, mu1, R1, &
    deltaMu, R0, deltaR

    real (dp), dimension( 3, 0:Floor( (mu1 - mu0)/deltaMu ), &
    0:Floor( (R1 - R0)/deltaR ))    :: list

    integer                            :: imax, jmax, i, j

    imax = Floor( (mu1 - mu0)/deltaMu ); jmax = Floor( (R1 - R0)/deltaR )

    list = 0

    do i = 0, imax

      list(1,i,:) = mu0 + i * deltaMu

      do j = 0, jmax
        list(2,i,j) = R0 + j * deltaR
        list(3,i,j) = self%MassFitter(iter, charm, n, order, list(1,i,j), &
        list(2,i,j), mUpsilon, lambda, method)
      end do

    end do

  end function MassList

!ccccccccccccccc

  ! function NRQCDList(self, iter, charm, n, order, mu0, mu1, deltaMu, R0, R1, &
  ! deltaR, mUpsilon, lambda, method) result(list)
  !   class (NRQCD)      , intent(inout) :: self
  !   character (len = *), intent(in)    :: method, charm, iter
  !   integer            , intent(in)    :: order, n
  !   real (dp)          , intent(in)    :: lambda, mUpsilon, mu0, mu1, R1, &
  !   deltaMu, R0, deltaR
  !
  !   real (dp), dimension( 3, 0:Floor( (mu1 - mu0)/deltaMu ), &
  !   0:Floor( (R1 - R0)/deltaR ))    :: list
  !
  !   integer                            :: imax, jmax, i, j
  !
  !   imax = Floor( (mu1 - mu0)/deltaMu ); jmax = Floor( (R1 - R0)/deltaR )
  !
  !   list = 0
  !
  !   do i = 0, imax
  !
  !     list(1,i,:) = mu0 + i * deltaMu
  !
  !     do j = 0, jmax
  !       list(2,i,j) = R0 + j * deltaR
  !       list(3,i,j) = self%MassFitter(iter, charm, n, order, list(1,i,j), &
  !       list(2,i,j), mUpsilon, lambda, method)
  !     end do
  !
  !   end do
  !
  ! end function NRQCDList

!ccccccccccccccc

  function MassError(self, iter, charm, n, order, mu0, mu1, deltaMu, R0, R1, &
  deltaR, x, mUpsilon, lambda, method) result(list)
    class (NRQCD)      , intent(inout) :: self
    character (len = *), intent(in)    :: method, charm, iter
    integer            , intent(in)    :: order, n
    real (dp)          , intent(in)    :: lambda, mUpsilon, mu0, mu1, R1, x, &
    deltaMu, R0, deltaR
    real (dp), dimension(2)            :: list
    real (dp)                          :: mu, R, mass, xinv, rat
    integer                            :: imax, jmax, i, j

    list = [0._dp, mUpsilon]; xinv = 1/x

    imax = Floor( (mu1 - mu0)/deltaMu ); jmax = Floor( (R1 - R0)/deltaR )

    do i = 0, imax

      mu = mu0 + i * deltaMu

      if ( self%scheme(:3) == 'MSR' ) then

        do j = 0, jmax
          R = R0 + j * deltaR
          rat = mu/R; if ( rat > x .or.  rat < xinv ) cycle
          mass = self%MassFitter(iter, charm, n, order, mu, R, mUpsilon, &
          lambda, method)
          if ( mass > list(1) ) list(1) = mass
          if ( mass < list(2) ) list(2) = mass
        end do

      else

        mass = self%MassFitter(iter, charm, n, order, mu, mu, mUpsilon, &
        lambda, method)
        if ( mass > list(1) ) list(1) = mass
        if ( mass < list(2) ) list(2) = mass

      end if

    end do

    list = [ sum(list), list(1) - list(2) ]/2

  end function MassError

!ccccccccccccccc

  function EnError(self, iter, charm, n, mu0, mu1, deltaMu, R0, R1, &
  deltaR, x, mUpsilon, lambda, method) result(list)
    class (NRQCD)      , intent(inout) :: self
    character (len = *), intent(in)    :: method, charm, iter
    integer            , intent(in)    :: n
    real (dp)          , intent(in)    :: lambda, mu0, mu1, R1, x, deltaR, R0, &
    deltaMu, mUpsilon
    real (dp), dimension(2,0:4)        :: list
    real (dp)                          :: mu, R, xinv, rat, upMass
    real (dp), dimension(0:4)          :: mass
    integer                            :: imax, jmax, i, j, k

    list(1,:) = 0; list(2,:) = 3 * self%mH; xinv = 1/x; mass = 0

    imax = Floor( (mu1 - mu0)/deltaMu ); jmax = Floor( (R1 - R0)/deltaR )

    do i = 0, imax

      mu = mu0 + i * deltaMu

      if ( self%scheme(:3) == 'MSR' ) then

        do j = 0, jmax

          R = R0 + j * deltaR; upMass = 0
          rat = mu/R; if ( rat > x .or.  rat < xinv ) cycle

          if ( iter(:10) == 'FixedOrder') then
            mass = self%En(charm, n, mu, R, lambda, method)
          else if ( iter(:8) == 'expanded') then
            mass = self%EnExpand(charm, n, mu, R, mUpsilon, lambda, method)
          else if ( iter(:9) == 'iterative') then
            mass = self%MassIter(charm, n, mu, R, mUpsilon, lambda, method)
          end if

          do k = 0, 4
            upMass = upMass + mass(k)
            if ( upMass > list(1,k) ) list(1,k) = upMass
            if ( upMass < list(2,k) ) list(2,k) = upMass
          end do

        end do

      else

        mass = self%En(charm, n, mu, mu, lambda, method); upMass = 0

        do k = 0, 4
          upMass = upMass + mass(k)
          if ( upMass > list(1,k) ) list(1,k) = upMass
          if ( upMass < list(2,k) ) list(2,k) = upMass
        end do

      end if

    end do

    do k = 0, 4
      list(:,k) = [  sum( list(:,k) ), list(1,k) - list(2,k)  ]/2
    end do

  end function EnError

!ccccccccccccccc

  real (dp) function MassFitter(self, iter, charm, n, order, mu, R, mUpsilon, &
  lambda, method)
    class (NRQCD)      , intent(inout) :: self
    character (len = *), intent(in)    :: method, charm, iter
    integer            , intent(in)    :: order, n
    real (dp)          , intent(in)    :: mu, R, lambda, mUpsilon
    real (dp)                          :: a

    MassFitter = FindRoot(mUpsilon/2)

    do
      a = FindRoot(MassFitter);  if ( abs(a - MassFitter) < 1e-10_dp ) exit
      MassFitter = a
    end do

  contains

    real (dp) function FindRoot(mass)
      real (dp), intent(in)     :: mass
      real (dp), dimension(0:7) :: list

      call self%SetMass(mass, self%rat * mass)

      if ( iter(:9) == 'iterative' ) then
        list(:4) = self%MassIter(charm, order, mu, R, mUpsilon, lambda, method)
        FindRoot = sum( list(:n) )
      else if ( iter(:10) == 'FixedOrder' ) then
        list = self%EnInv(charm, order, mu, R, lambda, method)
        list(2:3) = list(2:3) + list(6:7)
        FindRoot = mUpsilon/sum( list(:n) )/2 - list(5)
      else if ( iter(:8) == 'expanded' ) then
        list(:4) = self%EnExpand(charm, order, mu, R, mUpsilon, lambda, method)
        FindRoot = sum( list(:n) )
      end if

    end function FindRoot

  end function MassFitter

!ccccccccccccccc

  function En(self, charm, order, mu, R, lambda, method) result(list)
    class (NRQCD)   , intent(inout) :: self
    character (len = *), intent(in) :: method, charm
    integer            , intent(in) :: order
    real (dp)          , intent(in) :: mu, R, lambda
    real (dp), dimension(0:4)       :: list
    real (dp), dimension(0:7)       :: listA

    listA = self%EnInv(charm, order, mu, R, lambda, method)
    listA(2) = listA(2) + listA(6)
    listA(3) = listA(3) + listA(7)
    list  = 2 * ( self%mH + listA(5) ) * listA(:4)

  end function En

!ccccccccccccccc

  function EnExpand(self, charm, order, mu, R, mUpsilon, lambda, method) result(list)
    class (NRQCD)   , intent(inout) :: self
    character (len = *), intent(in) :: method, charm
    integer            , intent(in) :: order
    real (dp)          , intent(in) :: mu, R, lambda, mUpsilon
    real (dp), dimension(0:4)       :: list
    real (dp), dimension(0:7)       :: listA
    integer                         :: n

    listA = self%EnInv(charm, order, mu, R, lambda, method)

    list(0) = 1

    do n = 0, 3
      list(n + 1) = - sum( list(:n) * listA(n + 1:1:-1) )
    end do

    list(2) = list(2) - listA(6)
    list(3) = list(3) - 2 * listA(6) * list(1) - listA(7)

    list = mUpsilon * list/2

    list(0) = list(0) - listA(5)

  end function EnExpand

!ccccccccccccccc

  function EnInv(self, charm, order, mu, R, lambda, method) result(list)
    class (NRQCD)   , intent(inout) :: self
    character (len = *), intent(in) :: method, charm
    integer            , intent(in) :: order
    real (dp)          , intent(in) :: mu, R, lambda
    real (dp), dimension(0:4)       :: alphaList
    real (dp), dimension(0:7)       :: list
    real (dp), dimension(0:3)       :: logList
    real (dp), dimension(0:4)       :: delta
    real (dp), dimension(4)         :: lgmList
    real (dp), dimension(0:4,0:4)   :: deltaLog  ! (power, order)
    real (dp), dimension(0:4,4)     :: coefMSR
    real (dp), dimension(4,0:3)     :: c
    real (dp), dimension(2)         :: deltaM, deltaCharm
    integer                         :: i, j, k, l
    real (dp)                       :: alp, Rmass, mass, factor, lg, rat, delta2

    list = 0; list(0) = 1 ; alp = self%alphaMass%alphaQCD(mu); coefMSR = 0
    alphaList(0) = 1; alphaList(1:) = PowList(alp/Pi,4); delta(0) = 1; delta2 = 0
    factor = - 2 * alp**2/9/self%n**2 ; logList(0) = 1 ; deltaM = 0

    if ( self%scheme(:4) == 'pole' ) then
      delta(1:) = 0; Rmass = 0; mass = self%mH
    else if ( self%scheme(:5) == 'MSbar' ) then
      coefMSR(0,:) = self%mH * self%Andim%MSRDelta('MSRp')
      Rmass = self%mH; mass = self%mH
    else if ( self%scheme(:3) == 'MSR' ) then
      coefMSR(0,:) = R * self%andim%MSRDelta(self%scheme); Rmass = R
      mass = self%MSR%MSRmass(self%up, self%scheme, order, R, lambda, method)
    end if

    logList(1:) = PowList( log(3 * self%n * mu / 4 / alp / mass) + self%harm, 3 )

    if ( self%scheme(:4) == 'pole' ) then

      list(1:4) = matmul( self%c(:3,:), logList )
      list(4) = list(4) + self%c(4,0) * log(alp)
      list(1:4) = factor * alphaList(:3) * list(1:4)

    else

      call self%andim%expandAlpha(coefMSR); lgmList = PowList( log(mu/Rmass), 4 )
      call AddAlpha( coefMSR, alphaList(1:) )     ; c = 0; deltaLog = 0
      delta(1:) = DeltaComputer(coefMSR, lgmList, 0)/mass; deltaLog(0,0) = 1
      deltaLog(1,1:2) = delta(1:2) - [ 0._dp, delta(1)**2/2 ]

      if ( self%scheme(:3) == 'MSR' .or. self%up(:2) == 'up' ) then

        deltaM(1) = self%MSR%DeltaM(self%up, Rmass)
        deltaM(2) = self%MSR%DeltaM2(self%scheme, self%up, Rmass) + &
        self%beta(0) * log(mu/Rmass) * deltaM(1)

      else if ( self%scheme(:5) == 'MSbar' .and. self%up(:4) == 'down' ) then

        rat = self%mC/self%mH; lg = -log(rat);  deltaM(1) = deltaCharm2(rat)

        deltaM(2) = deltaCharm3(self%nf, 1, rat) + 2 * lg * deltaM(1)/3 +       &
        4 * lg**2/27 - 27.51152614489051_dp + 1.3053814981630874_dp * self%nf + &
        lg * (11.073375282398949_dp - 0.694244607447754_dp * self%nf)

        deltaM(1) = 4 * log(rat)/9 + deltaM(1) - 71._dp/144 - pi2/18
        deltaM(2) = deltaM(2) + self%beta(0) * log(mu/self%mH) * deltaM(1)

      end if

      deltaM = Rmass * alphaList(2:3) * deltaM /mass

      do i = 1, 4
        do k = 0, i - 1

          do j = k, i - 1
            c(i,k) = c(i,k) + self%c(i - 1,j) * self%Binomial(j,k) * logList(j - k)
          end do

          c(i,k) = (-1)**k * c(i,k)
          if (i == 4 .and. k == 0) c(i,k) = c(i,k) + self%c(4,0) * log(alp)

        end do

        c(i,:) = factor * alphaList(i - 1) * c(i,:)

      end do

      do j = 1, 4
        do i = 1, j
          l = min(i - 1, j - i)
          list(j) = list(j) + sum( deltaLog(:l,j - i) * c(i,:l) )
        end do
      end do

      do j = 4, 1, -1
        list(j) = sum( delta(:j) * list(j:0:-1) )
      end do

    end if

    deltaCharm = factor * alphaList(1:2) * [ self%DeltaCharmBin(alp, mass), &
    self%DeltaCharmExact('exact', mu, alp, mass) ]

    list(6:7) = deltaM + deltaCharm(1:2)

    list(7) = list(7) + factor * deltaM(1) + deltaCharm(1) * delta(1) - &
    factor * alphaList(1) * delta(1) * self%DeltaCharmDerBin(alp, mass)

    list(5) = mass - self%mH

  end function EnInv

!ccccccccccccccc

  function MassIter(self, charm, order, mu, R, mUpsilon, lambda, method) result(list)
    class (NRQCD)   , intent(inout) :: self
    character (len = *), intent(in) :: method, charm
    integer            , intent(in) :: order
    real (dp)          , intent(in) :: mu, R, lambda, mUpsilon
    real (dp), dimension(0:4)       :: list, listInv, alphaList
    real (dp), dimension(0:3)       :: logList
    real (dp), dimension(0:4)       :: delta, deltaInv
    real (dp), dimension(4)         :: lgmList
    real (dp), dimension(0:4,4)     :: coefMSR
    real (dp), dimension(2)         :: deltaM, deltaCharm
    integer                         :: n
    real (dp)                       :: alp, Rmass, mass, factor, mTree, lg, rat

    list = 0; list(0) = 1 ; alp = self%alphaMass%alphaQCD(mu); coefMSR = 0
    alphaList(0) = 1; alphaList(1:) = PowList(alp/Pi,4); delta(0) = 1
    factor = - 2 * alp**2/9/self%n**2 ; logList(0) = 1; mTree = mUpsilon/2
    deltaM = 0; listInv = 0; listInv(0) = 1

    if ( self%scheme(:4) == 'pole' ) then
      delta(1:) = 0; Rmass = 0; mass = self%mH
    else if ( self%scheme(:5) == 'MSbar' ) then
      coefMSR(0,:) = self%Andim%MSRDelta('MSRp')
      Rmass = mTree; mass = self%mH
    else if ( self%scheme(:3) == 'MSR' ) then
      coefMSR(0,:) = self%andim%MSRDelta(self%scheme); Rmass = R
      mass = self%MSR%MSRmass(self%up, self%scheme, order, R, lambda, method)
    end if

    logList(1:) = PowList( log(3 * self%n * mu / 4 / alp / mTree) + self%harm, 3 )

    listInv(1:) = matmul( self%c(:3,:), logList )
    listInv(4) = listInv(4) + self%c(4,0) * log(alp)
    listInv(1:) = factor * alphaList(:3) * listInv(1:)

    listInv(3) = listInv(3) + alphaList(1) * self%beta(0) * factor**2
    listInv(4) = listInv(4) - alphaList(1) * self%beta(0) * factor**3/2 + &
    alphaList(2) * factor**2 * (  5 * self%beta(0)**2 * logList(1)/2 + &
    ( 10 * self%beta(0) * self%c(1,0) - 2 * self%beta(0)**2 + self%beta(1) )/4  )

    do n = 0, 3
      list(n + 1) = - sum( list(:n) * listInv(n + 1:1:-1) )
    end do

    if ( self%scheme(:4) /= 'pole' ) then

      if ( self%scheme(:3) == 'MSR' .or. self%up(:2) == 'up' ) then
        deltaM(1) = self%MSR%DeltaM(self%up, Rmass)
        deltaM(2) = self%MSR%DeltaM2(self%scheme, self%up, Rmass) + &
        self%beta(0) * log(mu/Rmass) * deltaM(1)
      else if ( self%scheme(:5) == 'MSbar' .and. self%up(:4) == 'down' ) then

        rat = self%mC/mTree; lg = - log(rat);  deltaM(1) = deltaCharm2(rat)

        deltaM(2) = deltaCharm3(self%nf, 1, rat) + 2 * lg * deltaM(1)/3 +       &
        4 * lg**2/27 - 27.51152614489051_dp + 1.3053814981630874_dp * self%nf + &
        lg * (11.073375282398949_dp - 0.694244607447754_dp * self%nf)

        deltaM(1) = 4 * log(rat)/9 + deltaM(1) - 71._dp/144 - pi2/18
        deltaM(2) = deltaM(2) + self%beta(0) * log(mu/mTree) * deltaM(1)

      end if

      deltaM = Rmass * alphaList(2:3) * deltaM/mTree

      call self%andim%expandAlpha(coefMSR); lgmList = PowList( log(mu/Rmass), 4 )
      call AddAlpha( coefMSR, alphaList(1:) )
      delta(1:) = DeltaComputer(coefMSR, lgmList, 0)

      if ( self%scheme(:3) == 'MSR' ) then
        list(1:) = list(1:) - delta(1:) * Rmass/mTree
      else

        delta(3) = delta(3) + coefMSR(1,2) * ( delta(1) - list(1) )

        delta(4) = delta(4) + coefMSR(1,2) * ( list(1)**2 - 2 * list(2) - &
        delta(1)**2 + 2 * coefMSR(0,2) )/2 + coefMSR(1,3) * ( delta(1) &
         - list(1) ) + lgmList(1) * (  coefMSR(1,2)**2 + 2 * coefMSR(2,3) * &
        ( delta(1) - list(1) )  )

        deltaInv(0) = 1

        do n = 0, 3
          deltaInv(n + 1) = - sum( deltaInv(:n) * delta(n + 1:1:-1) )
        end do

        do n = 4, 1, -1
          list(n) = sum( deltaInv(:n) * list(n:0:-1) )
        end do

      end if

    end if

    deltaCharm = factor * alphaList(1:2) * [ self%DeltaCharmBin(alp, mTree), &
    self%DeltaCharmExact('exact', mu, alp, mTree) ]

    list(2:3) = list(2:3) - deltaM - deltaCharm

    list(3) = list(3) + factor**2 * alphaList(1) * &
    ( 2 * self%DeltaCharmBin(alp, mTree) - self%DeltaCharmDerBin(alp, mass) )

    if ( self%up(:2) == 'up' .and. self%mC > tiny(1._dp) ) then

      if ( self%scheme(:5) == 'MSbar' ) then

        rat = self%mC/mTree

        list(3) = list(3) + deltaM(1) * ( 2 * delta(1) - list(1) ) + &
        alphaList(2) * ( list(1) - delta(1) ) * rat * DeltaCharm2Der(rat)

        if ( self%up(:4) == 'down' ) list(3) = list(3) + 4._dp/9

      end if

    end if

    list = mTree * list; list(0) = list(0) + self%mH - mass

  end function MassIter

!ccccccccccccccc

  real (dp) function DeltaCharmBin3(self, mu, ln, mC) result(delta2)
    class (NRQCD), intent(in) :: self
    real (dp)    , intent(in) :: mu, ln, mC
    real (dp)                 :: lg

   lg = log(mu/mC)

    delta2 = self%cnl(2) + lg**2/3 - 2 * (ln**2 + self%h2) * (self%nf - 17)/3 + lg * &
    ( 13._dp/18 - 2._dp * self%nf/9 - self%cnl(1) ) + ln * ( lg * &
    (2 * self%nf - 35)/3 + (8 * self%nf - 79)/18._dp + (35/2._dp - self%nf) * &
    self%cnl(1) + (self%nf - 33/2._dp) * self%cnf(1) ) - self%cnf(2)

    if ( self%scheme(:4) == 'pole' ) then
      delta2 = delta2 - 7._dp/12
    else
      delta2 = delta2 + 11._dp/36
    end if

  end function DeltaCharmBin3

!ccccccccccccccc

  real (dp) function DeltaCharmExact(self, type, mu, alp, mass)
    class (NRQCD)   , intent(inout) :: self
    character (len = *), intent(in) :: type
    real (dp)          , intent(in) :: mu, alp, mass
    real (dp)                       :: lg, x, App1, App2, App3, gamma, beta, a, &
    b, e1, e2, mC
    real (dp), parameter            :: c1 = - 0.8316040946513316_dp, &
    c2 = 0.470_dp, d2 = 1.120_dp, d1 = 1.8316040946513317_dp

    DeltaCharmExact = 0; if ( self%mC <= tiny(1._dp) ) return

    mC = 1._dp/sqrt(1._dp * self%n)

    if ( type(:5) == 'exact' ) then

      gamma = 2 * mass * alp/3; x = self%mC/gamma

      if ( self%n == 1 ) then

        if (x < 7) then

          lg = log(x);  beta = 11 - 2 * self%nf/3._dp

          App1 = - 24.865811264136_dp * x + 11.983864345921031_dp * x**2 - &
          26.71170106052555_dp * x**3 + 15.542365205643007_dp * x**4 + &
          7.9194151299077715_dp * x**5 + 0.3053516028457726_dp * x**6 + &
          0.04307631370176248_dp * x * lg + 0.23973826692206646_dp * x**2 * lg - &
          12.19353255229631_dp * x**3 * lg - 17.71098150670699_dp * x**4 * lg - &
          2.864456892317382_dp * x**5 * lg - 0.05715945405398001_dp * x**6 * lg

          App2 = 1.9477460481502156_dp * x + 0.1527091423940376_dp * x**2 - &
          2.862610556473993_dp * x**3 + 2.474572078659877_dp * x**4 + &
          0.8542022081305681_dp * x**5 + 0.02932979796074642_dp * x**6 - &
          2.3473257597558925_dp * x * lg + 0.32395085128272094_dp * x**2 * lg - &
          2.4754591830262753_dp * x**3 * lg - 2.181651709563445_dp * x**4 * lg - &
          0.2964274725676275_dp * x**5 * lg - 0.005426744077104387_dp * x**6 * lg

          App3 = 1.553742421397914_dp * x - 0.8004938532161119_dp * x**2 + &
          1.380778522466775_dp * x**3 - 0.6923566636506896_dp * x**4 - &
          0.3823712766542248_dp * x**5 - 0.014981803878737314_dp * x**6 - &
          0.003436674677916459_dp * x * lg - 0.06597592537344886_dp * x**2 * lg + &
          0.4712708827971536_dp * x**3 * lg + 0.8364190131891226_dp * x**4 * lg + &
          0.13911559154554373_dp * x**5 * lg + 0.0028092361619528835_dp * x**6 * lg

          DeltaCharmExact = 2 * (  App1 + 57 * ( c1 * c2 * x/(1 + c2 * x) + &
          d1 * d2 * x/(1 + d2 * x) + c1 * log(1 + c2 * x) + d1 * log(1 + d2 * x) )/4 &
          + 3 * beta * ( App2 + 3 * App3 * ( log(mu/2/gamma) + 5/6._dp )/2 )  )/9

          if ( self%scheme(:4) /= 'pole' ) then
            DeltaCharmExact = DeltaCharmExact + 4 * self%DeltaCharmDerBin(alp, mass)/3
          end if

        else

          lg = log(3 * self%n * mu / 4 / alp / mass) + self%harm
          DeltaCharmExact = self%DeltaCharmBin3(mu, lg, self%mC)

        end if

      else

        lg = log(3 * self%n * mu / 4 / alp / mass) + self%harm

        if ( self%mc < mC) then

          e1 = self%DeltaCharmBin3( mu, lg, mC )

          e2 = ( self%DeltaCharmBin3(mu, lg, mC + 0.01_dp) - &
          self%DeltaCharmBin3(mu, lg, mC - 0.01_dp) )/0.02_dp

          a = ( e1 + (e1 - e2 * mC) * Log(mC) )/mC; b = e2 - e1/mC

          DeltaCharmExact = self%mC * ( a + b * log(self%mC) )

        else

          DeltaCharmExact = self%DeltaCharmBin3(mu, lg, self%mC)

        end if

      end if

    else

      lg = log(3 * self%n * mu / 4 / alp / mass) + self%harm

      DeltaCharmExact = self%DeltaCharmBin3(mu, lg, self%mC)

    end if

    if ( self%up(:4) == 'down' ) then
      lg = log(3 * self%n * mu / 4 / alp / mass) + self%harm
      DeltaCharmExact = DeltaCharmExact - self%DeltaCharmBin3(mu, lg, self%mC)
    end if

  end function DeltaCharmExact

!ccccccccccccccc

  pure real (dp) function c2(nl, n, l, j, s)
    integer, intent(in) :: nl, n, l, j, s
    integer             :: nl2, ss, ss2, jj, jj2

    c2 = 0; nl2 = nl**2; ss = s * (1 + s); ss2 = ss**2; jj = j * (j + 1)
    jj2 = jj**2

    if (n == 1) then

      if (l == 0) then
        c2 = 220.34174177286914_dp - 19.798187714015196_dp * nl + &
        0.496190504426009_dp * nl2 - 11.69730891980961_dp * ss
      end if

    else if (n == 2) then

      if (l == 0) then
        c2 = 234.57814959979672_dp - 24.946605003852223_dp * nl + &
        0.6522031495725855_dp * nl2 - 5.848654459904805_dp * ss
      else if (l == 1) then
        c2 = 146.39649245990387_dp - 1.7545963379714415_dp * jj   + &
        0.2193245422464302_dp  * jj2 - 16.951844980548632_dp * nl + &
        0.40993769432096167_dp * nl2 + (1.169730891980961_dp      - &
        0.4386490844928604_dp * jj) * ss + 0.2193245422464302_dp * ss2
      end if

    else if (n == 3) then

      if (l == 0) then
        c2 = 243.61209732248133_dp - 27.345022293689244_dp * nl + &
        0.7248824613858287_dp * nl2 - 3.8991029732698697_dp * ss
      else if (l == 1) then
        c2 = 170.30854367434102_dp - 1.169730891980961_dp * jj    +  &
        0.14621636149762013_dp * jj2 - 20.255693524026842_dp * nl +  &
        0.5100543168506044_dp * nl2 + (0.7798205946539739_dp      -  &
        0.29243272299524026_dp * jj) * ss + 0.14621636149762013_dp * ss2
      else if (l == 2) then
        c2 = 132.27874873495742_dp - 0.2228058841868497_dp * jj     + &
        0.0069626838808390535_dp * jj2 - 16.128767492567004_dp * nl + &
        0.3849959522609123_dp * nl2 + (0.1671044131401373_dp        - &
        0.013925367761678107_dp * jj) * ss + 0.0069626838808390535_dp * ss2
      end if

    else if (n == 4) then
      if (l == 0) then
        c2 = 249.6485479478922_dp - 28.758871682291705_dp * nl + &
        0.7677263822525698_dp * nl2 - 2.9243272299524024_dp * ss
      else if (l == 1) then
        c2 = 186.45187350545456_dp - 0.8772981689857208_dp * jj  + &
        0.1096622711232151_dp * jj2 - 22.445631564883403_dp * nl + &
        0.5764160756644395_dp * nl2 + (0.5848654459904805_dp     - &
        0.2193245422464302_dp * jj) * ss + 0.1096622711232151_dp * ss2
      else if (l == 2) then
        c2 = 150.96237118757924_dp - 0.1671044131401373_dp * jj    + &
        0.0052220129106292906_dp * jj2 - 18.50593326997487_dp * nl + &
        0.45703127884902944_dp * nl2 + (0.12532830985510296_dp     - &
        0.010444025821258581_dp * jj) * ss + 0.0052220129106292906_dp * ss2
      else if (l == 3) then
        c2 = 126.83861491499124_dp - 0.059182812987131954_dp * jj   + &
        0.0008703354851048817_dp * jj2 - 15.738818473520091_dp * nl + &
        0.37317931532009674_dp * nl2 + (0.04525744522545385_dp      - &
        0.0017406709702097634_dp * jj) * ss + 0.0008703354851048817_dp * ss2
      end if

    end if

  end function c2

!ccccccccccccccc

  pure real (dp) function c3log(nl, n, l, j, s)
    integer, intent(in) :: nl, n, l, j, s
    integer             :: nl2, ss, ss2, jj, jj2

    c3log = 0; nl2 = nl**2; ss = s * (1 + s); ss2 = ss**2; jj = j * (j + 1)
    jj2 = jj**2

    if (n == 1) then

      if (l == 0) then
        c3log = 597.1110662659062_dp - 61.41087182900045_dp * ss
      end if

    else if (n == 2) then

      if (l == 0) then
        c3log = 329.5351247252613_dp - 30.705435914500224_dp * ss
      else if (l == 1) then
        c3log = 114.45085696226215_dp - 4.1671663026821735_dp * jj + &
        0.6579736267392905_dp * jj2 + (2.4125699647107317_dp     - &
        1.315947253478581_dp  * jj) * ss + 0.6579736267392905_dp * ss2
      end if

    else if (n == 3) then

      if (l == 0) then
        c3log = 236.44404123844322_dp - 20.470290609666815_dp * ss
      else if (l == 1) then
        c3log = 93.05452939644374_dp - 2.7781108684547826_dp * jj + &
        0.4386490844928604_dp * jj2 + (1.6083799764738214_dp    - &
        0.8772981689857208_dp * jj) * ss + 0.4386490844928604_dp * ss2
      else if (l == 2) then
        c3log = 72.13862701840323_dp - 0.522201291062929_dp * jj + &
        0.020888051642517162_dp * jj2 + (0.3550968779227917_dp - &
        0.041776103285034324_dp * jj) * ss + 0.020888051642517162_dp * ss2
      end if

    else if (n == 4) then
      if (l == 0) then
        c3log = 189.16741768754605_dp - 15.352717957250112_dp * ss
      else if (l == 1) then
        c3log = 81.62528380604644_dp - 2.0835831513410867_dp * jj + &
        0.32898681336964525_dp * jj2 + (1.2062849823553659_dp   - &
        0.6579736267392905_dp * jj) * ss + 0.32898681336964525_dp * ss2
      else if (l == 2) then
        c3log = 65.93835702251604_dp - 0.3916509682971967_dp * jj + &
        0.01566603873188787_dp * jj2 + (0.2663226584420938_dp   - &
        0.03133207746377574_dp * jj) * ss + 0.01566603873188787_dp * ss2
      else if (l == 3) then
        c3log = 59.17062829034048_dp - 0.1383833421316762_dp * jj  + &
        0.0026110064553146453_dp * jj2 + (0.09660723884664187_dp - &
        0.0052220129106292906_dp * jj) * ss + 0.0026110064553146453_dp * ss2
      end if

    end if

  end function c3log

!ccccccccccccccc

  pure real (dp) function c3(nl, n, l, j, s)
    integer, intent(in) :: nl, n, l, j, s
    integer             :: nl3, nl2, ss, ss2, jj, jj2

    c3 = 0; nl2 = nl**2; ss = s * (1 + s); ss2 = ss**2; jj = j * (j + 1)
    jj2 = jj**2; nl3 = nl**3

    if (n == 1) then

      if (l == 0) then
        c3 = 1928.7628660728874_dp - 418.00274061279185_dp * nl   + &
        27.350849570543684_dp * nl2 - 0.4478789305261393_dp * nl3 + &
        (218.58922102668993_dp - 11.52783363067168_dp * nl) * ss
      end if

    else if (n == 2) then

      if (l == 0) then
        c3 = 1555.6647767522752_dp - 427.2864081082992_dp * nl    + &
        29.077705021060194_dp * nl2 - 0.4700406103132774_dp * nl3 + &
        (189.25038089189587_dp - 10.715520511240458_dp * nl) * ss
      else if (l == 1) then
        c3 = 1919.4334788948734_dp + jj * (32.43347193767492_dp    - &
        1.224489010765641_dp * nl) + jj2 * (-3.9574555630864614_dp + &
        0.16829199733504055_dp * nl) - 412.57539959999093_dp * nl  + &
        25.345102512389868_dp * nl2 - 0.41382274015759046_dp * nl3 + &
        (-22.002104070692376_dp + jj * (7.914911126172923_dp       - &
        0.3365839946700811_dp * nl) + 0.9381729750917773_dp * nl) * ss + &
        (-3.9574555630864614_dp + 0.16829199733504055_dp * nl) * ss2
      end if

    else if (n == 3) then

      if (l == 0) then
        c3 = 1419.3458052875883_dp - 418.4769180394629_dp * nl   + &
        28.6078892939055_dp * nl2 - 0.45420057781924034_dp * nl3 + &
        (158.96019573759003_dp - 9.145048480340092_dp * nl) * ss
      else if (l == 1) then
        c3 = 1988.787252623633_dp + jj * (30.640977403131618_dp      - &
        1.3788285034575731_dp  * nl) + jj2 * (-3.7569974635218113_dp + &
        0.18250747692508693_dp * nl) - 444.90038297813294_dp * nl    + &
        27.73818853998447_dp * nl2 - 0.4544693728045063_dp * nl3     + &
        (-20.703548812349908_dp + jj * (7.513994927043625_dp         - &
        0.36501495385017385_dp * nl) + 1.0004503142481713_dp * nl) * ss + &
        (-3.7569974635218113_dp + 0.18250747692508693_dp * nl) * ss2
      else if (l == 2) then
        c3 = 1901.5017267203552_dp + jj * (4.079182167775893_dp      - &
        0.11282798880462161_dp * nl) + jj2 * (-0.1299105555452488_dp + &
        0.0040335703497889385_dp * nl) - 402.34691146347046_dp * nl  + &
        24.612549396122763_dp * nl2 - 0.4008723708453693_dp * nl3    + &
        (-3.0561439858025263_dp + jj * (0.25982111109049766_dp       - &
        0.008067140699577877_dp * nl) + 0.10222110919114273_dp * nl) * ss + &
        (-0.1299105555452488_dp + 0.0040335703497889385_dp * nl) * ss2
      end if

    else if (n == 4) then
      if (l == 0) then
        c3 = 1334.1657153067044_dp - 407.28551220371247_dp * nl  + &
        27.83594776731537_dp * nl2 - 0.4347159974378718_dp * nl3 + &
        (137.000933835361_dp - 7.940117353183779_dp * nl) * ss
      else if (l == 1) then
        c3 = 1995.6832394075916_dp + jj * (27.958016166615735_dp       - &
        1.3404873242204196_dp * nl) + jj2 * (-3.4373507564943924_dp    + &
        0.17517635102222015_dp * nl) - 457.0465561083301_dp * nl       + &
        28.689641476076723_dp * nl2 - 0.4683744068804726_dp * nl3      + &
        (-18.85267096658803_dp + jj * (6.874701512988785_dp            - &
        0.3503527020444403_dp * nl) + 0.9545817001042881_dp * nl) * ss + &
        (-3.4373507564943924_dp + 0.17517635102222015_dp * nl) * ss2
      else if (l == 2) then
        c3 = 1997.2069835963655_dp + jj * (3.958438263203731_dp          - &
        0.14118158483048285_dp * nl) + jj2 * (-0.12522904650665853_dp    + &
        0.004792696300685975_dp * nl) - 431.6105518485525_dp * nl        + &
        26.621747666930165_dp * nl2 - 0.4360214587004526_dp * nl3        + &
        (-2.968790587941932_dp + jj * (0.25045809301331706_dp            - &
        0.00958539260137195_dp * nl) + 0.11908627681361951_dp * nl) * ss + &
        (-0.12522904650665853_dp + 0.004792696300685975_dp * nl) * ss2
      else if (l == 3) then
        c3 = 1882.185489900684_dp + jj * (1.1136805133484151_dp             - &
        0.02558450978851483_dp * nl) + jj2 * (-0.017109601029735064_dp      + &
        0.00044023804726528275_dp * nl) - 396.62888749856256_dp * nl        + &
        24.241424872934203_dp * nl2 - 0.3942467526859063_dp * nl3           + &
        (-0.844278574298179_dp + jj * (0.03421920205947013_dp               - &
        0.0008804760945305655_dp * nl) + 0.024342937599636186_dp * nl) * ss + &
        (-0.017109601029735064_dp + 0.00044023804726528275_dp * nl) * ss2
      end if

    end if

  end function c3

!ccccccccccccccc

   real (dp) function Binomial(self, i, j)
     class (NRQCD), intent(in) :: self
     integer      , intent(in) :: i, j
     integer                   :: k

    select case(j)
    case(0)  ;  Binomial = 1
    case(:-1);  Binomial = 0
    case default

      Binomial = i

      do k = 1, j - 1
        Binomial = Binomial * (i - k)
      end do

      Binomial = Binomial/self%Listfact(j)

    end select

   end function Binomial

!ccccccccccccccc

  pure function factList(n) result(list)
    integer, intent(in)     :: n
    integer, dimension(0:n) :: list
    integer                 :: i

    list(0) = 1

    do i = 1, n
      list(i) = list(i - 1) * i
    end do

  end function factList

!ccccccccccccccc

  real (dp) function DeltaCharmBin(self, alpha, mb)
    class (NRQCD), intent(inout)   :: self
    real (dp)    , intent(in)      :: mb, alpha
    real (dp)                      :: r, lg
    integer                        :: l

    if ( self%up(:4) == 'down' ) then
      DeltaCharmBin = self%DeltaCharm(alpha, mb); return
    end if

    r = 3 * self%n * self%mc/2/mb/alpha; DeltaCharmBin = 0

    if (100 * r > 1 ) then
      DeltaCharmBin = self%DeltaCharm(alpha, mb) - self%ZeroBin(alpha, mb)
      return
    else if ( r <= tiny(1._dp) ) then
      DeltaCharmBin = 0; return
    end if

    lg = log(r/2)

    if ( self%average(:3) == 'yes' ) then

      do l = 0, self%n - 1
        self%l = l
        DeltaCharmBin = DeltaCharmBin + (2 * l + 1) * res(r)
      end do

      DeltaCharmBin = DeltaCharmBin/self%n**2

    else

      DeltaCharmBin = res(r)

    end if

  contains

    real (dp) function res(r)
      real (dp), intent(in) :: r

      if (self%n == 1) then

        if (self%l == 0) then
          res = - 3 * Pi * r/2 + 9 * r**2/2 - 2 * Pi * r**3 + 3 * r**4/16 &
          - 5 * r**6/6 - 369 * r**8/512 - 401 * r**10/640 + lg * ( - 15 * r**4/4&
          - 7 * r**6/4 - 81 * r**8/64 - 33 * r**10/32)
        end if

      else if (self%n == 2) then

        if (self%l == 0) then

          res = - 3 * Pi * r + 18 * r**2 - 14 * Pi * r**3 + 237 * r**4/16 &
          - 713 * r**6/24 - 19161 * r**8/512 - 14693 * r**10/310 + lg * ( - 165 * r**4/4&
          - 77 * r**6/2 - 2997 * r**8/64 - 231 * r**10/4)

        else if (self%l == 1) then
          res = - 3 * Pi * r + 15 * r**2 - 10 * Pi * r**3 + 109 * r**4/16 &
          - 359 * r**6/24 - 9033 * r**8/512 - 1647 * r**10/80 + lg * ( - 105 * r**4/4&
          - 21 * r**6 - 1485 * r**8/64 - 429 * r**10/16)
        end if

      else if (self%n == 3) then

        if (self%l == 0) then
          res = - 9 * Pi * r/2 + 81 * r**2/2 - 46 * Pi * r**3 + 2097 * r**4/16 &
          - 7937 * r**6/24 - 323289 * r**8/512 - 707151 * r**10/640 + lg * ( - 765 * r**4/4&
          - 1281 * r**6/4 - 39933 * r**8/64 - 36333 * r**10/32)
        else if (self%l == 1) then
          res = - 9 * Pi * r/2 + 75 * r**2/2 - 40 * Pi * r**3 + 1599 * r**4/16 &
          - 11533 * r**6/48 - 222117 * r**8/512 - 232803 * r**10/320 + lg * ( - 315 * r**4/4&
          - 483 * r**6/2 - 28215 * r**8/64 - 24453 * r**10/32)
        else if (self%l == 2) then
          res = - 9 * Pi * r/2 + 63 * r**2/2 - 28 * Pi * r**3 + 3747 * r**4/80 &
          - 25037 * r**6/240 - 426537 * r**8/2560 - 16323 * r**10/64 + lg * ( - 189 * r**4/2&
          - 231 * r**6/2 - 11583 * r**8/64 - 9009 * r**10/32)
        end if

      else if (self%n == 4) then
        if (self%l == 0) then
          res = - 6 * Pi * r + 72 * r**2 - 108 * Pi * r**3 + 8817 * r**4/16 &
          - 15885 * r**6/8 - 2867661 * r**8/512 - 4399953 * r**10/320 + lg * ( - 585 * r**4&
          - 1582 * r**6 - 297999 * r**8/64 - 96657 * r**10/8)
        else if (self%l == 1) then
          res = - 6 * Pi * r + 69 * r**2 - 100 * Pi * r**3 + 7661 * r**4/16 &
          - 99209 * r**6/60 - 2290317 * r**8/512 - 679527 * r**10/64 + lg * ( - 525 * r**4&
          - 1344 * r**6 - 242055 * r**8/64 - 151437 * r**10/16)
        else if (self%l == 2) then
          res = - 6 * Pi * r + 63 * r**2 - 84 * Pi * r**3 + 27577 * r**4/80 &
          - 65089 * r**6/60 - 6857721 * r**8/2560 - 378641 * r**10/64 + lg * ( - 819 * r**4/2&
          - 924 * r**6 - 150579 * r**8/64 - 87087 * r**10/16)
        else if (self%l == 3) then
          res = - 6 * Pi * r + 54 * r**2 - 60 * Pi * r**3 + 19031 * r**4/112 &
          - 383231 * r**6/840 - 3433635 * r**8/3584 - 839075 * r**10/448 + lg * ( - 495 * r**4/2&
          - 429 * r**6 - 57915 * r**8/64 - 7293 * r**10/4)
        end if

      end if

      res = - res/3

    end function res

  end function DeltaCharmBin

!ccccccccccccccc

  real (dp) function DeltaCharm(self, alpha, mb)
    class (NRQCD), intent(inout)   :: self
    real (dp)    , intent(in)      :: mb, alpha
    real (dp)                      :: r, root, ArTan
    real (dp), dimension(4*self%n) :: rho
    integer                        :: l

    DeltaCharm = 0
    r = 3 * self%n * self%mc/2/mb/alpha; rho = powList(r,4*self%n)

    if ( r > 1 ) then
      ArTan = Atan(  sqrt( (r - 1)/(r + 1) )  ); root = sqrt( rho(2) - 1 )
    else
      root = sqrt( 1 - rho(2) ); ArTan = log( (1 + root)/r  )/2
    end if

    if ( self%average(:3) == 'yes' ) then

      do l = 0, self%n - 1
        self%l = l
        DeltaCharm = DeltaCharm + (2 * l + 1) * res(r)
      end do

      DeltaCharm = DeltaCharm/self%n**2

    else

      DeltaCharm = res(r)

    end if

  contains

    real (dp) function res(r)
      real (dp), intent(in) :: r

      if (self%n == 1) then

        if (self%l == 0) then
          res = 11._dp/3 - 3 * pi * r/2 + 4 * rho(2) - 2 * pi * rho(3) &
          - 2 * ( 2 - rho(2) - 4 * rho(4) )/root * ArTan
        end if

      else if (self%n == 2) then

        if (self%l == 0) then

          res = - pi * r * ( 3 + 14 * rho(2) ) + ArTan * &
          ( 10 * rho(2) - 4 + 75 * rho(4) - 128 * rho(6) + 56 * rho(8) )/root**5 + &
          ( 28 + 49 * rho(2) - 272 * rho(4) + 168 * rho(6) )/6/root**4

        else if (self%l == 1) then
          res = - pi * r * ( 3 + 10 * rho(2) ) + ArTan * &
          ( 10 * rho(2) - 4 + 45 * rho(4) - 88 * rho(6) + 40 * rho(8) )/root**5 + &
          ( 32 + 23 * rho(2) - 184 * rho(4) + 120 * rho(6) )/6/root**4
        end if

      else if (self%n == 3) then

        if (self%l == 0) then
          res = - pi * r * ( 9 + 92 * rho(2) )/2 + ArTan * &
          ( 36 * rho(2) + 702 * rho(4) - 2109 * rho(6) + 2736 * rho(8) &
          - 1620 * rho(10) + 368 * rho(12) - 8 )/2/root**9 + ( 64 + &
          224 * rho(2) - 3111 * rho(4) + 5528 * rho(6) - 4124 * rho(8) + &
          1104 * rho(10) )/12/root**8
        else if (self%l == 1) then
          res = - pi * r * ( 9 + 80 * rho(2) )/2 + ArTan * &
          ( 72 * rho(2) + 1134 * rho(4) - 3633 * rho(6) + 4716 * rho(8) &
          - 2808 * rho(10) + 640 * rho(12) - 16 )/4/root**9 + ( 140 + &
          328 * rho(2) - 5115 * rho(4) + 9556 * rho(6) - 7144 * rho(8) + &
          1920 * rho(10) )/24/root**8
        else if (self%l == 2) then
          res = - pi * r * ( 9 + 56 * rho(2) )/2 + ArTan * &
          ( 72 * rho(2) + 630 * rho(4) - 2373 * rho(6) + 3204 * rho(8) &
          - 1944 * rho(10) + 448 * rho(12) - 16 )/4/root**9 + ( 748 + &
          728 * rho(2) - 16035 * rho(4) + 32204 * rho(6) - 24680 * rho(8) + &
          6720 * rho(10) )/120/root**8
        end if

      else if (self%n == 4) then
        if (self%l == 0) then
          res = - 6 * pi * r * ( 1 + 18 * rho(2) ) + ArTan * ( 208 * rho(2) &
          + 8788 * rho(4) - 34670 * rho(6) + 76531 * rho(8) - 89072 * rho(10) + &
          60528 * rho(12) - 22272 * rho(14) + 3456 * rho(16) - 32 )/8/root**13 + &
          ( 280 + 1752 * rho(2) - 42852 * rho(4) + 116345 * rho(6) - 178800 * rho(8)&
          + 142416 * rho(10) - 59904 * rho(12) + 10368 * rho(14) )/48/root**12
        else if (self%l == 1) then
          res = - 2 * pi * r * ( 3 + 50 * rho(2) ) + ArTan * ( 208 * rho(2) &
          + 7828 * rho(4) - 32238 * rho(6) + 70137 * rho(8) - 82368 * rho(10) + &
          55952 * rho(12) - 20608 * rho(14) + 3200 * rho(16) - 32 )/8/root**13 + &
          ( 1496 + 7464 * rho(2) - 191160 * rho(4) + 541079 * rho(6) - 819744 * rho(8)&
          + 658800 * rho(10) - 277120 * rho(12) + 48000 * rho(14) )/240/root**12
        else if (self%l == 2) then
          res = - 6 * pi * r * ( 1 + 14 * rho(2) ) + ArTan * ( 208 * rho(2) &
          + 5980 * rho(4) - 26946 * rho(6) + 57915 * rho(8) - 68640 * rho(10) + &
          46800 * rho(12) - 17280 * rho(14) + 2688 * rho(16) - 32 )/8/root**13 + &
          ( 1576 + 5544 * rho(2) - 149136 * rho(4) + 454325 * rho(6) - 681408 * rho(8)&
          + 550704 * rho(10) - 232320 * rho(12) + 40320 * rho(14) )/240/root**12
        else if (self%l == 3) then
          res = - 6 * pi * r * ( 1 + 10 * rho(2) ) + ArTan * ( 208 * rho(2) &
          + 3388 * rho(4) - 18018 * rho(6) + 39897 * rho(8) - 48048 * rho(10) + &
          33072 * rho(12) - 12288 * rho(14) + 1920 * rho(16) - 32 )/8/root**13 + &
          ( 11512 + 20808 * rho(2) - 652380 * rho(4) + 2169953 * rho(6) - 3325968 * rho(8)&
          + 2719920 * rho(10) - 1155840 * rho(12) + 201600 * rho(14) )/1680/root**12
        end if

      end if

      res = - res/3

    end function res

  end function DeltaCharm

!ccccccccccccccc

  real (dp) function DeltaCharmDerBin(self, alpha, mb)
    class (NRQCD), intent(inout) :: self
    real (dp)    , intent(in)    :: mb, alpha
    real (dp)                    :: lg, r
    integer                      :: l

    DeltaCharmDerBin = 0

    if ( self%up(:4) == 'down' ) then
      DeltaCharmDerBin = self%DeltaCharmDer(alpha, mb); return
    end if

    r = 3 * self%n * self%mc/2/mb/alpha

    if (100 * r > 1 ) then
      DeltaCharmDerBin = self%DeltaCharmDer(alpha, mb) + 2._dp/3
      return
    else if ( r <= tiny(1._dp) ) then
      DeltaCharmDerBin = 0; return
    end if

    lg = log(r/2)

    if ( self%average(:3) == 'yes' ) then

      do l = 0, self%n - 1
        self%l = l
        DeltaCharmDerBin = DeltaCharmDerBin + (2 * l + 1) * res(r)
      end do

      DeltaCharmDerBin = DeltaCharmDerBin/self%n**2

    else

      DeltaCharmDerBin = res(r)

    end if

  contains

    real (dp) function res(r)
      real (dp), intent(in) :: r

      if (self%n == 1) then

        if (self%l == 0) then
          res = 2 - 3 * pi * r/2 + 9 * r**2 - 6 * pi * r**3 + &
          - (3 + 15 * lg) * r**4 - (27._dp/4 + 21 * lg/2) * r**6 - &
          (225._dp/32 + 81 * lg/8) * r**8 - (467._dp/64 + 165 * lg/16) * r**10
        end if

      else if (self%n == 2) then

        if (self%l == 0) then

          res = 2 - 3 * pi * r + 36 * r**2 - 42 * pi * r**3  - &
          (201._dp/2 + 165 * lg) * r**4 - (867._dp/4 + 231 * lg) * r**6 - &
          (11079._dp/32 + 2997 * lg/8) * r**8 - (16541._dp/32 + 1155 * lg/2) * r**10

        else if (self%l == 1) then
          res = 2 - 3 * pi * r + 30 * r**2 - 30 * pi * r**3  - &
          (107._dp/2 + 105 * lg) * r**4 - (443._dp/4 + 126 * lg) * r**6 - &
          (5259._dp/32 + 1485 * lg/8) * r**8 - (3723._dp/16 + 2145 * lg/8) * r**10
        end if

      else if (self%n == 3) then

        if (self%l == 0) then
          res = 2 - 9 * pi * r/2 + 81 * r**2 - 138 * pi * r**3 - &
          (1431._dp/2 + 765 * lg) * r**4 - (4609._dp/2 + 3843 * lg/2) * r**6  - &
          (181611._dp/32 + 39933 * lg/8) * r**8 - (779817._dp/64 + &
          181665 * lg/16) * r**10
        else if (self%l == 1) then
          res = 2 - 9 * pi * r/2 + 75 * r**2 - 120 * pi * r**3 - &
          (2229._dp/4 + 630 * lg) * r**4 - (13465._dp/8 + 1449 * lg) * r**6 - &
          (62583._dp/16 + 28215 * lg/8) * r**8 - (32157._dp/4 - 122265 * lg/16) * r**10
        else if (self%l == 2) then
          res = 2 - 9 * pi * r/2 + 63 * r**2 - 84 * pi * r**3 - &
          (5637._dp/20 + 378 * lg) * r**4 - (29657._dp/40 + 693 * lg) * r**6 - &
          (121113._dp/80 - 11583 * lg/8) * r**8 - (2832 + 45045 * lg/16) * r**10
        end if

      else if (self%n == 4) then
        if (self%l == 0) then
          res = 2 - 6 * pi * r + 144 * r**2 - 324 * pi * r**3 - &
          (11157._dp/4 + 2340 * lg) * r**4 - (53983._dp/4 + 9492 * lg) * r**6 - &
          (791415._dp/16 + 297999 * lg/8) * r**8 - (4786581._dp/32 + &
          483285 * lg/4) * r**10
        else if (self%l == 1) then
          res = 2 - 6 * pi * r + 138 * r**2 - 300 * pi * r**3 - &
          (9761._dp/4 + 2100 * lg) * r**4 - (112649._dp/10 + 8064 * lg) * r**6 - &
          (633093._dp/16 + 242055 * lg/8) * r**8 - (3700509._dp/32 + &
          757185 * lg/8) * r**10
        else if (self%l == 2) then
          res = 2 - 6 * pi * r + 126 * r**2 - 252 * pi * r**3 - &
          (35767._dp/20 + 1638 * lg) * r**4 - (74329._dp/10 + 5544 * lg) * r**6 &
          - (951327._dp/40 + 150579 * lg/8) * r**8 - (2067379._dp/32 + &
          435435 * lg/8) * r**10
        else if (self%l == 3) then
          res = 2 - 6 * pi * r + 108 * r**2 - 180 * pi * r**3 - &
          (25961._dp/28 + 990 * lg) * r**4 - (443291._dp/140 + 2574 * lg) * r**6 &
          - (59985._dp/7 + 57915 * lg/8) * r**8 - (4603783._dp/224 + &
          36465 * lg/2) * r**10
        end if

      end if

      res = - res/3

    end function res

  end function DeltaCharmDerBin

!ccccccccccccccc

  real (dp) function DeltaCharmDer(self, alpha, mb)
    class (NRQCD), intent(inout)   :: self
    real (dp)    , intent(in)      :: mb, alpha
    real (dp)                      :: r, root, ArTan
    real (dp), dimension(4*self%n) :: rho
    integer                        :: l

    DeltaCharmDer = 0

    r = 3 * self%n * self%mc/2/mb/alpha; rho = powList(r,4*self%n)

    if ( r > 1 ) then
      ArTan = Atan(  sqrt( (r - 1)/(r + 1) )  ); root = sqrt( rho(2) - 1 )
    else
      root = sqrt( 1 - rho(2) ); ArTan = - log( (1 + root)/r  )/2
    end if

    if ( self%average(:3) == 'yes' ) then

      do l = 0, self%n - 1
        self%l = l
        DeltaCharmDer = DeltaCharmDer + (2 * l + 1) * res(r)
      end do

      DeltaCharmDer = DeltaCharmDer/self%n**2

    else

      DeltaCharmDer = res(r)

    end if

  contains

    real (dp) function res(r)
      real (dp), intent(in) :: r

      if (self%n == 1) then

        if (self%l == 0) then
          res = (  4 + 14 * rho(2) - 24 * rho(4) + 3 * pi * r * &
          ( 4 * rho(4) - 1 - 3 * rho(2) )  )/( 2 - 2 * rho(2) ) + ( 6 * ArTan *  &
          rho(4) * ( 4 * rho(2) - 5 ) )/root**3
        end if

      else if (self%n == 2) then

        if (self%l == 0) then
          res = - 3 * pi * ( r + 14 * rho(3) ) + (  3 * ArTan * rho(4) &
          * (  231 * rho(2) - 110 - 192 * rho(4) + 56 * rho(6) )  )/root**7 + &
          ( 405 * rho(4) - 464 * rho(6) + 168 * rho(8) - 4 - 60 * rho(2) )/2/root**6
        else if (self%l == 1) then
          res = - 3 * pi * ( r + 10 * rho(3) ) + (  3 * ArTan * rho(4) * &
          ( 161 * rho(2) - 136 * rho(4) + 40 * rho(6) - 70 )  )/root**7 + &
          ( 275 * rho(4) - 328 * rho(6) + 120 * rho(8) - 4 - 48 * rho(2) )/2/root**6
        end if

      else if (self%n == 3) then

        if (self%l == 0) then
          res = - 3 * pi * r * ( 3 + 92 * rho(2) )/2 + (  3 * ArTan * &
          rho(4) *  ( 3048 * rho(2) - 1020 - 5187 * rho(4) + 4488 * rho(6) -&
          2012 * rho(8) + 368 * rho(10) )  )/2/root**11 + ( 4402 * rho(4) - 8 &
          - 284 * rho(2) - 9017 * rho(6) + 10048 * rho(8) - 5300 * rho(10) + &
          1104 * rho(12) )/4/root**10
        else if (self%l == 1) then
          res = - 3 * pi * r * ( 3 + 80 * rho(2) )/2 + (  3 * ArTan * rho(4)&
          * (  5376 * rho(2) - 1680 - 8943 * rho(4) + 7788 * rho(6) - 3496 * rho(8) &
          + 640 * rho(10) )  )/4/root**11 + ( 7298 * rho(4) - 16 - 520 * rho(2) - &
          15925 * rho(6) + 17396 * rho(8) - 9208 * rho(10) + 1920 * rho(12) )/&
          8/root**10
        else if (self%l == 2) then
          res = - 3 * pi * r * ( 3 + 56 * rho(2) )/2 + (  3 * ArTan * &
          rho(4) * ( 3696 * rho(2) - 1008 - 6171 * rho(4) + 5412 * rho(6) -  &
          2440 * rho(8) + 448 * rho(10) )  )/4/root**11 + ( 23074 * rho(4) - 80 &
          - 2120 * rho(2) - 54893 * rho(6) + 60364 * rho(8) - 32120 * rho(10) + &
          6720 * rho(12) )/40/root**10
        end if

      else if (self%n == 4) then
        if (self%l == 0) then
          res = - 6 * pi * ( r + 54 * rho(3) ) + (  3 * ArTan * rho(4) &
          * (  42976 * rho(2) - 12480 - 123186 * rho(4) + 169355 * rho(6) - 153040 * &
          rho(8) + 83760 * rho(10) - 25856 * rho(12) + 3456 * rho(14) )  )/&
          8/root**15 + ( 60084 * rho(4) - 32 - 2080 * rho(2) - 153088 * rho(6) + &
          320641 * rho(8) - 325392 * rho(10) + 205200 * rho(12) - 70656 * rho(14)&
          + 10368 * rho(16) )/16/root**14
        else if (self%l == 1) then
          res = - 6 * pi * ( r + 50 * rho(3) ) + (  3 * ArTan * rho(4) * &
          ( 40992 * rho(2) - 11200 - 111810 * rho(4) + 157665 * rho(6) -         &
          141440 * rho(8) + 77520 * rho(10) - 23936 * rho(12) + 3200 * rho(14) )  )&
          /8/root**15 + ( 269140 * rho(4) - 160 - 9920 * rho(2)- 733588 * rho(6) &
          + 1454511 * rho(8) - 1514848 * rho(10) + 949040 * rho(12) - &
          327040 * rho(14) + 48000 * rho(16) )/80/root**14
        else if (self%l == 2) then
          res = - 6 * pi * ( r + 42 * rho(3) ) + (  3 * ArTan * rho(4) * &
          ( 35952 * rho(2) - 8736 - 91566 * rho(4) + 132275 * rho(6) - &
          118560 * rho(8) + 65040 * rho(10) - 20096 * rho(12) + 2688 * rho(14) )  )&
          /8/root**15 + ( 210268 * rho(4) - 160 - 8960 * rho(2) - 645684 * rho(6) + &
          1198013 * rho(8) - 1270336 * rho(10) + 796144 * rho(12) - &
          274560 * rho(14) + 40320 * rho(16) )/80/root**14
        else if (self%l == 3) then
          res = - 6 * pi * ( r + 30 * rho(3) ) + (  3 * ArTan * rho(4) * &
          ( 25872 * rho(2) - 5280 - 64350 * rho(4) + 93665 * rho(6) - 84240 * rho(8)&
          + 46320 * rho(10) - 14336 * rho(12) + 1920 * rho(14) )  )/8/root**15 + &
          ( 919060 * rho(4) - 3230856 * rho(6) - 1120 - 52640 * rho(2) + &
          5925737 * rho(8) - 6313456 * rho(10) + 3967600 * rho(12)  &
          - 1370880 * rho(14) + 201600 * rho(16) )/560/root**14
        end if

      end if

      res = - res/3

    end function res

  end function DeltaCharmDer

!ccccccccccccccc

  real (dp) function ZeroBin(self, alpha, mb)
    class (NRQCD), intent(in) :: self
    real (dp)    , intent(in) :: mb, alpha
    real (dp)                 :: rho, coef

    rho = 3 * self%n * self%mc/2/mb/alpha

    if ( self%average(:3) /= 'yes' ) then
      coef = 2 * Harmonic(self%l + self%n)
    else
      coef = 1 - 1._dp/self%n + 2 * Harmonic(self%n)
    end if

    ZeroBin = - ( 5._dp/3 + coef + 2 * log(rho/2) )/3

  end function ZeroBin

!ccccccccccccccc

  pure real (dp) function Harmonic(n)
    integer, intent(in) :: n
    integer             :: i

    Harmonic = 0

    do i = 1, n
      Harmonic = Harmonic + 1._dp/i
    end do

  end function Harmonic

!ccccccccccccccc

end module NRQCDClass
