
Module ToppikClass
  use Constants; use QuadPack; use hyperGeo; implicit none; private

  integer  , parameter, public :: nmax = 360
  real (dp), parameter :: wgamma = 2.07_dp, zmass = 91.187_dp, eps = 1e-4_dp, &
  bmass = 4.7_dp, GFERMI = 1.16637e-5_dp, wmass = 80.33_dp

!ccccccccccccccc

  type, public :: Toppik

  integer      :: kinflg, gcflg, vflag, ILFLAG, np, nf
  real (dp)    :: alamb5, alphas, tmass, tgamma, xp, xpmax, dcut, energy, dd, &
  kincom, kincoa, kincov, cplas, scale, c0, c1, c2, cdeltc, cdeltl, cfullc,   &
  cfulll, crm2, QCUT, QMAT1, ALR, a55, b5, c5, alfmt, alphrc, vzero, xim, xdi

  complex (dp), dimension(nmax) :: zvfct
  real (dp)   , dimension(nmax) :: xpp, ww, dsdp

  contains

    procedure, pass(self), public  :: CrossSec, pGrid, Weight, Sigma, Kernel
    procedure, pass(self), private :: g0c, sae, a, g0, vhat

  end type Toppik

!ccccccccccccccc

  interface Toppik
    module procedure ToppikIn
  end interface Toppik

  contains

! *********************************************************************
!
! Working version with all the different original potentials
!  like (p^2+q^2)/|p-q|^2, not transformed in terms of delta and 1/r^2;
! accuracy eps=1.d-4 possible (only), but should be save, 13.8.'98, tt.
!
! *********************************************************************
!
! -- Calculation of the Green function in momentum space by solving the
!     Lippmann-Schwinger equation
!     G(p) = G_0(p) + G_0(p) int_0^cutn V(p,q) G(q) dq
!
! -- Written by Thomas Teubner, Hamburg, November 1998
!     Rewritten in fortran 95 by Vicent Mateu, July 2017, for hardwired npot = 5
!     * Based on TOPPIK Version 1.1
!        from M. Jezabek and TT, Karlsruhe, June 1992
!     * Version originally for non-constant top-width
!     * Constant width supplied here
!     * No generator included
!
! -- Use of double precision everywhere
!
! -- All masses, momenta, energies, widths in GeV
!
! -- Input parameters:
!
!    energy  :  E = Sqrt[s] - 2 * topmass
!    tm      :  top mass (in the Pole scheme)
!    tg      :  top width
!    alphas  :  alpha_s^{MSbar,n_f=5}(scale)
!    scale   :  soft scale  mu_{soft}
!    cutn    :  numerical UV cutoff on all momenta
!                (UV cutoff of the Gauss-Legendre grid)
!    cutv    :  renormalization cutoff on the
!                 delta-, the (p^2+q^2)/(p-q)^2-, and the
!                  1/r^2-[1/|p-q|]-potential:
!                 if (max(p,q).ge.cutv) then the three potentials
!                  are set to zero in the Lippmann-Schwinger equation
!    c0      :  0th order coefficient for the Coulomb potential,
!                 see calling example above
!    c1      :  1st order coefficient for the Coulomb potential
!    c2      :  2nd order coefficient for the Coulomb potential
!    cdeltc  :  constant of the delta(r)-
!                 [= constant in momentum space-] potential
!    cdeltl  :  constant for the additional log(q^2/mu^2)-part of the
!                 delta-potential:
!                  cdeltc*1 + cdeltl*log(q^2/mu^2)
!    cfullc  :  constant of the (p^2+q^2)/(p-q)^2-potential
!    cfulll  :  constant for the additional log(q^2/mu^2)-part of the
!                 (p^2+q^2)/(p-q)^2-potential
!    crm2    :  constant of the 1/r^2-[1/|p-q|]-potential
!    kincm   :  } kinetic corrections in the 0th order Green-function:
!    kinca   :  }  G_0(p):=1/[E+iGamma_t-p^2/m_t]*(1+kincm)+kinca
!     !!! WATCH THE SIGN IN G_0 !!!
!    jknflg   :  flag for these kinetic corrections:
!                 0 : no kinetic corrections applied
!                 1 : kinetic corrections applied with cutoff cutv
!                      for  kinca  only
!                 2 : kinetic corrections applied with cutoff cutv
!                      for  kinca  AND  kincm
!    jgcflg   :  flag for G_0(p) in the LS equation:
!                 0 (standard choice) : G_0(p) as given above
!                 1 (for TIPT)        : G_0(p) = G_c^{0}(p) the 0th
!                                        order Coulomb-Green-function
!                                        in analytical form; not for
!                                        momenta  p > 1000*topmass
!    xkincv   :  additional kinematic vertex correction in G_0, see below:
!    jvflg    :  flag for the additional vertex correction  xkincv  in the
!                 ``zeroth order'' G_0(p) in the LS-equation:
!                 0 : no correction, means  G = G_0 + G_0 int V G
!                      with G_0=1/[E+iGamma_t-p^2/m_t]*(1+kincm)+kinca
!                 1 : apply the correction in the LS equation as
!                      G = G_0 + xkincv*p^2/m_t^2/[E+iGamma_t-p^2/m_t] +
!                          G_0 int V G
!                     and correct the integral over Im[G(p)] to get sigma_tot
!                     from the optical theorem by the same factor.
!                     The cutoff  cutv  is applied for these corrections.
!
! -- Output:
!
!    xim      :  R_{ttbar} from the imaginary part of the green
!                 function
!    xdi      :  R_{ttbar} form the integral over the momentum
!                 distribution (no cutoff but the numerical one here!!)
!    np       :  number of points used for the grid; fixed in tttoppik
!    xpp      :  1-dim array (max. 400 elements) giving the momenta of
!                 the Gauss-Legendre grid (pp(i) in the code)
!    ww      :  1-dim array (max. 400 elements) giving the corresponding
!                 Gauss-Legendre weights for the grid
!    dsdp    :  1-dim array (max. 400 elements) giving the
!                 momentum distribution of top: d\sigma/dp,
!                  normalized to R,
!                  at the momenta of the Gauss-Legendre grid xpp(i)
!    zvfct    :  1-dim array (max. 400 elements) of COMPLEX*16 numbers
!                 giving the vertex function K(p), G(p)=K(p)*G_0(p)
!                 at the momenta of the grid
!
! *********************************************************************

  type (Toppik) function ToppikIn(nf, energy, tm, tg, alphas, scale, cutn, &
  cutv, c0, c1, c2, cdeltc, cdeltl, cfullc, cfulll, crm2, kincm, kinca, jknflg,&
  jgcflg, kincv, jvflg)

    integer  , intent(in) :: jknflg, jgcflg, jvflg, nf
    real (dp), intent(in) :: energy, tm, tg, alphas, scale, cutn, cutv, &
    c0, c1, c2, cdeltc, cdeltl, cfullc, cfulll, crm2, kincm, kinca, kincv

    integer   :: n, i
    real (dp) :: etot, consde, const, sig1, sig2, critp

    complex (dp), dimension(nmax)   :: bb, gg, a1
    real (dp)   , dimension(nmax)   :: pp, w1
    real (dp)   , dimension(nmax/3) :: xx, w2

    ToppikIn%nf = nf

! Number of points to evaluate on the integral equation (<=400 and n mod 3 = 0 !!):

    n = nmax;  ToppikIn%np = nmax

! Input:

    ToppikIn%energy = energy; ToppikIn%tmass = tm; ToppikIn%tgamma = tg
    ToppikIn%cplas = alphas; ToppikIn%c0 = c0; ToppikIn%c1 = c1
    ToppikIn%c2 = c2; ToppikIn%cdeltc = cdeltc; ToppikIn%cdeltl = cdeltl
    ToppikIn%cfullc = cfullc; ToppikIn%cfulll = cfulll; ToppikIn%crm2 = crm2
    ToppikIn%kincom = kincm; ToppikIn%kincoa = kinca; ToppikIn%kincov = kincv
    ToppikIn%kinflg = jknflg; ToppikIn%gcflg = jgcflg; ToppikIn%vflag = jvflg
    ToppikIn%alphas = alphas; ToppikIn%scale = scale

! Cut for divergent potential-terms for large momenta in the function vhat and in the integrals a(p):

    ToppikIn%dcut = cutv

! Numerical Cutoff of all momenta (maximal momenta of the grid):

    ToppikIn%xpmax = cutn

    if (cutv > cutn) then
      print*, ' dcut > xpmax  makes no sense! Stop.';  stop
    endif

! Not needed for the fixed order potentials:

    ToppikIn%alamb5 = 0.2_dp;  etot = 2 * tm + energy

! For pure coulomb and fixed order potentials there is no delta-part:

    consde = 0

! Delta-part of potential is absorbed by subtracting vzero from the
! original energy (shift from the potential to the free Hamiltonian):

    ToppikIn%vzero = consde / (2 * pi)**3

! Or different (simpler) method, good for V_JKT:

    if (energy <= 0) then
      critp = tm/3
    else
      critp = max( tm/3, 2 * sqrt(energy * tm) )
    endif

    call gauleg( 0._dp, critp, pp(:2*n/3), w1(:2*n/3), 2*n/3 )
    call gauleg(1/ToppikIn%xpmax, 1/critp, pp(2*n/3 + 1:n), w1(2*n/3 + 1:n), n/3)

! Do substitution p => 1/p for the last interval explicitly:

    pp(2 * n/3 + 1:) = 1/pp(2 * n/3 + 1:)

! Reorder the arrays for the third interval:

    xx = pp(2 * n/3 + 1:)     ; w2 = w1(2 * n/3 + 1:)
    pp(n:2 * n/3 + 1:-1) = xx ; w1(n:2 * n/3 + 1:-1) = w2

! Calculate the integrals a(p) for the given momenta pp(i)
!  and store weights and momenta for the output arrays:

    do  i = 1, n
      a1(i) = ToppikIn%a( pp(i) )
    end do

    ToppikIn%xpp(:n) = pp   ;  ToppikIn%ww(:n) = w1
    ToppikIn%xpp(n + 1:) = 0; ToppikIn%ww(n + 1:) = 0

! Solve the integral-equation by solving a system of algebraic equations:

    call ToppikIn%sae(pp, w1, bb, a1, n)

! (The substitution for the integration to infinity  pp => 1/pp is done already.)

    do  i = 1, n
      ToppikIn%zvfct(i) = bb(i);  gg(i) = bb(i) * ToppikIn%g0c( pp(i) )
    end do

! Normalisation on R:

    const = 8 * pi/tm**2

! Proove of the optical theorem for the output values of sae:
!  Simply check if sig1 = sig2.

    sig1 = 0; sig2 = 0

    do  i = 1, 2 * n/3

      if  ( pp(i) < ToppikIn%dcut .and. ToppikIn%vflag == 1 ) then
        sig1 = sig1 + w1(i) * pp(i)**2 * imagpart(   gg(i) * &
        (  1 + ToppikIn%kincov * ToppikIn%g0( pp(i) ) * (pp(i)/tm)**2/ToppikIn%g0c( pp(i) )  )   )
      else
        sig1 = sig1 + w1(i) * pp(i)**2 * imagpart( gg(i) )
      endif

      if  ( pp(i) < ToppikIn%dcut .and. ToppikIn%kinflg /= 0 ) then

        sig2 = sig2 + w1(i) * pp(i)**2 * abs( gg(i) )**2 * tg &
         * ( 1 - pp(i)**2/2/tm**2 )

        ToppikIn%dsdp(i) = pp(i)**2 * abs( gg(i) )**2 * tg &
        * ( 1 - pp(i)**2/2/tm**2 )/(2 * Pi2) * const

      else

        sig2 = sig2 + w1(i) * pp(i)**2 * abs( gg(i) )**2 * tg

        ToppikIn%dsdp(i) = pp(i)**2 * abs( gg(i) )**2 * tg/(2 * Pi2) * const

      endif

    end do

! '*p**2' because of substitution p => 1/p in the integration of p**2*G(p) to infinity

    do  i = 2 * n/3 + 1, n

      if  ( pp(i) < ToppikIn%dcut .and. ToppikIn%vflag == 1 ) then
        sig1 = sig1 + w1(i) * pp(i)**4 * imagpart(   gg(i) * (  1 + ToppikIn%kincov * &
        ToppikIn%g0( pp(i) ) * ( pp(i)/tm )**2/ToppikIn%g0c( pp(i) )  )   )
      else
         sig1 = sig1 + w1(i) * pp(i)**4 * imagpart( gg(i) )
      endif

      if  ( pp(i) < ToppikIn%dcut .and. ToppikIn%kinflg /= 0 ) then

        sig2 = sig2 + w1(i)*pp(i)**4 * abs( gg(i) )**2 * tg &
         * (1 - pp(i)**2/2/tm**2)

         ToppikIn%dsdp(i) = pp(i)**2 * abs( gg(i) )**2 * tg &
         * ( 1 - pp(i)**2/2/tm**2 )/(2 * Pi2) * const

      else

        sig2 = sig2 + w1(i) * pp(i)**4 * abs( gg(i) )**2 * tg

        ToppikIn%dsdp(i) = pp(i)**2 * abs( gg(i) )**2 * tg/(2 * Pi2) * const

      endif

    end do

    ToppikIn%dsdp(n + 1:) = 0; ToppikIn%zvfct(n + 1:) = 0

! Normalisation on R:

    sig1 = sig1 / (2 * Pi2) * const;  sig2 = sig2 / (2 * Pi2) * const

! The results from the momentum space approach finally are:

    ToppikIn%xim = - sig1;  ToppikIn%xdi = sig2

    DO i = 1, n
      IF(  ABS( a1(i) ) > 100  ) THEN
        ToppikIn%xim = 1e4_dp;  ToppikIn%xdi = 1e4_dp
      ENDIF
    ENDDO

  end

! cccccccccccccccccccc

  complex (dp) function a(self, p)
    class (Toppik), intent(inout) :: self
    real (dp), intent(in) :: p
    integer               :: ier, neval
    real(dp)              :: xb1, xb2, ddcut, a1, a2, a3, a4, a5, a6, buf, abserr

    self%xp = p; buf = 0;  a1 = 0; a2 = 0; a3 = 0; a4 = 0; a5 = 0; a6 = 0

    if (self%gcflg == 0) then
      ddcut = self%xpmax
    else if (self%gcflg == 1) then
      ddcut = self%dcut
    else
      write(*,*) ' gcflg wrong! Stop.';  stop
    endif

    if (2 * self%xp < ddcut) then

      xb1 = self%xp; xb2 = 2 * self%xp

! More stable for logarithmically divergent fixed order potentials:

      call qags(fretil1,buf,xb1, eps, eps, a1, abserr, neval, ier)
      call qags(fimtil1,buf,xb1, eps, eps, a2, abserr, neval, ier)
      call qags(fretil1,xb1,xb2, eps, eps, a3, abserr, neval, ier)
      call qags(fimtil1, xb1,xb2, eps, eps, a4, abserr, neval, ier)
      call qags(fretil2, 1/ddcut, 1/xb2, eps, eps, a5, abserr, neval, ier)
      call qags(fimtil2, 1/ddcut, 1/xb2, eps, eps, a6, abserr, neval, ier)

    else if (self%xp < ddcut) then

      xb1 = self%xp;  xb2 = ddcut

      call qags(fretil1, buf, xb1, eps, eps, a1, abserr, neval, ier)
      call qags(fimtil1, buf, xb1, eps, eps, a2, abserr, neval, ier)
      call qags(fretil1, xb1,xb2, eps, eps, a3, abserr, neval, ier)
      call qags(fimtil1, xb1,xb2, eps, eps, a4, abserr, neval, ier)

    else if (ddcut <= self%xp) then

    else
      write(*,*) ' Constellation not possible! Stop.';  stop
    endif

    IF(a1 + a2 + a3 + a4 + a5 + a6 >= 1e4_dp) THEN
      a = 1e4_dp
    ELSE
      a  = 1/(4 * pi2) * cmplx(a1 + a3 + a5, a2 + a4 + a6, kind = dp)
    ENDIF

  contains

! cccccccccccccccccccc

    real (dp) function fretil1(xk)
      real (dp), intent(in) :: xk
      fretil1 = freal(xk)
    end

! cccccccccccccccccccc

    real (dp) function fretil2(xk)
      real (dp), intent(in) :: xk
      fretil2 = freal(1/xk)/xk**2
    end

! cccccccccccccccccccc

    real (dp) function fimtil1(xk)
      real (dp), intent(in) :: xk
      fimtil1 = fim(xk)
    end

! cccccccccccccccccccc

    real (dp) function fimtil2(xk)
      real (dp), intent(in) :: xk
      fimtil2 = fim(1/xk)/xk**2
    end

! cccccccccccccccccccc

    real (dp) function freal(xk)
      real (dp), intent(in) :: xk
       freal = realpart( self%g0c(xk) * self%vhat(self%xp, xk) )
    end

! cccccccccccccccccccc

    real (dp) function fim(xk)
      real (dp), intent(in) :: xk
       fim = imagpart( self%g0c(xk) * self%vhat(self%xp, xk) )
    end

end function a

! cccccccccccccccccccc

  complex (dp) function vhat(self, p, xk)
    class (Toppik), intent(in) :: self
    real (dp)     , intent(in) :: p, xk

    real (dp) :: cnspot, phiint, pm, xkm, a1, a2, b0, b1, xkpln1st, xkpln2nd, xkpln3rd

    complex (dp), parameter :: zi = (0._dp,1._dp)
    real (dp)   , parameter :: cf = 4._dp/3, ca = 3._dp, tf = 0.5_dp

    b0 = 11 - 2._dp/3 * self%nf;  b1 = 102 - 38._dp/3 * self%nf

    a1 = 31 * ca/9 - 20 * tf * self%nf/9
    a2 = (4343._dp/162 + 4 * pi2 - pi2**2/4 + 22 * zeta3/3) * ca**2 - &
    (1798._dp/81 + 56 * zeta3/3) * ca * tf * self%nf - (55._dp/3 - 16 * zeta3) * &
    cf * tf * self%nf + (20 * tf * self%nf/9)**2

    pm = p; xkm = xk; cnspot = - 16 * pi/3

    if (p/xk <= 1.e-5_dp .and. p <= 1.e-5_dp) then
      xkpln1st = 2;  xkpln2nd = - 4 * log(self%scale/xk)
      xkpln3rd = - 6 * log(self%scale/xk)**2
    else if (xk/p <= 1.e-5_dp .and. xk <= 1.e-5_dp) then
      xkpln1st = 2 * (xk/p)**2; xkpln2nd = - 4 * (xk/p)**2 * log(self%scale/p)
      xkpln3rd = - 6 * (xk/p)**2 * log(self%scale/p)**2
    else
     xkpln1st = xk/p * (  log(p + xk) - log( abs(p - xk) )  )
     xkpln2nd = - xk/p * (  log( self%scale/(p + xk) )**2 - log( self%scale/abs(p - xk) )**2  )
     xkpln3rd = - 4 * xk/p * (  log( self%scale/(p + xk) )**3 - log( self%scale/abs(p - xk) )**3  )/3
    end if

    phiint = cnspot * self%alphas

    if ( max(xk,p) < self%dcut ) then
!Coulomb + first + second order corrections:
       vhat = phiint * (self%c0 * xkpln1st + self%c1 * xkpln2nd + self%c2 * xkpln3rd) &
! All other potentials:
       + self%cdeltc * 2 * xk**2 + self%cdeltl * xk/p/2 * ( &
       (p + xk)**2 * (   log(  ( (p + xk)/self%scale )**2  ) - 1   ) - &
       (p - xk)**2 * (   log(  ( (p - xk)/self%scale )**2  ) - 1 ) ) + &
       self%cfullc * (p**2 + xk**2) * xkpln1st + self%cfulll * (p**2 + xk**2) * xk/p/4 * &
       (   log(  ( (p + xk)/self%scale )**2  )**2 - log(  ( (p - xk)/self%scale )**2  )**2   ) &
       + self%crm2 * xk/p * ( p + xk - abs(xk - p) )
    else
     vhat = phiint * (self%c0 * xkpln1st + self%c1 * xkpln2nd + self%c2 * xkpln3rd)
    endif

  end function vhat

! cccccccccccccccccccc

  subroutine sae(self, pp, w1, bb, a1, n)
    class (Toppik), intent(in)              :: self
    real (dp)   , intent(in) , dimension(n) :: pp, w1
    complex (dp), intent(in) , dimension(n) :: a1
    complex (dp), intent(out), dimension(n) :: bb
    integer     , intent(in)                :: n

     complex (dp) :: ff(n,n), cw(n), svw
     real (dp)    :: d
     integer      :: i, j, indx(n)

     do i = 1, 2 * n/3
        cw(i) = w1(i) * self%g0c( pp(i) ) / ( 4 * pi2 )
     end do

     do i = 2 * n/3 + 1, n
        cw(i) = w1(i) * self%g0c( pp(i) ) * pp(i)**2 / (4 * pi2)
     end do

     do i = 1, n

        if ( pp(i) < self%dcut .and. self%vflag == 1 ) then
            bb(i) = 1 + self%kincov * self%g0( pp(i) ) * ( pp(i)/self%tmass )**2/self%g0c( pp(i) )
        else
          bb(i) = (1,0)
        endif

        svw = (0,0)

        do j = 1, n
          if (i /= j) then
             ff(i,j) = - self%vhat( pp(i), pp(j) ) * cw(j);  svw = svw + ff(i,j)
          endif
        end do

        ff(i,i) = 1 - a1(i) - svw

     end do

     call zldcmp(ff, n, n, indx, d);  call zlbksb(ff, n, n, indx, bb)

  end subroutine sae

! cccccccccccccccccccc

  complex (dp) function g0(self, p)
    class (Toppik), intent(in) :: self
    real(dp)      , intent(in) :: p
      g0 = 1/cmplx(self%energy - self%vzero - p**2/self%tmass, self%tgamma, kind = dp)
  end function g0

! cccccccccccccccccccc

  complex (dp) function g0c(self, p)
    class (Toppik), intent(in) :: self
    real (dp)     , intent(in) :: p

    complex (dp) :: green, zk, zi, amd2k, aa, bb, cc, zzp, zzm, hypp, hypm

    if (self%gcflg == 0) then

      if (self%kinflg == 0) then
        g0c = self%g0(p)
      else if (self%kinflg == 1 .and. p <= self%dcut) then
        g0c = self%g0(p) * (1 + self%kincom) + self%kincoa
      else if (self%kinflg == 1 .and. p >= self%dcut) then
        g0c = self%g0(p) * (1 + self%kincom)
      else if (self%kinflg == 2 .and. p <= self%dcut) then
        g0c = self%g0(p) * (1 + self%kincom) + self%kincoa
      else if (self%kinflg == 2 .and. p >= self%dcut) then
        g0c = self%g0(p)
      else
        write(*,*) ' kinflg wrong! Stop.';  stop
      endif

    else if (self%gcflg == 1) then

      zi = (0,1)
      zk = - self%tmass * cmplx( self%energy, self%tgamma, kind = dp )
      zk = sqrt(zk);  amd2k = 4 * self%alphas * self%tmass/6/zk; aa = (2,0)
      bb = (1, 0) ; cc = 2 - amd2k;  zzp = (1 + zi * p/zk)/2
      zzm = (1 - zi * p/zk)/2

      if ( abs(zzp) >  20 ) then
        hypp = (1 - zzp)**(-aa) * hypgeo( aa, cc - bb, cc, zzp/(zzp - 1) )
      else
        hypp = hypgeo(aa,bb,cc,zzp)
      end if

      if ( abs(zzm) >  20 ) then
        hypm = (1 - zzm)**(-aa) * hypgeo( aa, cc - bb, cc, zzm/(zzm - 1) )
      else
        hypm = hypgeo(aa, bb, cc, zzm)
      endif

      green = - zi * self%tmass/(4 * p * zk)/(1 - amd2k) * (hypp - hypm)

! VZ anders herum als in Andres Konvention, da bei ihm G_0=1/[-E-i G+p^2/m]:

      g0c = - green;  if (p >  1.e3_dp * self%tmass) then
      write(*,*) ' g0cana = ', g0c,' not reliable. Stop.';  stop
    end if

    else
      write(*,*) ' gcflg wrong! Stop.';  stop
    end if

  end function g0c

! cccccccccccccccccccc

  function CrossSec(self) result(res)
    class (Toppik), intent(in) :: self
    real (dp), dimension(2)    :: res

    res = [self%xim, self%xdi]

  end function CrossSec

! cccccccccccccccccccc

  function pGrid(self) result(res)
    class (Toppik), intent(in) :: self
    real (dp), dimension(nmax)    :: res

    res = self%xpp

  end function pGrid

! cccccccccccccccccccc

  function Weight(self) result(res)
    class (Toppik), intent(in) :: self
    real (dp), dimension(nmax) :: res

    res = self%ww

  end function Weight

! cccccccccccccccccccc

  function Sigma(self) result(res)
    class (Toppik), intent(in) :: self
    real (dp), dimension(nmax) :: res

    res = self%dsdp

  end function Sigma

! cccccccccccccccccccc

  function Kernel(self) result(res)
    class (Toppik), intent(in)    :: self
    complex (dp), dimension(nmax) :: res

    res = self%zvfct

  end function Kernel

end module ToppikClass
