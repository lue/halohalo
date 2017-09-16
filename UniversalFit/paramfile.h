c=======================================================================
c
c                      INCLUDE FILE mah.h
c
c=======================================================================
c  Written by: Arthur Fangzhou Jiang
c=======================================================================

      IMPLICIT NONE

c=======================================================================
c                        PARAMETERS
c=======================================================================


c---the constant pi

      REAL  pi
      PARAMETER (pi=3.141592654)

c---the constant SQRT(2*pi)

      REAL  sqrt2pi
      PARAMETER (sqrt2pi=2.506628275)

c---the constant SQRT(2/pi)

      REAL  sqrt2divpi
      PARAMETER (sqrt2divpi=0.79788456)

c---the constant Gamma(3/2) = sqrt(pi)/2

      REAL  gamma32
      PARAMETER (gamma32=0.886226925)

c---the constant ln(10)

      REAL  lnten
      PARAMETER (lnten=2.302585093)

c---a small number

      REAL*8 deps
      PARAMETER (deps=1.0E-10)

c---The gravitational constant in M_sun^-1 Mpc (km/s)^2

      REAL  gee
      PARAMETER (gee = 4.2994E-9)

c---The binsize (dex) of subhalo mass used for plotting

      REAL binsize
      PARAMETER (binsize=0.1)     

c---Number of interpolation points for mass variance

      INTEGER Nsigma
      REAL    Mminvar, Mmaxvar
      PARAMETER (Nsigma = 2000)
      PARAMETER (Mminvar = 4.0)
      PARAMETER (Mmaxvar =15.5)

c---parameters needed for quadpack integration routines [Arthur]

      INTEGER   N_MaxSubInt,N_MaxSubInt4
      REAL      AbsAcc,RelAcc
      PARAMETER (N_MaxSubInt=10000)
      PARAMETER (N_MaxSubInt4=41000) ! no less than N_MaxSubInt*4
      PARAMETER (AbsAcc=1.0E-4)
      PARAMETER (RelAcc=1.0E-4)

c---parameters needed for quadpack integration routines  [Frank]

      INTEGER Nlimit,Nlenw
      PARAMETER (Nlimit=10000)
      PARAMETER (Nlenw =41000)

c---the null character

      CHARACTER null*1
      PARAMETER (null='0')

c=======================================================================
c                        COMMON BLOCKS
c=======================================================================

c---cosmology

      REAL  omega_0,omega_lambda,sigma8,xhubble
      REAL  nspec,omega_b_h2,c8,Mf_WDM
      COMMON /cosmo_model/ omega_0,omega_lambda,sigma8,xhubble,
     &                     nspec,omega_b_h2,c8,Mf_WDM

      REAL  omega_b,f_bar,rho_aver_com
      COMMON /cosmo_baryon/ omega_b,f_bar,rho_aver_com

      INTEGER  ivar
      REAL     xH_0,xH_0_recip,tff0,t0,Gamma_cosmo,cstlbt1,cstlbt2
      REAL     rho_crit_0,sigma8_norm,deltacrit0
      LOGICAL  BBKS,EBW,EISHU
      COMMON /cosmo_param/ xH_0,xH_0_recip,tff0,t0,Gamma_cosmo,
     &             cstlbt1,cstlbt2,ivar,rho_crit_0,sigma8_norm,
     &             deltacrit0,BBKS,EBW,EISHU

c---Eisenstein & Hu power spectrum

      REAL     sEH,bnode,ksilk,keq,alpha_c,alpha_b,beta_c,beta_b
      COMMON /ehpar/ sEH,bnode,ksilk,keq,alpha_c,alpha_b,beta_c,beta_b

c---mass variance

      REAL   vectorM(Nsigma),vectorS(Nsigma)
      REAL   vectorMM(Nsigma),vectorSS(Nsigma)
      REAL   vectorM2(Nsigma),vectorS2(Nsigma)
      REAL   vectorD(Nsigma),vectorZ(Nsigma),vectorZ2(Nsigma)
      COMMON /spl_sigma/ vectorM,vectorS,vectorSS,vectorMM,vectorS2,
     &                   vectorM2,vectorD,vectorZ,vectorZ2


c---parameters regarding host halo
 
      REAL    xlgMhalo,Mhalo,Shalo
      REAL    znow,Wnow,tnow
      COMMON /halopar/ xlgMhalo,Mhalo,Shalo,znow,Wnow,tnow
      
c---parameters in the universal fitting function
c   dN/dln(x) = gamma*(a*x)^alpha * exp[-beta*(a*x)^omega]

      REAL    p(5)
      COMMON /UniversalFitParams/ p

c---parameters for quadpack integration routines         [Arthur]

      INTEGER N_eval,ErrFlag,N_SubInt,WorkArr1(N_MaxSubInt)
      REAL*8  WorkArr2(N_MaxSubInt4)
      COMMON /forQuadPack/ WorkArr2,WorkArr1,N_eval,N_SubInt,ErrFlag
      
c---work space required for quadpack integration routines [Frank]

      INTEGER Neval,ierr,last,iwork(Nlimit)
      REAL*8  work(Nlenw)
      COMMON /quadpackint/ work,iwork,Neval,ierr,last

c---Directories etc.

      CHARACTER moddir*60
      COMMON /modeldir/ moddir


c=======================================================================
c                             END
c=======================================================================





