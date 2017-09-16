c**********************************************************************
c                SUBROUTINES WITH INITIALIZATIONS
c
c
c  Fangzhou Arthur Jiang                              March 2014
c**********************************************************************

      SUBROUTINE read_global_param
c---------------------------------------------------------------------------
c
c  Subroutine to read in global parameters
c
c---------------------------------------------------------------------------
 
      INCLUDE 'paramfile.h'

      INTEGER      ipow

      WRITE(*,*)'User inputs are:'
      WRITE(*,*)'--------------------------------------------------'
      WRITE(*,*)'Cosmology :'

c---read parameters of cosmological parameters

      READ(*,*)omega_0
      WRITE(*,*)'-> Omega_m,0 = ',omega_0

      omega_lambda = 1.0 - omega_0
      WRITE(*,*)'-> Omega_Lambda,0 = ',omega_lambda

      READ(*,*)xhubble
      WRITE(*,*)'-> h (=H_0/100) = ',xhubble

      READ(*,*)sigma8
      WRITE(*,*)'-> sigma_8 = ',sigma8

      READ(*,*)nspec
      WRITE(*,*)'-> n_spec = ',nspec
      WRITE(*,*)'   [Harrisson-Zeldovich = 1.0]'

      READ(*,*)omega_b_h2
      WRITE(*,*)'-> Omega_b*h^2 = ',omega_b_h2
   
c---decide which power spectrum fitting function to use
c     BBKS = Bardeen, Bond, Kaiser & Szalay, 1986, ApJ, 304, 15
c     EBW  = Efstathiou, Bond & White, 1992, MNRAS, 258, 1
c     EH   = Eisenstein & Hu, 1998, ApJ, 496, 605

      READ(*,*)ipow
      WRITE(*,*)'-> power spectrum = ',ipow
      WRITE(*,*)'   [choice of the power spectrum' 
      WRITE(*,*)'   fitting function T(K) to use:' 
      WRITE(*,*)'   BBKS (1), EBW (2), EH (3)]'

      IF (ipow.NE.1.and.ipow.NE.2.and.ipow.NE.3) THEN
        CALL Terminate('Invalid ipow')
      END IF
 
      BBKS  = .FALSE.
      EBW   = .FALSE.
      EISHU = .FALSE.
      IF (ipow.EQ.1) BBKS  = .TRUE.
      IF (ipow.EQ.2) EBW   = .TRUE.
      IF (ipow.EQ.3) EISHU = .TRUE.

c---allow for the possibility of WDM by specifying the filter mass

      READ(*,*)Mf_WDM
      WRITE(*,*)'-> Mf_WDM = ',Mf_WDM
      WRITE(*,*)'   [WDM filter mass (0.0 = CDM)]'

c---indicate whether you want to compute the mass variance numerically
c   (slow) or using the fitting function (fast)

      READ(*,*)ivar
      WRITE(*,*)'-> method = ',ivar
      WRITE(*,*)'   [method of computing mass variance (1=num, 2=fit)]'

      IF (ivar.NE.1.AND.ivar.NE.2) THEN
        CALL Terminate('Invalid ivar')
      END IF

      RETURN
      END

c**********************************************************************

      SUBROUTINE read_halo_param
c---------------------------------------------------------------------------
c
c  Subroutine to read in parameters specific to the merger tree(s)
c
c---------------------------------------------------------------------------
 
      INCLUDE 'paramfile.h'

      WRITE(*,*)'Halo info:'

c---read present day halo mass

      READ(*,*)xlgMhalo
      WRITE(*,*)'-> Halo mass M0 = 10^',xlgMhalo,'[Msun/h]'
      Mhalo = 10.0**xlgMhalo

c---read redshift that corresponds to `present'

      READ(*,*)znow
      WRITE(*,*)'-> Obervation redshift z0 = ',znow
      
c---give directory to which output should be written

      WRITE(*,*)'output directory:'
      READ(*,'(A)')moddir
      WRITE(*,'(A)')moddir
      WRITE(*,*)'--------------------------------------------------'
      WRITE(*,*)' '
          
      RETURN
      END

c**********************************************************************

      SUBROUTINE init_cosmo    
c---------------------------------------------------------------------------
c
c  Subroutine to initialize cosmology related stuff
c
c---------------------------------------------------------------------------
 
      INCLUDE 'paramfile.h'

c---parameter that sets photon temperature (needed for EH transfer function)

      REAL   theta
      PARAMETER (theta = 1.0093)

      REAL     sTH,sGS,sSK
      REAL     f,bb1,bb2,zeq,zd,y,Fy,Gy,Req,Rd,a1,a2,b1,b2

      REAL     time,Delta_crit,XEXP
      EXTERNAL time,Delta_crit,XEXP
     
c---

c---Calculate H0 in km/s/Mpc and its reciprocal in Gyr. From the latter
c   we then calculate the age of the Universe, t0, in Gyr

      xH_0 = 100.0 * xhubble
      xH_0_recip = 1.0 / (xH_0 * 1.023E-3)            
  
c---calculate free-fall time at z=0 for an overdensity of 200
 
      tff0 = (pi/SQRT(800.0 * omega_0)) * xH_0_recip

c---calculate parameters used to speed up calculation of (lookback)time

      cstlbt1 = SQRT((1.0-omega_0)/omega_0)
      cstlbt2 = 2.0 / (3.0*SQRT(1.0-omega_0))

      t0 = time(0.0) * xH_0_recip

c---set critical density at z=0 and comoving mass density [h^2 Msun Mpc^{-3}]

      rho_crit_0 = 3.0E+4 / (8.0 * pi * gee)
      rho_aver_com = rho_crit_0 * omega_0

c---calculate critical density for collapse at z=0

      deltacrit0 = Delta_crit(0.0)

c---calculate baryon density and baryonic mass fraction

      omega_b = omega_b_h2/(xhubble**2.0)
      f_bar = omega_b/omega_0

c---define Gamma: two options here; i) go with the standard definition
c   (Gamma = Omega_0 h), which is for negligible baryon mass, or use
c   the baryonic correction from Sugiyama (1995).

c---with baryonic correction

c      Gamma_cosmo = omega_0 * xhubble * 
c     &              XEXP(-omega_b - SQRT(2.0*xhubble) * f_bar)

c---without baryonic correction

      Gamma_cosmo = omega_0 * xhubble

c---define a number of parameters needed to compute the Eisenstein & Hu
c   power specrum

      f = omega_0 * xhubble**2

      bb1 = 0.313 * f**(-0.419) * (1.0 + 0.607*f**0.674)
      bb2 = 0.238 * f**(0.223)
 
      bnode = 8.41*f**0.435

      keq = 7.46E-2 * f / theta**2.0
      ksilk = 1.6 * (omega_b_h2)**0.52 * f**0.73 * 
     &       (1.0 + (10.4*f)**(-0.95))

      zeq = 2.5E+4 * f / theta**4.0
      zd = 1291.0 * ((f**0.251)/(1.0 + 0.659*f**0.828)) *
     &               (1.0 + bb1 * omega_b_h2**bb2)

      y = ((1.0+zeq)/(1.0+zd))
      Fy = ALOG((SQRT(1.0+y) + 1.0) / (SQRT(1.0+y) - 1.0))
      Gy = y * (-6.0 * SQRT(1.0+y) + (2.0+3.0*y) * Fy)

      Req = 31.5 * omega_b_h2 * (1000.0/zeq) / theta**4.0
      Rd  = 31.5 * omega_b_h2 * (1000.0/zd) / theta**4.0

      sEH = (2.0/(3.0*keq)) * SQRT(6.0/Req) *
     &      ALOG((SQRT(1.0+Rd) + SQRT(Rd+Req))/(1.0+SQRT(Req)))
        
      a1 = (46.9*f)**(0.670) * (1.0 + (32.1*f)**(-0.532))
      a2 = (12.0*f)**(0.424) * (1.0 + (45.0*f)**(-0.582))
      b1 = 0.944 / (1.0+(458.0*f)**(-0.708))
      b2 = (0.395*f)**(-0.0266)

      alpha_c = a1**(-f_bar) * a2**(-(f_bar**3))
      beta_c = 1.0 + b1*(((omega_0-omega_b)/omega_0)**b2 - 1.0)
      beta_c = 1.0/beta_c

      alpha_b = 2.07 * keq * sEH * (1.0+Rd)**(-0.75) * Gy
      beta_b = 0.5 + f_bar + 
     &   (3.0-2.0*f_bar) * SQRT((17.2*f)**2 + 1.0)

      RETURN
      END

c**********************************************************************

      SUBROUTINE init_variance
c---------------------------------------------------------------------------
c
c  Subroutine to initialize mass variance
c
c---------------------------------------------------------------------------
 
      INCLUDE 'paramfile.h'

      INTEGER  i,ivtemp
      REAL     xM,z,yp1,ypn,xp1,xpn

      REAL     variance,Delta_collapse
      EXTERNAL variance,Delta_collapse

c---

c---compute un-normalized rms mass variance inside spherical
c   shell with radius Rf = 8 h^{-1} Msun
c   This is used to normalize the power-spectrum to sigma8.

      WRITE(*,*)'Computing Mass Variance ... ...'

      ivtemp = ivar
      ivar = 1

      sigma8_norm = sigma8
      sigma8_norm = variance(5.9543E+14 * omega_0)
      c8 = 1.0

      ivar = ivtemp

c---from now, `variance' is normalized to sigma8

      WRITE(*,*)'... ... Mass variance is now normalized to sigma_8.'
      WRITE(*,*)' '
      
c---set up the mass variance on a grid. The grid is a one-D vector,
c   for which we compute the mass variance numerically. The grid
c   consistes of Nsigma points with 5 <= log(M) <= 18.0

      ivtemp = ivar
      ivar = 1

      DO i=1,Nsigma
        vectorM(i) = Mminvar + 
     &       FLOAT(i-1)/FLOAT(Nsigma-1) * (Mmaxvar-Mminvar)
        vectorZ(i) = FLOAT(i-1)/FLOAT(Nsigma-1) * 100.0
        xM = 10.0**vectorM(i)
        z = vectorZ(i)
        vectorS(i) = variance(xM)
        vectorD(i) = Delta_collapse(z)
      END DO      

c---compute the derivatives at the two ends of one-D grids

      yp1 = (vectorS(2) - vectorS(1)) / 
     &      (vectorM(2) - vectorM(1))
      ypn = (vectorS(Nsigma) - vectorS(Nsigma-1)) / 
     &      (vectorM(Nsigma) - vectorM(Nsigma-1))

      xp1 = (vectorZ(2) - vectorZ(1)) /
     &      (vectorD(2) - vectorD(1))
      xpn = (vectorZ(Nsigma) - vectorZ(Nsigma-1)) / 
     &      (vectorD(Nsigma) - vectorD(Nsigma-1))

c---and compute the spline coefficients, to be used for spline interpolation
c   note that we compute the spline coefficients both ways!

      DO i=1,Nsigma
        vectorSS(i) = vectorS(Nsigma+1-i)
        vectorMM(i) = vectorM(Nsigma+1-i)
      END DO

      CALL spline(vectorM,vectorS,Nsigma,yp1,ypn,vectorS2)
      CALL spline(vectorSS,vectorMM,Nsigma,2.0E+30,2.0E+30,vectorM2)
      CALL spline(vectorD,vectorZ,Nsigma,xp1,xpn,vectorZ2)
      ivar = ivtemp

      RETURN
      END

c**********************************************************************

      SUBROUTINE init_param
c---------------------------------------------------------------------------
c
c  Subroutine to initialize parameters and variables
c
c---------------------------------------------------------------------------
 
      INCLUDE 'paramfile.h'
      
      REAL     variance,Delta_collapse,lookbacktime
      EXTERNAL variance,Delta_collapse,lookbacktime
c---  

      Shalo = variance(Mhalo)**2.0
      Wnow = Delta_collapse(znow)
      tnow = lookbacktime(znow)

      RETURN
      END

c**********************************************************************


