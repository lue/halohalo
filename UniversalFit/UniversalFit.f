      PROGRAM UniversalFit
c**********************************************************************
c  + Version 1.0: March 18, 2014
c----------------------------------------------------------------------
c
c Compute universal fitting functions for SHMF and SHVF
c
c----------------------------------------------------------------------
c  Author: Arthur Fangzhou Jiang                        Yale University  
c**********************************************************************

      INCLUDE 'paramfile.h'

      INTEGER   j

      REAL      z004,z025,z050,t050
      REAL      Ntau,Ntau_err
      
      REAL      fsub,a,cvir,Vvir,Vmax,fac_cvir
      REAL      x,xv
      
      INTEGER   k1,k2,k3,k4,k5,k6,k7
      CHARACTER id1,id2,id3,id4,id5,id6,id7
      CHARACTER outfil1*30,outfil2*30,outfil3*30,outfil4*30

      REAL      Delta_collapse,variance,zriddr
      EXTERNAL  Delta_collapse,variance,zriddr
      
      INTEGER   lblnk
      REAL      xH, Delta_crit,lookbacktime,time
      EXTERNAL  xH, Delta_crit,lookbacktime,time,lblnk
      
      REAL      zf,f_Ntau,dNdlog
      EXTERNAL  zf,f_Ntau,dNdlog

c**********************************************************************

c---read global and cosmological parameters

      CALL read_global_param

c---read parameters specific to the halo

      CALL read_halo_param

c---initialize cosmological parameters to be used throughout

      CALL init_cosmo

c---initialize and store mass variance

      CALL init_variance   

c---initialize and define parameters & variables to be used
c   (note that this subroutine must be called after init_variance)

      CALL init_param

c---compute formation redshifts z_f, where f=0.04, 0.25, 0.5

      z004 = zf(0.04)
      z025 = zf(0.25)
      z050 = zf(0.5)

      
c---compute halo dynamical age N_tau

      t050 = lookbacktime(z050)
      CALL qags(f_Ntau,tnow,t050,AbsAcc,RelAcc,Ntau,Ntau_err,N_eval,
     &   ErrFlag,N_MaxSubInt,N_MaxSubInt4,N_SubInt,WorkArr1,WorkArr2)
      
c---compute subhalo mass fraction fsub[(m/M0)>10^-4]
c   using Eq.(26) in Jiang & van den Bosch (2014) 
c   'Statistics of Dark Matter Substructure I'
      
      fsub = 0.3563/Ntau**0.6-0.075

c---compute the rescaling factor "a" in the fitting function
c   dN/dln(x) = gamma*(a*x)^alpha * exp[-beta*(a*x)^omega]

      Vvir = 159.43 * ((Mhalo/40.0/1.0E12)*(xH(z025)/xH_0))**(1.0/3.0)
      Vvir = Vvir * (Delta_crit(z025)/178.0)**(1.0/6.0)
      fac_cvir = time(z025) / (3.75 * time(z004))
      cvir = 4.0 * (1.0 + fac_cvir**8.4)**0.125
      Vmax = 0.465 * Vvir * SQRT(cvir/(ALOG(1.0+cvir)-cvir/(1.0+cvir)))
      a = 1.536 * Vvir/Vmax    

c---output useful info on screen

      WRITE(*,*)'Model predictions:'
      WRITE(*,*)'--------------------------------------------------'
      WRITE(*,*)'Halo properties:'
      WRITE(*,*)'Formation redshift z_f by which main progenitor'
      WRITE(*,*)'has assembled fM0'
      WRITE(*,*)'-> z_0.04 =',z004
      WRITE(*,*)'-> z_0.25 =',z025
      WRITE(*,*)'-> z_0.50 =',z050
      WRITE(*,*)'Dynamical Age:'
      WRITE(*,*)'-> N_tau =',Ntau
      WRITE(*,*)'   [number of dynamical times elapsed since z_0.50]'
      WRITE(*,*)'--------------------------------------------------'
      WRITE(*,*)'Parameters p1,p2,p3,p4,p5 of the fitting function' 
      WRITE(*,*)'dN/dln(x) = p1*(p5*x)^p2*exp[-p3*(p5*x)^p4],'
      WRITE(*,*)'where x = macc/M0, m/M0, Vacc/Vvir0, Vmax/Vvir0,'
      WRITE(*,*)'for USMF, ESMF, USVF, and ESVF, respectively:'
      WRITE(*,51)0.22,-0.91,6.0,3.0,1.0
      WRITE(*,52)0.31*fsub,-0.82,50.0,4.0,1.0
      WRITE(*,53)2.05,-3.2,2.2,13.0,a
      WRITE(*,54)5.45*fsub**1.4,-2.6,4.0,15.0,a
      WRITE(*,*)'--------------------------------------------------'
      WRITE(*,*)''
      
            
c---output the model functions for the unevolved/evolved SHMF/SHVF

      k1 = INT(xlgMhalo*100/1000)
      k2 = INT(MOD(xlgMhalo*100,1000)/100)
      k3 = INT(MOD(xlgMhalo*100,100)/10)
      k4 = MOD(xlgMhalo*100,10)
      id1 = CHAR(k1+ICHAR(null))
      id2 = CHAR(k2+ICHAR(null))
      id3 = CHAR(k3+ICHAR(null))
      id4 = CHAR(k4+ICHAR(null))
      
      k5 = INT(znow*100/100)
      k6 = INT(MOD(znow*100,100)/10)
      k7 = MOD(znow*100,10)
      id5 = CHAR(k5+ICHAR(null))
      id6 = CHAR(k6+ICHAR(null))
      id7 = CHAR(k7+ICHAR(null))

c---Unevolved Subhalo Mass Function (USMF)
      
      outfil1='/USMF_M'//id1//id2//id3//id4//'_z'//id5//id6//id7//'.dat'
      OPEN(10,file=moddir(1:lblnk(moddir))//outfil1(1:lblnk(outfil1)),
     &     status='UNKNOWN')

      p(1) = 0.22
      p(2) = -0.91
      p(3) = 6.0
      p(4) = 3.0
      p(5) = 1.0
      DO j=1,50
        x = binsize * FLOAT(j) - 5.0 - (binsize/2.0)
        x = 10.0**x
        WRITE(10,61)ALOG10(x),dNdlog(x)
      END DO

c---Evolved Subhalo Mass Function (ESMF)
      
      outfil2='/ESMF_M'//id1//id2//id3//id4//'_z'//id5//id6//id7//'.dat'
      OPEN(20,file=moddir(1:lblnk(moddir))//outfil2(1:lblnk(outfil2)),
     &     status='UNKNOWN')

      p(1) = 0.31*fsub
      p(2) = -0.82
      p(3) = 50.0
      p(4) = 4.0
      p(5) = 1.0
      DO j=1,50
        x = binsize * FLOAT(j) - 5.0 - (binsize/2.0)
        x = 10.0**x
        WRITE(20,61)ALOG10(x),dNdlog(x)
      END DO

c---Unevolved Subhalo Velocity Function (USVF)
      
      outfil3='/USVF_M'//id1//id2//id3//id4//'_z'//id5//id6//id7//'.dat'

      OPEN(30,file=moddir(1:lblnk(moddir))//outfil3(1:lblnk(outfil3)),
     &     status='UNKNOWN')
      p(1) = 2.05
      p(2) = -3.2
      p(3) = 2.2
      p(4) = 13.0
      p(5) = a
      DO j=1,50
        xv = (binsize/3.0) * FLOAT(j) - 1.0 - (binsize/6.0)
        xv = 10.0**xv
        WRITE(30,61)ALOG10(xv),dNdlog(xv)
      END DO  

c---Evolved Subhalo Velocity Function (ESVF)
      
      outfil4='/ESVF_M'//id1//id2//id3//id4//'_z'//id5//id6//id7//'.dat'

      OPEN(40,file=moddir(1:lblnk(moddir))//outfil4(1:lblnk(outfil4)),
     &     status='UNKNOWN')
      p(1) = 5.45*fsub**1.4
      p(2) = -2.6
      p(3) = 4.0
      p(4) = 15.0
      p(5) = a
      DO j=1,50
        xv = (binsize/3.0) * FLOAT(j) - 1.0 - (binsize/6.0)
        xv = 10.0**xv
        WRITE(40,61)ALOG10(xv),dNdlog(xv)
      END DO                        
      
c---close files

      CLOSE(10)
      CLOSE(20)
      CLOSE(30)
      CLOSE(40)

 51   FORMAT(' -> USMF:',4(F7.2,1X),F7.3)
 52   FORMAT(' -> ESMF:',4(F7.2,1X),F7.3) 
 53   FORMAT(' -> USVF:',4(F7.2,1X),F7.3)
 54   FORMAT(' -> ESVF:',4(F7.2,1X),F7.3)
 61   FORMAT(F7.4,2X,E12.5)

      STOP
      END


c**********************************************************************

      REAL FUNCTION zf(f)
c----------------------------------------------------------------------
c
c Calculate formation redshift zf by which the main progenitor has 
c assembled a fraction f of its present-day mass M0.
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'

      REAL     f,wtildef,alphaf,Sf,Wf
      
      REAL     Watz
      COMMON /zatW/ Watz
      
      REAL     Delta_collapse,variance,z_at_W,zriddr
      EXTERNAL Delta_collapse,variance,z_at_W,zriddr

c---
      alphaf = 0.815*EXP(-2.0*f**3.0)/f**0.707      
      wtildef = SQRT(2.0*ALOG(1.0+alphaf))
      Sf = variance(f*Mhalo)**2.0
      Wf = Wnow + wtildef*SQRT(Sf-Shalo)

c---find the redshift (zf) corresponding to Wf

      Watz = Wf
      zf = zriddr(z_at_W,znow,20.0,1.0E-4)

      END

c**********************************************************************

      REAL FUNCTION f_Ntau(t)
c----------------------------------------------------------------------
c
c f(t) = 1/tdyn(t) -- the integrand for Ntau,
c where tdyn [Gyr] is the dynamical timescale at lookback time t.
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'
      
      REAL     t,z,tdyn0,tdyn,fac
      
      REAL     tatz
      COMMON/zatt/ tatz
      
      REAL     Delta_crit,xH,z_at_t,zriddr
      EXTERNAL Delta_crit,xH,z_at_t,zriddr

c---
      tatz = t
      z = zriddr(z_at_t,znow,20.0,1.0E-4)
      tdyn0 = (pi/1.4142)*(deltacrit0**(-0.5))*xH_0_recip ![Gyr]
      fac = SQRT(Delta_crit(z)/deltacrit0) * xH(z)/xH_0
      tdyn = tdyn0 / fac
      
      f_Ntau = 1.0/tdyn

      END

c**********************************************************************

      REAL FUNCTION dNdlog(x)
c----------------------------------------------------------------------
c
c universal fitting function:
c dN/dlog(x) = ln(10) * gamma*(a*x)^alpha * exp[-beta*(a*x)^omega],
c where x is macc/M0, m/M0, Vacc/Vvir0, or Vmax/Vvir0 for 
c unevolved/evolved SHMF/SHVF respectively.
c
c----------------------------------------------------------------------

      INCLUDE 'paramfile.h'
      
      REAL     x

c---

      dNdlog = lnten*p(1)*(p(5)*x)**p(2)*EXP(-p(3)*(p(5)*x)**p(4))

      END

c**********************************************************************