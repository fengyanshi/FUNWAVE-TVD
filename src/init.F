!------------------------------------------------------------------------------------
!
!      FILE init.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!    ALLOCATE_VARIABLES is a subroutine to allocate variables
!
!    HISTORY:
!    05/01/2010 Fengyan Shi 
!    09/26/2013 Babak Tehranirad, added spacially varying cd
!    08/06/2015 Choi added V4xL, U4yL, V4xL, U4yL, WAVEMAKER_VIS
!
!-------------------------------------------------------------------------------------
SUBROUTINE ALLOCATE_VARIABLES
     USE GLOBAL

! coordinate for cartesian only
     ALLOCATE (Xco(Mloc),Yco(Nloc))
! allocate variables
     ALLOCATE(DelxU(Mloc,Nloc),DelxHU(Mloc,Nloc),DelxV(Mloc,Nloc),DelxEtar(Mloc,Nloc),&
              DelyU(Mloc,Nloc),DelyHV(Mloc,Nloc),DelyV(Mloc,Nloc),DelyEtar(Mloc,Nloc),&
              DelxHV(Mloc,Nloc),DelyHU(Mloc,Nloc), &
! U V HU H in x-direction
              UxL(Mloc1,Nloc),UxR(Mloc1,Nloc),VxL(Mloc1,Nloc),VxR(Mloc1,Nloc),&
              HUxL(Mloc1,Nloc),HUxR(Mloc1,Nloc),HVxL(Mloc1,Nloc),HVxR(Mloc1,Nloc), &
              HxL(Mloc1,Nloc),HxR(Mloc1,Nloc), &
! U V HV H in y-direction
              UyL(Mloc,Nloc1),UyR(Mloc,Nloc1),VyL(Mloc,Nloc1),VyR(Mloc,Nloc1),&
              HVyL(Mloc,Nloc1),HVyR(Mloc,Nloc1),HUyL(Mloc,Nloc1),HUyR(Mloc,Nloc1), &
              HyL(Mloc,Nloc1),HyR(Mloc,Nloc1), &
! cross-derivatives
              Uxy(Mloc,Nloc),Vxy(Mloc,Nloc),DUxy(Mloc,Nloc),DVxy(Mloc,Nloc), &
! second-derivatives
              Uxx(Mloc,Nloc),Vyy(Mloc,Nloc),DUxx(Mloc,Nloc),DVyy(Mloc,Nloc), &
! 1st-derivatives
              Ux(Mloc,Nloc),Uy(Mloc,Nloc),Vx(Mloc,Nloc),Vy(Mloc,Nloc), &
              DUx(Mloc,Nloc),DUy(Mloc,Nloc),DVx(Mloc,Nloc),DVy(Mloc,Nloc), &
              ETAT(Mloc,Nloc),ETAx(Mloc,Nloc),ETAy(Mloc,Nloc), &
              ETATx(Mloc,Nloc),ETATy(Mloc,Nloc), &
! time-derivatives
              U0(Mloc,Nloc),V0(Mloc,Nloc),Ut(Mloc,Nloc),Vt(Mloc,Nloc),&
              Utx(Mloc,Nloc),Vty(Mloc,Nloc),Utxx(Mloc,Nloc),Utxy(Mloc,Nloc),&
              Vtxy(Mloc,Nloc),Vtyy(Mloc,Nloc),&
              DUtxx(Mloc,Nloc),DUtxy(Mloc,Nloc),&
              DVtxy(Mloc,Nloc),DVtyy(Mloc,Nloc),DUtx(Mloc,Nloc),DVty(Mloc,Nloc),&
! P Q Eta, Fx, Fy
              PL(Mloc1,Nloc),PR(Mloc1,Nloc),QL(Mloc,Nloc1),QR(Mloc,Nloc1), &
              FxL(Mloc1,Nloc),FxR(Mloc1,Nloc),FyL(Mloc,Nloc1),FyR(Mloc,Nloc1), &
              GxL(Mloc1,Nloc),GxR(Mloc1,Nloc),GyL(Mloc,Nloc1),GyR(Mloc,Nloc1), &
              EtaRxL(Mloc1,Nloc),EtaRxR(Mloc1,Nloc), &
              EtaRyL(Mloc,Nloc1),EtaRyR(Mloc,Nloc1), &
! sponge
              SPONGE(Mloc,Nloc), SpongeMaker(Mloc,Nloc), &
! original variables at notes
              Fx(Mloc1,Nloc),Fy(Mloc,Nloc1),&
              U(Mloc,Nloc),V(Mloc,Nloc), HU(Mloc,Nloc),HV(Mloc,Nloc),&
              Gx(Mloc1,Nloc),Gy(Mloc,Nloc1), &
              P(Mloc1,Nloc),Q(Mloc,Nloc1), &
              SxL(Mloc1,Nloc),SxR(Mloc1,Nloc), &
              SyL(Mloc,Nloc1),SyR(Mloc,Nloc1),SourceX(Mloc,Nloc), &
              SourceY(Mloc,Nloc), &
! others
              Umean(Mloc,Nloc),Vmean(Mloc,Nloc),ETAmean(Mloc,Nloc),&
              Usum(Mloc,Nloc),Vsum(Mloc,Nloc),ETAsum(Mloc,Nloc), &
              nu_smg(Mloc,Nloc), &
              Num_Zero_Up(Mloc,Nloc), &
              WaveHeightRMS(Mloc,Nloc),  &
              WaveHeightAve(Mloc,Nloc),  &
              Emax(Mloc,Nloc),  &
              Emin(Mloc,Nloc), &
              HrmsSum(Mloc,Nloc), &
              HavgSum(Mloc,Nloc), &
              !ykchoi
              ETA2sum(Mloc,Nloc), ETA2mean(Mloc,Nloc), &
                  SigWaveHeight(Mloc,Nloc),  &
# if defined (CARTESIAN)
              U4xL(Mloc1,Nloc),U4xR(Mloc1,Nloc),&
              V4yL(Mloc,Nloc1),V4yR(Mloc,Nloc1), &
      ! ykchoi added V4xL and U4yL (08/06/15)
                  V4xL(Mloc1,Nloc),V4xR(Mloc1,Nloc),&  
                  U4yL(Mloc,Nloc1),U4yR(Mloc,Nloc1) & 
# else
# if defined (ZALPHA)
              U4xL(Mloc1,Nloc),U4xR(Mloc1,Nloc),&
              V4yL(Mloc,Nloc1),V4yR(Mloc,Nloc1), &
        ! ykchoi added V4xL and U4yL (08/06/15)
            V4xL(Mloc1,Nloc),V4xR(Mloc1,Nloc),&  
            U4yL(Mloc,Nloc1),U4yR(Mloc,Nloc1), &  
# endif
              Dx(Mloc,Nloc),Dy(Mloc,Nloc), Coriolis(Mloc,Nloc), &
              SlopeX(Mloc,Nloc),SlopeY(Mloc,Nloc), &
              Lat_theta(Mloc,Nloc) &
# endif
              )
      ALLOCATE(Depth(Mloc,Nloc),H(Mloc,Nloc),&
               Depthx(Mloc1,Nloc),Depthy(Mloc,Nloc1), &
               MASK(Mloc,Nloc),DepthNode(Mloc1,Nloc1), &
               MASK_STRUC(Mloc,Nloc),MASK9(Mloc,Nloc), &
               tmp4preview(Mloc,Nloc),Int2Flo(Mloc,Nloc),&
               Cd(Mloc,Nloc),CD_breakwater(Mloc,Nloc) &
              )
! updating variables
      ALLOCATE(Eta(Mloc,Nloc),Eta0(Mloc,Nloc), &
               Ubar0(Mloc,Nloc),Vbar0(Mloc,Nloc),&
               Ubar(Mloc,Nloc),Vbar(Mloc,Nloc))
# if defined(UseEtaScreen)
! wetting and drying
      ALLOCATE(EtaScreen(Mloc,Nloc))
# endif
! dispersion updating variables
# if defined (CARTESIAN)
      ALLOCATE(U4(Mloc,Nloc),V4(Mloc,Nloc),U1p(Mloc,Nloc), & 
               V1p(Mloc,Nloc),U1pp(Mloc,Nloc),V1pp(Mloc,Nloc),&
               U2(Mloc,Nloc),V2(Mloc,Nloc),U3(Mloc,Nloc),V3(Mloc,Nloc))
# else
# if defined (ZALPHA)
      ALLOCATE(U4(Mloc,Nloc),V4(Mloc,Nloc) & 
               )
# endif
      ALLOCATE(U1p(Mloc,Nloc), & 
               V1p(Mloc,Nloc) )
# endif
# if defined (ITERATION)
      ALLOCATE(Ui(Mloc,Nloc), Vi(Mloc,Nloc),ETAi(Mloc,Nloc),&
        UbarOld(Mloc,Nloc),VbarOld(Mloc,Nloc),EtaOld(Mloc,Nloc) )
# endif

  ! HeightMax will be used not only in output but also meteo module
        ALLOCATE(HeightMax(Mloc,Nloc))
        HeightMax=ZERO

      ALLOCATE(WaveMaker_Mass(Mloc,Nloc))
      WaveMaker_Mass = ZERO

      IF(WAVEMAKER(1:7)=='WK_TIME')THEN
        ALLOCATE(WAVE_COMP(NumWaveComp,3),Beta_genS(NumWaveComp),D_genS(NumWaveComp) )
      ENDIF

      IF(VISCOSITY_BREAKING.OR.SHOW_BREAKING)THEN
       ALLOCATE(AGE_BREAKING(Mloc,Nloc))
       ALLOCATE(nu_break(Mloc,Nloc))
       nu_break=nu_bkg
      ENDIF

      ALLOCATE(ROLLER_FLUX(Mloc,Nloc),UNDERTOW_U(Mloc,Nloc),UNDERTOW_V(Mloc,Nloc))

      IF(WAVEMAKER_VIS)THEN
       ALLOCATE(nu_break(Mloc,Nloc))
       nu_break=ZERO
      ENDIF

      IF(DIFFUSION_SPONGE)THEN
       ALLOCATE(nu_sponge(Mloc,Nloc))
       nu_sponge=ZERO
      ENDIF

      IF(OUT_Hmin)THEN
        ALLOCATE(HeightMin(Mloc,Nloc))
        HeightMin=ZERO
      ENDIF
      IF(OUT_Umax)THEN
        ALLOCATE(VelocityMax(Mloc,Nloc))
        VelocityMax=ZERO
      ENDIF
      IF(OUT_VORmax)THEN
        ALLOCATE(VorticityMax(Mloc,Nloc))
        VorticityMax=ZERO
      ENDIF
      IF(OUT_MFmax)THEN
        ALLOCATE(MomentumFluxMax(Mloc,Nloc))
        MomentumFluxMax=ZERO
      ENDIF
      IF(OUT_Time)THEN
        ALLOCATE(ARRTIME(Mloc,Nloc))
        ARRTIME=ZERO
      ENDIF      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Oct, 2021
!Salatin, R., Chen, Q., Bak, A. S., Shi, F., & Brandt, S. R. (2021). Effects of
!wave coherence on longshore variability of nearshore wave processes. Journal
!of Geophysical Research: Oceans,126, e2021JC017641.
!https://doi.org/10.1029/2021JC017641

! allocate some variables

    ALLOCATE(xmk_wk(Mloc),ymk_wk(Nloc))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE ALLOCATE_VARIABLES

!-------------------------------------------------------------------------------------
!
!    INITIALIZATION is subroutine for initialization
!
!    HISTORY:
!    05/01/2010 Fengyan Shi
!    09/26/2013 Babak Tehranirad, added varying cd
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIALIZATION
     USE GLOBAL
     USE INPUT_READ
     USE BATHY_CORRECTION_MODULE
# if defined(INTEL)
     USE IFPORT
# endif
     IMPLICIT NONE
     INTEGER :: VTYPE
# if defined (PARALLEL)
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGlob
     REAL(SP) :: myvar_tmp
# endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(SP) :: DXg,DYg                                   ! Salatin et al. 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if !defined (CARTESIAN)
     REAL(SP)::Dphi_r,Dtheta_r

! for spherical coordinates and Coriolis parameter

   IF(StretchGrid) THEN
# if defined (PARALLEL)
       CALL GetFile(DX_FILE,Dx)
       CALL GetFile(DY_FILE,Dy)
       CALL GetFile(Coriolis_FILE,Coriolis)
# else
      OPEN(1,FILE=TRIM(DX_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Dx(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(DY_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Dy(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(Coriolis_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Coriolis(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        Dx(I,J)=Dx(I,Jbeg)
        Dy(I,J)=Dy(I,Jbeg)
        Coriolis(I,J)=Coriolis(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        Dx(I,J)=Dx(I,Jend)
        Dy(I,J)=Dy(I,Jend)
        Coriolis(I,J)=Coriolis(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        Dx(I,J)=Dx(Ibeg,J)
        Dy(I,J)=Dy(Ibeg,J)
        Coriolis(I,J)=Coriolis(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        Dx(I,J)=Dx(Iend,J)
        Dy(I,J)=Dy(Iend,J)
        Coriolis(I,J)=Coriolis(Iend,J)
       ENDDO
     ENDDO  

# endif

    ELSE ! for stretch grid

     Dphi_r=Dphi*pi/180.0_SP
     Dtheta_r=Dtheta*pi/180.0_SP

     Do J=1,Nloc
     Do I=1,Mloc
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!       Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
!                       +(npy*Nglob/py+J-1)*Dtheta_r
       Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                       +((jjsta - 1)+J-1)*Dtheta_r
!---ykchoi Jan/23/2018]
# else
       Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                       +(J-1)*Dtheta_r
# endif
       Dx(I,J) = R_earth*Dphi_r*COS(Lat_theta(I,J))
       Dy(I,J) = R_earth*Dtheta_r
       Coriolis(I,J) = pi*SIN(Lat_theta(I,J)) / 21600.0_SP
     ENDDO
     ENDDO
    ENDIF ! stretchgrid

# endif

! parameter kappa for order of MUSCL
     IF(HIGH_ORDER(1:3)=='SEC')THEN
      Kappa = -1.0_SP
     ELSE
      Kappa = 1.0_SP/3.0_SP
     ENDIF

! set zeros
     T_sum = ZERO
     Umean = ZERO
     Vmean = ZERO
     ETAmean = ZERO
     nu_smg = ZERO
     Num_Zero_Up = 0
     WaveHeightRMS = ZERO
     WaveHeightAve =ZERO 
     Emax = ZERO
     Emin = ZERO
     HrmsSum = ZERO
     HavgSum = ZERO
     DelxU=0.0_SP
     DelxHU=0.0_SP
     DelxV=0.0_SP
     DelxEtar=0.0_SP
     DelyU=0.0_SP
     DelyHV=0.0_SP
     DelyV=0.0_SP
     DelyEtar=0.0_SP
     DelxHV=0.0_SP
     DelyHU=0.0_SP
     UxL=0.0_SP
     UxR=0.0_SP
     VxL=0.0_SP
     VxR=0.0_SP
     HUxL=0.0_SP
     HUxR=0.0_SP
     HVxL=0.0_SP
     HVxR=0.0_SP
     HxL=0.0_SP
     HxR=0.0_SP
     UyL=0.0_SP
     UyR=0.0_SP
     VyL=0.0_SP
     VyR=0.0_SP
     HVyL=0.0_SP
     HVyR=0.0_SP
     HUyL=0.0_SP
     HUyR=0.0_SP
     HyL=0.0_SP
     HyR=0.0_SP
# if defined (CARTESIAN)
     U4xL=ZERO
     U4xR=ZERO
     V4yL=ZERO
     V4yR=ZERO
# endif
     Uxy=ZERO
     Vxy=ZERO
     DUxy=ZERO
     DVxy=ZERO
     Uxx=ZERO
     Vyy=ZERO
     DUxx=ZERO
     DVyy=ZERO 
     U0=ZERO
     V0=ZERO
     Ut=ZERO
     Vt=ZERO
     Utx=ZERO
     Vty=ZERO
     Utxx=ZERO
     Utxy=ZERO
     Vtxy=ZERO
     Vtyy=ZERO
     DUtxx=ZERO
     DUtxy=ZERO
     DVtxy=ZERO
     DVtyy=ZERO
     DUtx=ZERO
     DVty=ZERO    
     PL=0.0_SP
     PR=0.0_SP
     QL=0.0_SP
     QR=0.0_SP
     FxL=0.0_SP
     FxR=0.0_SP
     FyL=0.0_SP
     FyR=0.0_SP
     GxL=0.0_SP
     GxR=0.0_SP
     GyL=0.0_SP
     GyR=0.0_SP
     SxL=0.0_SP
     SxR=0.0_SP
     SyL=0.0_SP
     SyR=0.0_SP
! original variables
     Ubar=0.0_SP
     Vbar=0.0_SP
     Ubar0=0.0_SP
     Vbar0=0.0_SP
     U=0.0_SP
     V=0.0_SP
     HU=0.0_SP
     HV=0.0_SP
     Fx=0.0_SP
     Fy=0.0_SP
     Gx=0.0_SP
     Gy=0.0_SP
     P=0.0_SP
     Q=0.0_SP
     U1p=ZERO
     V1p=ZERO
# if defined (CARTESIAN)
     U4=ZERO
     V4=ZERO
     U1pp=ZERO
     V1pp=ZERO
     U2=ZERO
     V2=ZERO
     U3=ZERO
     V3=ZERO
# endif
     Depth=10.0_SP
     DepthNode=10.0_SP
     H=0.0_SP
     Eta=0.0_SP
     SourceX=0.0_SP
     SourceY=0.0_SP
     PLOT_COUNT=0.0_SP
     PLOT_COUNT_STATION=0.0_SP
     HOTSTART_COUNT=ZERO
     MASK=1
     MASK_STRUC=1
     SCREEN_COUNT=ZERO
     SPONGE=1.0_SP
     SpongeMaker=1.0_SP

! coordinate for cartesian only
! Xco, and Yco

# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Xco(Ibeg) = npx*(Mloc-2*Nghost)*DX
# if defined (CARTESIAN)
     Xco(Ibeg) = (iista-1)*DX
# else
    Xco(Ibeg) = (iista-1)*DX(1,1)
# endif 
! cartesian
!---ykchoi Jan/23/2018]
# else
     Xco(Ibeg) = ZERO
# endif

     DO I = Ibeg+1,Mloc
# if defined (CARTESIAN)
       Xco(I) = Xco(I-1)+DX
# else
       Xco(I) = Xco(I-1)+DX(1,1)
# endif
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
# if defined (CARTESIAN)
       Xco(I) = Xco(I+1)-DX
# else
       Xco(I) = Xco(I+1)-DX(1,1)
# endif
     ENDDO

# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Yco(Jbeg) = npy*(Nloc-2*Nghost)*DY
# if defined (CARTESIAN)
     Yco(Jbeg) = (jjsta-1)*DY
# else
     Yco(Jbeg) = (jjsta-1)*DY(1,1)
#endif
!---ykchoi Jan/23/2018]
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
# if defined (CARTESIAN)
       Yco(J) = Yco(J-1)+DY
# else
       Yco(J) = Yco(J-1)+DY(1,1)
# endif
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
# if defined (CARTESIAN)
       Yco(J) = Yco(J+1)-DY
# else
       Yco(J) = Yco(J+1)-DY(1,1)
# endif
     ENDDO

          
     IF(SHOW_BREAKING)THEN
     AGE_BREAKING = ZERO
     ENDIF

     ROLLER_FLUX = ZERO
     UNDERTOW_U =ZERO
     UNDERTOW_V =ZERO

     PLOT_COUNT=PLOT_INTV
     PLOT_COUNT_STATION=PLOT_INTV_STATION

     SCREEN_COUNT=SCREEN_INTV
   
     DT=ZERO

     IF(Gamma3 > ZERO) THEN
     ELSE
!      gamma3 = 0 means linear shallow water equation
      Gamma1 =ZERO
# if defined (CARTESIAN)
      Gamma2 = ZERO
# endif
     ENDIF

     IF(DISPERSION)THEN
        ! make sure gamma1 and gamma2 are right
     ELSE
       Gamma1=ZERO
# if defined (CARTESIAN)
       Gamma2=ZERO
# endif
     ENDIF

# if defined (CARTESIAN)
! physics - below are fully nonlinear Boussinesq for reference
      a1=Beta_ref*Beta_ref/2.0_SP - 1.0_SP/6.0_SP
      a2=Beta_ref + 1.0_SP/2.0_SP
      b1=Beta_ref*Beta_ref
      b2=Beta_ref
! kennedy equation
      Beta_1=Beta_ref+1.0_SP
      Beta_2=(1.0_SP/5.0_SP)**2/1.0_SP
# else
! for weakly nonlinear equations on sphere (kirby et al 2009)
       a1=ZERO
       a2=ZERO
       b1=1.0_SP/3.0_SP
       b2=-1.0_SP/2.0_SP
# if defined (ZALPHA)
! redefine a1 a2 and b1 and b2 used in trid just for connivence 
! although they can be replaced by beta_1 and beta_2
      a1=Beta_ref*Beta_ref/2.0_SP - 1.0_SP/6.0_SP
      a2=Beta_ref + 1.0_SP/2.0_SP
      b1=Beta_ref*Beta_ref
      b2=Beta_ref
      Beta_1=Beta_ref+1.0_SP
# endif
# endif

! bathymetry

  IF(DEPTH_TYPE(1:3)=='DAT')THEN

! check existing

 INQUIRE(FILE=TRIM(DEPTH_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(DEPTH_FILE), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(DEPTH_FILE), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF  ! exist

# if defined (PARALLEL)
     call GetFile (DEPTH_FILE,Depth)
# else
     OPEN(1,FILE=TRIM(DEPTH_FILE))

       DO J=Jbeg,Jend
        READ(1,*)(Depth(I,J),I=Ibeg,Iend)
       ENDDO

     CLOSE(1)
# endif
  ENDIF

  IF(DEPTH_TYPE(1:3)=='FLA') THEN
    DO J=1,Nloc
     DO I=1,Mloc
      Depth(I,J) = Depth_FLat
     ENDDO
    ENDDO

  ENDIF

# if defined (PARALLEL)

  IF(DEPTH_TYPE(1:3)=='SLO') THEN
    IF(.NOT.ALLOCATED(VarGlob)) ALLOCATE (VarGlob(Mglob,Nglob))
    DO J=1,Nglob
     DO I=1,Mglob
      VarGlob(I,J) = Depth_FLat
     ENDDO
# if defined (CARTESIAN)
     DO I=INT(Xslp/DX)+1,Mglob
      VarGlob(I,J) = Depth_Flat-SLP*(I-(INT(Xslp/DX)+1))*DX
     ENDDO
# else
     DO I=INT(Xslp/DX(1,1))+1,Mglob
      VarGlob(I,J) = Depth_Flat-SLP*(I-(INT(Xslp/DX(1,1))+1))*DX(1,1)
     ENDDO
# endif
    ENDDO ! end of J

    CALL DISTRIBUTE_VarGlob (VarGlob,Depth)

    DEALLOCATE (VarGlob)

  ENDIF

# else

  IF(DEPTH_TYPE(1:3)=='SLO') THEN
# if defined (CARTESIAN)
! dont need this
!    IF(INT(Xslp/DX)>Mloc.OR.INT(Xslp/DX)<1)THEN
# else
!    IF(INT(Xslp/DX(1,1))>Mloc.OR.INT(Xslp/DX(1,1))<1)THEN
# endif
!      WRITE(*,*)'Please reset Xslp, STOP!'
!      STOP
!    ENDIF

    DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      Depth(I,J) = Depth_FLat
     ENDDO

# if defined (CARTESIAN)
     DO I=INT(Xslp/DX)+1+Ibeg,Iend
      Depth(I,J) = Depth_Flat-SLP*(I-(INT(Xslp/DX)+1+Ibeg))*DX
     ENDDO
# else
     DO I=INT(Xslp/DX(1,1))+1+Ibeg,Iend
      Depth(I,J) = Depth_Flat-SLP*(I-(INT(Xslp/DX(1,1))+1+Ibeg))*DX(1,1)
     ENDDO
# endif
    ENDDO
  ENDIF
# endif

! depth at ghost cells and re-construct depth at x,y-interfaces
! use mirror instead of continuous used before, date: 02/24/2019
! the reason for using mirror is eta bc is mirror at open bc
! it is important for complex bathy/topo at bc

    VTYPE=1
# if defined (CARTESIAN)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Depth,VTYPE,PERIODIC)
# else
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Depth,VTYPE)
# endif

! correction bathymetry
    IF (BATHY_CORRECTION)THEN
      CALL CORRECTION 
    ENDIF

! add waterlevel 03/29/2016

      Depth = Depth + WaterLevel
      DEP_WK = DEP_WK + WaterLevel
      IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
        Dep_Ser = Dep_Ser + WaterLevel
      ENDIF

! re-construct Depth

     DO J=1,Nloc
     DO I=2,Mloc
      DepthX(I,J)=0.5_SP*(Depth(I-1,J)+Depth(I,J))
     ENDDO
     ENDDO
     DO J=1,Nloc
      DepthX(1,J)=0.5_SP*(3.0_SP*Depth(1,J)-Depth(2,J))
      DepthX(Mloc1,J)=0.5_SP*(3.0_SP*Depth(Mloc,J)-Depth(Mloc-1,J))
     ENDDO

     DO J=2,Nloc
     DO I=1,Mloc
      DepthY(I,J)=0.5_SP*(Depth(I,J-1)+Depth(I,J))
     ENDDO
     ENDDO
     DO I=1,Mloc
      DepthY(I,1)=0.5_SP*(3.0_SP*Depth(I,1)-Depth(I,2))
      DepthY(I,Nloc1)=0.5_SP*(3.0_SP*Depth(I,Nloc)-Depth(I,Nloc-1))
     ENDDO

!   dont need to re-construct depth in terms of using well-balanced scheme
!   01/21/2012
!   initially we tried to reconstruct depth in order to keep consistency
!   between depth_ele and depth_nod. It turned out it is easy to cause
!   errors if trying to make artifial high walls for dry points 
!   (results in large slopes). We decide to remove the option of reading depth_nod
!   and the reconstruction. fyshi (11/14/2016)

!     DO J=1,Nloc
!     DO I=1,Mloc
!       Depth(I,J)=0.25_SP*(Depthx(I,J)+Depthx(I+1,J)+Depthy(I,J)+Depthy(I,J+1))
!     ENDDO
!     ENDDO


! blowup threshold
      EtaBlowVal=100.0*MAXVAL( abs(Depth(Ibeg:Iend,Jbeg:Jend)) )
# if defined (PARALLEL)
      CALL MPI_ALLREDUCE(EtaBlowVal,myvar_tmp,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
      EtaBlowVal = myvar_tmp
# endif     
 
!friction
     IF(IN_Cd) THEN
# if defined (PARALLEL)
      call GetFile(CD_FILE,Cd)
# else
      OPEN(1,FILE=TRIM(CD_FILE))
       DO J=Jbeg,Jend
        READ(1,*)(Cd(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)
! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        Cd(I,J)=Cd(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        Cd(I,J)=Cd(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        Cd(I,J)=Cd(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        Cd(I,J)=Cd(Iend,J)
       ENDDO
     ENDDO  


# endif
      ELSE
          DO J=1,Nloc
            DO I=1,Mloc
              Cd(I,J) = Cd_fixed
          ENDDO
          ENDDO
      ENDIF



! initial eta u and v for deforming
     IF(INI_UVZ) THEN
        CALL INITIAL_UVZ
       IF(BED_DEFORMATION)THEN
         DO J=1,Nloc
         DO I=1,Mloc
           Depth(I,J)=Depth(I,J)-ETA(I,J)
         ENDDO
         ENDDO         
       ENDIF
     ENDIF 

! initial solitary wave
     IF(WaveMaker(1:7)=='INI_SOL') THEN
       CALL INITIAL_SOLITARY_WAVE(Mloc,Nloc, DX,Xwavemaker,& 
          AMP_SOLI,Dep_Soli,Beta_ref,U,V,Eta,SolitaryPositiveDirection)
     ENDIF



! initial N wave
     IF(WaveMaker(1:6)=='N_WAVE') THEN
       CALL INITIAL_N_WAVE(Mloc,Nloc, DX,x1_Nwave,& 
          x2_Nwave,a0_Nwave,gamma_Nwave,dep_Nwave,U,V,Eta)
     ENDIF



! initial rectangular hump
     IF(WaveMaker(1:7)=='INI_REC') THEN
# if defined (CARTESIAN)
       CALL INITIAL_RECTANGULAR(Mloc,Nloc,Nghost,DX,DY,Xc,Yc,WID,AMP_SOLI, &
                      Eta)
# else
       CALL INITIAL_RECTANGULAR(Mloc,Nloc,Nghost,Dphi,Dtheta,Lon_West,Lat_South, &
             Xc,Yc,WID,AMP_SOLI,Eta)
# endif
     ENDIF


! initial gausian hump
     IF(WaveMaker(1:7)=='INI_GAU') THEN
# if defined (CARTESIAN)
       CALL INITIAL_GAUSIAN(Mloc,Nloc,Nghost,DX,DY,Xc,Yc,AMP_SOLI, WID,&
                      Eta)
# else
       CALL INITIAL_GAUSIAN(Mloc,Nloc,Nghost,Dphi,Dtheta,Lon_West,Lat_South, &
             Xc,Yc,AMP_SOLI,WID,Eta)
# endif
     ENDIF
! initial dipole from x-derivative of gausian hump
     IF(WaveMaker(1:7)=='INI_DIP') THEN
# if defined (CARTESIAN)
       CALL INITIAL_DIPOLE(Mloc,Nloc,Nghost,DX,DY,Xc,Yc,AMP_SOLI, WID,&
                      Eta)
# else
       CALL INITIAL_DIPOLE(Mloc,Nloc,Nghost,Dphi,Dtheta,Lon_West,Lat_South, &
             Xc,Yc,AMP_SOLI,WID,Eta)
# endif
     ENDIF



! initial wave surface
     IF(WaveMaker(1:7)=='INI_OTH') THEN
       CALL INITIAL_WAVE
     ENDIF


     CALL WAVEMAKER_INITIALIZATION
   
   
     IF(DIRECT_SPONGE)THEN
       CALL CALCULATE_SPONGE(Mloc,Nloc,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,SPONGE)
     ENDIF

     IF(DIFFUSION_SPONGE)THEN
       CALL CALCULATE_DIFFUSION_SPONGE(Mloc,Nloc,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,nu_sponge)
     ENDIF

     IF(FRICTION_SPONGE)THEN

       IF(.NOT.ALLOCATED(CD_4_SPONGE)) ALLOCATE(CD_4_SPONGE(Mloc,Nloc))
       CD_4_SPONGE = ZERO

       CALL CALCULATE_FRICTION_SPONGE(Mloc,Nloc,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,CD_4_SPONGE)
     ENDIF


     IF(WaveMaker(1:3)=='ABS')THEN
       CALL CALCULATE_SPONGE_MAKER(Mloc,Nloc,Nghost,DX,DY,&
                            WidthWaveMaker, &
                            R_sponge_wavemaker,A_sponge_wavemaker,SpongeMaker)
     ENDIF

! get Eta and H

    IF(NO_MASK_FILE)THEN
     DO J=1,Nloc
     DO I=1,Mloc
      IF(Eta(I,J)<-DEPTH(I,J))THEN
       MASK(I,J)=0
       Eta(I,J)=MinDepth-Depth(I,J)
# if defined(UseEtaScreen)
       EtaScreen(I,J)=Eta(I,J)
# endif
      ELSE
       MASK(I,J)=1
      ENDIF
     ENDDO
     ENDDO
    ENDIF    

     H=MAX(Eta*Gamma3+Depth,MinDepthFrc)
     HU=H*U
     HV=H*V
    
     IF(DISPERSION)THEN
       CALL CAL_DISPERSION
     ENDIF
    
     Ubar=HU+gamma1*U1p*H
     Vbar=HV+gamma1*V1p*H

  
! read obstacle structures 
     IF(OBSTACLE)THEN

  INQUIRE(FILE=TRIM(OBSTACLE_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(OBSTACLE_FILE), ' specified in input.txt but CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(OBSTACLE_FILE), ' specified in input.txt but CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF
     
# if defined (PARALLEL)
     IF(.NOT.ALLOCATED(VarGlob)) ALLOCATE (VarGlob(Mloc,Nloc)) ! use local here 
     
     call GetFile ( OBSTACLE_FILE, VarGlob )    

     MASK_STRUC = INT(VarGlob)
     DEALLOCATE(VarGlob)
# else
     OPEN(1,FILE=TRIM(OBSTACLE_FILE))
       DO J=Jbeg,Jend
        READ(1,*)(MASK_STRUC(I,J),I=Ibeg,Iend)
       ENDDO
     CLOSE(1)
# endif
     ENDIF

! read breakwater 
     IF(BREAKWATER)THEN

  INQUIRE(FILE=TRIM(BREAKWATER_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(BREAKWATER_FILE), ' specified in input.txt but CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(BREAKWATER_FILE), ' specified in input.txt but CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! we have to use global for width 
     ALLOCATE(BreakWaterWidth(Mglob+2*Nghost,Nglob+2*Nghost))
     BreakWaterWidth = 0.0_SP     

# if defined (PARALLEL)
     if (myid.eq.0) then
        OPEN(1,FILE=TRIM(BREAKWATER_FILE))
        DO J=Nghost+1,NGlob+NGhost
           READ(1,*)(BreakWaterWidth(I,J),I=Nghost+1,MGlob+Nghost)
        ENDDO
        CLOSE(1)
! ghost cells
        DO I=Nghost+1,MGlob+Nghost
           DO J=1,Nghost
              BreakWaterWidth(I,J)=BreakWaterWidth(I,Nghost+1)
           ENDDO
           DO J=NGlob+Nghost+1,NGlob+2*Nghost
              BreakWaterWidth(I,J)=BreakWaterWidth(I,NGlob+Nghost)
           ENDDO
        ENDDO
        DO J=1,NGlob+2*Nghost
           DO I=1,Nghost
              BreakWaterWidth(I,J)=BreakWaterWidth(Nghost+1,J)
           ENDDO
           DO I=MGlob+Nghost+1,MGlob+2*Nghost
              BreakWaterWidth(I,J)=BreakWaterWidth(MGlob+Nghost,J)
           ENDDO
        ENDDO
     endif

# else
     OPEN(1,FILE=TRIM(BREAKWATER_FILE))

       DO J=Jbeg,Jend
        READ(1,*)(BreakWaterWidth(I,J),I=Ibeg,Iend)
       ENDDO

     CLOSE(1)


! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        BreakWaterWidth(I,J)=BreakWaterWidth(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        BreakWaterWidth(I,J)=BreakWaterWidth(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        BreakWaterWidth(I,J)=BreakWaterWidth(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        BreakWaterWidth(I,J)=BreakWaterWidth(Iend,J)
       ENDDO
     ENDDO  

# endif
     CALL CALCULATE_CD_BREAKWATER

     DEALLOCATE(BreakWaterWidth)

     ENDIF  ! endif breakwater
! end breakwater

     DO J=1,Nloc
     DO I=1,Mloc
      IF(MASK_STRUC(I,J)==0)THEN
        Depth(I,J)=-LARGE
      ENDIF
     ENDDO
     ENDDO

     MASK=MASK*MASK_STRUC

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      MASK9(I,J)=MASK(I,J)*MASK(I-1,J)*MASK(I+1,J)  &
                *MASK(I+1,J+1)*MASK(I,J+1)*MASK(I-1,J+1) &
                *MASK(I+1,J-1)*MASK(I,J-1)*MASK(I-1,J-1) 
     ENDDO
     ENDDO

# if defined(PARALLEL)
     CALL PHI_INT_EXCH(MASK)
     CALL PHI_INT_EXCH(MASK9)
# endif



# if defined (COUPLING)

     CALL COUPLING_INITIAL

# endif
 

# if defined(CARTESIAN)
# else
    ! use slope control in caution, should set false unless for a spherical
    ! ocean basin domain and a huge slope at coastline. 
    SLOPE_CTR = .FALSE.
    MAX_SLOPE = 0.1
    IF(SLOPE_CTR)THEN
      WRITE(*,*)'SLOPE CONTROL IS ON'
    ENDIF

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
      SlopeX(I,J) = (Depthx(I+1,J)-Depthx(I,J))/DX(I,J)
      SlopeY(I,J) = (Depthy(I,J+1)-Depthy(I,J))/DY(I,J)
      IF(SLOPE_CTR)THEN
        IF(ABS(SlopeX(I,J))>MAX_SLOPE)SlopeX(I,J)=SIGN(MAX_SLOPE,SlopeX(I,J))
        IF(ABS(SlopeY(I,J))>MAX_SLOPE)SlopeY(I,J)=SIGN(MAX_SLOPE,SlopeY(I,J))
      ENDIF
    ENDDO
    ENDDO

# endif

! deal with masks, this is great for an extremely large bed slope
! at edges of mask point, depth at cell interface can cause unreasonable large 
! depth gradient, even if the slope cap is on. Making depth locally flat can avoid
! this happening. This scheme is also used in update_mask subroutine
! 01/21/2012  
! HOWEVER, this truncation can affect the model accuracy as pointed by 
! Choi (private communication, in Thacker bowl test case, 07/03/2016). 
! In the following, I keep an option to use more accurate solution. 
! In Makefile, define -DIGNORE_BIG_SLOPE to ignore big slopes
# if defined (IGNORE_BIG_SLOPE)
    ! do not truncate depthx and depthy
# else
      DO J=2,Nloc-1
      DO I=2,Mloc-1
        IF(MASK(I,J)<1)THEN
         DepthX(I,J)=Depth(I-1,J)
         DepthX(I+1,J)=Depth(I+1,J)
         DepthY(I,J)=Depth(I,J-1)
         DepthY(I,J+1)=Depth(I,J+1)
        ENDIF
      ENDDO
      ENDDO
# endif

      IF(SHOW_BREAKING)THEN
       T_brk=20.0_SP
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Oct, 2021
!Salatin, R., Chen, Q., Bak, A. S., Shi, F., & Brandt, S. R. (2021). Effects of
!wave coherence on longshore variability of nearshore wave processes. Journal
!of Geophysical Research: Oceans,126, e2021JC017641.
!https://doi.org/10.1029/2021JC017641

    ! This is used inside source.F to speed up the calculations
# if defined (CARTESIAN)
     DXg=DX
     DYg=DY
# else
     DXg=DX(1,1)
     DYg=DY(1,1)
# endif
    ilo = Iend+1
    ihi = Ibeg-1
    jlo = Jend+1
    jhi = Jbeg-1
    DO J=Jbeg,Jend
# if defined (PARALLEL)
        ymk_wk(J)=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
        ymk_wk(J)=(J-Jbeg)*DYg
# endif
        IF(ABS(ymk_wk(J)-Yc_WK)<Ywidth_WK/2.0_SP) THEN
            IF(J<jlo) jlo = J
            IF(J>jhi) jhi = J
        ENDIF
    ENDDO
    DO I=Ibeg,Iend
# if defined (PARALLEL)
        xmk_wk(I)=(I-Ibeg)*DXg + (iista-1)*DXg
# else
        xmk_wk(I)=(I-Ibeg)*DXg
# endif
        IF(ABS(xmk_wk(I)-Xc_WK)<Width_WK) THEN
            IF(I<ilo) ilo = I
            IF(I>ihi) ihi = I
        ENDIF
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE INITIALIZATION


!-------------------------------------------------------------------------------------
!
!    INITIAL_UVZ is subroutine of given initial u v and eta 
!
!    HISTORY:
!    02/03/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIAL_UVZ
      USE GLOBAL
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mloc,Nloc) :: rMASK
      LOGICAL :: FILE_EXIST

IF(.NOT.NO_UV_FILE)THEN

  INQUIRE(FILE=TRIM(U_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0) WRITE(*,*) TRIM(U_FILE), ' specified in input.txt but does not exist. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(U_FILE), ' specified in input.txt but does not exist. STOP'
    STOP
# endif
  ENDIF


  INQUIRE(FILE=TRIM(V_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0) WRITE(*,*) TRIM(V_FILE), ' specified in input.txt but does not exist. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(V_FILE), ' specified in input.txt but does not exist. STOP'
    STOP
# endif
  ENDIF

ENDIF

  INQUIRE(FILE=TRIM(ETA_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0) WRITE(*,*) TRIM(ETA_FILE), ' specified in input.txt but does not exist. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(ETA_FILE), ' specified in input.txt but does not exist. STOP'
    STOP
# endif
  ENDIF

# if defined (PARALLEL)
IF(.NOT.NO_UV_FILE)THEN
      CALL GetFile(U_FILE,U)
      CALL GetFile(V_FILE,V)
ELSE
      U=ZERO
      V=ZERO
ENDIF

      CALL GetFile(ETA_FILE,ETA)

      IF(.NOT.NO_MASK_FILE)THEN
      CALL GetFile(MASK_FILE,rMASK)
      ENDIF

# else

IF(.NOT.NO_UV_FILE)THEN
      OPEN(1,FILE=TRIM(U_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(U(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(V_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(V(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)
ELSE
      U=ZERO
      V=ZERO
ENDIF

      OPEN(1,FILE=TRIM(ETA_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(ETA(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)

     IF(.NOT.NO_MASK_FILE)THEN
      OPEN(1,FILE=TRIM(MASK_FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(rMASK(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)
     ENDIF

# endif

     IF(.NOT.NO_MASK_FILE)THEN
      MASK=INT(rMASK)
     ELSE
      MASK = 1
     ENDIF

END SUBROUTINE INITIAL_UVZ




