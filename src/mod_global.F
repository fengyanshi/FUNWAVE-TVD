!------------------------------------------------------------------------------------
!
!      FILE mod_global.F
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
!    GLOBAL is the module to define all global variables, parameters
!    
!    HISTORY: 05/01/2010 Fengyan Shi
!             the module is updated corresponding to modifications in subroutines
!
! --------------------------------------------------
MODULE GLOBAL
       USE PARAM
       IMPLICIT NONE
       SAVE

# if defined (PARALLEL)
! MPI variables
       INTEGER :: myid,ier
       INTEGER :: comm2d
       INTEGER :: n_west, n_east, n_suth, n_nrth
       INTEGER :: npx,npy
       INTEGER :: ndims=2
       INTEGER :: NumberProcessor
       INTEGER,DIMENSION(:,:),ALLOCATABLE :: ProcessorID 
       INTEGER, DIMENSION(2) :: dims, coords
       LOGICAL, DIMENSION(2) :: periods
       LOGICAL :: reorder = .true.
       INTEGER :: nprocs                 !ykchoi(04/May/2017)
       INTEGER :: iista, iiend, jjsta, jjend !ykchoi(04/May/2017)
# endif
       INTEGER :: px,py

! station data
       INTEGER :: NumberStations
       INTEGER,DIMENSION(:),ALLOCATABLE :: ista,jsta,nsta
       REAL(SP):: PLOT_INTV_STATION,PLOT_COUNT_STATION
       INTEGER :: StationOutputBuffer

! timer
       REAL(8) :: tbegin,tend
       INTEGER :: ISTAGE

# if defined(FILTERING)
       INTEGER :: filter_count=0, filter_interval = 1
# endif

! define parameters

       ! [mayhl] Added global input file name for modules
       CHARACTER(LEN=80) INPUT_FILE_NAME
       CHARACTER(LEN=80) TITLE
       CHARACTER(LEN=80) ReadType
       CHARACTER(LEN=80) DEPTH_TYPE
       CHARACTER(LEN=80) DEPTH_FILE
       CHARACTER(LEN=80) ETA_FILE
       CHARACTER(LEN=80) U_FILE
       CHARACTER(LEN=80) V_FILE
       CHARACTER(LEN=80) MASK_FILE
       CHARACTER(LEN=80) DX_FILE
       CHARACTER(LEN=80) DY_FILE
       CHARACTER(LEN=80) Coriolis_FILE
       CHARACTER(LEN=80) OBSTACLE_FILE
       CHARACTER(LEN=80) BREAKWATER_FILE
       CHARACTER(LEN=80) RESULT_FOLDER
       CHARACTER(LEN=80) STATIONS_FILE
       CHARACTER(LEN=80) WaveMaker
       CHARACTER(LEN=80) Time_Scheme
       CHARACTER(LEN=80) CONSTR
       CHARACTER(LEN=80) HIGH_ORDER
       CHARACTER(LEN=80) FIELD_IO_TYPE !mayhl 17/06/12
       ![ykchoi(14.12.24)
       REAL(SP) :: HotStartTime
       !ykchoi(14.12.24)]

       LOGICAL :: NO_MASK_FILE = .TRUE.
       LOGICAL :: NO_UV_FILE = .TRUE.

        CHARACTER(LEN=16)::FORMAT_LEN=' '


       REAL(SP),PARAMETER,DIMENSION(3)::alpha=(/0.0_SP,3.0_SP/4.0_SP,1.0_SP/3.0_SP/)
        REAL(SP),PARAMETER,DIMENSION(3)::beta=(/1.0_SP,1.0_SP/4.0_SP,2.0_SP/3.0_SP/)
       REAL(SP)::Kappa   ! parameter for order of muscl

! some global variables
! Mloc1=Mloc+1, Nloc1=Nloc+1
       INTEGER :: Mglob,Nglob,Mloc,Nloc,Mloc1,Nloc1
       INTEGER, PARAMETER :: Nghost = 3
       INTEGER :: Ibeg,Iend,Jbeg,Jend,Iend1,Jend1
# if defined (CARTESIAN)
       REAL(SP):: DX,DY
# else
       REAL(SP) :: Lon_West,Lat_South,Dphi,Dtheta
# endif
       REAL(SP)::  DT,TIME,TOTAL_TIME,PLOT_INTV,PLOT_COUNT,&
                  SCREEN_INTV,SCREEN_COUNT, &
                   DT_fixed,PLOT_START_TIME
       REAL(SP) :: HOTSTART_INTV,HOTSTART_COUNT
       INTEGER :: icount=-1  ! for output file number
       INTEGER :: icount_mean=0  ! for output mean file number
       INTEGER :: icount_hotstart=0
       INTEGER :: FileNumber_HOTSTART
       
       REAL(SP) :: MinDepth=0.001
       REAL(SP) :: MinDepthFrc=0.5
       REAL(SP) :: ArrTimeMin=0.001
       REAL(SP) :: CFL=0.15_SP
       REAL(SP) :: FroudeCap=10.0_SP
       REAL(SP) :: SWE_ETA_DEP=0.7_SP

       LOGICAL :: BATHY_CORRECTION = .FALSE.
       LOGICAL :: BREAKWATER = .FALSE.

       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
                BreakWaterWidth,CD_breakwater
       REAL(SP) :: BreakWaterAbsorbCoef

! coordinates for Cartesian only
      REAL(SP),DIMENSION(:),ALLOCATABLE ::  Xco,Yco

! switch for dispersion
!   gamma1 is for extra M term
       REAL(SP) :: gamma1
# if defined (CARTESIAN)
       REAL(SP) :: gamma2
# endif
! gamma3 is for linear shallow water equation
       REAL(SP) :: gamma3
       REAL(SP) :: Beta_ref=-0.531_SP
! kennedys equations
       REAL(SP) :: Beta_1,Beta_2
! a1=beta_ref^2/2-1/6, a2=beta_ref+1/2, b1=beta_ref^2, b2=beta_ref
       REAL(SP) :: a1,a2,b1,b2
       LOGICAL :: DISPERSION=.FALSE.
       LOGICAL :: DISP_TIME_LEFT=.FALSE. 
       LOGICAL :: StretchGrid = .FALSE.
                 ! put time derivative dispersion term on left
       LOGICAL :: FIXED_DT = .FALSE.

! some local variables
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
# if defined (CARTESIAN)
          U4xL,U4xR,V4yL,V4yR, &
          V4xL,V4xR,U4yL,U4yR, &   !ykchoi (15.08.06.)
          U4,V4,U1p,V1p, &
          U1pp,V1pp, &
          U2,V2,U3,V3, &
# else
# if defined(ZALPHA)
          U4xL,U4xR,V4yL,V4yR, &
          V4xL,V4xR,U4yL,U4yR, &   !ykchoi (15.08.06.)
          U4,V4,&
# endif
          Dx,Dy, Coriolis,Lat_theta, &
          U1p,V1p, SlopeX,SlopeY,&
# endif
          DelxU,DelxHU,DelxV,DelxEtar,&
          DelxHV, DelyHU, &
          DelyU,DelyHV,DelyV,DelyEtar,&
          UxL,UxR,VxL,VxR,&
          HUxL,HUxR,HUyL,HUyR,HxL,HxR, &
          EtaRxL,EtaRxR,&
          UyL,UyR,VyL,VyR,&
          HVxL,HVxR,HVyL,HVyR,HyL,HyR, &
          EtaRyL,EtaRyR, &
          PL,PR,QL,QR, &
          FxL,FxR,FyL,FyR, &
          GxL,GxR,GyL,GyR, &
          SxL,SxR,SyL,SyR, &
! cross-derivatives 
          Vxy,DVxy,Uxy,DUxy, &
! second-derivatives
          Uxx,DUxx,Vyy,DVyy, &
! first-derivatives
          Ux,Vx,Uy,Vy,DUx,DUy,DVx,DVy, &
          ETAx,ETAy, ETAT, ETATx,ETATy, &
! time-derivative
          U0,V0,Ut,Vt,Utx,Vty,Utxx,Utxy,Vtxy,Vtyy,&
          DUtxx,DUtxy,DVtxy,DVtyy,DUtx,DVty,&
! original variables
          Fx,Fy,U,V,HU,HV,&
          Gx,Gy,P,Q,SourceX,SourceY,Int2Flo, &
          tmp4preview,HeightMax,HeightMin,VelocityMax,&
          MomentumFluxMax,VorticityMax,ARRTIME
# if defined (ITERATION)
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Ui,Vi,ETAi,&
                        UbarOld,VbarOld,EtaOld
       INTEGER :: Ki
       LOGICAL :: ACCURATE = .FALSE.
# endif
! 
! wetting and drying
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: MASK,MASK_STRUC,MASK9
        REAL(SP) :: Dmass,WetArea,DwetEta
# if defined(UseEtaScreen)
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: EtaScreen
# endif
! wave maker
        LOGICAL :: EqualEnergy = .FALSE.
        LOGICAL :: WaveMakerCurrentBalance = .FALSE.
        REAL(SP) :: WaveMakerCd
        REAL(SP)::AMP_SOLI,DEP_SOLI,LAG_SOLI, CPH_SOLI,XWAVEMAKER, &
                  Xc,Yc, WID,Xc_WK,Tperiod,AMP_WK,DEP_WK,Theta_WK, &
                  rlamda,Time_ramp,D_gen,Beta_gen,Width_WK,Delta_WK,&
                  Ywidth_WK,Yc_WK
        LOGICAL :: SolitaryPositiveDirection = .TRUE.
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: D_gen_ir,rlamda_ir,phase_ir
        REAL(SP),DIMENSION(:),ALLOCATABLE :: Beta_gen_ir,omgn_ir
        REAL(SP),DIMENSION(:),ALLOCATABLE :: omgn2D
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Wavemaker_Mass
        REAL(SP) :: FreqMin,FreqMax,FreqPeak,GammaTMA,Hmo,ThetaPeak,&
                    Sigma_Theta
        REAL(SP),DIMENSION(:,:,:),ALLOCATABLE ::Cm,Sm
        INTEGER :: Nfreq,Ntheta 
        REAL(SP)::x1_Nwave = 5.0, &
                  x2_Nwave = 5.0, &
                  a0_Nwave = 1.0, &
                  gamma_Nwave = -3.0, &
                  dep_Nwave = 1.0
!      for measure time series
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: WAVE_COMP
       REAL(SP),DIMENSION(:),ALLOCATABLE :: Beta_genS,D_genS
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Beta_gen2D,D_gen2D,rlamda2D, &
                                              Phase2D
       REAL(SP),DIMENSION(:),ALLOCATABLE :: Freq,Dire
       REAL(SP) :: PeakPeriod
       INTEGER :: NumWaveComp,NumFreq,NumDir
       CHARACTER(LEN=80) WaveCompFile
! friction
        LOGICAL :: IN_Cd=.FALSE.
        REAL(SP):: Cd_fixed
        CHARACTER(LEN=80) CD_FILE
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Cd

        LOGICAL :: ROLLER=.FALSE.
        REAL(SP) :: ROLLER_SWITCH 

! sponge
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: SPONGE,SpongeMaker
        REAL(SP)::Sponge_west_width,Sponge_east_width, &
                  Sponge_south_width,Sponge_north_width, &
                  R_sponge,A_sponge


! smagorinsky and wave height
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Umean,Vmean,&
                  ETAmean,Usum,Vsum,ETAsum, nu_smg
      REAL(SP)::T_INTV_mean = 20.0,T_sum=0.0,C_smg=0.25
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
                WaveHeightRMS,WaveHeightAve,Emax,Emin,& 
                HrmsSum,HavgSum
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: Num_Zero_Up

      !ykchoi
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: ETA2sum,ETA2mean,SigWaveHeight

! depth H=Eta+Depth, 
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Depth,H,&
            DepthNode,Depthx,Depthy
       REAL(SP)::Depth_Flat, SLP,Xslp

! updating variables
       REAL(SP),DIMENSION(:,:),ALLOCATABLE::Ubar0,Vbar0,Eta0,&
                               Ubar,Vbar,Eta

! water level, fyshi 03/29/2016
       REAL(SP) :: WaterLevel = 0.0

# if defined(COUPLING)
! coupling
       CHARACTER(LEN=80) COUPLING_FILE
       INTEGER :: icount_coupling = 1
       INTEGER :: N_COUPLING_EAST,&
                  N_COUPLING_WEST,&
                  N_COUPLING_SOUTH,&
                  N_COUPLING_NORTH,&
                  N_COUPLING_DATA,&
                  J_START_EAST,&
                  J_START_WEST,&
                  I_START_SOUTH,&
                  I_START_NORTH
       REAL(SP),DIMENSION(:,:),ALLOCATABLE::U_COUPLING_EAST,&
                  V_COUPLING_EAST,Z_COUPLING_EAST,&
                  U_COUPLING_WEST,V_COUPLING_WEST,Z_COUPLING_WEST,&
                  U_COUPLING_SOUTH,V_COUPLING_SOUTH,Z_COUPLING_SOUTH,&
                  U_COUPLING_NORTH,V_COUPLING_NORTH,Z_COUPLING_NORTH
       REAL(SP) :: TIME_COUPLING_1,TIME_COUPLING_2
       INTEGER :: Kstart_EAST,Kend_EAST,Kshift_EAST,&
                  Kstart_WEST,Kend_WEST,Kshift_WEST,&
                  Kstart_SOUTH,Kend_SOUTH,Kshift_SOUTH,&
                  Kstart_NORTH,Kend_NORTH,Kshift_NORTH
       LOGICAL :: IN_DOMAIN_EAST,IN_DOMAIN_WEST,&
                  IN_DOMAIN_SOUTH,IN_DOMAIN_NORTH
# endif

! output logical parameters
       INTEGER :: OUTPUT_RES = 1
       LOGICAL :: OUT_U=.FALSE., OUT_V=.FALSE.,OUT_ETA=.FALSE., &
# if defined(UseEtaScreen)
                  OUT_ETASCREEN=.FALSE., &
# endif
                  OUT_MASK=.FALSE.,OUT_SXL=.FALSE.,OUT_SXR=.FALSE.,&
                  OUT_SYL=.FALSE., OUT_SYR=.FALSE.,&
                  OUT_SourceX=.FALSE., OUT_SourceY=.FALSE., &
                  OUT_P=.FALSE., OUT_Q=.FALSE., &
                  OUT_Fx=.FALSE., OUT_Fy=.FALSE.,&
                  OUT_Gx=.FALSE., OUT_Gy=.FALSE.,&
                  OUT_MASK9=.FALSE., OUT_DEPTH=.FALSE., &
                  OUT_TMP=.FALSE., OUT_AGE=.FALSE.,OUT_NU=.FALSE., &
                  OUT_Hmax=.FALSE., &
                  OUT_Hmin=.FALSE., &
                  OUT_Umax=.FALSE.
       LOGICAL :: OUT_MFmax=.FALSE., &
                  OUT_VORmax=.FALSE.,OUT_Time, &
                  OUT_ROLLER = .FALSE., &
                  OUT_UNDERTOW = .FALSE.

       LOGICAL :: OUT_Umean=.FALSE.,OUT_Vmean=.FALSE.,&
                  OUT_ETAmean=.FALSE.,&
                  OUT_WaveHeight = .FALSE.,&
                  OUT_WaveHeightSig = .FALSE.    !ykchoi



!# if defined (CARTESIAN)
! periodic boundary conditions
       LOGICAL :: PERIODIC=.FALSE.
!# endif
! dispersion control
       LOGICAL :: OBSTACLE=.FALSE., HOT_START=.FALSE.
! sponge
       LOGICAL :: SPONGE_ON=.FALSE.
! breaking
       LOGICAL :: SHOW_BREAKING=.TRUE.

! slope control, use it in caution, should set false unless for a spherical
! ocean basin domain and slope > 1:5
       LOGICAL :: SLOPE_CTR = .FALSE.
       REAL(SP) :: MAX_SLOPE

! absorbing-generating wavemaker
       REAL(SP) :: Dep_Ser,WidthWaveMaker,&
                   R_sponge_wavemaker,A_sponge_wavemaker

       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Amp_Ser,Phase_LEFT
       REAL(SP),DIMENSION(:),ALLOCATABLE :: Per_Ser,Phase_Ser,&
                                            Theta_Ser
       REAL(SP),DIMENSION(:),ALLOCATABLE :: Segma_Ser,Wave_Number_Ser
       REAL(SP) :: Stokes_Drift_Ser
       REAL(SP),DIMENSION(:,:,:),ALLOCATABLE ::Cm_eta,Sm_eta, &
                Cm_u,Sm_u,Cm_v,Sm_v
!       LOGICAL :: WAVE_DATA = .FALSE.
       CHARACTER(LEN=80) WAVE_DATA_TYPE

! eddy viscosity breaking
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: AGE_BREAKING
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: ROLLER_FLUX,UNDERTOW_U,UNDERTOW_V
        REAL(SP) :: Cbrk1=0.65,Cbrk2=0.35,T_brk
        ! use T_brk to judge breakers
        LOGICAL :: INI_UVZ=.FALSE.
        LOGICAL :: BED_DEFORMATION = .FALSE.

       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: nu_break,nu_sponge
       REAL(SP) :: nu_bkg
       LOGICAL :: DIFFUSION_SPONGE = .FALSE.
       LOGICAL :: DIRECT_SPONGE = .FALSE.
       LOGICAL :: FRICTION_SPONGE = .FALSE.
       LOGICAL :: VISCOSITY_BREAKING = .FALSE.
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: CD_4_SPONGE
       REAL(SP) :: Csp = 0.15
       REAL(SP) :: CDsponge = 0.0
       REAL(SP) :: WAVEMAKER_Cbrk
       LOGICAL :: ETA_LIMITER = .FALSE.
       REAL(SP) :: CrestLimit, TroughLimit
       !ykchoi 
       REAL(SP) :: EtaBlowVal
       REAL(SP) :: STEADY_TIME
       !ykchoi(08.18.2015) :: for viscosity of wavemaker
       REAL(SP) :: visbrk, WAVEMAKER_visbrk
       LOGICAL :: WAVEMAKER_VIS=.FALSE.
!      REAL(SP) :: PLOT_SMALLINTV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Oct, 2021
!Salatin, R., Chen, Q., Bak, A. S., Shi, F., & Brandt, S. R. (2021). Effects of
!wave coherence on longshore variability of nearshore wave processes. Journal
!of Geophysical Research: Oceans,126, e2021JC017641.
!https://doi.org/10.1029/2021JC017641

! introduce some global variables

    INTEGER,DIMENSION(:),ALLOCATABLE :: loop_index
    INTEGER :: jlo,jhi,ilo,ihi
    REAL(SP), DIMENSION(:),ALLOCATABLE :: xmk_wk, ymk_wk
    REAL(SP)::alpha_c ! coherence percentage for WK_NEW_IRR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE GLOBAL
