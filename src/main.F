!------------------------------------------------------------------------------------
!
!      FILE main.F
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
Program FUNWAVE_TVD
!-------------------------------------------------------------------------------------
!    VERSION 3.0
!
!    MAIN - READ_INPUT
!         - INDEX
!         - ALLOCATE_VARIABLES
!         - INITIALIZATION 
!         DO LOOP
!             - VARIABLE_UPDATING
!             - EXCHANGE
!             - ESTIMATE_DT
!             - RUNGE-KUTTA
!             - DISPERSION
!             - FLUXES          
!             - SourceTerms      
!             - ESTIMATE_HUV
!             - EXCHANGE 
!             - FLUXES                 
!             - SourceTerms again
!
!             - STATISTICS
!             - PREVIEW
!          ENDDO LOOP
!-------------------------------------------------------------------------------------
! ** OPEN FILES **
!  (1): read input, (2): output, (3): log, (4): !write/read hotstart
!-------------------------------------------------------------------------------------
! ** HOT START DATA **
!   NOTE: read input.txt first, if HOT_START, then read  
        ! -- dimension
! Mloc,Nloc,Mloc1,Nloc1
! Nghost
! Ibeg,Iend,Jbeg,Jend,Iend1,Jend1
!   NOTE: need to confirm if the saved data is consistent with input.txt
        ! -- time
! TIME
! TOTAL_TIME
! PLOT_INTV
! PLOT_COUNT
! SCREEN_INTV
! SCREEN_COUNT
! HOTSTART_INTV
! ICOUNT
        ! spacing
! DX,DY
        ! -- physics
! DISPERSION
! Gamma1
! a1,a2,b1,b2
! SWE_ETA_DEP
        ! -- numerics
! Time_Scheme
! HIGH_ORDER
! CONSTR
! CFL
! FroudeCap
! DISP_TIME_LEFT
        ! -- wet-dry
! MinDepth,MinDepthfrc

        ! -- depth
! DEPTH
! DEPTHx
! DEPTHy
        ! variables
! U
! V
! if (.NOT.DISP_TIME_LEFT)THEN
! U0
! V0
! endif
! Ubar
! Vbar
! ETA 
! H
! MASK
! MASK9
! MAST_STRUC
!
       ! -- wavemaker
! if (WAVEMAKER is WK_IRR)
! turns out the data for Cm Sm too large, calculate it when hotstart
!
! if (WAVEMAKER is WK_REG)
! D_gen
! Beta_gen
! rlamda
! 
!
!-------------------------------------------------------------------------------------
     USE GLOBAL
# if defined (VESSEL)
     USE VESSEL_MODULE
# endif
# if defined (SEDIMENT)
     USE SEDIMENT_MODULE
# endif

# if defined (METEO)
     USE METEO_MODULE
# endif
# if defined (TRACKING)
     USE TRACER
# endif

     USE TIDE_MODULE

     IMPLICIT NONE

!     INTEGER::ISTAGE ! moved to mod_global 09/12/2017
!     REAL(SP) :: tbegin,tend  ! moved to mod_global 07/29/2016

# if defined (PARALLEL)
     CALL MPI_INIT ( ier )
# endif

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'MPI_INIT'     
#    endif
# endif

     CALL READ_INPUT

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'READ_INPUT'     
#    else   
       WRITE(*,*) 'READ_INPUT'
#    endif
# endif

      ![ykchoi(14.12.24.)
      IF(INI_UVZ)THEN
        TIME=HotStartTime
      ENDIF
      !ykchoi(14.12.24.)]

     CALL INDEX

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'INDEX'     
#    else   
       WRITE(*,*) 'INDEX'
#    endif
# endif

! allocate variables
     CALL ALLOCATE_VARIABLES

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'ALLOCATION'     
#    else   
       WRITE(*,*) 'ALLOCATION'
#    endif
# endif

     CALL INITIALIZATION

     CALL TIDE_INITIAL
      
# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'INITIALIZATION'     
#    else   
       WRITE(*,*) 'INITIALIZATION'
#    endif
# endif

# if defined (COUPLING)
     TIME=TIME_COUPLING_1
# endif

# if defined (VESSEL)
     CALL VESSEL_INITIAL
# endif

# if defined (METEO)
     CALL METEO_INITIAL
# endif

# if defined (SEDIMENT)
     CALL SEDIMENT_INITIAL
# endif

#if defined (TRACKING)
     CALL TRACER_INITIAL
# endif

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'INITIAL VESSEL, METEO, OR SEDIMENT if included'     
#    else   
       WRITE(*,*) 'INITIAL VESSEL, METEO, OR SEDIMENT if included'
#    endif
# endif

! time integration

     ! record wall time
# if defined (PARALLEL)
     if(myid == 0) tbegin = MPI_Wtime( )
# else
     CALL CPU_TIME(tbegin)
# endif


   DO WHILE (TIME<TOTAL_TIME)

!     move output here to get the initial condition 11/27/2018

      CALL OUTPUT

     IF(WaveMaker(1:7)=='LEF_SOL')THEN
       CALL SOLITARY_WAVE_LEFT_BOUNDARY
     ENDIF   

! update three variables
     Eta0=Eta
     Ubar0=Ubar
     Vbar0=Vbar  
# if defined (SEDIMENT)
     CHH0=CHH
# endif

! previous update_mask
   

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'UPDATE Eta Ubar Vbar, MASK'     
#    else   
       WRITE(*,*) 'UPDATE Eta Ubar Vbar, MASK'
#    endif
# endif

     CALL EXCHANGE

# if defined (COUPLING)
         CALL OneWayCoupling
# endif

 
     CALL ESTIMATE_DT(Mloc,Nloc,DX,DY,U,V,H,MinDepthFrc,DT,CFL,TIME)

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'ESTIMATE_DT'     
#    else   
       WRITE(*,*) 'ESTIMATE_DT'
#    endif
# endif
  
! U0, V0 are moved to following part due to computation of Ut, Vt.
      U0=U   !ykchoi(15. 08. 06.)
      V0=V   !ykchoi

# if defined (VESSEL)
       CALL VESSEL_FORCING
# endif

# if defined (METEO)
       CALL METEO_FORCING
# endif

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'vessel, meteo forcing, if included'     
#    else   
       WRITE(*,*) 'vessel, meteo forcing, if included'
#    endif
# endif

     ! 3-ORDER RUNGE-KUTTA TIME STEPPING
     DO ISTAGE=1,3

       IF(DISPERSION)THEN
         CALL Cal_Dispersion
       ENDIF 

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'Cal_Dispersion'     
#    else   
       WRITE(*,*) 'Cal_Dispersion'
#    endif
# endif

       CALL FLUXES

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'FLUXES'     
#    else   
       WRITE(*,*) 'FLUXES'
#    endif
# endif

       CALL SourceTerms   ! put sourceterms after fluxes in order to get eta_t

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'SourceTerms'     
#    else   
       WRITE(*,*) 'SourceTerms'
#    endif
# endif

       CALL ESTIMATE_HUV(ISTAGE) 

! etascreen was added, update_mask was moved here from outside RK 
       CALL UPDATE_MASK

       IF(TIDAL_BC)THEN
         CALL TIDE_BC
       ENDIF

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'ESTIMATE_HUV'     
#    else   
       WRITE(*,*) 'ESTIMATE_HUV'
#    endif
# endif   

# if defined (SEDIMENT)
     CALL SEDIMENT_ADVECTION_DIFFUSION(ISTAGE) 
# endif

       CALL WAVE_BREAKING

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'WAVE_BREAKING'     
#    else   
       WRITE(*,*) 'WAVE_BREAKING'
#    endif
# endif
       
       CALL EXCHANGE

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'EXCHANGE'     
#    else   
       WRITE(*,*) 'EXCHANGE'
#    endif
# endif
       
# if defined (COUPLING)
       CALL OneWayCoupling
# endif
         

       IF(WaveMaker(1:3)=='ABS') THEN
         CALL ABSORBING_GENERATING_BC

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'ABS wavemaker'     
#    else   
       WRITE(*,*) 'ABS wavemaker'
#    endif
# endif

       ENDIF

       IF(WaveMaker(1:11)=='LEFT_BC_IRR') THEN
         CALL IRREGULAR_LEFT_BC
       ENDIF

       IF(DIRECT_SPONGE)THEN
           CALL SPONGE_DAMPING
       ENDIF

     ENDDO

# if defined (SEDIMENT)
!   moved out of the RK loop 07/02/2019
     IF(Bed_Change)THEN
       CALL MORPHOLOGICAL_CHANGE
     ENDIF
# endif

     CALL MIXING_STUFF

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'MIXING_STUFF'     
#    else   
       WRITE(*,*) 'MIXING_STUFF'
#    endif
# endif

# if defined (TRACKING)
     CALL TRACK_XY
# endif


!  find maximum eta velocity 

      IF (OUT_Hmax.OR.OUT_Hmin.OR.OUT_Umax.OR.OUT_MFmax.OR.OUT_VORmax.OR.OUT_Time)THEN
        CALL MAX_MIN_PROPERTY
      ENDIF        

      CALL CHECK_BLOWUP

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'CHECK_BLOWUP'     
#    else   
       WRITE(*,*) 'CHECK_BLOWUP'
#    endif
# endif

# if defined (DEBUG)
#    if defined (PARALLEL)
       IF(myid==0) WRITE(*,*) 'CHECK_BLOWUP'     
#    else   
       WRITE(*,*) 'CHECK_BLOWUP'
#    endif
# endif
  
   END DO

   IF(NumberStations>0)THEN
# if defined (SPHERICAL_IJ_STATION)
       ! use IJ for coupling purpose
       CALL STATIONS_SPHERICAL_IJ
# else
       CALL STATIONS
# endif
   ENDIF


# if defined (COUPLING)
     CLOSE(11)
# endif

# if defined (VESSEL)

DO Kves = 201,200+NumVessel
    CLOSE(Kves)
ENDDO
# endif

     ! record wall time at the end
# if defined (PARALLEL)
     if(myid.eq.0) tend = MPI_Wtime( )
# else
     CALL CPU_TIME (tend)
# endif

# if defined (PARALLEL)
     if(myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if(myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     if (myid.eq.0) WRITE(*,*)'Normal Termination!'
     if (myid.eq.0) WRITE(3,*)'Normal Termination!'
# else
     write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     WRITE(*,*)'Normal Termination!'
     WRITE(3,*)'Normal Termination!'
# endif

# if defined (PARALLEL)
     call MPI_FINALIZE ( ier )
# endif

END PROGRAM FUNWAVE_TVD


!-------------------------------------------------------------------------------------
! This part is not subroutines
!  DEFINITIONS OF VARIABLES
! 
!    Last Update: 02/18/2016 Fengyan Shi
!-------------------------------------------------------------------------------------
!
! Depth(): still water depth at element point
! DepthNode(): still water depth at node
! DepthX(): still water depth at x-interface
! DepthY(): still water depth at y-interface
! Eta():   surface elevation
! Eta0(): Eta at previous time level
!  for dry point, Eta() = MinDepth+Z()
! MASK(): 1 - wet
!         0 - dry
! MASK_STRUC(): 0 - permanent dry point
! MASK9: mask for itself and 8 elements around
! 
! U():  depth-averaged u or u at the reference level (u_alpha) at element
! V():  depth-averaged v or v at the reference level (v_alpha) at element
! HU(): (dep+eta)*u at element
! HV(): (dep+eta)*v at element
! P(): HU + dispersion at x-interface
! Q(): HV + dispersion at y-interface
! Fx(): F at x-interface
! Fy(): F at y-interface
! Gx(): G at x-interface
! Gy(): G at y-interface
! Ubar(:,:,:): Ubar
! Vbar(:,:,:): Vbar

! dispersion
! U1p(:,:): x-component of V1p
! V1p(:,:): y-component of V1p

! 
! EtaRxL(): Eta Left value at x-interface
! EtaRxR(): Eta Right value at x-interface
! EtaRyL(): Eta Left value at y-interface
! EtaRyR(): Eta Right value at y-interface
! HxL():   total depth  Left value at x-interface
! HxR():   total depth  Right value at x-interface
! HyL():   total depth  Left value at y-interface
! HyR():   total depth  Right value at y-interface

! HUxL(): HU Left value at x-interface
! HUxR(): HU Right value at x-interface
! HUyL(): HV Left value at y-interface
! HUyR(): HV Right value at y-interface

! PL(): HU + dispersion, Left value at x-interface
! PR(): HU + dispersion, Right value at x-interface
! QL(): HV + dispersion, Left value at y-interface
! QR(): HV + dispersion, Right value at y-interface

! FxL = HUxL*UxL + 1/2*g*(EtaRxL^2 + 2*EtaRxL*Depthx)
! FxR = HUxR*UxR + 1/2*g*(EtaRxR^2 + 2*EtaRxR*Depthx)
! FyL = HyL*UyL*VyL
! FyR = HyR*UyR*VyR

! GxL = HxL*UxL*VxL
! GxR = HxR*UxR*VxR
! GyL = HVyL*VyL + 1/2*g*(EtaRyL^2 + 2*EtaRyL*Depthy)
! GyR = HVyR*VyR + 1/2*g*(EtaRyR^2 + 2*EtaRyR*Depthy) 





