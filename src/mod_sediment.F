!------------------------------------------------------------------------------------
!
!      FILE mod_sediment.F
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
!  SEDIMENT is a module to model sediment transport and morpho changes    
!
!  HISTORY :
!    03/14/2016 Babak Tehranirad
!    06/09/2017 Fengyan Shi
!    02/05/2018 put in the sediment module 
!-------------------------------------------------------------------------------------
#define nop(x) associate( x => x ); end associate

# if defined (SEDIMENT)

MODULE SEDIMENT_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Mloc1,Nloc1,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,U,V,P,Q,MinDepth,MASK,DT,Gamma3,Depth,tmp4preview, &
                    ALPHA,BETA,MASK9,DepthX,DepthY,PERIODIC, &
                    UNDERTOW_U, UNDERTOW_V, ROLLER_SWITCH,ROLLER
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

    CHARACTER(LEN=80) Mask_s_FILE
    CHARACTER(LEN=80) Sed_Scheme

    LOGICAL :: Bed_Change=.FALSE.
    LOGICAL :: Avalanche=.FALSE.
    LOGICAL :: IN_Mask_s=.FALSE.
    LOGICAL :: BedLoad = .FALSE.

    INTEGER :: Counter_s=0, Counter_ava=0
    INTEGER :: Morph_step,c_dum,Morph_factor

    REAL(SP):: PLOT_INTV_SEDIMENT,PLOT_COUNT_SEDIMENT

    REAL(SP) :: D50,Sdensity,n_porosity,DT1,Counter,WS,Kappa1,Kappa2, &
                MinDepthPickup
! cohesive
    REAL(SP) :: a_coh,b_coh,n_coh,m_coh,alpha_coh,E_coh,k_coh,Tau_crd_coh,&
                Tau_cr_coh
    LOGICAL :: SoftBed = .TRUE.
    LOGICAL :: CohesiveSediment = .FALSE.
    
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Zb,P_ave,D_ave,Mask_s,D,Pickup,Delta_c
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Zb_avalanching
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::CH,CHX,CHY,CHH,CHHxR,CHHxL,CHXxL,CHXxR
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::CHHyR,CHHyL,CHYyL,CHYyR,ZbOld,dZb
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::FsL,FsR,GsL,GsR,Fs,Gs,CH0,Zs,ava
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::BedFluxX,BedFluxY
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Hpo,CHH0,TotalSuspendLoad,TotalBedLoad
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Depth_ini
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Tau_xy


    REAL(SP) :: Shields_cr,A_sed,Tau_cr,D_sed,htt,k_s,u_c,ustar_c,ec_max
    REAL(SP) :: Shields_cr_bedload,Tau_cr_bedload
    REAL(SP) :: beta_d,reduction,c_a,c_b,ec_max_w,ec_max_c
    REAL(SP) :: u_c1,u_c2,u_c3,u_c4,ustar_c1,ustar_c2,ustar_c3,ustar_c4
    REAL(SP) :: k1,k2,k3,k4,F1,F2,F3,F4,tan_phi,Dstar,viscosity
    REAL(SP) :: Satuation

#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif    


CONTAINS
  
SUBROUTINE SEDIMENT_INITIAL
  USE GLOBAL,ONLY : SMALL, INPUT_FILE_NAME
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  INTEGER :: ierr
# if !defined (PARALLEL)
  INTEGER :: I,J
# endif

! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME

! input is based on Babaks original but make default values
      CALL READ_STRING(Sed_Scheme,FILE_NAME,'Sed_Scheme',ierr)

      IF(ierr==1)THEN
        Sed_Scheme = 'Upwinding'
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
         WRITE(3,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
      endif
# else
         WRITE(*,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
         WRITE(3,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
# endif
      ENDIF

      CALL READ_LOGICAL(CohesiveSediment,FILE_NAME,'CohesiveSediment',ierr)

      IF(ierr == 1)THEN
       CohesiveSediment = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'CohesiveSediment not defined, Default: False'
       WRITE(*,'(A40)')'CohesiveSediment not defined, Default: False'
      endif
# else
      WRITE(3,'(A40)')'CohesiveSediment not defined, Default: False'
      WRITE(*,'(A40)')'CohesiveSediment not defined, Default: False'
# endif 
      ENDIF

!   for cohesive sediment: a_coh,b_coh,n_coh,m_coh,alpha_coh,E_coh     

    IF(CohesiveSediment)THEN

! SoftBed
      CALL READ_LOGICAL(SoftBed,FILE_NAME,'SoftBed',ierr)

      IF(ierr == 1)THEN
       SoftBed = .TRUE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'SoftBed not defined, Default: True'
       WRITE(*,'(A40)')'SoftBed not defined, Default: True'
      endif
# else
      WRITE(3,'(A40)')'SoftBed not defined, Default: True'
      WRITE(*,'(A40)')'SoftBed not defined, Default: True'
# endif 
      ENDIF

! tau_crd_coh, default 0.001

      CALL READ_Float(Tau_crd_coh,FILE_NAME,'Tau_crd_coh',ierr)

      IF(ierr==1)THEN
        Tau_crd_coh = 0.001_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify Tau_crd_coh, use default: 0.001'
         WRITE(3,'(A40)')'You dont specify Tau_crd_coh, use default: 0.001'
      endif
# else
         WRITE(*,'(A40)')'You dont specify Tau_crd_coh, use default: 0.001'
         WRITE(3,'(A40)')'You dont specify Tau_crd_coh, use default: 0.001'
# endif
      ENDIF 

! tau_cr_coh, default 0.001

      CALL READ_Float(Tau_cr_coh,FILE_NAME,'Tau_cr_coh',ierr)

      IF(ierr==1)THEN
        Tau_cr_coh = 0.001_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify Tau_cr_coh, use default: 0.001'
         WRITE(3,'(A40)')'You dont specify Tau_cr_coh, use default: 0.001'
      endif
# else
         WRITE(*,'(A40)')'You dont specify Tau_cr_coh, use default: 0.001'
         WRITE(3,'(A40)')'You dont specify Tau_cr_coh, use default: 0.001'
# endif
      ENDIF 

! k_coh, default: 0.000001
      CALL READ_Float(k_coh,FILE_NAME,'k_coh',ierr)

      IF(ierr==1)THEN
        k_coh = 0.000001_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify k_coh, use default: 10E-6'
         WRITE(3,'(A40)')'You dont specify k_coh, use default: 10E-6'
      endif
# else
         WRITE(*,'(A40)')'You dont specify k_coh, use default: 10E-6'
         WRITE(3,'(A40)')'You dont specify k_coh, use default: 10E-6'
# endif
      ENDIF 

! a_coh
      CALL READ_Float(a_coh,FILE_NAME,'a_coh',ierr)

      IF(ierr==1)THEN
        a_coh = 0.1_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify a_coh, use default: 0.1'
         WRITE(3,'(A40)')'You dont specify a_coh, use default: 0.1'
      endif
# else
         WRITE(*,'(A40)')'You dont specify a_coh, use default: 0.1'
         WRITE(3,'(A40)')'You dont specify a_coh, use default: 0.1'
# endif
      ENDIF 

      CALL READ_Float(b_coh,FILE_NAME,'b_coh',ierr)

      IF(ierr==1)THEN
        b_coh = 2.0_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify b_coh, use default: 2'
         WRITE(3,'(A40)')'You dont specify b_coh, use default: 2'
      endif
# else
         WRITE(*,'(A40)')'You dont specify b_coh, use default: 2'
         WRITE(3,'(A40)')'You dont specify b_coh, use default: 2'
# endif
      ENDIF 

      CALL READ_Float(n_coh,FILE_NAME,'n_coh',ierr)

      IF(ierr==1)THEN
        n_coh = 0.5_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify n_coh, use default: 0.5'
         WRITE(3,'(A40)')'You dont specify n_coh, use default: 0.5'
      endif
# else
         WRITE(*,'(A40)')'You dont specify n_coh, use default: 0.5'
         WRITE(3,'(A40)')'You dont specify n_coh, use default: 0.5'
# endif
      ENDIF 

      CALL READ_Float(m_coh,FILE_NAME,'m_coh',ierr)

      IF(ierr==1)THEN
        m_coh = 1.5_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify m_coh, use default: 1.5'
         WRITE(3,'(A40)')'You dont specify m_coh, use default: 1.5'
      endif
# else
         WRITE(*,'(A40)')'You dont specify m_coh, use default: 1.5'
         WRITE(3,'(A40)')'You dont specify m_coh, use default: 1.5'
# endif
      ENDIF

      CALL READ_Float(E_coh,FILE_NAME,'E_coh',ierr)

      IF(ierr==1)THEN
        E_coh = 0.0001_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify E_coh, use default: 0.0001'
         WRITE(3,'(A40)')'You dont specify E_coh, use default: 0.0001'
      endif
# else
         WRITE(*,'(A40)')'You dont specify E_coh, use default: 0.0001'
         WRITE(3,'(A40)')'You dont specify E_coh, use default: 0.0001'
# endif
      ENDIF

      CALL READ_Float(alpha_coh,FILE_NAME,'alpha_coh',ierr)

      IF(ierr==1)THEN
        alpha_coh = 1.0_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify alpha_coh, use default: 1.0'
         WRITE(3,'(A40)')'You dont specify alpha_coh, use default: 1.0'
      endif
# else
         WRITE(*,'(A40)')'You dont specify alpha_coh, use default: 1.0'
         WRITE(3,'(A40)')'You dont specify alpha_coh, use default: 1.0'
# endif
      ENDIF

    ENDIF ! end cohesive

      CALL READ_Float(D50,FILE_NAME,'D50',ierr)

      IF(ierr==1)THEN
       IF(CohesiveSediment)THEN
        D50 = 0.000005_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify D50, use default: 5nm'
         WRITE(3,'(A40)')'You dont specify D50, use default: 5nm'
      endif
# else
         WRITE(*,'(A40)')'You dont specify D50, use default: 5nm'
         WRITE(3,'(A40)')'You dont specify D50, use default: 5nm'
# endif
       ELSE
        D50 = 0.0005_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify D50, use default: 0.5mm'
         WRITE(3,'(A40)')'You dont specify D50, use default: 0.5mm'
      endif
# else
         WRITE(*,'(A40)')'You dont specify D50, use default: 0.5mm'
         WRITE(3,'(A40)')'You dont specify D50, use default: 0.5mm'
# endif
      ENDIF ! cohesive
     ENDIF ! ierr

      CALL READ_Float(Sdensity,FILE_NAME,'Sdensity',ierr)

      IF(ierr==1)THEN
        Sdensity = 2.68_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Sdensity, use default: 2.68'
         WRITE(3,'(A60)')'You dont specify Sdensity, use default: 2.68'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Sdensity, use default: 2.68'
         WRITE(3,'(A60)')'You dont specify Sdensity, use default: 2.68'
# endif
      ENDIF 

      CALL READ_LOGICAL(BED_CHANGE,FILE_NAME,'Bed_Change',ierr)

      IF(ierr == 1)THEN
       BED_CHANGE = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'Bed_Change not defined, Default: False'
       WRITE(*,'(A40)')'Bed_Change not defined, Default: False'
      endif
# else
      WRITE(3,'(A40)')'Bed_Change not defined, Default: False'
      WRITE(*,'(A40)')'Bed_Change not defined, Default: False'
# endif 
      ENDIF

      IF(CohesiveSediment)THEN
       BedLoad = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'BedLoad turned off for Cohesive sed'
       WRITE(*,'(A40)')'BedLoad turned off for Cohesive sed'
      endif
# else
      WRITE(3,'(A40)')'BedLoad turned off for Cohesive sed'
      WRITE(*,'(A40)')'BedLoad turned off for Cohesive sed'
# endif 

      ELSE ! non-cohesive

      CALL READ_LOGICAL(BedLoad,FILE_NAME,'BedLoad',ierr)

      IF(ierr == 1)THEN
       BedLoad = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'BedLoad not defined, Default: False'
       WRITE(*,'(A40)')'BedLoad not defined, Default: False'
      endif
# else
      WRITE(3,'(A40)')'BedLoad not defined, Default: False'
      WRITE(*,'(A40)')'BedLoad not defined, Default: False'
# endif 
      ENDIF

     ENDIF ! end cohesive

      CALL READ_Float(MinDepthPickup,FILE_NAME,'MinDepthPickup',ierr)

      IF(ierr == 1)THEN
       MinDepthPickup = 0.1_SP 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'MinDepthPickup not defined, Default: 0.1m'
       WRITE(*,'(A40)')'MinDepthPickup not defined, Default: 0.1m'
      endif
# else
      WRITE(3,'(A40)')'MinDepthPickup not defined, Default: 0.1m'
      WRITE(*,'(A40)')'MinDepthPickup not defined, Default: 0.1m'
# endif 
      ENDIF

      CALL READ_Float(n_porosity,FILE_NAME,'n_porosity',ierr)

      IF(ierr==1)THEN
        n_porosity = 0.47_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify n_porosity, use default: 0.47'
         WRITE(3,'(A60)')'You dont specify n_porosity, use default: 0.47'
      endif
# else
         WRITE(*,'(A60)')'You dont specify n_porosity, use default: 0.47'
         WRITE(3,'(A60)')'You dont specify n_porosity, use default: 0.47'
# endif
      ENDIF 

      CALL READ_Float(WS,FILE_NAME,'WS',ierr)

      IF(ierr==1)THEN
        IF(CohesiveSediment)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify WS, use formula for cohesive sed'
         WRITE(3,'(A60)')'You dont specify WS, use formula for cohesive sed'
      endif
# else
         WRITE(*,'(A60)')'You dont specify WS, use formula for cohesive sed'
         WRITE(3,'(A60)')'You dont specify WS, use formula for cohesive sed'
# endif
        ELSE ! non-cohesive
        WS = 0.0125_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify WS, use default: 0.0125m/s'
         WRITE(3,'(A60)')'You dont specify WS, use default: 0.0125m/s'
      endif
# else
         WRITE(*,'(A60)')'You dont specify WS, use default: 0.0125m/s'
         WRITE(3,'(A60)')'You dont specify WS, use default: 0.0125m/s'
# endif
        ENDIF ! end if or not cohesive
      ENDIF

      CALL READ_Float(Shields_cr,FILE_NAME,'Shields_cr',ierr)

      IF(ierr==1)THEN
        Shields_cr = 0.055_SP
       IF(.NOT.CohesiveSediment)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Shields_cr, use default: 0.055'
         WRITE(3,'(A60)')'You dont specify Shields_cr, use default: 0.055'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Shields_cr, use default: 0.055'
         WRITE(3,'(A60)')'You dont specify Shields_cr, use default: 0.055'
# endif
       ENDIF ! cohesive
      ENDIF

      CALL READ_Float(Shields_cr_bedload,FILE_NAME,'Shields_cr_bedload',ierr)

      IF(ierr==1)THEN
        Shields_cr_bedload = Shields_cr
       IF(.NOT.CohesiveSediment)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify Shields_cr_bedload, use Shields_cr'
         WRITE(3,'(A80)')'You dont specify Shields_cr_bedload, use Shields_cr'
      endif
# else
         WRITE(*,'(A80)')'You dont specify Shields_cr_bedload, use Shields_cr'
         WRITE(3,'(A80)')'You dont specify Shields_cr_bedload, use Shields_cr'
# endif
       ENDIF
      ENDIF

      CALL READ_INTEGER(Morph_step,FILE_NAME,'Morph_step',ierr)

      IF(ierr==1)THEN
        Morph_step = 25
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Morph_step, use default: 25'
         WRITE(3,'(A60)')'You dont specify Morph_step, use default: 25'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Morph_step, use default: 25'
         WRITE(3,'(A60)')'You dont specify Morph_step, use default: 25'
# endif
      ENDIF

      CALL READ_INTEGER(Morph_factor,FILE_NAME,'Morph_factor',ierr)

      IF(ierr==1)THEN
        Morph_factor = 1
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Morph_factor, use default: 1'
         WRITE(3,'(A60)')'You dont specify Morph_factor, use default: 1'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Morph_factor, use default: 1'
         WRITE(3,'(A60)')'You dont specify Morph_factor, use default: 1'
# endif
      ENDIF

      CALL READ_LOGICAL(Avalanche,FILE_NAME,'Avalanche',ierr)
      IF(ierr == 1)THEN
       Avalanche = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'Avalanche not defined, Default: False'
       WRITE(*,'(A40)')'Avalanche not defined, Default: False'
      endif
# else
      WRITE(3,'(A40)')'Avalanche not defined, Default: False'
      WRITE(*,'(A40)')'Avalanche not defined, Default: False'
# endif 
      ENDIF

      CALL READ_Float(tan_phi,FILE_NAME,'Tan_phi',ierr)

      IF(ierr==1)THEN
        Tan_phi = 0.7_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Tan_phi, use default: 0.7'
         WRITE(3,'(A60)')'You dont specify Tan_phi, use default: 0.7'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Tan_phi, use default: 0.7'
         WRITE(3,'(A60)')'You dont specify Tan_phi, use default: 0.7'
# endif
      ENDIF

      CALL READ_LOGICAL(IN_Mask_s,FILE_NAME,'Hard_bottom',ierr)

      IF(ierr==1)THEN
        IN_Mask_s = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Hard_bottom, use default: False'
         WRITE(3,'(A60)')'You dont specify Hard_bottom, use default: False'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Hard_bottom, use default: False'
         WRITE(3,'(A60)')'You dont specify Hard_bottom, use default: False'
# endif
      ENDIF

     IF(IN_Mask_s)THEN
      CALL READ_STRING(Mask_s_File,FILE_NAME,'Hard_bottom_file',ierr) 

      IF(ierr==1)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
         WRITE(3,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
      endif
       call MPI_FINALIZE ( ier )
# else
         WRITE(*,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
         WRITE(3,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
# endif
        STOP
      ENDIF ! end ierr

     ENDIF ! end IN_Mask_s


      CALL READ_Float(Kappa1,FILE_NAME,'Kappa1',ierr)

      IF(ierr==1)THEN
        Kappa1 = 0.3333_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
         WRITE(3,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
      endif
# else
         WRITE(*,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
         WRITE(3,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
# endif
      ENDIF

      CALL READ_Float(Kappa2,FILE_NAME,'Kappa2',ierr)
      IF(ierr==1)THEN
        Kappa2 = 1.0_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
         WRITE(3,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
         WRITE(3,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
# endif
      ENDIF

! plot sediment intitial
     PLOT_COUNT_SEDIMENT = 0
     CALL READ_FLOAT(PLOT_INTV_SEDIMENT,FILE_NAME,'PLOT_INTV_SEDIMENT',ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A80)')'PLOT_INTV_SEDIMENT not specified, use SAME'
# else
      WRITE(3,'(A50)')'PLOT_INTV_SEDIMENT not specified, use SAME'
# endif
       PLOT_INTV_SEDIMENT = SMALL
     ENDIF

! summary of sediment parameters

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) '----- Sediment Parameters'
    IF(CohesiveSediment)THEN
      WRITE(3,*) 'Cohesive Sediment'
      WRITE(3,*) 'SoftBed', SoftBed
      WRITE(3,*) 'Critical deposition stress: Tau_crd_coh ', Tau_crd_coh
      WRITE(3,*) 'Critical pickup stress: Tau_cr_coh ', Tau_cr_coh
      WRITE(3,*) 'Diffusion coefficient: k_coh ', k_coh
      WRITE(3,*) 'flocs settling parameter: a_coh ', a_coh
      WRITE(3,*) 'flocs settling parameter: b_coh ', b_coh
      WRITE(3,*) 'flocs settling parameter: m_coh ', m_coh
      WRITE(3,*) 'flocs settling parameter: n_coh ', n_coh
      WRITE(3,*) 'Erosion rate: E_coh ', E_coh
      WRITE(3,*) 'alpha for softbed: alpha_coh ', alpha_coh

      WRITE(3,*) 'D50:  ', D50
      WRITE(3,*) 'specific gravity of sediment, Sdensity ', Sdensity
      WRITE(3,*) 'Bed_Change:  ', BED_CHANGE
      WRITE(3,*) 'Minimum Depth  for Pickup: MinDepthPickup  ', MinDepthPickup
      WRITE(3,*) 'BedLoad:  ', BedLoad
      WRITE(3,*) 'Porosity:  ', n_porosity
      WRITE(3,*) 'Settling velocity: calculated based on formula '
      WRITE(3,*) 'Morph_step:  ', Morph_step
      WRITE(3,*) 'Morph_factor:  ', Morph_factor
      WRITE(3,*) 'Avalanche:  ', Avalanche
      WRITE(3,*) 'Tan_phi for Avalanche:  ', Tan_phi
      WRITE(3,*) 'Hard_bottom:  ', IN_Mask_s
      WRITE(3,*) '----- end sedment parameters  '
 
    ELSE

      WRITE(3,*) 'Non-cohesive Sediment'
      WRITE(3,*) 'D50:  ', D50
      WRITE(3,*) 'specific gravity of sediment, Sdensity ', Sdensity
      WRITE(3,*) 'Bed_Change:  ', BED_CHANGE
      WRITE(3,*) 'Minimum Depth  for Pickup: MinDepthPickup  ', MinDepthPickup
      WRITE(3,*) 'BedLoad:  ', BedLoad
      WRITE(3,*) 'Porosity:  ', n_porosity
      WRITE(3,*) 'Settling velocity Ws: ', WS
      WRITE(3,*) 'Critical Shields number:  ', Shields_cr
      WRITE(3,*) 'Critical Shields number for bedload:  ', Shields_cr_bedload
      WRITE(3,*) 'Morph_step:  ', Morph_step
      WRITE(3,*) 'Morph_factor:  ', Morph_factor
      WRITE(3,*) 'Avalanche:  ', Avalanche
      WRITE(3,*) 'Tan_phi for Avalanche:  ', Tan_phi
      WRITE(3,*) 'Hard_bottom:  ', IN_Mask_s
      WRITE(3,*) '----- end sedment parameters  '
 
    ENDIF ! end cohesive 

   ENDIF ! end myid
# else
 
   WRITE(3,*) '----- Sediment Parameters'
    IF(CohesiveSediment)THEN
      WRITE(3,*) 'Cohesive Sediment'
      WRITE(3,*) 'SoftBed', SoftBed
      WRITE(3,*) 'Critical deposition stress: Tau_crd_coh ', Tau_crd_coh
      WRITE(3,*) 'Critical pickup stress: Tau_cr_coh ', Tau_cr_coh
      WRITE(3,*) 'Diffusion coefficient: k_coh ', k_coh
      WRITE(3,*) 'flocs settling parameter: a_coh ', a_coh
      WRITE(3,*) 'flocs settling parameter: b_coh ', b_coh
      WRITE(3,*) 'flocs settling parameter: m_coh ', m_coh
      WRITE(3,*) 'flocs settling parameter: n_coh ', n_coh
      WRITE(3,*) 'Erosion rate: E_coh ', E_coh
      WRITE(3,*) 'alpha for softbed: alpha_coh ', alpha_coh

      WRITE(3,*) 'D50:  ', D50
      WRITE(3,*) 'specific gravity of sediment, Sdensity ', Sdensity
      WRITE(3,*) 'Bed_Change:  ', BED_CHANGE
      WRITE(3,*) 'Minimum Depth  for Pickup: MinDepthPickup  ', MinDepthPickup
      WRITE(3,*) 'BedLoad:  ', BedLoad
      WRITE(3,*) 'Porosity:  ', n_porosity
      WRITE(3,*) 'Settling velocity: calculated based on formula '
      WRITE(3,*) 'Morph_step:  ', Morph_step
      WRITE(3,*) 'Morph_factor:  ', Morph_factor
      WRITE(3,*) 'Avalanche:  ', Avalanche
      WRITE(3,*) 'Tan_phi for Avalanche:  ', Tan_phi
      WRITE(3,*) 'Hard_bottom:  ', IN_Mask_s
      WRITE(3,*) '----- end sedment parameters  '
 
    ELSE

      WRITE(3,*) 'Non-cohesive Sediment'
      WRITE(3,*) 'D50:  ', D50
      WRITE(3,*) 'specific gravity of sediment, Sdensity ', Sdensity
      WRITE(3,*) 'Bed_Change:  ', BED_CHANGE
      WRITE(3,*) 'Minimum Depth  for Pickup: MinDepthPickup  ', MinDepthPickup
      WRITE(3,*) 'BedLoad:  ', BedLoad
      WRITE(3,*) 'Porosity:  ', n_porosity
      WRITE(3,*) 'Settling velocity Ws: ', WS
      WRITE(3,*) 'Critical Shields number:  ', Shields_cr
      WRITE(3,*) 'Critical Shields number for bedload:  ', Shields_cr_bedload
      WRITE(3,*) 'Morph_step:  ', Morph_step
      WRITE(3,*) 'Morph_factor:  ', Morph_factor
      WRITE(3,*) 'Avalanche:  ', Avalanche
      WRITE(3,*) 'Tan_phi for Avalanche:  ', Tan_phi
      WRITE(3,*) 'Hard_bottom:  ', IN_Mask_s
      WRITE(3,*) '----- end sedment parameters  '
 
    ENDIF ! end cohesive 

# endif

! ALLOCATION

        ALLOCATE(Zb(Mloc,Nloc),Zb_avalanching(Mloc,Nloc))
        ALLOCATE(CHH(Mloc,Nloc),CHH0(Mloc,Nloc))
        ALLOCATE(CHHxR(Mloc1,Nloc),CHHxL(Mloc1,Nloc),  &
                 CHHyR(Mloc,Nloc1),CHHyL(Mloc,Nloc1))
        ALLOCATE(CHY(Mloc,Nloc),CH(Mloc,Nloc), &
                 CHXxL(Mloc1,Nloc),CHXxR(Mloc1,Nloc))
        ALLOCATE(CHX(Mloc,Nloc),CHYyL(Mloc,Nloc1),CHYyR(Mloc,Nloc1)) 
        ALLOCATE(Mask_s(Mloc,Nloc))
        ALLOCATE(D(Mloc,Nloc),Pickup(Mloc,Nloc))
        ALLOCATE(Delta_c(Mloc,Nloc),D_ave(Mloc,Nloc),P_ave(Mloc,Nloc))
        ALLOCATE(FsL(Mloc1,Nloc),FsR(Mloc1,Nloc), &
                 GsL(Mloc,Nloc1),GsR(Mloc,Nloc1))
        ALLOCATE(Gs(Mloc,Nloc1),Fs(Mloc1,Nloc),CH0(Mloc,Nloc))
        ALLOCATE(ZbOld(Mloc,Nloc),dZb(Mloc,Nloc), &
                 Zs(Mloc,Nloc),ava(Mloc,Nloc),Hpo(Mloc,Nloc), &
                 TotalSuspendLoad(Mloc,Nloc),TotalBedLoad(Mloc,Nloc))
        ALLOCATE(BedFluxX(Mloc,Nloc),BedFluxY(Mloc,Nloc),Depth_ini(Mloc,Nloc))
        ALLOCATE(Tau_xy(Mloc,Nloc))

        Depth_ini = Depth
        ava=ZERO
        Zs=LARGE
        ZbOld=ZERO
        dZb=ZERO
        CH0=ZERO
        CHXxR=ZERO
        CHXxL=ZERO
        CH=ZERO
        CHX=ZERO
        CHXxR=ZERO
        CHXxL=ZERO
        CHYyR=ZERO
        CHYyL=ZERO   
        CHY=ZERO
        CHH=ZERO
        CHHxR=ZERO
        CHHxL=ZERO
        CHHyR=ZERO
        CHHyL=ZERO
        FsL=ZERO
        FsR=ZERO
        GsL=ZERO
        GsR=ZERO
        Fs=ZERO
        Gs=ZERO

        D=ZERO
        Pickup=ZERO
        D_ave=ZERO
        P_ave=ZERO
        Delta_c=ZERO
        Mask_s=ZERO
        Zb=ZERO
        Zb_avalanching = ZERO
        TotalSuspendLoad = ZERO
        TotalBedLoad = ZERO
        BedFluxX = ZERO
        BedFluxY = ZERO

! sediment parameters 
     IF (CohesiveSediment) THEN
      viscosity=k_coh
     ELSE
      viscosity=0.000001_SP
     ENDIF  ! end cohesive

! the following parameters are for non-cohesive
      Dstar = D50*((Sdensity-1.0_SP)*grav/viscosity**2.0_SP)**(1.0_SP/3.0_SP)
      Tau_cr=(Sdensity-1.0_SP)*grav*D50*Shields_cr
      Tau_cr_bedload=(Sdensity-1.0_SP)*grav*D50*Shields_cr_bedload
      k_s = 2.5_SP*D50 
      Satuation = 0.001

     IF(IN_Mask_s)THEN
! check existing

 INQUIRE(FILE=TRIM(Mask_s_File),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(Mask_s_File), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(Mask_s_File), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF  ! exist

# if defined (PARALLEL)
     call GetFile (Mask_s_File,Zs)
# else
     OPEN(1,FILE=TRIM(Mask_s_File))

       DO J=Jbeg,Jend
        READ(1,*)(Zs(I,J),I=Ibeg,Iend)
       ENDDO

     CLOSE(1)
# endif

     ENDIF ! end IN_Mask_s



END SUBROUTINE SEDIMENT_INITIAL

SUBROUTINE SEDIMENT_ADVECTION_DIFFUSION(ISTEP)
  USE GLOBAL,ONLY : SMALL
  IMPLICIT NONE
  INTEGER::ISTEP,IVAR
  REAL(SP) :: gamma_cao
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Scalx,Scaly
  REAL(SP),DIMENSION(Ibeg:Iend,Jbeg:Jend)::R1
  REAL(SP) :: Weight0 = 0.9_SP, Weight1=0.1_SP
  REAL(SP) :: Angle_cur,bedf,fr,Pd,Cb
  INTEGER :: I,J
  REAL(SP) :: tmp1,tmp2

    ALLOCATE(Scalx(Mloc1,Nloc),Scaly(Mloc,Nloc1))

    H=Eta*Gamma3+Depth
    DO J=1,Nloc
    DO I=1,Mloc
      Hpo(I,J)=MAX(H(I,J),MinDepth)
    ENDDO
    ENDDO

 
  IF(Sed_Scheme(1:3)=='Upw')THEN

! advection in x direction
     Scalx = ZERO
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend+1
      IF(ROLLER) THEN
       tmp1=P(I,J)+0.5_SP*(UNDERTOW_U(I-1,J)+UNDERTOW_U(I,J))
      ELSE
       tmp1=P(I,J)
      ENDIF
      IF(tmp1>=0.0)THEN
        IF(MASK(I-1,J)==0)THEN
          Scalx(I,J) = ZERO
        ELSE
          Scalx(I,J) = tmp1*CH(I-1,J)
        ENDIF
      ELSE
        IF(MASK(I,J)==0)THEN
          Scalx(I,J) = ZERO
        ELSE
          Scalx(I,J) = tmp1*CH(I,J)
        ENDIF
      ENDIF
     ENDDO
     ENDDO

! advection in y direction
     Scaly = ZERO
     DO J=Jbeg,Jend+1
     DO I=Ibeg,Iend
      IF(ROLLER) THEN
        tmp2=Q(I,J)+0.5_SP*(UNDERTOW_V(I,J-1)+UNDERTOW_V(I,J))
      ELSE
        tmp2=Q(I,J)
      ENDIF
      IF(tmp2>=0.0)THEN
        IF(MASK(I,J-1)==0)THEN
          Scaly(I,J) = ZERO
        ELSE
          Scaly(I,J) = tmp2*CH(I,J-1)
        ENDIF
      ELSE
        IF(MASK(I,J)==0)THEN
          Scaly(I,J) = ZERO
        ELSE
          Scaly(I,J) = tmp2*CH(I,J)
        ENDIF
      ENDIF
     ENDDO
     ENDDO

  ELSE  ! TVD, minmod

! advection in x direction
     Scalx = ZERO
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend+1
      IF(ROLLER)THEN
        tmp1=P(I,J)+0.5_SP*(UNDERTOW_U(I-1,J)+UNDERTOW_U(I,J))
      ELSE
        tmp1=P(I,J)
      ENDIF
      IF(tmp1>=ZERO)THEN
        IF(MASK(I-1,J)==0)THEN
          Scalx(I,J) = ZERO
        ELSE
          Scalx(I,J) = tmp1*(Weight0*CH(I-1,J)+Weight1*CH(I,J))
        ENDIF
      ELSE
        IF(MASK(I,J)==0)THEN
          Scalx(I,J) = ZERO
        ELSE
          Scalx(I,J) = tmp1*(Weight0*CH(I,J)+Weight1*CH(I-1,J))
        ENDIF
      ENDIF
     ENDDO
     ENDDO

! advection in y direction
     Scaly = ZERO
     DO J=Jbeg,Jend+1
     DO I=Ibeg,Iend
      IF(ROLLER)THEN
        tmp2=Q(I,J)+0.5_SP*(UNDERTOW_V(I,J-1)+UNDERTOW_V(I,J))
      ELSE
        tmp2=Q(I,J)
      ENDIF
      IF(tmp2>=ZERO)THEN
        IF(MASK(I,J-1)==0)THEN
          Scaly(I,J) = ZERO
        ELSE
          Scaly(I,J) = tmp2*(Weight0*CH(I,J-1)+Weight1*CH(I,J))
        ENDIF
      ELSE
        IF(MASK(I,J)==0)THEN
          Scaly(I,J) = ZERO
        ELSE
          Scaly(I,J) = tmp2*(Weight0*CH(I,J)+Weight1*CH(I,J-1))

        ENDIF
      ENDIF
     ENDDO
     ENDDO

  ENDIF

! diffusion (roller effect is not included in diffusion)

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend

       IF(MASK(I,J).GT.0)THEN
            u_c =SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
            ustar_c = 0.4_SP*u_c/(-1.0_SP+log(30.0_SP*Hpo(I,J)/k_s))
! diffusion x
         IF(MASK(I-1,J).GT.0)THEN
            u_c2=SQRT(U(I-1,J)*U(I-1,J)+V(I-1,J)*V(I-1,J))
            ustar_c2= 0.4_SP*u_c2/(-1.0_SP+log(30.0_SP*Hpo(I-1,J)/k_s))
            k2=5.93_SP*(ustar_c2+ustar_c)*(Hpo(I-1,J)+Hpo(I,J))/4.0_SP
# if defined (CARTESIAN)
         Scalx(I,J)=Scalx(I,J) &
              - k2*(Hpo(I-1,J)+Hpo(I,J))*(CH(I,J)-CH(I-1,J))/2.0_SP/DX
# else
         Scalx(I,J)=Scalx(I,J) &
              - k2*(Hpo(I-1,J)+Hpo(I,J))*(CH(I,J)-CH(I-1,J))/2.0_SP/DX(I,J)
# endif

         ENDIF ! mask x

         IF(MASK(I,J-1).GT.0)THEN
            u_c4=SQRT(U(I,J-1)*U(I,J-1)+V(I,J-1)*V(I,J-1))
            ustar_c4= 0.4_SP*u_c4/(-1.0_SP+log(30.0_SP*Hpo(I,J-1)/k_s))
            k4=5.93_SP*(ustar_c4+ustar_c)*(Hpo(I,J)+Hpo(I,J-1))/4.0_SP
# if defined (CARTESIAN)
         Scaly(I,J)=Scaly(I,J) &
              - k4*(Hpo(I,J-1)+Hpo(I,J))*(CH(I,J)-CH(I,J-1))/2.0_SP/DY
# else
         Scaly(I,J)=Scaly(I,J) &
              - k4*(Hpo(I,J-1)+Hpo(I,J))*(CH(I,J)-CH(I,J-1))/2.0_SP/DY(I,J)
# endif

         ENDIF  ! mask y
       ENDIF !  mask

    ENDDO
    ENDDO

! bc
     IVAR = 1
     CALL FLUX_SCALAR_BC(IVAR,Scalx,Scaly)

! solver

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      
       IF(MASK(I,J).GT.0)THEN

# if defined (CARTESIAN)
       R1(I,J)= - ((Scalx(i+1,j)-Scalx(i,j))/DX  &
                            + (Scaly(i,j+1)-Scaly(i,j))/DY) &
                + Pickup(I,J)-D(I,J)
# else
       R1(I,J)= - ((Scalx(i+1,j)-Scalx(i,j))/DX(I,J)  &
                            + (Scaly(i,j+1)-Scaly(i,j))/DY(I,J)) &
                + Pickup(I,J)-D(I,J)
# endif

        CHH(I,J)=ALPHA(ISTEP)*CHH0(I,J)+BETA(ISTEP)*(CHH(I,J)+DT*R1(I,J))

        IF(CHH(I,J)<ZERO) CHH(I,J)=ZERO

        CH(I,J)=CHH(I,J)/Hpo(I,J)

      ENDIF ! mask

     ENDDO
     ENDDO

# if defined (PARALLEL)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,CH,1,PERIODIC)
# endif


! sources

      BedFluxX=ZERO
      BedFluxY=ZERO

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
    IF(MASK(I,J).GT.0.AND.H(I,J)>MinDepthPickup)THEN
          u_c=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))

! avoid sheet flow condition
          fr=SQRT(H(I,J)*grav)*0.5
          IF(u_c > fr) u_c = fr

    ! note rho_w is removed, also removed in tau_cr, 

          Tau_xy(I,J)= 0.16_SP/(1.0+log(k_s/(30.0_SP*Hpo(I+1,J))))**2*(u_c**2.)

       IF (CohesiveSediment)THEN

          IF (Tau_xy(I,J).GT.Tau_cr_coh) THEN
            IF(SoftBed)THEN
              Pickup(I,J)=E_coh*EXP(alpha_coh*SQRT(ABS(Tau_xy(I,J)-tau_cr_coh)))
            ELSE                           
              Pickup(I,J)=E_coh*(Tau_xy(I,J)/Max(Tau_cr_coh,SMALL)-1.0_SP)
            ENDIF
          ELSE 
            Pickup(I,J)=0.0_SP
          ENDIF ! end tau_cr

       ELSE ! non-cohesive

!  suspended load
          IF (Tau_xy(I,J).GT.Tau_cr) THEN
            c_b=0.015_SP*(((Tau_xy(I,J)-Tau_cr)/Tau_cr)**1.5_SP)*Dstar**(-0.3_SP)
            reduction=MIN(1.0_SP,0.65_SP/c_b)
            c_a=reduction*c_b*D50/(MAX(0.01_SP*Hpo(I,J),0.01_SP))
            Pickup(I,J)=MAX(0.0_SP,c_a*WS)
! satuation
            IF(Pickup(I,J)*DT*0.3333_SP/Hpo(I,J)>=Satuation)  &
               Pickup(I,J)=Satuation*Hpo(I,J)*3.0_SP/DT

          ELSE
            Pickup(I,J)=0.0_SP
          ENDIF

         IF(BedLoad)THEN
!  bedload

          IF (Tau_xy(I,J).GT.Tau_cr_bedload) THEN
             Angle_cur = ATAN2(V(I,J),U(I,J))
             bedf = 8.0_SP*(Tau_xy(I,J)-tau_cr_bedload)**1.5/grav/(Sdensity-1.0_SP)
             BedFluxX(I,J) = bedf*COS(Angle_cur)
             BedFluxY(I,J) = bedf*SIN(Angle_cur)
          ELSE
             BedFluxX(I,J) = ZERO
             BedFluxY(I,J) = ZERO
          ENDIF
         ENDIF ! end bedload

       ENDIF ! end cohesive

    ELSE ! dry point
          Pickup(I,J)=0.0_SP
    ENDIF ! end dry point

        IF(IN_Mask_s) THEN
          IF(Zb(I,J).GE.(Zs(I,J)-0.001_SP)) THEN
             Pickup(I,J)=0.0_SP
             IF(BedLoad)THEN
               BedFluxX(I,J) = ZERO
               BedFluxY(I,J) = ZERO
             ENDIF
          ENDIF

        ENDIF

    ENDDO
    ENDDO

# if defined (PARALLEL)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,BedFluxX,1,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,BedFluxY,1,PERIODIC)
# endif

! Calculate the Deposition Rate D(I,J)
      DO J=Jbeg,Jend+1
      DO I=Ibeg,Iend+1
        IF(MASK(I,J).GT.0)THEN

         IF(CohesiveSediment)THEN
           Cb=CH(I,J)*Sdensity*RHO_WATER
           WS = a_coh*Cb**n_coh/MAX(SMALL, (Cb**2+b_coh**2)**m_coh)
           Pd = 1.0_SP - Tau_xy(I,J)/MAX(SMALL,Tau_crd_coh)
           D(I,J) = WS*CH(I,J)*Pd

         ELSE ! non-cohesive
!  Cao(2004)
          gamma_cao = MIN(2.0_SP,(1.0_SP-n_porosity)/MAX(SMALL,CH(I,J)))
!
!         D(I,J) = ZERO      
          D(I,J)=gamma_cao*CH(I,J)*WS* &
                 (1-gamma_cao*CH(I,J))**2.0_SP

         ENDIF ! end non-cohesive

        ELSE
          D(I,J)=0.0_SP
        ENDIF ! mask
      ENDDO
      ENDDO

      DEALLOCATE(Scalx,Scaly)

END SUBROUTINE SEDIMENT_ADVECTION_DIFFUSION

SUBROUTINE FLUX_SCALAR_BC(IVAR,Scalx,Scaly)
!--------------------------------------------------------
!   Specify boundary conditions for scalar convection
! 
!   fyshi 02/07/2018
!-------------------------------------------------------
    USE GLOBAL, ONLY : Mloc1,Mloc,Nloc1,Nloc,Ibeg,Iend,  &
                       Iend1,Jbeg,Jend,Jend1
# if defined (PARALLEL)
    USE GLOBAL, ONLY : n_west,n_east,n_nrth,n_suth
# endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IVAR
    REAL(SP), DIMENSION(Mloc1,Nloc), INTENT(INOUT) :: Scalx
    REAL(SP), DIMENSION(Mloc,Nloc1), INTENT(INOUT) :: Scaly
    INTEGER :: I,J
    nop(IVAR)

    ! west
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     DO J = Jbeg,Jend
         Scalx(Ibeg,J) = ZERO
     ENDDO
# if defined (PARALLEL)
     endif
# endif

    ! east
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     DO J = Jbeg,Jend
         Scalx(Iend1,J) = ZERO
     ENDDO
# if defined (PARALLEL)
     endif
# endif

    ! south
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     DO I = Ibeg,Iend
         Scaly(I,Jbeg) = ZERO
     ENDDO
# if defined (PARALLEL)
     endif
# endif

    ! north
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     DO I = Ibeg,Iend
         Scaly(I,Jend1) = ZERO
     ENDDO
# if defined (PARALLEL)
     endif
# endif

! mask
    DO J = Jbeg,Jend
    DO I = Ibeg,Iend
      IF(MASK(I,J)==0) THEN
        Scalx(I,J) = ZERO
        Scalx(I+1,J) = ZERO
        Scaly(I,J) = ZERO
        Scaly(I,J+1) = ZERO
      ENDIF
    ENDDO
    ENDDO

    RETURN
END SUBROUTINE FLUX_SCALAR_BC

SUBROUTINE MORPHOLOGICAL_CHANGE
  IMPLICIT NONE
  INTEGER :: VTYPE,i4,i4_record
  REAL(SP) :: slope_max,dh
  REAL(SP),DIMENSION(4) :: slope4
  INTEGER :: I,J

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
      TotalSuspendLoad(I,J) = TotalSuspendLoad(I,J)+(-Pickup(I,J)+D(I,J))*DT
# if defined (CARTESIAN)
      TotalBedLoad(I,J) = TotalBedLoad(I,J)  &
                    -(BedFluxX(I+1,J)-BedFluxX(I-1,J))/2.0_SP/DX*DT  &
                    -(BedFluxY(I,J+1)-BedFluxY(I,J-1))/2.0_SP/DY*DT
# else
      TotalBedLoad(I,J) = TotalBedLoad(I,J)  &
                    -(BedFluxX(I+1,J)-BedFluxX(I-1,J))/2.0_SP/DX(I,J)*DT  &
                    -(BedFluxY(I,J+1)-BedFluxY(I,J-1))/2.0_SP/DY(I,J)*DT
# endif

           ! zb: positive for erosion
      Zb(I,J) = -(TotalSuspendLoad(I,J) + TotalBedLoad(I,J))/(1.0_SP-n_porosity)

      IF(Zb(I,J)>Zs(I,J)) Zb(I,J)=Zs(I,J)
 
    ENDDO
    ENDDO

! avalanching
    IF(Avalanche)THEN
    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined(CARTESIAN)
        slope4(1)=(Depth(I,J)-Depth(I-1,J))/DX
        slope4(2)=(Depth(I,J)-Depth(I+1,J))/DX
        slope4(3)=(Depth(I,J)-Depth(I,J-1))/DY
        slope4(4)=(Depth(I,J)-Depth(I,J+1))/DY
# else
        slope4(1)=(Depth(I,J)-Depth(I-1,J))/DX(I,J)
        slope4(2)=(Depth(I,J)-Depth(I+1,J))/DX(I,J)
        slope4(3)=(Depth(I,J)-Depth(I,J-1))/DY(I,J)
        slope4(4)=(Depth(I,J)-Depth(I,J+1))/DY(I,J)
# endif
        slope_max=0.0_SP
        dh=0.0_SP
        DO i4=1,4
          IF(slope4(i4)>slope_max)THEN
           slope_max=slope4(i4)
           i4_record = i4
          ENDIF
        ENDDO
        
        IF(slope_max>Tan_phi)THEN 
          IF(i4_record.EQ.1) THEN
            IF(Zb(I-1,J)<Zs(I-1,J))THEN
              dh=0.5_SP*(Depth(I,J)-Depth(I-1,J))
              Zb_avalanching(I,J)=Zb_avalanching(I,J)+dh
              Zb_avalanching(I-1,J)=Zb_avalanching(I-1,J)-dh
            ENDIF
          ENDIF
          IF(i4_record.EQ.2) THEN
            IF(Zb(I+1,J)<Zs(I+1,J))THEN
              dh=0.5_SP*(Depth(I,J)-Depth(I+1,J))
              Zb_avalanching(I,J)=Zb_avalanching(I,J)+dh
              Zb_avalanching(I+1,J)=Zb_avalanching(I+1,J)-dh
            ENDIF
          ENDIF
          IF(i4_record.EQ.3) THEN
            IF(Zb(I,J-1)<Zs(I,J-1))THEN
              dh=0.5_SP*(Depth(I,J)-Depth(I,J-1))
              Zb_avalanching(I,J)=Zb_avalanching(I,J)+dh
              Zb_avalanching(I,J)=Zb_avalanching(I,J-1)-dh
            ENDIF
          ENDIF
          IF(i4_record.EQ.4) THEN
            IF(Zb(I,J+1)<Zs(I,J+1))THEN
              dh=0.5_SP*(Depth(I,J)-Depth(I,J+1))
              Zb_avalanching(I,J)=Zb_avalanching(I,J)+dh
              Zb_avalanching(I,J)=Zb_avalanching(I,J+1)-dh
            ENDIF
          ENDIF 
        ENDIF ! > tan_phi
    ENDDO
    ENDDO
    ENDIF
! end avalanche

! update depth

    Depth = Depth_ini + (Zb+Zb_avalanching)*Morph_factor

    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Depth,VTYPE,PERIODIC)

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


END SUBROUTINE MORPHOLOGICAL_CHANGE


END MODULE SEDIMENT_MODULE

#endif 
! end sediment

