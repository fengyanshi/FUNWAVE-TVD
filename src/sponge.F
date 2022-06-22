!------------------------------------------------------------------------------------
!
!      FILE sponge.F
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
!    SPONGE_DAMPING is subroutine for dampping waves using Larsen-Dancy(1983)
!    type sponge layer 
!    
!    HISTORY: 10/27/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
#define nop(x) associate( x => x ); end associate

SUBROUTINE SPONGE_DAMPING
     USE GLOBAL
     IMPLICIT NONE
     INTEGER :: I,J

     DO J=1,Nloc
     DO I=1,Mloc
      IF(MASK(I,J)>ZERO)THEN
       ETA(I,J)=ETA(I,J)/SPONGE(I,J)
      ENDIF
       U(I,J)=U(I,J)/SPONGE(I,J)
       V(I,J)=V(I,J)/SPONGE(I,J)
     ENDDO
     ENDDO

END SUBROUTINE SPONGE_DAMPING

!-------------------------------------------------------------------------------------
!
!    ABSORBING_GENERATING_BC is subroutine for wave generation using Larsen-Dancy(1983)
!      type wave-maker
!    
!    HISTORY: 
!      05/01/2014 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE ABSORBING_GENERATING_BC
     USE GLOBAL, ONLY : ETA,U,V,SPONGEMAKER,Mloc,Nloc,ZERO
     USE GLOBAL, ONLY : SP,PI,TIME,NumFreq,Phase_Ser,&
                       Segma_Ser,Cm_eta,Sm_eta,&
                       Cm_u,Sm_u,Cm_v,Sm_v,tmp4preview, &
                       HU, HV,Depth

     IMPLICIT NONE
     INTEGER :: I,J,KK
     real(SP),DIMENSION(Mloc,Nloc) :: Ein2D
     real(SP),DIMENSION(Mloc,Nloc) :: Uin2D,Vin2D
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(SP),DIMENSION(NumFreq) :: BB,CC                  ! Salatin et al. 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     Ein2D = ZERO
     Uin2D = ZERO
     Vin2D = ZERO
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Oct, 2021
!Salatin, R., Chen, Q., Bak, A. S., Shi, F., & Brandt, S. R. (2021). Effects of
!wave coherence on longshore variability of nearshore wave processes. Journal
!of Geophysical Research: Oceans,126, e2021JC017641.
!https://doi.org/10.1029/2021JC017641

! Eliminate redundant operations to speed up!

    ! eliminate sine and cosine from redundant calculations
    DO KK = 1, NumFreq
        BB(KK) = COS(pi/2.0_SP+Segma_Ser(KK)*TIME+Phase_Ser(KK))
        CC(KK) = SIN(pi/2.0_SP+Segma_Ser(KK)*TIME+Phase_Ser(KK))
    ENDDO
    ! use constant values to speed up
    DO J = 1,Nloc
        DO I = 1,Mloc
            DO KK = 1, NumFreq
                Ein2D(I,J)=Ein2D(I,J)+Cm_eta(I,J,KK)*BB(KK)+Sm_eta(I,J,KK)*&
                CC(KK)
                Uin2D(I,J)=Uin2D(I,J)+Cm_u(I,J,KK)*BB(KK)+Sm_u(I,J,KK)*CC(KK)
                Vin2D(I,J)=Vin2D(I,J)+Cm_v(I,J,KK)*BB(KK)+Sm_v(I,J,KK)*CC(KK)
            ENDDO
            !Ein2D(I,J)=Ein2D(I,J)/SPONGE(I,J)**4
            !Uin2D(I,J)=Uin2D(I,J)/SPONGE(I,J)**4
            !Vin2D(I,J)=Vin2D(I,J)/SPONGE(I,J)**4
            Eta(I,J) = Ein2D(I,J)+ (Eta(I,J)-Ein2D(I,J))/SpongeMaker(I,J)
            !U(I,J) = Uin2D(I,J)+(U(I,J)-Uin2D(I,J))/SpongeMaker(I,J)
            !V(I,J) = Vin2D(I,J)+(V(I,J)-Vin2D(I,J))/SpongeMaker(I,J)
            HU(I,J)=(Depth(I,J)+ETA(I,J))*U(I,J)
            HV(I,J)=(Depth(I,J)+ETA(I,J))*V(I,J)
        ENDDO
    ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       tmp4preview = Uin2D

END SUBROUTINE ABSORBING_GENERATING_BC

!-------------------------------------------------------------------------------------
!    CALCULATE_SPONGE_MAKER is subroutine for calculation of 
!      spong layer coefficient
!    
!    HISTORY: 
!      05/01/2014 Fengyan Shi
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_SPONGE_MAKER(M,N,Nghost,DX,DY,&
                            Sponge_west_width, &
                            R_sponge,A_sponge,SPONGE)
     USE PARAM
# if defined (PARALLEL)
     USE GLOBAL, ONLY : iista   !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN)::DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY
# endif
     REAL(SP),INTENT(IN) :: &
                          Sponge_west_width, &
                          R_sponge,A_sponge
     REAL(SP),DIMENSION(M,N),INTENT(INOUT)::SPONGE
     REAL(SP)::ri,lim
     INTEGER::I,J,Iwidth
     nop(DX)
     nop(DY)


# if defined (CARTESIAN)
     Iwidth=INT(Sponge_west_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_west_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,Iwidth
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri = R_Sponge**(50*(i+npx*Mglob/px-1)/(Iwidth-1))
       ri = R_Sponge**(50*(i + (iista - 1) -1)/(Iwidth-1))
!---ykchoi Jan/23/2018]
# else
       ri=R_sponge**(50*(I-1)/(Iwidth-1))
# endif
       Sponge(i,j)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO

END SUBROUTINE CALCULATE_SPONGE_MAKER

!-------------------------------------------------------------------------------------
!
!    CALCULATE_FRICTION_SPONGE is subroutine for coefficient of 
!      friction type sponge layer
!    
!    HISTORY: 05/01/2014 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_FRICTION_SPONGE(M,N,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,SPONGE)
     USE PARAM
     USE GLOBAL, ONLY : CDsponge
                        
# if defined (PARALLEL)
     USE GLOBAL, ONLY : px,py,npx,npy,&
                        iista, iiend, jjsta, jjend    !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN)::DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY
# endif
     REAL(SP),INTENT(IN) :: &
                          Sponge_west_width,Sponge_east_width,&
                          Sponge_south_width,Sponge_north_width, &
                          R_sponge,A_sponge
     REAL(SP),DIMENSION(M,N),INTENT(OUT)::SPONGE
     REAL(SP)::ri,lim
     INTEGER::I,J,Iwidth
     REAL(SP) :: DXg,DYg
     REAL(SP),DIMENSION(M,N) :: tmp_2d_1,tmp_2d_2
     nop(A_sponge)
     nop(R_sponge)


! note that I used SPONGE represent CD from input, do not mess up 
! with the sponge defined in the direction sponge 

# if defined (CARTESIAN)
     DXg=DX
     DYg=DY
# else
     DXg=DX(1,1)
     DYg=DY(1,1)
# endif

     SPONGE = ZERO
     tmp_2d_2 = ZERO
     tmp_2d_1 = ZERO

! west

     IF(Sponge_west_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_west_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_west_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,M
         lim=0.0_SP
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!        ri=MAX(0.0,REAL(Iwidth-i-npx*Mglob/px))
        ri=MAX(0.0,REAL( Iwidth-i-(iista-1) ))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-I))
# endif
       tmp_2d_1(i,j)=MAX(CDsponge*TANH(ri/10.0),lim)

     ENDDO
     ENDDO
     ENDIF

! east

     IF(Sponge_east_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_east_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_east_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,M
         lim=0.0_SP
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-M+I-(px-npx-1)*Mglob/px))
       ri=MAX(0.0,REAL(Iwidth-M+I-(px-npx-1)*(iiend - iista + 1) ))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-M+I))
# endif
       tmp_2d_2(i,j)=MAX(CDsponge*TANH(ri/10.0),lim)
     ENDDO
     ENDDO
     ENDIF

      DO J=1,N
      DO I=1,M
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
      ENDDO
      ENDDO

! south

     IF(Sponge_south_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_south_width/DY)+Nghost
# else
     Iwidth=INT(Sponge_south_width/DY(1,1))+Nghost
# endif
     DO I=1,M
     DO J=1,N
       IF(SPONGE(I,J)>0.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=0.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-J-npy*Nglob/py))
       ri=MAX(0.0,REAL(Iwidth-J-(jjsta - 1)))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-J))
# endif
       tmp_2d_1(i,j)=MAX(CDsponge*TANH(ri/10.0),lim)
     ENDDO
     ENDDO
     ENDIF


! north

     IF(Sponge_north_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_north_width/DY)+Nghost
# else
     Iwidth=INT(Sponge_north_width/DY(1,1))+Nghost
# endif
     DO I=1,M
     DO J=1,N
       IF(SPONGE(I,J)>0.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=0.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-N+J-(py-npy-1)*Nglob/py))
       ri=MAX(0.0,REAL(Iwidth-N+J-(py-npy-1)*(jjend - jjsta + 1) ))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-N+J))
# endif
       tmp_2d_2(i,j)=MAX(CDsponge*TANH(ri/10.0),lim)
     ENDDO
     ENDDO
     ENDIF

      DO J=1,N
      DO I=1,M
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
      ENDDO
      ENDDO


END SUBROUTINE CALCULATE_FRICTION_SPONGE

!-------------------------------------------------------------------------------------
!
!    CALCULATE_DIFFUSION_SPONGE is subroutine for coefficient of 
!      diffusion type sponge layer
!    
!    HISTORY: 05/01/2014 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_DIFFUSION_SPONGE(M,N,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,SPONGE)
     USE PARAM
     USE GLOBAL, ONLY : Csp
                        
# if defined (PARALLEL)
     USE GLOBAL, ONLY : px,py,npx,npy,&
                        iista, iiend, jjsta, jjend    !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN)::DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY
# endif
     REAL(SP),INTENT(IN) :: &
                          Sponge_west_width,Sponge_east_width,&
                          Sponge_south_width,Sponge_north_width, &
                          R_sponge,A_sponge
     REAL(SP),DIMENSION(M,N),INTENT(INOUT)::SPONGE
     REAL(SP)::ri,lim
     INTEGER::I,J,Iwidth
     REAL(SP) :: DXg,DYg
     REAL(SP),DIMENSION(M,N) :: tmp_2d_1,tmp_2d_2
     nop(A_sponge)
     nop(R_sponge)


     tmp_2d_2 = ZERO
     tmp_2d_1 = ZERO

# if defined (CARTESIAN)
     DXg=DX
     DYg=DY
# else
     DXg=DX(1,1)
     DYg=DY(1,1)
# endif

! west

     IF(Sponge_west_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_west_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_west_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,M
       IF(SPONGE(I,J)>0.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=0.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-i-npx*Mglob/px))
       ri=MAX(0.0,REAL(Iwidth-i-(iista - 1)))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-I))
# endif
       tmp_2d_1(i,j)=MAX(Csp*TANH(ri/10.0),lim)

     ENDDO
     ENDDO
     ENDIF

! east

     IF(Sponge_east_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_east_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_east_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,M
       IF(SPONGE(I,J)>0.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=0.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-M+I-(px-npx-1)*Mglob/px))
       ri=MAX(0.0,REAL(Iwidth-M+I-(px-npx-1)*(iiend - iista + 1) ))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-M+I))
# endif
       tmp_2d_2(i,j)=MAX(Csp*TANH(ri/10.0),lim)
     ENDDO
     ENDDO
     ENDIF

      DO J=1,N
      DO I=1,M
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
      ENDDO
      ENDDO

! south

     IF(Sponge_south_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_south_width/DY)+Nghost
# else
     Iwidth=INT(Sponge_south_width/DY(1,1))+Nghost
# endif
     DO I=1,M
     DO J=1,N
       IF(SPONGE(I,J)>0.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=0.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-J-npy*Nglob/py))
       ri=MAX(0.0,REAL(Iwidth-J-(jjsta - 1)))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-J))
# endif
       tmp_2d_1(i,j)=MAX(Csp*TANH(ri/10.0),lim)
     ENDDO
     ENDDO
     ENDIF


! north

     IF(Sponge_north_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_north_width/DY)+Nghost
# else
     Iwidth=INT(Sponge_north_width/DY(1,1))+Nghost
# endif
     DO I=1,M
     DO J=1,N
       IF(SPONGE(I,J)>0.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=0.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=MAX(0.0,REAL(Iwidth-N+J-(py-npy-1)*Nglob/py))
       ri=MAX(0.0,REAL(Iwidth-N+J-(py-npy-1)*(jjend - jjsta + 1) ))
!---ykchoi Jan/23/2018]
# else
       ri=MAX(0.0,REAL(Iwidth-N+J))
# endif
       tmp_2d_2(i,j)=MAX(Csp*TANH(ri/10.0),lim)
     ENDDO
     ENDDO
     ENDIF

      DO J=1,N
      DO I=1,M
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
      ENDDO
      ENDDO

END SUBROUTINE CALCULATE_DIFFUSION_SPONGE

!-------------------------------------------------------------------------------------
!
!    CALCULATE_SPONGE is subroutine of sponge layer to get coefficient
!
!    HISTORY: 
!     10/27/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------

SUBROUTINE CALCULATE_SPONGE(M,N,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,SPONGE)

     USE PARAM
# if defined (PARALLEL)
     USE GLOBAL, ONLY : px,py,npx,npy,&
                        iista, iiend, jjsta, jjend    !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN)::DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY
# endif
     REAL(SP),INTENT(IN) :: &
                          Sponge_west_width,Sponge_east_width,&
                          Sponge_south_width,Sponge_north_width, &
                          R_sponge,A_sponge

     REAL(SP),DIMENSION(M,N),INTENT(INOUT)::SPONGE
     REAL(SP)::ri,lim
     INTEGER::I,J,Iwidth
     REAL(SP) :: DXg,DYg
     REAL(SP),DIMENSION(M,N) :: tmp_2d_1,tmp_2d_2 

     tmp_2d_2 = ZERO
     tmp_2d_1 = ZERO

# if defined (CARTESIAN)
     DXg=DX
     DYg=DY
# else
     DXg=DX(1,1)
     DYg=DY(1,1)
# endif

! west

     IF(Sponge_west_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_west_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_west_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,M
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri = R_Sponge**(50*(i+npx*Mglob/px-1)/(Iwidth-1))
       ri = R_Sponge**(50*(i+(iista-1)-1)/(Iwidth-1))
!---ykchoi Jan/23/2018]
# else
       ri=R_sponge**(50*(I-1)/(Iwidth-1))
# endif
       tmp_2d_1(i,j)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF

! east

     IF(Sponge_east_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_east_width/DX)+Nghost
# else
     Iwidth=INT(Sponge_east_width/DX(1,1))+Nghost
# endif
     DO J=1,N
     DO I=1,M
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri = R_Sponge**(50*(M-i+(px-npx-1)*Mglob/px)/(Iwidth-1))
       ri = R_Sponge**(50*(M-i+(px-npx-1)*(iiend - iista + 1))/(Iwidth-1))
!---ykchoi Jan/23/2018]
# else
       ri=R_sponge**(50*(M-I)/(Iwidth-1))
# endif
       tmp_2d_2(i,j)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF

      DO J=1,N
      DO I=1,M
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
        IF(Sponge(I,J)<A_Sponge**(R_Sponge**50))Sponge(I,J)=1.0_SP
      ENDDO
      ENDDO

! south

     IF(Sponge_south_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_south_width/DY)+Nghost
# else
     Iwidth=INT(Sponge_south_width/DY(1,1))+Nghost
# endif
     DO I=1,M
     DO J=1,N
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri=R_sponge**(50*(J+npy*Nglob/py-1)/(Iwidth-1))
       ri=R_sponge**(50*(J+ (jjsta - 1) -1)/(Iwidth-1))
!---ykchoi Jan/23/2018]
# else
       ri=R_sponge**(50*(J-1)/(Iwidth-1))
# endif
       tmp_2d_1(i,j)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF


! north

     IF(Sponge_north_width>ZERO)THEN
# if defined (CARTESIAN)
     Iwidth=INT(Sponge_north_width/DY)+Nghost
# else
     Iwidth=INT(Sponge_north_width/DY(1,1))+Nghost
# endif
     DO I=1,M
     DO J=1,N
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
# if defined(PARALLEL)
![---ykchoi Jan/23/2018
!       ri = R_Sponge**(50*(N-J+(py-npy-1)*Nglob/py)/(Iwidth-1))
       ri = R_Sponge**(50*(N-J+(py-npy-1)*(jjend - jjsta + 1) )/(Iwidth-1))
!---ykchoi Jan/23/2018]
# else
       ri=R_sponge**(50*(N-J)/(Iwidth-1))
# endif
       tmp_2d_2(i,j)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF

      DO J=1,N
      DO I=1,M
        IF(tmp_2d_1(I,J)>tmp_2d_2(I,J)) THEN
          Sponge(I,J)=tmp_2d_1(I,J)
        ELSE
          Sponge(I,J)=tmp_2d_2(I,J)
        ENDIF
        IF(Sponge(I,J)<A_Sponge**(R_Sponge**50))Sponge(I,J)=1.0_SP
      ENDDO
      ENDDO

END SUBROUTINE CALCULATE_SPONGE


!-------------------------------------------------------------------------------------
!
!    CALCULATE_CD_BREAKWATER is subroutine for coefficient of 
!      friction type breakwater
!    
!    HISTORY: 04/27/2017 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_CD_BREAKWATER
     USE PARAM
     USE GLOBAL, ONLY : CD_breakwater,BreakWaterWidth, &
                        Mglob,Nglob,Nghost,DX,DY,&
                        BreakWaterAbsorbCoef
# if defined (PARALLEL)
     USE GLOBAL, ONLY : Mloc,Nloc,px,py,myid,ier,NumberProcessor, &
                        iista, iiend, jjsta, jjend    !ykchoi (28/Jan/2018)   
# endif                 
     IMPLICIT NONE
# if defined (PARALLEL)
![-------ykchoi (28/Jan/2018)      
      !INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
      !REAL(SP),DIMENSION(NumberProcessor) :: xx
      !INTEGER :: l
      INTEGER :: Nista, Niend, Njsta, Njend
      INTEGER, ALLOCATABLE :: Nistas(:), Niends(:), Njstas(:), Njends(:)
      INTEGER :: irank, lenx, leny, lenxy, ireq
      INTEGER :: istanum, iendnum, jstanum, jendnum
      INTEGER :: istatus(mpi_status_size)
      REAL(SP), ALLOCATABLE :: xx(:,:)
!-------ykchoi (28/Jan/2018)]
# endif

     REAL(SP)::ri
     INTEGER::I,J,Iwidth,Jwidth,Ib,Jb
     REAL(SP) :: DXg,DYg,tmp_2d
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: tmp_glob
 
     ALLOCATE (tmp_glob(Mglob+2*Nghost,Nglob+2*Nghost))

# if defined (PARALLEL)
     if (myid.eq.0) then
# endif

# if defined (CARTESIAN)
     DXg=DX
     DYg=DY
# else
     DXg=DX(1,1)
     DYg=DY(1,1)
# endif

     CD_breakwater = ZERO

! global

    tmp_glob = ZERO
    DO J=1,Nglob+2*Nghost
    DO I=1,Mglob+2*Nghost
      IF(BreakWaterWidth(I,J).GT.ZERO)THEN
        Iwidth=INT(BreakWaterWidth(I,J)/DXg)
        Jwidth=INT(BreakWaterWidth(I,J)/DYg)
    ! fixed bug from gangfeng 02/26/2021, should add 2*Nghost
        DO Jb=MAX(1,J-Jwidth),Min(Nglob+2*Nghost,J+Jwidth)
        DO Ib=MAX(1,I-Iwidth),Min(Mglob+2*Nghost,I+Iwidth)
          ri=SQRT(((Ib-I)*DXg)**2+((Jb-J)*DYg)**2)
          IF(ri.LE.BreakWaterWidth(I,J))THEN
           tmp_2d=BreakWaterAbsorbCoef* &
             TANH( (BreakWaterWidth(I,J)-ri) / (0.3_SP*BreakWaterWidth(I,J)) )

!if (J==20+Nghost.AND.I==825+Nghost)then
!print*,Ib,Jb,I,J,ri,tmp_2d,tmp_glob(I,J)
!endif
           IF(tmp_2d.GT.tmp_glob(Ib,Jb)) tmp_glob(Ib,Jb)=tmp_2d

          ENDIF
        ENDDO
        ENDDO
      ENDIF
    ENDDO
    ENDDO

# if defined (PARALLEL)
     endif ! end myid
# endif

! distribute to CD_breakwater 

# if defined (PARALLEL)
![-------ykchoi (28/Jan/2018)
!     call MPI_Gather(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,&
!          0,MPI_COMM_WORLD,ier)
!     call MPI_Gather(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,&
!          0,MPI_COMM_WORLD,ier)

!     do i=1,Mloc
!     do j=1,Nloc
!        if (myid.eq.0) then
!           do l=1,px*py
!              xx(l) = tmp_glob(i+npxs(l)*(Iend-Ibeg+1),&
!                   j+npys(l)*(Jend-Jbeg+1))
!           enddo
!        endif
!        call MPI_Scatter(xx,1,MPI_SP,&
!             CD_breakwater(i,j),1,MPI_SP,0,MPI_COMM_WORLD,ier)
!     enddo
!     enddo

     Nista = iista + Nghost;
     Niend = iiend + Nghost;     
     Njsta = jjsta + Nghost;     
     Njend = jjend + Nghost;

     allocate( Nistas(NumberProcessor), Niends(NumberProcessor),   &
               Njstas(NumberProcessor), Njends(NumberProcessor) )
          
     call MPI_Gather( Nista, 1, MPI_INTEGER, Nistas, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( Niend, 1, MPI_INTEGER, Niends, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( Njsta, 1, MPI_INTEGER, Njstas, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( Njend, 1, MPI_INTEGER, Njends, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )

     if( myid == 0 )then
           CD_breakwater = tmp_glob( 1:Mloc, 1:Nloc )
     endif          

     do irank=1, px*py-1
          if( myid == 0 ) then
          istanum = Nistas(irank+1) - Nghost
          iendnum = Niends(irank+1) + Nghost
          jstanum = Njstas(irank+1) - Nghost
          jendnum = Njends(irank+1) + Nghost

          lenx = iendnum - istanum + 1
          leny = jendnum - jstanum + 1
          lenxy = lenx*leny
          allocate( xx(lenx, leny) )
          
          xx = tmp_glob( istanum:iendnum, jstanum:jendnum )
          call mpi_isend( xx, lenxy, mpi_sp, irank, 1, mpi_comm_world, ireq, ier )
          call mpi_wait( ireq, istatus, ier )
          deallocate( xx )
          
        elseif( myid == irank ) then

          lenx = Niend-Nista+1+2*Nghost
          leny = Njend-Njsta+1+2*Nghost
          lenxy = lenx*leny
          
          call mpi_irecv( CD_breakwater, lenxy, mpi_sp, 0, 1, mpi_comm_world, ireq, ier )
          call mpi_wait( ireq, istatus, ier )
          
        endif
     enddo

     deallocate( Nistas, Niends, Njstas, Njends )
!-------ykchoi (28/Jan/2018)]

# else
     CD_breakwater = tmp_glob
# endif
    
    DEALLOCATE (tmp_glob)

END SUBROUTINE CALCULATE_CD_BREAKWATER



