!------------------------------------------------------------------------------------
!
!      FILE etauv_solver.F
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
!    ESTIMATE_HUV is subroutine to calculate eta, ubar and vbar
!      using 3rd-order LK scheme 
!
!    HISTORY: 
!      05/12/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE ESTIMATE_HUV(ISTEP)
     USE GLOBAL
# if defined (VESSEL)
     USE VESSEL_MODULE
# endif
     IMPLICIT NONE
     INTEGER,INTENT(IN)::ISTEP
     REAL(SP),PARAMETER::n_left=-1.0_SP,n_right=1.0_SP,n_bottom=-1.0_SP,n_top=1.0_SP
     REAL(SP)::F_left,F_right,F_bottom,F_top
     REAL(SP),DIMENSION(Ibeg:Iend,Jbeg:Jend)::R1,R2,R3
! now work for spherical # if defined (CARTESIAN)
     REAL(SP)::xmk,ymk
! now work for spherical # endif
     REAL(SP)::DXg,DYg
     INTEGER :: I,J

! MUSCL-Hancock, Zhou et al., p. 7

# if defined (CARTESIAN)
     DXg=DX
     DYg=DY
# else
! only for wavemaker
     DXg=DX(1,1)
     DYg=DY(1,1)
# endif

! solve eta
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=P(I,J)
      F_right=P(I+1,J)
      F_bottom=Q(I,J)
      F_top=Q(I,J+1)
! now work for spherical # if defined (CARTESIAN)
        IF(WaveMaker(1:6)=='WK_IRR'.OR.WaveMaker(1:6)=='TMA_1D'  &
           .OR.WaveMaker(1:6)=='JON_1D'.OR.WaveMaker(1:6)=='JON_2D' &
           .OR.WAVEMAKER(1:10)=='WK_NEW_IRR'&             ! Salatin et al. 2021
           .OR.WAVEMAKER(1:13)=='WK_NEW_DATA2D')THEN      ! Salatin et al. 2021
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
!            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg
            xmk=(I-Ibeg)*DXg + (iista-1)*DXg
            ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
!---ykchoi Jan/23/2018]
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN

! WaveMaker_Mass was done in sources.F 02/27/2017

!          WK_Source=ZERO
!          DO kf=1,Nfreq
!           WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm(I,J,kf) &
!                       *COS(OMGN_IR(KF)*TIME) &
!                       +Sm(I,J,kf)*SIN(OMGN_IR(KF)*TIME))
!          ENDDO

          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
        ! wavemaker
                 +WaveMaker_Mass(I,J)
!                +WK_Source      
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF
       ELSEIF(WAVEMAKER(1:6)=='WK_REG')THEN
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
!            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg
            xmk=(I-Ibeg)*DXg + (iista-1)*DXg
            ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
!---ykchoi Jan/23/2018]
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
        ! wavemaker 
                 +WaveMaker_Mass(I,J)    
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF
       ELSEIF(WAVEMAKER(1:7)=='WK_TIME')THEN
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
!            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg
            xmk=(I-Ibeg)*DXg + (iista-1)*DXg
            ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
!---ykchoi Jan/23/2018]
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
                +WaveMaker_Mass(I,J)      
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF   

       ELSEIF(WAVEMAKER(1:9)=='WK_DATA2D')THEN
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
!            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg
            xmk=(I-Ibeg)*DXg + (iista-1)*DXg
            ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
!---ykchoi Jan/23/2018]
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN!
          
          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
                +WaveMaker_Mass(I,J)      
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF 
   
      ELSE ! no wk_wavemaker, theres bug in version 1.1 Dxg,Dyg should be 
           ! replaced by Dxg() and Dy()
# if defined (CARTESIAN)
        R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
# if defined (VESSEL)
        R1(I,J)=R1(I,J) +VesselFluxGradient(I,J)
# endif
  ! end vessel

# else
        R1(I,J)=-1.0_SP/DX(I,J)*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DY(I,J)*(F_top*n_top+F_bottom*n_bottom)
# endif
      ENDIF
# if defined (CARTESIAN)
! do nothing
# else

# if defined (ZALPHA)
! extra terms for z_alpha: 1/R tan theta (HV+hV4)
        R1(I,J)=R1(I,J)+1.0_SP/R_earth*TAN(Lat_theta(I,J))*(HV(I,J)+ &
                Depth(I,J)*V4(I,J)*Gamma1)
# else
! extra terms for depth-averaged u: 1/R tan theta HV
        R1(I,J)=R1(I,J)+1.0_SP/R_earth*TAN(Lat_theta(I,J))*HV(I,J)
# endif

# endif

# if defined (ITERATION)
      Eta(I,J)=ALPHA(ISTEP)*Eta0(I,J)+BETA(ISTEP)*(EtaOld(I,J)+DT*R1(I,J))
# else
      Eta(I,J)=ALPHA(ISTEP)*Eta0(I,J)+BETA(ISTEP)*(Eta(I,J)+DT*R1(I,J))

! eta_limiter is used for the case as wave touches the seabed within wavemaker
      IF(ETA_LIMITER)THEN
        IF(Eta(I,J)<TroughLimit)Eta(I,J)=TroughLimit
        IF(Eta(I,J)>CrestLimit)Eta(I,J)=CrestLimit
      ENDIF ! end eta_limiter

# endif
     ENDDO
     ENDDO

! solve ubar
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=Fx(I,J)
      F_right=Fx(I+1,J)
      F_bottom=Fy(I,J)
      F_top=Fy(I,J+1)
# if defined (CARTESIAN)
      R2(I,J)=-1.0_SP/DX*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY*(F_top*n_top+F_bottom*n_bottom) &
                        +SourceX(I,J)
# else
      R2(I,J)=-1.0_SP/DX(I,J)*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY(I,J)*(F_top*n_top+F_bottom*n_bottom) &
                        +SourceX(I,J)
# endif
# if defined (ITERATION)
      Ubar(I,J)=ALPHA(ISTEP)*Ubar0(I,J)+BETA(ISTEP)*(UbarOld(I,J)+DT*R2(I,J))
# else
      Ubar(I,J)=ALPHA(ISTEP)*Ubar0(I,J)+BETA(ISTEP)*(Ubar(I,J)+DT*R2(I,J))
# endif
     ENDDO
     ENDDO

! solve vbar
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=Gx(I,J)
      F_right=Gx(I+1,J)
      F_bottom=Gy(I,J)
      F_top=Gy(I,J+1)
# if defined (CARTESIAN)
      R3(I,J)=-1.0_SP/DX*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY*(F_top*n_top+F_bottom*n_bottom) &
                       +SourceY(I,J)
# else
      R3(I,J)=-1.0_SP/DX(I,J)*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY(I,J)*(F_top*n_top+F_bottom*n_bottom) &
                       +SourceY(I,J)
# endif
# if defined (ITERATION)
      Vbar(I,J)=ALPHA(ISTEP)*Vbar0(I,J)+BETA(ISTEP)*(VbarOld(I,J)+DT*R3(I,J))
# else
      Vbar(I,J)=ALPHA(ISTEP)*Vbar0(I,J)+BETA(ISTEP)*(Vbar(I,J)+DT*R3(I,J))
# endif
     ENDDO
     ENDDO

     CALL GET_Eta_U_V_HU_HV

END SUBROUTINE ESTIMATE_HUV

!-------------------------------------------------------------------------------------
!
!    GET_Eta_U_V_HU_HV is subroutine to obtain Eta, u,v,hu,hv
!
!  HISTORY: 
!       09/17/2010 Fengyan Shi
!       10/14/2012 Fengyan Shi, added nesting bc
!       08/06/2015 Young-Kwang Choi, modified U0 V0 shift 
!       01/27/2016 Fengyan Shi, made serial/parallel codes consistent
!                               added parallel code of periodic bc in y
!
!-------------------------------------------------------------------------------------
SUBROUTINE GET_Eta_U_V_HU_HV
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::Fr,Utotal,Utheta,dep,depl,depr
     REAL(SP),DIMENSION(Mloc,Nloc) :: myA,myC,myD,myF
# if defined (CARTESIAN)
     REAL(SP)::reta,retal,retar
# endif
     INTEGER :: I,J
     REAL(SP) :: tmp1,tmp2,tmp3,tmp4

! calculate etar, u and vetar, HU, HV
     H=Eta*Gamma3+Depth


!   tridiagonal coefficient
! x direction

! shift U and V
!     U0=U    !ykchoi (15. 08. 06.) 
!     V0=V    !ykchoi
! ykchoi : U0, V0 are moved to the above part of "Do istage=1,3"

   IF(DISPERSION)THEN

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       dep=Max(Depth(I,J),MinDepthFrc)
       depl=Max(Depth(I-1,J),MinDepthFrc)
       depr=Max(Depth(I+1,J),MinDepthFrc)

# if defined (CARTESIAN)
       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DX/DX*dep*dep + b2/DX/DX*depl*dep)
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DX/DX*dep*dep-2.0_SP*b2/DX/DX*dep*dep)
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DX/DX*dep*dep + b2/DX/DX*dep*depr)
       tmp4=Ubar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
            + Gamma1*MASK9(I,J)*( -b1/2.0_SP*dep*dep*Vxy(I,J)-b2*dep*DVxy(I,J))
# else
       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DX(I,J)/DX(I,J)*dep*dep + b2/DX(I,J)/DX(I,J)*depl*dep)
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DX(I,J)/DX(I,J)*dep*dep-2.0_SP*b2/DX(I,J)/DX(I,J)*dep*dep)
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DX(I,J)/DX(I,J)*dep*dep + b2/DX(I,J)/DX(I,J)*dep*depr)
       tmp4=Ubar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
            + Gamma1*MASK9(I,J)*( -b1/2.0_SP*dep*dep*Vxy(I,J)-b2*dep*DVxy(I,J)) 
! remember to document this part for spherical version 
# endif

! I added coupling condition 10/14/2012
!  added west_bc wavemaker 09/12/2017

# if defined (COUPLING)
     IF(IN_DOMAIN_WEST)THEN
       IF(I.eq.Ibeg)THEN
         tmp4=tmp4-tmp1*U(I-1,J)
       ENDIF
     ELSEIF(IN_DOMAIN_EAST)THEN
       IF(I.eq.Iend)THEN
         tmp4=tmp4-tmp3*U(I+1,J)
       ENDIF
     ENDIF
# endif

!  left_bc wavemaker
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif

    IF (WaveMaker(1:11)=='LEFT_BC_IRR')THEN
       IF(I.eq.Ibeg)THEN
         tmp4=tmp4-tmp1*U(I-1,J)
       ENDIF
    ENDIF

# if defined (PARALLEL)
    endif  
# endif


       IF(tmp2.NE.0.0_SP.OR.MASK(I,J).GT.0)THEN
          myA(I,J)=tmp1/tmp2
          myC(I,J)=tmp3/tmp2
          myD(I,J)=tmp4/tmp2
       ELSE
          myA(I,J)=ZERO
          myC(I,J)=ZERO
          myD(I,J)=ZERO
       ENDIF
     ENDDO
     ENDDO

# if defined (PARALLEL)
     call TRIDx(myA,myC,myD,myF)
     U(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)
# else 
     call TRIDx_ser(myA,myC,myD,myF)
     U(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)
# endif

! y direction
# if defined (PARALLEL)
     myA=ZERO
     myC=ZERO
     myD=ZERO
# endif
     DO I=Ibeg,Iend
     DO J=Jbeg,Jend
       dep=Max(Depth(I,J),MinDepthFrc)
       depl=Max(Depth(I,J-1),MinDepthFrc)
       depr=Max(Depth(I,J+1),MinDepthFrc)

# if defined (CARTESIAN)
     IF(DISP_TIME_LEFT)THEN
       reta=Eta(I,J)
       retal=Eta(I,J-1)
       retar=Eta(I,J+1)
       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY/DY*dep*dep + b2/DY/DY*depl*dep) &
             -Gamma2*MASK9(I,J)*((reta+retal)*depl/2.0_SP/DY/DY+(retal+reta)*(retal+reta)/8.0_SP/DY/DY)
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DY/DY*dep*dep-2.0_SP*b2/DY/DY*dep*dep) &
             +Gamma2*MASK9(I,J)*((retar+retal+2.0_SP*reta)/2.0_SP/DY/DY &
                        +(retar*retar+2.0_SP*reta*reta &
                          +2.0_SP*reta*retar+2.0_SP*retal*reta+retal*retal)/8.0_SP/DY/DY)
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY/DY*dep*dep + b2/DY/DY*dep*depr) &
             -Gamma2*MASK9(I,J)*((reta+retar)*depr/2.0_SP/DY/DY+(retar+reta)*(retar+reta)/8.0_SP/DY/DY)
       tmp4=Vbar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
             + Gamma1*MASK9(I,J)*(-b1/2.0_SP*dep*dep*Uxy(I,J)-b2*dep*DUxy(I,J)) &
            + Gamma2*MASK9(I,J)*(reta*reta/2.0_SP*Uxy(I,J)+reta*DUxy(I,J) &
                             + reta*ETAy(I,J)*Ux(I,J) + ETAy(I,J)*DUx(I,J) )
     ELSE
       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY/DY*dep*dep + b2/DY/DY*depl*dep) 
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DY/DY*dep*dep-2.0_SP*b2/DY/DY*dep*dep) 
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY/DY*dep*dep + b2/DY/DY*dep*depr)
       tmp4=Vbar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
             + Gamma1*MASK9(I,J)*(-b1/2.0_SP*dep*dep*Uxy(I,J)-b2*dep*DUxy(I,J))
     ENDIF  
# else
       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY(I,J)/DY(I,J)*dep*dep + b2/DY(I,J)/DY(I,J)*depl*dep) 
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DY(I,J)/DY(I,J)*dep*dep-2.0_SP*b2/DY(I,J)/DY(I,J)*dep*dep) 
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY(I,J)/DY(I,J)*dep*dep + b2/DY(I,J)/DY(I,J)*dep*depr)
       tmp4=Vbar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
             + Gamma1*MASK9(I,J)*(-b1/2.0_SP*dep*dep*Uxy(I,J)-b2*dep*DUxy(I,J)) 
! remember to document this part for spherical version
# endif

# if defined (COUPLING)
     IF(IN_DOMAIN_SOUTH)THEN
       IF(J.eq.Jbeg)THEN
         tmp4=tmp4-tmp1*V(I,J-1)
       ENDIF
     ELSEIF(IN_DOMAIN_SOUTH)THEN
       IF(J.eq.Jend)THEN
         tmp4=tmp4-tmp3*V(I,J+1)
       ENDIF
     ENDIF
# endif
  
       IF(tmp2.NE.0.0_SP.OR.MASK(I,J).GT.0)THEN
         myA(I,J)=tmp1/tmp2
         myC(I,J)=tmp3/tmp2
         myD(I,J)=tmp4/tmp2
       ELSE
         myA(I,J)=ZERO
         myC(I,J)=ZERO
         myD(I,J)=ZERO
       ENDIF
     ENDDO
     ENDDO ! end I

     IF(PERIODIC) THEN

!  Sherman-Morrison algorithm 01/27/2016 fyshi

      CALL TRIDy_periodic (myA,myC,myD,myF)  
       V(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)
     ELSE ! no periodic
# if defined (PARALLEL)
       CALL TRIDy(myA,myC,myD,myF)
       V(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)
# else
       CALL TRIDy_ser(myA,myC,myD,myF)
       V(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)
# endif

     ENDIF ! end if periodic

   ELSE  ! if no dispersion
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend  
        U(I,J)=Ubar(I,J)/Max(H(I,J),MinDepthFrc)
        V(I,J)=Vbar(I,J)/Max(H(I,J),MinDepthFrc)
     ENDDO
     ENDDO   

   ENDIF  ! end dispersion

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend   
       IF(MASK(I,J)<1)THEN
        Ubar(I,J)=ZERO
        Vbar(I,J)=ZERO
        U(I,J)=ZERO
        V(I,J)=ZERO
        HU(I,J)=ZERO
        HV(I,J)=ZERO
       ELSE
        HU(I,J)=Max(H(I,J),MinDepthFrc)*U(I,J)
        HV(I,J)=Max(H(I,J),MinDepthFrc)*V(I,J)
! apply Froude cap
        Utotal=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
        Fr=SQRT(GRAV*Max(H(I,J),MinDepthFrc))
        IF(Utotal/Fr.gt.FroudeCap)THEN
          Utheta=ATAN2(V(I,J),U(I,J))
          U(I,J)=FroudeCap*Fr*COS(Utheta)
          V(I,J)=FroudeCap*Fr*SIN(Utheta)
          HU(I,J)=U(I,J)*Max(H(I,J),MinDepthFrc)
          HV(I,J)=V(I,J)*Max(H(I,J),MinDepthFrc)
        ENDIF
! end Froude cap
       ENDIF
     ENDDO
     ENDDO

!------------ykchoi 07/26/2016

     IF( .NOT. DISPERSION) THEN

        Ubar = HU
        Vbar = HV

     ENDIF

END SUBROUTINE GET_Eta_U_V_HU_HV
