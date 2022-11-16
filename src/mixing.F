!------------------------------------------------------------------------------------
!
!      FILE mixing.F
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
!-------------------------------------------------------------------------------------
!
!    MIXING_STUFF is subroutine to calculate mixing related, time-averaged properties
!    mean eta is also calculated.
!    
!    HISTORY: 05/02/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE MIXING_STUFF
     USE GLOBAL
     IMPLICIT NONE

! calculate mean for smagorinsky s mixing and wave height
	!ykchoi (be careful)
	!I think Umean, Vmean is not using in other routine.
      IF( time >= STEADY_TIME )THEN    
      CALL CALCULATE_MEAN(T_INTV_mean,T_sum,DT,DXg,DYg,Mloc,Nloc,U,V,Wsurf,ETA,ETA0,&
           Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,WaveHeightRMS, &
           WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
           HrmsSum,HavgSum, &
	     !ykchoi
	     ETA2sum, ETA2mean, SigWaveHeight, &
                 UUsum,UUmean,UVsum,UVmean,VVsum,VVmean, &
                 WWsum,WWmean,FRCXsum,FRCXmean,FRCYsum,FRCYmean, &
                 DxSxx,DySxy,DySyy,DxSxy,PgrdX,PgrdY,DxUUH,DyUVH,DyVVH,DxUVH,&
                 Cd,Depth,U_davg,V_davg,U_davg_sum,V_davg_sum, &
                  U_davg_mean,V_davg_mean,P_center,Q_center,P_sum,Q_sum, &
                  P_mean,Q_mean,MinDepthFrc)
      ENDIF    !ykchoi

END SUBROUTINE MIXING_STUFF


!-------------------------------------------------------------------------------------
!
!    CALCULATE_MEAN is subroutine to calculate mean u v required by 
!      smagorinsky mixing and wave height
!      mean eta is also calculated.
!    
!    HISTORY: 
!      05/02/2011 Fengyan Shi
!                 Young-Kwang Choi added some time-averaging stuff
!
!-------------------------------------------------------------------------------------
      SUBROUTINE CALCULATE_MEAN(T_INTV_mean,T_sum,DT,DXg,DYg,M,N,U,V,Wsurf,ETA,ETA0,&
                Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,&
                WaveHeightRMS, &
                WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
                HrmsSum,HavgSum, &
	          !ykchoi
			  ETA2sum,ETA2mean,SigWaveHeight, &
                  ! for radiation stress
                 UUsum,UUmean,UVsum,UVmean,VVsum,VVmean, &
                 WWsum,WWmean,FRCXsum,FRCXmean,FRCYsum,FRCYmean, &
                 DxSxx,DySxy,DySyy,DxSxy,PgrdX,PgrdY,DxUUH,DyUVH,DyVVH,DxUVH, &
                 Cd,Depth,U_davg,V_davg,U_davg_sum,V_davg_sum, &
                  U_davg_mean,V_davg_mean,P_center,Q_center,P_sum,Q_sum, &
                  P_mean,Q_mean,MinDepthFrc)
! calculate mean for smagorinsky s mixing and wave height

      USE PARAM
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
      REAL(SP),DIMENSION(M,N),INTENT(IN)::U,V,ETA,ETA0,Wsurf,Cd,Depth
      REAL(SP),INTENT(IN) :: T_INTV_mean,DT,DXg,DYg,MinDepthFrc
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Umean,Vmean
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: WaveHeightRMS,WaveHeightAve
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) ::ETAmean
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Usum,Vsum,ETAsum
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: HrmsSum,HavgSum
      REAL(SP),INTENT(OUT) :: T_sum
      REAL(SP)::Tmpe,Tmp_0
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Emax,Emin
      INTEGER,DIMENSION(M,N),INTENT(INOUT) :: Num_Zero_Up
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: UUmean,UUsum,UVmean,UVsum,VVmean,VVsum
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: WWmean,WWsum
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: FRCXmean,FRCXsum
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: FRCYmean,FRCYsum
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: DxSxx,DySxy,DySyy,DxSxy, &
                                            PgrdX,PgrdY,DxUUH,DyUVH, &
                                            DyVVH,DxUVH	
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: U_davg_sum,V_davg_sum, &
                                             U_davg_mean,V_davg_mean,&
                                             P_sum,Q_sum, &
                                             P_mean,Q_mean
      REAL(SP),DIMENSION(M,N),INTENT(IN) :: U_davg,V_davg,P_center,Q_center

	!ykchoi
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: ETA2sum,ETA2mean,SigWaveHeight
      
      T_sum=T_sum+DT
      IF(T_sum.GE.T_INTV_mean)THEN

	ETA2sum = (Eta-ETAmean)*(Eta-ETAmean)*DT + ETA2sum   !ykchoi
	ETA2mean = ETA2sum/T_sum

        Usum=U*DT+Usum
        Vsum=V*DT+Vsum
        ETAsum=ETA*DT+ETAsum
        Umean=Usum/T_sum
        Vmean=Vsum/T_sum
        ETAmean=ETAsum/T_sum

        U_davg_sum=U_davg*DT+U_davg_sum
        V_davg_sum=V_davg*DT+V_davg_sum
        P_sum=P_center*DT+P_sum
        Q_sum=Q_center*DT+Q_sum
        UUsum = (U_davg-U_davg_mean)*(U_davg-U_davg_mean)*(Eta-ETAmean+Depth)*DT + UUsum
        UVsum = (U_davg-U_davg_mean)*(V_davg-V_davg_mean)*(Eta-ETAmean+Depth)*DT + UVsum
        VVsum = (V_davg-V_davg_mean)*(V_davg-V_davg_mean)*(Eta-ETAmean+Depth)*DT + VVsum
        WWsum = 0.25*Wsurf*Wsurf*(Eta-ETAmean+Depth)*DT + WWsum
        FRCXsum = U*SQRT(U*U+V*V)*DT + FRCXsum
        FRCYsum = V*SQRT(U*U+V*V)*DT + FRCYsum

        U_davg_mean = U_davg_sum/T_sum
        V_davg_mean = V_davg_sum/T_sum
        P_mean = P_sum/T_sum
        Q_mean = Q_sum/T_sum
        UUmean = UUsum/T_sum
        UVmean = UVsum/T_sum
        VVmean = VVsum/T_sum
        WWmean = WWsum/T_sum
        FRCXmean = FRCXsum/T_sum
        FRCYmean = FRCYsum/T_sum

! radiation stresses,
      DO J=Jbeg,Jend
      DO I=Ibeg,Iend  
        DxSxx(I,J) = (UUmean(I+1,J)-UUmean(I-1,J))/2.0/DXg &
                    -(WWmean(I+1,J)-WWmean(I-1,J))/2.0/DXg &
                    +0.5*9.80*(ETA2mean(I+1,J)-ETA2mean(I-1,J))/2.0/DXg
        DySxy(I,J) = (UVmean(I,J+1)-UVmean(I,J-1))/2.0/DYg
        DySyy(I,J) = (VVmean(I,J+1)-VVmean(I,J-1))/2.0/DYg &
                    -(WWmean(I,J+1)-WWmean(I,J-1))/2.0/DYg &
                    +0.5*9.80*(ETA2mean(I,J+1)-ETA2mean(I,J-1))/2.0/DYg
        DxSxy(I,J) = (UVmean(I+1,J)-UVmean(I-1,J))/2.0/DXg
        PgrdX(I,J) = 9.80*(Depth(I,J)+ETAmean(I,J))* &
                     (ETAmean(I+1,J)-ETAmean(I-1,J))/2.0/DXg
        PgrdY(I,J) = 9.80*(Depth(I,J)+ETAmean(I,J))* &
                     (ETAmean(I,J+1)-ETAmean(I,J-1))/2.0/DYg
        DxUUH(I,J) = (P_mean(I+1,J)*P_mean(I+1,J)/ &
                     Max(Depth(I+1,J)+ETAmean(I+1,J),MinDepthFrc) - &
                     P_mean(I-1,J)*P_mean(I-1,J)/ &
                     Max(Depth(I-1,J)+ETAmean(I-1,J),MinDepthFrc))/2.0/DXg
        DyUVH(I,J) = (P_mean(I,J+1)*Q_mean(I,J+1)/ &
                     Max(Depth(I,J+1)+ETAmean(I,J+1),MinDepthFrc) - &
                     P_mean(I,J-1)*Q_mean(I,J-1)/ &
                     Max(Depth(I,J-1)+ETAmean(I,J-1),MinDepthFrc))/2.0/DYg
        DyVVH(I,J) = (Q_mean(I,J+1)*Q_mean(I,J+1)/ &
                     Max(Depth(I,J+1)+ETAmean(I,J+1),MinDepthFrc) - &
                     Q_mean(I,J-1)*Q_mean(I,J-1)/ &
                     Max(Depth(I,J-1)+ETAmean(I,J-1),MinDepthFrc))/2.0/DYg
        DxUVH(I,J) = (P_mean(I+1,J)*Q_mean(I+1,J)/ &
                     Max(Depth(I+1,J)+ETAmean(I+1,J),MinDepthFrc) - &
                     P_mean(I-1,J)*Q_mean(I-1,J)/ &
                     Max(Depth(I-1,J)+ETAmean(I-1,J),MinDepthFrc))/2.0/DXg
        FRCXmean(I,J) = FRCXmean(I,J)*Cd(I,J)
        FRCYmean(I,J) = FRCYmean(I,J)*Cd(I,J)
      ENDDO
      ENDDO       
	  
	  !ykchoi includes ETAmean, fyshi added 03/22/2016, ykchoi move the 
          ! two statements right after IF because ETA2sum should use the previsou
          ! ETAmean. 04/20/2016
	  !ETA2sum = (Eta-ETAmean)*(Eta-ETAmean)*DT + ETA2sum 
	  !ETA2mean = ETA2sum/T_sum

        T_sum=T_sum-T_INTV_mean   ! T_sum=ZERO? (ykchoi)
        Usum=ZERO
        Vsum=ZERO
        ETAsum=ZERO
	ETA2sum=ZERO     !ykchoi


        U_davg_sum = ZERO
        V_davg_sum = ZERO
        P_sum = ZERO
        Q_sum = ZERO
        UUsum = ZERO
        UVsum = ZERO
        VVsum = ZERO
        WWsum = ZERO
        FRCXsum = ZERO
        FRCYsum = ZERO

	  SigWaveHeight = 4.004*SQRT( ETA2mean )  !ykcho

! wave height
       DO J=1,N
       DO I=1,M
        IF(Num_Zero_Up(I,J)>=2)THEN
          WaveHeightAve(I,J)=HavgSum(I,J)/(Num_Zero_Up(I,J)-1)
          WaveHeightRMS(I,J)=SQRT(HrmsSum(I,J)/(Num_Zero_Up(I,J)-1))
        ENDIF
!        Num_Zero_Up(I,J)=0
!        HavgSum(I,J)=ZERO
!        HrmsSum(I,J)=ZERO
       ENDDO
       ENDDO

        CALL PREVIEW_MEAN

      ELSE

        Usum=U*DT+Usum
        Vsum=V*DT+Vsum
        ETAsum=ETA*DT+ETAsum
	  !ykchoi, fyshi added ETAmean 03/22/2016
	ETA2sum = (Eta-ETAmean)*(Eta-ETAmean)*DT + ETA2sum

          ! fyshi added radiation calc 06/07/2022   
        U_davg_sum=U_davg*DT+U_davg_sum
        V_davg_sum=V_davg*DT+V_davg_sum
        P_sum=P_center*DT+P_sum
        Q_sum=Q_center*DT+Q_sum
        UUsum = (U_davg-U_davg_mean)*(U_davg-U_davg_mean)*(Eta-ETAmean+Depth)*DT + UUsum
        UVsum = (U_davg-U_davg_mean)*(V_davg-V_davg_mean)*(Eta-ETAmean+Depth)*DT + UVsum
        VVsum = (V_davg-V_davg_mean)*(V_davg-V_davg_mean)*(Eta-ETAmean+Depth)*DT + VVsum
        WWsum = 0.25*Wsurf*Wsurf*(Eta-ETAmean+Depth)*DT + WWsum
        FRCXsum = U*SQRT(U*U+V*V)*DT + FRCXsum
        FRCYsum = V*SQRT(U*U+V*V)*DT + FRCYsum

! wave height
       DO J=1,N
       DO I=1,M
         if(Eta(i,j)>Emax(i,j)) Emax(i,j) = Eta(i,j)
         if(Eta(i,j)<Emin(i,j)) Emin(i,j) = Eta(i,j)
         Tmpe = Eta(i,j)-ETAmean(i,j)
         Tmp_0 = Eta0(i,j)-ETAmean(i,j)
         if(Tmpe>Tmp_0.and.Tmpe*Tmp_0<=Zero) then
           Num_Zero_Up(i,j) = Num_Zero_Up(i,j)+1
           if(Num_Zero_Up(i,j)>=2) then
               HavgSum(i,j) = HavgSum(i,j)+Emax(i,j)-Emin(i,j)
               HrmsSum(i,j) = HrmsSum(i,j)+(Emax(i,j)-Emin(i,j))**2
           endif
           ! reset Emax and Emin to find next wave
           Emax(i,j) = -1000.
           Emin(i,j) = 1000.
         endif  
       ENDDO
       ENDDO

      ENDIF  ! end average time

      END SUBROUTINE CALCULATE_MEAN
