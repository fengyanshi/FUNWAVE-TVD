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
#define nop(x) associate( x => x ); end associate

SUBROUTINE MIXING_STUFF
     USE GLOBAL
     IMPLICIT NONE

! calculate mean for smagorinsky s mixing and wave height
      !ykchoi (be careful)
      !I think Umean, Vmean is not using in other routine.
      IF( time >= STEADY_TIME )THEN    
      CALL CALCULATE_MEAN(T_INTV_mean,T_sum,DT,Mloc,Nloc,U,V,ETA,ETA0,&
           Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,WaveHeightRMS, &
           WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
           HrmsSum,HavgSum, &
           !ykchoi
           ETA2sum, ETA2mean, SigWaveHeight)
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
      SUBROUTINE CALCULATE_MEAN(T_INTV_mean,T_sum,DT,M,N,U,V,ETA,ETA0,&
                Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,&
                WaveHeightRMS, &
                WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
                HrmsSum,HavgSum, &
                !ykchoi
                    ETA2sum,ETA2mean,SigWaveHeight)
! calculate mean for smagorinsky s mixing and wave height

      USE PARAM
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
      REAL(SP),DIMENSION(M,N),INTENT(IN)::U,V,ETA,ETA0
      REAL(SP),INTENT(IN) :: T_INTV_mean,DT
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Umean,Vmean
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: WaveHeightRMS,WaveHeightAve
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) ::ETAmean
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Usum,Vsum,ETAsum
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: HrmsSum,HavgSum
      REAL(SP),INTENT(OUT) :: T_sum
      REAL(SP)::Tmpe,Tmp_0
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Emax,Emin
      INTEGER,DIMENSION(M,N),INTENT(INOUT) :: Num_Zero_Up

      !ykchoi
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: ETA2sum,ETA2mean,SigWaveHeight
      INTEGER :: I,J
      nop(Ibeg)
      nop(Iend)
      nop(Jbeg)
      nop(Jend)
      
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

        SigWaveHeight = 4.004*SQRT( ETA2mean )  !ykchoi

! wave height
       DO J=1,N
       DO I=1,M
        IF(Num_Zero_Up(I,J)>=2)THEN
          WaveHeightAve(I,J)=HavgSum(I,J)/Num_Zero_Up(I,J)
          WaveHeightRMS(I,J)=SQRT(HrmsSum(I,J)/Num_Zero_Up(I,J))
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
