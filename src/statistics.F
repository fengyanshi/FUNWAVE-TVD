!------------------------------------------------------------------------------------
!
!      FILE statistics.F
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
!    STATISTICS is subroutine to calculate statistics
!
!    HISTORY: 
!       05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE STATISTICS
     USE GLOBAL
     IMPLICIT NONE

     REAL(SP)::MassVolume=ZERO,Energy=ZERO,MaxEta=ZERO,MinEta=ZERO, &
              MaxU=ZERO,MaxV=ZERO,Fr=ZERO,UTotal=ZERO,UTotalMax=ZERO
# if defined (PARALLEL)
     REAL(SP)::myvar
# endif
!


     MassVolume=ZERO
     Energy=ZERO
     UTotalMax=ZERO

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend

! Vol=SUM(Eta*dx*dy), reference is at z=0
! Energy=SUM(1/2*g*H^2*dx*dy+0.5*u^2*H*dx*dy)

# if defined (CARTESIAN)
       IF(MASK(I,J)>0)THEN
         MassVolume=MassVolume+Eta(I,J)*DX*DY
       ENDIF

       Energy=Energy+0.5_SP*H(I,J)*H(I,J)*GRAV*DX*DY &
             +0.5_SP*U(I,J)*U(I,J)*H(I,J)*DX*DY &
             +0.5_SP*V(I,J)*V(I,J)*H(I,J)*DX*DY
# else
       IF(MASK(I,J)>0)THEN
         MassVolume=MassVolume+Eta(I,J)*DX(I,J)*DY(I,J)
       ENDIF
       Energy=Energy+0.5_SP*H(I,J)*H(I,J)*GRAV*DX(I,J)*DY(I,J) &
             +0.5_SP*U(I,J)*U(I,J)*H(I,J)*DX(I,J)*DY(I,J) &
             +0.5_SP*V(I,J)*V(I,J)*H(I,J)*DX(I,J)*DY(I,J)
# endif      
     ENDDO
     ENDDO

     MaxEta=MAXVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MinEta=MINVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MaxU=MAXVAL(ABS(U(Ibeg:Iend,Jbeg:Jend)))
     MaxV=MAXVAL(ABS(V(Ibeg:Iend,Jbeg:Jend)))

! found Froude vs. max speed
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      IF(MASK(I,J)>ZERO)THEN
       Utotal=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
       IF(Utotal.gt.UtotalMax)THEN
         UtotalMax=Utotal
         Fr=SQRT(GRAV*Max(H(I,J),MinDepthfrc))
       ENDIF
      ENDIF
     ENDDO
     ENDDO
     IF(Fr==ZERO)Fr=SQRT(GRAV*MinDepthfrc)

# if defined (PARALLEL)
     call MPI_ALLREDUCE(MassVolume,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     MassVolume = myvar
     call MPI_ALLREDUCE(Energy,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     Energy = myvar
     call MPI_ALLREDUCE(MaxEta,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxEta = myvar
     call MPI_ALLREDUCE(MinEta,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     MinEta = myvar
     call MPI_ALLREDUCE(MaxU,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxU = myvar
     call MPI_ALLREDUCE(MaxV,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxV = myvar
     call MPI_ALLREDUCE(UTotalMax,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     UTotalMax = myvar
     call MPI_ALLREDUCE(Fr,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     Fr = myvar
# endif


# if defined (PARALLEL)
     if (myid.eq.0) then
# endif

! print screen
     WRITE(*,*) '----------------- STATISTICS ----------------'
     WRITE(*,*) ' TIME        DT'
     WRITE(*,101) Time, DT
     WRITE(*,*) ' MassVolume  Energy      MaxEta      MinEta      Max U       Max V '
     WRITE(*,101)  MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV
     WRITE(*,*) ' MaxTotalU   PhaseS      Froude      WetDryMass'
     WRITE(*,101) UTotalMax, Fr, UTotalMax/Fr, Dmass
! print log file
     WRITE(3,*) '----------------- STATISTICS ----------------'
     WRITE(3,*) ' TIME        DT'
     WRITE(3,101) Time, DT
     WRITE(3,*) ' MassVolume  Energy      MaxEta      MinEta      Max U       Max V '
     WRITE(3,101)  MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV
     WRITE(3,*) ' MaxTotalU   PhaseS      Froude      WetDryMass'
     WRITE(3,101)  UTotalMax, Fr, UTotalMax/Fr, Dmass
# if defined (PARALLEL)
     endif
# endif

101  FORMAT(6E12.4)

# if defined (PARALLEL)
     if (myid.eq.0) then
# endif

     CALL CHECK_STATISTICS(MassVolume)
     CALL CHECK_STATISTICS(Energy)
     CALL CHECK_STATISTICS(MaxEta)
     CALL CHECK_STATISTICS(MinEta)
     CALL CHECK_STATISTICS(MaxU)
     CALL CHECK_STATISTICS(MaxV)
     CALL CHECK_STATISTICS(UTotalMax)
     CALL CHECK_STATISTICS(Fr)
     CALL CHECK_STATISTICS(UTotalMax/Fr)
     CALL CHECK_STATISTICS(Dmass)
# if defined (PARALLEL)
     endif
# endif


END SUBROUTINE STATISTICS

SUBROUTINE CHECK_STATISTICS( val )

    USE GLOBAL
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: val


     IF (ISNAN(val)) THEN
       WRITE(*,*) 'NaN detected in statistics, stopping simulation!'
       WRITE(3,*) 'NaN detected in statistics, stopping simulation!'

#if defined (PARALLEL)
       CALL MPI_ABORT( MPI_COMM_WORLD , 9 , ier)
#endif
       STOP

     END IF

 END SUBROUTINE CHECK_STATISTICS





