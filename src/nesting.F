!------------------------------------------------------------------------------------
!
!      FILE nesting.F
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
# if defined (COUPLING)
!-------------------------------------------------------------------------------------
!
!    COUPLING_INITIAL is the subroutine to initialize coupling/nesting conditions
!    
!    HISTORY: 05/01/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE COUPLING_INITIAL
    USE GLOBAL
    USE INPUT_READ
    IMPLICIT NONE
    INTEGER :: I
  
 INQUIRE(FILE=TRIM(COUPLING_FILE),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(COUPLING_FILE), ' specified in input.txt but CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(COUPLING_FILE), ' specified in input.txt but CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

       OPEN(11,FILE=TRIM(COUPLING_FILE))


         READ(11,*)  ! title
         READ(11,*)  ! boundary info
! boundary basic info including point number of coupling, start point, etc
! east
         READ(11,*)  ! east
         READ(11,*) N_COUPLING_EAST,J_START_EAST
! west 
         READ(11,*)  ! west
         READ(11,*) N_COUPLING_WEST,J_START_WEST
! south 
         READ(11,*)  ! south
         READ(11,*) N_COUPLING_SOUTH,I_START_SOUTH
! north 
         READ(11,*)  ! north
         READ(11,*) N_COUPLING_NORTH,I_START_NORTH

! read time and variable at the first level

         READ(11,*) ! time start title
         READ(11,*) TIME_COUPLING_1 
! initialize time_2
         TIME_COUPLING_2 = TIME_COUPLING_1

! east
         IF(N_COUPLING_EAST.GT.0)THEN
           ALLOCATE(U_COUPLING_EAST(N_COUPLING_EAST,2),&
               V_COUPLING_EAST(N_COUPLING_EAST,2),&
               Z_COUPLING_EAST(N_COUPLING_EAST,2))
             READ(11,*)   ! east
             READ(11,119)(U_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119)(V_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
!   initialize first step
             U_COUPLING_EAST(:,1)=U_COUPLING_EAST(:,2)
             V_COUPLING_EAST(:,1)=V_COUPLING_EAST(:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_east
119      FORMAT(5E16.6)

! west
         IF(N_COUPLING_WEST.GT.0)THEN
           ALLOCATE(U_COUPLING_WEST(N_COUPLING_WEST,2),&
               V_COUPLING_WEST(N_COUPLING_WEST,2),&
               Z_COUPLING_WEST(N_COUPLING_WEST,2))
             READ(11,*)   ! west
             READ(11,119)(U_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119)(V_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
!   initialize first step
             U_COUPLING_WEST(:,1)=U_COUPLING_WEST(:,2)
             V_COUPLING_WEST(:,1)=V_COUPLING_WEST(:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_west
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
           ALLOCATE(U_COUPLING_SOUTH(N_COUPLING_SOUTH,2),&
               V_COUPLING_SOUTH(N_COUPLING_SOUTH,2),&
               Z_COUPLING_SOUTH(N_COUPLING_SOUTH,2))
             READ(11,*)   ! south
             READ(11,119)(U_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119)(V_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
!   initialize first step
             U_COUPLING_SOUTH(:,1)=U_COUPLING_SOUTH(:,2)
             V_COUPLING_SOUTH(:,1)=V_COUPLING_SOUTH(:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_south
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
           ALLOCATE(U_COUPLING_NORTH(N_COUPLING_NORTH,2),&
               V_COUPLING_NORTH(N_COUPLING_NORTH,2),&
               Z_COUPLING_NORTH(N_COUPLING_NORTH,2))
             READ(11,*)   ! north
             READ(11,119)(U_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119)(V_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)

!   initialize first step
             U_COUPLING_NORTH(:,1)=U_COUPLING_NORTH(:,2)
             V_COUPLING_NORTH(:,1)=V_COUPLING_NORTH(:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_north


! specify boundary start points

! west boundary
   IF(N_COUPLING_WEST>0)THEN
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
![---ykchoi Jan/23/2018
!      Kstart_WEST=J_START_WEST+Nghost -npy*Nglob/py
!      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1 -npy*Nglob/py
      Kstart_WEST=J_START_WEST+Nghost - (jjsta - 1)
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1 - (jjsta - 1)
!---ykchoi Jan/23/2018]
      IF((Kstart_WEST>Nghost.AND.Kstart_WEST<Nloc-Nghost+1).OR.&
     (Kend_WEST<Nloc-Nghost+1.AND.Kend_WEST>Nghost))THEN
       IF(Kstart_WEST<Nghost+1)THEN
         Kshift_WEST= -Kstart_WEST+Nghost+1
         Kstart_WEST=Nghost+1
       ELSE
         Kshift_WEST=-(Kstart_WEST-Nghost)+1
       ENDIF
       IF(Kend_WEST>Nloc-Nghost)THEN
         Kend_WEST=Nloc-Nghost
       ENDIF
       IN_DOMAIN_WEST=.TRUE.
      ELSE
       IF(Kstart_WEST<=Nghost.AND.Kend_WEST>=Nloc-Nghost+1)THEN
         Kshift_WEST=-Kstart_WEST+Nghost+1
         Kstart_WEST = Nghost+1
         Kend_WEST = Nloc-Nghost
         IN_DOMAIN_WEST=.TRUE.
       ELSE
         IN_DOMAIN_WEST=.FALSE.
       ENDIF
      ENDIF
! check print*,myid,Kshift_WEST,Kstart_WEST,Kend_WEST
     endif
# else
      Kstart_WEST=J_START_WEST+Nghost
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1
      Kshift_WEST = -(Kstart_WEST-Nghost)+1
      IN_DOMAIN_WEST = .TRUE.
# endif

   ENDIF

! east boundary
   IF(N_COUPLING_EAST>0)THEN
# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
![---ykchoi Jan/23/2018
!      Kstart_EAST=J_START_EAST+Nghost -npy*Nglob/py
!      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1 -npy*Nglob/py
      Kstart_EAST=J_START_EAST+Nghost - (jjsta - 1)
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1 - (jjsta - 1)
!---ykchoi Jan/23/2018]
      IF((Kstart_EAST>Nghost.AND.Kstart_EAST<Nloc-Nghost+1).OR.&
     (Kend_EAST<Nloc-Nghost+1.AND.Kend_EAST>Nghost))THEN
       IF(Kstart_EAST<Nghost+1)THEN
         Kshift_EAST= -Kstart_EAST+Nghost+1
         Kstart_EAST=Nghost+1
       ELSE
         Kshift_EAST=-(Kstart_EAST-Nghost)+1
       ENDIF
       IF(Kend_EAST>Nloc-Nghost)THEN
         Kend_EAST=Nloc-Nghost
       ENDIF
       IN_DOMAIN_EAST=.TRUE.
      ELSE

       IF(Kstart_EAST<=Nghost.AND.Kend_EAST>=Nloc-Nghost+1)THEN
         Kshift_EAST=-Kstart_EAST+Nghost+1
         Kstart_EAST = Nghost+1
         Kend_EAST = Nloc-Nghost
         IN_DOMAIN_EAST=.TRUE.
       ELSE
         IN_DOMAIN_EAST=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_EAST=J_START_EAST+Nghost
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1
      Kshift_EAST = -(Kstart_EAST-Nghost)+1
      IN_DOMAIN_EAST = .TRUE.
# endif
    ENDIF

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
![---ykchoi Jan/23/2018
!      Kstart_SOUTH=I_START_SOUTH+Nghost -npx*Mglob/px
!      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1 -npx*Mglob/px
      Kstart_SOUTH=I_START_SOUTH+Nghost -(iista - 1)
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1 -(iista - 1)
!---ykchoi Jan/23/2018]
      IF((Kstart_SOUTH>Nghost.AND.Kstart_SOUTH<Mloc-Nghost+1).OR.&
     (Kend_SOUTH<Mloc-Nghost+1.AND.Kend_SOUTH>Nghost))THEN
       IF(Kstart_SOUTH<Nghost+1)THEN
         Kshift_SOUTH= -Kstart_SOUTH+Nghost+1
         Kstart_SOUTH=Nghost+1
       ELSE
         Kshift_SOUTH=-(Kstart_SOUTH-Nghost)+1
       ENDIF
       IF(Kend_SOUTH>Mloc-Nghost)THEN
         Kend_SOUTH=Mloc-Nghost
       ENDIF
       IN_DOMAIN_SOUTH=.TRUE.
      ELSE

       IF(Kstart_SOUTH<=Nghost.AND.Kend_SOUTH>=Mloc-Nghost+1)THEN
         Kshift_SOUTH=-Kstart_SOUTH+Nghost+1
         Kstart_SOUTH = Nghost+1
         Kend_SOUTH = Mloc-Nghost
         IN_DOMAIN_SOUTH=.TRUE.
       ELSE
         IN_DOMAIN_SOUTH=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_SOUTH=I_START_SOUTH+Nghost
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1
      Kshift_SOUTH = -(Kstart_SOUTH-Nghost)+1
      IN_DOMAIN_SOUTH = .TRUE.
# endif
   ENDIF

! north boundary
   IF(N_COUPLING_NORTH>0)THEN
# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
![---ykchoi Jan/23/2018
!      Kstart_NORTH=I_START_NORTH+Nghost -npx*Mglob/px
!      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1 -npx*Mglob/px
      Kstart_NORTH=I_START_NORTH+Nghost - (iista - 1)
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1 - (iista - 1)
!---ykchoi Jan/23/2018]
      IF((Kstart_NORTH>Nghost.AND.Kstart_NORTH<Mloc-Nghost+1).OR.&
     (Kend_NORTH<Mloc-Nghost+1.AND.Kend_NORTH>Nghost))THEN
       IF(Kstart_NORTH<Nghost+1)THEN
         Kshift_NORTH= -Kstart_NORTH+Nghost+1
         Kstart_NORTH=Nghost+1
       ELSE
         Kshift_NORTH=-(Kstart_NORTH-Nghost)+1
       ENDIF
       IF(Kend_NORTH>Mloc-Nghost)THEN
         Kend_NORTH=Mloc-Nghost
       ENDIF
       IN_DOMAIN_NORTH=.TRUE.
      ELSE

       IF(Kstart_NORTH<=Nghost.AND.Kend_NORTH>=Mloc-Nghost+1)THEN
         Kshift_NORTH=-Kstart_NORTH+Nghost+1
         Kstart_NORTH = Nghost+1
         Kend_NORTH = Mloc-Nghost
         IN_DOMAIN_NORTH=.TRUE.
       ELSE
         IN_DOMAIN_NORTH=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_NORTH=I_START_NORTH+Nghost
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1
      Kshift_NORTH = -(Kstart_NORTH-Nghost)+1
      IN_DOMAIN_NORTH = .TRUE.
# endif
   ENDIF

END SUBROUTINE COUPLING_INITIAL
# endif 


# if defined (COUPLING)
!-------------------------------------------------------------------------------------
!
!   OneWayCoupling subroutine is used to pass coupling variables into ghost cells                                                         
!
!   HISTORY: 04/28/2012 Fengyan Shi   
!   1) used format I/O  2) read through read(11,119)                                      
!
!-------------------------------------------------------------------------------------
SUBROUTINE OneWayCoupling
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: I,J
    REAL(SP) :: tmp1,tmp2

119      FORMAT(5E16.6)  ! this is a fixed format for I/O
 
! determine time slot

    IF(TIME>TIME_COUPLING_1.AND.TIME>TIME_COUPLING_2) THEN
         TIME_COUPLING_1=TIME_COUPLING_2

         READ(11,*,END=120) TIME_COUPLING_2 
! east
         IF(N_COUPLING_EAST.GT.0)THEN
             READ(11,*,END=120)   ! east

             U_COUPLING_EAST(:,1)=U_COUPLING_EAST(:,2)
             V_COUPLING_EAST(:,1)=V_COUPLING_EAST(:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)

             READ(11,119,END=120)(U_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119,END=120)(V_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119,END=120)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
         ELSE
             READ(11,*,END=120)   ! east            
         ENDIF
! west
         IF(N_COUPLING_WEST.GT.0)THEN
             READ(11,*,END=120)   ! west

             U_COUPLING_WEST(:,1)=U_COUPLING_WEST(:,2)
             V_COUPLING_WEST(:,1)=V_COUPLING_WEST(:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)


             READ(11,119,END=120)(U_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119,END=120)(V_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119,END=120)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
         ELSE
             READ(11,*,END=120)   ! west            
         ENDIF
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
             READ(11,*,END=120)   ! south

             U_COUPLING_SOUTH(:,1)=U_COUPLING_SOUTH(:,2)
             V_COUPLING_SOUTH(:,1)=V_COUPLING_SOUTH(:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)


             READ(11,119,END=120)(U_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119,END=120)(V_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119,END=120)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
         ELSE
             READ(11,*,END=120)   ! south            
         ENDIF
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
             READ(11,*,END=120)   ! north

             U_COUPLING_NORTH(:,1)=U_COUPLING_NORTH(:,2)
             V_COUPLING_NORTH(:,1)=V_COUPLING_NORTH(:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)

             READ(11,119,END=120)(U_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119,END=120)(V_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119,END=120)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
         ELSE
             READ(11,*,END=120)   ! north            
         ENDIF


    ENDIF  ! time>time_2 and time_1

120 CONTINUE

    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TIME_COUPLING_1)THEN
      IF(TIME_COUPLING_1.EQ.TIME_COUPLING_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TIME_COUPLING_2-TIME) &
            /MAX(SMALL, ABS(TIME_COUPLING_2-TIME_COUPLING_1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

! west boundary
   IF(N_COUPLING_WEST>0)THEN
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_WEST)THEN

      DO J=Kstart_WEST,Kend_WEST 
      DO I=1,Nghost
        ETA(I,J)=Z_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +Z_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2
        U(I,J)=U_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +U_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2
        V(I,J)=V_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +V_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2
        HU(I,J)=(Depth(I,J)+ETA(I,J))*U(I,J)
        HV(I,J)=(Depth(I,J)+ETA(I,J))*V(I,J)
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_west>0

! east boundary
   IF(N_COUPLING_EAST>0)THEN
# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_EAST)THEN
      DO J=Kstart_EAST,Kend_EAST  
      DO I=Iend+1,Iend+Nghost
        ETA(I,J)=Z_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +Z_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2
        U(I,J)=U_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +U_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2
        V(I,J)=V_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +V_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2
        HU(I,J)=(Depth(I,J)+ETA(I,J))*U(I,J)
        HV(I,J)=(Depth(I,J)+ETA(I,J))*V(I,J)
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_east>0

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_SOUTH)THEN
      DO I=Kstart_SOUTH,Kend_SOUTH  
      DO J=1,Nghost
        ETA(I,J)=Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2
        U(I,J)=U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2
        V(I,J)=V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2
        HU(I,J)=(Depth(I,J)+ETA(I,J))*U(I,J)
        HV(I,J)=(Depth(I,J)+ETA(I,J))*V(I,J)
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_south>0

! north boundary
   IF(N_COUPLING_NORTH>0)THEN
# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_NORTH)THEN
      DO I=Kstart_NORTH,Kend_NORTH  
      DO J=Jend+1,Jend+Nghost
        ETA(I,J)=Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2
        U(I,J)=U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2
        V(I,J)=V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2
        HU(I,J)=(Depth(I,J)+ETA(I,J))*U(I,J)
        HV(I,J)=(Depth(I,J)+ETA(I,J))*V(I,J)
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_north>0

END SUBROUTINE OneWayCoupling
# endif 
  
