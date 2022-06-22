!------------------------------------------------------------------------------------
!
!      FILE mod_tide.F
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
!  TIDE_MODULE is a module to add tide/surge boundary conditions into wave simulation    
!
!  HISTORY :
!    11/26/2019 Fengyan Shi
!-------------------------------------------------------------------------------------

MODULE TIDE_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Mloc1,Nloc1,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,U,V,P,Q,MinDepth,MASK,DT,Gamma3,Depth,tmp4preview, &
                    ALPHA,BETA,MASK9,DepthX,DepthY,PERIODIC, &
                    UNDERTOW_U, UNDERTOW_V, ROLLER_SWITCH,ROLLER,Mglob,Nglob, &
                    TIME,Sponge_west_width,Sponge_east_width,Sponge_south_width, &
                    Sponge_north_width
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

  CHARACTER(LEN=80) TideBcType
  LOGICAL :: TIDAL_BC=.FALSE.
  LOGICAL :: TideEast=.TRUE.
  LOGICAL :: TideWest=.TRUE.
  LOGICAL :: TideNorth=.TRUE.
  LOGICAL :: TideSouth=.TRUE.
  REAL(SP) :: TideEast_U,TideEast_V,TideEast_ETA
  REAL(SP) :: TideWest_U,TideWest_V,TideWest_ETA
  REAL(SP) :: TideNorth_U,TideNorth_V,TideNorth_ETA
  REAL(SP) :: TideSouth_U,TideSouth_V,TideSouth_ETA
  REAL(SP),DIMENSION(:,:),ALLOCATABLE :: SPONGE_TIDE_WEST, &
              SPONGE_TIDE_EAST,SPONGE_TIDE_SOUTH,SPONGE_TIDE_NORTH
  INTEGER :: Iwidth
  

CONTAINS
  
SUBROUTINE TIDE_INITIAL
  USE GLOBAL,ONLY : INPUT_FILE_NAME
  USE INPUT_READ
  IMPLICIT NONE

  CHARACTER(LEN=80)::FILE_NAME=' '
  INTEGER :: ierr


! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME

      CALL READ_LOGICAL(TIDAL_BC,FILE_NAME,'TIDAL_BC',ierr)
      IF(ierr == 1)THEN
       TIDAL_BC = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A80)')'TIDAL_BC not defined, Default: False'
       WRITE(*,'(A80)')'TIDAL_BC not defined, Default: False'
      endif
# else
      WRITE(3,'(A80)')'TIDAL_BC not defined, Default: False'
      WRITE(*,'(A80)')'TIDAL_BC not defined, Default: False'
# endif 
      ENDIF

    IF (TIDAL_BC) THEN

# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A80)')'TIDAL_BC is defined'
       WRITE(*,'(A80)')'TIDAL_BC is defined'
      endif
# else
      WRITE(3,'(A80)')'TIDAL_BC is defined'
      WRITE(*,'(A80)')'TIDAL_BC is defined'
# endif 

      ALLOCATE(SPONGE_TIDE_WEST(Mloc,Nloc),SPONGE_TIDE_EAST(Mloc,Nloc), &
                SPONGE_TIDE_SOUTH(Mloc,Nloc),SPONGE_TIDE_NORTH(Mloc,Nloc))

      CALL TIDE_SPONGE

! input tide details
      CALL READ_STRING(TideBcType,FILE_NAME,'TideBcType',ierr)
      IF(ierr==1)THEN
        TideBcType = 'CONSTANT'
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You dont specify TideBcType, use CONSTANT.'
         WRITE(3,'(A50)')'You dont specify TideBcType, use CONSTANT.'
      endif
# else
         WRITE(*,'(A50)')'You dont specify TideBcType, use CONSTANT.'
         WRITE(3,'(A50)')'You dont specify TideBcType, use CONSTANT.'
# endif
      ENDIF

    IF(TideBcType(1:4)=='CONS')THEN

!  west ---

      CALL READ_Float(TideWest_ETA,FILE_NAME,'TideWest_ETA',ierr)

      IF(ierr==1)THEN
        TideWest = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideWest_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideWest_ETA, use default: FALSE'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideWest_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideWest_ETA, use default: FALSE'
# endif
      ENDIF

      CALL READ_Float(TideWest_U,FILE_NAME,'TideWest_U',ierr)

      IF(ierr==1)THEN
        TideWest_U = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideWest_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideWest_U, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideWest_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideWest_U, use default: 0'
# endif
      ENDIF

      CALL READ_Float(TideWest_V,FILE_NAME,'TideWest_V',ierr)

      IF(ierr==1)THEN
        TideWest_V = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideWest_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideWest_V, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideWest_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideWest_V, use default: 0'
# endif
      ENDIF


!  east ---

      CALL READ_Float(TideEast_ETA,FILE_NAME,'TideEast_ETA',ierr)

      IF(ierr==1)THEN
        TideEast = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideEast_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideEast_ETA, use default: FALSE'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideEast_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideEast_ETA, use default: FALSE'
# endif
      ENDIF

      CALL READ_Float(TideEast_U,FILE_NAME,'TideEast_U',ierr)

      IF(ierr==1)THEN
        TideEast_U = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideEast_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideEast_U, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideEast_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideEast_U, use default: 0'
# endif
      ENDIF

      CALL READ_Float(TideEast_V,FILE_NAME,'TideEast_V',ierr)

      IF(ierr==1)THEN
        TideEast_V = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideEast_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideEast_V, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideEast_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideEast_V, use default: 0'
# endif
      ENDIF


!  south ---

      CALL READ_Float(TideSouth_ETA,FILE_NAME,'TideSouth_ETA',ierr)

      IF(ierr==1)THEN
        TideSouth = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideSouth_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideSouth_ETA, use default: FALSE'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideSouth_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideSouth_ETA, use default: FALSE'
# endif
      ENDIF

      CALL READ_Float(TideSouth_U,FILE_NAME,'TideSouth_U',ierr)

      IF(ierr==1)THEN
        TideSouth_U = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideSouth_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideSouth_U, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideSouth_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideSouth_U, use default: 0'
# endif
      ENDIF

      CALL READ_Float(TideSouth_V,FILE_NAME,'TideSouth_V',ierr)

      IF(ierr==1)THEN
        TideSouth_V = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideSouth_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideSouth_V, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideSouth_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideSouth_V, use default: 0'
# endif
      ENDIF

!  north ---

      CALL READ_Float(TideNorth_ETA,FILE_NAME,'TideNorth_ETA',ierr)

      IF(ierr==1)THEN
        TideNorth = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideNorth_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideNorth_ETA, use default: FALSE'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideNorth_ETA, use default: FALSE'
         WRITE(3,'(A80)')'You dont specify TideNorth_ETA, use default: FALSE'
# endif
      ENDIF

      CALL READ_Float(TideNorth_U,FILE_NAME,'TideNorth_U',ierr)

      IF(ierr==1)THEN
        TideNorth_U = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideNorth_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideNorth_U, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideNorth_U, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideNorth_U, use default: 0'
# endif
      ENDIF

      CALL READ_Float(TideNorth_V,FILE_NAME,'TideNorth_V',ierr)

      IF(ierr==1)THEN
        TideNorth_V = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify TideNorth_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideNorth_V, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify TideNorth_V, use default: 0'
         WRITE(3,'(A80)')'You dont specify TideNorth_V, use default: 0'
# endif
      ENDIF

      IF(TideWest)THEN

        Sponge_west_width = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You sed TideWest BC,  Sponge_west is invalid '
         WRITE(3,'(A80)')'You sed TideWest BC,  Sponge_west is invalid '
      endif
# else
         WRITE(*,'(A80)')'You sed TideWest BC,  Sponge_west is invalid '
         WRITE(3,'(A80)')'You sed TideWest BC,  Sponge_west is invalid '
# endif         
      ENDIF

      IF(TideEast)THEN

        Sponge_east_width = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You sed TideEast BC,  Sponge_east is invalid '
         WRITE(3,'(A80)')'You sed TideEast BC,  Sponge_east is invalid '
      endif
# else
         WRITE(*,'(A80)')'You sed TideEast BC,  Sponge_east is invalid '
         WRITE(3,'(A80)')'You sed TideEast BC,  Sponge_east is invalid '
# endif         
      ENDIF

      IF(TideSouth)THEN

        Sponge_south_width = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You sed TideSouth BC,  Sponge_south is invalid '
         WRITE(3,'(A80)')'You sed TideSouth BC,  Sponge_south is invalid '
      endif
# else
         WRITE(*,'(A80)')'You sed TideSouth BC,  Sponge_south is invalid '
         WRITE(3,'(A80)')'You sed TideSouth BC,  Sponge_south is invalid '
# endif         
      ENDIF

      IF(TideNorth)THEN

        Sponge_North_width = ZERO
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You sed TideNorth BC,  Sponge_North is invalid '
         WRITE(3,'(A80)')'You sed TideNorth BC,  Sponge_North is invalid '
      endif
# else
         WRITE(*,'(A80)')'You sed TideNorth BC,  Sponge_North is invalid '
         WRITE(3,'(A80)')'You sed TideNorth BC,  Sponge_North is invalid '
# endif         
      ENDIF


      IF((.NOT.TideWest).AND.(.NOT.TideEast)  &
                        .AND.(.NOT.TideEast).AND.(.NOT.TideEast)) THEN
        TIDAL_BC = .FALSE.

# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You specify Tide but no BCs, set tide as false'
         WRITE(3,'(A80)')'You dont specify TideNorth_V, use default: 0'
      endif
# else
         WRITE(*,'(A80)')'You specify Tide but no BCs, set tide as false'
         WRITE(3,'(A80)')'You specify Tide but no BCs, set tide as false'
# endif      

      ENDIF ! end all false
 
     ENDIF ! end tide constant

    ENDIF ! end tide


END SUBROUTINE TIDE_INITIAL

SUBROUTINE TIDE_SPONGE
  IMPLICIT NONE
  REAL(SP) :: ri,lim,R_sp,A_sp
  INTEGER :: i,j

  Iwidth = 30
  lim = 1.0_SP
  R_sp = 0.85_SP
  A_sp = 10.0_SP

! west

       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(i+npx*Mglob/px-1)/(Iwidth-1))
# else
         ri = R_sp**(50*(i-1)/(Iwidth-1))
# endif
         SPONGE_TIDE_WEST(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! east
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(Mloc-i+(px-npx-1)*Mglob/px)/(Iwidth-1))
# else
         ri = R_sp**(50*(Mloc-i)/(Iwidth-1))
# endif
         SPONGE_TIDE_EAST(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! south

       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(j+npy*Nglob/py-1)/(Iwidth-1))
# else
         ri = R_sp**(50*(j-1)/(Iwidth-1))
# endif
         SPONGE_TIDE_SOUTH(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

! north
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_sp**(50*(Nloc-j+(py-npy-1)*Nglob/py)/(Iwidth-1))
# else
         ri = R_sp**(50*(Nloc-i)/(Iwidth-1))
# endif
         SPONGE_TIDE_NORTH(i,j) = max(A_sp**ri,lim)
       enddo
       enddo

       SPONGE_TIDE_WEST = 1.0_SP/SPONGE_TIDE_WEST
       SPONGE_TIDE_EAST = 1.0_SP/SPONGE_TIDE_EAST
       SPONGE_TIDE_SOUTH = 1.0_SP/SPONGE_TIDE_SOUTH
       SPONGE_TIDE_NORTH = 1.0_SP/SPONGE_TIDE_NORTH



END SUBROUTINE TIDE_SPONGE

SUBROUTINE TIDE_BC
  IMPLICIT NONE
  INTEGER :: I,J

   IF(TideWest)THEN
     DO J=1,Nloc
     DO I=1,Iwidth
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideWest_ETA +(ETA(I,J)-TideWest_ETA)*SPONGE_TIDE_WEST(I,J)
         U(I,J)=TideWest_U +(U(I,J)-TideWest_U)*SPONGE_TIDE_WEST(I,J)
         V(I,J)=TideWest_V +(V(I,J)-TideWest_V)*SPONGE_TIDE_WEST(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(TideEast)THEN
     DO J=1,Nloc
     DO I=Mloc-Iwidth+1,Mloc
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideEast_ETA +(ETA(I,J)-TideEast_ETA)*SPONGE_TIDE_EAST(I,J)
         U(I,J)=TideEast_U +(U(I,J)-TideEast_U)*SPONGE_TIDE_EAST(I,J)
         V(I,J)=TideEast_V +(V(I,J)-TideEast_V)*SPONGE_TIDE_EAST(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(TideSouth)THEN
     DO J=1,Iwidth
     DO I=1,Mloc
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideSouth_ETA +(ETA(I,J)-TideSouth_ETA)*SPONGE_TIDE_SOUTH(I,J)
         U(I,J)=TideSouth_U +(U(I,J)-TideSouth_U)*SPONGE_TIDE_SOUTH(I,J)
         V(I,J)=TideSouth_V +(V(I,J)-TideSouth_V)*SPONGE_TIDE_SOUTH(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(TideNorth)THEN
     DO J=Nloc-Iwidth+1,Nloc
     DO I=1,Mloc
       IF(MASK(I,J)==1)THEN
         ETA(I,J)=TideNorth_ETA +(ETA(I,J)-TideNorth_ETA)*SPONGE_TIDE_NORTH(I,J)
         U(I,J)=TideNorth_U +(U(I,J)-TideNorth_U)*SPONGE_TIDE_NORTH(I,J)
         V(I,J)=TideNorth_V +(V(I,J)-TideNorth_V)*SPONGE_TIDE_NORTH(I,J)
       ENDIF
     ENDDO
     ENDDO
   ENDIF

END SUBROUTINE TIDE_BC

END MODULE TIDE_MODULE



    



