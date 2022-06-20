!------------------------------------------------------------------------------------
!
!      FILE mod_bathy_correction.F
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
!  bathy_correction_module is a module to correct bathy/topo data if slope exceeds
!  some threahold such as slope > 1.0 (default)     
!
!  HISTORY :
!    04/24/2017  Fengyan Shi
!
!-------------------------------------------------------------------------------------

MODULE BATHY_CORRECTION_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY,Depth, INPUT_FILE_NAME
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Depth0, Depth1,Gradx0,Grady0, &
                            Gradx,Grady,Depth_org
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: MASK0
    REAL(SP) :: SmoothBelowDepth, SlopeCap

#if defined (PARALLEL)
    REAL(SP) :: myvar
    INTEGER :: myint
# endif

CONTAINS

SUBROUTINE CORRECTION
  USE GLOBAL,ONLY : LARGE,tmp1,MinDepthFrc

  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  INTEGER :: ierr,pointcount
  REAL(SP) :: change_percentage, threshold

     pointcount=0
     threshold = 0.05_SP

! read vessel number and folder from input.txt
      FILE_NAME=INPUT_FILE_NAME


      CALL READ_FLOAT(SmoothBelowDepth,FILE_NAME,'SmoothBelowDepth',ierr)
      IF(ierr == 1)THEN
       SmoothBelowDepth = -LARGE 
      ENDIF

      CALL READ_FLOAT(SlopeCap,FILE_NAME,'SlopeCap',ierr)
      IF(ierr == 1)THEN
       SlopeCap = 1.0_SP 
      ENDIF

# if defined (PARALLEL)
        IF (myid.eq.0) THEN
          write(3,*) 'Bathymetry correction ...'
          write(3,'(A25,F12.2)') 'SmoothBelowDepth:', SmoothBelowDepth
          write(3,'(A25,F12.2)') 'SlopeCap:', SlopeCap
          write(*,*) 'Bathymetry correction ...'
          write(*,'(A25,F12.2)') 'SmoothBelowDepth:', SmoothBelowDepth
          write(*,'(A25,F12.2)') 'SlopeCap:', SlopeCap
        ENDIF
# else
          write(3,*) 'Bathymetry correction ...'
          write(3,'(A25,F12.2)') 'SmoothBelowDepth:', SmoothBelowDepth
          write(3,'(A25,F12.2)') 'SlopeCap:', SlopeCap
          write(*,*) 'Bathymetry correction ...'
          write(*,'(A25,F12.2)') 'SmoothBelowDepth:', SmoothBelowDepth
          write(*,'(A25,F12.2)') 'SlopeCap:', SlopeCap
# endif


        ALLOCATE (Depth0(Mloc,Nloc), Depth1(Mloc,Nloc), Depth_org(Mloc,Nloc), &
                  Gradx0(Mloc,Nloc),Grady0(Mloc,Nloc), &
                  Gradx(Mloc,Nloc),Grady(Mloc,Nloc), &
                  MASK0(Mloc,Nloc))

!     mask smooth area
      MASK0=1
      Depth_org=Depth
      Depth0=Depth
      Depth1=Depth

      DO J=Jbeg-1,Jend+1
      DO I=Ibeg-1,Iend+1
        IF(Depth0(I,J)<SmoothBelowDepth) THEN
          MASK0(I,J) = 0
          MASK0(I+1,J) = 0
          MASK0(I-1,J) = 0
          MASK0(I,J+1) = 0
          MASK0(I,J-1) = 0
        ENDIF
      ENDDO
      ENDDO

# if defined (PARALLEL)
     call phi_int_exch(MASK0)
# endif

      Gradx0=0.0
      Grady0=0.0

!     initial slope
      DO J=Jbeg,Jend
      DO I=Ibeg,Iend
        IF(MASK0(I,J)==1)THEN
# if defined (CARTESIAN)
          Gradx0(I,J)=ABS(Depth0(I+1,J)-Depth0(I-1,J))/DX
          Grady0(I,J)=ABS(Depth0(I,J+1)-Depth0(I,J-1))/DY
# else
          Gradx0(I,J)=ABS(Depth0(I+1,J)-Depth0(I-1,J))/DX(I,J)
          Grady0(I,J)=ABS(Depth0(I,J+1)-Depth0(I,J-1))/DY(I,J)
# endif
        ENDIF
      ENDDO
      ENDDO


      change_percentage = 1.0_SP
      DO WHILE (change_percentage.GT.threshold .AND. pointcount .LE.1000)

# if defined (PARALLEL)
       IF(myid == 0)THEN
        write(3,'(A12,I4,A25,F6.2)') 'iteration: ', pointcount, & 
                   'convegence percentage: ', change_percentage
        write(*,'(A12,I4,A25,F6.2)') 'iteration: ', pointcount, & 
                   'convegence percentage: ', change_percentage
       ENDIF
# else
        write(3,'(A12,I4,A25,F6.2)') 'iteration: ', pointcount, & 
                   'convegence percentage: ', change_percentage
        write(*,'(A12,I4,A25,F6.2)') 'iteration: ', pointcount, & 
                   'convegence percentage: ', change_percentage
# endif

        change_percentage = 0.0_SP
        DO J=Jbeg,Jend
        DO I=Ibeg,Iend
        IF(MASK0(I,J)==1)THEN
# if defined (CARTESIAN)
          Gradx(I,J)=ABS(Depth0(I+1,J)-Depth0(I-1,J))/DX
          Grady(I,J)=ABS(Depth0(I,J+1)-Depth0(I,J-1))/DY
# else
          Gradx(I,J)=ABS(Depth0(I+1,J)-Depth0(I-1,J))/DX(I,J)
          Grady(I,J)=ABS(Depth0(I,J+1)-Depth0(I,J-1))/DY(I,J)
# endif

          IF(MAX(Gradx(I,J),Grady(I,J)).GT.SlopeCap) THEN
            Depth1(I,J)=0.4_SP*Depth0(I,J)   &
                +0.15_SP*(Depth0(I+1,J)+Depth0(I-1,J)+Depth0(I,J+1)+Depth0(I,J-1))
            tmp1=ABS(Depth0(I,J)-Depth1(I,J))/MAX(MinDepthFrc*10_SP,Depth0(I,J))
            IF(tmp1 .GT. change_percentage) change_percentage = tmp1

          ENDIF

        ENDIF ! end mask0
        ENDDO
        ENDDO    

# if defined (PARALLEL)
     call phi_exch(Depth1)
# endif

        Depth0 = Depth1

        PointCount = PointCount +1

# if defined (PARALLEL)
        call MPI_ALLREDUCE(change_percentage,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
        change_percentage = myvar
# endif

      ENDDO  ! do while


!        call MPI_WAITALL( nreq, req, status, ier )

# if defined (PARALLEL)

        IF (myid.eq.0) THEN

        write(3,'(A25,I10)') 'total iteration: ', pointcount
        write(*,'(A25,I10)') 'total iteration: ', pointcount
        ENDIF
# else 
        write(3,'(A25,I10)') 'total iteration: ', pointcount
        write(*,'(A25,I10)') 'total iteration: ', pointcount
# endif


      Depth = Depth0

# if defined (PARALLEL)
      call phi_exch (Depth)
      call phi_exch (Gradx)
      call phi_exch (Grady)
      call phi_exch (Gradx0)
      call phi_exch (Grady0)
# endif

      CALL OUTPUT_CORRECTION

END SUBROUTINE CORRECTION

SUBROUTINE OUTPUT_CORRECTION
     USE GLOBAL, ONLY : RESULT_FOLDER

     IMPLICIT NONE

     CHARACTER(LEN=80)::TMP_NAME=' '
     CHARACTER(LEN=80)::FDIR=' '

     FDIR=TRIM(RESULT_FOLDER)

        TMP_NAME = TRIM(FDIR)//'depth_org.txt'
        call PutFile(TMP_NAME,Depth_org)
        TMP_NAME = TRIM(FDIR)//'depth_change.txt'
        call PutFile(TMP_NAME,Depth_org-Depth0)
        TMP_NAME = TRIM(FDIR)//'gradx0.txt'
        call PutFile(TMP_NAME,Gradx0)
        TMP_NAME = TRIM(FDIR)//'grady0.txt'
        call PutFile(TMP_NAME,Grady0)
        TMP_NAME = TRIM(FDIR)//'gradx.txt'
        call PutFile(TMP_NAME,Gradx)
        TMP_NAME = TRIM(FDIR)//'grady.txt'
        call PutFile(TMP_NAME,Grady)

END SUBROUTINE OUTPUT_CORRECTION

END MODULE BATHY_CORRECTION_MODULE






  
