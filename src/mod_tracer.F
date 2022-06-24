!------------------------------------------------------------------------------------
!
!      FILE mod_tracer.F
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
!  VESSEL is a module to model ship-wakes    
!
!  HISTORY :
!    09/19/2018  Fengyan Shi
!
!-------------------------------------------------------------------------------------
#if defined (TRACKING)

MODULE TRACER
     USE PARAM
     USE GLOBAL, ONLY : TIME, Jbeg,Jend,Ibeg,Iend,Mloc,Nloc,&
                 Mloc1,Nloc1,MASK,U,V,DT,&
                 DX,DY, Xco,Yco, RESULT_FOLDER
     USE INPUT_READ

# if defined (PARALLEL)
     USE MPI
     USE GLOBAL, ONLY : myid,ier,npx,npy,NumberProcessor,PX,PY,MPI_SP,&
                        n_east,n_west,n_suth,n_nrth
# endif
     IMPLICIT NONE


     REAL(SP),SAVE :: TRACK_TIME_COUNT = 0.0_SP
     REAL(SP),SAVE,DIMENSION(:),ALLOCATABLE :: X_TRACK,Y_TRACK,Sc,S1,S2,S3,&
                   U_TRACER,V_TRACER,X_TRACK_PRE,Y_TRACK_PRE, &
                   Sc_pre,S1_pre,S2_pre,S3_pre
                       
     INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: nx1,ny1,nx2,ny2,nx3,ny3,&
        nx1_pre,ny1_pre,nx2_pre,ny2_pre,nx3_pre,ny3_pre
     INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: LAYER_012
     LOGICAL,SAVE :: FirstCallTrack = .TRUE.
     REAL(SP),SAVE,DIMENSION(:),ALLOCATABLE :: TRACKER_START
     REAL(SP),SAVE,DIMENSION(:,:),ALLOCATABLE :: Usurf,Vsurf,Ubott,Vbott
     LOGICAL,SAVE,DIMENSION(:),ALLOCATABLE :: IN_CELL,STOP_SEARCH
     INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: FOUND_IN_DOMAIN
     REAL(SP), SAVE :: PLOT_INTV_TRACKING,PLOT_COUNT_TRACKING

# if defined (PARALLEL)
     INTEGER::myint
# endif
     
     INTEGER,SAVE :: NumTracker=1

     CONTAINS

SUBROUTINE TRACER_INITIAL

      USE GLOBAL, ONLY : INPUT_FILE_NAME
      IMPLICIT NONE
      INTEGER :: I_tmp, ierr
      CHARACTER(LEN=80) :: WHAT
      CHARACTER(LEN=80) :: TRACER_FILE,FILE_NAME

! read  from input.txt
      FILE_NAME=INPUT_FILE_NAME

      CALL READ_STRING(TRACER_FILE,FILE_NAME,'TRACER_FILE',ierr) 

       IF(ierr==1)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'You use tracking ', 'but no TRACER_FILE found, STOP'
         WRITE(3,'(A40,A40)')'You use tracking ', 'but no TRACER_FILE found, STOP'
      endif
       call MPI_FINALIZE ( ier )
# else
         WRITE(*,'(A40,A40)')'You use tracking ', 'but no TRACER_FILE found, STOP'
         WRITE(3,'(A40,A40)')'You use tracking ', 'but no TRACER_FILE found, STOP'
# endif
        STOP
      ENDIF


      OPEN(1,FILE=TRIM(TRACER_FILE))
        READ(1,*)WHAT
        READ(1,*)NumTracker
        READ(1,*)PLOT_INTV_TRACKING

       ALLOCATE(X_TRACK(NumTracker), &
                Y_TRACK(NumTracker),TRACKER_START(NumTracker), &
                LAYER_012(NumTracker) )
       ALLOCATE(Usurf(Mloc,Nloc),Vsurf(Mloc,Nloc),Ubott(Mloc,Nloc),Vbott(Mloc,Nloc))

        DO I_tmp=1,NumTracker
           READ(1,*)X_TRACK(I_tmp),Y_TRACK(I_tmp),TRACKER_START(I_tmp),&
               LAYER_012(I_tmp)
        ENDDO
      CLOSE(1)

     IF(.NOT.ALLOCATED(Sc)) THEN
       ALLOCATE(  &
                Sc(NumTracker), &
                S1(NumTracker), &
                S2(NumTracker), &
                S3(NumTracker), &
                nx1(NumTracker), &
                nx2(NumTracker), &
                nx3(NumTracker), &
                ny1(NumTracker), &
                ny2(NumTracker), &
                ny3(NumTracker), &
                Sc_pre(NumTracker), &
                S1_pre(NumTracker), &
                S2_pre(NumTracker), &
                S3_pre(NumTracker), &
                nx1_pre(NumTracker), &
                nx2_pre(NumTracker), &
                nx3_pre(NumTracker), &
                ny1_pre(NumTracker), &
                ny2_pre(NumTracker), &
                ny3_pre(NumTracker), &
                IN_CELL(NumTracker), &
                FOUND_IN_DOMAIN(NumTracker), &
                STOP_SEARCH(NumTracker), &
                X_TRACK_PRE(NumTracker), &
                Y_TRACK_PRE(NumTracker), &
                U_TRACER(NumTracker), &
                V_TRACER(NumTracker) &
               )
        IN_CELL = .FALSE.
        FOUND_IN_DOMAIN = 0
        STOP_SEARCH = .FALSE.
     ENDIF

!   print*,myid,x_track(1),y_track(1),x_track(2),y_track(2)

       CALL GET_XY_POSITION
       PLOT_COUNT_TRACKING = 0
     
END SUBROUTINE TRACER_INITIAL

SUBROUTINE GET_XY_POSITION
      IMPLICIT NONE
      REAL(SP) :: x1,x2,x3,y1,y2,y3,&
                  area1,area2,area3
      INTEGER :: I,I_start,I_end,J_start,J_end,search_count
      INTEGER :: II,JJ

! ---  areas of the four triangles, area will be negative
!      if an order is clockwise
!       Sc -- triangle 1,2,3
!       S1 -- triangle 2,3,c
!       S2 -- triangle 3,1,c
!       S3 -- triangle 1,2,c

! --- find the triangle includes the points in grid2

  
        DO I=1,NumTracker
          x1=X_TRACK(I)
          y1=Y_TRACK(I)
          
          search_count=0
300       CONTINUE  ! re-search
          search_count=search_count+1

         IF(.NOT.STOP_SEARCH(I)) THEN

          IF(IN_CELL(I))THEN
           I_start=MAX(nx1(I)-2,Ibeg)
           I_end=MIN(nx1(I)+2,Iend)
           J_start=MAX(ny2(I)-2,Jbeg)
           J_end=MIN(ny2(I)+2,Jend)
          ELSEIF(FOUND_IN_DOMAIN(I)<1)THEN
# if defined (PARALLEL)
! use Iend-1 is a wrong idea for parallelization,
! since theres no x or y connected between processor interfaces
           I_start=Ibeg
           if ( n_east .eq. MPI_PROC_NULL ) then
            I_end=Iend-1
           else
            I_end=Iend
           endif
           J_start=Jbeg
           if ( n_nrth .eq. MPI_PROC_NULL ) then
             J_end=Jend-1
           else
             J_end=Jend
           endif
# else
           I_start=Ibeg
           I_end=Iend-1
           J_start=Jbeg
           J_end=Jend-1
# endif

          ELSE  ! for processors dont contain tracers
           I_start=Ibeg
           I_end=Ibeg
           J_start=Jbeg
           J_end=Jbeg 
          ENDIF

          do JJ=J_start,J_end
          do II=I_start,I_end

            x2=Xco(ii+1)
            y2=Yco(jj)
            x3=Xco(ii)
            y3=Yco(jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=Xco(ii)
            y2=Yco(jj+1)
            x3=Xco(ii)
            y3=Yco(jj)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=Xco(ii)
            y2=Yco(jj)
            x3=Xco(ii+1)
            y3=Yco(jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              IN_CELL(I) = .TRUE.
              FOUND_IN_DOMAIN(I) = 1
              nx1(I)=ii
              ny1(I)=jj
              nx2(I)=ii+1
              ny2(I)=jj
              nx3(I)=ii
              ny3(I)=jj+1
              S1(I)=area1
              S2(I)=area2
              S3(I)=area3
               
              x1=Xco(nx1(I))
              y1=Yco(ny1(I))
              x2=Xco(nx2(I))
              y2=Yco(ny2(I))
              x3=Xco(nx3(I))
              y3=Yco(ny3(I))
              Sc(I)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 120
            endif

            x2=Xco(ii+1)
            y2=Yco(jj)
            x3=Xco(ii+1)
            y3=Yco(jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=Xco(ii+1)
            y2=Yco(jj+1)
            x3=Xco(ii)
            y3=Yco(jj+1)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=Xco(ii)
            y2=Yco(jj+1)
            x3=Xco(ii+1)
            y3=Yco(jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              IN_CELL(I) = .TRUE.
              FOUND_IN_DOMAIN(I) = 1
              nx1(I)=ii
              ny1(I)=jj+1
              nx2(I)=ii+1
              ny2(I)=jj
              nx3(I)=ii+1
              ny3(I)=jj+1
              S1(I)=area1
              S2(I)=area2
              S3(I)=area3
               
              x1=Xco(nx1(I))
              y1=Yco(ny1(I))
              x2=Xco(nx2(I))
              y2=Yco(ny2(I))
              x3=Xco(nx3(I))
              y3=Yco(ny3(I))
              Sc(I)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 120
            endif

          enddo
          enddo

            IN_CELL(I) = .FALSE.
            FOUND_IN_DOMAIN(I) = 0

120       continue


# if defined (PARALLEL)
!  gather found_in_domain
     call MPI_ALLREDUCE(FOUND_IN_DOMAIN(I),myint,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier) 
          FOUND_IN_DOMAIN(I) = myint  

          IF(FOUND_IN_DOMAIN(I)<1)THEN
!            WRITE(*,*)'myid',myid,'TRACER',I,'LOST, RE-SEARCH'
            IF(search_count>1) THEN
!              WRITE(*,*)'myid',myid,'TRACER',I,'LOST permanently'
              STOP_SEARCH(I)=.TRUE.
            ELSE
              STOP_SEARCH(I)=.FALSE.
              GOTO 300
            ENDIF
          ENDIF ! couldnt find
   
# else
          IF(FOUND_IN_DOMAIN(I)<1)THEN
!            WRITE(*,*)'TRACER',I,'LOST, RE-SEARCH'
            IF(search_count>1) THEN
!              WRITE(*,*)'TRACER',I,'LOST permanently'
              STOP_SEARCH(I)=.TRUE.
            ELSE
              STOP_SEARCH(I)=.FALSE.
              GOTO 300
            ENDIF
          ENDIF ! couldnt find
# endif

        ENDIF ! end stop search

        ENDDO ! end numtracker

! print*,myid,in_cell(2),Y_TRACK(2),STOP_SEARCH(2),FOUND_IN_DOMAIN(2)
!      CALL SWEXITMPI  
        
END SUBROUTINE GET_XY_POSITION

SUBROUTINE TRACK_XY
     IMPLICIT NONE
     INTEGER :: I
       X_TRACK_PRE = X_TRACK
       Y_TRACK_PRE = Y_TRACK
       nx1_pre=nx1
       ny1_pre=ny1
       nx2_pre=nx2
       ny2_pre=ny2
       nx3_pre=nx3
       ny3_pre=ny3
       S1_pre=S1
       S2_pre=S2
       S3_pre=S3
       Sc_pre=Sc

       CALL MOVE_TRACER

# if defined (PARALLEL)
       CALL BROADCAST_XY
# endif

       CALL GET_XY_POSITION

       DO I=1,NumTracker
        IF(IN_CELL(I))THEN
         IF(MASK(nx1(I),ny1(I))+MASK(nx2(I),ny2(I))+MASK(nx2(I),ny2(I))<1)THEN
          nx1(I)=nx1_pre(I)
          ny1(I)=ny1_pre(I)
          nx2(I)=nx2_pre(I)
          ny2(I)=ny2_pre(I)
          nx3(I)=nx3_pre(I)
          ny3(I)=ny3_pre(I)
          S1(I)=S1_pre(I)
          S2(I)=S2_pre(I)
          S3(I)=S3_pre(I)
          Sc(I)=Sc_pre(I)
          X_TRACK(I)=X_TRACK_PRE(I)
          Y_TRACK(I)=Y_TRACK_PRE(I)
         ENDIF
        ENDIF ! end in cell
       ENDDO

END SUBROUTINE TRACK_XY


# if defined (PARALLEL)
SUBROUTINE BROADCAST_XY
     IMPLICIT NONE
     INTEGER :: I,l
     INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
     REAL(SP),DIMENSION(NumberProcessor) :: xx
     LOGICAL,DIMENSION(NumberProcessor) :: yy
       REAL(SP),DIMENSION(NumTracker,NumberProcessor) :: X_TRACK_ALL,Y_TRACK_ALL
       LOGICAL,DIMENSION(NumTracker,NumberProcessor) :: IN_CELL_ALL

     call MPI_Gather(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
     call MPI_Gather(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)

! gather to myid=0
       DO I=1,NumTracker

! x
        call MPI_Gather(X_TRACK(I),1,MPI_SP,&
             xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)

!        if (I.eq.1) call MPI_Barrier(MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
               X_TRACK_ALL(I,l)= xx(l)
           enddo
        endif 

! y       
        call MPI_Gather(Y_TRACK(I),1,MPI_SP,&
             xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
               Y_TRACK_ALL(I,l)= xx(l)
           enddo
        endif 

! in_cell
        call MPI_Gather(IN_CELL(I),1,MPI_LOGICAL,&
             yy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
               IN_CELL_ALL(I,l)= yy(l)
           enddo
        endif 

        ENDDO

       DO I=1,NumTracker
         DO l=1,px*py
          IF(IN_CELL_ALL(I,l))THEN
             X_TRACK(I)=X_TRACK_ALL(I,l)
             Y_TRACK(I)=Y_TRACK_ALL(I,l)
          ENDIF
         ENDDO
       ENDDO        

! send x,y to all processors


       DO I=1,NumTracker
! send x
        if (myid.eq.0) then
           do l=1,px*py
              xx(l) =  X_TRACK(I)
           enddo
        endif
        call MPI_Scatter(xx,1,MPI_SP,&
             X_TRACK(I),1,MPI_SP,0,MPI_COMM_WORLD,ier)
! send y
        if (myid.eq.0) then
           do l=1,px*py
              xx(l) =  Y_TRACK(I)
           enddo
        endif
        call MPI_Scatter(xx,1,MPI_SP,&
             Y_TRACK(I),1,MPI_SP,0,MPI_COMM_WORLD,ier)
       ENDDO


END SUBROUTINE BROADCAST_XY
# endif

SUBROUTINE MOVE_TRACER
     IMPLICIT NONE
     INTEGER :: I

        DO I=1,NumTracker
           IF(IN_CELL(I).AND.TIME>TRACKER_START(I))THEN            
              IF(Sc(I)==0.0)THEN
                write(*,*) 'Sc is ZERO for TRACER', I
              ELSE
                IF(LAYER_012(I)==1)THEN  ! surface
                  U_TRACER(I)=(S1(I)*Usurf(nx1(I),ny1(I)) &
                        +S2(I)*Usurf(nx2(I),ny2(I)) &
                        +S3(I)*Usurf(nx3(I),ny3(I))) &
                        /Sc(I) 
                  V_TRACER(I)=(S1(I)*Vsurf(nx1(I),ny1(I)) &
                        +S2(I)*Vsurf(nx2(I),ny2(I)) &
                        +S3(I)*Vsurf(nx3(I),ny3(I))) &
                        /Sc(I) 
                ELSEIF(LAYER_012(I)==2)THEN ! bottom
                  U_TRACER(I)=(S1(I)*Ubott(nx1(I),ny1(I)) &
                        +S2(I)*Ubott(nx2(I),ny2(I)) &
                        +S3(I)*Ubott(nx3(I),ny3(I))) &
                        /Sc(I) 
                  V_TRACER(I)=(S1(I)*Vbott(nx1(I),ny1(I)) &
                        +S2(I)*Vbott(nx2(I),ny2(I)) &
                        +S3(I)*Vbott(nx3(I),ny3(I))) &
                        /Sc(I)
                ELSE        ! depth averaged
                  U_TRACER(I)=(S1(I)*U(nx1(I),ny1(I)) &
                        +S2(I)*U(nx2(I),ny2(I)) &
                        +S3(I)*U(nx3(I),ny3(I))) &
                        /Sc(I) 
                  V_TRACER(I)=(S1(I)*V(nx1(I),ny1(I)) &
                        +S2(I)*V(nx2(I),ny2(I)) &
                        +S3(I)*V(nx3(I),ny3(I))) &
                        /Sc(I)
                ENDIF
                X_TRACK(I)=X_TRACK(I)+U_TRACER(I)*DT
                Y_TRACK(I)=Y_TRACK(I)+V_TRACER(I)*DT
              ENDIF 
            ENDIF ! in_cell    
        ENDDO

END SUBROUTINE MOVE_TRACER 
 
SUBROUTINE OUTPUT_TRACKING
     IMPLICIT NONE
     LOGICAL, SAVE :: FirstCallTrackingOutput = .TRUE.
     CHARACTER(LEN=14)::FORMAT_LEN_TRACKING=' '
     INTEGER :: I

# if defined (PARALLEL)
    if(myid==0)then
# endif
     IF(FirstCallTrackingOutput)THEN
       FirstCallTrackingOutput = .FALSE.

! format length
        write(FORMAT_LEN_TRACKING(1:1),'(A1)') '('
        write(FORMAT_LEN_TRACKING(2:8),'(I7)') NumTracker+1
        write(FORMAT_LEN_TRACKING(9:13),'(A5)') 'E16.6'
        write(FORMAT_LEN_TRACKING(14:14),'(A1)') ')'

       OPEN(8,FILE=TRIM(RESULT_FOLDER)//'tk_x.txt')
       OPEN(9,FILE=TRIM(RESULT_FOLDER)//'tk_y.txt')
     ENDIF

     WRITE(8,FORMAT_LEN_TRACKING)TIME,(X_TRACK(I),I=1,NumTracker)
     WRITE(9,FORMAT_LEN_TRACKING)TIME,(Y_TRACK(I),I=1,NumTracker)
# if defined (PARALLEL)
    endif
# endif

END SUBROUTINE OUTPUT_TRACKING


END MODULE TRACER
#endif
! end tracking

