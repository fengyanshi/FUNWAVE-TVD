!------------------------------------------------------------------------------------
!
!      FILE tridiagnal.F
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
!    TRIDy_periodic is subroutine to
!      solve the periodic, quasi-tridiagonal problem in the y direction
!
!      This subroutine is based on the Sherman-Morrison algorithm as described in
!      Thomas, J. W., 1995, Numerical Partial Differential Equations, Springer,
!      Section 5.6.1
!
!    HISTORY:
!      03/07/2014  Jim Kirby
!      02/29/2016  Fengyan Shi, finished and tested 
!
!-----------------------------------------------------------------------------------

    SUBROUTINE TRIDy_periodic ( a, c, d, f)
    USE PARAM

# if defined (PARALLEL)
    USE GLOBAL, ONLY : ier,myid,PX,PY,&
                       ProcessorID
    USE MPI
# endif
    USE GLOBAL, ONLY : Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    IMPLICIT NONE

# if defined (PARALLEL)
    INTEGER,DIMENSION(1) :: req
    INTEGER :: len,l,II
!    INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    REAL(SP),DIMENSION(Mloc) :: xx
    REAL(SP), dimension(Mloc) :: y2_end,y1_end
# endif

    REAL(SP), dimension(Mloc,Nloc) :: a, c, d, f
    REAL(SP), dimension(Mloc,Nloc) :: c1, y1, y2
    REAL(SP), dimension(Mloc) :: beta, a_beg,c_end
    INTEGER :: I,J

!
!  off diagonal corner values are stored in a(i,Jbeg) and c(i, Jend)
!
!  First, construct the diagonalized matrix including renormalization of the first and last rows to retain 
!  normalization of main diagonal at 1. 
!

# if defined (PARALLEL)
    IF (PY>1)THEN

      len=Mloc

! south
    DO II = 1,PX

    if(myid==ProcessorID(II,1)) then
      do i = 1,Mloc
        a_beg(i)=a(i,Jbeg)
      enddo

! send from master
        call MPI_SEND(a_beg,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,ier)

    endif ! end myid

        if(myid==ProcessorID(II,PY))then
        call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
          DO I=1,Mloc
            a_beg(I)=xx(I)
          ENDDO
        endif

! north
    if(myid==ProcessorID(II,PY)) then
      do i = 1,Mloc
        c_end(i)=c(i,Jend)
      enddo
! send from master
        call MPI_SEND(c_end,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
    endif ! end myid

        if(myid==ProcessorID(II,1))then
        call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
          DO I=1,Mloc
            c_end(I)=xx(I)
          ENDDO
        endif
 
    ENDDO  ! end PX

    ELSE  ! PY = 1
      DO I = Ibeg,Iend
       A_beg(I) = a(I,Jbeg)
       C_end(I) = c(I,Jend)
      ENDDO
    ENDIF

# else
    DO I = Ibeg,Iend
     A_beg(I) = a(I,Jbeg)
     C_end(I) = c(I,Jend)
    ENDDO
# endif

# if defined (PARALLEL)
   DO II=1,PX
    if(myid == ProcessorID(II,1)) then
      do i= Ibeg, Iend
            c(i,Jbeg)=c(i,Jbeg)/(1+c_end(I))
            d(i,Jbeg)=d(i,Jbeg)/(1+c_end(I))
      enddo
    endif
    if(myid == ProcessorID(II,PY)) then
      do i= Ibeg, Iend
            a(i,Jend)=a(i,Jend)/(1+a_beg(I))
            d(i,Jend)=d(i,Jend)/(1+a_beg(I))
      enddo
    endif
   ENDDO  ! end PX

# else
    do i= Ibeg, Iend
            c(i,Jbeg)=c(i,Jbeg)/(1+c_end(I))
            d(i,Jbeg)=d(i,Jbeg)/(1+c_end(I))
            a(i,Jend)=a(i,Jend)/(1+a_beg(I))
            d(i,Jend)=d(i,Jend)/(1+a_beg(I))
    enddo
# endif

!
!  Diagonal matrix now stored in a,c.  c and d will get trashed in first call to TRIDy.
!  We need c again but not d
!
      c1=c
!
!  Solution pass one:  Solve B y1 = d
!
# if defined(PARALLEL)
        CALL TRIDy ( a, c1, d, y1)
# else
      call TRIDy_ser( a, c1, d, y1)
# endif
!
!  Now, the RHS for the second problem, which is ones in the first and last positions (also renormalized)
!
# if defined (PARALLEL)

!  initialize d
    DO J=Jbeg, Jend
    DO I=Ibeg, Iend
       d(i,j)=ZERO
    ENDDO
    ENDDO

    DO II=1,PX
        if(myid==ProcessorID(II,1))then
          do i=Ibeg, Iend               
                  d(i,Jbeg)=1.0_SP/(1.0_SP+c_end(I))
          enddo
        endif

        if(myid==ProcessorID(II,PY))then
          do i=Ibeg, Iend
              d(i,Jend)=-1.0_SP/(1.0_SP+a_beg(I))
          enddo
      endif
    ENDDO ! end PX
# else 
      do i=Ibeg, Iend
            d(i,Jbeg)=1.0_SP/(1.0_SP+c_end(I))
            do j = Jbeg+1, Jend-1
                  d(i,j)=ZERO
            enddo
            d(i,Jend)=-1.0_SP/(1.0_SP+a_beg(I))
      enddo
# endif

!
!  Solution pass two:  Solve By2=w
!

# if defined (PARALLEL)
        CALL TRIDy ( a, c, d, y2)
# else
      call TRIDy_ser( a, c, d, y2)
# endif
!
!  Compute scale factor
!

# if defined (PARALLEL)

     IF(PY>1)THEN

! transfer y2 y1 end
       DO II=1,PX
! y2
    if(myid==ProcessorID(II,PY)) then
      do i = 1,Mloc
        y2_end(i)=y2(i,Jend)
      enddo
! send from master
        call MPI_SEND(y2_end,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
    endif ! end myid

        if(myid==ProcessorID(II,1))then
        call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
          DO I=1,Mloc
            y2_end(I)=xx(I)
          ENDDO
        endif
! y1
    if(myid==ProcessorID(II,PY)) then
      do i = 1,Mloc
        y1_end(i)=y1(i,Jend)
      enddo
! send from master
        call MPI_SEND(y1_end,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
    endif ! end myid

        if(myid==ProcessorID(II,1))then
        call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
          DO I=1,Mloc
            y1_end(I)=xx(I)
          ENDDO
        endif

! calculate beta
        if(myid==ProcessorID(II,1))then
         DO I=Ibeg,Iend
            beta(i) = ( c_end(I)*y1(i,Jbeg) - a_beg(I)*y1_end(I)) / &
                  ( 1.0_SP - ( c_end(I)*y2(i,Jbeg) - a_beg(I)*y2_end(I) ) )

         ENDDO
        endif ! end myid
! transfer to other II
        if(myid==ProcessorID(II,1))then
          DO l = 2,PY
            call MPI_SEND(beta,len,MPI_SP,ProcessorID(II,l),0,MPI_COMM_WORLD,ier)
          ENDDO
        endif  
        DO l = 2,PY
          if(myid==ProcessorID(II,l))then
           call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req(1),ier)
           call MPI_WAITALL(1,req,status,ier)
           DO I=1,Mloc
             beta(I) = xx(I)
           ENDDO
          endif
        ENDDO      ! end l

       ENDDO ! end PX
     
     ELSE  ! PY=1

      do i=Ibeg, Iend
            beta(i) = ( c_end(I)*y1(i,Jbeg) - a_beg(I)*y1(i,Jend)) / &
                  ( 1.0_SP - ( c_end(I)*y2(i,Jbeg) - a_beg(I)*y2(i,Jend) ) )
        enddo

     ENDIF ! if PY>1 
                 
# else
      do i=Ibeg, Iend
            beta(i) = ( c_end(I)*y1(i,Jbeg) - a_beg(I)*y1(i,Jend)) / &
                  ( 1.0_SP - ( c_end(I)*y2(i,Jbeg) - a_beg(I)*y2(i,Jend) ) )
        enddo
# endif
      do i=Ibeg, Iend        
            do j = Jbeg, Jend
                  f(i,j)=y1(i,j)+beta(i)*y2(i,j)
            enddo
      enddo

      RETURN

      END SUBROUTINE TRIDy_periodic

# ifndef PARALLEL
!-----------------------------------------------------------------------------------
!
!    TRIDy_ser is serial version of the subroutine to
!      solve tridiagonal problem in the y direction based on Thomas algorithm. 
!
!    HISTORY:
!    03/07/2014  Jim Kirby
!    02/29/2016  Fengyan Shi, finished and tested 
!
!-----------------------------------------------------------------------------------
SUBROUTINE TRIDy_ser(A,C,D,Z)
        USE PARAM
        USE GLOBAL, ONLY : Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
        IMPLICIT NONE
        INTEGER :: II,JJ
        REAL(SP),DIMENSION(Mloc,Nloc),INTENT(INOUT) :: A,C,D
        REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: Z

      DO II = Ibeg,Iend
        DO JJ=Jbeg+1,Jend
          IF(A(II,JJ).NE.ZERO)THEN
            C(II,JJ)=C(II,JJ)/A(II,JJ)/(1.0_SP/A(II,JJ)-C(II,JJ-1))
            D(II,JJ)=(D(II,JJ)/A(II,JJ)-D(II,JJ-1))/(1.0_SP/A(II,JJ)-C(II,JJ-1))
          ENDIF
        ENDDO

        Z(II,Jend)=D(II,Jend)

        DO JJ=Jend-1,Jbeg,-1
          Z(II,JJ)=D(II,JJ)-C(II,JJ)*Z(II,JJ+1)
        ENDDO
      ENDDO  ! end II

END SUBROUTINE TRIDy_ser
# endif

# ifndef PARALLEL
!-----------------------------------------------------------------------------------
!
!    TRIDx_ser is serial version of the subroutine to
!      solve tridiagonal problem in the x direction based on Thomas algorithm. 
!
!    HISTORY:
!      02/11/2010  Fengyan Shi
!      02/29/2016  Fengyan Shi, modified to 2D arrays
!
!-----------------------------------------------------------------------------------
SUBROUTINE TRIDx_ser(A,C,D,Z)
        USE PARAM
        USE GLOBAL, ONLY : Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
        IMPLICIT NONE
        INTEGER :: II,JJ
        REAL(SP),DIMENSION(Mloc,Nloc),INTENT(INOUT) :: A,C,D
        REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: Z

       DO JJ=Jbeg,Jend
        DO II=Ibeg+1,Iend
          IF(A(II,JJ).NE.ZERO)THEN
            C(II,JJ)=C(II,JJ)/A(II,JJ)/(1.0_SP/A(II,JJ)-C(II-1,JJ))
            D(II,JJ)=(D(II,JJ)/A(II,JJ)-D(II-1,JJ))/(1.0_SP/A(II,JJ)-C(II-1,JJ))
          ENDIF
        ENDDO

        Z(Iend,JJ)=D(Iend,JJ)

        DO II=Iend-1,Ibeg,-1
          Z(II,JJ)=D(II,JJ)-C(II,JJ)*Z(II+1,JJ)
        ENDDO
       ENDDO

END SUBROUTINE TRIDx_ser
# endif

# if defined (PARALLEL)
!---------------------------------------------------------------------------------
!
!    TRIDx is parallel version of the subroutine to
!
!      solve tridiagonal problem in the x direction based on Thomas algorithm. 
!
!    HISTORY:
!      02/14/2010  Jeff Harris
!      05/01/2011  Fengyan Shi, tested the TVD code
!
!-----------------------------------------------------------------------------------
SUBROUTINE TRIDx ( a, c, d, f)
        USE PARAM
        USE GLOBAL, ONLY : n_east, n_west, comm2d, ier, Mloc, Nloc, &
                           Ibeg,Iend,Jbeg,Jend
        IMPLICIT NONE

        REAL(SP), dimension(Mloc,Nloc) :: a, c, d, f
        REAL(SP), dimension(Nloc,2) :: rmsg, smsg
        integer status(MPI_STATUS_SIZE), req
        integer :: i,j
        
! forward sweep

        if ( n_west .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Nloc, MPI_SP,&
                n_west, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           do j = Jbeg, Jend
           if (a(Ibeg,j).ne.ZERO) then
              c(Ibeg,j) = c(Ibeg,j) / a(Ibeg,j) / (1.0_SP/a(Ibeg,j) - rmsg(j,2))
              d(Ibeg,j) = (d(Ibeg,j) / a(Ibeg,j) - rmsg(j,1)) &
                   / (1.0_SP/a(Ibeg,j) - rmsg(j,2))
           endif
           enddo
        endif   

        do j = Jbeg, Jend
        do i = Ibeg+1, Iend
           if (a(i,j).ne.ZERO) then
              c(i,j) = c(i,j) / a(i,j) / (1.0_SP/a(i,j) - c(i-1,j))
              d(i,j) = (d(i,j) / a(i,j) - d(i-1,j)) &
                   / (1.0_SP/a(i,j) - c(i-1,j))
           endif
        enddo
        enddo

        if ( n_east .ne. MPI_PROC_NULL ) then
           do j = Jbeg, Jend
              smsg(j,1) = d(Iend,j)
              smsg(j,2) = c(Iend,j)
           enddo
           call MPI_ISEND( smsg, 2*Nloc, MPI_SP,&
                n_east, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

! back substitution

        if ( n_east .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Nloc, MPI_SP,&
                n_east, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           do j = Jbeg, Jend
              f(Iend,j) = d(Iend,j) - c(Iend,j) * rmsg(j,1)
           enddo
        else
           do j = Jbeg, Jend
              f(Iend,j) = d(Iend,j)
           enddo
        endif   

        do j = Jbeg, Jend
        do i = Iend-1, Ibeg, -1
           f(i,j) = d(i,j) - c(i,j) * f(i+1,j)
        enddo
        enddo

        if ( n_west .ne. MPI_PROC_NULL ) then
           do j = Jbeg, Jend
              smsg(j,1) = f(Ibeg,j)
           enddo
           call MPI_ISEND( smsg, 2*Nloc, MPI_SP,&
                n_west, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

END SUBROUTINE TRIDx
# endif

# if defined (PARALLEL)
!---------------------------------------------------------------------------------
!
!    TRIDx is parallel version of the subroutine to
!      solve tridiagonal problem in the y direction based on Thomas algorithm. 
!
!    HISTORY:
!      02/14/2010  Jeff Harris
!      05/01/2011  Fengyan Shi, tested the TVD code
!
!-----------------------------------------------------------------------------------
SUBROUTINE TRIDy ( a, c, d, f)
        USE PARAM
        USE GLOBAL, ONLY : n_suth, n_nrth, comm2d, ier, Mloc, Nloc, &
                           Ibeg,Iend,Jbeg,Jend
        IMPLICIT NONE 

        REAL(SP), dimension(Mloc,Nloc) :: a, c, d, f
        REAL(SP), dimension(Mloc,2) :: rmsg, smsg
        integer status(MPI_STATUS_SIZE), req
        integer :: i, j
        
! forward sweep

        if ( n_suth .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Mloc, MPI_SP,&
                n_suth, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           do i = Ibeg, Iend
           if (a(i,Jbeg).ne.ZERO) then
              c(i,Jbeg) = c(i,Jbeg) / a(i,Jbeg) / (1.0_SP/a(i,Jbeg) - rmsg(i,2))
              d(i,Jbeg) = (d(i,Jbeg) / a(i,Jbeg) - rmsg(i,1)) &
                   / (1.0_SP/a(i,Jbeg) - rmsg(i,2))
           endif
           enddo
        endif   

        do j = Jbeg+1, Jend
        do i = Ibeg, Iend
           if (a(i,j).ne.ZERO) then
              c(i,j) = c(i,j) / a(i,j) / (1.0_SP/a(i,j) - c(i,j-1))
              d(i,j) = (d(i,j) / a(i,j) - d(i,j-1)) &
                   / (1.0_SP/a(i,j) - c(i,j-1))
           endif
        enddo
        enddo

        if ( n_nrth .ne. MPI_PROC_NULL ) then
           do i = Ibeg, Iend
              smsg(i,1) = d(i,Jend)
              smsg(i,2) = c(i,Jend)
           enddo
           call MPI_ISEND( smsg, 2*Mloc, MPI_SP,&
                n_nrth, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

! back substitution

        if ( n_nrth .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Mloc, MPI_SP,&
                n_nrth, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           do i = Ibeg, Iend
              f(i,Jend) = d(i,Jend) - c(i,Jend) * rmsg(i,1)
           enddo
        else
           do i = Ibeg, Iend
              f(i,Jend) = d(i,Jend)
           enddo
        endif   

        do j = Jend-1, Jbeg, -1
        do i = Ibeg, Iend
           f(i,j) = d(i,j) - c(i,j) * f(i,j+1)
        enddo
        enddo

        if ( n_suth .ne. MPI_PROC_NULL ) then
           do i = Ibeg, Iend
              smsg(i,1) = f(i,Jbeg)
           enddo
           call MPI_ISEND( smsg, 2*Mloc, MPI_SP,&
                n_suth, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

END SUBROUTINE TRIDy
# endif
