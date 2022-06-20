!------------------------------------------------------------------------------------
!
!      FILE misc.F
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
!    INDEX is subroutine to index for MPI
!
!    HISTORY: 
!    05/06/2010 Fengyan Shi
!    01/27/2016 Fengyan Shi, added ProcessorID for periodic bc
!
!-------------------------------------------------------------------------------------

SUBROUTINE INDEX
    USE GLOBAL
    IMPLICIT NONE

![--------------ykchoi 02/May/2017
# if defined (PARALLEL)
    INTEGER :: icp
    INTEGER :: myidi, myidj, irank
    INTEGER, ALLOCATABLE :: iistas(:), iiends(:), jjstas(:), jjends(:)
# endif
!--------------ykchoi 02/May/2017]

! TEMP
# if defined (PARALLEL)

    ALLOCATE(ProcessorID(PX,PY)) 

    NumberProcessor = px*py
    dims(1) = px
    dims(2) = py
    periods(1) = .false.
    periods(2) = .false.

! dont use periodic because our TRID solver is not based on
!  the cyclic topology  
!    IF(PERIODIC) periods(2) = .true.  

    coords(1) = 0
    coords(2) = 0

![--------------ykchoi 04/May/2017
   if( nprocs /= NumberProcessor) then
     if( myid==0 ) then
       print *, '======================================================='
       print *, '*** STOP :: Number of processors(',nprocs,') in MPIRUN /= Px*Py(',PX*PY,') ***'
       print *, '======================================================='
     endif
     call MPI_FINALIZE ( ier )
   endif   
!--------------ykchoi 04/May/2017]

    call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, &
         periods, reorder, comm2d, ier )
    call MPI_CART_COORDS( comm2d, myid, 2, coords, ier)

    npx = coords(1)
    npy = coords(2)

    call MPI_Cart_shift( comm2d, 0, 1, n_west, n_east, ier )
    call MPI_Cart_shift( comm2d, 1, 1, n_suth, n_nrth, ier )

    icp=0
    DO I=1,PX
    DO J=1,PY
      ProcessorID(I,J) = icp
      ![----ykchoi(04/May/2017) 
      if( myid == icp ) then      
        myidi = I-1; myidj = J-1
      endif
      !----ykchoi]
      icp=icp+1
    ENDDO
    ENDDO

! check
! print*,myid, n_west,n_east,n_suth,n_nrth,ProcessorID(1,1),ProcessorID(1,2)
!      print*,myid,ProcessorID(1,1),ProcessorID(2,1),ProcessorID(1,2),ProcessorID(2,2)
!      call MPI_FINALIZE ( ier )

# else
  px=1
  py=1
# endif

! now for serial code
![--------------ykchoi 04/May/2017
    !Mloc=Mglob/px+2*Nghost
    !Nloc=Nglob/py+2*Nghost

# if defined (PARALLEL)
    call grid_range_per_procs( 1, Mglob, px, myidi, iista, iiend )
    call grid_range_per_procs( 1, Nglob, py, myidj, jjsta, jjend )

    Mloc = (iiend-iista+1) + 2*Nghost
    Nloc = (jjend-jjsta+1) + 2*Nghost 
    
    !===========for check
    allocate( iistas(nprocs), iiends(nprocs), jjstas(nprocs), jjends(nprocs) ) 
    call mpi_gather( iista, 1, mpi_integer, iistas, 1, mpi_integer, &
                     0, mpi_comm_world, ier )
    call mpi_gather( iiend, 1, mpi_integer, iiends, 1, mpi_integer, &
                     0, mpi_comm_world, ier )
    call mpi_gather( jjsta, 1, mpi_integer, jjstas, 1, mpi_integer, &
                     0, mpi_comm_world, ier )
    call mpi_gather( jjend, 1, mpi_integer, jjends, 1, mpi_integer, &
                     0, mpi_comm_world, ier )

    if( myid == 0 ) then
      open(500, file='Grid_Range.out', status='unknown')
      write(500,'(A)')'Core#  istart#  iend#  jstart#  jend#'
      do irank=0,nprocs-1
         write(500,'(5(i6,1x))') irank, iistas(irank+1), iiends(irank+1), &
                                        jjstas(irank+1), jjends(irank+1)
      enddo
      close(500)
    endif
    deallocate( iistas, iiends, jjstas, jjends ) 
    !=============================

# else
    Mloc = Mglob + 2*Nghost
    Nloc = Nglob + 2*Nghost

# endif
!--------------ykchoi 04/May/2017]

    Mloc1=Mloc+1
    Nloc1=Nloc+1

    Ibeg=Nghost+1
    Iend=Mloc-Nghost
    Iend1=Mloc1-Nghost
    Jbeg=Nghost+1
    Jend=Nloc-Nghost
    Jend1=Nloc1-Nghost

END SUBROUTINE INDEX

![--------------ykchoi 04/May/2017
subroutine grid_range_per_procs( n1, n2, nprocs, myid, stag, endg )
!-------------------------------------------------------------------------------------
!
!    Grid_range_per_procs is subroutine to get iiend etc
!
!    HISTORY: 
!       05/10/2017 Fengyan Shi copied from Chois codes
!
!-------------------------------------------------------------------------------------

    implicit none
    integer, intent(in) :: n1, n2, nprocs, myid
    integer, intent(out) :: stag, endg

    integer :: inum1, inum2

    inum1 = int( ( n2 - n1 + 1 )/nprocs )
    inum2 = mod( n2-n1+1, nprocs )

    stag = myid*inum1 + n1 + min(myid, inum2)
    endg = stag + inum1 - 1

    if( inum2 > myid ) endg=endg+1
endsubroutine grid_range_per_procs
!--------------ykchoi 04/May/2017]

!-------------------------------------------------------------------------------------
!
!    ESTIMATE_DT is subroutine evaluate dt based in CFL
!
!    HISTORY: 
!       05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------

SUBROUTINE ESTIMATE_DT(M,N,DX,DY,U,V,H,MinDepthFrc,DT,CFL,TIME)
     USE PARAM
     USE GLOBAL, ONLY : DT_fixed, FIXED_DT
# if defined (PARALLEL)
     USE GLOBAL, ONLY : ier
# endif
     IMPLICIT NONE
     INTEGER,INTENT(IN)::M,N
# if defined (PARALLEL)
     REAL(SP) :: myvar
# endif

# if defined (CARTESIAN)
     REAL(SP),INTENT(IN)::DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY
# endif
     REAL(SP),INTENT(IN),DIMENSION(M,N)::U,V,H
     REAL(SP),INTENT(IN)::CFL,MinDepthFrc
     REAL(SP),INTENT(OUT)::DT
     REAL(SP),INTENT(INOUT)::TIME
     REAL(SP) :: DT_tmp

     TMP3=LARGE
     DO J=1,N
     DO I=1,M
! x direction
      TMP1=ABS(U(I,J))+SQRT(GRAV*MAX(H(I,J),MinDepthFrc))
      IF(TMP1<SMALL)THEN
# if defined (CARTESIAN)
       TMP2=DX/SMALL
# else
       TMP2=DX(I,J)/SMALL
# endif
      ELSE
# if defined (CARTESIAN)
       TMP2=DX/TMP1
# else
       TMP2=DX(I,J)/TMP1
# endif
      ENDIF
      IF(TMP2<TMP3)TMP3=TMP2
! y direction
      TMP1=ABS(V(I,J))+SQRT(GRAV*MAX(H(I,J),MinDepthFrc))
      IF(TMP1<SMALL)THEN
# if defined (CARTESIAN)
       TMP2=DY/SMALL
# else
       TMP2=DY(I,J)/SMALL
# endif
      ELSE
# if defined (CARTESIAN)
       TMP2=DY/TMP1
# else
       TMP2=DY(I,J)/TMP1
# endif
      ENDIF
      IF(TMP2<TMP3)TMP3=TMP2      
     ENDDO
     ENDDO
# if defined (PARALLEL)
     call MPI_ALLREDUCE (TMP3,myvar,1,MPI_SP,MPI_MIN,&
          MPI_COMM_WORLD,ier)
     TMP3 = myvar
# endif
     DT_tmp=CFL*TMP3

    IF(FIXED_DT)THEN

      DT = DT_fixed
      DO WHILE (DT > DT_tmp)
        DT=DT/2.0_SP
      ENDDO
    ELSE
     DT = DT_tmp
    ENDIF

! TEMP
     TIME=TIME+DT

END SUBROUTINE ESTIMATE_DT

!-------------------------------------------------------------------------------------
!
!    MAX_MIN_PROPERTY is subroutine to calculate Max and Min properties 
!        based on Chen et al., 2004
!
!    HISTORY: 
!        02/10/2016 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE MAX_MIN_PROPERTY
    USE GLOBAL
# if defined (METEO)
    USE METEO_MODULE
# endif
    IMPLICIT NONE
    REAL(SP) :: maxv
# if !defined (CARTESIAN)
    REAL(SP) :: omega_0
# endif

      DO J=1,Nloc
      DO I=1,Mloc
# if defined (METEO)
       IF(OUT_Hmax.OR.WindForce)THEN
# else
       IF(OUT_Hmax)THEN
# endif
        IF(MASK(I,J).GT.0)THEN
        IF(Eta(I,J).GT.HeightMax(I,J)) HeightMax(I,J)=Eta(I,J)
        ENDIF
       ENDIF

       IF(OUT_Hmin)THEN
        IF(MASK(I,J).GT.0)THEN
        IF(Eta(I,J).LT.HeightMin(I,J)) HeightMin(I,J)=Eta(I,J)
        ENDIF
       ENDIF

       IF(OUT_Umax)THEN
        IF(MASK(I,J).GT.0)THEN
          maxv=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
          IF(maxV.GT.VelocityMax(I,J)) VelocityMax(I,J)=maxV
        ENDIF
       ENDIF

       IF(OUT_MFmax)THEN
        IF(MASK(I,J).GT.0)THEN
          maxv=(U(I,J)*U(I,J)+V(I,J)*V(I,J))*(H(I,J))
          IF(maxv.GT.MomentumFluxMax(I,J)) MomentumFluxMax(I,J)=maxv
        ENDIF
       ENDIF

       !Lauren Schambach 3/2/2020 Matrix of Arrival Time
       IF(OUT_Time)THEN
         IF(MASK(I,J).GT.0)THEN !Check if wet
           IF(ARRTIME(I,J).EQ.0)THEN !Only record time if a value doesnt already exist
             IF(ABS(Eta(I,J)).GT.ArrTimeMin) ARRTIME(I,J) = TIME !If eta is greater than threshold, record arrival time
           ENDIF
         ENDIF
       ENDIF

      ENDDO
      ENDDO

      IF(OUT_VORmax) THEN
# if defined (CARTESIAN)
   ! the vorticity has been calculated in dispersion.F
# else
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,V,Vx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,U,Uy)
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
        omega_0=ABS(Vx(I,J)-Uy(I,J))
        IF(omega_0.GT.VorticityMax(I,J)) THEN
          VorticityMax(I,J)=omega_0
        ENDIF
     ENDDO
     ENDDO
# endif
      ENDIF ! max vorticity

# if defined (PARALLEL)
      IF(OUT_VORmax) THEN
       CALL phi_exch(VorticityMax)
      ENDIF
# endif   


END SUBROUTINE MAX_MIN_PROPERTY

!-------------------------------------------------------------------------------------
!
!    CHECK_BLOWUP is subroutine to check numerical stability 
!
!    HISTORY: 
!        01/23/2015 Young-Kwang Choi
!        02/15/2016 Fengyan Shi, added the threshold to 100*max_depth
!
!-------------------------------------------------------------------------------------
SUBROUTINE CHECK_BLOWUP
    USE GLOBAL
    IMPLICIT NONE
![ykchoi 15.01.23.
     REAL(SP) :: MaxAbsEta
# if defined (PARALLEL)
     REAL(SP)::myvar_tmp
# endif

      MaxAbsEta=MAXVAL( abs(Eta(Ibeg:Iend,Jbeg:Jend)) )
# if defined (PARALLEL)
      CALL MPI_ALLREDUCE(MaxAbsEta,myvar_tmp,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
      MaxAbsEta = myvar_tmp
# endif

      if (MaxAbsEta > EtaBlowVal) then
# if defined (PARALLEL)
         if (myid.eq.0) then
# endif
            WRITE(*,*) "========================================="
              WRITE(*,*) "BlowUp Time, MaxAbsEta=", Time, MaxAbsEta
            WRITE(*,*) "========================================="
# if defined (PARALLEL)
         endif
# endif
         ICOUNT=99998;
         CALL PREVIEW

# if defined (PARALLEL)
     ! 2019/09/04 mayhl: Switched MPI_FINALIZE to MPI_ABORT 
     !call MPI_FINALIZE ( ier )
     call MPI_ABORT( MPI_COMM_WORLD , 9 , ier )
# else
     STOP
# endif

      endif


END SUBROUTINE CHECK_BLOWUP

!-------------------------------------------------------------------------------------
!
!   wall_time_secs is used to calculate current wall time
!
!   HISTORY: 
!   Gangfeng Ma, 09/12/2011
!
!-------------------------------------------------------------------------------------
    SUBROUTINE wall_time_secs(tcurrent)
    IMPLICIT NONE
    INTEGER, dimension(8) :: walltime
    real, INTENT(OUT) :: tcurrent
    real :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &                                                                                   
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
    end SUBROUTINE wall_time_secs

#if defined (FILTERING)
SUBROUTINE SHAPIRO_3(M,N,Nghost,MASK,VAR_IN_OUT)  
  USE PARAM,ONLY : SP
  IMPLICIT NONE
  REAL(SP),PARAMETER :: a1=0.68750_SP    ! 44/64
  REAL(SP),PARAMETER :: a2=0.234375_SP   ! 15/64
  REAL(SP),PARAMETER :: a3=-0.09375_SP   ! -6/64
  REAL(SP),PARAMETER :: a4=0.015625_SP        ! 1/64
  INTEGER,          INTENT(IN) :: M,N,Nghost
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: VAR_IN_OUT
  REAL(SP),DIMENSION(M,N) :: VAR
  INTEGER, DIMENSION(M,N),INTENT(IN) :: MASK
  REAL(SP) :: p0,pp1,pp2,pp3,pn1,pn2,pn3
  INTEGER :: I,J


! x direction
  DO J=1,N
  DO I=1+Nghost,M-Nghost
    IF(MASK(I,J)==1) THEN
      p0=VAR_IN_OUT(I,J)
      pp1=VAR_IN_OUT(I+1,J)*MASK(I+1,J)+(1.0-MASK(I+1,J))*p0
      pp2=VAR_IN_OUT(I+2,J)*MASK(I+2,J)+(1.0-MASK(I+2,J))*p0
      pp3=VAR_IN_OUT(I+3,J)*MASK(I+3,J)+(1.0-MASK(I+3,J))*p0
      pn1=VAR_IN_OUT(I-1,J)*MASK(I-1,J)+(1.0-MASK(I-1,J))*p0
      pn2=VAR_IN_OUT(I-2,J)*MASK(I-2,J)+(1.0-MASK(I-2,J))*p0
      pn3=VAR_IN_OUT(I-3,J)*MASK(I-3,J)+(1.0-MASK(I-3,J))*p0

      VAR(I,J) = a1*p0+a2*(pp1+pn1)+a3*(pp2+pn2)+a4*(pp3+pn3)

    ENDIF
  ENDDO
  ENDDO


! y direction

  DO J=1+Nghost,N-Nghost
  DO I=1,M
    IF(MASK(I,J)==1) THEN
      p0=VAR(I,J)
      pp1=VAR(I,J+1)*MASK(I,J+1)+(1.0-MASK(I,J+1))*p0
      pp2=VAR(I,J+2)*MASK(I,J+2)+(1.0-MASK(I,J+2))*p0
      pp3=VAR(I,J+3)*MASK(I,J+3)+(1.0-MASK(I,J+3))*p0
      pn1=VAR(I,J-1)*MASK(I,J-1)+(1.0-MASK(I,J-1))*p0
      pn2=VAR(I,J-2)*MASK(I,J-2)+(1.0-MASK(I,J-2))*p0
      pn3=VAR(I,J-3)*MASK(I,J-3)+(1.0-MASK(I,J-3))*p0

      VAR_IN_OUT(I,J) = a1*p0+a2*(pp1+pn1)+a3*(pp2+pn2)+a4*(pp3+pn3)

    ENDIF
  ENDDO
  ENDDO


END SUBROUTINE SHAPIRO_3
# endif

