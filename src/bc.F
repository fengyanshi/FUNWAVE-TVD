!------------------------------------------------------------------------------------
!
!      FILE bc.F
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
!    BOUNDARY_CONDITION is subroutine to provide boundary conditions at edges of domain
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!      04/05/2011  Jeff Harris, corrected bugs in do-loop
!
! -----------------------------------------------------------------------------------
#define nop(x) associate( x => x ); end associate

SUBROUTINE BOUNDARY_CONDITION
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::Xi,Deps
     INTEGER :: I,J


! four sides of computational domain

# if defined (PARALLEL)
        if ( n_west .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
   IF(IN_DOMAIN_WEST)THEN
     DO J=Jbeg,Kstart_WEST-1
      P(Ibeg,J)=ZERO
      Xi=EtaRxR(Ibeg,J)
      Deps=Depthx(Ibeg,J)
      Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Ibeg,J)=ZERO
      ENDDO
     DO J=Kend_WEST+1,Jend
      P(Ibeg,J)=ZERO
      Xi=EtaRxR(Ibeg,J)
      Deps=Depthx(Ibeg,J)
      Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Ibeg,J)=ZERO
      ENDDO
   ENDIF
# else
   IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
     ! do nothing
   ELSE
     DO J=Jbeg,Jend
      P(Ibeg,J)=ZERO
      Xi=EtaRxR(Ibeg,J)
      Deps=Depthx(Ibeg,J)
      Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Ibeg,J)=ZERO
      ENDDO
   ENDIF ! left bc wavemaker
# endif 

# if defined (PARALLEL)
      endif
# endif


# if defined (PARALLEL)
        if ( n_east .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
   IF(IN_DOMAIN_EAST)THEN
     DO J=Jbeg,Kstart_EAST-1
      P(Iend1,J)=ZERO
      Xi=EtaRxL(Iend1,J)
      Deps=Depthx(Iend1,J)
      Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Iend1,J)=ZERO
     ENDDO
     DO J=Kend_EAST+1,Jend
      P(Iend1,J)=ZERO
      Xi=EtaRxL(Iend1,J)
      Deps=Depthx(Iend1,J)
      Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Iend1,J)=ZERO
     ENDDO
   ENDIF
# else
     DO J=Jbeg,Jend
      P(Iend1,J)=ZERO
      Xi=EtaRxL(Iend1,J)
      Deps=Depthx(Iend1,J)
      Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Iend1,J)=ZERO
     ENDDO
# endif 

# if defined (PARALLEL)
      endif
# endif

# if defined (CARTESIAN)
! y direction
   IF(PERIODIC)THEN
!   do nothing
   ELSE
# endif

# if defined (PARALLEL)
      if ( n_suth .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
   IF(IN_DOMAIN_SOUTH)THEN
     DO I=Ibeg,Kstart_SOUTH-1
      Q(I,Jbeg)=ZERO
      Fy(I,Jbeg)=ZERO
      Xi=EtaRyR(I,Jbeg)
      Deps=Depthy(I,Jbeg)
      Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO
     DO I=Kend_SOUTH+1,Iend
      Q(I,Jbeg)=ZERO
      Fy(I,Jbeg)=ZERO
      Xi=EtaRyR(I,Jbeg)
      Deps=Depthy(I,Jbeg)
      Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO
   ENDIF
# else
     DO I=Ibeg,Iend
      Q(I,Jbeg)=ZERO
      Fy(I,Jbeg)=ZERO
      Xi=EtaRyR(I,Jbeg)
      Deps=Depthy(I,Jbeg)
      Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      ENDDO
# endif  

# if defined (PARALLEL)
      endif
# endif
# if defined (PARALLEL)
      if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
   IF(IN_DOMAIN_NORTH)THEN
     DO I=Ibeg,Kstart_NORTH-1
      Q(I,Jend1)=ZERO
      Fy(I,Jend1)=ZERO
      Xi=EtaRyL(I,Jend1)
      Deps=Depthy(I,Jend1)
      Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO
     DO I=Kend_NORTH+1,Iend
      Q(I,Jend1)=ZERO
      Fy(I,Jend1)=ZERO
      Xi=EtaRyL(I,Jend1)
      Deps=Depthy(I,Jend1)
      Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO
   ENDIF
# else
     DO I=Ibeg,Iend
      Q(I,Jend1)=ZERO
      Fy(I,Jend1)=ZERO
      Xi=EtaRyL(I,Jend1)
      Deps=Depthy(I,Jend1)
      Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO
# endif 
# if defined (PARALLEL)
     endif
# endif

# if defined (CARTESIAN)
    ENDIF
# endif

! mask points
! Jeff pointed out the loop should be Jbeg-1, Jend+1
! The problem is that the fluxes on the inter-processor boundaries may be
!modified if the point next to the boundary (e.g., in the ghost cells,
!managed by a different processor) is land, but as is the routine doesnt
!check for this. 

     DO j=Jbeg-1,Jend+1
     DO i=Ibeg-1,Iend+1
      IF(MASK(I,J)<1)THEN
        P(I,J)=ZERO
! Jeff reported a bug here for parallel version
# if defined(PARALLEL)
        IF((I/=Ibeg).or.(n_west.ne.MPI_PROC_NULL))THEN
# else
        IF(I/=Ibeg)THEN
# endif
!         Fx(I,J)=0.5_SP*GRAV*HxL(I,J)*HxL(I,J)*MASK(I-1,J)
!new splitting method
      Xi=EtaRxL(I,J)
      Deps=Depthx(I,J)
         Fx(I,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I-1,J)
        ELSE
         Fx(I,J)=ZERO
        ENDIF
        Gx(I,J)=ZERO

        P(I+1,J)=ZERO
! Jeff also here
# if defined(PARALLEL)
        IF((I/=Iend).or.(n_east.ne.MPI_PROC_NULL))THEN
# else
        IF(I/=Iend)THEN
# endif
!         Fx(I+1,J)=0.5_SP*GRAV*HxR(I+1,J)*HxR(I+1,J)*MASK(I+1,J)
! new splitting method
      Xi=EtaRxR(I+1,J)
      Deps=Depthx(I+1,J)
         Fx(I+1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I+1,J)
        ELSE
         Fx(I+1,J)=ZERO
        ENDIF
        Gx(I+1,J)=ZERO

        Q(I,J)=ZERO
        Fy(I,J)=ZERO
! Jeff also here
# if defined(PARALLEL)
        IF((J/=Jbeg).or.(n_suth.ne.MPI_PROC_NULL))THEN
# else
        IF(J/=Jbeg)THEN
# endif
!         Gy(I,J)=0.5_SP*GRAV*HyL(I,J)*HyL(I,J)*MASK(I,J-1)
! new splitting method
      Xi=EtaRyL(I,J)
      Deps=Depthy(I,J)
         Gy(I,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I,J-1)
        ELSE
         Gy(I,J)=ZERO
        ENDIF

        Q(I,J+1)=ZERO
        Fy(I,J+1)=ZERO
! Jeff also here
# if defined(PARALLEL)
        IF((J/=Jend).or.(n_nrth.ne.MPI_PROC_NULL))THEN
# else
        IF(J/=Jend)THEN
# endif
!         Gy(I,J+1)=0.5_SP*GRAV*HyR(I,J+1)*HyR(I,J+1)*MASK(I,J+1)
! new splitting method
      Xi=EtaRyR(I,J+1)
      Deps=Depthy(I,J+1)
         Gy(I,J+1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I,J+1)
        ELSE
         Gy(I,J+1)=ZERO
        ENDIF
      ENDIF
     ENDDO
     ENDDO

END SUBROUTINE BOUNDARY_CONDITION


!-------------------------------------------------------------------
!
!   This subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi, use dummy variables 2) add vtype=3
!
!-------------------------------------------------------------------
SUBROUTINE EXCHANGE_DISPERSION
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: VTYPE

# if defined (CARTESIAN)
    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxx,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxx,VTYPE,PERIODIC)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vyy,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVyy,VTYPE,PERIODIC)

    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxy,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxy,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vxy,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVxy,VTYPE,PERIODIC)

    IF(Gamma2>ZERO)THEN

      IF(DISP_TIME_LEFT)THEN
        VTYPE=1 ! symetric in both direction
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAT,VTYPE,PERIODIC)
        VTYPE=2  ! like u
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETATx,VTYPE,PERIODIC)
        VTYPE=3  ! like v
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETATy,VTYPE,PERIODIC) 
      ELSE
        VTYPE=2
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ut,VTYPE,PERIODIC)
        VTYPE=3
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vt,VTYPE,PERIODIC)

        VTYPE=1
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utx,VTYPE,PERIODIC)
        VTYPE=1
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vty,VTYPE,PERIODIC)

        VTYPE=2
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utxx,VTYPE,PERIODIC)
        VTYPE=3
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vtyy,VTYPE,PERIODIC)

        VTYPE=1
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utxy,VTYPE,PERIODIC) 
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vtxy,VTYPE,PERIODIC) 

        VTYPE=2
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUtxx,VTYPE,PERIODIC)
        VTYPE=3
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVtyy,VTYPE,PERIODIC)

        VTYPE=1
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUtxy,VTYPE,PERIODIC) 
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVtxy,VTYPE,PERIODIC) 
    
      ENDIF

      VTYPE=1  ! symetric in both direction
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ux,VTYPE,PERIODIC)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUx,VTYPE,PERIODIC)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vy,VTYPE,PERIODIC)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVy,VTYPE,PERIODIC)
      VTYPE=3  !like v
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uy,VTYPE,PERIODIC)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUy,VTYPE,PERIODIC)
      Vtype=2  !like u
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vx,VTYPE,PERIODIC)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVx,VTYPE,PERIODIC)
      VTYPE=2  ! like u
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAx,VTYPE,PERIODIC)
      VTYPE=3  ! like v
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAy,VTYPE,PERIODIC)   

    ENDIF
# else
    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxx,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxx,VTYPE)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vyy,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVyy,VTYPE)

    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxy,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxy,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vxy,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVxy,VTYPE)

      VTYPE=1  ! symetric in both direction
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ux,VTYPE)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUx,VTYPE)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vy,VTYPE)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVy,VTYPE)
      VTYPE=3  !like v
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uy,VTYPE)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUy,VTYPE)
      Vtype=2  !like u
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vx,VTYPE)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVx,VTYPE)
      VTYPE=2  ! like u
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAx,VTYPE)
      VTYPE=3  ! like v
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAy,VTYPE)   
# endif
   
    
END SUBROUTINE EXCHANGE_DISPERSION

!---------------------------------------------------------------------------------------
!
!   EXCHANGE subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi 
!     1) use dummy variables 2) add vtype=3
!   08/19/2015 Choi, corrected segmentation fault, maybe memory leaking                                       
!
!---------------------------------------------------------------------------------------
SUBROUTINE EXCHANGE
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: VTYPE
# if defined (CARTESIAN)
    REAL(SP),DIMENSION(Mloc,Nloc) :: rMASK
#   if defined (FILTERING)
      INTEGER :: I,J
#   endif
# endif

# if defined (CARTESIAN)
    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta,VTYPE,PERIODIC)

    IF(VISCOSITY_BREAKING) THEN
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,AGE_BREAKING,VTYPE,PERIODIC)
    ENDIF
    
    !ykchoi (08.19.2015) :: new variable (WAVEMAKER_VIS)
    IF(VISCOSITY_BREAKING .OR. WAVEMAKER_VIS) THEN
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,nu_break,VTYPE,PERIODIC)
    ENDIF

    IF(VISCOSITY_BREAKING .OR. DIFFUSION_SPONGE .OR. WAVEMAKER_VIS)THEN
      CALL PHI_COLL_VARIABLE_LENGTH(Mloc1,Nloc,Ibeg,Iend1,Jbeg,Jend,Nghost,P,VTYPE)
      CALL PHI_COLL_VARIABLE_LENGTH(Mloc,Nloc1,Ibeg,Iend,Jbeg,Jend1,Nghost,Q,VTYPE)
    ENDIF

    rMASK = MASK ! for periodic boundary condition
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,rMASK,VTYPE,PERIODIC)  
    MASK = INT(rMASK)

    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU,VTYPE,PERIODIC)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV,VTYPE,PERIODIC)
# else
    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta,VTYPE)

    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU,VTYPE)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V,VTYPE)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV,VTYPE)
# endif

! etaR x mask is a wrong idea
!    Eta=Eta*MASK

    U=U*MASK
    V=V*MASK
    HU=HU*MASK
    HV=HV*MASK

# if defined (FILTERING)
     filter_count = filter_count + 1 
# if defined (CARTESIAN)
     IF (filter_count > filter_interval) THEN
       filter_count =0

       CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK,Eta)

       CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK,U)
       CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK,V)

    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta,VTYPE,PERIODIC)
    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U,VTYPE,PERIODIC)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V,VTYPE,PERIODIC)

    H=Eta*Gamma3+Depth

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend   
       IF(MASK(I,J)==1)THEN
        HU(I,J)=Max(H(I,J),MinDepthFrc)*U(I,J)
        HV(I,J)=Max(H(I,J),MinDepthFrc)*V(I,J)
       ENDIF
     ENDDO
     ENDDO

     ENDIF ! end filter_count
# endif
# endif
    
END SUBROUTINE EXCHANGE

!-----------------------------------------------------------------------------------
!
!   PHI_COLL_VARIABLE_LENGTH subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi                                      
!
!-----------------------------------------------------------------------------------
SUBROUTINE PHI_COLL_VARIABLE_LENGTH(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE)

    USE PARAM
    USE GLOBAL,ONLY : WaveMaker

# if defined (PARALLEL)
    USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth
# endif

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)
    INTEGER :: I,J,K
    nop(VTYPE)

      ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif

    IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF

# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif

      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO

# if defined (PARALLEL)
    endif
# endif


      ! y-direction and corners
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif

      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO

# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif

      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO

# if defined (PARALLEL)
    endif
# endif


END SUBROUTINE PHI_COLL_VARIABLE_LENGTH

!-------------------------------------------------------------------------------------
!
!   PHI_COLL subroutine is used to collect data into ghost cells
!
!   HISTORY:
!     05/01/2010  Fengyan Shi
!     09/07/2010 Fengyan Shi, fix:
!       1) u v symmetric problem, 2) remove use global 3) fix bug
!     05/27/2010 Gangfeng Ma, corrected some bugs
!
!-------------------------------------------------------------------------------------
# if defined (CARTESIAN)
SUBROUTINE PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,PERIODIC)
# else
SUBROUTINE PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE)
# endif
    USE PARAM
    USE GLOBAL,ONLY : WaveMaker
# if defined (COUPLING)
    USE GLOBAL,ONLY : IN_DOMAIN_WEST,IN_DOMAIN_EAST,IN_DOMAIN_SOUTH,IN_DOMAIN_NORTH
# endif
# if defined (PARALLEL)
# if defined (COUPLING)
    USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth,&
        IN_DOMAIN_WEST,IN_DOMAIN_EAST,IN_DOMAIN_SOUTH,IN_DOMAIN_NORTH
# else
    USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth
# endif
# if defined (CARTESIAN)
    USE GLOBAL, ONLY : ier,myid,PX,PY,ProcessorID
# endif
# endif
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)
# if defined (CARTESIAN)
    LOGICAL :: PERIODIC
# endif
# if defined (PARALLEL)
# if defined (CARTESIAN)
    INTEGER,DIMENSION(1) :: req
    INTEGER :: len,II
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    REAL(SP),DIMENSION(Mloc,Nghost) :: xx,send2d
# endif
# endif
    INTEGER I,J,K

! periodic first because it is not related to VTYPE
# if defined (CARTESIAN)
      IF(PERIODIC)THEN
# if defined (PARALLEL)
!   _____________________________ exchange
    IF (PY>1)THEN

      len=Mloc*Nghost
! south
    DO II = 1,PX

    if(myid==ProcessorID(II,1)) then
      DO I = 1,Mloc
      DO J = 1,Nghost
        send2d(I,J)=PHI(I,J+Nghost)
      ENDDO
      ENDDO

! send from master
        call MPI_SEND(send2d,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,ier)

    endif ! end myid

        if(myid==ProcessorID(II,PY))then
        call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
          DO I=1,Mloc
          DO J=1,Nghost
            PHI(I,Jend+J)=xx(I,J)
          ENDDO
          ENDDO
        endif

! north
    if(myid==ProcessorID(II,PY)) then
      DO I = 1,Mloc
      DO J = 1,Nghost
        send2d(I,J)=PHI(I,Jend-Nghost+J)
      ENDDO
      ENDDO
! send from master
        call MPI_SEND(send2d,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
    endif ! end myid

        if(myid==ProcessorID(II,1))then
        call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
          DO I=1,Mloc
          DO J=1,Nghost
            PHI(I,J)=xx(I,J)
          ENDDO
          ENDDO
        endif
 
    ENDDO  ! end PX

    ELSE  ! PY = 1
      DO I = Ibeg,Iend
      DO J = 1,Nghost
       PHI(I,J) = PHI(I,Jend-Nghost+J)
       PHI(I,Jend+J) = PHI(I,J+Nghost)
      ENDDO
      ENDDO
    ENDIF

!   ----------------------------- end exchange
# else 
!    not parallel
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=PHI(I,Jend-Nghost+K)
          PHI(I,Jend+K)=PHI(I,Nghost+K)
        ENDDO
        ENDDO        
# endif 
! end parallel

    ENDIF ! end periodic 

# endif 
! end cartesian

! I added coupling condition 10/14/2012
!   add left_bc wavemaker 09/13/2017

! for Eta
    IF(VTYPE==1) THEN  ! for eta
      ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF  ! end left bc wavemaker

# if defined (COUPLING)
    ENDIF
# endif

# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif


      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=1

! for u
    IF(VTYPE==2) THEN  ! for u (x-mirror condition)
      ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF ! end left_bc wavemaker

# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif


      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 

# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=2

    IF(VTYPE==3) THEN ! for v (y-mirror condition)
! for v
      ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF ! end left_bc wavemaker

# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif


      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 
  ! end cartesian
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif
# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=3

! for cross-derivatives
    IF(VTYPE==4) THEN ! VTYPE==4 for u and v cross-mirror
     ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=ZERO
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=ZERO
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

      ! y-direction and corners, this one is not an exact solution
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=ZERO
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=ZERO
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif


     ENDIF ! end vtype=4

! for symmetric
    IF(VTYPE==5)THEN
      ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
       ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

      ! y-direction and corners

# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

     ENDIF ! end vtype=5

! for anti-symmetric
     IF(VTYPE==6)THEN
      ! x-direction
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO 
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO 
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif

      ! y-direction and corners
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO   
# if defined (COUPLING)
    ENDIF
# endif  
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO     
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif
# endif


    ENDIF ! end vtype=6

# if defined (PARALLEL)
    call phi_exch (PHI)
# endif

END SUBROUTINE PHI_COLL



