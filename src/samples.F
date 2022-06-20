!------------------------------------------------------------------------------------
!
!      FILE sample.F
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
!    SOLITARY_WAVE_LEFT_BOUNDARY is subroutine to provide 
!    solitary wave at left boundary
!    it can be specified in input.txt giving 'SOL'
!
!    HISTORY: 
!      05/28/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
#define nop(x) associate( x => x ); end associate

SUBROUTINE SOLITARY_WAVE_LEFT_BOUNDARY
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP):: aa,h00,c1,tex,tlag,zz
     INTEGER::Iwavemaker

# if defined (PARALLEL)
       if ( n_west .eq. MPI_PROC_NULL ) then
# endif
       Iwavemaker=Ibeg
       aa=AMP_SOLI
       h00=DEP_SOLI
       c1=sqrt(GRAV*h00*(1.0_SP+aa/h00))
       DO J=1,Nloc
         tex=sqrt(0.75_SP*aa/h00**3)
         tlag=4.0_SP*h00/sqrt(aa/h00)
         zz=aa/COSH(tex*(Lag_soli-c1*TIME))**2
         Eta(Iwavemaker,J)=zz
         H(Iwavemaker,J)=Eta(Iwavemaker,J)+Depth(Iwavemaker,J) 
! note: can not provide u and hu at boundary for dispersive equations!
!         U(Iwavemaker,J)= SQRT(grav/h00)*zz
!         HU(Iwavemaker,J)=h00*U(Iwavemaker,J)       
       enddo   
# if defined (PARALLEL)
       endif
# endif
     
END SUBROUTINE SOLITARY_WAVE_LEFT_BOUNDARY


!-------------------------------------------------------------------------------------
!
!    INITIAL_WAVE is subroutine to provide initial wave condition
!    it can be specified in input.txt giving 'INI'
!
!    HISTORY: 
!       10/01/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIAL_WAVE
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP),Dimension(Mloc,Nloc)::XX,YY
     REAL(SP):: sigma,x_c,y_c,a
     INTEGER :: ii1,ii2

! the initial wave domain includes ghost cells

     do j=1,Nloc
     do i=1,Mloc
# if defined (CARTESIAN)
       xx(i,j)=(i-1.)*dx
       yy(i,j)=(j-1.)*dy
# else
       xx(i,j)=(i-1.)*dx(I,J)
       yy(i,j)=(j-1.)*dy(I,J)
# endif
     enddo
     enddo

     sigma=0.5
     x_c=10.0
     y_c=10.0
     a=0.1

     do j=1,Nloc
     do i=1,Mloc
     goto 200
! box
     tmp1=21+Nghost
     tmp2=31+Nghost
     if(i>tmp1.and.i<tmp2.and.j>tmp1.and.j<tmp2)then
      Eta(i,j)=1.0_SP
     else
      Eta(i,j)=zero
     endif

200  continue

     goto 100
! dam break
        if(i<100+Nghost)then
         Eta(i,j)=5.0_SP
        else
         Eta(i,j)=0.0_SP
        endif
100   continue
! gausian
!         r=sqrt((xx(i,j)-x_c)**2+(yy(i,j)-y_c)**2)
!         Eta(i,j)=a*exp(-r**2/2./sigma**2)
     enddo
     enddo

! alongshore crest
     goto 213
     ii1=21+Nghost
     ii2=25+Nghost
     DO J=1,Nloc
       DO I=ii1,ii2
         Eta(I,J)=1.0_SP
       ENDDO
     ENDDO
213  continue

END SUBROUTINE INITIAL_WAVE


!-------------------------------------------------------------------------------------
!
!    INITIAL_RECTANGULAR is subroutine to provide initial rectangular hump 
!    it can be specified in input.txt giving 'INI_REC'
!
!    HISTORY: 
!      10/11/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
# if defined (CARTESIAN)
SUBROUTINE INITIAL_RECTANGULAR(M,N,Nghost,DX,DY,Xc,Yc,WID,AMP,Eta)
# else
SUBROUTINE INITIAL_RECTANGULAR(M,N,Nghost,Dphi,Dtheta,Lon_West,Lat_South, &
             Xc,Yc,WID,AMP,Eta)
# endif
     USE PARAM
# if defined (PARALLEL)
     USE GLOBAL, ONLY : iista, jjsta   !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX,DY,Xc,Yc,WID,AMP
# else
     REAL(SP),INTENT(IN) :: Xc,Yc,WID,AMP,Lon_West,Lat_South,Dphi,Dtheta
# endif
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Eta
     INTEGER :: Il,Jl,Ir,Jr

     Eta = ZERO
# if defined (PARALLEL)

# if defined (CARTESIAN)
![---ykchoi Jan/23/2018
!     Il=Nghost+Xc/DX-WID/DX +1 - npx*Mglob/px
!     Ir=Nghost+Xc/DX+WID/DX +1 - npx*Mglob/px
!     Jl=Nghost+Yc/DY-WID/DY +1 - npy*Nglob/py
!     Jr=Nghost+Yc/DY+WID/DY +1 - npy*Nglob/py
     Il=INT(Nghost+Xc/DX-WID/DX) +1 - (iista - 1)
     Ir=INT(Nghost+Xc/DX+WID/DX) +1 - (iista - 1)
     Jl=INT(Nghost+Yc/DY-WID/DY) +1 - (jjsta - 1)
     Jr=INT(Nghost+Yc/DY+WID/DY) +1 - (jjsta - 1)
!---ykchoi Jan/23/2018]
# else
![---ykchoi Jan/23/2018
!     Il=Nghost+(Xc-Lon_West)/dphi-WID/dphi +1 - npx*Mglob/px
!     Ir=Nghost+(Xc-Lon_West)/dphi+WID/dphi +1 - npx*Mglob/px
!     Jl=Nghost+(Yc-Lat_South)/dtheta-WID/dtheta +1 - npy*Nglob/py
!     Jr=Nghost+(Yc-Lat_South)/dtheta+WID/dtheta +1 - npy*Nglob/py
     Il=INT(Nghost+(Xc-Lon_West)/dphi-WID/dphi) +1 - (iista - 1)
     Ir=INT(Nghost+(Xc-Lon_West)/dphi+WID/dphi) +1 - (iista - 1)
     Jl=INT(Nghost+(Yc-Lat_South)/dtheta-WID/dtheta) +1 - (jjsta - 1)
     Jr=INT(Nghost+(Yc-Lat_South)/dtheta+WID/dtheta) +1 - (jjsta - 1)
!---ykchoi Jan/23/2018]
# endif

     IF(Il>M.or.Ir<1.or.Jl>N.or.Jr<1)THEN
     ELSE
      IF(Il<1)Il=1
      IF(Ir>M)Ir=M
      IF(Jl<1)Jl=1
      IF(Jr>N)Jr=N
      Eta(Il:Ir,Jl:Jr) = AMP
     ENDIF
# else

# if defined (CARTESIAN)
     Il=INT(Xc/DX+1 - WID/DX) + Nghost
     Ir=INT(Xc/DX+1 + WID/DX) + Nghost
     Jl=INT(Yc/DY+1 - WID/DY) + Nghost
     Jr=INT(Yc/DY+1 + WID/DY) + Nghost
# else
     Il=INT(Nghost+(Xc-Lon_West)/dphi-WID/dphi) +1
     Ir=INT(Nghost+(Xc-Lon_West)/dphi+WID/dphi) +1
     Jl=INT(Nghost+(Yc-Lat_South)/dtheta-WID/dtheta) +1
     Jr=INT(Nghost+(Yc-Lat_South)/dtheta+WID/dtheta) +1
# endif

     Eta(Il:Ir,Jl:Jr) = AMP
# endif
 
END SUBROUTINE INITIAL_RECTANGULAR


!-------------------------------------------------------------------------------------
!
!    INITIAL_GAUSIAN is subroutine to provide initial gausian hump 
!    it can be specified in input.txt giving 'INI_REC'
!
!    HISTORY:
!       10/11/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
# if defined (CARTESIAN)
SUBROUTINE INITIAL_GAUSIAN(M,N,Nghost,DX,DY,Xc,Yc,AMP,WID,Eta)
# else
SUBROUTINE INITIAL_GAUSIAN(M,N,Nghost,Dphi,Dtheta,Lon_West,Lat_South, &
             Xc,Yc,AMP,WID,Eta)
# endif
     USE PARAM
# if defined (PARALLEL)
     USE GLOBAL, ONLY : iista, jjsta   !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX,DY,Xc,Yc,WID,AMP
# else
     REAL(SP),INTENT(IN) :: Xc,Yc,WID,AMP,Lon_West,Lat_South,Dphi,Dtheta
# endif
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Eta
     REAL(SP) :: r,xr,yr
     
     Eta = ZERO
# if defined (PARALLEL)

# if defined (CARTESIAN)
     DO J=1,N
     DO I=1,M
![---ykchoi Jan/23/2018
!       xr=(I-1)*DX-Nghost*DX + (npx*Mglob/px)*DX
!       yr=(J-1)*DY-Nghost*DY + (npy*Mglob/py)*DY
       xr=(I-1)*DX-Nghost*DX + (iista - 1)*DX
       yr=(J-1)*DY-Nghost*DY + (jjsta - 1)*DY
!---ykchoi Jan/23/2018]
       r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       Eta(I,J) = AMP*EXP(-r**2/1.0_SP/WID**2)
     ENDDO
     ENDDO
# else
     DO J=1,N
     DO I=1,M
![---ykchoi Jan/23/2018
!      xr=(I-1)*dphi-Nghost*dphi + (npx*Mglob/px)*dphi+Lon_West
!      yr=(J-1)*dtheta-Nghost*dtheta + (npy*Nglob/py)*dtheta+Lat_South
      xr=(I-1)*dphi-Nghost*dphi + (iista - 1)*dphi+Lon_West
      yr=(J-1)*dtheta-Nghost*dtheta + (jjsta - 1)*dtheta+Lat_South
!---ykchoi Jan/23/2018]
      r=SQRT((xr-Xc)**2+(yr-Yc)**2)
      Eta(I,J) = AMP*EXP(-r**2/1.0_SP/WID**2)
     ENDDO
     ENDDO
# endif

# else

# if defined (CARTESIAN)
     DO J=1,N
     DO I=1,M
       xr=(I-1)*DX-Nghost*DX 
       yr=(J-1)*DY-Nghost*DY
       r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       Eta(I,J) = AMP*EXP(-r**2/1.0_SP/WID**2)
     ENDDO
     ENDDO
# else
     DO J=1,N
     DO I=1,M
      xr=(I-1)*dphi-Nghost*dphi +Lon_West
      yr=(J-1)*dtheta-Nghost*dtheta +Lat_South
      r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       Eta(I,J) = AMP*EXP(-r**2/1.0_SP/WID**2)
     ENDDO
     ENDDO
# endif

# endif
 
END SUBROUTINE INITIAL_GAUSIAN


!-------------------------------------------------------------------------------------
!
!    INITIAL_DIPOLE is subroutine to provide initial dipole 
!    it can be specified in input.txt giving 'INI_DIP'
!
!    HISTORY: 
!        08/24/2012 Fengyan Shi
!
!-------------------------------------------------------------------------------------
# if defined (CARTESIAN)
SUBROUTINE INITIAL_DIPOLE(M,N,Nghost,DX,DY,Xc,Yc,AMP,WID,Eta)
# else
SUBROUTINE INITIAL_DIPOLE(M,N,Nghost,Dphi,Dtheta,Lon_West,Lat_South, &
             Xc,Yc,AMP,WID,Eta)
# endif
     USE PARAM
# if defined (PARALLEL)
     USE GLOBAL, ONLY : iista, jjsta   !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Nghost
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX,DY,Xc,Yc,WID,AMP
# else
     REAL(SP),INTENT(IN) :: Xc,Yc,WID,AMP,Lon_West,Lat_South,Dphi,Dtheta
# endif
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Eta
     REAL(SP) :: r,xr,yr,uni
     
     Eta = ZERO
# if defined (PARALLEL)

# if defined (CARTESIAN)
     DO J=1,N
     DO I=1,M
![---ykchoi Jan/23/2018
!       xr=(I-1)*DX-Nghost*DX + (npx*Mglob/px)*DX
!       yr=(J-1)*DY-Nghost*DY + (npy*Nglob/py)*DY
       xr=(I-1)*DX-Nghost*DX + (iista - 1)*DX
       yr=(J-1)*DY-Nghost*DY + (jjsta - 1)*DY
!---ykchoi Jan/23/2018]
       r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       uni=-1.414_SP/WID*(xr-Xc)
       Eta(I,J) = AMP*uni*EXP(-r**2/1.0_SP/WID**2+0.5_SP)
     ENDDO
     ENDDO
# else
     DO J=1,N
     DO I=1,M
![---ykchoi Jan/23/2018
!      xr=(I-1)*dphi-Nghost*dphi + (npx*Mglob/px)*dphi+Lon_West
!      yr=(J-1)*dtheta-Nghost*dtheta + (npy*Nglob/py)*dtheta+Lat_South
      xr=(I-1)*dphi-Nghost*dphi + (iista - 1)*dphi+Lon_West
      yr=(J-1)*dtheta-Nghost*dtheta + (jjsta - 1)*dtheta+Lat_South
!---ykchoi Jan/23/2018]
      r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       uni=-1.414_SP/WID*(xr-Xc)
       Eta(I,J) = AMP*uni*EXP(-r**2/1.0_SP/WID**2+0.5_SP)
     ENDDO
     ENDDO
# endif

# else

# if defined (CARTESIAN)
     DO J=1,N
     DO I=1,M
       xr=(I-1)*DX-Nghost*DX 
       yr=(J-1)*DY-Nghost*DY
       r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       uni=-1.414_SP/WID*(xr-Xc)
       Eta(I,J) = AMP*uni*EXP(-r**2/1.0_SP/WID**2+0.5_SP)
     ENDDO
     ENDDO
# else
     DO J=1,N
     DO I=1,M
      xr=(I-1)*dphi-Nghost*dphi+Lon_West
      yr=(J-1)*dtheta-Nghost*dtheta+Lat_South
      r=SQRT((xr-Xc)**2+(yr-Yc)**2)
       uni=-1.414_SP/WID*(xr-Xc)
       Eta(I,J) = AMP*uni*EXP(-r**2/1.0_SP/WID**2+0.5_SP)
     ENDDO
     ENDDO
# endif

# endif
 
END SUBROUTINE INITIAL_DIPOLE


!-------------------------------------------------------------------------------------
!
!    INITIAL_N_WAVE is subroutine to provide initial N wave solution 
!    it can be specified in input.txt giving N_WAVE
!
!    HISTORY: 
!        11/30/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIAL_N_WAVE(M,N,DX,x1,x2,Ha,gamma_s,dep,U,V,Eta)
# if defined (PARALLEL)
     USE GLOBAL, only: iista   !ykchoi Jan/23/2018
# endif
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N
     REAL(SP),INTENT(IN) :: x1,x2,Ha,gamma_s,dep
     REAL(SP) :: a0,gamma
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX
# endif
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: U,V,Eta
     nop(gamma_s)
     nop(x2)
     IF(dep==ZERO)THEN
      WRITE(*,*)'depth for Nwave should not be 0, stop!'
      STOP
     ENDIF
! recalculate a0 and gamma based on Synolakis etal (A15)
! note non-dimensional formula in Tadepalli and Synolakis we use dimensional one
! like in xi zhao etal icce 2010
     a0=3.0_SP*SQRT(3.0_SP)*Ha/2.0_SP
     gamma=3.0_SP*0.5_SP/dep*SQRT(SQRT(3.0_SP/4.0_SP)*Ha/dep)
!print*,gamma

     DO J=1,N
     DO I=1,M
# if defined (CARTESIAN)
# if defined (PARALLEL)
! start from ghost cell
![---ykchoi Jan/23/2018       
!     tmp1=(npx*Mglob/px+I-1)*DX
      tmp1=( (iista - 1) +I-1)*DX
!---ykchoi Jan/23/2018]
# else
       tmp1=(I-1)*DX
# endif
# else 
# if defined (PARALLEL)
! start from ghost cell
![---ykchoi Jan/23/2018       
!       tmp1=(npx*Mglob/px+I-1)*DX(I,J)
      tmp1=( (iista - 1) +I-1)*DX(I,J)
!---ykchoi Jan/23/2018]
# else
       tmp1=(I-1)*DX(I,J)
# endif
# endif
! general n-wave
!       Eta(I,J)=a0*(tmp1-x2)/COSH(gamma*(tmp1-x1))
! isosceles n-wave
       Eta(I,J)=a0*COSH(gamma*(tmp1-x1))**(-2)*tanh(gamma*(tmp1-x1))
       U(I,J) = SQRT(GRAV/dep)*Eta(I,J)
       V(I,J) = ZERO
     ENDDO
     ENDDO


END SUBROUTINE INITIAL_N_WAVE


!-------------------------------------------------------------------------------------
!
!    INITIAL_SOLITARY_WAVE is subroutine to provide initial solitary wave solution 
!    it can be specified in input.txt giving SOL
!
!    HISTORY: 
!        10/01/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE INITIAL_SOLITARY_WAVE(M,N,DX,Xwavemaker,&
           AMP_Soli,Dep_Soli,Beta,U,V,Eta,SolitaryPositiveDirection)
     USE PARAM
# if defined (PARALLEL)
     USE GLOBAL, only: iista    !ykchoi Jan/23/2018
# endif
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX 
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX 
# endif
     REAL(SP),INTENT(IN) :: Dep_Soli,Xwavemaker,Beta
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: U,V,Eta
     REAL(SP),INTENT(IN) :: AMP_Soli
     REAL(SP) :: A,B,A1,A2,alpha,SC,Cph,C_ph
     LOGICAL,INTENT(IN) :: SolitaryPositiveDirection
     nop(Beta)

!     alpha=0.5_SP*Beta*Beta+Beta
   ! something wrong here. if beta=-0.531, then alpha should be -0.39
   ! the solitary wave shape does not match the equation. but alpha = -0.10 seems OK
     alpha=-0.39

      CALL SUB_SLTRY(AMP_Soli,Dep_Soli,alpha,Cph,B,A1,A2,A,C_ph)

     DO J=1,N
     DO I=1,M
# if defined (CARTESIAN)
# if defined (PARALLEL)
! start from ghost cell
![---ykchoi Jan/23/2018
!       SC=1.0_SP/COSH(B*(npx*Mglob/px+I-Xwavemaker/DX-1)*DX)
       SC=1.0_SP/COSH(B*( (iista - 1) +I-Xwavemaker/DX-1)*DX)
!---ykchoi Jan/23/2018]
# else
       SC=1.0_SP/COSH(B*(I-Xwavemaker/DX-1)*DX)
# endif
# else
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!       SC=1.0_SP/COSH(B*(npx*Mglob/px+I-Xwavemaker/DX(I,J)-1)*DX(I,J))
       SC=1.0_SP/COSH(B*( (iista - 1) +I-Xwavemaker/DX(I,J)-1)*DX(I,J))
!---ykchoi Jan/23/2018]
# else
       SC=1.0_SP/COSH(B*(I-Xwavemaker/DX(I,J)-1)*DX(I,J))
# endif
# endif
       Eta(I,J)=A1*SC*SC+A2*SC*SC*SC*SC
       IF(SolitaryPositiveDirection)THEN
         U(I,J)=A*SC*SC
       ELSE
         U(I,J)=-A*SC*SC
       ENDIF
       V(I,J)=ZERO
     ENDDO
     ENDDO


END SUBROUTINE INITIAL_SOLITARY_WAVE


!-------------------------------------------------------------------------------------
!
!    SUB_SLTRY is subroutine to provide solitary wave solution of Nogu equations
!
!    HISTORY: 
!        10/05/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE SUB_SLTRY(a0, h1, alp, r1, bue, ae1, ae2, au, C_ph)
    USE PARAM
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a0,h1,alp
    REAL(SP), INTENT(OUT) :: r1,bue,ae1,ae2,au,C_ph
    REAL(SP) :: alp2,p,q,r,x,fx, fpx,rx,cph,eps,zr,g
    INTEGER :: ite

!
!--------coefficients for third order polynomial equations
!               " x**3+p*x**2+q*x+r=0  "

         alp2 = alp + 1.0_SP/3.0_SP
         eps  = a0/h1
         g=GRAV

         p = -(alp2+2.0_SP*alp*(1.0_SP+eps))/(2.0_SP*alp)
         q = eps*alp2/alp
         r = alp2/(2.0_SP*alp)

!--------Newton-Rapsons method to solve x ( >1 )for the above equation
         ite = 0
         x=1.2_SP
1        fx  = r+x*(q+x*(p+x))
         fpx = q+x*(2.0_SP*p+3.0_SP*x)
         x = x-fx/fpx

         ite = ite+1
         if (ite.gt.10) then
            !write(*,*) 'no solitary wave solution (check eps = a0/h1)'
            stop
         endif
         if (abs(fx).ge.1e-5) goto 1
         rx = sqrt(x)
         cph = sqrt(g*h1)
         r1 = rx*cph
         C_ph=rx*cph

         !write (*,*) rx, r1


!---------coefficients for solitary solutions :
!         "   u =  au/(cosh(bue*xi))**2    "
!         "  et = ae1/(cosh(bue*xi)**2+ae2/(cosh(bue*xi)**4   "

          au =  (x-1.0_SP)/(eps*rx)*cph*eps
         bue =  sqrt((x-1.0_SP)/(4.0_SP*(alp2-alp*x)))/h1
         ae1 =  (x-1.0_SP)/(eps*3.0_SP*(alp2-alp*x))*a0
         ae2 = -(x-1.0_SP)/(2.0_SP*eps)*(x-1.0_SP)*(2.0_SP*alp*x+alp2) &
                                      /(x*(alp2-alp*x))*a0
         zr  = bue*rx*cph

END SUBROUTINE SUB_SLTRY
