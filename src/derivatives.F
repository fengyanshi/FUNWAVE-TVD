!------------------------------------------------------------------------------------
!
!      FILE derivatives.F
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
!    DERIVATIVE_Y_High is subroutine to calculation first-derivative y with higher-order
!
!    HISTORY:
!      10/11/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_Y_High(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DY,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DY
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J)= (Uin(I,J+2)+2.0_SP*Uin(I,J+1)    &
                 -Uin(I,J-2)-2.0_SP*Uin(I,J-1))/DY/8.0_SP*MASK(I,J)
# else
       Uout(I,J)= (Uin(I,J+2)+2.0_SP*Uin(I,J+1)    &
                 -Uin(I,J-2)-2.0_SP*Uin(I,J-1))/DY(I,J)/8.0_SP*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_Y_High

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_X_High is subroutine to calculation 
!    the first-derivative x with higher-order
!
!    HISTORY: 
!      10/11/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_X_High(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J)= (Uin(I+2,J)+2.0_SP*Uin(I+1,J)    &
                 -Uin(I-2,J)-2.0_SP*Uin(I-1,J))/DX/8.0_SP*MASK(I,J)
# else
       Uout(I,J)= (Uin(I+2,J)+2.0_SP*Uin(I+1,J)    &
                 -Uin(I-2,J)-2.0_SP*Uin(I-1,J))/DX(I,J)/8.0_SP*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_X_High

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_Y is subroutine to calculation first-derivative y
!
!    HISTORY: 
!      10/11/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_Y(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DY,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DY
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J)= (Uin(I,J+1)   &
                 -Uin(I,J-1))/DY/2.0_SP*MASK(I,J)
# else
       Uout(I,J)= (Uin(I,J+1)   &
                 -Uin(I,J-1))/DY(I,J)/2.0_SP*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_Y

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_X is subroutine to calculation first-derivative x
!
!    HISTORY: 
!      10/11/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_X(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J)= (Uin(I+1,J)    &
                 -Uin(I-1,J))/DX/2.0_SP*MASK(I,J)
# else
       Uout(I,J)= (Uin(I+1,J)    &
                 -Uin(I-1,J))/DX(I,J)/2.0_SP*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_X

!-------------------------------------------------------------------------------------
!    DERIVATIVE_YY is subroutine to calculation 2nd-derivative yy
!
!    HISTORY: 
!       09/21/2010 Fengyan Shi, University of Delaware
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_YY(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DY,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DY
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J)= (Uin(I,J+1)-2.0_SP*Uin(I,J) & 
                 +Uin(I,J-1))/DY/DY*MASK(I,J)
# else
       Uout(I,J)= (Uin(I,J+1)-2.0_SP*Uin(I,J) & 
                 +Uin(I,J-1))/DY(I,J)/DY(I,J)*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_YY

!-------------------------------------------------------------------------------------
!    DERIVATIVE_XY is subroutine to calculation 2nd-derivative xy
!
!    HISTORY: 
!       09/21/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_XY(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,DY,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX,DY
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J
     REAL(SP) :: tmp1,tmp2

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       tmp1=(Uin(I+1,J+1)-Uin(I+1,J-1))/2.0_SP/DY
       tmp2=(Uin(I-1,J+1)-Uin(I-1,J-1))/2.0_SP/DY
       Uout(I,J)= (tmp1-tmp2)/2.0_SP/DX*MASK(I,J)
# else
       tmp1=(Uin(I+1,J+1)-Uin(I+1,J-1))/2.0_SP/DY(I,J)
       tmp2=(Uin(I-1,J+1)-Uin(I-1,J-1))/2.0_SP/DY(I,J)
       Uout(I,J)= (tmp1-tmp2)/2.0_SP/DX(I,J)*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_XY

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_XX is subroutine to calculation 2nd-derivative xx
!
!    HISTORY: 
!       09/21/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_XX(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

! I assume no 2nd derivative 
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J)= (Uin(I+1,J)-2.0_SP*Uin(I,J) & 
                 +Uin(I-1,J))/DX/DX*MASK(I,J)
# else
       Uout(I,J)= (Uin(I+1,J)-2.0_SP*Uin(I,J) & 
                 +Uin(I-1,J))/DX(I,J)/DX(I,J)*MASK(I,J)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_XX

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_XY_HIGH is subroutine to calculation 4th-derivative xy
!
!    HISTORY:
!       09/21/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_XY_HIGH(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,DY,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX,DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX,DY
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J
     REAL(SP) :: tmp1,tmp2,tmp3,tmp4

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       tmp1=1.0_SP/12.0_SP/DY*(Uin(I-2,J-2)-8.0_SP*Uin(I-2,J-1) &
                             +8.0_SP*Uin(I-2,J+1)-Uin(I-2,J+2))
       tmp2=1.0_SP/12.0_SP/DY*(Uin(I-1,J-2)-8.0_SP*Uin(I-1,J-1) &
                             +8.0_SP*Uin(I-1,J+1)-Uin(I-1,J+2))
       tmp3=1.0_SP/12.0_SP/DY*(Uin(I+1,J-2)-8.0_SP*Uin(I+1,J-1) &
                             +8.0_SP*Uin(I+1,J+1)-Uin(I+1,J+2))
       tmp4=1.0_SP/12.0_SP/DY*(Uin(I+2,J-2)-8.0_SP*Uin(I+2,J-1) &
                             +8.0_SP*Uin(I+2,J+1)-Uin(I+2,J+2))
       Uout(I,J)=MASK(I,J)/12.0_SP/DX*(tmp1-8.0_SP*tmp2 &
                             +8.0_SP*tmp3-tmp4)
# else
! to assure symmetric, use locally constant dx and dy
       tmp1=1.0_SP/12.0_SP/DY(I,J)*(Uin(I-2,J-2)-8.0_SP*Uin(I-2,J-1) &
                             +8.0_SP*Uin(I-2,J+1)-Uin(I-2,J+2))
       tmp2=1.0_SP/12.0_SP/DY(I,J)*(Uin(I-1,J-2)-8.0_SP*Uin(I-1,J-1) &
                             +8.0_SP*Uin(I-1,J+1)-Uin(I-1,J+2))
       tmp3=1.0_SP/12.0_SP/DY(I,J)*(Uin(I+1,J-2)-8.0_SP*Uin(I+1,J-1) &
                             +8.0_SP*Uin(I+1,J+1)-Uin(I+1,J+2))
       tmp4=1.0_SP/12.0_SP/DY(I,J)*(Uin(I+2,J-2)-8.0_SP*Uin(I+2,J-1) &
                             +8.0_SP*Uin(I+2,J+1)-Uin(I+2,J+2))
       Uout(I,J)=MASK(I,J)/12.0_SP/DX(I,J)*(tmp1-8.0_SP*tmp2 &
                             +8.0_SP*tmp3-tmp4)
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_XY_HIGH

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_XX_HIGH is subroutine to calculation 4th-derivative xx
!
!    HISTORY: 
!       05/30/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DERIVATIVE_XX_HIGH(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DX
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DX
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

! I assume no 2nd derivative 
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J) = MASK(I,J)*1.0_SP/12.0_SP/DX/DX*(-Uin(I-2,J)+16.0_SP*Uin(I-1,J)   &
                -30.0_SP*Uin(I,J)+16.0_SP*Uin(I+1,J)-Uin(I+2,J))
# else
       Uout(I,J) = MASK(I,J)*1.0_SP/12.0_SP/DX(I,J)/DX(I,J)*(-Uin(I-2,J)+16.0_SP*Uin(I-1,J)   &
                -30.0_SP*Uin(I,J)+16.0_SP*Uin(I+1,J)-Uin(I+2,J))
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_XX_HIGH

!-------------------------------------------------------------------------------------
!
!    DERIVATIVE_YY_HIGH is subroutine to calculation 4th-derivative yy
!
!    HISTORY: 
!    05/30/2011 Fengyan Shi
!
! --------------------------------------------------
SUBROUTINE DERIVATIVE_YY_HIGH(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DY,Uin,Uout)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: Uin
# if defined (CARTESIAN)
     REAL(SP),INTENT(IN) :: DY
# else
     REAL(SP),DIMENSION(M,N),INTENT(IN) :: DY
# endif
     INTEGER,DIMENSION(M,N),INTENT(IN) :: MASK
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Uout
     INTEGER :: I,J

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
# if defined (CARTESIAN)
       Uout(I,J) = MASK(I,J)*1.0_SP/12.0_SP/DY/DY*(-Uin(I,J-2)+16.0_SP*Uin(I,J-1)   &
                -30.0_SP*Uin(I,J)+16.0_SP*Uin(I,J+1)-Uin(I,J+2))
# else
       Uout(I,J) = MASK(I,J)*1.0_SP/12.0_SP/DY(I,J)/DY(I,J)*(-Uin(I,J-2)+16.0_SP*Uin(I,J-1)   &
                -30.0_SP*Uin(I,J)+16.0_SP*Uin(I,J+1)-Uin(I+2,J+2))
# endif
     ENDDO
     ENDDO

END SUBROUTINE DERIVATIVE_YY_HIGH

