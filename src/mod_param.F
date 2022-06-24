!------------------------------------------------------------------------------------
!
!      FILE mod_param.F
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
!    PARAM is the module to define all fixed parameters
!    
!    HISTORY: 
!      05/01/2010 Fengyan Shi
!          the module is updated corresponding to modifications in subroutines
!-------------------------------------------------------------------------------------

MODULE PARAM  
# if defined PARALLEL
       USE MPI
# endif
       IMPLICIT NONE
# if defined DOUBLE_PRECISION 
       INTEGER, PARAMETER::SP=8
# if defined PARALLEL
       INTEGER, PARAMETER::MPI_SP=MPI_DOUBLE_PRECISION
# endif
# else
       INTEGER, PARAMETER::SP=SELECTED_REAL_KIND(6,30)
# if defined PARALLEL
       INTEGER, PARAMETER::MPI_SP=MPI_REAL
# endif
# endif

! define parameters
       REAL(SP), PARAMETER::pi=3.141592653
       REAL(SP), PARAMETER::R_earth = 6371000.0_SP
       REAL(SP), PARAMETER::SMALL=0.000001_SP
       REAL(SP), PARAMETER::LARGE=999999.0_SP
       REAL(SP), PARAMETER:: grav=9.81_SP
       REAL(SP), PARAMETER:: zero = 0.0_SP
       REAL(SP), PARAMETER:: RHO_AW = 0.0012041_SP  ! relative to water
       REAL(SP), PARAMETER:: RHO_AIR = 1.15_SP  ! absolute value
       REAL(SP), PARAMETER:: RHO_WATER = 1000.0_SP
       REAL(SP), PARAMETER:: DEG2RAD = 0.0175_SP

# if defined (ITERATION)
! give a maximum iteration number of 2 like in Nwogu 
       INTEGER :: ITER_NUM = 2
       REAL(SP) :: ERRORi = 0.0001_SP
# endif

END MODULE PARAM
