!Crown Copyright 2014 AWE.
!
! This file is part of Fortran LULESH.
!
! Fortran LULESH is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! Fortran LULESH is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! Fortran LULESH. If not, see http://www.gnu.org/licenses/.
!
!
! Authors:
!   Duncan Harris
!   Andy Herdman
!
!
!
!Copyright (c) 2010.
!Lawrence Livermore National Security, LLC.
!Produced at the Lawrence Livermore National Laboratory.
!LLNL-CODE-461231
!All rights reserved.
!
!This file is part of LULESH, Version 1.0.
!Please also read this link -- http://www.opensource.org/licenses/index.php
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions
!are met:
!
!* Redistributions of source code must retain the above copyright
!notice, this list of conditions and the disclaimer below.
!
!* Redistributions in binary form must reproduce the above copyright
!notice, this list of conditions and the disclaimer (as noted below)
!in the documentation and/or other materials provided with the
!distribution.
!
!* Neither the name of the LLNS/LLNL nor the names of its contributors
!may be used to endorse or promote products derived from this software
!without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
!THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
!OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
!EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!Additional BSD Notice

!1. This notice is required to be provided under our contract with the U.S.
!Department of Energy (DOE). This work was produced at Lawrence Livermore
!National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
!
!2. Neither the United States Government nor Lawrence Livermore National
!Security, LLC nor any of their employees, makes any warranty, express
!or implied, or assumes any liability or responsibility for the accuracy,
!completeness, or usefulness of any information, apparatus, product, or
!process disclosed, or represents that its use would not infringe
!privately-owned rights.
!
!3. Also, reference herein to any specific commercial products, process, or
!services by trade name, trademark, manufacturer or otherwise does not
!necessarily constitute or imply its endorsement, recommendation, or
!favoring by the United States Government or Lawrence Livermore National
!Security, LLC. The views and opinions of authors expressed herein do not
!necessarily state or reflect those of the United States Government or
!Lawrence Livermore National Security, LLC, and shall not be used for
!advertising or product endorsement purposes.



PROGRAM lulesh

USE lulesh_domain


IMPLICIT NONE

INTEGER(KIND=4), PARAMETER :: VolumeError = -1
INTEGER(KIND=4), PARAMETER :: QStopError  = -2

INTEGER(KIND=4), PARAMETER :: RLK = 8





! Start of main
TYPE(domain_type) :: domain
INTEGER :: edgeElems 
INTEGER :: edgeNodes
REAL(KIND=8) :: tx, ty, tz 
INTEGER :: nidx, zidx 
INTEGER :: domElems 

INTEGER :: plane, row, col, i, j, k
INTEGER :: planeInc, rowInc

REAL(KIND=8),DIMENSION(0:7) :: x_local, y_local, z_local

INTEGER :: gnode, lnode, idx
INTEGER(KIND=4), DIMENSION(:), POINTER :: localNode => NULL()

REAL(KIND=8) :: volume
REAL(KIND=8) :: starttim, endtim
REAL(KIND=8) :: elapsed_time


! Stuff needed for boundary conditions
! 2 BCs on each of 6 hexahedral faces (12 bits)
INTEGER, PARAMETER :: XI_M        = z'003' ! 0x003
INTEGER, PARAMETER :: XI_M_SYMM   = z'001' ! 0x001
INTEGER, PARAMETER :: XI_M_FREE   = z'002' ! 0x002

INTEGER, PARAMETER :: XI_P        = z'00c' ! 0x00c
INTEGER, PARAMETER :: XI_P_SYMM   = z'004' ! 0x004
INTEGER, PARAMETER :: XI_P_FREE   = z'008' ! 0x008

INTEGER, PARAMETER :: ETA_M       = z'030' ! 0x030
INTEGER, PARAMETER :: ETA_M_SYMM  = z'010' ! 0x010
INTEGER, PARAMETER :: ETA_M_FREE  = z'020' ! 0x020

INTEGER, PARAMETER :: ETA_P       = z'0c0' ! 0x0c0
INTEGER, PARAMETER :: ETA_P_SYMM  = z'040' ! 0x040
INTEGER, PARAMETER :: ETA_P_FREE  = z'080' ! 0x080

INTEGER, PARAMETER :: ZETA_M      = z'300' ! 0x300
INTEGER, PARAMETER :: ZETA_M_SYMM = z'100' ! 0x100
INTEGER, PARAMETER :: ZETA_M_FREE = z'200' ! 0x200

INTEGER, PARAMETER :: ZETA_P      = z'c00' ! 0xc00
INTEGER, PARAMETER :: ZETA_P_SYMM = z'400' ! 0x400
INTEGER, PARAMETER :: ZETA_P_FREE = z'800' ! 0x800

CHARACTER(len=10) :: arg

INTEGER ::  ElemId = 0

REAL(KIND=8) :: MaxAbsDiff   = 0.0_RLK
REAL(KIND=8) :: TotalAbsDiff = 0.0_RLK
REAL(KIND=8) :: MaxRelDiff   = 0.0_RLK

REAL(KIND=8) :: AbsDiff, RelDiff


CALL GETARG(1, arg)
READ(arg,*) edgeElems

edgeNodes = edgeElems+1

! get run options to measure various metrics 

!****************************
!*   Initialize Sedov Mesh  *
!****************************

! construct a uniform box for this processor

domain%m_sizeX   = edgeElems 
domain%m_sizeY   = edgeElems 
domain%m_sizeZ   = edgeElems 
domain%m_numElem = edgeElems*edgeElems*edgeElems 
domain%m_numNode = edgeNodes*edgeNodes*edgeNodes 

domElems = domain%m_numElem 


! allocate field memory

CALL AllocateElemPersistent(domain%m_numElem)
CALL AllocateElemTemporary (domain%m_numElem) 

CALL AllocateNodalPersistent(domain%m_numNode) 
CALL AllocateNodesets(edgeNodes*edgeNodes) 


! initialize nodal coordinates 

nidx = 0
tz = 0.0_RLK

DO plane=0, edgeNodes-1
   ty = 0.0_RLK
   DO row=0, edgeNodes-1
      tx = 0.0_RLK
      DO col=0, edgeNodes-1
         domain%m_x(nidx) = tx
         domain%m_y(nidx) = ty
         domain%m_z(nidx) = tz
         nidx = nidx+1

         tx = (1.125_RLK*REAL((col+1),8))/REAL(edgeElems,8)
      END DO

      ty = 1.125_RLK*REAL((row+1),8)/REAL(edgeElems,8)
   END DO

   tz = 1.125_RLK*REAL((plane+1),8)/REAL(edgeElems,8)
END DO


! embed hexehedral elements in nodal point lattice

nidx = 0
zidx = 0

DO plane=0, edgeElems-1
   DO row=0, edgeElems-1
      DO col=0, edgeElems-1
         localNode => domain%m_nodelist(zidx*8:)
!  fortran pointer index starts from 1
         localNode(1) = nidx
         localNode(2) = nidx                                   + 1
         localNode(3) = nidx                       + edgeNodes + 1
         localNode(4) = nidx                       + edgeNodes
         localNode(5) = nidx + edgeNodes*edgeNodes
         localNode(6) = nidx + edgeNodes*edgeNodes             + 1
         localNode(7) = nidx + edgeNodes*edgeNodes + edgeNodes + 1
         localNode(8) = nidx + edgeNodes*edgeNodes + edgeNodes
         zidx = zidx + 1
         nidx = nidx + 1
      END DO
      nidx = nidx + 1
   END DO
   nidx = nidx + edgeNodes
END DO

NULLIFY(localNode)

CALL AllocateNodeElemIndexes()

!Create a material IndexSet (entire domain same material for now)
DO i=0, domElems-1
   domain%m_matElemlist(i) = i
END DO

  
! initialize material parameters
  domain%m_dtfixed         = -1.0e-7_RLK
  domain%m_deltatime       =  1.0e-7_RLK
  domain%m_deltatimemultlb =  1.1_RLK
  domain%m_deltatimemultub =  1.2_RLK
  domain%m_stoptime        =  1.0e-2_RLK
  domain%m_dtcourant       =  1.0e+20_RLK
  domain%m_dthydro         =  1.0e+20_RLK
  domain%m_dtmax           =  1.0e-2_RLK
  domain%m_time            =  0.0_RLK
  domain%m_cycle           =  0
  
  domain%m_e_cut = 1.0e-7_RLK
  domain%m_p_cut = 1.0e-7_RLK
  domain%m_q_cut = 1.0e-7_RLK
  domain%m_u_cut = 1.0e-7_RLK
  domain%m_v_cut = 1.0e-10_RLK
  
  domain%m_hgcoef  = 3.0_RLK
  domain%m_ss4o3   =(4.0_RLK)/(3.0_RLK)
  
  domain%m_qstop              = 1.0e+12_RLK
  domain%m_monoq_max_slope    = 1.0_RLK
  domain%m_monoq_limiter_mult = 2.0_RLK
  domain%m_qlc_monoq          = 0.5_RLK
  domain%m_qqc_monoq          = (2.0_RLK)/(3.0_RLK)
  domain%m_qqc                = 2.0_RLK
  
  domain%m_pmin =  0.0_RLK
  domain%m_emin = -1.0e+15_RLK
  
  domain%m_dvovmax =  0.1_RLK
  
  domain%m_eosvmax =  1.0e+9_RLK
  domain%m_eosvmin =  1.0e-9_RLK
  
  domain%m_refdens =  1.0_RLK

 ! initialize field data
 DO i=0, domElems-1
    DO lnode=0,7
       gnode = domain%m_nodelist(i*8+lnode)
       x_local(lnode) = domain%m_x(gnode)
       y_local(lnode) = domain%m_y(gnode)
       z_local(lnode) = domain%m_z(gnode)
    END DO
   
    ! volume calculations
    volume = CalcElemVolume(x_local, y_local, z_local)
    domain%m_volo(i) = volume
    domain%m_elemMass(i) = volume
    DO j=0, 7
       idx = domain%m_nodelist(i*8+j)
       domain%m_nodalMass(idx) =  domain%m_nodalMass(idx) + ( volume / 8.0_RLK)
    END DO
 END DO
   

 ! deposit energy 
 domain%m_e(0) = 3.948746e+7
  
 ! set up symmetry nodesets
 nidx = 0
 
 DO i=0,edgeNodes-1
    planeInc = i*edgeNodes*edgeNodes
    rowInc   = i*edgeNodes
    DO j=0,edgeNodes-1
       domain%m_symmX(nidx) = planeInc + j*edgeNodes
       domain%m_symmY(nidx) = planeInc + j
       domain%m_symmZ(nidx) = rowInc   + j
       nidx=nidx+1
    END DO
 END DO

 ! set up elemement connectivity information
 domain%m_lxim(0) = 0
 DO i=1,domElems-1
    domain%m_lxim(i)   = i-1
    domain%m_lxip(i-1) = i
 END DO

 domain%m_lxip(domElems-1) = domElems

 DO i=0, edgeElems-1
    domain%m_letam(i)=i
    domain%m_letap(domElems-edgeElems+i) = domElems-edgeElems+i
 END DO

DO i=edgeElems,domElems-1
   domain%m_letam(i) = i-edgeElems
   domain%m_letap(i-edgeElems) = i
END DO

DO i=0,edgeElems*edgeElems-1
   domain%m_lzetam(i) = i
   domain%m_lzetap(domElems-edgeElems*edgeElems+i) = domElems-edgeElems*edgeElems+i
END DO

DO i=(edgeElems*edgeElems), domElems-1
   domain%m_lzetam(i) = i - edgeElems*edgeElems
   domain%m_lzetap(i-edgeElems*edgeElems) = i
END DO

! set up boundary condition information
domain%m_elemBC = 0 ! clear BCs by default

! faces on "external" boundaries will be
! symmetry plane or free surface BCs
DO i=0,edgeElems-1
   planeInc = i*edgeElems*edgeElems
   rowInc   = i*edgeElems
   DO j=0,edgeElems-1
      domain%m_elemBC(planeInc+j*edgeElems)                 = IOR(domain%m_elemBC(planeInc+(j)*edgeElems), XI_M_SYMM)
      domain%m_elemBC(planeInc+(j)*edgeElems+edgeElems-1)     = IOR(domain%m_elemBC(planeInc+(j)*edgeElems+edgeElems-1), XI_P_FREE)
      domain%m_elemBC(planeInc+j)                               = IOR(domain%m_elemBC(planeInc+j),ETA_M_SYMM)
      domain%m_elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) = IOR(domain%m_elemBC(planeInc+j+edgeElems*edgeElems-edgeElems),ETA_P_FREE) 
      domain%m_elemBC(rowInc+j)                                 = IOR(domain%m_elemBC(rowInc+j),ZETA_M_SYMM) 
      domain%m_elemBC(rowInc+j+domElems-edgeElems*edgeElems)    = IOR(domain%m_elemBC(rowInc+j+domElems-edgeElems*edgeElems),ZETA_P_FREE)
   END DO
END DO


! timestep to solution
!!$ timeval start, end
!!$ gettimeofday(&start, NULL)
CALL CPU_TIME(starttim)

DO
   call TimeIncrement
   CALL LagrangeLeapFrog

#ifdef LULESH_SHOW_PROGRESS
   PRINT *,"time = ", domain%m_time, " dt=",domain%m_deltatime
#endif

   IF(domain%m_time >= domain%m_stoptime) EXIT
END DO


CALL CPU_TIME(endtim)
!!$  gettimeofday(&end, NULL);

elapsed_time = endtim - starttim
!!$  double elapsed_time = double(end.tv_sec - start.tv_sec) + double(end.tv_usec - start.tv_usec) *1e-6;
!!$  

PRINT *,""
PRINT *,""
PRINT '("Elapsed time = ", E13.6)', elapsed_time


ElemId = 0

PRINT *,"Run completed:"
PRINT '("   Problem size        = ", I8)',    edgeElems
PRINT '("   Iteration count     = ", I8)',    domain%m_cycle
PRINT '("   Final Origin Energy = ", e13.6)', domain%m_e(ElemId)
PRINT *,""

  
 MaxAbsDiff = 0.0_RLK
 TotalAbsDiff = 0.0_RLK
 MaxRelDiff = 0.0_RLK


! MIGHT WANT TO DOUBLE CHECK THESE LOOPS
 DO j=0, edgeElems-1
    DO k=j+1, edgeElems-1
       AbsDiff = ABS(domain%m_e(j*edgeElems+k) - domain%m_e(k*edgeElems+j))
       TotalAbsDiff  = TotalAbsDiff+AbsDiff

       if (MaxAbsDiff <AbsDiff) MaxAbsDiff = AbsDiff

       RelDiff = AbsDiff / domain%m_e(k*edgeElems+j)

       if (MaxRelDiff <RelDiff)  MaxRelDiff = RelDiff

    END DO
 END DO

 PRINT *,"  Testing Plane 0 of Energy Array:"
 PRINT '("        MaxAbsDiff   = ", E13.6)', MaxAbsDiff
 PRINT '("        TotalAbsDiff = ", E13.6)', TotalAbsDiff
 PRINT '("        MaxRelDiff   = ", E13.6)', MaxRelDiff
 PRINT *,""

CONTAINS



REAL(KIND=8) FUNCTION CBRT(dat)

  IMPLICIT NONE
  REAL(KIND=8) :: dat

  CBRT = dat**(1.0_RLK/3.0_RLK)

END FUNCTION CBRT


REAL(KIND=8) FUNCTION CalcElemVolume( x, y, z )

  IMPLICIT NONE
  REAL(KIND=8) ,DIMENSION(0:7) :: x, y, z

  REAL(KIND=8)  :: volume=0.0_RLK

!!$REAL(KIND=8)  :: x0, x1, x2, x3, x4, x5, x6, x7
!!$REAL(KIND=8)  :: y0, y1, y2, y3, y4, y5, y6, y7
!!$REAL(KIND=8)  :: z0, z1, z2, z3, z4, z5, z6, z7

  REAL(KIND=8) :: twelveth = (1.0_RLK)/(12.0_RLK)

  REAL(KIND=8) :: dx61
  REAL(KIND=8) :: dy61
  REAL(KIND=8) :: dz61

  REAL(KIND=8) :: dx70
  REAL(KIND=8) :: dy70
  REAL(KIND=8) :: dz70
  
  REAL(KIND=8) :: dx63
  REAL(KIND=8) :: dy63 
  REAL(KIND=8) :: dz63

  REAL(KIND=8) :: dx20 
  REAL(KIND=8) :: dy20
  REAL(KIND=8) :: dz20

  REAL(KIND=8) :: dx50 
  REAL(KIND=8) :: dy50
  REAL(KIND=8) :: dz50

  REAL(KIND=8) :: dx64
  REAL(KIND=8) :: dy64
  REAL(KIND=8) :: dz64

  REAL(KIND=8) :: dx31
  REAL(KIND=8) :: dy31 
  REAL(KIND=8) :: dz31

  REAL(KIND=8) :: dx72
  REAL(KIND=8) :: dy72
  REAL(KIND=8) :: dz72

  REAL(KIND=8) :: dx43 
  REAL(KIND=8) :: dy43
  REAL(KIND=8) :: dz43

  REAL(KIND=8) :: dx57
  REAL(KIND=8) :: dy57
  REAL(KIND=8) :: dz57

  REAL(KIND=8) :: dx14
  REAL(KIND=8) :: dy14 
  REAL(KIND=8) :: dz14

  REAL(KIND=8) :: dx25
  REAL(KIND=8) :: dy25 
  REAL(KIND=8) :: dz25


  dx61 = x(6) - x(1)
  dy61 = y(6) - y(1)
  dz61 = z(6) - z(1)

  dx70 = x(7) - x(0)
  dy70 = y(7) - y(0)
  dz70 = z(7) - z(0)

  dx63 = x(6) - x(3)
  dy63 = y(6) - y(3)
  dz63 = z(6) - z(3)

  dx20 = x(2) - x(0)
  dy20 = y(2) - y(0)
  dz20 = z(2) - z(0)

  dx50 = x(5) - x(0)
  dy50 = y(5) - y(0)
  dz50 = z(5) - z(0)

  dx64 = x(6) - x(4)
  dy64 = y(6) - y(4)
  dz64 = z(6) - z(4)

  dx31 = x(3) - x(1)
  dy31 = y(3) - y(1)
  dz31 = z(3) - z(1)

  dx72 = x(7) - x(2)
  dy72 = y(7) - y(2)
  dz72 = z(7) - z(2)

  dx43 = x(4) - x(3)
  dy43 = y(4) - y(3)
  dz43 = z(4) - z(3)

  dx57 = x(5) - x(7)
  dy57 = y(5) - y(7)
  dz57 = z(5) - z(7)

  dx14 = x(1) - x(4)
  dy14 = y(1) - y(4)
  dz14 = z(1) - z(4)

  dx25 = x(2) - x(5)
  dy25 = y(2) - y(5)
  dz25 = z(2) - z(5)

  volume =  TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,   &
                           dy31 + dy72, dy63, dy20,   &
                           dz31 + dz72, dz63, dz20) + &
            TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,   &
                           dy43 + dy57, dy64, dy70,   &
                           dz43 + dz57, dz64, dz70) + &
            TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,   &
                           dy14 + dy25, dy61, dy50,   &
                           dz14 + dz25, dz61, dz50)

  volume = volume*twelveth


  CalcElemVolume=volume


END FUNCTION CalcElemVolume


REAL(KIND=8) FUNCTION TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3)

  REAL(KIND=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3

  TRIPLE_PRODUCT = ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3)  &
                  - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  RETURN

END FUNCTION TRIPLE_PRODUCT


FUNCTION AreaFace( x0, x1, x2, x3,  &
                   y0, y1, y2, y3,  &
                   z0, z1, z2, z3  ) RESULT(area)


  IMPLICIT NONE
  REAL(KIND=8)  :: x0, x1, x2, x3
  REAL(KIND=8)  :: y0, y1, y2, y3
  REAL(KIND=8)  :: z0, z1, z2, z3

  REAL(KIND=8) :: fx, fy, fz
  REAL(KIND=8) :: gx, gy, gz
  REAL(KIND=8) :: area

  fx = (x2 - x0) - (x3 - x1)
  fy = (y2 - y0) - (y3 - y1)
  fz = (z2 - z0) - (z3 - z1)
  gx = (x2 - x0) + (x3 - x1)
  gy = (y2 - y0) + (y3 - y1)
  gz = (z2 - z0) + (z3 - z1)

  area =                             &
    (fx * fx + fy * fy + fz * fz) *  &
    (gx * gx + gy * gy + gz * gz) -  &
    (fx * gx + fy * gy + fz * gz) *  &
    (fx * gx + fy * gy + fz * gz)

  RETURN

END FUNCTION AreaFace



REAL(KIND=8) FUNCTION SUM4(a, b, c, d)

  IMPLICIT NONE 
  REAL(KIND=8) :: a, b, c, d

  SUM4 = a + b + c + d

  RETURN

END FUNCTION SUM4


END PROGRAM lulesh
