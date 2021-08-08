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

! Import the communcation layer
use lulesh_comm

IMPLICIT NONE

INTEGER(KIND=4), PARAMETER :: VolumeError = -1
INTEGER(KIND=4), PARAMETER :: QStopError  = -2

INTEGER(KIND=4), PARAMETER :: RLK = 8

TYPE domain_type


  
  !****************
  ! Implementation 
  !****************
  
  ! Node-centered 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_x   ! coordinates 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_y 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_z 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_xd  ! velocities 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_yd 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_zd 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_xdd  ! accelerations 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_ydd 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_zdd 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_fx   ! forces 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_fy 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_fz 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_nodalMass   ! mass 
  
  INTEGER,     DIMENSION(:), ALLOCATABLE ::m_symmX   ! symmetry plane nodesets 
  INTEGER,     DIMENSION(:), ALLOCATABLE ::m_symmY 
  INTEGER,     DIMENSION(:), ALLOCATABLE ::m_symmZ 
  
  INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemCount 
  INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemStart 
  !   INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemList 
  INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemCornerList 
  
  ! Element-centered 
  
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_matElemlist   ! material indexset 
  INTEGER,     DIMENSION(:), POINTER     :: m_nodelist => NULL()  ! elemToNode connectivity 
  
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lxim   ! element connectivity across each face 
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lxip 
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_letam 
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_letap 
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lzetam 
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lzetap 
  
  INTEGER,     DIMENSION(:), ALLOCATABLE :: m_elemBC   ! symmetry/free-surface flags for each elem face 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_dxx   ! principal strains -- temporary 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_dyy 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_dzz 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv_xi     ! velocity gradient -- temporary 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv_eta 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv_zeta 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delx_xi     ! coordinate gradient -- temporary 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delx_eta 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delx_zeta 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_e    ! energy 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_p    ! pressure 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_q    ! q 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_ql   ! linear term for q 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_qq   ! quadratic term for q 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_v      ! relative volume 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_volo   ! reference volume 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_vnew   ! new relative volume -- temporary 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv   ! m_vnew - m_v 
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_vdov   ! volume derivative over volume 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_arealg   ! characteristic length of an element 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_ss       ! "sound speed" 
  
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_elemMass   ! mass 

  ! Parameters 
  
  REAL(KIND=8)       ::   m_dtfixed            ! fixed time increment 
  REAL(KIND=8)       ::   m_time               ! current time 
  REAL(KIND=8)       ::   m_deltatime          ! variable time increment 
  REAL(KIND=8)       ::   m_deltatimemultlb 
  REAL(KIND=8)       ::   m_deltatimemultub 
  REAL(KIND=8)       ::   m_stoptime           ! end time for simulation 
  
  REAL(KIND=8)       ::   m_u_cut              ! velocity tolerance 
  REAL(KIND=8)       ::   m_hgcoef             ! hourglass control 
  REAL(KIND=8)       ::   m_qstop              ! excessive q indicator 
  REAL(KIND=8)       ::   m_monoq_max_slope 
  REAL(KIND=8)       ::   m_monoq_limiter_mult 
  REAL(KIND=8)       ::   m_e_cut              ! energy tolerance 
  REAL(KIND=8)       ::   m_p_cut              ! pressure tolerance 
  REAL(KIND=8)       ::   m_ss4o3 
  REAL(KIND=8)       ::   m_q_cut              ! q tolerance 
  REAL(KIND=8)       ::   m_v_cut              ! relative volume tolerance 
  REAL(KIND=8)       ::   m_qlc_monoq          ! linear term coef for q 
  REAL(KIND=8)       ::   m_qqc_monoq          ! quadratic term coef for q 
  REAL(KIND=8)       ::   m_qqc 
  REAL(KIND=8)       ::   m_eosvmax 
  REAL(KIND=8)       ::   m_eosvmin 
  REAL(KIND=8)       ::   m_pmin               ! pressure floor 
  REAL(KIND=8)       ::   m_emin               ! energy floor 
  REAL(KIND=8)       ::   m_dvovmax            ! maximum allowable volume change 
  REAL(KIND=8)       ::   m_refdens            ! reference density 
  
  REAL(KIND=8)       ::   m_dtcourant          ! courant constraint 
  REAL(KIND=8)       ::   m_dthydro            ! volume change constraint 
  REAL(KIND=8)       ::   m_dtmax              ! maximum allowable time increment 
  
  INTEGER ::   m_cycle              ! iteration count for simulation 
  
  INTEGER ::   m_sizeX            ! X,Y,Z extent of this block 
  INTEGER ::   m_sizeY 
  INTEGER ::   m_sizeZ 
  
  INTEGER ::   m_numElem          ! Elements/Nodes in this domain 
  INTEGER ::   m_numNode 

  
END TYPE domain_type



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

SUBROUTINE AllocateNodalPersistent(size)

  IMPLICIT NONE 
  INTEGER :: size

  ALLOCATE(domain%m_x(0:size-1))
  ALLOCATE(domain%m_y(0:size-1)) 
  ALLOCATE(domain%m_z(0:size-1)) 
  
  ALLOCATE(domain%m_xd(0:size-1)) 
  ALLOCATE(domain%m_yd(0:size-1)) 
  ALLOCATE(domain%m_zd(0:size-1)) 
  domain%m_xd = 0.0_RLK
  domain%m_yd = 0.0_RLK
  domain%m_zd = 0.0_RLK

  ALLOCATE(domain%m_xdd(0:size-1)) 
  ALLOCATE(domain%m_ydd(0:size-1)) 
  ALLOCATE(domain%m_zdd(0:size-1)) 
  domain%m_xdd = 0.0_RLK
  domain%m_ydd = 0.0_RLK
  domain%m_zdd = 0.0_RLK
  
  ALLOCATE(domain%m_fx(0:size-1)) 
  ALLOCATE(domain%m_fy(0:size-1)) 
  ALLOCATE(domain%m_fz(0:size-1)) 
  
  ALLOCATE(domain%m_nodalMass(0:size-1))
  domain%m_nodalMass = 0.0_RLK

END SUBROUTINE AllocateNodalPersistent




SUBROUTINE AllocateElemPersistent(size)
  IMPLICIT NONE
  INTEGER :: size 
  
  ALLOCATE(domain%m_matElemlist(0:size-1)) 
  ALLOCATE(domain%m_nodelist(0:8*size-1)) 
  
  ALLOCATE(domain%m_lxim(0:size-1)) 
  ALLOCATE(domain%m_lxip(0:size-1)) 
  ALLOCATE(domain%m_letam(0:size-1)) 
  ALLOCATE(domain%m_letap(0:size-1)) 
  ALLOCATE(domain%m_lzetam(0:size-1)) 
  ALLOCATE(domain%m_lzetap(0:size-1)) 
    
  ALLOCATE(domain%m_elemBC(0:size-1)) 
  
  ALLOCATE(domain%m_e(0:size-1)) 
  domain%m_e = 0.0_RLK
  ALLOCATE(domain%m_p(0:size-1)) 
  domain%m_p = 0.0_RLK
  ALLOCATE(domain%m_q(0:size-1)) 
  ALLOCATE(domain%m_ql(0:size-1)) 
  ALLOCATE(domain%m_qq(0:size-1)) 
  
  ALLOCATE(domain%m_v(0:size-1)) 
  domain%m_v = 1.0_RLK
  ALLOCATE(domain%m_volo(0:size-1)) 
  ALLOCATE(domain%m_delv(0:size-1)) 
  ALLOCATE(domain%m_vdov(0:size-1)) 
  
  ALLOCATE(domain%m_arealg(0:size-1)) 
  
  ALLOCATE(domain%m_ss(0:size-1)) 
  
  ALLOCATE(domain%m_elemMass(0:size-1)) 
  
END SUBROUTINE AllocateElemPersistent





!Temporaries should not be initialized in bulk but
!this is a runnable placeholder for now
SUBROUTINE AllocateElemTemporary( size)
  IMPLICIT NONE 
  INTEGER :: size
  
  ALLOCATE(domain%m_dxx(0:size-1))
  ALLOCATE(domain%m_dyy(0:size-1))
  ALLOCATE(domain%m_dzz(0:size-1))
  
  ALLOCATE(domain%m_delv_xi(0:size-1))
  ALLOCATE(domain%m_delv_eta(0:size-1))
  ALLOCATE(domain%m_delv_zeta(0:size-1))
  
  ALLOCATE(domain%m_delx_xi(0:size-1))
  ALLOCATE(domain%m_delx_eta(0:size-1))
  ALLOCATE(domain%m_delx_zeta(0:size-1))
  
  ALLOCATE(domain%m_vnew(0:size-1))
  
END SUBROUTINE AllocateElemTemporary


SUBROUTINE AllocateNodesets( size)
  IMPLICIT NONE 
  INTEGER :: size

  ALLOCATE(domain%m_symmX(0:size-1))
  ALLOCATE(domain%m_symmY(0:size-1))
  ALLOCATE(domain%m_symmZ(0:size-1))

END SUBROUTINE AllocateNodesets





SUBROUTINE AllocateNodeElemIndexes()
  IMPLICIT NONE 

  INTEGER :: m
  INTEGER :: i,j,k
  INTEGER :: offset
  INTEGER :: clSize, clv
  INTEGER :: numElem 
  INTEGER :: numNode
  INTEGER :: nodelist_len

  numElem = domain%m_numElem
  numNode = domain%m_numNode

  ! set up node-centered indexing of elements
  ALLOCATE(domain%m_nodeElemCount(0:numNode-1))
  !  m_nodeElemCount.resize(numNode);

  domain%m_nodeElemCount=0

  DO i=0, SIZE(domain%m_nodelist)-1  
     domain%m_nodeElemCount(domain%m_nodelist(i))=   &
          domain%m_nodeElemCount(domain%m_nodelist(i))+1
  END DO

  ALLOCATE(domain%m_nodeElemStart(0:numNode-1))  

  domain%m_nodeElemStart=0

  DO i=1,numNode-1
     domain%m_nodeElemStart(i) = domain%m_nodeElemStart(i-1) + domain%m_nodeElemCount(i-1)
  END DO

  ALLOCATE(domain%m_nodeElemCornerList(0:(domain%m_nodeElemStart(numNode-1) +  &
                                          domain%m_nodeElemCount(numNode-1))))

  domain%m_nodeElemCount=0  


  !CHECK THIS LOOP - THINK IT'S OK ????

  DO i=0, SIZE(domain%m_nodelist)-1
        m=domain%m_nodelist(i)
        offset = domain%m_nodeElemStart(m)+domain%m_nodeElemCount(m)
        domain%m_nodeElemCornerList(offset) = i
        domain%m_nodeElemCount(m) = domain%m_nodeElemCount(m) + 1
  END DO

  clSize = SIZE(domain%m_nodeElemCornerList)
  DO i=0, clSize-1
     clv=domain%m_nodeElemCornerList(i)
     IF ((clv.LT.0).OR.(clv.GT.numElem*8))THEN
        PRINT*,"AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!"
        CALL luabort(1)
     END IF
  END DO


END SUBROUTINE AllocateNodeElemIndexes


SUBROUTINE TimeIncrement()
  IMPLICIT NONE 

  REAL(KIND=8) :: targetdt
  REAL(KIND=8) :: ratio, olddt, newdt

  targetdt = domain%m_stoptime - domain%m_time

  IF (( domain%m_dtfixed <= 0.0_RLK) .AND. (domain%m_cycle /= 0)) THEN

     olddt = domain%m_deltatime
     
     ! This will require a reduction in parallel
     newdt = 1.0e+20

     IF (domain%m_dtcourant < newdt) THEN
        newdt = domain%m_dtcourant / 2.0_RLK
     END IF

     IF (domain%m_dthydro < newdt) THEN
        newdt = domain%m_dthydro * (2.0_RLK/3.0_RLK)
     END IF

     ratio = newdt / olddt

     IF (ratio >= 1.0_RLK) THEN
        IF (ratio < domain%m_deltatimemultlb) THEN
           newdt = olddt
        ELSE IF (ratio > domain%m_deltatimemultub) THEN
           newdt = olddt*domain%m_deltatimemultub
        END IF
     END IF


     IF (newdt > domain%m_dtmax) THEN
        newdt = domain%m_dtmax
     END IF

     domain%m_deltatime = newdt

  END IF


  ! TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE
  IF ((targetdt > domain%m_deltatime) .AND. (targetdt < 4.0_RLK * domain%m_deltatime / 3.0_RLK)) THEN
     targetdt = 2.0_RLK * domain%m_deltatime / 3.0_RLK
  END IF

  IF (targetdt < domain%m_deltatime) THEN
     domain%m_deltatime = targetdt
  END IF

  domain%m_time = domain%m_time+domain%m_deltatime

  domain%m_cycle=domain%m_cycle+1

END SUBROUTINE TimeIncrement







SUBROUTINE InitStressTermsForElems( numElem,  sigxx, sigyy, sigzz)

  IMPLICIT NONE

  INTEGER         :: numElem
  REAL(KIND=8), DIMENSION(0:) :: sigxx
  REAL(KIND=8), DIMENSION(0:) :: sigyy
  REAL(KIND=8), DIMENSION(0:) :: sigzz
  INTEGER(KIND=4) :: ii

  DO ii = 0, numElem-1
    sigxx(ii) =  - domain%m_p(ii) - domain%m_q(ii)
    sigyy(ii) =  - domain%m_p(ii) - domain%m_q(ii)
    sigzz(ii) =  - domain%m_p(ii) - domain%m_q(ii)
  ENDDO

END SUBROUTINE InitStressTermsForElems





SUBROUTINE CalcElemShapeFunctionDerivatives( x, y, z,   &
                                             b,         &
                                             el_volume   )
  IMPLICIT NONE 

  REAL(KIND=8), DIMENSION(0:7)  :: x, y, z
  REAL(KIND=8), DIMENSION(0:7,0:2) :: b ! alloc 2nd dim to 8 or 0:7
  REAL(KIND=8), INTENT(INOUT) :: el_volume
  REAL(KIND=8)  :: x0,x1,x2,x3,x4,x5,x6,x7
  REAL(KIND=8)  :: y0,y1,y2,y3,y4,y5,y6,y7
  REAL(KIND=8)  :: z0,z1,z2,z3,z4,z5,z6,z7

  REAL(KIND=8)  :: fjxxi, fjxet, fjxze
  REAL(KIND=8)  :: fjyxi, fjyet, fjyze
  REAL(KIND=8)  :: fjzxi, fjzet, fjzze
  REAL(KIND=8)  :: cjxxi, cjxet, cjxze
  REAL(KIND=8)  :: cjyxi, cjyet, cjyze
  REAL(KIND=8)  :: cjzxi, cjzet, cjzze

  x0 = x(0)
  x1 = x(1)
  x2 = x(2)
  x3 = x(3)
  x4 = x(4)
  x5 = x(5)
  x6 = x(6)
  x7 = x(7)

  y0 = y(0)
  y1 = y(1)
  y2 = y(2)
  y3 = y(3)
  y4 = y(4)
  y5 = y(5)
  y6 = y(6)
  y7 = y(7)

  z0 = z(0)
  z1 = z(1)
  z2 = z(2)
  z3 = z(3)
  z4 = z(4)
  z5 = z(5)
  z6 = z(6)
  z7 = z(7)

  fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) )
  fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) )
  fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) )

  fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) )
  fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) )
  fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) )

  fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) )
  fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) )
  fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) )

! compute cofactors
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze)
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze)
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet)

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze)
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze)
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet)

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze)
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze)
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet)

! calculate partials :
!     this need only be done for l = 0,1,2,3   since , by symmetry ,
!     (6,7,4,5) = - (0,1,2,3) .
  b(0,0) =   -  cjxxi  -  cjxet  -  cjxze
  b(1,0) =      cjxxi  -  cjxet  -  cjxze
  b(2,0) =      cjxxi  +  cjxet  -  cjxze
  b(3,0) =   -  cjxxi  +  cjxet  -  cjxze
  b(4,0) = -b(2,0)
  b(5,0) = -b(3,0)
  b(6,0) = -b(0,0)
  b(7,0) = -b(1,0)

  b(0,1) =   -  cjyxi  -  cjyet  -  cjyze
  b(1,1) =      cjyxi  -  cjyet  -  cjyze
  b(2,1) =      cjyxi  +  cjyet  -  cjyze
  b(3,1) =   -  cjyxi  +  cjyet  -  cjyze
  b(4,1) = -b(2,1)
  b(5,1) = -b(3,1)
  b(6,1) = -b(0,1)
  b(7,1) = -b(1,1)

  b(0,2) =   -  cjzxi  -  cjzet  -  cjzze
  b(1,2) =      cjzxi  -  cjzet  -  cjzze
  b(2,2) =      cjzxi  +  cjzet  -  cjzze
  b(3,2) =   -  cjzxi  +  cjzet  -  cjzze
  b(4,2) = -b(2,2)
  b(5,2) = -b(2,3)
  b(6,2) = -b(2,0)
  b(7,2) = -b(2,1)

! calculate jacobian determinant (volume)
  el_volume = 8.0_RLK * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet)

END SUBROUTINE CalcElemShapeFunctionDerivatives





SUBROUTINE SumElemFaceNormal(normalX0, normalY0, normalZ0, &
                             normalX1, normalY1, normalZ1, &
                             normalX2, normalY2, normalZ2, &
                             normalX3, normalY3, normalZ3, &
                              x0,  y0,  z0,    &
                              x1,  y1,  z1,    &
                              x2,  y2,  z2,    &
                              x3,  y3,  z3     )

  IMPLICIT NONE

  REAL(KIND=8) :: normalX0,normalY0,normalZ0
  REAL(KIND=8) :: normalX1,normalY1,normalZ1
  REAL(KIND=8) :: normalX2,normalY2,normalZ2
  REAL(KIND=8) :: normalX3,normalY3,normalZ3
  REAL(KIND=8) :: x0,y0,z0
  REAL(KIND=8) :: x1,y1,z1
  REAL(KIND=8) :: x2,y2,z2  
  REAL(KIND=8) :: x3,y3,z3
 
  REAL(KIND=8) :: bisectX0
  REAL(KIND=8) :: bisectY0
  REAL(KIND=8) :: bisectZ0
  REAL(KIND=8) :: bisectX1
  REAL(KIND=8) :: bisectY1
  REAL(KIND=8) :: bisectZ1
  REAL(KIND=8) :: areaX
  REAL(KIND=8) :: areaY
  REAL(KIND=8) :: areaZ
  REAL(KIND=8), PARAMETER :: RHALF = 0.5_RLK
  REAL(KIND=8), PARAMETER :: RQTR  = 0.25_RLK

  bisectX0 = RHALF * (x3 + x2 - x1 - x0)
  bisectY0 = RHALF * (y3 + y2 - y1 - y0)
  bisectZ0 = RHALF * (z3 + z2 - z1 - z0)
  bisectX1 = RHALF * (x2 + x1 - x3 - x0)
  bisectY1 = RHALF * (y2 + y1 - y3 - y0)
  bisectZ1 = RHALF * (z2 + z1 - z3 - z0)
  areaX = RQTR * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1)
  areaY = RQTR * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1)
  areaZ = RQTR * (bisectX0 * bisectY1 - bisectY0 * bisectX1)

  normalX0 = normalX0 + areaX
  normalX1 = normalX1 + areaX
  normalX2 = normalX2 + areaX
  normalX3 = normalX3 + areaX

  normalY0 = normalY0 + areaY
  normalY1 = normalY1 + areaY
  normalY2 = normalY2 + areaY
  normalY3 = normalY3 + areaY

  normalZ0 = normalZ0 + areaZ
  normalZ1 = normalZ1 + areaZ
  normalZ2 = normalZ2 + areaZ
  normalZ3 = normalZ3 + areaZ

END SUBROUTINE SumElemFaceNormal







SUBROUTINE CalcElemNodeNormals( pfx,pfy, pfz, x, y, z  )

  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:) :: pfx,pfy,pfz
  REAL(KIND=8), DIMENSION(0:) :: x, y, z 

  DO i = 0, 7
    pfx(i) = 0.0_RLK
    pfy(i) = 0.0_RLK
    pfz(i) = 0.0_RLK
  ENDDO

! evaluate face one: nodes 0, 1, 2, 3
  CALL SumElemFaceNormal(pfx(0), pfy(0), pfz(0),              &
                         pfx(1), pfy(1), pfz(1),              &
                         pfx(2), pfy(2), pfz(2),              &
                         pfx(3), pfy(3), pfz(3),              &
                         x(0), y(0), z(0), x(1), y(1), z(1),  &
                         x(2), y(2), z(2), x(3), y(3), z(3))
! evaluate face two: nodes 0, 4, 5, 1 */
  CALL SumElemFaceNormal(pfx(0), pfy(0), pfz(0),              &
                         pfx(4), pfy(4), pfz(4),              &
                         pfx(5), pfy(5), pfz(5),              &
                         pfx(1), pfy(1), pfz(1),              &
                         x(0), y(0), z(0), x(4), y(4), z(4),  &
                         x(5), y(5), z(5), x(1), y(1), z(1))
! evaluate face three: nodes 1, 5, 6, 2 */
  CALL SumElemFaceNormal(pfx(1), pfy(1), pfz(1),              &
                         pfx(5), pfy(5), pfz(5),              &
                         pfx(6), pfy(6), pfz(6),              &
                         pfx(2), pfy(2), pfz(2),              &
                         x(1), y(1), z(1), x(5), y(5), z(5),  &
                         x(6), y(6), z(6), x(2), y(2), z(2))
! evaluate face four: nodes 2, 6, 7, 3 */
  CALL SumElemFaceNormal(pfx(2), pfy(2), pfz(2),              &
                         pfx(6), pfy(6), pfz(6),              &
                         pfx(7), pfy(7), pfz(7),              &
                         pfx(3), pfy(3), pfz(3),              &
                         x(2), y(2), z(2), x(6), y(6), z(6),  &
                         x(7), y(7), z(7), x(3), y(3), z(3))
! evaluate face five: nodes 3, 7, 4, 0 */
  CALL SumElemFaceNormal(pfx(3), pfy(3), pfz(3),              &
                         pfx(7), pfy(7), pfz(7),              &
                         pfx(4), pfy(4), pfz(4),              &
                         pfx(0), pfy(0), pfz(0),              &
                         x(3), y(3), z(3), x(7), y(7), z(7),  &
                         x(4), y(4), z(4), x(0), y(0), z(0))
! evaluate face six: nodes 4, 7, 6, 5 */
  CALL SumElemFaceNormal(pfx(4), pfy(4), pfz(4),              &
                         pfx(7), pfy(7), pfz(7),              &
                         pfx(6), pfy(6), pfz(6),              &
                         pfx(5), pfy(5), pfz(5),              &
                         x(4), y(4), z(4), x(7), y(7), z(7),  &
                         x(6), y(6), z(6), x(5), y(5), z(5))

END SUBROUTINE CalcElemNodeNormals




SUBROUTINE SumElemStressesToNodeForces(B, stress_xx, stress_yy, stress_zz,  fx,  fy,  fz)

  IMPLICIT NONE 
  REAL(KIND=8) ,DIMENSION(0:7,0:2) :: B ! alloc 2nd dim to 8 or 0:7
  REAL(KIND=8) :: stress_xx, stress_yy, stress_zz
  REAL(KIND=8), DIMENSION(0:7) ::  fx,  fy,  fz 

  REAL(KIND=8) :: pfx0, pfx1, pfx2, pfx3, pfx4, pfx5, pfx6, pfx7
  REAL(KIND=8) :: pfy0, pfy1, pfy2, pfy3, pfy4, pfy5, pfy6, pfy7
  REAL(KIND=8) :: pfz0, pfz1, pfz2, pfz3, pfz4, pfz5, pfz6, pfz7

  pfx0 = B(0,0)
  pfx1 = B(1,0)
  pfx2 = B(2,0)
  pfx3 = B(3,0)
  pfx4 = B(4,0)
  pfx5 = B(5,0)
  pfx6 = B(6,0)
  pfx7 = B(7,0)

  pfy0 = B(0,1)
  pfy1 = B(1,1)
  pfy2 = B(2,1)
  pfy3 = B(3,1)
  pfy4 = B(4,1)
  pfy5 = B(5,1)
  pfy6 = B(6,1)
  pfy7 = B(7,1)

  pfz0 = B(0,2)
  pfz1 = B(1,2)
  pfz2 = B(2,2)
  pfz3 = B(3,2)
  pfz4 = B(4,2)
  pfz5 = B(5,2)
  pfz6 = B(6,2)
  pfz7 = B(7,2)

  fx(0) = -( stress_xx * pfx0 )
  fx(1) = -( stress_xx * pfx1 )
  fx(2) = -( stress_xx * pfx2 )
  fx(3) = -( stress_xx * pfx3 )
  fx(4) = -( stress_xx * pfx4 )
  fx(5) = -( stress_xx * pfx5 )
  fx(6) = -( stress_xx * pfx6 )
  fx(7) = -( stress_xx * pfx7 )
  
  fy(0) = -( stress_yy * pfy0  )
  fy(1) = -( stress_yy * pfy1  )
  fy(2) = -( stress_yy * pfy2  )
  fy(3) = -( stress_yy * pfy3  )
  fy(4) = -( stress_yy * pfy4  )
  fy(5) = -( stress_yy * pfy5  )
  fy(6) = -( stress_yy * pfy6  )
  fy(7) = -( stress_yy * pfy7  )
  
  fz(0) = -( stress_zz * pfz0 )
  fz(1) = -( stress_zz * pfz1 )
  fz(2) = -( stress_zz * pfz2 )
  fz(3) = -( stress_zz * pfz3 )
  fz(4) = -( stress_zz * pfz4 )
  fz(5) = -( stress_zz * pfz5 )
  fz(6) = -( stress_zz * pfz6 )
  fz(7) = -( stress_zz * pfz7 )

END SUBROUTINE SumElemStressesToNodeForces






SUBROUTINE IntegrateStressForElems(numElem, sigxx, sigyy, sigzz, determ)
  IMPLICIT NONE

  INTEGER      :: numElem
  REAL(KIND=8),DIMENSION(0:) :: sigxx, sigyy, sigzz
  REAL(KIND=8),DIMENSION(0:), INTENT(INOUT) :: determ  

  REAL(KIND=8) :: fx, fy, fz
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: fx_elem, fy_elem, fz_elem
  REAL(KIND=8),DIMENSION(0:7,0:2) :: B   ! shape function derivatives
  REAL(KIND=8),DIMENSION(0:7)   :: x_local
  REAL(KIND=8),DIMENSION(0:7)   :: y_local
  REAL(KIND=8),DIMENSION(0:7)   :: z_local
  INTEGER(KIND=4), DIMENSION(:), POINTER :: elemNodes => NULL()
  INTEGER      :: lnode, gnode, count, start, elem, kk
  INTEGER      :: numNode, numElem8

  numElem8 = numElem * 8
  ALLOCATE(fx_elem(0:numElem8-1))
  ALLOCATE(fy_elem(0:numElem8-1))
  ALLOCATE(fz_elem(0:numElem8-1))
  
! loop over all elements
  DO kk=0, numElem-1
    elemNodes => domain%m_nodelist(kk*8:)

!   get nodal coordinates from global arrays and copy into local arrays.
    DO lnode=0, 7
      gnode = elemNodes(lnode+1)
      x_local(lnode) = domain%m_x(gnode)
      y_local(lnode) = domain%m_y(gnode)
      z_local(lnode) = domain%m_z(gnode)
    ENDDO

!   Volume calculation involves extra work for numerical consistency.
    CALL CalcElemShapeFunctionDerivatives(x_local, y_local, z_local, &
                                          B, determ(kk))

    CALL CalcElemNodeNormals( B(:,0) , B(:,1), B(:,2), x_local, y_local, z_local )

    CALL SumElemStressesToNodeForces( B, sigxx(kk), sigyy(kk), sigzz(kk),  &
                                      fx_elem(kk*8), fy_elem(kk*8), fz_elem(kk*8) )

#if 0
!   copy nodal force contributions to global force arrray.
    DO lnode=0, 7
      node = elemNodes(lnode+1)
      domain%m_fx(gnode) = domain%m_fx(gnode) + fx_local(lnode)
      domain%m_fy(gnode) = domain%m_fy(gnode) + fy_local(lnode)
      domain%m_fz(gnode) = domain%m_fz(gnode) + fz_local(lnode)
    ENDDO
#endif
  ENDDO

  numNode = domain%m_numNode

  DO gnode=0, numNode-1
    count = domain%m_nodeElemCount(gnode)
    start = domain%m_nodeElemStart(gnode)
    fx = (0.0_RLK)
    fy = (0.0_RLK)
    fz = (0.0_RLK)
    DO i=0, count-1
      elem = domain%m_nodeElemCornerList(start+i)
      fx = fx + fx_elem(elem)
      fy = fy + fy_elem(elem)
      fz = fz + fz_elem(elem)
    ENDDO
    domain%m_fx(gnode) = fx
    domain%m_fy(gnode) = fy
    domain%m_fz(gnode) = fz
  ENDDO

  DEALLOCATE(fz_elem)
  DEALLOCATE(fy_elem)
  DEALLOCATE(fx_elem)

END SUBROUTINE IntegrateStressForElems



SUBROUTINE CollectDomainNodesToElemNodes(elemToNode, elemX, elemY, elemZ)

  IMPLICIT NONE 

  INTEGER, DIMENSION(:), POINTER :: elemToNode
  REAL(KIND=8),DIMENSION(0:7)    :: elemX, elemY, elemZ

  INTEGER(KIND=4) :: nd0i, nd1i, nd2i, nd3i
  INTEGER(KIND=4) :: nd4i, nd5i, nd6i, nd7i

  nd0i = elemToNode(1)
  nd1i = elemToNode(2)
  nd2i = elemToNode(3)
  nd3i = elemToNode(4)
  nd4i = elemToNode(5)
  nd5i = elemToNode(6)
  nd6i = elemToNode(7)
  nd7i = elemToNode(8)

  elemX(0) = domain%m_x(nd0i)
  elemX(1) = domain%m_x(nd1i)
  elemX(2) = domain%m_x(nd2i)
  elemX(3) = domain%m_x(nd3i)
  elemX(4) = domain%m_x(nd4i)
  elemX(5) = domain%m_x(nd5i)
  elemX(6) = domain%m_x(nd6i)
  elemX(7) = domain%m_x(nd7i)

  elemY(0) = domain%m_y(nd0i)
  elemY(1) = domain%m_y(nd1i)
  elemY(2) = domain%m_y(nd2i)
  elemY(3) = domain%m_y(nd3i)
  elemY(4) = domain%m_y(nd4i)
  elemY(5) = domain%m_y(nd5i)
  elemY(6) = domain%m_y(nd6i)
  elemY(7) = domain%m_y(nd7i)

  elemZ(0) = domain%m_z(nd0i)
  elemZ(1) = domain%m_z(nd1i)
  elemZ(2) = domain%m_z(nd2i)
  elemZ(3) = domain%m_z(nd3i)
  elemZ(4) = domain%m_z(nd4i)
  elemZ(5) = domain%m_z(nd5i)
  elemZ(6) = domain%m_z(nd6i)
  elemZ(7) = domain%m_z(nd7i)

END SUBROUTINE CollectDomainNodesToElemNodes









SUBROUTINE VoluDer(x0, x1, x2,      &
                   x3, x4, x5,      &
                   y0, y1, y2,      &
                   y3, y4, y5,      &
                   z0, z1, z2,      &
                   z3, z4, z5,      &
                   dvdx, dvdy, dvdz )
  IMPLICIT NONE
  REAL(KIND=8) :: x0, x1, x2, x3, x4, x5
  REAL(KIND=8) :: y0, y1, y2, y3, y4, y5
  REAL(KIND=8) :: z0, z1, z2, z3, z4, z5
  REAL(KIND=8) :: dvdx, dvdy, dvdz

  REAL(KIND=8), PARAMETER :: twelfth = 1.0_RLK / 12.0_RLK

  dvdx =                                              &
    (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +   &
    (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -   &
    (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5)

  dvdy =                                              &
    - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) - &
    (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +   &
    (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5)

  dvdz =                                              &
    - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) - &
    (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +   &
    (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5)

  dvdx = dvdx * twelfth
  dvdy = dvdy * twelfth
  dvdz = dvdz * twelfth

END SUBROUTINE VoluDer



SUBROUTINE CalcElemVolumeDerivative(dvdx,dvdy,dvdz, x, y, z)


  IMPLICIT NONE
  REAL(KIND=8),DIMENSION(0:7) :: dvdx, dvdy, dvdz
  REAL(KIND=8),DIMENSION(0:7) :: x, y, z

  CALL VoluDer(x(1), x(2), x(3), x(4), x(5), x(7),  &
               y(1), y(2), y(3), y(4), y(5), y(7),  &
               z(1), z(2), z(3), z(4), z(5), z(7),  &
               dvdx(0), dvdy(0), dvdz(0))
  CALL VoluDer(x(0), x(1), x(2), x(7), x(4), x(6),  &
               y(0), y(1), y(2), y(7), y(4), y(6),  &
               z(0), z(1), z(2), z(7), z(4), z(6),  &
               dvdx(3), dvdy(3), dvdz(3))
  CALL VoluDer(x(3), x(0), x(1), x(6), x(7), x(5),  &
               y(3), y(0), y(1), y(6), y(7), y(5),  &
               z(3), z(0), z(1), z(6), z(7), z(5),  &
               dvdx(2), dvdy(2), dvdz(2))
  CALL VoluDer(x(2), x(3), x(0), x(5), x(6), x(4),  &
               y(2), y(3), y(0), y(5), y(6), y(4),  &
               z(2), z(3), z(0), z(5), z(6), z(4),  &
               dvdx(1), dvdy(1), dvdz(1))
  CALL VoluDer(x(7), x(6), x(5), x(0), x(3), x(1),  &
               y(7), y(6), y(5), y(0), y(3), y(1),  &
               z(7), z(6), z(5), z(0), z(3), z(1),  &
               dvdx(4), dvdy(4), dvdz(4))
  CALL VoluDer(x(4), x(7), x(6), x(1), x(0), x(2),  &
               y(4), y(7), y(6), y(1), y(0), y(2),  &
               z(4), z(7), z(6), z(1), z(0), z(2),  &
               dvdx(5), dvdy(5), dvdz(5))
  CALL VoluDer(x(5), x(4), x(7), x(2), x(1), x(3),  &
               y(5), y(4), y(7), y(2), y(1), y(3),  &
               z(5), z(4), z(7), z(2), z(1), z(3),  &
               dvdx(6), dvdy(6), dvdz(6))
  CALL VoluDer(x(6), x(5), x(4), x(3), x(2), x(0),  &
               y(6), y(5), y(4), y(3), y(2), y(0),  &
               z(6), z(5), z(4), z(3), z(2), z(0),  &
               dvdx(7), dvdy(7), dvdz(7))

END SUBROUTINE CalcElemVolumeDerivative



SUBROUTINE CalcElemFBHourglassForce(xd, yd, zd, &
                                    hourgam0, hourgam1, &
                                    hourgam2, hourgam3, &
                                    hourgam4, hourgam5, &
                                    hourgam6, hourgam7, &
                                    coefficient, hgfx,  &
                                    hgfy, hgfz          )
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:7) :: xd,yd,zd
  REAL(KIND=8), DIMENSION(0:3) :: hourgam0, hourgam1,  &
                                  hourgam2, hourgam3,  &
                                  hourgam4, hourgam5,  &
                                  hourgam6, hourgam7
  REAL(KIND=8) :: coefficient
  REAL(KIND=8),DIMENSION(0:7) :: hgfx,hgfy,hgfz
  REAL(KIND=8) :: h00,h01,h02,h03

  INTEGER(KIND=4),PARAMETER :: i00 = 0_4
  INTEGER(KIND=4),PARAMETER :: i01 = 1_4
  INTEGER(KIND=4),PARAMETER :: i02 = 2_4
  INTEGER(KIND=4),PARAMETER :: i03 = 3_4

  h00 =                                             &
    hourgam0(i00) * xd(0) + hourgam1(i00) * xd(1) + &
    hourgam2(i00) * xd(2) + hourgam3(i00) * xd(3) + &
    hourgam4(i00) * xd(4) + hourgam5(i00) * xd(5) + &
    hourgam6(i00) * xd(6) + hourgam7(i00) * xd(7)
  
  h01 =                                             &
    hourgam0(i01) * xd(0) + hourgam1(i01) * xd(1) + &
    hourgam2(i01) * xd(2) + hourgam3(i01) * xd(3) + &
    hourgam4(i01) * xd(4) + hourgam5(i01) * xd(5) + &
    hourgam6(i01) * xd(6) + hourgam7(i01) * xd(7)
  
  h02 =                                             &
    hourgam0(i02) * xd(0) + hourgam1(i02) * xd(1) + &
    hourgam2(i02) * xd(2) + hourgam3(i02) * xd(3) + &
    hourgam4(i02) * xd(4) + hourgam5(i02) * xd(5) + &
    hourgam6(i02) * xd(6) + hourgam7(i02) * xd(7)
  
  h03 =                                             &
    hourgam0(i03) * xd(0) + hourgam1(i03) * xd(1) + &
    hourgam2(i03) * xd(2) + hourgam3(i03) * xd(3) + &
    hourgam4(i03) * xd(4) + hourgam5(i03) * xd(5) + &
    hourgam6(i03) * xd(6) + hourgam7(i03) * xd(7)
  
  hgfx(0) = coefficient *                       &
   (hourgam0(i00) * h00 + hourgam0(i01) * h01 + &
    hourgam0(i02) * h02 + hourgam0(i03) * h03)
  
  hgfx(1) = coefficient *                       &
   (hourgam1(i00) * h00 + hourgam1(i01) * h01 + &
    hourgam1(i02) * h02 + hourgam1(i03) * h03)
  
  hgfx(2) = coefficient *                       &
   (hourgam2(i00) * h00 + hourgam2(i01) * h01 + &
    hourgam2(i02) * h02 + hourgam2(i03) * h03)
  
  hgfx(3) = coefficient *                       &
   (hourgam3(i00) * h00 + hourgam3(i01) * h01 + &
    hourgam3(i02) * h02 + hourgam3(i03) * h03)
  
  hgfx(4) = coefficient *                       &
   (hourgam4(i00) * h00 + hourgam4(i01) * h01 + &
    hourgam4(i02) * h02 + hourgam4(i03) * h03)
  
  hgfx(5) = coefficient *                       &
   (hourgam5(i00) * h00 + hourgam5(i01) * h01 + &
    hourgam5(i02) * h02 + hourgam5(i03) * h03)
  
  hgfx(6) = coefficient *                       &
   (hourgam6(i00) * h00 + hourgam6(i01) * h01 + &
    hourgam6(i02) * h02 + hourgam6(i03) * h03)
  
  hgfx(7) = coefficient *                       &
   (hourgam7(i00) * h00 + hourgam7(i01) * h01 + &
    hourgam7(i02) * h02 + hourgam7(i03) * h03)
  
  h00 =                                             &
    hourgam0(i00) * yd(0) + hourgam1(i00) * yd(1) + &
    hourgam2(i00) * yd(2) + hourgam3(i00) * yd(3) + &
    hourgam4(i00) * yd(4) + hourgam5(i00) * yd(5) + &
    hourgam6(i00) * yd(6) + hourgam7(i00) * yd(7)
  
  h01 =                                             &
    hourgam0(i01) * yd(0) + hourgam1(i01) * yd(1) + &
    hourgam2(i01) * yd(2) + hourgam3(i01) * yd(3) + &
    hourgam4(i01) * yd(4) + hourgam5(i01) * yd(5) + &
    hourgam6(i01) * yd(6) + hourgam7(i01) * yd(7)
  
  h02 =                                            &
    hourgam0(i02) * yd(0) + hourgam1(i02) * yd(1)+ &
    hourgam2(i02) * yd(2) + hourgam3(i02) * yd(3)+ &
    hourgam4(i02) * yd(4) + hourgam5(i02) * yd(5)+ &
    hourgam6(i02) * yd(6) + hourgam7(i02) * yd(7)
  
  h03 =                                             &
    hourgam0(i03) * yd(0) + hourgam1(i03) * yd(1) + &
    hourgam2(i03) * yd(2) + hourgam3(i03) * yd(3) + &
    hourgam4(i03) * yd(4) + hourgam5(i03) * yd(5) + &
    hourgam6(i03) * yd(6) + hourgam7(i03) * yd(7)
  
  
  hgfy(0) = coefficient *                       &
   (hourgam0(i00) * h00 + hourgam0(i01) * h01 + &
    hourgam0(i02) * h02 + hourgam0(i03) * h03)
  
  hgfy(1) = coefficient *                       &
   (hourgam1(i00) * h00 + hourgam1(i01) * h01 + &
    hourgam1(i02) * h02 + hourgam1(i03) * h03)
  
  hgfy(2) = coefficient *                       &
   (hourgam2(i00) * h00 + hourgam2(i01) * h01 + &
    hourgam2(i02) * h02 + hourgam2(i03) * h03)
  
  hgfy(3) = coefficient *                       &
   (hourgam3(i00) * h00 + hourgam3(i01) * h01 + &
    hourgam3(i02) * h02 + hourgam3(i03) * h03)
  
  hgfy(4) = coefficient *                       &
   (hourgam4(i00) * h00 + hourgam4(i01) * h01 + &
    hourgam4(i02) * h02 + hourgam4(i03) * h03)
  
  hgfy(5) = coefficient *                       &
   (hourgam5(i00) * h00 + hourgam5(i01) * h01 + &
    hourgam5(i02) * h02 + hourgam5(i03) * h03)
  
  hgfy(6) = coefficient *                       &
   (hourgam6(i00) * h00 + hourgam6(i01) * h01 + &
    hourgam6(i02) * h02 + hourgam6(i03) * h03)
  
  hgfy(7) = coefficient *                       &
   (hourgam7(i00) * h00 + hourgam7(i01) * h01 + &
    hourgam7(i02) * h02 + hourgam7(i03) * h03)
  
  h00 =                                              &
    hourgam0(i00) * zd(0) + hourgam1(i00) * zd(1) +  &
    hourgam2(i00) * zd(2) + hourgam3(i00) * zd(3) +  &
    hourgam4(i00) * zd(4) + hourgam5(i00) * zd(5) +  &
    hourgam6(i00) * zd(6) + hourgam7(i00) * zd(7)
  
  h01 =                                              &
    hourgam0(i01) * zd(0) + hourgam1(i01) * zd(1) +  &
    hourgam2(i01) * zd(2) + hourgam3(i01) * zd(3) +  &
    hourgam4(i01) * zd(4) + hourgam5(i01) * zd(5) +  &
    hourgam6(i01) * zd(6) + hourgam7(i01) * zd(7)
  
  h02 =                                              &
    hourgam0(i02) * zd(0) + hourgam1(i02) * zd(1)+   &
    hourgam2(i02) * zd(2) + hourgam3(i02) * zd(3)+   &
    hourgam4(i02) * zd(4) + hourgam5(i02) * zd(5)+   &
    hourgam6(i02) * zd(6) + hourgam7(i02) * zd(7)
  
  h03 =                                              &
    hourgam0(i03) * zd(0) + hourgam1(i03) * zd(1) +  &
    hourgam2(i03) * zd(2) + hourgam3(i03) * zd(3) +  &
    hourgam4(i03) * zd(4) + hourgam5(i03) * zd(5) +  &
    hourgam6(i03) * zd(6) + hourgam7(i03) * zd(7)
  
  
  hgfz(0) = coefficient *                        &
   (hourgam0(i00) * h00 + hourgam0(i01) * h01 +  &
    hourgam0(i02) * h02 + hourgam0(i03) * h03)
  
  hgfz(1) = coefficient *                        &
   (hourgam1(i00) * h00 + hourgam1(i01) * h01 +  &
    hourgam1(i02) * h02 + hourgam1(i03) * h03)
  
  hgfz(2) = coefficient *                        &
   (hourgam2(i00) * h00 + hourgam2(i01) * h01 +  &
    hourgam2(i02) * h02 + hourgam2(i03) * h03)
  
  hgfz(3) = coefficient *                        &
   (hourgam3(i00) * h00 + hourgam3(i01) * h01 +  &
    hourgam3(i02) * h02 + hourgam3(i03) * h03)
  
  hgfz(4) = coefficient *                        &
   (hourgam4(i00) * h00 + hourgam4(i01) * h01 +  &
    hourgam4(i02) * h02 + hourgam4(i03) * h03)
  
  hgfz(5) = coefficient *                        &
   (hourgam5(i00) * h00 + hourgam5(i01) * h01 +  &
    hourgam5(i02) * h02 + hourgam5(i03) * h03)
  
  hgfz(6) = coefficient *                        &
   (hourgam6(i00) * h00 + hourgam6(i01) * h01 +  &
    hourgam6(i02) * h02 + hourgam6(i03) * h03)
  
  hgfz(7) = coefficient *                        &
   (hourgam7(i00) * h00 + hourgam7(i01) * h01 +  &
    hourgam7(i02) * h02 + hourgam7(i03) * h03)

END SUBROUTINE CalcElemFBHourglassForce



REAL(KIND=8) FUNCTION CBRT(dat)

  IMPLICIT NONE
  REAL(KIND=8) :: dat

  CBRT = dat**(1.0_RLK/3.0_RLK)

END FUNCTION CBRT




SUBROUTINE CalcFBHourglassForceForElems(determ,           &
                                        x8n, y8n, z8n,    &
                                        dvdx, dvdy, dvdz, &
                                        hourg             )

! *************************************************
! *
! *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
! *               force.
! *
! *************************************************


  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:) :: determ
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x8n, y8n, z8n
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dvdx, dvdy, dvdz
  REAL(KIND=8) :: hourg

  REAL(KIND=8) :: coefficient, volinv, ss1, mass1, volume13
  REAL(KIND=8) :: hourmodx, hourmody, hourmodz
  REAL(KIND=8) :: fx, fy, fz
  REAL(KIND=8), DIMENSION(0:7) :: hgfx, hgfy, hgfz
  REAL(KIND=8), DIMENSION(0:7) :: xd1, yd1, zd1
  REAL(KIND=8), DIMENSION(0:3) :: hourgam0, hourgam1, hourgam2, hourgam3
  REAL(KIND=8), DIMENSION(0:3) :: hourgam4, hourgam5, hourgam6, hourgam7
  REAL(KIND=8), DIMENSION(0:7,0:3) :: gamma
  REAL(KIND=8), DIMENSION(:), POINTER :: fx_elem, fy_elem, fz_elem
  REAL(KIND=8), DIMENSION(:), POINTER :: fx_local, fy_local, fz_local
  INTEGER(KIND=4) :: numElem, numElem8, i, i2, i3, i1
  INTEGER(KIND=4) :: numNode, gnode, elem, count, start
  INTEGER(KIND=4) :: n0si2, n1si2, n2si2, n3si2, n4si2, n5si2, n6si2, n7si2
  INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()

  NULLIFY(fx_local, fy_local, fz_local)
  NULLIFY(fx_elem, fy_elem, fz_elem)
  numElem = domain%m_numElem
  numElem8 = numElem * 8
  ALLOCATE(fx_elem(0:numElem8-1))
  ALLOCATE(fy_elem(0:numElem8-1))
  ALLOCATE(fz_elem(0:numElem8-1))

  gamma(0,0) = ( 1.0_RLK)
  gamma(1,0) = ( 1.0_RLK)
  gamma(2,0) = (-1.0_RLK)
  gamma(3,0) = (-1.0_RLK)
  gamma(4,0) = (-1.0_RLK)
  gamma(5,0) = (-1.0_RLK)
  gamma(6,0) = ( 1.0_RLK)
  gamma(7,0) = ( 1.0_RLK)
  gamma(0,1) = ( 1.0_RLK)
  gamma(1,1) = (-1.0_RLK)
  gamma(2,1) = (-1.0_RLK)
  gamma(3,1) = ( 1.0_RLK)
  gamma(4,1) = (-1.0_RLK)
  gamma(5,1) = ( 1.0_RLK)
  gamma(6,1) = ( 1.0_RLK)
  gamma(7,1) = (-1.0_RLK)
  gamma(0,2) = ( 1.0_RLK)
  gamma(1,2) = (-1.0_RLK)
  gamma(2,2) = ( 1.0_RLK)
  gamma(3,2) = (-1.0_RLK)
  gamma(4,2) = ( 1.0_RLK)
  gamma(5,2) = (-1.0_RLK)
  gamma(6,2) = ( 1.0_RLK)
  gamma(7,2) = (-1.0_RLK)
  gamma(0,3) = (-1.0_RLK)
  gamma(1,3) = ( 1.0_RLK)
  gamma(2,3) = (-1.0_RLK)
  gamma(3,3) = ( 1.0_RLK)
  gamma(4,3) = ( 1.0_RLK)
  gamma(5,3) = (-1.0_RLK)
  gamma(6,3) = ( 1.0_RLK)
  gamma(7,3) = (-1.0_RLK)
  
! *************************************************
! compute the hourglass modes
  
  
  DO i2=0, numElem-1

    elemToNode => domain%m_nodelist(i2*8:)

    i3=8*i2
    volinv= (1.0_RLK)/determ(i2)

    DO i1=0, 3

      hourmodx =                                             &
        x8n(i3)   * gamma(0,i1) + x8n(i3+1) * gamma(1,i1) +  &
        x8n(i3+2) * gamma(2,i1) + x8n(i3+3) * gamma(3,i1) +  &
        x8n(i3+4) * gamma(4,i1) + x8n(i3+5) * gamma(5,i1) +  &
        x8n(i3+6) * gamma(6,i1) + x8n(i3+7) * gamma(7,i1)

      hourmody =                                             &
        y8n(i3)   * gamma(0,i1) + y8n(i3+1) * gamma(1,i1) +  &
        y8n(i3+2) * gamma(2,i1) + y8n(i3+3) * gamma(3,i1) +  &
        y8n(i3+4) * gamma(4,i1) + y8n(i3+5) * gamma(5,i1) +  &
        y8n(i3+6) * gamma(6,i1) + y8n(i3+7) * gamma(7,i1)

      hourmodz =                                             &
        z8n(i3)   * gamma(0,i1) + z8n(i3+1) * gamma(1,i1) +  &
        z8n(i3+2) * gamma(2,i1) + z8n(i3+3) * gamma(3,i1) +  &
        z8n(i3+4) * gamma(4,i1) + z8n(i3+5) * gamma(5,i1) +  &
        z8n(i3+6) * gamma(6,i1) + z8n(i3+7) * gamma(7,i1)

      hourgam0(i1) = gamma(0,i1) -  volinv*(dvdx(i3  ) * hourmodx +  &
                    dvdy(i3  ) * hourmody + dvdz(i3  ) * hourmodz )

      hourgam1(i1) = gamma(1,i1) -  volinv*(dvdx(i3+1) * hourmodx +  &
                    dvdy(i3+1) * hourmody + dvdz(i3+1) * hourmodz )

      hourgam2(i1) = gamma(2,i1) -  volinv*(dvdx(i3+2) * hourmodx +  &
                    dvdy(i3+2) * hourmody + dvdz(i3+2) * hourmodz )

      hourgam3(i1) = gamma(3,i1) -  volinv*(dvdx(i3+3) * hourmodx +  &
                    dvdy(i3+3) * hourmody + dvdz(i3+3) * hourmodz )

      hourgam4(i1) = gamma(4,i1) -  volinv*(dvdx(i3+4) * hourmodx +  &
                    dvdy(i3+4) * hourmody + dvdz(i3+4) * hourmodz )

      hourgam5(i1) = gamma(5,i1) -  volinv*(dvdx(i3+5) * hourmodx +  &
                    dvdy(i3+5) * hourmody + dvdz(i3+5) * hourmodz )

      hourgam6(i1) = gamma(6,i1) -  volinv*(dvdx(i3+6) * hourmodx +  &
                    dvdy(i3+6) * hourmody + dvdz(i3+6) * hourmodz )

      hourgam7(i1) = gamma(7,i1) -  volinv*(dvdx(i3+7) * hourmodx +  &
                    dvdy(i3+7) * hourmody + dvdz(i3+7) * hourmodz )

    ENDDO

!   compute forces
!   store forces into h arrays (force arrays)

    ss1=domain%m_ss(i2)
    mass1=domain%m_elemMass(i2)
    volume13=CBRT(determ(i2))

    n0si2 = elemToNode(1)
    n1si2 = elemToNode(2)
    n2si2 = elemToNode(3)
    n3si2 = elemToNode(4)
    n4si2 = elemToNode(5)
    n5si2 = elemToNode(6)
    n6si2 = elemToNode(7)
    n7si2 = elemToNode(8)

    xd1(0) = domain%m_xd(n0si2)
    xd1(1) = domain%m_xd(n1si2)
    xd1(2) = domain%m_xd(n2si2)
    xd1(3) = domain%m_xd(n3si2)
    xd1(4) = domain%m_xd(n4si2)
    xd1(5) = domain%m_xd(n5si2)
    xd1(6) = domain%m_xd(n6si2)
    xd1(7) = domain%m_xd(n7si2)

    yd1(0) = domain%m_yd(n0si2)
    yd1(1) = domain%m_yd(n1si2)
    yd1(2) = domain%m_yd(n2si2)
    yd1(3) = domain%m_yd(n3si2)
    yd1(4) = domain%m_yd(n4si2)
    yd1(5) = domain%m_yd(n5si2)
    yd1(6) = domain%m_yd(n6si2)
    yd1(7) = domain%m_yd(n7si2)

    zd1(0) = domain%m_zd(n0si2)
    zd1(1) = domain%m_zd(n1si2)
    zd1(2) = domain%m_zd(n2si2)
    zd1(3) = domain%m_zd(n3si2)
    zd1(4) = domain%m_zd(n4si2)
    zd1(5) = domain%m_zd(n5si2)
    zd1(6) = domain%m_zd(n6si2)
    zd1(7) = domain%m_zd(n7si2)

    coefficient = - hourg * (0.01_RLK) * ss1 * mass1 / volume13

    CALL CalcElemFBHourglassForce(xd1,yd1,zd1,                          &
                                  hourgam0,hourgam1,hourgam2,hourgam3,  &
                                  hourgam4,hourgam5,hourgam6,hourgam7,  &
                                  coefficient, hgfx, hgfy, hgfz)

    fx_local(0:) => fx_elem(i3:)
    fx_local(0) = hgfx(0)
    fx_local(1) = hgfx(1)
    fx_local(2) = hgfx(2)
    fx_local(3) = hgfx(3)
    fx_local(4) = hgfx(4)
    fx_local(5) = hgfx(5)
    fx_local(6) = hgfx(6)
    fx_local(7) = hgfx(7)

    fy_local(0:) => fy_elem(i3:)
    fy_local(0) = hgfy(0)
    fy_local(1) = hgfy(1)
    fy_local(2) = hgfy(2)
    fy_local(3) = hgfy(3)
    fy_local(4) = hgfy(4)
    fy_local(5) = hgfy(5)
    fy_local(6) = hgfy(6)
    fy_local(7) = hgfy(7)

    fz_local(0:) => fz_elem(i3:)
    fz_local(0) = hgfz(0)
    fz_local(1) = hgfz(1)
    fz_local(2) = hgfz(2)
    fz_local(3) = hgfz(3)
    fz_local(4) = hgfz(4)
    fz_local(5) = hgfz(5)
    fz_local(6) = hgfz(6)
    fz_local(7) = hgfz(7)

#if 0
    domain%m_fx(n0si2) = domain%m_fx(n0si2) + hgfx(0)
    domain%m_fy(n0si2) = domain%m_fy(n0si2) + hgfy(0)
    domain%m_fz(n0si2) = domain%m_fz(n0si2) + hgfz(0)

    domain%m_fx(n1si2) = domain%m_fx(n1si2) + hgfx(1)
    domain%m_fy(n1si2) = domain%m_fy(n1si2) + hgfy(1)
    domain%m_fz(n1si2) = domain%m_fz(n1si2) + hgfz(1)

    domain%m_fx(n2si2) = domain%m_fx(n2si2) + hgfx(2)
    domain%m_fy(n2si2) = domain%m_fy(n2si2) + hgfy(2)
    domain%m_fz(n2si2) = domain%m_fz(n2si2) + hgfz(2)

    domain%m_fx(n3si2) = domain%m_fx(n3si2) + hgfx(3)
    domain%m_fy(n3si2) = domain%m_fy(n3si2) + hgfy(3)
    domain%m_fz(n3si2) = domain%m_fz(n3si2) + hgfz(3)

    domain%m_fx(n4si2) = domain%m_fx(n4si2) + hgfx(4)
    domain%m_fy(n4si2) = domain%m_fy(n4si2) + hgfy(4)
    domain%m_fz(n4si2) = domain%m_fz(n4si2) + hgfz(4)

    domain%m_fx(n5si2) = domain%m_fx(n5si2) + hgfx(5)
    domain%m_fy(n5si2) = domain%m_fy(n5si2) + hgfy(5)
    domain%m_fz(n5si2) = domain%m_fz(n5si2) + hgfz(5)

    domain%m_fx(n6si2) = domain%m_fx(n6si2) + hgfx(6)
    domain%m_fy(n6si2) = domain%m_fy(n6si2) + hgfy(6)
    domain%m_fz(n6si2) = domain%m_fz(n6si2) + hgfz(6)

    domain%m_fx(n7si2) = domain%m_fx(n7si2) + hgfx(7)
    domain%m_fy(n7si2) = domain%m_fy(n7si2) + hgfy(7)
    domain%m_fz(n7si2) = domain%m_fz(n7si2) + hgfz(7)
#endif
  ENDDO

  numNode = domain%m_numNode

  DO gnode=0, numNode-1

    count = domain%m_nodeElemCount(gnode)
    start = domain%m_nodeElemStart(gnode)
    fx = (0.0_RLK)
    fy = (0.0_RLK)
    fz = (0.0_RLK)
    DO i = 0, count-1
      elem = domain%m_nodeElemCornerList(start+i)
      fx = fx + fx_elem(elem)
      fy = fy + fy_elem(elem)
      fz = fz + fz_elem(elem)
    ENDDO
    domain%m_fx(gnode) = domain%m_fx(gnode) + fx
    domain%m_fy(gnode) = domain%m_fy(gnode) + fy
    domain%m_fz(gnode) = domain%m_fz(gnode) + fz
  ENDDO

  DEALLOCATE(fz_elem)
  DEALLOCATE(fy_elem)
  DEALLOCATE(fx_elem)


END SUBROUTINE CalcFBHourglassForceForElems





SUBROUTINE CalcHourglassControlForElems(determ, hgcoef)

  IMPLICIT NONE
  REAL(KIND=8),DIMENSION(0:) :: determ
  REAL(KIND=8) :: hgcoef

  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: dvdx
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: dvdy
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: dvdz
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: x8n
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: y8n
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: z8n
  REAL(KIND=8),DIMENSION(0:7) :: x1, y1, z1
  REAL(KIND=8),DIMENSION(0:7) :: pfx, pfy, pfz
  INTEGER(KIND=4) :: numElem, numElem8, i, ii, jj
  INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()

  numElem = domain%m_numElem
  numElem8 = numElem * 8
  ALLOCATE(dvdx(0:numElem8-1))
  ALLOCATE(dvdy(0:numElem8-1))
  ALLOCATE(dvdz(0:numElem8-1))
  ALLOCATE(x8n(0:numElem8-1))
  ALLOCATE(y8n(0:numElem8-1))
  ALLOCATE(z8n(0:numElem8-1))

! start loop over elements
  DO i=0, numElem-1
    elemToNode => domain%m_nodelist(i*8:)
    CALL CollectDomainNodesToElemNodes(elemToNode, x1, y1, z1)
    CALL CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1)

!   load into temporary storage for FB Hour Glass control
    DO ii=0, 7
      jj=8*i+ii

      dvdx(jj) = pfx(ii)
      dvdy(jj) = pfy(ii)
      dvdz(jj) = pfz(ii)
      x8n(jj)  = x1(ii)
      y8n(jj)  = y1(ii)
      z8n(jj)  = z1(ii)
    ENDDO

    determ(i) = domain%m_volo(i) * domain%m_v(i)

!   Do a check for negative volumes
    IF ( domain%m_v(i) <= (0.0_RLK) ) THEN
      CALL luabort(VolumeError)
    ENDIF
  ENDDO

  IF ( hgcoef > (0.0_RLK) ) THEN
    CALL CalcFBHourglassForceForElems(determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef)
  ENDIF

  DEALLOCATE(z8n)
  DEALLOCATE(y8n)
  DEALLOCATE(x8n)
  DEALLOCATE(dvdz)
  DEALLOCATE(dvdy)
  DEALLOCATE(dvdx)

  RETURN

END SUBROUTINE CalcHourglassControlForElems



SUBROUTINE CalcVolumeForceForElems()

  IMPLICIT NONE
  INTEGER(KIND=4) :: numElem
  INTEGER(KIND=4) :: k
  REAL(KIND=8) :: hgcoef
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sigxx
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sigyy
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sigzz
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: determ


  numElem = domain%m_numElem
  IF (numElem /= 0) THEN
    hgcoef = domain%m_hgcoef
    ALLOCATE(sigxx(0:numElem-1))
    ALLOCATE(sigyy(0:numElem-1))
    ALLOCATE(sigzz(0:numElem-1))
    ALLOCATE(determ(0:numElem-1))

!   Sum contributions to total stress tensor
    CALL InitStressTermsForElems(numElem, sigxx, sigyy, sigzz)

!   call elemlib stress integration loop to produce nodal forces from
!   material stresses.
    CALL IntegrateStressForElems( numElem, sigxx, sigyy, sigzz, determ)

!   check for negative element volume and abort if found
    DO k=0, numElem-1
       IF (determ(k) <= 0.0_RLK) THEN
         CALL luabort(VolumeError)
       ENDIF
    ENDDO

    CALL CalcHourglassControlForElems(determ, hgcoef)

    DEALLOCATE(determ)
    DEALLOCATE(sigzz)
    DEALLOCATE(sigyy)
    DEALLOCATE(sigxx)
  ENDIF

END SUBROUTINE CalcVolumeForceForElems


SUBROUTINE CalcForceForNodes()

  IMPLICIT NONE 
  INTEGER(KIND=4) :: numNode
  INTEGER(KIND=4) :: i

  numNode = domain%m_numNode
  DO i=0, numNode-1
    domain%m_fx(i) = 0.0_RLK
    domain%m_fy(i) = 0.0_RLK
    domain%m_fz(i) = 0.0_RLK
  ENDDO

! Calcforce calls partial, force, hourq
  CALL CalcVolumeForceForElems()

! Calculate Nodal Forces at domain boundaries
! problem->commSBN->Transfer(CommSBN::forces)

END SUBROUTINE CalcForceForNodes


SUBROUTINE CalcAccelerationForNodes()

  IMPLICIT NONE 
  INTEGER(KIND=4) :: numNode
  INTEGER(KIND=4) :: i

  numNode = domain%m_numNode
  DO i=0, numNode-1
    domain%m_xdd(i) = domain%m_fx(i) / domain%m_nodalMass(i)
    domain%m_ydd(i) = domain%m_fy(i) / domain%m_nodalMass(i)
    domain%m_zdd(i) = domain%m_fz(i) / domain%m_nodalMass(i)
  ENDDO

END SUBROUTINE CalcAccelerationForNodes


SUBROUTINE ApplyAccelerationBoundaryConditionsForNodes()

  IMPLICIT NONE 
  INTEGER(KIND=4) :: numNodeBC
  INTEGER(KIND=4) :: i

  numNodeBC = (domain%m_sizeX+1)*(domain%m_sizeX+1)

  DO i=0, numNodeBC-1
    domain%m_xdd(domain%m_symmX(i)) = 0.0_RLK
  ENDDO
  DO i=0, numNodeBC-1
    domain%m_ydd(domain%m_symmY(i)) = 0.0_RLK
  ENDDO
  DO i=0, numNodeBC-1
    domain%m_zdd(domain%m_symmZ(i)) = 0.0_RLK
  ENDDO

END SUBROUTINE ApplyAccelerationBoundaryConditionsForNodes




SUBROUTINE CalcVelocityForNodes( dt, u_cut)

  IMPLICIT NONE 
  REAL(KIND=8)    :: dt, u_cut
  INTEGER(KIND=4) :: numNode
  INTEGER(KIND=4) :: i
  REAL(KIND=8)    :: xdtmp, ydtmp, zdtmp

  numNode = domain%m_numNode

  DO i = 0, numNode-1

    xdtmp = domain%m_xd(i) + domain%m_xdd(i) * dt
    if( ABS(xdtmp) < u_cut ) xdtmp = 0.0_RLK
    domain%m_xd(i) = xdtmp

    ydtmp = domain%m_yd(i) + domain%m_ydd(i) * dt
    if( ABS(ydtmp) < u_cut ) ydtmp = 0.0_RLK
    domain%m_yd(i) = ydtmp

    zdtmp = domain%m_zd(i) + domain%m_zdd(i) * dt
    if( ABS(zdtmp) < u_cut ) zdtmp = 0.0_RLK
    domain%m_zd(i) = zdtmp
  ENDDO


END SUBROUTINE CalcVelocityForNodes






SUBROUTINE CalcPositionForNodes(dt)

  IMPLICIT NONE 
  REAL(KIND=8)    :: dt
  INTEGER(KIND=4) :: numNode
  INTEGER(KIND=4) :: i

  numNode = domain%m_numNode

  DO i = 0, numNode-1
    domain%m_x(i) = domain%m_x(i) + domain%m_xd(i) * dt
    domain%m_y(i) = domain%m_y(i) + domain%m_yd(i) * dt
    domain%m_z(i) = domain%m_z(i) + domain%m_zd(i) * dt
  ENDDO

END SUBROUTINE CalcPositionForNodes





SUBROUTINE LagrangeNodal()

  IMPLICIT NONE 
  REAL(KIND=8) :: delt
  REAL(KIND=8) :: u_cut

  delt  = domain%m_deltatime
  u_cut = domain%m_u_cut

! time of boundary condition evaluation is beginning of step for force and
! acceleration boundary conditions.
  CALL CalcForceForNodes()

  CALL CalcAccelerationForNodes()

  CALL ApplyAccelerationBoundaryConditionsForNodes()

  CALL CalcVelocityForNodes( delt, u_cut )

  CALL CalcPositionForNodes( delt )

END SUBROUTINE LagrangeNodal



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






FUNCTION CalcElemCharacteristicLength( x, y, z, volume) RESULT(charLength)

  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:7) :: x, y, z
  REAL(KIND=8) :: volume
  REAL(KIND=8) :: a
  REAL(KIND=8) :: charLength

  charLength = 0.0_RLK

  a = AreaFace(x(0),x(1),x(2),x(3),  &
               y(0),y(1),y(2),y(3),  &
               z(0),z(1),z(2),z(3))
  charLength = MAX(a,charLength)
  
  a = AreaFace(x(4),x(5),x(6),x(7),  &
               y(4),y(5),y(6),y(7),  &
               z(4),z(5),z(6),z(7))
  charLength = MAX(a,charLength)
  
  a = AreaFace(x(0),x(1),x(5),x(4),  &
               y(0),y(1),y(5),y(4),  &
               z(0),z(1),z(5),z(4))
  charLength = MAX(a,charLength)

  a = AreaFace(x(1),x(2),x(6),x(5),  &
               y(1),y(2),y(6),y(5),  &
               z(1),z(2),z(6),z(5))
  charLength = MAX(a,charLength)
  
  a = AreaFace(x(2),x(3),x(7),x(6),  &
               y(2),y(3),y(7),y(6),  &
               z(2),z(3),z(7),z(6))
  charLength = MAX(a,charLength)
  
  a = AreaFace(x(3),x(0),x(4),x(7),  &
               y(3),y(0),y(4),y(7),  &
               z(3),z(0),z(4),z(7))
  charLength = MAX(a,charLength)
  
  charLength = (4.0_RLK) * volume / SQRT(charLength);

  RETURN

END FUNCTION CalcElemCharacteristicLength










SUBROUTINE CalcElemVelocityGrandient( xvel, yvel, zvel, &
                                      b, detJ, d )

  IMPLICIT NONE 

  REAL(KIND=8), DIMENSION(0:7),     INTENT(IN)  :: xvel, yvel, zvel
  REAL(KIND=8), DIMENSION(0:7,0:2), INTENT(IN)  :: b   ![3,8]
  REAL(KIND=8),                     INTENT(IN)  :: detJ
  REAL(KIND=8), DIMENSION(0:5),     INTENT(OUT) :: d

  REAL(KIND=8) :: dyddx, dxddy, dzddx, dxddz, dzddy, dyddz
  REAL(KIND=8) :: inv_detJ
  REAL(KIND=8), DIMENSION(0:7) :: pfx
  REAL(KIND=8), DIMENSION(0:7) :: pfy
  REAL(KIND=8), DIMENSION(0:7) :: pfz

  inv_detJ = (1.0_RLK) / detJ
  pfx = b(:,0)
  pfy = b(:,1)
  pfz = b(:,2)

  d(0) = inv_detJ * ( pfx(0) * (xvel(0)-xvel(6))   &
                    + pfx(1) * (xvel(1)-xvel(7))   &
                    + pfx(2) * (xvel(2)-xvel(4))   &
                    + pfx(3) * (xvel(3)-xvel(5)) )
  
  d(1) = inv_detJ * ( pfy(0) * (yvel(0)-yvel(6))   &
                    + pfy(1) * (yvel(1)-yvel(7))   &
                    + pfy(2) * (yvel(2)-yvel(4))   &
                    + pfy(3) * (yvel(3)-yvel(5)) )
  
  d(2) = inv_detJ * ( pfz(0) * (zvel(0)-zvel(6))   &
                    + pfz(1) * (zvel(1)-zvel(7))   &
                    + pfz(2) * (zvel(2)-zvel(4))   &
                    + pfz(3) * (zvel(3)-zvel(5)) )

  dyddx = inv_detJ * ( pfx(0) * (yvel(0)-yvel(6))  &
                     + pfx(1) * (yvel(1)-yvel(7))  &
                     + pfx(2) * (yvel(2)-yvel(4))  &
                     + pfx(3) * (yvel(3)-yvel(5)) )
  
  dxddy = inv_detJ * ( pfy(0) * (xvel(0)-xvel(6))  &
                     + pfy(1) * (xvel(1)-xvel(7))  &
                     + pfy(2) * (xvel(2)-xvel(4))  &
                     + pfy(3) * (xvel(3)-xvel(5)) )
  
  dzddx = inv_detJ * ( pfx(0) * (zvel(0)-zvel(6))  &
                     + pfx(1) * (zvel(1)-zvel(7))  &
                     + pfx(2) * (zvel(2)-zvel(4))  &
                     + pfx(3) * (zvel(3)-zvel(5)) )
  
  dxddz = inv_detJ * ( pfz(0) * (xvel(0)-xvel(6))  &
                     + pfz(1) * (xvel(1)-xvel(7))  &
                     + pfz(2) * (xvel(2)-xvel(4))  &
                     + pfz(3) * (xvel(3)-xvel(5)) )
  
  dzddy = inv_detJ * ( pfy(0) * (zvel(0)-zvel(6))  &
                     + pfy(1) * (zvel(1)-zvel(7))  &
                     + pfy(2) * (zvel(2)-zvel(4))  &
                     + pfy(3) * (zvel(3)-zvel(5)) )
  
  dyddz = inv_detJ * ( pfz(0) * (yvel(0)-yvel(6))  &
                     + pfz(1) * (yvel(1)-yvel(7))  &
                     + pfz(2) * (yvel(2)-yvel(4))  &
                     + pfz(3) * (yvel(3)-yvel(5)) )

  d(5) = (0.5_RLK) * ( dxddy + dyddx )
  d(4) = (0.5_RLK) * ( dxddz + dzddx )
  d(3) = (0.5_RLK) * ( dzddy + dyddz )

END SUBROUTINE CalcElemVelocityGrandient








SUBROUTINE CalcKinematicsForElems( numElem, dt )

  IMPLICIT NONE
  INTEGER      :: numElem
  INTEGER      :: k, lnode, gnode, j
  REAL(KIND=8) :: dt

  REAL(KIND=8), DIMENSION(0:7,0:2) :: B  ! shape function derivatives
  REAL(KIND=8), DIMENSION(0:5):: D
  REAL(KIND=8), DIMENSION(0:7) :: x_local
  REAL(KIND=8), DIMENSION(0:7) :: y_local
  REAL(KIND=8), DIMENSION(0:7) :: z_local
  REAL(KIND=8), DIMENSION(0:7) :: xd_local
  REAL(KIND=8), DIMENSION(0:7) :: yd_local
  REAL(KIND=8), DIMENSION(0:7) :: zd_local
  REAL(KIND=8) :: detJ, volume, relativeVolume,dt2
  INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()

  detJ = 0.0_RLK

! loop over all elements
  DO k = 0, numElem-1
    elemToNode => domain%m_nodelist(k*8:)

!   get nodal coordinates from global arrays and copy into local arrays
    DO lnode=0, 7
      gnode = elemToNode(lnode+1)
      x_local(lnode) = domain%m_x(gnode)
      y_local(lnode) = domain%m_y(gnode)
      z_local(lnode) = domain%m_z(gnode)
    ENDDO

!   volume calculations
    volume = CalcElemVolume(x_local, y_local, z_local )
    relativeVolume = volume / domain%m_volo(k)
    domain%m_vnew(k) = relativeVolume
    domain%m_delv(k) = relativeVolume - domain%m_v(k)

!   set characteristic length
    domain%m_arealg(k) = CalcElemCharacteristicLength(x_local, y_local,  &
                                                      z_local, volume)

!   get nodal velocities from global array and copy into local arrays.
    DO lnode=0, 7
      gnode = elemToNode(lnode+1);
      xd_local(lnode) = domain%m_xd(gnode)
      yd_local(lnode) = domain%m_yd(gnode)
      zd_local(lnode) = domain%m_zd(gnode)
    ENDDO
    
    dt2 = (0.5_RLK) * dt
    DO j=0, 7
      x_local(j) = x_local(j) - dt2 * xd_local(j)
      y_local(j) = y_local(j) - dt2 * yd_local(j)
      z_local(j) = z_local(j) - dt2 * zd_local(j)
    ENDDO
    
    CALL CalcElemShapeFunctionDerivatives( x_local, y_local, z_local,  &
                                           B, detJ )
    
    CALL CalcElemVelocityGrandient( xd_local, yd_local, zd_local,  &
                                    B, detJ, D )
    
!   put velocity gradient quantities into their global arrays.
    domain%m_dxx(k) = D(0);
    domain%m_dyy(k) = D(1);
    domain%m_dzz(k) = D(2);
  ENDDO

END SUBROUTINE CalcKinematicsForElems



SUBROUTINE CalcLagrangeElements( deltatime)

  IMPLICIT NONE 
  REAL(KIND=8)    :: deltatime
  REAL(KIND=8)    :: vdov, vdovthird
  INTEGER(KIND=4) :: numElem, k

  numElem = domain%m_numElem
  IF (numElem > 0) THEN
    CALL CalcKinematicsForElems(numElem, deltatime)

!   element loop to do some stuff not included in the elemlib function.

    DO k=0, numElem-1
!     calc strain rate and apply as constraint (only done in FB element)
      vdov = domain%m_dxx(k) + domain%m_dyy(k) + domain%m_dzz(k)
      vdovthird = vdov/(3.0_RLK)

!     make the rate of deformation tensor deviatoric
      domain%m_vdov(k) = vdov
      domain%m_dxx(k) = domain%m_dxx(k) - vdovthird
      domain%m_dyy(k) = domain%m_dyy(k) - vdovthird
      domain%m_dzz(k) = domain%m_dzz(k) - vdovthird

!     See if any volumes are negative, and take appropriate action.
      IF (domain%m_vnew(k) <= (0.0_RLK)) THEN
        call luabort(VolumeError)
      ENDIF
    ENDDO
  ENDIF

END SUBROUTINE CalcLagrangeElements




REAL(KIND=8) FUNCTION SUM4(a, b, c, d)

  IMPLICIT NONE 
  REAL(KIND=8) :: a, b, c, d

  SUM4 = a + b + c + d

  RETURN

END FUNCTION SUM4




SUBROUTINE  CalcMonotonicQGradientsForElems()

  IMPLICIT NONE 
  REAL(KIND=8), PARAMETER :: ptiny = 1.e-36_RLK
  REAL(KIND=8)            :: ax,ay,az,dxv,dyv,dzv
  REAL(KIND=8)            :: x0,x1,x2,x3,x4,x5,x6,x7
  REAL(KIND=8)            :: y0,y1,y2,y3,y4,y5,y6,y7
  REAL(KIND=8)            :: z0,z1,z2,z3,z4,z5,z6,z7
  REAL(KIND=8)            :: xv0,xv1,xv2,xv3,xv4,xv5,xv6,xv7
  REAL(KIND=8)            :: yv0,yv1,yv2,yv3,yv4,yv5,yv6,yv7
  REAL(KIND=8)            :: zv0,zv1,zv2,zv3,zv4,zv5,zv6,zv7
  REAL(KIND=8)            :: vol, norm
  REAL(KIND=8)            :: dxi,dxj,dxk
  REAL(KIND=8)            :: dyi,dyj,dyk
  REAL(KIND=8)            :: dzi,dzj,dzk
  INTEGER(KIND=4)         :: numElem, i
  INTEGER(KIND=4)         :: n0,n1,n2,n3,n4,n5,n6,n7
  INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()

  numElem = domain%m_numElem

  DO i=0, numElem-1
    
    elemToNode => domain%m_nodelist(i*8:)
    n0 = elemToNode(1)
    n1 = elemToNode(2)
    n2 = elemToNode(3)
    n3 = elemToNode(4)
    n4 = elemToNode(5)
    n5 = elemToNode(6)
    n6 = elemToNode(7)
    n7 = elemToNode(8)

    x0 = domain%m_x(n0)
    x1 = domain%m_x(n1)
    x2 = domain%m_x(n2)
    x3 = domain%m_x(n3)
    x4 = domain%m_x(n4)
    x5 = domain%m_x(n5)
    x6 = domain%m_x(n6)
    x7 = domain%m_x(n7)
    
    y0 = domain%m_y(n0)
    y1 = domain%m_y(n1)
    y2 = domain%m_y(n2)
    y3 = domain%m_y(n3)
    y4 = domain%m_y(n4)
    y5 = domain%m_y(n5)
    y6 = domain%m_y(n6)
    y7 = domain%m_y(n7)
    
    z0 = domain%m_z(n0)
    z1 = domain%m_z(n1)
    z2 = domain%m_z(n2)
    z3 = domain%m_z(n3)
    z4 = domain%m_z(n4)
    z5 = domain%m_z(n5)
    z6 = domain%m_z(n6)
    z7 = domain%m_z(n7)
    
    xv0 = domain%m_xd(n0)
    xv1 = domain%m_xd(n1)
    xv2 = domain%m_xd(n2)
    xv3 = domain%m_xd(n3)
    xv4 = domain%m_xd(n4)
    xv5 = domain%m_xd(n5)
    xv6 = domain%m_xd(n6)
    xv7 = domain%m_xd(n7)
    
    yv0 = domain%m_yd(n0)
    yv1 = domain%m_yd(n1)
    yv2 = domain%m_yd(n2)
    yv3 = domain%m_yd(n3)
    yv4 = domain%m_yd(n4)
    yv5 = domain%m_yd(n5)
    yv6 = domain%m_yd(n6)
    yv7 = domain%m_yd(n7)
    
    zv0 = domain%m_zd(n0)
    zv1 = domain%m_zd(n1)
    zv2 = domain%m_zd(n2)
    zv3 = domain%m_zd(n3)
    zv4 = domain%m_zd(n4)
    zv5 = domain%m_zd(n5)
    zv6 = domain%m_zd(n6)
    zv7 = domain%m_zd(n7)
    
    vol = domain%m_volo(i)*domain%m_vnew(i)
    norm = (1.0_RLK) / ( vol + ptiny )
    
    dxj = (-0.25_RLK)*(SUM4(x0,x1,x5,x4) - SUM4(x3,x2,x6,x7))
    dyj = (-0.25_RLK)*(SUM4(y0,y1,y5,y4) - SUM4(y3,y2,y6,y7))
    dzj = (-0.25_RLK)*(SUM4(z0,z1,z5,z4) - SUM4(z3,z2,z6,z7))
    
    dxi = ( 0.25_RLK)*(SUM4(x1,x2,x6,x5) - SUM4(x0,x3,x7,x4))
    dyi = ( 0.25_RLK)*(SUM4(y1,y2,y6,y5) - SUM4(y0,y3,y7,y4))
    dzi = ( 0.25_RLK)*(SUM4(z1,z2,z6,z5) - SUM4(z0,z3,z7,z4))
    
    dxk = ( 0.25_RLK)*(SUM4(x4,x5,x6,x7) - SUM4(x0,x1,x2,x3))
    dyk = ( 0.25_RLK)*(SUM4(y4,y5,y6,y7) - SUM4(y0,y1,y2,y3))
    dzk = ( 0.25_RLK)*(SUM4(z4,z5,z6,z7) - SUM4(z0,z1,z2,z3))
    
!   find delvk and delxk ( i cross j )
    
    ax = dyi*dzj - dzi*dyj
    ay = dzi*dxj - dxi*dzj
    az = dxi*dyj - dyi*dxj
    
    domain%m_delx_zeta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny)
    
    ax = ax * norm
    ay = ay * norm
    az = az * norm
    
    dxv = (0.25_RLK)*(SUM4(xv4,xv5,xv6,xv7) - SUM4(xv0,xv1,xv2,xv3))
    dyv = (0.25_RLK)*(SUM4(yv4,yv5,yv6,yv7) - SUM4(yv0,yv1,yv2,yv3))
    dzv = (0.25_RLK)*(SUM4(zv4,zv5,zv6,zv7) - SUM4(zv0,zv1,zv2,zv3))
    
    domain%m_delv_zeta(i) = ax*dxv + ay*dyv + az*dzv
    
!   find delxi and delvi ( j cross k )
    
    ax = dyj*dzk - dzj*dyk ;
    ay = dzj*dxk - dxj*dzk ;
    az = dxj*dyk - dyj*dxk ;
    
    domain%m_delx_xi(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
    
    ax = ax * norm
    ay = ay * norm
    az = az * norm
    
    dxv = (0.25_RLK)*(SUM4(xv1,xv2,xv6,xv5) - SUM4(xv0,xv3,xv7,xv4)) ;
    dyv = (0.25_RLK)*(SUM4(yv1,yv2,yv6,yv5) - SUM4(yv0,yv3,yv7,yv4)) ;
    dzv = (0.25_RLK)*(SUM4(zv1,zv2,zv6,zv5) - SUM4(zv0,zv3,zv7,zv4)) ;
    
    domain%m_delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;
    
!   find delxj and delvj ( k cross i )
    
    ax = dyk*dzi - dzk*dyi ;
    ay = dzk*dxi - dxk*dzi ;
    az = dxk*dyi - dyk*dxi ;
    
    domain%m_delx_eta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
    
    ax = ax * norm
    ay = ay * norm
    az = az * norm
    
    dxv = (-0.25_RLK)*(SUM4(xv0,xv1,xv5,xv4) - SUM4(xv3,xv2,xv6,xv7)) ;
    dyv = (-0.25_RLK)*(SUM4(yv0,yv1,yv5,yv4) - SUM4(yv3,yv2,yv6,yv7)) ;
    dzv = (-0.25_RLK)*(SUM4(zv0,zv1,zv5,zv4) - SUM4(zv3,zv2,zv6,zv7)) ;
    
    domain%m_delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
  ENDDO

END SUBROUTINE CalcMonotonicQGradientsForElems




SUBROUTINE CalcMonotonicQRegionForElems( qlc_monoq, qqc_monoq,                &
                                         monoq_limiter_mult, monoq_max_slope, &
                                         ptiny,                               &
                                         elength                              ) 

  IMPLICIT NONE 
  REAL(KIND=8) :: qlc_monoq,  qqc_monoq
  REAL(KIND=8) :: monoq_limiter_mult,  monoq_max_slope
  REAL(KIND=8) :: ptiny
  INTEGER(KIND=4) :: elength     ! the elementset length
  INTEGER(KIND=4) :: ielem, i, bcMask
  REAL(KIND=8) :: qlin, qquad, phixi, phieta, phizeta, delvm, delvp
  REAL(KIND=8) :: norm, delvxxi, delvxeta, delvxzeta, rho


  DO ielem = 0, elength-1
    i = domain%m_matElemlist(ielem)
    bcMask = domain%m_elemBC(i)

!   phixi
    norm = (1.0_RLK) / ( domain%m_delv_xi(i) + ptiny )

    SELECT CASE(IAND(bcMask, XI_M))
      CASE (0)
        delvm = domain%m_delv_xi(domain%m_lxim(i))
      CASE (XI_M_SYMM)
        delvm = domain%m_delv_xi(i)
      CASE (XI_M_FREE)
        delvm = (0.0_RLK)
      CASE DEFAULT
!       ERROR
    END SELECT

    SELECT CASE(IAND(bcMask, XI_P))
      CASE (0)
        delvp = domain%m_delv_xi(domain%m_lxip(i))
      CASE (XI_P_SYMM)
        delvp = domain%m_delv_xi(i)
      CASE (XI_P_FREE)
        delvp = (0.0_RLK)
      CASE DEFAULT
!       ERROR 
    END SELECT

    delvm = delvm * norm
    delvp = delvp * norm

    phixi = (0.5_RLK) * ( delvm + delvp )
    
    delvm = delvm * monoq_limiter_mult
    delvp = delvp * monoq_limiter_mult
    
    if ( delvm < phixi ) phixi = delvm
    if ( delvp < phixi ) phixi = delvp
    if ( phixi < 0.0_RLK ) phixi = (0.0_RLK)
    if ( phixi > monoq_max_slope) phixi = monoq_max_slope
    
    
!   phieta
    norm = (1.0_RLK) / ( domain%m_delv_eta(i) + ptiny )
    
    SELECT CASE(IAND(bcMask, ETA_M))
      CASE (0)
        delvm = domain%m_delv_eta(domain%m_letam(i))
      CASE (ETA_M_SYMM)
        delvm = domain%m_delv_eta(i)
      CASE (ETA_M_FREE)
        delvm = 0.0_RLK
      CASE DEFAULT
!       ERROR
    END SELECT
    SELECT CASE(IAND(bcMask, ETA_P))
      CASE (0)
        delvp = domain%m_delv_eta(domain%m_letap(i))
      CASE (ETA_P_SYMM)
        delvp = domain%m_delv_eta(i)
      CASE (ETA_P_FREE)
        delvp = (0.0_RLK)
      CASE DEFAULT
!       ERROR
    END SELECT
    
    delvm = delvm * norm
    delvp = delvp * norm
    
    phieta = (0.5_RLK) * ( delvm + delvp )
    
    delvm = delvm * monoq_limiter_mult
    delvp = delvp * monoq_limiter_mult
    
    if ( delvm  < phieta ) phieta = delvm
    if ( delvp  < phieta ) phieta = delvp
    if ( phieta < (0.0_RLK)) phieta = (0.0_RLK)
    if ( phieta > monoq_max_slope)  phieta = monoq_max_slope
    
!   phizeta
    norm = (1.0_RLK) / ( domain%m_delv_zeta(i) + ptiny ) ;
    
    SELECT CASE(IAND(bcMask, ZETA_M))
      CASE (0)
        delvm = domain%m_delv_zeta(domain%m_lzetam(i))
      CASE (ZETA_M_SYMM)
        delvm = domain%m_delv_zeta(i)
      CASE (ZETA_M_FREE)
        delvm = (0.0_RLK)
      CASE DEFAULT
!       ERROR
    END SELECT
    SELECT CASE(IAND(bcMask, ZETA_P))
      CASE (0)
        delvp = domain%m_delv_zeta(domain%m_lzetap(i))
      CASE (ZETA_P_SYMM)
        delvp = domain%m_delv_zeta(i)
      CASE (ZETA_P_FREE)
        delvp = (0.0_RLK)
      CASE DEFAULT
!       ERROR
    END SELECT

    delvm = delvm * norm
    delvp = delvp * norm

    phizeta = (0.5_RLK) * ( delvm + delvp )

    delvm = delvm * monoq_limiter_mult
    delvp = delvp * monoq_limiter_mult

    IF ( delvm   < phizeta ) phizeta = delvm
    IF ( delvp   < phizeta ) phizeta = delvp
    IF ( phizeta < (0.0_RLK) ) phizeta = (0.0_RLK)
    IF ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope

!   Remove length scale

    IF ( domain%m_vdov(i) > (0.0_RLK) ) THEN
      qlin  = (0.0_RLK)
      qquad = (0.0_RLK)
    ELSE
      delvxxi   = domain%m_delv_xi(i)   * domain%m_delx_xi(i)
      delvxeta  = domain%m_delv_eta(i)  * domain%m_delx_eta(i)
      delvxzeta = domain%m_delv_zeta(i) * domain%m_delx_zeta(i)

      IF ( delvxxi   > (0.0_RLK) ) delvxxi   = (0.0_RLK)
      IF ( delvxeta  > (0.0_RLK) ) delvxeta  = (0.0_RLK)
      IF ( delvxzeta > (0.0_RLK) ) delvxzeta = (0.0_RLK)

      rho = domain%m_elemMass(i) / (domain%m_volo(i) * domain%m_vnew(i))

      qlin = -qlc_monoq * rho *                      &
             (  delvxxi   * ((1.0_RLK) - phixi)  +     &
                delvxeta  * ((1.0_RLK) - phieta) +     &
                delvxzeta * ((1.0_RLK) - phizeta)  )

      qquad = qqc_monoq * rho *                                       &
             (  delvxxi*delvxxi     * ((1.0_RLK) - phixi*phixi)   +     &
                delvxeta*delvxeta   * ((1.0_RLK) - phieta*phieta) +     &
                delvxzeta*delvxzeta * ((1.0_RLK) - phizeta*phizeta)  )
    ENDIF

    domain%m_qq(i) = qquad
    domain%m_ql(i) = qlin
  ENDDO

END SUBROUTINE CalcMonotonicQRegionForElems




SUBROUTINE CalcMonotonicQForElems()

  IMPLICIT NONE 
  REAL(KIND=8), PARAMETER :: ptiny = 1.e-36_RLK
  REAL(KIND=8) :: monoq_max_slope
  REAL(KIND=8) :: monoq_limiter_mult
  REAL(KIND=8) :: qlc_monoq
  REAL(KIND=8) :: qqc_monoq
  INTEGER(KIND=4) :: elength     ! the elementset length
!
! initialize parameters
!
  monoq_max_slope    = domain%m_monoq_max_slope
  monoq_limiter_mult = domain%m_monoq_limiter_mult

!
! calculate the monotonic q for pure regions
!
  elength = domain%m_numElem
  IF (elength > 0) THEN
    qlc_monoq = domain%m_qlc_monoq
    qqc_monoq = domain%m_qqc_monoq
    CALL CalcMonotonicQRegionForElems( qlc_monoq, qqc_monoq, &
                                       monoq_limiter_mult,   &
                                       monoq_max_slope,      &
                                       ptiny, elength )
  ENDIF

END SUBROUTINE CalcMonotonicQForElems






SUBROUTINE CalcQForElems()

  IMPLICIT NONE 
  REAL(KIND=8)    :: qstop
  INTEGER(KIND=4) :: numElem, idx, i

  qstop = domain%m_qstop
  numElem = domain%m_numElem
!
! MONOTONIC Q option
!

! Calculate velocity gradients
  CALL CalcMonotonicQGradientsForElems()

! Transfer veloctiy gradients in the first order elements
! problem->commElements->Transfer(CommElements::monoQ)
  CALL CalcMonotonicQForElems()

! Don't allow excessive artificial viscosity
  IF (numElem /= 0) THEN
    idx = -1
    DO i = 0, numElem-1
      IF ( domain%m_q(i) > qstop ) THEN
        idx = i
        EXIT
      ENDIF
    ENDDO

    IF (idx >= 0) THEN
      CALL luabort(QStopError)
    ENDIF
  ENDIF

END SUBROUTINE CalcQForElems






SUBROUTINE CalcPressureForElems( p_new, bvc,         &
                                 pbvc, e_old,        &
                                 compression, vnewc, &
                                 pmin,               &
                                 p_cut,eosvmax,      &
                                 length              )

  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:) :: p_new, bvc, pbvc, e_old
  REAL(KIND=8), DIMENSION(0:) ::  compression
  REAL(KIND=8), DIMENSION(0:) :: vnewc
  REAL(KIND=8)    :: pmin
  REAL(KIND=8)    :: p_cut
  REAL(KIND=8)    :: eosvmax
  INTEGER(KIND=4) :: length 

  INTEGER(KIND=4) :: i
  REAL(KIND=8), PARAMETER :: c1s = (2.0_RLK)/(3.0_RLK)

  DO i = 0, length-1
    bvc(i) = c1s * (compression(i) + (1.0_RLK))
    pbvc(i) = c1s
  ENDDO

  DO i = 0, length-1
    p_new(i) = bvc(i) * e_old(i)

    IF (ABS(p_new(i)) < p_cut) THEN
      p_new(i) = (0.0_RLK)
    ENDIF

    IF ( vnewc(i) >= eosvmax ) THEN  ! impossible condition here?
      p_new(i) = (0.0_RLK)
    ENDIF

    IF (p_new(i) < pmin) THEN
      p_new(i) = pmin
    ENDIF
  ENDDO

END SUBROUTINE CalcPressureForElems




SUBROUTINE  CalcEnergyForElems( p_new,  e_new,  q_new,          &
                                bvc,  pbvc,                     &
                                p_old,  e_old,  q_old,          &
                                compression,  compHalfStep,     &
                                vnewc,  work,  delvc,  pmin,    &
                                p_cut,   e_cut,  q_cut,  emin,  &
                                qq,  ql,                        &
                                rho0,                           &
                                eosvmax,                        &
                                length                          ) 

  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:) :: p_new, e_new, q_new
  REAL(KIND=8), DIMENSION(0:) :: bvc,  pbvc
  REAL(KIND=8), DIMENSION(0:) :: p_old, e_old, q_old
  REAL(KIND=8), DIMENSION(0:) :: compression, compHalfStep
  REAL(KIND=8), DIMENSION(0:) :: vnewc, work, delvc
  REAL(KIND=8)    :: pmin, p_cut,  e_cut, q_cut, emin
  REAL(KIND=8), DIMENSION(0:) :: qq, ql
  REAL(KIND=8)    :: rho0
  REAL(KIND=8)    :: eosvmax
  INTEGER(KIND=4) :: length

  INTEGER(KIND=4) :: i
  REAL(KIND=8)    :: vhalf, ssc, q_tilde
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: pHalfStep
  REAL(KIND=8), PARAMETER :: TINY1 = 0.111111e-36_RLK
  REAL(KIND=8), PARAMETER :: TINY3 = 0.333333e-18_RLK
  REAL(KIND=8), PARAMETER :: SIXTH = (1.0_RLK) / (6.0_RLK)


  ALLOCATE(pHalfStep(0:length-1))

  DO i = 0, length-1
    e_new(i) = e_old(i) - (0.5_RLK) * delvc(i) * (p_old(i) + q_old(i))  &
             + (0.5_RLK) * work(i)

    IF (e_new(i)  < emin ) THEN
      e_new(i) = emin
    ENDIF
  ENDDO

  CALL CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep,  &
                            vnewc, pmin, p_cut, eosvmax, length)

  DO i = 0, length-1
    vhalf = (1.0_RLK) / ((1.0_RLK) + compHalfStep(i))
    
    IF ( delvc(i) > (0.0_RLK) ) THEN
!      q_new(i) /* = qq(i) = ql(i) */ = Real_t(0.) ;
      q_new(i) = (0.0_RLK)
    ELSE
      ssc = ( pbvc(i) * e_new(i)   &
          + vhalf * vhalf * bvc(i) * pHalfStep(i) ) / rho0

      IF ( ssc <= TINY1 ) THEN
        ssc = TINY3
      ELSE
        ssc = SQRT(ssc)
      ENDIF

      q_new(i) = (ssc*ql(i) + qq(i))
    ENDIF

    e_new(i) = e_new(i) + (0.5_RLK) * delvc(i)   &
       * (  (3.0_RLK)*(p_old(i)     + q_old(i))  &
          - (4.0_RLK)*(pHalfStep(i) + q_new(i)))
  ENDDO
  
  DO i = 0, length-1
    e_new(i) = e_new(i) + (0.5_RLK) * work(i)

    IF (ABS(e_new(i)) < e_cut) THEN
      e_new(i) = (0.0_RLK)
    ENDIF
    IF (e_new(i)  < emin ) THEN
      e_new(i) = emin
    ENDIF
  ENDDO

  CALL CalcPressureForElems(p_new, bvc, pbvc, e_new, compression,  &
                            vnewc, pmin, p_cut, eosvmax, length)

  DO i = 0, length-1
    IF (delvc(i) > (0.0_RLK)) THEN
      q_tilde = (0.0_RLK)
    ELSE
      ssc = ( pbvc(i) * e_new(i)     &
          + vnewc(i) * vnewc(i) * bvc(i) * p_new(i) ) / rho0

      IF ( ssc <= TINY1 ) THEN
        ssc = TINY3
      ELSE
        ssc = SQRT(ssc)
      ENDIF

      q_tilde = (ssc*ql(i) + qq(i))
    ENDIF

    e_new(i) = e_new(i) - (  (7.0_RLK)*(p_old(i)     + q_old(i))   &
                        -    (8.0_RLK)*(pHalfStep(i) + q_new(i))   &
                        + (p_new(i) + q_tilde)) * delvc(i)*SIXTH

    IF (ABS(e_new(i)) < e_cut) THEN
      e_new(i) = (0.0_RLK)
    ENDIF
    IF ( e_new(i)  < emin ) THEN
      e_new(i) = emin
    ENDIF
  ENDDO

  CALL CalcPressureForElems(p_new, bvc, pbvc, e_new, compression,  &
                            vnewc, pmin, p_cut, eosvmax, length)

  DO i = 0, length-1

    IF ( delvc(i) <= (0.0_RLK) ) THEN
      ssc = ( pbvc(i) * e_new(i)  &
          + vnewc(i) * vnewc(i) * bvc(i) * p_new(i) ) / rho0

      IF ( ssc <= TINY1 ) THEN
        ssc = TINY3
      ELSE
        ssc = SQRT(ssc)
      ENDIF

      q_new(i) = (ssc*ql(i) + qq(i))

      if (ABS(q_new(i)) < q_cut) q_new(i) = (0.0_RLK)
    ENDIF
  ENDDO

  DEALLOCATE(pHalfStep)

END SUBROUTINE CalcEnergyForElems




SUBROUTINE CalcSoundSpeedForElems(vnewc,  rho0, enewc, &
                                  pnewc, pbvc,         &
                                  bvc, ss4o3, nz       )
  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:) :: vnewc, enewc
  REAL(KIND=8), DIMENSION(0:) :: pnewc, pbvc
  REAL(KIND=8), DIMENSION(0:) :: bvc
  REAL(KIND=8) :: rho0
  REAL(KIND=8) :: ss4o3
  INTEGER      :: nz
  REAL(KIND=8), PARAMETER :: TINY1 = 0.111111e-36_RLK
  REAL(KIND=8), PARAMETER :: TINY3 = 0.333333e-18_RLK
  REAL(KIND=8) :: ssTmp
  INTEGER      :: i, iz

  DO i = 0, nz - 1
    iz = domain%m_matElemlist(i)
    ssTmp = (pbvc(i) * enewc(i) + vnewc(i) * vnewc(i) *  &
                         bvc(i) * pnewc(i)) / rho0
    IF (ssTmp <= TINY1) THEN
      ssTmp = TINY3
    ELSE
      ssTmp = SQRT(ssTmp)
    ENDIF
    domain%m_ss(iz) = ssTmp
  ENDDO

END SUBROUTINE CalcSoundSpeedForElems



SUBROUTINE EvalEOSForElems(vnewc, length)

  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:) :: vnewc
  INTEGER :: length

  REAL(KIND=8) :: e_cut, p_cut, ss4o3, q_cut
  REAL(KIND=8) :: eosvmax, eosvmin, pmin, emin, rho0
  REAL(KIND=8) :: vchalf
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: e_old, &
                delvc, p_old, q_old, compression,   &
                compHalfStep, qq, ql, work, p_new,  &
                e_new, q_new, bvc, pbvc
  INTEGER      :: i, zidx

  e_cut = domain%m_e_cut
  p_cut = domain%m_p_cut
  ss4o3 = domain%m_ss4o3
  q_cut = domain%m_q_cut

  eosvmax = domain%m_eosvmax
  eosvmin = domain%m_eosvmin
  pmin    = domain%m_pmin
  emin    = domain%m_emin
  rho0    = domain%m_refdens

  ALLOCATE(e_old(0:length-1))
  ALLOCATE(delvc(0:length-1))
  ALLOCATE(p_old(0:length-1))
  ALLOCATE(q_old(0:length-1))
  ALLOCATE(compression(0:length-1))
  ALLOCATE(compHalfStep(0:length-1))
  ALLOCATE(qq(0:length-1))
  ALLOCATE(ql(0:length-1))
  ALLOCATE(work(0:length-1))
  ALLOCATE(p_new(0:length-1))
  ALLOCATE(e_new(0:length-1))
  ALLOCATE(q_new(0:length-1))
  ALLOCATE(bvc(0:length-1))
  ALLOCATE(pbvc(0:length-1))

! compress data, minimal set
  DO i = 0, length-1
    zidx = domain%m_matElemlist(i) ;
    e_old(i) = domain%m_e(zidx)
  ENDDO
  
  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    delvc(i) = domain%m_delv(zidx)
  ENDDO
  
  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    p_old(i) = domain%m_p(zidx)
  ENDDO
  
  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    q_old(i) = domain%m_q(zidx)
  ENDDO
  
  DO i = 0, length-1
    compression(i) = (1.0_RLK) / vnewc(i) - (1.0_RLK)
    vchalf = vnewc(i) - delvc(i) * (0.5_RLK)
    compHalfStep(i) = (1.0_RLK) / vchalf - (1.0_RLK)
  ENDDO
  
! Check for v > eosvmax or v < eosvmin
  IF ( eosvmin /= (0.0_RLK) ) THEN
    DO i = 0, length-1
      IF (vnewc(i) <= eosvmin) THEN  ! impossible due to calling func?
        compHalfStep(i) = compression(i)
      ENDIF
    ENDDO
  ENDIF
  IF ( eosvmax /= (0.0_RLK) ) THEN
    DO i = 0, length-1
      IF (vnewc(i) >= eosvmax) THEN ! impossible due to calling func? 
        p_old(i)        = (0.0_RLK)
        compression(i)  = (0.0_RLK)
        compHalfStep(i) = (0.0_RLK)
      ENDIF
    ENDDO
  ENDIF

  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    qq(i) = domain%m_qq(zidx)
    ql(i) = domain%m_ql(zidx)
    work(i) = (0.0_RLK)
  ENDDO

  CALL CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,          &
                          p_old, e_old,  q_old, compression,       &
                          compHalfStep, vnewc, work,  delvc, pmin, &
                          p_cut, e_cut, q_cut, emin,               &
                          qq, ql, rho0, eosvmax, length)


  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    domain%m_p(zidx) = p_new(i)
  ENDDO

  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    domain%m_e(zidx) = e_new(i)
  ENDDO

  DO i = 0, length-1
    zidx = domain%m_matElemlist(i)
    domain%m_q(zidx) = q_new(i)
  ENDDO


  CALL CalcSoundSpeedForElems(vnewc, rho0, e_new, p_new,  &
                              pbvc, bvc, ss4o3, length)

  DEALLOCATE(pbvc)
  DEALLOCATE(bvc)
  DEALLOCATE(q_new)
  DEALLOCATE(e_new)
  DEALLOCATE(p_new)
  DEALLOCATE(work)
  DEALLOCATE(ql)
  DEALLOCATE(qq)
  DEALLOCATE(compHalfStep)
  DEALLOCATE(compression)
  DEALLOCATE(q_old)
  DEALLOCATE(p_old)
  DEALLOCATE(delvc)
  DEALLOCATE(e_old)


END SUBROUTINE EvalEOSForElems


SUBROUTINE ApplyMaterialPropertiesForElems()

  IMPLICIT NONE 
  REAL(KIND=8)    :: eosvmin, eosvmax, vc
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vnewc
  INTEGER(KIND=4) :: length, zn

  length = domain%m_numElem

  IF (length /= 0) THEN
!   Expose all of the variables needed for material evaluation
    eosvmin = domain%m_eosvmin
    eosvmax = domain%m_eosvmax
    ALLOCATE(vnewc(0:length-1))

    DO i = 0, length-1
      zn = domain%m_matElemlist(i)
      vnewc(i) = domain%m_vnew(zn)
    ENDDO

    IF (eosvmin /= (0.0_RLK)) THEN
      DO i = 0, length-1
        IF (vnewc(i) < eosvmin) vnewc(i) = eosvmin
      ENDDO
    ENDIF

    IF (eosvmax /= (0.0_RLK)) THEN
      DO i = 0, length-1
        IF (vnewc(i) > eosvmax) vnewc(i) = eosvmax
      ENDDO
    ENDIF

    DO i = 0, length-1
      zn = domain%m_matElemlist(i)
      vc = domain%m_v(zn)
      IF (eosvmin /= (0.0_RLK)) THEN
        IF (vc < eosvmin) vc = eosvmin
      ENDIF
      IF (eosvmax /= (0.0_RLK)) THEN
        IF (vc > eosvmax) vc = eosvmax
      ENDIF
      IF (vc <= 0.0_RLK) THEN
        CALL luabort(VolumeError)
      ENDIF
    ENDDO

    CALL EvalEOSForElems(vnewc, length)

    DEALLOCATE(vnewc)

  ENDIF

END SUBROUTINE ApplyMaterialPropertiesForElems


SUBROUTINE UpdateVolumesForElems()

  IMPLICIT NONE 
  ReAL(KIND=8)    :: v_cut, tmpV
  INTEGER(KIND=4) :: numElem, i

  numElem = domain%m_numElem

  IF (numElem /= 0) THEN
    v_cut = domain%m_v_cut

    DO i = 0, numElem - 1
      tmpV = domain%m_vnew(i)

      IF ( ABS(tmpV - (1.0_RLK)) < v_cut ) THEN
        tmpV = (1.0_RLK)
      ENDIF
      domain%m_v(i) = tmpV
    ENDDO
  ENDIF

END SUBROUTINE UpdateVolumesForElems



SUBROUTINE LagrangeElements()

  IMPLICIT NONE 
  REAL(KIND=8) :: deltatime

  deltatime = domain%m_deltatime

  CALL CalcLagrangeElements(deltatime)

! Calculate Q.  (Monotonic q option requires communication)
  CALL CalcQForElems()

  CALL ApplyMaterialPropertiesForElems()

  CALL UpdateVolumesForElems()

END SUBROUTINE LagrangeElements



SUBROUTINE CalcCourantConstraintForElems()

  IMPLICIT NONE
  REAL(KIND=8)    :: dtcourant
  INTEGER(KIND=4) :: COURANT_ELEM

  REAL(KIND=8) :: qqc, qqc2, dtf
  REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: dtcourant_per_thread
  INTEGER(KIND=4) :: length, threads, i, indx, thread_num
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: courant_elem_per_thread

  dtcourant    = 1.0e+20_RLK
  COURANT_ELEM = -1

  qqc = domain%m_qqc
  length = domain%m_numElem

  qqc2 = (64.0_RLK) * qqc * qqc

! For sequential run (non OpenMP, replace line above with next line.  
  threads = 1_4

  ALLOCATE(dtcourant_per_thread(0:threads-1))
  ALLOCATE(courant_elem_per_thread(0:threads-1))

  DO i = 0, threads-1
    courant_elem_per_thread(i) = -1
    dtcourant_per_thread(i) =  (1.0e+20_RLK)
  ENDDO


  DO i = 0, length-1
    indx = domain%m_matElemlist(i)

    dtf = domain%m_ss(indx) * domain%m_ss(indx)

    IF ( domain%m_vdov(indx) < (0.0_RLK) ) THEN

      dtf = dtf + qqc2 * domain%m_arealg(indx) * domain%m_arealg(indx)  &
                * domain%m_vdov(indx) * domain%m_vdov(indx)
    ENDIF

    dtf = SQRT(dtf)

    dtf = domain%m_arealg(indx) / dtf

!   determine minimum timestep with its corresponding elem
    IF (domain%m_vdov(indx) /= (0.0_RLK)) THEN

!     For sequential run (non OpenMP, replace line above with next line.
      thread_num = (0_4)


      IF ( dtf < dtcourant_per_thread(thread_num) ) THEN

        dtcourant_per_thread(thread_num) = dtf
        courant_elem_per_thread(thread_num) = indx
      ENDIF
    ENDIF

  ENDDO

  DO i = 0, threads-1
    IF(dtcourant_per_thread(i) < dtcourant) THEN
      dtcourant = dtcourant_per_thread(i)
      courant_elem =  courant_elem_per_thread(i)
    ENDIF
  ENDDO


! Don't try to register a time constraint if none of the elements
! were active
  IF (courant_elem /= -1) THEN
    domain%m_dtcourant = dtcourant
  ENDIF

  DEALLOCATE(dtcourant_per_thread)
  DEALLOCATE(courant_elem_per_thread)

  RETURN

END SUBROUTINE CalcCourantConstraintForElems

SUBROUTINE CalcHydroConstraintForElems()

  IMPLICIT NONE 

  REAL(KIND=8) :: dthydro
  REAL(KIND=8) :: dvovmax, dtdvov
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dthydro_per_thread
  INTEGER(KIND=4) :: hydro_elem
  INTEGER(KIND=4) :: threads, i, length, indx, thread_num

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: hydro_elem_per_thread

  dthydro = 1.0e+20_RLK
  hydro_elem = -1
  dvovmax = domain%m_dvovmax
  length = domain%m_numElem

! For sequential run (non OpenMP, replace line above with next line.
  threads = (1_4)

  ALLOCATE(dthydro_per_thread(0:threads-1))
  ALLOCATE(hydro_elem_per_thread(0:threads-1))

  DO i = 0, threads-1
    hydro_elem_per_thread(i) = hydro_elem
    dthydro_per_thread(i) = dthydro
  ENDDO


  DO i = 0, length-1
    indx = domain%m_matElemlist(i) ;

    IF (domain%m_vdov(indx) /= (0.0_RLK)) THEN
      dtdvov = dvovmax / (ABS(domain%m_vdov(indx))+(1.e-20_RLK))

!     For sequential run (non OpenMP, replace line above with next line.
      thread_num = (0)

      IF ( dthydro_per_thread(thread_num) > dtdvov ) THEN
        dthydro_per_thread(thread_num) = dtdvov
        hydro_elem_per_thread(thread_num) = indx
      ENDIF
    ENDIF
  ENDDO


  DO i = 0, threads-1
    IF (dthydro_per_thread(i) < dthydro) THEN
      dthydro = dthydro_per_thread(i)
      hydro_elem =  hydro_elem_per_thread(i)
    ENDIF
  ENDDO

  IF (hydro_elem /= -1) THEN
    domain%m_dthydro = dthydro
  ENDIF

  RETURN

END SUBROUTINE CalcHydroConstraintForElems



SUBROUTINE CalcTimeConstraintsForElems()

  IMPLICIT NONE

! evaluate time constraint
  CALL CalcCourantConstraintForElems()

! check hydro constraint
  CALL CalcHydroConstraintForElems()

END SUBROUTINE CalcTimeConstraintsForElems



SUBROUTINE LagrangeLeapFrog()

  IMPLICIT NONE

! calculate nodal forces, accelerations, velocities, positions, with
! applied boundary conditions and slide surface considerations

  CALL LagrangeNodal()

! calculate element quantities (i.e. velocity gradient & q), and update
! material states

  CALL LagrangeElements()
  CALL CalcTimeConstraintsForElems()

! CALL LagrangeRelease()  ! Creation/destruction of temps may be important to capture 

END SUBROUTINE LagrangeLeapFrog



SUBROUTINE luabort(errcode)
  IMPLICIT NONE
  INTEGER(KIND=4) :: errcode
  WRITE(6,*) "ERROR CODE: ", errcode
  WRITE(6,*) "ABORTING"
  STOP
END SUBROUTINE luabort


END PROGRAM lulesh


