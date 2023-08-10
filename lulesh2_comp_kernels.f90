MODULE lulesh2_comp_kernels

    IMPLICIT NONE
    PRIVATE

    ! Definition of the domain data structure
    TYPE domain_type
        ! Node-centered

        ! coordinates: m_coord[idx] in cpp
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_x
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_y
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_z

        ! velocities: m_vel[idx] in cpp
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_xd
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_yd
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_zd

        ! accelerations: m_acc[idx] in cpp
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_xdd
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_ydd
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_zdd

        ! forces: m_force[idx] in cpp
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_fx
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_fy
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_fz

        ! mass: m_nodalMass[idx] in cpp
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_nodalMass

        ! symmetry plane nodeset: m_symmX[idx] in cpp
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_symmX
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_symmY
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_symmZ
        LOGICAL :: m_symm_is_set

        ! Missing region information

        ! Size of region sets
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_regElemSize

        ! Region number per domain element
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_regNumList

        ! Region indexset
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_regElemList

        ! Keys to the slices of the region indexset
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_regElemKeys

        INTEGER, DIMENSION(:), ALLOCATABLE :: m_nodeElemCount
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_nodeElemStart
        INTEGER, DIMENSION(:, :), ALLOCATABLE :: m_nodeElemCornerList

        ! Element-centered

        ! Material indexset
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_matElemList
        ! elemToNode connectivity
        INTEGER, DIMENSION(:, :), POINTER :: m_nodeList

        ! Element connectivity across each face
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_lxim
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_lxip
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_letam
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_letap
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_lzetam
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_lzetap

        ! Symmetry/free-surface flags for each elem face
        INTEGER, DIMENSION(:), ALLOCATABLE :: m_elemBC

        ! Principal strains -- temporary
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_dxx
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_dyy
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_dzz

        ! Velocity gradient -- temporary
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delv_xi
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delv_eta
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delv_zeta

        ! Coordinate gradient -- temporary
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delx_xi
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delx_eta
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delx_zeta

        ! Energy
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_e

        ! Pressure
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_p

        ! q
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_q

        ! Linear term for q
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_ql

        ! Quadratic term for q
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_qq

        ! Relative volume
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_v

        ! Reference volume
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_volo

        ! New relative volume -- temporary
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_vnew

        ! Change in relative volume -- m_vnew - m_v
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_delv

        ! Volume derivative over volume
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_vdov

        ! Characteristic length of an element
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_arealg

        ! Speed of sound
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_ss

        ! Mass
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m_elemMass

        ! Parameters

        ! Fixed time increment
        REAL(KIND=8) :: m_dtfixed

        ! Current time
        REAL(KIND=8) :: m_time

        ! Variable time increment
        REAL(KIND=8) :: m_deltatime
        REAL(KIND=8) :: m_deltatimemultlb
        REAL(KIND=8) :: m_deltatimemultub

        ! End time for the simulation
        REAL(KIND=8) :: m_stoptime

        ! Velocity tolerance
        REAL(KIND=8) :: m_u_cut

        ! Hourglass control
        REAL(KIND=8) :: m_hgcoef

        ! Excessive q-indicator
        REAL(KIND=8) :: m_qstop
        REAL(KIND=8) :: m_monoq_max_slope
        REAL(KIND=8) :: m_monoq_limiter_mult

        ! Energy tolerance
        REAL(KIND=8) :: m_e_cut

        ! Pressure tolerance
        REAL(KIND=8) :: m_p_cut
        REAL(KIND=8) :: m_ss4o3

        ! Q tolerance
        REAL(KIND=8) :: m_q_cut

        ! Relative volume tolerance
        REAL(KIND=8) :: m_v_cut

        ! Linear term coefficient for q
        REAL(KIND=8) :: m_qlc_monoq

        ! Quadratic term coefficient for q
        REAL(KIND=8) :: m_qqc_monoq
        REAL(KIND=8) :: m_qqc
        REAL(KIND=8) :: m_eosvmax
        REAL(KIND=8) :: m_eosvmin

        ! Pressure floor
        REAL(KIND=8) :: m_pmin

        ! Energy floor
        REAL(KIND=8) :: m_emin

        ! Maximum allowable volume change
        REAL(KIND=8) :: m_dvovmax

        ! Reference density
        REAL(KIND=8) :: m_refdens

        ! Courant constraint
        REAL(KIND=8) :: m_dtcourant

        ! Volume change constraint
        REAL(KIND=8) :: m_dthydro

        ! Maximum allowable time increment
        REAL(KIND=8) :: m_dtmax

        ! Iterative count for simulation
        INTEGER :: m_cycle

        ! Number of regions
        INTEGER :: m_numReg

        ! Imbalance cost
        INTEGER :: m_cost

        ! X, Y, Z extent of this block
        INTEGER :: m_sizeX
        INTEGER :: m_sizeY
        INTEGER :: m_sizeXZ

        ! Elements/nodes in this domain
        INTEGER :: m_numEem
        INTEGER :: m_numNode

    END TYPE domain_type

    ! Publicly available subroutines, and Datastructs
    PUBLIC :: domain_type
    PUBLIC :: AllocateNodalPersistent
    PUBLIC :: AllocateElemPersistent
    PUBLIC :: AllocateGradients
    PUBLIC :: DeallocateGradients
    PUBLIC :: AllocateStrains
    PUBLIC :: DeallocateStrains
    PUBLIC :: AllocateNodesets
    PUBLIC :: AllocateNodeElemIndexes

    CONTAINS

    SUBROUTINE AllocateNodalPersistent(domain, numNode)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER :: numNode
        INTEGER(KIND=4), PARAMETER :: RLK = 8

        ! TODO: Implement the allocations
        ! Coordinates

        ! Velocities

        ! Acceleration

        ! Forces

    END SUBROUTINE AllocateNodalPersistent

    SUBROUTINE AllocateElemPersistent(domain, numElem)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER :: numElem
        INTEGER(KIND=4), PARAMETER :: RLK = 8

        ! TODO: Implement the persistent allocation

    END SUBROUTINE AllocateElemPersistent

    SUBROUTINE AllocateGradients(domain, numElem, allElem)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER :: numElem, allElem

        ! TODO: Implement gradient allocation

    END SUBROUTINE AllocateGradients

    SUBROUTINE DeallocateGradients(domain)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain

        ! TODO: Implement gradient deallocation

    END SUBROUTINE DeallocateGradients

    SUBROUTINE AllocateStrains(domain, numElem)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER :: numElem

        ! TODO: Implement strain allocation

    END SUBROUTINE AllocateStrains

    SUBROUTINE DeallocateStrains(domain)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain

        ! TODO: Implement strain deallocation

    END SUBROUTINE DeallocateStrains

    SUBROUTINE AllocateNodesets(domain, edgeNodes_sq)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER :: edgeNodes_sq

        ! TODO: Implement nodeset allocation -- the one kernel I am highly unsure about!

    END SUBROUTINE AllocateNodesets

    SUBROUTINE AllocateNodeElemIndexes(domain)
        IMPLICIT NONE

        TYPE(domain_type), INTENT(INOUT) :: domain
        INTEGER :: numElem
        INTEGER, DIMENSION(0:7) :: m
        INTEGER :: i, j, k
        INTEGER :: offset
        INTEGER :: clSize, clv
        INTEGER :: numNode
        INTEGER :: nodelist_len

        ! TODO: Implement NodeElemIndexes allocation -- one of the kernels I am more unsure about!

    END SUBROUTINE AllocateNodeElemIndexes

    SUBROUTINE InitMeshDecomp(numRanks, myRank, col, row, plane, side)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: numRanks
        INTEGER, INTENT(IN) :: myRank
        INTEGER, INTENT(OUT) :: col
        INTEGER, INTENT(OUT) :: row
        INTEGER, INTENT(OUT) :: plane
        INTEGER, INTENT(OUT) :: side

        INTEGER(KIND=4) :: testProcs
        INTEGER(KIND=4) :: dx, dy, dz
        INTEGER(KIND=4) :: myDom
        INTEGER(KIND=4) :: remainder
        INTEGER(KIND=4), PARAMETER :: RLK = 8

        ! TODO: Implement the initial mesh decomposition

    END SUBROUTINE InitMeshDecomp

END MODULE lulesh2_comp_kernels