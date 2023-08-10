PROGRAM lulesh2
! TODO: Create an action, to build both LULESH2's against each other --> should fail fow now.
    USE lulesh2_comp_kernels

    IMPLICIT NONE

    INTEGER(KIND=4), PARAMETER :: VolumeError = -1
    INTEGER(KIND=4), PARAMETER :: QStopError  = -2
    INTEGER(KIND=4), PARAMETER :: RLK = 8

    ! Begin main
    TYPE(domain_type) :: domain
    INTEGER :: edgeElems
    INTEGER :: edgeNodes
    INTEGER :: opts_its, opts_nx, opts_numReg, opts_showProg, &
            opts_quiet, opts_viz, opts_balance, opts_cost
    REAL(KIND=8) :: opts_numFiles
    REAL(KIND=8) :: tx, ty, tz
    INTEGER :: nidx, zidx
    INTEGER :: col, row, plane, side
    INTEGER :: domElems
    INTEGER(KIND=4) :: grad_domElems
    INTEGER :: i, j, k
    INTEGER :: planeInc, rowInc
    REAL(KIND=8), DIMENSION(0:7) :: x_local, y_local, z_local
    INTEGER :: gnode, lnode, idx
    REAL(KIND=8) :: volume
    REAL(KIND=8) :: starttime, endtime
    REAL(KIND=8) :: elapsed_time

    ! TODO: Unsure about this point definition
    INTEGER(KIND=4), DIMENSION(:), POINTER :: localNode => NULL()

    ! Initial energy to be deposited
    REAL(KIND=8), PARAMETER :: ebase = 3.948746e+7_RLK
    REAL(KIND=8) :: scale, einit

    ! Define the boundary conditions using hexahedral faces (12 bits)
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

    INTEGER :: ElemID

    REAL(KIND=8) :: MaxAbsDiff
    REAL(KIND=8) :: TotalAbsDiff
    REAL(KIND=8) :: MaxRelDiff

    REAL(KIND=8) :: AbsDiff, RelDiff

    INTEGER :: regionNum, regionVar, binSize, lastReg, &
            elements, runto, costDenominator
    INTEGER, DIMENSION(:), ALLOCATABLE :: regBinEnd

    ElemID = 0
    MaxAbsDiff = 0.0_RLK
    TotalAbsDiff = 0.0_RLK
    MaxRelDiff = 0.0_RLK

    ! Symmetry planes not set yet
    domain%m_symm_is_set = .FALSE.

    ! Debugging configuration
    ! TODO: To be replaced with the realistic configurations
    edgeElems = 2  ! Normally opts.nx
    edgeNodes = edgeElems + 1
    numRanks = 1  ! Serial execution
    myRank = 0  ! Rank of the executor

    ! Options as set in the C++ code
    ! Set defaults that can be overridden by command line args
    opts_its      = 9999999
    opts_nx       = edgeElems
    opts_numReg   = 11
    opts_numFiles = (numRanks + 10)/9
    opts_showProg = 0
    opts_quiet    = 0
    opts_viz      = 0
    opts_balance  = 1
    opts_cost     = 1

    ! TODO: Am I missing the domain initialization here??

    ! Set up mesh, and decompose
    CALL InitMeshDecomp(numRanks, myRank, col, row, plane, side)

    ! Construction of the domain
    edgeElems = opts_nx
    edgeNodes = edgeElems + 1

    ! Populate the domain data structure
    domain%m_tp = tp
    domain%m_numRanks = numRanks

    ! Initialize the Sedov mesh
    domain%m_sizeX = edgeElems
    domain%m_sizeY = edgeElems
    domain%m_sizeZ = edgeElems
    domain%m_numElem = edgeElems * edgeElems * edgeElems
    domain%m_numNode = edgeNodes * edgeNodes * edgeNodes

    ! Define the material indexset
    ALLOCATE(domain$m_regNumList(0:numElem-1))

    ! Element-centered
    CALL AllocateElemPersistent(domain, domain%m_numNode)
    CALL AllocateNodesets(domain, edgeNodes * edgeNodes)

    ! TODO: Debug to this point here, and make sure that all the data structure agree.

END PROGRAM lulesh2