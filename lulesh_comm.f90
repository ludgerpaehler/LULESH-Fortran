MODULE lulesh_comm

! if USE_MPI

    INCLUDE 'mpif.h'

! allow_unpacked_plane to be set to false
! allow_unpacked_row to be set to false
! allow_unpacked_col to be set to false

!--------------------------------------------------------------------
!   There are coherence issues for packing and unpacking message
!   buffers.  Ideally, you would like a lot of threads to 
!   cooperate in the assembly/dissassembly of each message.
!   To do that, each thread should really be operating in a
!   different coherence zone.
!
!   Let's assume we have three fields, f1 through f3, defined on
!   a 61x61x61 cube.  If we want to send the block boundary
!   information for each field to each neighbor processor across
!   each cube face, then we have three cases for the
!   memory layout/coherence of data on each of the six cube
!   boundaries:
!
!      (a) Two of the faces will be in contiguous memory blocks
!      (b) Two of the faces will be comprised of pencils of
!          contiguous memory.
!      (c) Two of the faces will have large strides between
!          every value living on the face.
!
!   How do you pack and unpack this data in buffers to
!   simultaneous achieve the best memory efficiency and
!   the most thread independence?
!
!   Do do you pack field f1 through f3 tighly to reduce message
!   size?  Do you align each field on a cache coherence boundary
!   within the message so that threads can pack and unpack each
!   field independently?  For case (b), do you align each
!   boundary pencil of each field separately?  This increases
!   the message size, but could improve cache coherence so
!   each pencil could be processed independently by a separate
!   thread with no conflicts.
!
!   Also, memory access for case (c) would best be done without
!   going through the cache (the stride is so large it just causes
!   a lot of useless cache evictions).  Is it worth creating
!   a special case version of the packing algorithm that uses
!   non-coherent load/store opcodes?
!--------------------------------------------------------------------



!    void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
!                  Index_t dx, Index_t dy, Index_t dz, bool doRecv, bool planeOnly) {


SUBROUTINE CommRecv(msgType, xferFields, dx, dy, dz, doRecv, planeOnly)
    IMPLICIT NONE

    ! Initialize the variables coming in
    INTEGER :: msgType
    INTEGER :: xferFields
    INTEGER :: dx, dy, dz
    LOGICAL :: doRecv, planeOnly

    IF (domain%num_ranks == 1) THEN
        RETURN

    !----------------------------------------------------------------
    ! Post receive buffers for all incoming messages
    !----------------------------------------------------------------

    INTEGER :: myRank
    INTEGER :: maxPlaneComm = xferFields * domain%maxPlaneSize()
    INTEGER :: maxEdgeComm  = xferFields * domain%maxEdgeSize()
    ! Plane communication messaging
    INTEGER :: pmsg = 0
    ! Edge communication messaging
    INTEGER :: emsg = 0
    ! Corner communication messaging
    INTEGER :: cmsg = 0
    TYPE(MPI_DOUBLE) :: baseType
    LOGICAL :: rowMin, rowMax, colMin, colMax, planeMin, planeMax

    ! Assuming communication to all 6 neighbors by default
    rowMin = .TRUE.
    rowMax = .TRUE.
    colMin = .TRUE.
    colMax = .TRUE.
    planeMin = .TRUE.
    planeMax = .TRUE.

    ! Unsure whether the addresses exist in the type and whether the indexing here is correct
    IF (domain%rowLoc() == 0) THEN
        rowMin = .FALSE.
    ENDIF
    IF (domain%rowLoc() == (domain%tp() - 1)) THEN
        rowMax = .FALSE.
    ENDIF
    IF (domain%colLoc() == 0) THEN
        colMin = .TRUE.
    ENDIF
    IF (domain%colLoc() == (domain%tp() - 1)) THEN
        colMax = .FALSE.
    ENDIF
    IF (domain%planeLoc() == 0) THEN
        planeMin = .FALSE.
    ENDIF
    IF (domain%planeLoc() == (domain%tp() - 1)) THEN
        planeMax = .FALSE.
    ENDIF

    DO i=1, 27  ! 1-based indexing as opposed to the 0-based indexing of C++
        domain%recvRequest[i] = MPI_REQUEST_NULL
    END DO

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank)

    !----------------------------------------------------------------
    ! Post Receives
    !----------------------------------------------------------------

    ! Receive data from neighboring domain faces
    IF (planeMin .AND. doRecv) THEN
        INTEGER :: fromRank = myRank - domain%tp() * domain%tp()
        INTEGER :: recvCount = dx * dy * xferFields
        CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm] &
                      , recvCount, baseType, fromRank, msgType  &
                      , MPI_COMM_WORLD, domain%recvRequest[pmsg])
        pmsg = pmsg + 1
    ENDIF

    IF (planeMax) THEN
        INTEGER :: fromRank = myRank + domain%tp() * domain%tp()
        INTEGER :: recvCount = dx * dy * xferFields
        CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm] &
                      , recvCount, baseType, fromRank, msgType  &
                      , MPI_COMM_WORLD, domain%recvRequest[pmsg])
        pmsg = pmsg + 1
    ENDIF

    IF (rowMin .AND. doRecv) THEN
        ! Semi-contiguous memory
        INTEGER :: fromRank = myRank - domain%tp()
        INTEGER :: recvCount = dx * dz * xferFields
        CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm] &
                      , recvCount, baseType, fromRank, msgType  &
                      , MPI_COMM_WORLD, domain%recvRequest[pmsg])
        pmsg = pmsg + 1
    ENDIF

    IF (rowMax) THEN
        ! Semi-contiguous memory
        INTEGER :: fromRank = myRank + domain%tp()
        INTEGER :: recvCount = dx * dz * xferFields
        CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm] &
                      , recvCount, baseType, fromRank, msgType  &
                      , MPI_COMM_WORLD, domain%recvRequest[pmsg])
        pmsg = pmsg + 1
    ENDIF

    IF (colMin .AND. doRecv) THEN
        ! Scattered memory
        INTEGER :: fromRank = myRank - 1
        INTEGER :: recvCount = dy * dz * xferFields
        CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm] &
                      , recvCount, baseType, fromRank, msgType  &
                      , MPI_COMM_WORLD, domain%recvRequest[pmsg])
        pmsg = pmsg + 1
    ENDIF

    IF (colMax) THEN
        ! Scattered memory
        INTEGER :: fromRank = myRank + 1
        INTEGER :: recvCount = dy * dz * xferFields
        CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm] &
                      , recvCount, baseType, fromRank, msgType  &
                      , MPI_COMM_WORLD, domain%recvRequest[pmsg])
        pmsg = pmsg + 1
    ENDIF


    IF (.NOT. planeOnly) THEN

        ! Receive data from domains connected only by an edge
        IF (rowMin .AND. colMin .AND. doRecv) THEN
            INTEGER :: fromRank = myRank - domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dz * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMin .AND. planeMin .AND. doRecv) THEN
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() - domain&tp()
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dx * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (colMin .AND. planeMin .AND. doRecv) THEN
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dy * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMax .AND. colMax) THEN
            INTEGER :: fromRank = myRank + domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dz * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMax .AND. planeMax) THEN
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() + domain%tp()
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dx * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (colMax .AND. planeMax) THEN
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() + 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dy * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMax .AND. colMin) THEN
            INTEGER :: fromRank = myRank + domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dz * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMin .AND. planeMax) THEN
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() - domain%tp()
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dx * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (colMin .AND. planeMax) THEN
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dy * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMin .AND. colMax .AND. doRecv) THEN
            INTEGER :: fromRank = myRank - domain%tp() + 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dz * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (rowMax .AND. planeMin .AND. doRecv) THEN
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() + domain%tp()
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dx * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        IF (colMax .AND. planeMin .AND. doRecv) THEN
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() + 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm] &
                          , dy * xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD &
                          , domain%recvRequest[pmsg+emsg])
            emsg = emsg + 1
        ENDIF

        !------------------------------------------------------------
        ! Receive data from domains connected only by a corner
        !------------------------------------------------------------
        IF (rowMin .AND. colMin .AND. planeMin .AND. doRecv) THEN
            ! Corner at domain logical coord (0, 0, 0)
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() - domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMin .AND. colMin .AND. planeMax) THEN
            ! Corner at domain logical coord (0, 0, 1)
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() - domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMin .AND. colMax .AND. planeMin .AND. doRecv) THEN
            ! Corner at domain logical coord (1, 0, 0)
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() - domain%tp() + 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMin .AND. colMax .AND. planeMax) THEN
            ! Corner at domain logical coord (1, 0, 1)
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() - domain%tp() + 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMax .AND. colMin .AND. planeMin .AND. doRecv) THEN
            ! Corner at domain logical coord (0, 1, 0)
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() + domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMax .AND. colMin .AND. planeMax) THEN
            ! Corner at domain logical coord (0, 1, 1)
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() + domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMax .AND. colMax .AND. planeMin .AND. doRecv) THEN
            ! Corner at domain logical coord (1, 1, 0)
            INTEGER :: fromRank = myRank - domain%tp() * domain%tp() + domain%tp() - 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

        IF (rowMax .AND. colMax .AND. planeMax) THEN
            ! Corner at domain logical coord (1, 1, 1)
            INTEGER :: fromRank = myRank + domain%tp() * domain%tp() + domain%tp() + 1
            CALL MPI_IRECV(domain%commDataRecv[pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL], xferFields, baseType, fromRank, msgType, MPI_COMM_WORLD, domain%recvRequest[pmsg+emsg+cmsg])
            cmsg = cmsg + 1
        ENDIF

    ENDIF

END SUBROUTINE CommRecv




SUBROUTINE CommSend
    IMPLICIT NONE


END SUBROUTINE CommSend




SUBROUTINE CommSBN
    IMPLICIT NONE


END SUBROUTINE CommSBN




SUBROUTINE CommSyncPosVel
    IMPLICIT NONE


END SUBROUTINE CommSyncPosVel




SUBROUTINE CommMonoQ


END SUBROUTINE CommMonoQ










END MODULE lulesh_comm


