!######################################################################
!#                               GATHER
!#          
!#          Assumes each element of type array has been calculated
!#          or completed proc_id = mod(n-1, proc_count) n=1,total
!#          
!#          
!######################################################################


!######################################################################
!#                               BROADCAST
!#          
!#          Broadcasts all the elements in the array from the root
!#          out to the workers.  The array must be there and allocated,
!#          but arrays that make up the type are allocated automatically.
!#          
!######################################################################


MODULE mpi_udgb

USE kinds
USE efs_type, ONLY: t_efs
USE mpi

IMPLICIT NONE


INTEGER(kind=StandardInteger), ALLOCATABLE, DIMENSION(:) :: bheader        ! HEADER DATA
INTEGER(kind=StandardInteger), ALLOCATABLE, DIMENSION(:) :: bint           ! INTEGERS
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:) ::         bdp            ! DOUBLES
LOGICAL, ALLOCATABLE, DIMENSION(:) ::                       blogical       ! LOGICALS
INTEGER(kind=StandardInteger) ::                            ecount         ! ARRAY ELEMENT COUNT
INTEGER(kind=StandardInteger) ::                            ecounter       ! ARRAY ELEMENT COUNT
INTEGER(kind=StandardInteger) ::                            kcount
INTEGER(kind=StandardInteger) ::                            intcount
INTEGER(kind=StandardInteger) ::                            dpcount
INTEGER(kind=StandardInteger) ::                            lcount
INTEGER(kind=StandardInteger) ::                            intoffset
INTEGER(kind=StandardInteger) ::                            kpos
INTEGER(kind=StandardInteger) ::                            ipos
INTEGER(kind=StandardInteger) ::                            dppos
INTEGER(kind=StandardInteger) ::                            lpos
INTEGER(kind=StandardInteger) ::                            proc_id
INTEGER(kind=StandardInteger) ::                            proc_count

! PACK/UNPACK
INTERFACE pup_g
MODULE PROCEDURE pup_g_int_1, pup_g_int_2, pup_g_int_3, &
                 pup_g_dp_1, pup_g_dp_2, pup_g_dp_3, &
                 pup_g_l_1, pup_g_l_2, pup_g_l_3
END INTERFACE pup_g

INTERFACE pup_b
MODULE PROCEDURE pup_b_int_1, pup_b_int_2, pup_b_int_3, &
                 pup_b_dp_1, pup_b_dp_2, pup_b_dp_3, &
                 pup_b_l_1, pup_b_l_2, pup_b_l_3
END INTERFACE pup_b


CONTAINS

SUBROUTINE run_gather_efs(arr)
!###########################################################################
TYPE(t_efs), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
!###########################################################################
INTEGER(kind=StandardInteger) :: n, m, loop
INTEGER(kind=StandardInteger) :: tag, error
INTEGER(kind=StandardInteger) :: worker_id
INTEGER(kind=StandardInteger), DIMENSION(MPI_STATUS_SIZE) :: status


!###########################################################################

! Get MPI details
Call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,error)
Call MPI_Comm_size(MPI_COMM_WORLD,proc_count,error)

! SET UP BHEADER - stores sizes of buffers
IF (ALLOCATED(bheader)) THEN
  DEALLOCATE(bheader)
END IF
ALLOCATE(bheader(1:3))

! LOOP THROUGH WORKERS
DO worker_id = 1, proc_count-1

  ! IF THIS WORKER OR ROOT
  IF((proc_id .EQ. 0) .OR. (proc_id .EQ. worker_id))THEN

    ! ROOT AND WORKER:  RESET BUFFER HEADER


    intcount = 0 
    dpcount = 0 
    lcount = 0 

    ! Loop 1 calc size of buffers
    ! Loop 2 pack and send from worker    
    ! Loop 3 receive by root       
    DO loop = 1,3 

      IF (loop .EQ. 1) THEN    
        bheader = 0  
        ecount = 0
        kcount = 0
      END IF  

      IF (loop .EQ. 2 .AND. proc_id .EQ. worker_id) THEN
        ! ARRAY SIZE
        bint(1) = SIZE(arr,1)   ! STORE TOTAL NUMBER OF ELEMENTS
        bint(2) = ecount        ! STORE TOTAL NUMBER OF ELEMENTS FROM THIS WORKER
        bint(3) = kcount        ! STORE TOTAL NUMBER OF ELEMENTS FROM THIS WORKER
        ! STORE WHICH ELEMENTS
        m = 1
        DO n = 1, size(arr, 1)
          IF (MOD((n - 1), proc_count) .EQ. proc_id) THEN
            bint(3 + m) = n
            m  = m + 1
          END IF
        END DO           
      END IF      

      ! ON WORKER:  
      IF (loop .EQ. 1 .AND. proc_id .EQ. worker_id) THEN
        ! CALC HOW MANY ELEMENTS ARE BEING SENT
        DO n = 1, size(arr, 1)
          IF (MOD((n - 1), proc_count) .EQ. proc_id) THEN
            ecount = ecount + 1
          END IF
        END DO   
      END IF

      IF((proc_id .EQ. 0) .AND. (loop .EQ. 3))THEN
        ecount = bint(2)  
        kcount = bint(3)        
      END IF

      ! ON ROOT: SET POSITIONS
      IF((proc_id .EQ. 0) .AND. (loop .EQ. 3))THEN
        kpos = 3 + ecount
      END IF

      ! LOOP OVER ELEMENTS IN ARR
      DO n = 1, size(arr, 1)
        IF (MOD((n - 1), proc_count) .EQ. worker_id) THEN
          IF(loop .EQ. 2)THEN
            ecounter = ecounter + 1
          END IF
!############################################################################################
!############################################################################################

          CALL pup_g(arr(n)%counter, loop) 
          CALL pup_g(arr(n)%proc_id, loop) 
          CALL pup_g(arr(n)%proc_count, loop) 
          CALL pup_g(arr(n)%e_on, loop) 
          CALL pup_g(arr(n)%f_on, loop) 
          CALL pup_g(arr(n)%s_on, loop) 
          CALL pup_g(arr(n)%density, loop) 
          CALL pup_g(arr(n)%energy, loop) 
          CALL pup_g(arr(n)%forces, loop) 
          CALL pup_g(arr(n)%stress, loop) 
          CALL pup_g(arr(n)%energy_pair, loop) 
          CALL pup_g(arr(n)%energy_embed, loop) 


!############################################################################################
!############################################################################################
        END IF
      END DO

      ! ON WORKER: ALLOCATE BUFFERS
      IF((proc_id .EQ. worker_id) .AND. (loop .EQ. 1))THEN
        CALL deallocate_arrays()
        kpos = 3 + ecount
        ipos = 3 + ecount + kcount
        dppos = 1
        lpos = 1
        intoffset = 3 + ecount + kcount
        bheader(1) = intoffset + intcount
        bheader(2) = dpcount
        bheader(3) = lcount
        ALLOCATE(bint(1:bheader(1)))
        ALLOCATE(bdp(1:bheader(2)))
        ALLOCATE(blogical(1:bheader(3)))
        bint = 0
        bdp = 0.0D0
        blogical = .TRUE.       
      END IF

      ! SEND BUFFER HEADER TO ROOT      
      IF(loop .EQ. 1)THEN
        tag = 100 * worker_id + 1
        IF (proc_id .GT. 0) THEN
          CALL MPI_SEND(bheader, 3, MPI_INTEGER, 0, tag, MPI_comm_world, error)   
        ELSE 
          CALL MPI_RECV(bheader, 3, MPI_INTEGER, worker_id, tag, MPI_comm_world, status, error)   
        END IF
      END IF

      ! ON ROOT: ALLOCATE BUFFERS
      IF((proc_id .EQ. 0) .AND. (loop .EQ. 1))THEN
        CALL deallocate_arrays()
        ALLOCATE(bint(1:bheader(1)))
        ALLOCATE(bdp(1:bheader(2)))
        ALLOCATE(blogical(1:bheader(3)))
        bint = 0
        bdp = 0.0D0
        blogical = .TRUE.
      END IF


      ! SEND BUFFERS TO ROOT      
      IF(loop .EQ. 2)THEN
        IF (proc_id .GT. 0) THEN
          tag = 100 * worker_id + 2
          CALL MPI_SEND(bint, bheader(1), MPI_INTEGER, 0, tag, MPI_comm_world, error)  
          tag = 100 * worker_id + 3 
          CALL MPI_SEND(bdp, bheader(2), MPI_DOUBLE_PRECISION, 0, tag, MPI_comm_world, error)  
          tag = 100 * worker_id + 4
          CALL MPI_SEND(blogical, bheader(3), MPI_LOGICAL, 0, tag, MPI_comm_world, error)  
        ELSE 
          tag = 100 * worker_id + 2
          CALL MPI_RECV(bint, bheader(1), MPI_INTEGER, worker_id, tag, MPI_comm_world, status, error)
          tag = 100 * worker_id + 3
          CALL MPI_RECV(bdp, bheader(2), MPI_DOUBLE_PRECISION, worker_id, tag, MPI_comm_world, status, error)
          tag = 100 * worker_id + 4
          CALL MPI_RECV(blogical, bheader(3), MPI_LOGICAL, worker_id, tag, MPI_comm_world, status, error)
        END IF
      END IF      
    END DO
  END IF
END DO

END SUBROUTINE run_gather_efs









SUBROUTINE run_broadcast_efs(arr)
!###########################################################################
TYPE(t_efs), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
!###########################################################################
INTEGER(kind=StandardInteger) :: n, m, loop
INTEGER(kind=StandardInteger) :: tag, error
INTEGER(kind=StandardInteger) :: worker_id
INTEGER(kind=StandardInteger), DIMENSION(MPI_STATUS_SIZE) :: status


!###########################################################################

! Get MPI details
Call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,error)
Call MPI_Comm_size(MPI_COMM_WORLD,proc_count,error)

! SET UP BHEADER - stores sizes of buffers
IF (ALLOCATED(bheader)) THEN
  DEALLOCATE(bheader)
END IF
ALLOCATE(bheader(1:3))


! ROOT AND WORKER:  RESET BUFFER HEADER

intcount = 0 
dpcount = 0 
lcount = 0 

! Loop 1 calc size of buffers
! Loop 2 pack and broadcast from root    
! Loop 3 receive and unpack by workers
DO loop = 1,3

  IF (loop .EQ. 1) THEN    
    bheader = 0  
    ecount = 0
    kcount = 0
  END IF  

  IF (loop .EQ. 2 .AND. proc_id .EQ. 0) THEN
    ! ARRAY SIZE
    bint(1) = SIZE(arr,1)   ! STORE TOTAL NUMBER OF ELEMENTS
    bint(2) = ecount        ! STORE TOTAL NUMBER OF ELEMENTS FROM THIS WORKER
    bint(3) = kcount        ! STORE TOTAL NUMBER OF ELEMENTS FROM THIS WORKER
    ! STORE WHICH ELEMENTS
    m = 1
    DO n = 1, size(arr, 1)
      bint(3 + m) = n
      m  = m + 1
    END DO           
  END IF      

  ! ON ROOT: COUNT ELEMENTS BEING SENT (ALL)  
  IF (loop .EQ. 1 .AND. proc_id .EQ. 0) THEN
    ! CALC HOW MANY ELEMENTS ARE BEING SENT
    DO n = 1, size(arr, 1)
      ecount = ecount + 1
    END DO   
  END IF

  ! ON WORKER: LOAD ECOUNT AND KCOUNT FROM BUFFER 
  IF((proc_id .NE. 0) .AND. (loop .EQ. 3))THEN
    ecount = bint(2)  
    kcount = bint(3)        
  END IF

  ! ON WORKER: SET KEY POSITIONS
  IF((proc_id .NE. 0) .AND. (loop .EQ. 3))THEN
    kpos = 3 + ecount
  END IF

  ! LOOP OVER ELEMENTS IN ARR
  print *,"SIZE ", size(arr, 1) 
      dppos = 1
    lpos = 1
  DO n = 1, size(arr, 1)
    IF(loop .EQ. 2)THEN
      ecounter = ecounter + 1
    END IF
!############################################################################################
!############################################################################################

          CALL pup_b(arr(n)%counter, loop) 
          CALL pup_b(arr(n)%proc_id, loop) 
          CALL pup_b(arr(n)%proc_count, loop) 
          CALL pup_b(arr(n)%e_on, loop) 
          CALL pup_b(arr(n)%f_on, loop) 
          CALL pup_b(arr(n)%s_on, loop) 
          CALL pup_b(arr(n)%density, loop) 
          CALL pup_b(arr(n)%energy, loop) 
          CALL pup_b(arr(n)%forces, loop) 
          CALL pup_b(arr(n)%stress, loop) 
          CALL pup_b(arr(n)%energy_pair, loop) 
          CALL pup_b(arr(n)%energy_embed, loop) 


!############################################################################################
!############################################################################################
  END DO




  ! ON ROOT: ALLOCATE BUFFERS
  IF((proc_id .EQ. 0) .AND. (loop .EQ. 1))THEN
    CALL deallocate_arrays()
    kpos = 3 + ecount
    ipos = 3 + ecount + kcount    
    dppos = 1
    lpos = 1
    intoffset = 3 + ecount + kcount
    bheader(1) = intoffset + intcount
    bheader(2) = dpcount
    bheader(3) = lcount
    ALLOCATE(bint(1:bheader(1)))
    ALLOCATE(bdp(1:bheader(2)))
    ALLOCATE(blogical(1:bheader(3)))
    bint = 0
    bdp = 0.0D0
    blogical = .TRUE.    
    print *, bheader(:)   
  END IF

  ! SEND BUFFER HEADER TO ROOT      
  IF(loop .EQ. 1)THEN
    CALL MPI_BCAST(bheader, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, error) 
  END IF

  ! ON WORKER: ALLOCATE BUFFERS
  IF((proc_id .NE. 0) .AND. (loop .EQ. 1))THEN
    print *, bheader(:)
    CALL deallocate_arrays()
    ALLOCATE(bint(1:bheader(1)))
    ALLOCATE(bdp(1:bheader(2)))
    ALLOCATE(blogical(1:bheader(3)))
    bint = 0
    bdp = 0.0D0
    blogical = .TRUE.
  END IF


  ! SEND BUFFERS TO ROOT      
  IF (loop .EQ. 2) THEN  
    CALL MPI_BCAST(bint, bheader(1), MPI_INTEGER, 0, MPI_COMM_WORLD, error) 
    CALL MPI_BCAST(bdp, bheader(2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error) 
    CALL MPI_BCAST(blogical, bheader(3), MPI_LOGICAL, 0, MPI_COMM_WORLD, error) 
  END IF

END DO

END SUBROUTINE run_broadcast_efs




SUBROUTINE pup_g_int_1(int_in, loop)
!###########################################################################
INTEGER(kind=StandardInteger), INTENT(INOUT) :: int_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  ! COUNT
  kcount = kcount + 2
  intcount = intcount + 1
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = ipos
  bint(kpos + 1) = ipos
  kpos = kpos + 2
  ! STORE DATA
  bint(ipos) = int_in
  ipos = ipos + 1
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS  
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  kpos = kpos + 2
  ! UNPACK
  int_in = bint(a)
END IF
END SUBROUTINE pup_g_int_1


SUBROUTINE pup_g_int_2(int_in, loop)
!###########################################################################
INTEGER(kind=StandardInteger), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: int_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b, rows
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 3
  intcount = intcount + SIZE(int_in, 1)
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = ipos                                    ! Key start
  bint(kpos + 1) = -1 * (ipos + SIZE(int_in, 1) - 1)   ! Key end
  bint(kpos + 2) = SIZE(int_in, 1)                     ! Dim size
  kpos = kpos + 3
  ! STORE DATA
  bint(ipos:ipos + SIZE(int_in, 1) - 1) = int_in(:)
  ipos = ipos + SIZE(int_in, 1)
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  kpos = kpos + 3
  ! ALLOCATE
  IF (ALLOCATED(int_in)) THEN
    DEALLOCATE(int_in)
  END IF
  ALLOCATE(int_in(1:rows))
  ! UNPACK
  int_in(:) = bint(a:b)  
END IF

END SUBROUTINE pup_g_int_2


SUBROUTINE pup_g_int_3(int_in, loop)
!###########################################################################
INTEGER(kind=StandardInteger), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: int_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: a, b, rows, cols
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 4
  intcount = intcount + SIZE(int_in, 1) * SIZE(int_in, 2)
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = -1 * ipos                                    ! Key start
  bint(kpos + 1) = (ipos + SIZE(int_in, 1) * SIZE(int_in, 2) - 1)   ! Key end
  bint(kpos + 2) = SIZE(int_in, 1)                     ! Dim size
  bint(kpos + 3) = SIZE(int_in, 2)                     ! Dim size
  kpos = kpos + 4
  ! STORE DATA
  DO n = 1, SIZE(int_in, 2)
    bint(ipos:ipos + SIZE(int_in, 1) - 1) = int_in(:, n)
    ipos = ipos + SIZE(int_in, 1)
  END DO  
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  cols = bint(kpos+3)
  kpos = kpos + 4
  ! ALLOCATE
  IF (ALLOCATED(int_in)) THEN
    DEALLOCATE(int_in)
  END IF
  ALLOCATE(int_in(1:rows, 1:cols))
  ! UNPACK
  DO n = 1,cols
    int_in(:,n) = bint(a:a+rows-1)
    a = a + rows
  END DO
END IF

END SUBROUTINE pup_g_int_3




SUBROUTINE pup_g_dp_1(dp_in, loop)
!###########################################################################
REAL(kind=DoubleReal), INTENT(INOUT) :: dp_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 2
  dpcount = dpcount + 1
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = dppos
  bint(kpos + 1) = dppos
  kpos = kpos + 2
  ! STORE DATA  
  bdp(dppos) = dp_in
  dppos = dppos + 1
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS  
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  kpos = kpos + 2
  ! UNPACK
  dp_in = bdp(a)  
END IF
END SUBROUTINE pup_g_dp_1


SUBROUTINE pup_g_dp_2(dp_in, loop)
!###########################################################################
REAL(kind=DoubleReal), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: dp_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b, rows
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 3
  dpcount = dpcount + SIZE(dp_in, 1)
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = dppos                                    ! Key start
  bint(kpos + 1) = -1 * (dppos + SIZE(dp_in, 1) - 1)   ! Key end
  bint(kpos + 2) = SIZE(dp_in, 1)                     ! Dim size
  kpos = kpos + 3
  ! STORE DATA
  bdp(dppos:dppos + SIZE(dp_in, 1) - 1) = dp_in(:)
  dppos = dppos + SIZE(dp_in, 1)  
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  kpos = kpos + 3
  ! ALLOCATE
  IF (ALLOCATED(dp_in)) THEN
    DEALLOCATE(dp_in)
  END IF
  ALLOCATE(dp_in(1:rows))
  ! UNPACK
  dp_in(1:rows) = bdp(a:b)  
END IF
END SUBROUTINE pup_g_dp_2


SUBROUTINE pup_g_dp_3(dp_in, loop)
!###########################################################################
REAL(kind=DoubleReal), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: dp_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: a, b, rows, cols
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 4
  dpcount = dpcount + SIZE(dp_in, 1) * SIZE(dp_in, 2)
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = -1 * dppos                                    ! Key start
  bint(kpos + 1) = (dppos + SIZE(dp_in, 1) * SIZE(dp_in, 2) - 1)   ! Key end
  bint(kpos + 2) = SIZE(dp_in, 1)                     ! Dim size
  bint(kpos + 3) = SIZE(dp_in, 2)                     ! Dim size
  kpos = kpos + 4
  ! STORE DATA
  DO n = 1, SIZE(dp_in, 1)
    bdp(dppos:dppos + SIZE(dp_in, 2) - 1) = dp_in(n, :)
    dppos = dppos + SIZE(dp_in, 2)
  END DO   
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  cols = bint(kpos+3)
  kpos = kpos + 4
  ! ALLOCATE
  IF (ALLOCATED(dp_in)) THEN
    DEALLOCATE(dp_in)
  END IF
  ALLOCATE(dp_in(1:rows, 1:cols))
  ! UNPACK
  DO n = 1,cols
    dp_in(:,n) = bdp(a:a+rows-1)
    a = a + rows
  END DO
END IF
END SUBROUTINE pup_g_dp_3





SUBROUTINE pup_g_l_1(l_in, loop)
!###########################################################################
LOGICAL, INTENT(INOUT) :: l_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 2
  lcount = lcount + 1
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = lpos
  bint(kpos + 1) = lpos
  kpos = kpos + 2
  ! STORE DATA  
  blogical(lpos) = l_in
  lpos = lpos + 1
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS  
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  kpos = kpos + 2
  ! UNPACK
  l_in = blogical(a)  
END IF
END SUBROUTINE pup_g_l_1


SUBROUTINE pup_g_l_2(l_in, loop)
!###########################################################################
LOGICAL, INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: l_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b, rows
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 3
  lcount = lcount + SIZE(l_in, 1)
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = lpos                                    ! Key start
  bint(kpos + 1) = -1 * (lpos + SIZE(l_in, 1) - 1)   ! Key end
  bint(kpos + 2) = SIZE(l_in, 1)                     ! Dim size
  kpos = kpos + 3
  ! STORE DATA
  blogical(lpos:lpos + SIZE(l_in, 1) - 1) = l_in(:)
  lpos = lpos + SIZE(l_in, 1) 
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  kpos = kpos + 3
  ! ALLOCATE
  IF (ALLOCATED(l_in)) THEN
    DEALLOCATE(l_in)
  END IF
  ALLOCATE(l_in(1:rows))
  ! UNPACK
  l_in(:) = blogical(a:b)  
END IF
END SUBROUTINE pup_g_l_2


SUBROUTINE pup_g_l_3(l_in, loop)
!###########################################################################
LOGICAL, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: l_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: a, b, rows, cols
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .NE. 0) THEN
  kcount = kcount + 4
  lcount = lcount + SIZE(l_in, 1) * SIZE(l_in, 2)
ELSE IF (loop .eq. 2 .AND. proc_id .NE. 0) THEN
  ! STORE KEY
  bint(kpos) = -1 * lpos                                    ! Key start
  bint(kpos + 1) = (lpos + SIZE(l_in, 1) * SIZE(l_in, 2) - 1)   ! Key end
  bint(kpos + 2) = SIZE(l_in, 1)                     ! Dim size
  bint(kpos + 3) = SIZE(l_in, 2)                     ! Dim size
  kpos = kpos + 4
  ! STORE DATA
  DO n = 1, SIZE(l_in, 1)
    blogical(lpos:lpos + SIZE(l_in, 2) - 1) = l_in(n, :)
    lpos = lpos + SIZE(l_in, 2)
  END DO   
ELSE IF (loop .eq. 3 .AND. proc_id .EQ. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  cols = bint(kpos+3)
  kpos = kpos + 4
  ! ALLOCATE
  IF (ALLOCATED(l_in)) THEN
    DEALLOCATE(l_in)
  END IF
  ALLOCATE(l_in(1:rows, 1:cols))
  ! UNPACK
  DO n = 1,cols
    l_in(:,n) = blogical(a:a+rows-1)
    a = a + rows
  END DO
END IF

END SUBROUTINE pup_g_l_3





SUBROUTINE pup_b_int_1(int_in, loop)
!###########################################################################
INTEGER(kind=StandardInteger), INTENT(INOUT) :: int_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  ! COUNT
  kcount = kcount + 2
  intcount = intcount + 1
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = ipos
  bint(kpos + 1) = ipos
  kpos = kpos + 2
  ! STORE DATA
  bint(ipos) = int_in
  ipos = ipos + 1
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS  
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  kpos = kpos + 2
  ! UNPACK
  int_in = bint(a)
END IF
END SUBROUTINE pup_b_int_1


SUBROUTINE pup_b_int_2(int_in, loop)
!###########################################################################
INTEGER(kind=StandardInteger), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: int_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b, rows
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 3
  intcount = intcount + SIZE(int_in, 1)
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = ipos                                    ! Key start
  bint(kpos + 1) = -1 * (ipos + SIZE(int_in, 1) - 1)   ! Key end
  bint(kpos + 2) = SIZE(int_in, 1)                     ! Dim size
  kpos = kpos + 3
  ! STORE DATA
  bint(ipos:ipos + SIZE(int_in, 1) - 1) = int_in(:)
  ipos = ipos + SIZE(int_in, 1)
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  kpos = kpos + 3
  ! ALLOCATE
  IF (ALLOCATED(int_in)) THEN
    DEALLOCATE(int_in)
  END IF
  ALLOCATE(int_in(1:rows))
  ! UNPACK
  int_in(:) = bint(a:b)  
END IF

END SUBROUTINE pup_b_int_2


SUBROUTINE pup_b_int_3(int_in, loop)
!###########################################################################
INTEGER(kind=StandardInteger), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: int_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: a, b, rows, cols
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 4
  intcount = intcount + SIZE(int_in, 1) * SIZE(int_in, 2)
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = -1 * ipos                                    ! Key start
  bint(kpos + 1) = (ipos + SIZE(int_in, 1) * SIZE(int_in, 2) - 1)   ! Key end
  bint(kpos + 2) = SIZE(int_in, 1)                     ! Dim size
  bint(kpos + 3) = SIZE(int_in, 2)                     ! Dim size
  kpos = kpos + 4
  ! STORE DATA
  DO n = 1, SIZE(int_in, 2)
    bint(ipos:ipos + SIZE(int_in, 1) - 1) = int_in(:, n)
    ipos = ipos + SIZE(int_in, 1)
  END DO  
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  cols = bint(kpos+3)
  kpos = kpos + 4
  ! ALLOCATE
  IF (ALLOCATED(int_in)) THEN
    DEALLOCATE(int_in)
  END IF
  ALLOCATE(int_in(1:rows, 1:cols))
  ! UNPACK
  DO n = 1,cols
    int_in(:,n) = bint(a:a+rows-1)
    a = a + rows
  END DO
END IF

END SUBROUTINE pup_b_int_3




SUBROUTINE pup_b_dp_1(dp_in, loop)
!###########################################################################
REAL(kind=DoubleReal), INTENT(INOUT) :: dp_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 2
  dpcount = dpcount + 1
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = dppos
  bint(kpos + 1) = dppos
  kpos = kpos + 2
  ! STORE DATA  
  bdp(dppos) = dp_in
  dppos = dppos + 1
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS  
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  kpos = kpos + 2
  ! UNPACK
  dp_in = bdp(a)  
END IF
END SUBROUTINE pup_b_dp_1


SUBROUTINE pup_b_dp_2(dp_in, loop)
!###########################################################################
REAL(kind=DoubleReal), INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: dp_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b, rows
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 3
  dpcount = dpcount + SIZE(dp_in, 1)
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = dppos                                    ! Key start
  bint(kpos + 1) = -1 * (dppos + SIZE(dp_in, 1) - 1)   ! Key end
  bint(kpos + 2) = SIZE(dp_in, 1)                     ! Dim size
  kpos = kpos + 3
  ! STORE DATA
  bdp(dppos:dppos + SIZE(dp_in, 1) - 1) = dp_in(:)
  dppos = dppos + SIZE(dp_in, 1)  
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  kpos = kpos + 3
  ! ALLOCATE
  IF (ALLOCATED(dp_in)) THEN
    DEALLOCATE(dp_in)
  END IF
  ALLOCATE(dp_in(1:rows))
  ! UNPACK
  dp_in(1:rows) = bdp(a:b)  
END IF
END SUBROUTINE pup_b_dp_2


SUBROUTINE pup_b_dp_3(dp_in, loop)
!###########################################################################
REAL(kind=DoubleReal), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: dp_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: a, b, rows, cols
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 4
  dpcount = dpcount + SIZE(dp_in, 1) * SIZE(dp_in, 2)
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = -1 * dppos                                    ! Key start
  bint(kpos + 1) = (dppos + SIZE(dp_in, 1) * SIZE(dp_in, 2) - 1)   ! Key end
  bint(kpos + 2) = SIZE(dp_in, 1)                     ! Dim size
  bint(kpos + 3) = SIZE(dp_in, 2)                     ! Dim size
  kpos = kpos + 4
  ! STORE DATA
  DO n = 1, SIZE(dp_in, 1)
    bdp(dppos:dppos + SIZE(dp_in, 2) - 1) = dp_in(n, :)
    dppos = dppos + SIZE(dp_in, 2)
  END DO   
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  cols = bint(kpos+3)
  kpos = kpos + 4
  ! ALLOCATE
  IF (ALLOCATED(dp_in)) THEN
    DEALLOCATE(dp_in)
  END IF
  ALLOCATE(dp_in(1:rows, 1:cols))
  ! UNPACK
  DO n = 1,cols
    dp_in(:,n) = bdp(a:a+rows-1)
    a = a + rows
  END DO
END IF
END SUBROUTINE pup_b_dp_3





SUBROUTINE pup_b_l_1(l_in, loop)
!###########################################################################
LOGICAL, INTENT(INOUT) :: l_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 2
  lcount = lcount + 1
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = lpos
  bint(kpos + 1) = lpos
  kpos = kpos + 2
  ! STORE DATA  
  blogical(lpos) = l_in
  lpos = lpos + 1
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS  
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  kpos = kpos + 2
  ! UNPACK
  l_in = blogical(a)  
END IF
END SUBROUTINE pup_b_l_1


SUBROUTINE pup_b_l_2(l_in, loop)
!###########################################################################
LOGICAL, INTENT(INOUT), ALLOCATABLE, DIMENSION(:) :: l_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: a, b, rows
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 3
  lcount = lcount + SIZE(l_in, 1)
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = lpos                                    ! Key start
  bint(kpos + 1) = -1 * (lpos + SIZE(l_in, 1) - 1)   ! Key end
  bint(kpos + 2) = SIZE(l_in, 1)                     ! Dim size
  kpos = kpos + 3
  ! STORE DATA
  blogical(lpos:lpos + SIZE(l_in, 1) - 1) = l_in(:)
  lpos = lpos + SIZE(l_in, 1) 
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  kpos = kpos + 3
  ! ALLOCATE
  IF (ALLOCATED(l_in)) THEN
    DEALLOCATE(l_in)
  END IF
  ALLOCATE(l_in(1:rows))
  ! UNPACK
  l_in(:) = blogical(a:b)  
END IF
END SUBROUTINE pup_b_l_2


SUBROUTINE pup_b_l_3(l_in, loop)
!###########################################################################
LOGICAL, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: l_in
INTEGER(kind=StandardInteger) :: loop
INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: a, b, rows, cols
!###########################################################################
IF (loop .eq. 1 .AND. proc_id .EQ. 0) THEN
  kcount = kcount + 4
  lcount = lcount + SIZE(l_in, 1) * SIZE(l_in, 2)
ELSE IF (loop .eq. 2 .AND. proc_id .EQ. 0) THEN
  ! STORE KEY
  bint(kpos) = -1 * lpos                                    ! Key start
  bint(kpos + 1) = (lpos + SIZE(l_in, 1) * SIZE(l_in, 2) - 1)   ! Key end
  bint(kpos + 2) = SIZE(l_in, 1)                     ! Dim size
  bint(kpos + 3) = SIZE(l_in, 2)                     ! Dim size
  kpos = kpos + 4
  ! STORE DATA
  DO n = 1, SIZE(l_in, 1)
    blogical(lpos:lpos + SIZE(l_in, 2) - 1) = l_in(n, :)
    lpos = lpos + SIZE(l_in, 2)
  END DO   
ELSE IF (loop .eq. 3 .AND. proc_id .NE. 0) THEN  
  ! KEYS
  a = abs(bint(kpos))
  b = abs(bint(kpos+1))
  rows = bint(kpos+2)
  cols = bint(kpos+3)
  kpos = kpos + 4
  ! ALLOCATE
  IF (ALLOCATED(l_in)) THEN
    DEALLOCATE(l_in)
  END IF
  ALLOCATE(l_in(1:rows, 1:cols))
  ! UNPACK
  DO n = 1,cols
    l_in(:,n) = blogical(a:a+rows-1)
    a = a + rows
  END DO
END IF

END SUBROUTINE pup_b_l_3





SUBROUTINE deallocate_arrays()
!###########################################################################
IF (ALLOCATED(bint)) THEN
  DEALLOCATE(bint)
END IF
IF (ALLOCATED(bdp)) THEN
  DEALLOCATE(bdp)
END IF
IF (ALLOCATED(blogical)) THEN
  DEALLOCATE(blogical)
END IF
END SUBROUTINE deallocate_arrays

END MODULE mpi_udgb



