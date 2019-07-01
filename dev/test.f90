PROGRAM test
! University of Birmingham
! Ben Palmer
USE kinds
USE calc, ONLY: run_calc
USE shape_type, ONLY: tshape
USE mpi

IMPLICIT NONE

CALL main()

CONTAINS

! Subroutines


SUBROUTINE main()
!###########################################
! PRIVATE VARIABLES
INTEGER(kind=StandardInteger) :: error
INTEGER(kind=StandardInteger) :: proc_id, proc_count
INTEGER(kind=StandardInteger) :: thread_id, thread_count
INTEGER(kind=StandardInteger) :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
INTEGER(kind=StandardInteger) :: test_int

INTEGER(kind=StandardInteger) :: n, m
REAL(kind=DoubleReal) :: r
INTEGER(kind=StandardInteger), DIMENSION(1:12) :: seed
INTEGER(kind=StandardInteger), DIMENSION(MPI_STATUS_SIZE) :: status


TYPE(tshape), ALLOCATABLE, DIMENSION(:) :: shapes

!###########################################
CALL MPI_Init(error)

Call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,error)
Call MPI_Comm_size(MPI_COMM_WORLD,proc_count,error)

! Randomise seed on all processes
seed(:) = proc_id
CALL RANDOM_SEED(put=seed)

! Create array on all processes
ALLOCATE(shapes(1:12))

! Run calculation
CALL run_calc(shapes)




!print *, proc_id, results(1:4)

CALL MPI_Finalize(error)
END SUBROUTINE main



END PROGRAM test
