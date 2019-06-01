MODULE calc

USE kinds
USE shape_type, ONLY: tshape
USE mpi_udgb, ONLY: run_gather, run_broadcast
USE store_shape, ONLY: store_before, store_after

IMPLICIT NONE

! Include MPI header
INCLUDE 'mpif.h'


CONTAINS



SUBROUTINE run_calc(shapes)
!##################################################
TYPE(tshape), ALLOCATABLE, DIMENSION(:) :: shapes
INTEGER(kind=StandardInteger) :: proc_id, proc_count, error
INTEGER(kind=StandardInteger), DIMENSION(MPI_STATUS_SIZE) :: status
!##################################################
INTEGER(kind=StandardInteger) :: n
!##################################################

! Get MPI details
CALL MPI_Comm_rank(MPI_COMM_WORLD,proc_id,error)
CALL MPI_Comm_size(MPI_COMM_WORLD,proc_count,error)

! Loop through 
DO n = 1,12
  IF (mod(n-1,proc_count) .EQ. proc_id) THEN
    CALL shape_calc(shapes(n))
  END IF
END DO



CALL store_before(shapes)
CALL run_gather(shapes)
CALL run_broadcast(shapes)
CALL store_after(shapes)


END SUBROUTINE run_calc





SUBROUTINE shape_calc(shape_inout)
!##################################################
TYPE(tshape), INTENT(INOUT) :: shape_inout
!##################################################
INTEGER(kind=StandardInteger) :: n, rows, cols
REAL(kind=DoubleReal) :: r
!##################################################
CALL RANDOM_NUMBER(r)

shape_inout%sides = floor(9 * r) + 1

ALLOCATE(shape_inout%lengths(1:shape_inout%sides))

shape_inout%perimeter = 0.0D0
DO n = 1,shape_inout%sides
  CALL RANDOM_NUMBER(r)
  shape_inout%lengths(n) = r
  shape_inout%perimeter = shape_inout%perimeter + r 
END DO


CALL RANDOM_NUMBER(r)
rows = 2 + floor(7 * r) + 1
ALLOCATE(shape_inout%int1d(1:rows))
shape_inout%int1d = 111

CALL RANDOM_NUMBER(r)
rows = 2 + floor(7 * r) + 1
CALL RANDOM_NUMBER(r)
cols = 3
ALLOCATE(shape_inout%int2d(1:rows, 1:cols))
shape_inout%int2d = 222



END SUBROUTINE shape_calc




END MODULE calc