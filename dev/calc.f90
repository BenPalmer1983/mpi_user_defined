MODULE calc

USE kinds
USE shape_type, ONLY: tshape
USE mpi_udgb, ONLY: run_gather, run_broadcast
USE store_shape, ONLY: store_before, store_after
USE mpi

IMPLICIT NONE


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

IF(proc_id .EQ. 0)THEN
print *, shapes(1)%coords1(1,:)
print *, shapes(1)%coords1(1024,:)
print *, shapes(2)%coords1(1,:)
print *, shapes(2)%coords1(1024,:)
END IF


CALL run_broadcast(shapes)
CALL store_after(shapes)




END SUBROUTINE run_calc





SUBROUTINE shape_calc(shape_inout)
!##################################################
TYPE(tshape), INTENT(INOUT) :: shape_inout
!##################################################
INTEGER(kind=StandardInteger) :: n, m, rows, cols
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



ALLOCATE(shape_inout%coords1(1:1024,1:3))
ALLOCATE(shape_inout%coords2(1:1024,1:3))
ALLOCATE(shape_inout%matrix1(1:10,1:10))
ALLOCATE(shape_inout%enddp(1:1,1:1))


DO n = 1,1024
  shape_inout%coords1(n, 1) = 3.0D0 * (n-1) + 1
  shape_inout%coords1(n, 2) = 3.0D0 * (n-1) + 2
  shape_inout%coords1(n, 3) = 3.0D0 * (n-1) + 3
END DO

DO n = 1,1024
  shape_inout%coords2(n, 1) = 10000.0D0 + 3.0D0 * (n-1) + 1
  shape_inout%coords2(n, 2) = 10000.0D0 + 3.0D0 * (n-1) + 2
  shape_inout%coords2(n, 3) = 10000.0D0 + 3.0D0 * (n-1) + 3
END DO

DO n = 1,10
  DO m = 1,10
    shape_inout%matrix1(n, m) = n * m
  END DO
END DO

shape_inout%enddp = -5.0D0

END SUBROUTINE shape_calc




END MODULE calc
