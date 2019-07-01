MODULE store_shape

USE kinds
USE shape_type, ONLY: tshape
USE mpi

IMPLICIT NONE

CONTAINS



SUBROUTINE store_before(shapes)
TYPE(tshape), ALLOCATABLE, DIMENSION(:) :: shapes
INTEGER(kind=StandardInteger) :: fileid, proc_id, proc_count, error
INTEGER(kind=StandardInteger), DIMENSION(MPI_STATUS_SIZE) :: status
INTEGER(kind=StandardInteger) :: n, m

Call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,error)
Call MPI_Comm_size(MPI_COMM_WORLD,proc_count,error)

IF(proc_id .EQ. 0)THEN
  fileid = 100
  OPEN(UNIT=fileid,FILE='out_before_0.txt')
ELSE IF(proc_id .EQ. 1)THEN
  fileid = 101
  OPEN(UNIT=fileid,FILE='out_before_1.txt')
ELSE IF(proc_id .EQ. 2)THEN
  fileid = 102
  OPEN(UNIT=fileid,FILE='out_before_2.txt')
ELSE IF(proc_id .EQ. 3)THEN
  fileid = 103
  OPEN(UNIT=fileid,FILE='out_before_3.txt')
END IF

DO n=1, size(shapes,1)
  WRITE(fileid,*) "========================================"
  WRITE(fileid,*) "      SHAPE ", n
  WRITE(fileid,*) "========================================"
  WRITE(fileid,*) "sides: ", shapes(n)%sides
  IF(ALLOCATED(shapes(n)%lengths))THEN
    DO m=1, size(shapes(n)%lengths,1)
      WRITE(fileid,*) "len ", m,": ", shapes(n)%lengths(m)
    END DO
  END IF
  WRITE(fileid,*) "perimeter: ", shapes(n)%perimeter
  IF(ALLOCATED(shapes(n)%coords1))THEN
    WRITE(fileid,*) "coords1_1:    ", shapes(n)%coords1(1,:)
    WRITE(fileid,*) "coords1_1000: ", shapes(n)%coords1(1000,:)    
  END IF
  IF(ALLOCATED(shapes(n)%coords2))THEN
    WRITE(fileid,*) "coords2_1:    ", shapes(n)%coords2(1,:)
    WRITE(fileid,*) "coords2_1000: ", shapes(n)%coords2(1000,:)    
  END IF
  IF(ALLOCATED(shapes(n)%enddp))THEN
    WRITE(fileid,*) "enddp:    ", shapes(n)%enddp(1,1)
  END IF
  
  WRITE(fileid,*) ""
  WRITE(fileid,*) ""
END DO


CLOSE(fileid)



END SUBROUTINE store_before



SUBROUTINE store_after(shapes)
TYPE(tshape), ALLOCATABLE, DIMENSION(:) :: shapes
INTEGER(kind=StandardInteger) :: fileid, proc_id, proc_count, error
INTEGER(kind=StandardInteger), DIMENSION(MPI_STATUS_SIZE) :: status
INTEGER(kind=StandardInteger) :: n, m

Call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,error)
Call MPI_Comm_size(MPI_COMM_WORLD,proc_count,error)

IF(proc_id .EQ. 0)THEN
  fileid = 100
  OPEN(UNIT=fileid,FILE='out_after_0.txt')
ELSE IF(proc_id .EQ. 1)THEN
  fileid = 101
  OPEN(UNIT=fileid,FILE='out_after_1.txt')
ELSE IF(proc_id .EQ. 2)THEN
  fileid = 102
  OPEN(UNIT=fileid,FILE='out_after_2.txt')
ELSE IF(proc_id .EQ. 3)THEN
  fileid = 103
  OPEN(UNIT=fileid,FILE='out_after_3.txt')
END IF

DO n=1, size(shapes,1)
  WRITE(fileid,*) "========================================"
  WRITE(fileid,*) "      SHAPE ", n
  WRITE(fileid,*) "========================================"
  WRITE(fileid,*) "sides: ", shapes(n)%sides
  IF(ALLOCATED(shapes(n)%lengths))THEN
    DO m=1, size(shapes(n)%lengths,1)
      WRITE(fileid,*) "len ", m,": ", shapes(n)%lengths(m)
    END DO
  END IF
  WRITE(fileid,*) "perimeter: ", shapes(n)%perimeter
  IF(ALLOCATED(shapes(n)%coords1))THEN
    WRITE(fileid,*) "coords1_1:    ", shapes(n)%coords1(1,:)
    WRITE(fileid,*) "coords1_1000: ", shapes(n)%coords1(1000,:)    
  END IF
  IF(ALLOCATED(shapes(n)%coords2))THEN
    WRITE(fileid,*) "coords2_1:    ", shapes(n)%coords2(1,:)
    WRITE(fileid,*) "coords2_1000: ", shapes(n)%coords2(1000,:)    
  END IF
  IF(ALLOCATED(shapes(n)%enddp))THEN
    WRITE(fileid,*) "enddp:    ", shapes(n)%enddp(1,1)
  END IF
  
  WRITE(fileid,*) ""
  WRITE(fileid,*) ""
END DO


CLOSE(fileid)


END SUBROUTINE store_after

END MODULE store_shape
