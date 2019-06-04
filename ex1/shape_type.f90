MODULE shape_type

USE kinds

IMPLICIT NONE

TYPE tshape
INTEGER(kind=StandardInteger) :: sides = 0
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:) :: lengths
INTEGER(kind=StandardInteger), ALLOCATABLE, DIMENSION(:) :: int1d
INTEGER(kind=StandardInteger), ALLOCATABLE, DIMENSION(:,:) :: int2d
REAL(kind=DoubleReal) :: perimeter
LOGICAL :: calculated = .FALSE.
END TYPE tshape

END MODULE shape_type
