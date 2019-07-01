mkdir -p mod
mpif90 -g \
-fbounds-check \
-mtune=native \
-fopenmp \
-O3 \
-ffast-math \
kinds.f90 \
shape_type.f90 \
store_shape.f90 \
mod.f90 \
calc.f90 \
test.f90 \
-J mod \
-o test_mpi.x \
&& export OMP_NUM_THREADS=2 \
&& mpirun -n 4 test_mpi.x
