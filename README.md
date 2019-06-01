A python3 script that creates a Fortran module to package up a user defined data type array and:

1. collect parts of array from other processes back to the root process where mod(n-1,process_count) = process_id
2. broadcast entire array from root to all processes

It can handle integers, doubles, logicals as scalars, 1D arrays and 2D arrays.  All arrays must be allocatable.




There's an example program showing how I've used the module.  I still need to test is thoroughly.
