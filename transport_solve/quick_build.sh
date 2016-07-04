gfortran -c transport_solve.f95 matrix_tools.f95 transport_tools.f95
gfortran transport_solve.o matrix_tools.o transport_tools.o