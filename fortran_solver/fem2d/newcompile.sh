gfortran -O2 -c mod_solid.f90
gfortran -O2 -c mod_io.f90
gfortran -O2 -c fem2d.f90
gfortran -O2 -c soft_particles.f90
gfortran -O2 -o ibmc ibmc.f90 mod_solid.o mod_io.o fem2d.o soft_particles.o
