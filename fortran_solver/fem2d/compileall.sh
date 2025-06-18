# Option # 1

gfortran -O2 -c mod_solid.f90
gfortran -O2 -c mod_io.f90
gfortran -O2 -c fem2d.f90
gfortran -O2 -c /usr/local/include/gmsh.f90
gfortran -O2 -c ../mesh_handling_fortran/matrix_writer.f90
gfortran -O2 -c ../mesh_handling_fortran/mesh_module.f90
gfortran -O2 -c soft_particles.f90  # This now can use mesh_module.mod
gfortran -O2 -c ibmc.f90

gfortran -O2 -o ibmc ibmc.o mod_solid.o mod_io.o fem2d.o soft_particles.o gmsh.o matrix_writer.o mesh_module.o -L/usr/local/lib -lgmsh

# Option # 2

#gfortran -O2 -c /usr/local/include/gmsh.f90 matrix_writer.f90 mesh_module.f90
#gfortran -O2 -c mod_solid.f90 mod_io.f90 fem2d.f90 soft_particles.f90 ibmc.f90
#gfortran -O2 -o ibmc *.o -L/usr/local/lib -lgmsh
