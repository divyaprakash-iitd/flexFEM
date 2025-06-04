# #gfortran -c mod_solid.f90 -o mod_solid.o
# #gfortran -c mod_io.f90 -o mod_io.o
# #gfortran -c mod_solid.o fem2d.f90 -o fem2d.o
# #gfortran -c -fPIC mod_io.o fem2d.o soft_particles.f90 -o soft_particles.o
# #gfortran -shared -o libsoftparticles.so soft_particles.o
# 
# # Compile individual modules
# gfortran -c mod_solid.f90 -o mod_solid.o
# gfortran -c mod_io.f90 -o mod_io.o
# 
# # Compile fem2d.f90 and mod_solid.f90 into object files
# gfortran -c fem2d.f90 -o fem2d.o
# 
# # Compile soft_particles.f90 and link with object files to create shared library
# gfortran -c -fPIC soft_particles.f90 -o soft_particles.o
# gfortran -shared -o libsoftparticles.so soft_particles.o fem2d.o mod_solid.o mod_io.o
# 
# gfortran ibmc.f90 -L. -lsoftparticles -o ibmc

# New instuctions
gfortran -c -fPIC mod_solid.f90 -o mod_solid.o
gfortran -c -fPIC mod_io.f90 -o mod_io.o
gfortran -c -fPIC fem2d.f90 -o fem2d.o
gfortran -c -fPIC soft_particles.f90 -o soft_particles.o
gfortran -shared -o libsoftparticles.so *.o
gfortran ibmc.f90 -L. -lsoftparticles -o ibmc

# Set library path (if needed)
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
./ibmc
