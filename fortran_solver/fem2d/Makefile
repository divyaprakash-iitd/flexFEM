# Makefile for flexFEM project
# Compiler and flags
FC = gfortran
FFLAGS = -O2
LDFLAGS = -L/usr/local/lib
LIBS = -lgmsh

# Directories
GMSH_INCLUDE = /usr/local/include

# Target executable
TARGET = main

# Object files
OBJS = main.o \
       fem2d.o \
       soft_particles.o \
       gmsh.o \
       matrix_writer.o \
       mesh_module.o

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# Compile local Fortran source files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Compile gmsh.f90 from include directory
gmsh.o: $(GMSH_INCLUDE)/gmsh.f90
	$(FC) $(FFLAGS) -c $<

# Dependencies
matrix_writer.o: gmsh.o
mesh_module.o: gmsh.o matrix_writer.o
fem2d.o: mesh_module.o
soft_particles.o: fem2d.o mesh_module.o matrix_writer.o
main.o: soft_particles.o

# Clean up generated files
clean:
	rm -f *.o *.mod $(TARGET)

# Phony targets
.PHONY: all clean

