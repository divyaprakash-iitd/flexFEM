# Makefile for IBMC project
# Compiler and flags
FC = gfortran
FFLAGS = -O2
LDFLAGS = -L/usr/local/lib
LIBS = -lgmsh

# Directories
MESH_DIR = ../mesh_handling_fortran
GMSH_INCLUDE = /usr/local/include

# Target executable
TARGET = ibmc

# Object files
OBJS = ibmc.o \
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

# Compile files from mesh handling directory
matrix_writer.o: $(MESH_DIR)/matrix_writer.f90
	$(FC) $(FFLAGS) -c $<

mesh_module.o: $(MESH_DIR)/mesh_module.f90
	$(FC) $(FFLAGS) -c $<

# Dependencies (based on compilation order in script)
# Basic modules compiled first
gmsh.o:
matrix_writer.o:
mesh_module.o:

# soft_particles.f90 uses fem2d, mesh_module, and matrix_writer
soft_particles.o: fem2d.o mesh_module.o matrix_writer.o

# ibmc.f90 uses soft_particles module
ibmc.o: soft_particles.o

# Clean up generated files
clean:
	rm -f *.o *.mod $(TARGET)

# Phony targets
.PHONY: all clean

# Additional useful targets
rebuild: clean all

install: $(TARGET)
	@echo "Add installation commands here if needed"

help:
	@echo "Available targets:"
	@echo "  all      - Build the $(TARGET) executable (default)"
	@echo "  clean    - Remove object files, modules, and executable"
	@echo "  rebuild  - Clean and build from scratch"
	@echo "  help     - Show this help message"
