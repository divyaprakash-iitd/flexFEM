# Mesh Handling Module

This directory contains Fortran modules and programs for handling mesh data, including reading mesh files, extracting physical group information, and writing matrix/vector data to files.

## Directory Structure

- **`test_mesh_module.f90`**: A test program that demonstrates the usage of `mesh_module` and `matrix_writer` modules.
- **`mesh_module.f90`**: Contains subroutines for reading mesh files and extracting node and physical group data using GMSH.
- **`matrix_writer.f90`**: Provides subroutines for writing vectors and matrices to text files.
- **`donut2d_mesh.msh`**: A sample GMSH mesh file used for testing.
- **`coordAll.txt`**: Output file containing node coordinates written by the test program.
- **`README.md`**: Instructions for compiling and understanding the directory contents.
- **`mesh_module.mod`**, **`matrix_writer.mod`**, **`gmsh.mod`**: Compiled module files generated during the build process.

## Compile Instructions

To compile the test program, use the following command:

```bash
gfortran -I/usr/local/include -o test_mesh_module gmsh.f90 matrix_writer.f90 mesh_module.f90 test_mesh_module.f90 -L/usr/local/lib -lgmsh
```

## Usage

1. Compile the program using the instructions above.
2. Run the `test_mesh_module` executable to test mesh handling functionality.
3. Check the output files (e.g., `coordAll.txt`) for results.

## Dependencies

- **GMSH**: The program uses GMSH's API for mesh handling. Ensure GMSH is installed and its libraries are accessible.

## Notes

- The `test_mesh_module.f90` program demonstrates reading mesh data from `donut2d_mesh.msh` and writing node coordinates to `coordAll.txt`.
- Modify the test program or mesh file as needed for your use case.
```