# README for Mesh Processing Codes

[Reference](https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/api/gmsh.f90?ref_type=heads)

## Overview

This repository contains Fortran codes for processing 2D mesh data using the GMSH API. The codes can:

1. Read GMSH mesh files (.msh)
2. Extract physical group information (boundaries and surfaces)
3. Retrieve node coordinates and connectivity data
4. Write mesh data to text files for further processing

## Files

- `donut2d_mesh.geo`: GMSH script to generate a 2D annular mesh (donut shape)
- `matrix_writer.f90`: Module for writing matrices and vectors to files
- `mesh_module.f90`: Main module for mesh processing using GMSH API
- `test_mesh_module.f90`: Test program demonstrating the functionality

## Compilation Instructions

To compile the code, run:

```bash
gfortran -o test_mesh_module /usr/local/include/gmsh.f90 matrix_writer.f90 mesh_module.f90 test_mesh_module.f90 -L/usr/local/lib -lgmsh
```

## Usage

1. First generate a mesh file using GMSH:
   ```bash
   gmsh donut2d_mesh.geo -2 -o donut2d_mesh.msh
   ```

2. Run the compiled program:
   ```bash
   ./test_mesh_module
   ```

The program will:
- Extract physical group information (inner and outer boundaries, surface)
- Retrieve all node coordinates and connectivity data
- Write coordinate and connectivity data to `coordAll.txt` and `connectivity.txt`

## Output Files

- `coordAll.txt`: Contains all node coordinates (x,y,z)
- `connectivity.txt`: Contains element connectivity information

## Dependencies

- GMSH (including Fortran API)
- gfortran compiler

## Notes

- The code is specifically tested with a 2D annular mesh but can be adapted for other 2D meshes
- Ensure GMSH is properly installed with Fortran API support
- The paths in the compilation command may need adjustment based on your GMSH installation
