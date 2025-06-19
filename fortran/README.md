# Implementation in Fortran

## Overview

This repository provides a Fortran-based code for Finite Element Method (FEM) simulations. For demonstartion purposes a
deformation of annular (donut-shaped) structures under applied forces is simulated. However the provided functions can be used
to simulate any type of 2D geometry.

## Requirements

To run this simulation, ensure you have the following installed:

- **GMSH library**: Compiled from source. Instructions are available [here](https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/fortran/README.txt).
- **Fortran compiler**: e.g., `gfortran`.
- **PyVista**: For generating VTK files for visualization (included in the provided `environment.yml`).
- **ParaView**: For viewing the generated VTK files.
- **Conda**: For managing the Python environment (optional but recommended).

The `environment.yml` file is provided to set up a Conda environment with PyVista and the GMSH Python bindings.

## Installation

1. **Compile GMSH from Source**:
   - Follow the [GMSH Fortran tutorial](https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/fortran/README.txt) to compile GMSH with Fortran support.
   - Note the installation path (e.g., `/path/to/gmsh`), as it will be needed for compilation.

2. **Set Up the Conda Environment** (Optional):
   - Create and activate the environment using the provided `environment.yml`:
     ```bash
     conda env create -f environment.yml
     conda activate mesh-env
     ```

3. **Compile the Fortran Code**:
   - Clone or download this repository.
   - Open the `Makefile` and ensure the `GMSH_PATH` variable matches your GMSH source installation path.
   - Run the following command to compile the executable:
     ```bash
     make
     ```
   - This generates the `main` executable in the root directory.

## Running the Simulation

1. **Generate the Mesh**:
   - Navigate to the `tutorials/donut2d` directory.
   - Use GMSH to create the mesh file from the provided geometry file:
     ```bash
     gmsh donut2d_mesh.geo -2 -o donut2d_mesh.msh
     ```
   - This generates `donut2d_mesh.msh`, the mesh file used by the simulation.

2. **Run the Simulation**:
   - Specify the simulation parameters in the `input_params.dat`.
   - Execute the simulation script:
     ```bash
     ./run.sh
     ```
    This deletes the previous simulation files, copies the executable and runs it.
   - The simulation runs for 50,000 iterations, applying forces for the first half and then relaxing the structure. Output files (`F_*.txt` for forces and `P_*.txt` for positions) are saved in the `images` directory every 200 iterations.

## Visualization

1. **Generate VTK Files**:
   - After running the simulation, use the provided Python script to convert output files to VTK format:
     ```bash
     python plot_mesh.py
     ```
   - Ensure the Conda environment (`fem2d_env`) is active if using Conda.
   - This generates VTK files in the `images` directory.

2. **View Results in ParaView**:
   - Open ParaView.
   - Load the VTK files from the `images` directory to visualize the deformation of the annular ring over time.

## Code Structure

The repository includes the following key files:

- **`fem2d.f90`**:
  - Implements the core FEM solver, including force calculations and position updates for a 2D mesh.
  - Defines the `festruct` type to store mesh data and simulation parameters.

- **`fem2d_interface.f90`**:
  - Provides a simplified interface to the FEM solver with subroutines like `calculateforces`, `getforces`, `getpositions`, and `updatepositions`.
  - Designed for C interoperability (using `bind(C)`), enabling future use as a shared library.

- **`mesh_module.f90`**:
  - Interfaces with the GMSH Fortran API to load mesh data (nodes and connectivity) from GMSH mesh files.
  - Includes subroutines like `get_nodes_connectivity` for mesh handling.

- **`matrix_writer.f90`**:
  - Utility module for writing simulation results (matrices and vectors) to text files.
  - Supports both real and integer data types.

- **`main.f90`**:
  - The main program that orchestrates the simulation.
  - Applies boundary conditions (forces on the annular ring’s periphery) and iterates through time steps.

Additional files:
- **`tutorials/donut2d/donut2d_mesh.geo`**: Geometry file for the annular ring mesh.
- **`tutorials/donut2d/run.sh`**: Script to execute the simulation.
- **`tutorials/donut2d/plot_mesh.py`**: Python script for generating VTK files.
- **`input_params.dat`**: Input file specifying simulation parameters (e.g., mesh file, stiffness `Kp`, and damping `Bp`).

## Simulation Details

- **Geometry**: A 2D annular ring with an outer radius of 1.0 unit.
- **Boundary Conditions**: Forces of magnitude 500 are applied horizontally (outward on the right, inward on the left) on the periphery between -45° and 45° relative to the y-axis for the first 25,000 iterations.
- **Time Stepping**: 50,000 iterations with a time step `dt = 0.001`.
- **Output**: Forces and positions are saved every 200 iterations in the `images` directory.

## Notes

- Ensure the GMSH library path in the `Makefile` matches your installation to avoid linking errors.
- The simulation duration depends on the number of iterations and mesh size; progress is printed every 1,000 iterations.
- Modify `input_params.dat` to adjust simulation parameters like `Kp` (stiffness) or `Bp` (material coefficient).
- If visualization fails, verify that PyVista is installed and the `images` directory contains the expected output files.
