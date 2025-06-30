# Implementation in MATLAB (3D)

## Overview

This repository provides a MATLAB-based library for 3D Finite Element Method (FEM) simulations. The demonstration simulates the deformation of a 3D ellipsoidal shell under applied stretching and twisting forces. The provided functions can be adapted to simulate other 3D geometries.

## Requirements

To use this solver, ensure you have the following installed:

- **MATLAB**: Version R2020b or later (required for the `pagemtimes` function).
- **GMSH**: For generating the mesh from `.geo` files.
- **ParaView** (optional): For advanced visualization of the generated VTK files.
- A computer with sufficient memory to handle a tetrahedral mesh and up to 15,000 iterations.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. **Install GMSH**:
   - Download and install GMSH from [http://gmsh.info/](http://gmsh.info/).
   - Ensure GMSH is accessible from the command line.

3. **Set Up MATLAB**:
   - Ensure MATLAB is installed.
   - Add the repository root directory to MATLAB's path:
     ```matlab
     addpath(genpath('<repository-directory>'));
     ```

## Running the Simulation

1. **Generate the Mesh**:
   - Navigate to the `tutorials/ellipsoidal_shell` directory.
   - Use GMSH to create the mesh file from the provided geometry file:
     ```bash
     gmsh ellipsoidal_shell_mesh.geo -3 -format m -o ellipsoidal_shell_mesh.m
     ```
   - This generates `ellipsoidal_shell_mesh.m`, the mesh file used by the simulation.

2. **Run the Simulation**:
   - Open MATLAB and navigate to the `tutorials/ellipsoidal_shell` directory.
   - Run the simulation script:
     ```matlab
     ellipsoidal_shell
     ```
   - The simulation runs for 15,000 iterations, applying stretching and twisting forces for the first 7,500 iterations, then relaxing the structure. Output images and VTK files are saved in the `images` and `vtk` directories every 50 iterations.

## Visualization

- The simulation generates images in `tutorials/ellipsoidal_shell/images/` showing:
  - The deformed tetrahedral mesh (blue elements).
  - Force vectors applied at selected nodes.
- VTK files are generated in `tutorials/ellipsoidal_shell/vtk/` for advanced visualization in ParaView.

## Code Structure

The repository includes the following key files:

- **`calculate_forces.m`**:
  - Computes nodal forces for each tetrahedral element based on deformation tensors, Jacobian derivatives, and material properties.
  - Implements force calculation methods for 3D FEM simulations.

- **`cofactor.m`**:
  - Calculates the cofactor matrix (signed minors) of a square matrix, used in shape function coefficient calculations.

- **`create_festruct.m`**:
  - Initializes a structure (`festruct`) to store mesh and simulation data, streamlining data passing in FEM computations.

- **`shapefunctioncoefficients.m`**:
  - Computes shape function coefficients and element volumes for tetrahedral finite elements.

- **`tutorials/ellipsoidal_shell/ellipsoidal_shell.m`**:
  - The main simulation script that orchestrates the simulation, applies boundary conditions, and iterates through time steps.

- **`tutorials/ellipsoidal_shell/ellipsoidal_shell_mesh.geo`**:
  - Geometry file for the 3D ellipsoidal shell mesh.

## Simulation Details

- **Geometry**: A 3D ellipsoidal shell.
- **Boundary Conditions**: 
  - Stretching forces of magnitude 500 are applied along the z-axis (upward on the top patch, downward on the bottom patch) for the first 7,500 iterations.
  - Twisting (azimuthal) forces of magnitude 500 are applied on the top and bottom patches in opposite directions for the first 7,500 iterations.
  - Forces are applied to nodes within a 45° cone relative to the z-axis.
- **Time Stepping**: 15,000 iterations with a time step `dt = 0.001`.
- **Output**: Mesh deformation and force vectors are visualized every 50 iterations in the `images` and `vtk` directories.

## Customizing the Simulation

- **Parameter Modification**: Adjust parameters in `ellipsoidal_shell.m`:
  - `fmag`: Force magnitude (default: 500).
  - `niter`: Number of iterations (default: 15,000).
  - `dt`: Time step (default: 0.001).
  - `mu`: Damping coefficient (default: 1000).
  - `K` and `co`: Material constants (default: 100,000 and 50,000 respectively).

- **Working with Custom Meshes**:
  - Create a custom `.geo` file for your 3D geometry.
  - Define physical groups for boundaries and domains of interest.
  - Generate the MATLAB mesh file using GMSH.
  - Adapt the boundary identification code (e.g., `get_patches`) to your physical groups.
  - Update the force application logic to match your new geometry.

## Notes

- Ensure MATLAB's path includes the repository root directory.
- The simulation duration depends on the number of iterations and mesh size; progress is visualized every 50 iterations.
- Modify `ellipsoidal_shell.m` to adjust simulation parameters or boundary conditions.
- VTK files can be visualized in ParaView for detailed analysis of the deformation and forces.
