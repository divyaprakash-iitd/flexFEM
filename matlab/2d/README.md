# Implementation in MATLAB (2D)

## Overview

This repository provides a MATLAB-based library for Finite Element Method (FEM) simulations. For demonstartion purposes a
deformation of annular (donut-shaped) structures under applied forces is simulated. However the provided functions can be used
to simulate any type of 2D geometry.

## Requirements

To use this solver, ensure you have the following installed:

- **MATLAB**: Version R2020b or later (required for the `pagemtimes` function).
- **GMSH**: For generating the mesh from `.geo` files.
- **ParaView** (optional): For advanced visualization of the generated VTK files.
- A computer with sufficient memory to handle a mesh with 428 nodes and up to 50,000 iterations.

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
   - Navigate to the `tutorials/donut2d` directory.
   - Use GMSH to create the mesh file from the provided geometry file:
     ```bash
     gmsh donut2d_mesh.geo -2 -format m -o donut2d_mesh.m
     ```
   - This generates `donut2d_mesh.m`, the mesh file used by the simulation.

2. **Run the Simulation**:
   - Open MATLAB and navigate to the `tutorials/donut2d` directory.
   - Run the simulation script:
     ```matlab
     donut2d
     ```
   - The simulation runs for 50,000 iterations, applying forces for the first half and then relaxing the structure. Output images are saved in the `images` directory every 500 iterations.

## Visualization

- The simulation generates images in `tutorials/donut2d/images/` showing:
  - The deformed mesh (red triangles)
  - Applied forces (black arrows) on boundary nodes
  - Total forces (blue arrows) at selected nodes
- For advanced visualization, you can use ParaView to view the generated VTK files (if generated).

## Code Structure

The repository includes the following key files:

- **`calculate_forces.m`**:
  - Computes nodal forces for each triangular element based on deformation tensors, Jacobian derivatives, and material properties.
  - Implements the force calculation methods described in the referenced papers.

- **`cofactor.m`**:
  - Calculates the cofactor matrix (signed minors) of a square matrix, used in shape function coefficient calculations.

- **`create_festruct.m`**:
  - Initializes a structure (`festruct`) to store mesh and simulation data, streamlining data passing in FEM computations.

- **`shapefunctioncoefficients.m`**:
  - Computes shape function coefficients and element areas for triangular finite elements.

- **`tutorials/donut2d/donut2d.m`**:
  - The main simulation script that orchestrates the simulation, applies boundary conditions, and iterates through time steps.

- **`tutorials/donut2d/donut2d_mesh.geo`**:
  - Geometry file for the annular ring mesh.

## Simulation Details

- **Geometry**: A 2D annular ring with an outer radius of 1.0 unit.
- **Boundary Conditions**: Forces of magnitude 500 are applied horizontally (outward on the right, inward on the left) on the periphery between -45° and 45° relative to the y-axis for the first 25,000 iterations.
- **Time Stepping**: 50,000 iterations with a time step `dt = 0.001`.
- **Output**: Mesh deformation and force vectors are visualized every 500 iterations in the `images` directory.

## Customizing the Simulation

- **Parameter Modification**: Adjust parameters in `donut2d.m`:
  - `Ftip`: Force magnitude (default: 500)
  - `niter`: Number of iterations (default: 50,000)
  - `dt`: Time step (default: 0.001)
  - `mu`: Damping coefficient (default: 1000)
  - `K` and `co`: Material constants (default: 10000 and 5000 respectively)

- **Working with Custom Meshes**:
  - Create a custom `.geo` file for your geometry.
  - Define physical groups for boundaries and domains of interest.
  - Generate the MATLAB mesh file using GMSH.
  - Adapt the boundary identification code to your physical groups.
  - Update the force application to match your new geometry.

## Notes

- Ensure MATLAB's path includes the repository root directory.
- The simulation duration depends on the number of iterations and mesh size; progress is printed every 1,000 iterations.
- Modify `donut2d.m` to adjust simulation parameters or boundary conditions.
