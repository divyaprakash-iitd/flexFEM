# 2D FEM Simulation Library for Annular Deformation

![Animation](donut_animation.webp)

This repository provides a MATLAB-based library for 2D Finite Element Method (FEM) simulations, focusing on the deformation of annular (donut-shaped) structures under applied forces. The implementation follows the formulations described in:

- Bhattacharya, A., & Balazs, A. C. (2013). Stiffness-modulated motion of soft microscopic particles over active adhesive cilia. Soft Matter, 9(15), 3945-3955. https://doi.org/10.1039/c3sm27445f
- Duki, S. F., Kolmakov, G. V., Yashin, V. V., Kowalewski, T., Matyjaszewski, K., & Balazs, A. C. (2011). Modeling the nanoscratching of self-healing materials. The Journal of Chemical Physics, 134(8), 084901. https://doi.org/10.1063/1.3554193


- A. Bhattacharya and A. C. Balazs, [Soft Matter, 2013, 9, 3945–3955](https://pubs.rsc.org/en/content/articlelanding/2013/sm/c3sm00028a)
- More detailed formulations in [Journal of Chemical Physics, 134, 084901 (2011)](https://pubs.aip.org/aip/jcp/article/134/8/084901/960252)

## Prerequisites

- **MATLAB**: Version R2020b or later (required for `pagemtimes` function).
- **Gmsh**: For generating the mesh from `.geo` files.
- A computer with sufficient memory to handle a mesh with 428 nodes and up to 50,000 iterations.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. **Install Gmsh**:
   - Download and install Gmsh from [http://gmsh.info/](http://gmsh.info/).
   - Ensure Gmsh is accessible from the command line.

3. **Set Up MATLAB**:
   - Ensure MATLAB is installed.
   - Add the repository root directory to MATLAB's path:
     ```matlab
     addpath(genpath('<repository-directory>'));
     ```

## Core Library Source Code

The library consists of modular source code files that implement the FEM formulations from the referenced papers:

### `calculate_forces.m`

**Purpose**: Computes nodal forces for each triangular element based on deformation tensors, Jacobian derivatives, and material properties.

**Input**: A structure (`festruct`) containing:
- `M`: Mesh connectivity
- `x`: Current node coordinates
- `xorg`: Original node coordinates
- `FN`: Force vectors
- `co`, `K`: Material constants
- `b`: Shape function coefficients
- `Aelem`: Element areas

**Output**: Updates `festruct.FN` with computed nodal forces.

**Key Operations**:
- Calculates deformation tensor using `pagemtimes`
- Computes invariant and Jacobian derivative terms
- Applies material properties to compute forces per element
- Sums forces at each node

This implementation closely follows the force calculation methods described in the papers, particularly the formulations for hyperelastic materials.

### `cofactor.m`

**Purpose**: Calculates the cofactor matrix (signed minors) of a square matrix, used in shape function coefficient calculations.

**Input**: A square matrix `A`.

**Output**: Matrix `C` containing cofactors of each element.

**Key Operations**:
- Iterates over each element, removes corresponding row and column, and computes the determinant with a sign factor `(-1)^(i+j)`
- Used by `shapefunctioncoefficients.m` for triangular elements

### `create_festruct.m`

**Purpose**: Initializes a structure (`festruct`) to store mesh and simulation data, streamlining data passing in FEM computations.

**Input**: 
- `M`: Mesh connectivity
- `x`, `y`: Node coordinates
- `FN`: Initial forces
- `co`, `K`: Material constants

**Output**: A structure with fields `M`, `x`, `xorg`, `FN`, `co`, `K`, `nNodes`, `nElem`, `b`, and `Aelem`.

**Key Operations**:
- Combines coordinates into `festruct.x` and `festruct.xorg`
- Computes shape function coefficients using `shapefunctioncoefficients.m`

### `shapefunctioncoefficients.m`

**Purpose**: Computes shape function coefficients and element areas for triangular finite elements.

**Input**: 
- `M`: Mesh connectivity
- `x`, `y`: Node coordinates

**Output**: Coefficients `a`, `b`, and element areas `Aelem`.

**Key Operations**:
- Constructs a 3x3 matrix for each element using node coordinates
- Uses `cofactor.m` to compute cofactors, which yield shape function coefficients
- Calculates element areas via determinant

## Understanding GMSH Meshes

When creating a mesh in GMSH (either through the GUI or scripting), you can define physical groups to identify specific regions in your geometry. These physical groups are essential for applying boundary conditions in FEM simulations.

### Mesh Structure and Physical Groups

In the case of our annular mesh, we define three physical groups:

1. **Outer circumference** (boundary): Tag 1
2. **Inner circumference** (boundary): Tag 2
3. **Annular surface** (domain): Tag 3

When the mesh is saved in MATLAB format (`.m` extension), it is loaded as a structure `msh` with these key fields:

| Field | Description |
|-------|-------------|
| `POS` | An `NP×3` matrix where `NP` is the number of points. Each row contains the X, Y, Z coordinates of a mesh node. |
| `LINES` | An `NL×3` matrix where `NL` is the number of lines. The first two columns contain indices of endpoints, and the third column contains the physical group tag. |
| `TRIANGLES` | An `NT×4` matrix where `NT` is the number of triangles. The first three columns contain vertex indices, and the fourth column contains the physical group tag. |
| `nbNod` | The total number of nodes in the mesh. |

### Extracting Boundary Nodes from the Mesh

To apply boundary conditions, we need to identify which nodes belong to each physical group:

```matlab
% Extract lines belonging to each boundary
ocLinesId = msh.LINES(:,end) == ocTag;  % Outer circumference
icLinesId = msh.LINES(:,end) == icTag;  % Inner circumference

% Get unique points on each boundary
ocPointsId = unique(msh.LINES(ocLinesId,1:2));
icPointsId = unique(msh.LINES(icLinesId,1:2));

% Create boolean masks for each boundary
nNodes = msh.nbNod;
ocMask = false(nNodes,1);
ocMask(ocPointsId) = true;  % True for nodes on outer circumference

icMask = false(nNodes,1);
icMask(icPointsId) = true;  % True for nodes on inner circumference
```

These masks can then be used to apply forces or constraints to specific boundaries.

## Tutorial: Annular Deformation Simulation

The `tutorials/donut2d` directory contains a tutorial demonstrating how to use the source code files to perform a 2D FEM simulation on an annular domain:

### Preparing the Mesh

1. Navigate to the `tutorials/donut2d` directory.
2. Generate the MATLAB mesh file using Gmsh:
   ```bash
   gmsh donut2d_mesh.geo -2 -format m -o donut2d_mesh.m
   ```

### Running the Simulation

The `donut2d.m` script simulates the deformation of an annular mesh under lateral forces:

```matlab
donut2d
```
### How the Tutorial Works

The tutorial `donut2d.m` demonstrates the complete simulation workflow:

1. **Setup Phase**: Loads the mesh, identifies boundaries using physical group tags, and initializes material parameters (bulk modulus `K=10000`, shear modulus `co=5000`).

2. **Boundary Selection**: Creates masks for specific portions of the outer boundary where forces will be applied, using angle and position criteria.

3. **Simulation Loop**: Runs 50,000 time steps with these key operations:
   - Applies lateral forces to opposite sides of the outer boundary for the first 25,000 iterations
   - Calculates internal element forces using the hyperelastic material model
   - Updates node positions using explicit time integration
   - Captures snapshots every 500 iterations

This simple workflow demonstrates how the library handles mesh processing, boundary condition application, and time evolution of the deformation, following the formulations described in the referenced papers.
### Output and Visualization

The simulation generates images in `tutorials/donut2d/images/`, showing:
- The deformed mesh (red triangles)
- Applied forces (black arrows) on boundary nodes
- Total forces (blue arrows) at selected nodes

## Customizing the Simulation

### Parameter Modification

You can customize the simulation by modifying parameters in `donut2d.m`:
- `Ftip`: Force magnitude (default: 500)
- `niter`: Number of iterations (default: 50,000)
- `dt`: Time step (default: 0.001)
- `mu`: Damping coefficient (default: 1000)
- `K` and `co`: Material constants (default: 10000 and 5000 respectively)

### Working with Custom Meshes

To use your own mesh:
1. Create a custom `.geo` file for your geometry
2. Define physical groups for boundaries and domains of interest
3. Generate the MATLAB mesh file using Gmsh
4. Adapt the boundary identification code to your physical groups
5. Update the force application to match your new geometry


## Author

- **Divyaprakash**
- **Email**: divyaprakash.poddar@gmail.com

## License

This project is licensed under the [MIT License](./LICENSE).  
Feel free to use, modify, and share it with proper attribution.

