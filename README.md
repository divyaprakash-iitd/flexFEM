
# flexFEM: Finite Element Simulation Toolkit

This repository provides a toolkit for performing **Finite Element Method (FEM)** simulations on 2D/3D geometries, with two complementary implementations:

* A **MATLAB-based version** for prototyping and accessibility.
* A high-performance **Fortran-based version** suitable for large-scale simulations and integration with open-source toolchains.

The code is demonstrated on a soft 2D/3D structures deforming under applied forces, but is general enough to be adapted for arbitrary 2D/3D triangulated domains.

    Documentation containing implementation details are in the `doc` branch. However it is very basic and a work in progress.

![Animation](matlab/3d/tutorial/ellipsoidal_shell/ellipsoidal_shell.gif)
---

## Motivation

While MATLAB offers ease of use, the Fortran implementation is designed for **efficiency and openness**, utilizing freely available tools for meshing (GMSH), solving (Fortran), and visualization (PyVista/ParaView). The Fortran version is approximately **60× faster** than its MATLAB counterpart, making it suitable for extensive simulations or future coupling with external solvers.

---

## Applications and References

The implementations follow the formulations described in the following works:

* Bhattacharya, A., & Balazs, A. C. (2013). *Stiffness-modulated motion of soft microscopic particles over active adhesive cilia*. **Soft Matter**, 9(15), 3945–3955.
* Duki, S. F., et al. (2011). *Modeling the nanoscratching of self-healing materials*. **J. Chem. Phys.**, 134(8), 084901.

These formulations can be adapted to model a variety of soft structures interacting with applied forces or embedded in dynamic environments.

---

## Repository Structure

```
flexFEM/
├── fortran/    # High-performance implementation using Fortran and GMSH
├── matlab/     # MATLAB implementation with educational focus
├── LICENSE
└── README.md   # (This file)
```

Each subdirectory includes a separate README file with detailed instructions for compilation, simulation, and visualization.

---

## Getting Started

### 1. Fortran Version

* Requires `gfortran`, GMSH (compiled from source), and optionally `conda` for setting up Python-based visualization tools (PyVista).
* Compilation and execution are handled via a `Makefile` and shell scripts.
* Outputs VTK files for ParaView-based visualization.

➡ See `fortran/README.md` for setup and usage instructions.

---

### 2. MATLAB Version

* Requires MATLAB R2020b or newer (for `pagemtimes`) and a working GMSH installation.
* Suitable for quick experimentation and understanding of FEM concepts.
* Outputs PNG plots and optionally VTK files.

➡ See `matlab/README.md` for setup and usage instructions.

---

## Licensing

This repository is released under the [MIT License](./LICENSE). You are free to use, modify, and distribute the code with proper attribution.

---

## Contributions

Contributions, bug reports, and feature suggestions are welcome. Please open an issue or submit a pull request.

---
