import numpy as np
import pyvista as pv
import os
import glob

images_dir = 'images/'

# Create images directory if it doesn't exist
if not os.path.exists(images_dir):
    os.makedirs(images_dir)

def create_vtk_mesh(coords, connectivity, point_data=None):
    n_cells = connectivity.shape[0]
    
    # Each tetrahedron has 4 points; PyVista expects [n_points_in_cell, id1, id2, id3, id4]
    cells = np.hstack([np.full((n_cells, 1), 4), connectivity]).flatten()
    celltypes = np.full(n_cells, pv.CellType.TETRA)

    grid = pv.UnstructuredGrid(cells, celltypes, coords)
    
    # Add point data if provided
    if point_data is not None:
        for name, data in point_data.items():
            grid.point_data[name] = data
    
    return grid

# Load connectivity and subtract 1 for 0-based indexing
conn = np.loadtxt("connectivity.txt", dtype=int)
conn = conn - 1  # assuming 1-based indices in your input

# Handle if shape is (4, M) instead of (M, 4)
if conn.shape[0] == 4:
    conn = conn.T

pfiles = sorted(glob.glob("P_*.txt"), key=lambda x: int(x.split('_')[1].split('.')[0]))
ffiles = sorted(glob.glob("F_*.txt"), key=lambda x: int(x.split('_')[1].split('.')[0]))

for i, (pfile, ffile) in enumerate(zip(pfiles, ffiles)):
    coords = np.loadtxt(pfile)  # shape (N, 3)
    forces = np.loadtxt(ffile)  # shape (N, 3)

    point_data = {
        'force': forces
    }

    mesh = create_vtk_mesh(coords, conn, point_data)
    mesh.save(f"{images_dir}mesh_{i:04d}.vtu")

