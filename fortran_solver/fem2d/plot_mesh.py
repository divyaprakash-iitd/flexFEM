import numpy as np
import pyvista as pv
import os
import glob

def load_data(timestep):
    coords = np.loadtxt(f"coords_t{timestep}.txt")  # shape (N, 2) or (N, 3)
    
    connectivity = np.loadtxt(f"connectivity_t{timestep}.txt", dtype=int)  # shape (M, 3)
    data = np.loadtxt(f"data_t{timestep}.txt")  # shape (N,) or (N, k)
    
    return coords, connectivity, data

def create_vtk_mesh(coords, connectivity, point_data=None):
    n_cells = connectivity.shape[0]
    
    # Each triangle is 3 points; PyVista expects [n_points_in_cell, id1, id2, id3]
    cells = np.hstack([np.full((n_cells, 1), 3), connectivity]).flatten()
    celltypes = np.full(n_cells, pv.CellType.TRIANGLE)

    grid = pv.UnstructuredGrid(cells, celltypes, coords)
    
    # Add point data if provided
    if point_data is not None:
        for name, data in point_data.items():
            grid.point_data[name] = data  # e.g., 'force', 'displacement'
    
    return grid

conn = np.loadtxt("connectivity.txt", dtype=int)  # shape (M, 3)
# Subtract 1 to get the python index starting from 0
conn = conn - 1
pfiles = glob.glob("P_*.txt")
ffiles = glob.glob("F_*.txt")
sorted_files_p = sorted(pfiles, key=lambda x: int(x.split('_')[1].split('.')[0]))
sorted_files_f = sorted(ffiles, key=lambda x: int(x.split('_')[1].split('.')[0]))
for i, (pfile, ffile) in enumerate(zip(sorted_files_p, sorted_files_f)):
    coords = np.loadtxt(pfile)  # shape (N, 2) or (N, 3)
    forces = np.loadtxt(ffile)  # shape (N, 2) or (N, 3)
    point_data = {
        'force': forces.T          # Shape (N, 2) or (N, 3)
        #'displacement': disp.T      # Optional additional fields
    }
    mesh = create_vtk_mesh(coords.T, conn.T, point_data)
    mesh.save(f"mesh_{i:04d}.vtu")

#plotter = pv.Plotter()
#for i in range(len(files)):
#    mesh = pv.read(f"mesh_{i+1:04d}.vtu")
#    plotter.clear()
#    plotter.add_mesh(mesh)
#    plotter.show(auto_close=False)

