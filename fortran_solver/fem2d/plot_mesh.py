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
    
    if point_data is not None:
        if point_data.ndim == 1:
            grid.point_data['data'] = point_data
        else:
            for i in range(point_data.shape[1]):
                grid.point_data[f'data_{i}'] = point_data[:, i]
    
    return grid

conn = np.loadtxt("connectivity.txt", dtype=int)  # shape (M, 3)
# Subtract 1 to get the python index starting from 0
conn = conn - 1
files = glob.glob("P_*.txt")
sorted_files = sorted(files, key=lambda x: int(x.split('_')[1].split('.')[0]))
for i, file in enumerate(sorted_files):
    coords = np.loadtxt(file)  # shape (N, 2) or (N, 3)
    mesh = create_vtk_mesh(coords.T, conn.T)
    mesh.save(f"mesh_{i:04d}.vtu")

#plotter = pv.Plotter()
#for i in range(len(files)):
#    mesh = pv.read(f"mesh_{i+1:04d}.vtu")
#    plotter.clear()
#    plotter.add_mesh(mesh)
#    plotter.show(auto_close=False)

