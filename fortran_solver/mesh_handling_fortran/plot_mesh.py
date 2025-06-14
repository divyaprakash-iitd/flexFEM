import numpy as np
import pyvista as pv
import os

def load_data(timestep):
    coords = np.loadtxt(f"coords_t{timestep}.txt")  # shape (N, 2) or (N, 3)
    if coords.shape[1] == 2:
        coords = np.hstack([coords, np.zeros((coords.shape[0], 1))])  # pad z=0
    
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


c = np.loadtxt('connectivity.txt',dtype=int)
# Subtract 1 to get the python index starting from 0
c = c - 1
x = np.loadtxt('coordAll.txt')

mesh = create_vtk_mesh(x.T,c.T)
mesh.save(f"mesh_1.vtu")
