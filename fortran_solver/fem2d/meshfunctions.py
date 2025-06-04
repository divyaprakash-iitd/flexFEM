import numpy as np
from scipy.spatial import Delaunay
from matplotlib.path import Path

# Generate points for a triangle


def rectanglecor(L, B, Nx, Ny):
    """Genrates the delaunay triangulation for a rectangle

    Args:
        L (float): Rectangle's Length
        B (float): Rectangle's Breadth
        Nx (int): Number of nodes along x-direction
        Ny (int): Number of nodes along y-direction

    Returns:
        Delaunay object: The Delaunay triangulation object
    """
    X = np.linspace(0, B, Nx)
    Y = np.linspace(0, L, Ny)
    X, Y = np.meshgrid(X, Y)

    X = X.reshape((X.size, -1))
    Y = Y.reshape((Y.size, -1))

    M = Delaunay(np.concatenate((X, Y), axis=1))

    return M


def cofactor(A):
    """Generates the cofactors (signed minors) of a given matrix

    Args:
        A (float): A 2D numpy array

    Returns:
        float: A 2D numpy array containig the cofactors
    """
    return np.linalg.inv(A).T * np.linalg.det(A)

# Calculate shape function coefficients


def shapefunctioncoefficients(M):
    """Generates the shape function coefficients for the given mesh

    Args:
        M (Delaunay triangulation object): Delaunay triangulation

    Returns:
        tuple: A tuple containing the coefficients and elemental area
    """
    nnode = M.simplices.shape[1]  # No. of nodes in the FE (triangle)
    ncor = M.ndim   # No. of coordinates
    nelem = M.nsimplex  # No. of elements in the FE mesh

    # Shape function coefficients
    a = np.zeros((nelem, 1, nnode))
    b = np.zeros((nelem, ncor, nnode))

    # Go over each of the elements
    # For each element a matrix is created, the cofactors of
    # which give us the coefficient of the shape function

    A = np.ones((3, 3))
    Aelem = np.zeros((nelem))
    for ielem in range(nelem):
        A[0:-1, :] = M.points[M.simplices[ielem]].T
        Aelem[ielem] = np.abs(np.linalg.det(A))  # Elemental area (+ve)
        # Go over each node of the element
        # Calculate the coefficients for each node (a and b)
        C = cofactor(A)
        b[ielem, :, :] = C[0:-1, :]
        a[ielem, :, :] = C[-1, :]

        # Divide by area of the element
        b[ielem, :, :] = b[ielem, :, :]/Aelem[ielem]
        a[ielem, :, :] = a[ielem, :, :]/Aelem[ielem]

    return a, b, Aelem

def ellipsecor(a,b,theta):
    # Generate a particle's points and boundary points

    px = a*np.cos(theta)
    py = b*np.sin(theta)
    ds = np.linalg.norm([px[0]-px[1], py[0]-py[1]])
    boundary = np.vstack((px, py)).T

    # Create a background canvas
    x = np.arange(-a+ds, a, ds)
    y = np.arange(-b+ds, b, ds)
    x, y = np.meshgrid(x, y)
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x, y)).T

    p = Path(boundary)
    grid = p.contains_points(points)

    inpoints = np.vstack((x[grid], y[grid])).T

    M = Delaunay(np.vstack((boundary, inpoints)))

    return M, boundary

def particlecor(radius, N):
    # Generate a particle's points and boundary points
    theta = np.linspace(0, 2*np.pi, N)

    px = radius*np.cos(theta)
    py = radius*np.sin(theta)
    ds = np.linalg.norm([px[0]-px[1], py[0]-py[1]])
    boundary = np.vstack((px, py)).T

    # Create a background canvas
    x = np.arange(-radius+ds, radius, ds)
    y = np.arange(-radius+ds, radius, ds)
    x, y = np.meshgrid(x, y)
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x, y)).T

    p = Path(boundary)
    grid = p.contains_points(points)

    inpoints = np.vstack((x[grid], y[grid])).T

    M = Delaunay(np.vstack((boundary, inpoints)))

    return M, boundary
