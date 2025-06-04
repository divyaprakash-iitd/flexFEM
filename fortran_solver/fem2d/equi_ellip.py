import numpy as np
import scipy as sp
import scipy.optimize
import f90nml
import meshfunctions as mf

def angles_in_ellipse(
        num,
        a,
        b):
    assert(num > 0)
    #assert(a < b)
    angles = 2 * np.pi * np.arange(num) / num
    if a != b:
        e2 = (1.0 - a ** 2.0 / b ** 2.0)
        tot_size = sp.special.ellipeinc(2.0 * np.pi, e2)
        arc_size = tot_size / num
        arcs = np.arange(num) * arc_size
        res = sp.optimize.root(
            lambda x: (sp.special.ellipeinc(x, e2) - arcs), angles)
        angles = res.x 
    return angles

nml = f90nml.read('input_params.dat')
a = nml['particleprops']['aa']
b = nml['particleprops']['bb']
n = nml['particleprops']['nvp']

phi = angles_in_ellipse(n, a, b)

# Define the file path
file_path = 'equidistant_ellipse_angles.txt'

# Write the array elements as a column to the file
with open(file_path, 'w') as file:
    for element in phi:
        file.write(f"{element}\n")

# Generate a particle
MP, bpoints = mf.ellipsecor(a,b,phi)
x = MP.points[:, 0]
y = MP.points[:, 1]
PP = np.zeros((MP.npoints, 2), order='F')
PP[:, 0] = x
PP[:, 1] = y
# boundary[:, :] = False
_, pb, paelem = mf.shapefunctioncoefficients(MP)

femdata = np.ones((2))
femdata[0] = len(x)
femdata[1] = len(paelem)
# femdata = (np.rint(femdata)).astype(int32)
femdata = femdata.astype(np.int32)

np.savetxt("bpoints.csv", bpoints)
np.savetxt("mpnumpy.csv", MP.simplices)
np.asfortranarray(femdata).tofile("femdata.bin")
np.asfortranarray(pb).T.tofile("pb.bin")
np.asfortranarray(paelem).T.tofile("paelem.bin")
np.asfortranarray((MP.simplices+1).astype(np.int32)).T.tofile("mp.bin")
PP.T.tofile("pp.bin")
## plotting
#import matplotlib.pyplot as plt
#fig = plt.figure()
#ax = fig.gca()
#ax.axes.set_aspect('equal')
#ax.scatter(b * np.sin(phi), a * np.cos(phi))
#
##phi = np.linspace(0,2*np.pi,n)
##ax.scatter(b * np.sin(phi), a * np.cos(phi))
#
#plt.show()
#
