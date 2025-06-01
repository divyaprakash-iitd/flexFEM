import gmsh
import numpy as np
import matplotlib.pyplot as plt

def plot_nodes(p, ax,m='*',ms=2):
    ax.scatter(p[:,0],p[:,1],marker=m,s=ms)

# Initialize Gmsh
gmsh.initialize()

# Read the mesh file
gmsh.open("donut2d_mesh.msh")

# Now you can access all mesh data
print("Mesh successfully loaded!")

# Don't forget to finalize
#gmsh.finalize()

# Figure out the total number of physical groups
# Physical groups that we are ineterested in is 1D
physicalGroups = gmsh.model.getPhysicalGroups()
nPhysicalGroups = len(physicalGroups)

fig, ax = plt.subplots()

for i, s in zip(range(nPhysicalGroups),[80,40,10]):
    p = gmsh.model.mesh.getNodesForPhysicalGroup(*physicalGroups[i])
    x = np.array(p[1]).reshape(-1,3)
    plot_nodes(x,ax,ms=s)

# # Get nodes according to physical groups
# p1 = gmsh.model.mesh.getNodesForPhysicalGroup(*physicalGroups[0])
# p2 = gmsh.model.mesh.getNodesForPhysicalGroup(*physicalGroups[1])
# p3 = gmsh.model.mesh.getNodesForPhysicalGroup(*physicalGroups[2])

# x1 = np.array(p1[1]).reshape(-1,3)
# x2 = np.array(p2[1]).reshape(-1,3)
# x3 = np.array(p3[1]).reshape(-1,3)

# # Plot the nodes
# plt.plot(x1[:,0],x1[:,1],'ko')
# plt.plot(x2[:,0],x2[:,1],'bo')
# plt.scatter(x3[:,0],x3[:,1],marker="*",s=2)

plt.show()
