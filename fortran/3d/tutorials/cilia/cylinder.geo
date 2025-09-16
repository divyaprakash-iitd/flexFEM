SetFactory("OpenCASCADE");

// Parameters
r = 0.1;     // radius
h = 1.0;     // height
lc = 0.05;    // mesh size

// Create the cylinder
Cylinder(1) = {0, 0, 0, 0, 0, h, r}; // base at (0,0,0), height along z, radius r

// Extract the surfaces of the cylinder
// These are automatically numbered: check Gmsh GUI or print for confirmation
SurfacesList[] = Boundary{ Volume{1}; };

// Usually for OpenCASCADE, the surfaces are returned in order:
// - one bottom disk
// - one top disk
// - one lateral surface

// Define physical groups (you might need to swap indices based on visualization)
Physical Surface("Bottom") = { SurfacesList[0] };
Physical Surface("Top") = { SurfacesList[1] };
Physical Surface("Side") = { SurfacesList[2] };

// Define the volume
Physical Volume("SolidCylinder") = {1};

// Set mesh size
Characteristic Length{ PointsOf{ Volume{1}; } } = lc;

