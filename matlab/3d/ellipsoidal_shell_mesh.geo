SetFactory("OpenCASCADE");

// Create inner sphere and outer sphere
Sphere(1) = {0, 0, 0, 0.8};
Dilate {{0, 0, 0}, {1, 1, 1.3}} { Volume{1}; }
Sphere(2) = {0, 0, 0, 1.0};
Dilate {{0, 0, 0}, {1, 1, 1.3}} { Volume{2}; }

// Perform Boolean Difference to create hollow shell
BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };

// Get the surfaces of the hollow sphere
SurfacesList[] = Boundary{ Volume{3}; };

// Assuming only two surfaces are returned: inner and outer
// You may need to check which is which
// For symmetric case, the order is usually: {inner, outer}

// Assign physical groups
Physical Surface("InnerSurface") = { SurfacesList[0] };
Physical Surface("OuterSurface") = { SurfacesList[1] };

// You can also define the volume if needed
Physical Volume("HollowSphere") = { 3 };

// Set mesh size
Characteristic Length{ PointsOf{ Volume{3}; } } = 0.15;

