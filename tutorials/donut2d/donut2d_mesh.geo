// annular_mesh.geo
// GMSH script to create an annular domain (circle with hole) meshed with triangles
// Both inner and outer circumferences are defined as named physical groups (patches)
// Exports mesh in MATLAB format for easy access to circumference points

// Parameters
outerRadius = 1.0;      // Outer circle radius
innerRadius = 0.4;      // Inner circle radius (hole)
meshSize = 0.1;         // Target mesh element size
centerX = 0.0;          // X-coordinate of circle center
centerY = 0.0;          // Y-coordinate of circle center

// Create the center point
Point(1) = {centerX, centerY, 0, meshSize};

// Create 4 points on the outer circle perimeter (90 degree intervals)
Point(2) = {centerX + outerRadius, centerY, 0, meshSize};         // Right
Point(3) = {centerX, centerY + outerRadius, 0, meshSize};         // Top
Point(4) = {centerX - outerRadius, centerY, 0, meshSize};         // Left
Point(5) = {centerX, centerY - outerRadius, 0, meshSize};         // Bottom

// Create 4 points on the inner circle perimeter (90 degree intervals)
Point(6) = {centerX + innerRadius, centerY, 0, meshSize};         // Right
Point(7) = {centerX, centerY + innerRadius, 0, meshSize};         // Top
Point(8) = {centerX - innerRadius, centerY, 0, meshSize};         // Left
Point(9) = {centerX, centerY - innerRadius, 0, meshSize};         // Bottom

// Create 4 arcs to form the outer circle
Circle(1) = {2, 1, 3};  // First quadrant arc (outer)
Circle(2) = {3, 1, 4};  // Second quadrant arc (outer)
Circle(3) = {4, 1, 5};  // Third quadrant arc (outer)
Circle(4) = {5, 1, 2};  // Fourth quadrant arc (outer)

// Create 4 arcs to form the inner circle (hole)
Circle(5) = {6, 1, 7};  // First quadrant arc (inner)
Circle(6) = {7, 1, 8};  // Second quadrant arc (inner)
Circle(7) = {8, 1, 9};  // Third quadrant arc (inner)
Circle(8) = {9, 1, 6};  // Fourth quadrant arc (inner)

// Create a curve loop for the outer circle
Curve Loop(1) = {1, 2, 3, 4};

// Create a curve loop for the inner circle (hole)
// Note: The inner curve must be in the opposite direction
// to create a hole, so we reverse the order
Curve Loop(2) = {5, 6, 7, 8};

// Create a plane surface with a hole
// First argument is the outer loop, second is the inner loop
Plane Surface(1) = {1, 2};

// Define the mesh type as triangles
Mesh.Algorithm = 5;  // Delaunay meshing algorithm

// Define the outer boundary as a physical group named "OuterCircumference"
Physical Curve("OuterCircumference") = {1, 2, 3, 4};

// Define the inner boundary as a physical group named "InnerCircumference"
Physical Curve("InnerCircumference") = {5, 6, 7, 8};

// Define the annular surface as a physical group named "AnnularSurface"
Physical Surface("AnnularSurface") = {1};

// Set some visualization options
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.MeshSizeMin = meshSize / 2;
Mesh.MeshSizeMax = meshSize;

// Generate the 2D mesh
Mesh 2;

