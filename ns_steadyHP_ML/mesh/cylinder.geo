// ---- 2D Circular Cylinder Gmsh Tutorial ----
// 2D_cylinder_tutorial.geo
// --------------------------------------------

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 3;
cl2 = 0.5;
cl3 = 5;

radius = 0.5;

// Exterior (bounding box) of mesh
Point(1) = {-30, -30, 0, cl1};
Point(2) = { 50, -30, 0, cl3};
Point(3) = { 50,  30, 0, cl3};
Point(4) = {-30,  30, 0, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Circle
Point(5) = {0,   0, 0, cl2};
Point(6) = {0,  radius, 0, cl2};
Point(7) = {0, -radius, 0, cl2};
Point(8) = {radius,  0, 0, cl2};
Point(9) = {-radius, 0, 0, cl2};

Circle(5) = {7, 5, 8};
Circle(6) = {6, 5, 9};
Circle(13) = {8, 5, 6};
Circle(14) = {9, 5, 7};

Line Loop(5) = {1, 2, 3, 4, 5, 6, 13, 14}; // Exterior


Physical Line(0) = {1,3};
Physical Line(1)  = {2};
Physical Line(2)    = {4};
Physical Line(3)   = {5,6,13,14};


Plane Surface(6) = {5};
Physical Surface(7) = {6};

Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
