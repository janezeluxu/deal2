Point(1) = {0, 0.5, 0, 1};
Point(2) = {2, 0.5, 0, 1};
Point(3) = {2, 0, 0, 1};
Point(4) = {8, 0, 0, 1};
Point(5) = {8, 1, 0, 1};
Point(6) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line Loop(5) = {1, 2, 3, 4,5,6};
Physical Line(0) = {1,2,3,5}; //top+bottom
Physical Line(1) = {6}; //inlet
Physical Line(2) = {4}; //outlet

Plane Surface(6) = {5};
Physical Surface(7) = {6};

Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.05;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";

