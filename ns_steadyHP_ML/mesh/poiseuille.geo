Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 0.2, 0, 1};
Point(4) = {0, 0.2, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Physical Line(0) = {1,3}; //top1+bottom1
Physical Line(1) = {2}; //inlet
Physical Line(2) = {4}; //outlet

Plane Surface(6) = {5};
Physical Surface(7) = {6};
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.02;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 40;
Show "*";

