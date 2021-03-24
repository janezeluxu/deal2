Point(1) = {0, 0, 0, 1};
Point(2) = {0.2, 0, 0, 1};
Point(3) = {1, 0, 0, 1};
Point(4) = {1, 0.2, 0, 1};
Point(5) = {0.2, 0.2, 0, 1};
Point(6) = {0, 0.2, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line Loop(5) = {1, 2, 3, 4,5,6};
Physical Line(0) = {1,5}; //top1+bottom1
Physical Line(1) = {2,4}; //top2+bottom2
Physical Line(2) = {6}; //inlet
Physical Line(3) = {3}; //outlet

Plane Surface(6) = {5};
Physical Surface(7) = {6};
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.02;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 40;
Show "*";

