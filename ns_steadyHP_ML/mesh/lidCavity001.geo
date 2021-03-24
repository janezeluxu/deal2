Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {0, 1, 0, 1};
Point(4) = {1, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line Loop(5) = {4, 1, 2, 3};
Physical Line(0) = {3};
Physical Line(1) = {4};
Physical Line(2) = {1};
Physical Line(3) = {2};

Plane Surface(6) = {5};
Physical Surface(7) = {6};
//Transfinite Surface {6} = {1, 2, 4, 3};
//Transfinite Line {4, 2, 3, 1} = 10 Using Progression 1;
//Recombine Surface {6};

Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.05;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
