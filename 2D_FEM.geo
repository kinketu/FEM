// chalacteristic length definition
lc = 0.1;

// Definition of points
Point(1) = {0, 0, 0, lc};
Point(2) = {0.5, 0, 0, lc};
Point(3) = {0.5, 0.5, 0, lc};
Point(4) = {0, 0.5, 0, lc};

//Definition of lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//Difinition of surface
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};