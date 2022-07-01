//+
Point(1) = {-0.5, 0.03125, 0, 1.0};
//+
Point(2) = {0.5, 0.03125, 0, 1.0};
//+
Point(3) = {0.5, -0.03125, 0, 1.0};
//+
Point(4) = {-0.5, -0.03125, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Physical Curve("top", 5) = {1};
//+
Physical Curve("bottom", 6) = {3};
//+
Physical Curve("left", 7) = {4};
//+
Physical Curve("right", 8) = {2};
//+
Plane Surface(1) = {1};
//+
Physical Surface("domain", 9) = {1};
//+
Physical Curve("edges", 10) = {1, 4, 3, 2};
