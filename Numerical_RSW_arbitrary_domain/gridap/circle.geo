//+
SetFactory("OpenCASCADE");
Circle(1) = {-0, -0, 0, 0.5, 0, 2*Pi};
//+
Physical Curve("edges", 2) = {1};
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
//+
Physical Surface("domain", 3) = {1};
