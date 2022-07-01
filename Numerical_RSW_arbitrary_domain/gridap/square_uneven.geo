// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 10
//
//  Mesh size fields
//
// -----------------------------------------------------------------------------

// In addition to specifying target mesh sizes at the points of the geometry
// (see `t1.geo') or using a background mesh (see `t7.geo'), you can use general
// mesh size "Fields".

// Let's create a simple rectangular geometry
lc = 0.03;
Point(1) = {-0.5,-0.5,0,lc}; Point(2) = { 0.5,-0.5,0,lc};
Point(3) = { 0.5, 0.5,0,lc}; Point(4) = {-0.5, 0.5,0,lc};
Point(5) = { 0.0,-0.5,0,lc}; Point(6) = { 0.0, 0.5,0,lc};

Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1}; Line(5) = {5,6};

Curve Loop(5) = {1,2,3,4}; Plane Surface(6) = {5};

// Say we would like to obtain mesh elements with size lc/30 near curve 2 and
// point 5, and size lc elsewhere. To achieve this, we can use two fields:
// "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
// points 5 and on curve 2. This field returns the distance to point 5 and to
// (100 equidistant points on) curve 2.
Field[1] = Distance;
Field[1].CurvesList = {1,2,3,4};
Field[1].NumPointsPerCurve = 64;


// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc / 2;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.05;
Field[2].DistMax = 0.20;

// Finally, let's use the minimum of all the fields as the background mesh size
// field
Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

// To determine the size of mesh elements, Gmsh locally computes the minimum of
//
// 1) the size of the model bounding box;
// 2) if `Mesh.MeshSizeFromPoints' is set, the mesh size specified at
//    geometrical points;
// 3) if `Mesh.MeshSizeFromCurvature' is positive, the mesh size based on
//    curvature (the value specifying the number of elements per 2 * pi rad);
// 4) the background mesh size field;
// 5) any per-entity mesh size constraint.
//
// This value is then constrained in the interval [`Mesh.MeshSizeMin',
// `Mesh.MeshSizeMax'] and multiplied by `Mesh.MeshSizeFactor'.  In addition,
// boundary mesh sizes (on curves or surfaces) are interpolated inside the
// enclosed entity (surface or volume, respectively) if the option
// `Mesh.MeshSizeExtendFromBoundary' is set (which is the case by default).
//
// When the element size is fully specified by a background mesh size field (as
// it is in this example), it is thus often desirable to set

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// This will prevent over-refinement due to small mesh sizes on the boundary.

Curve Loop(1) = {1, 2, 3, 4};
Physical Curve("top", 5) = {1};
Physical Curve("bottom", 6) = {3};
Physical Curve("left", 7) = {4};
Physical Curve("right", 8) = {2};
Plane Surface(1) = {1};
Physical Surface("domain", 9) = {1};
Physical Curve("edges", 10) = {1, 2, 3, 4};
