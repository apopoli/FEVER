SetFactory("OpenCASCADE");
Merge "GIUNTO CAD con superpian2.iges";

Dilate {{0, 0, 0}, {1/1000, 1/1000, 1/1000}} {
  Curve{:}; Surface{:};
}

Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Curve{:}; Surface{:}; 
}

Translate {-0.41121107898583415, +1.322018717455649960, 0} {
  Curve{:}; Surface{:}; 
}

Coherence;

// Transfinite Curve {35,43} = 5 Using Progression 1;
// Transfinite Curve {14} = 250 Using Progression 1;

Physical Surface("Cond_in", 54) = {3};
Physical Surface("Cond_out", 53) = {1};
Physical Surface("Insul_in_L", 55) = {2};
Physical Surface("Insul_in_R", 56) = {7};
Physical Surface("Semicon_in", 57) = {4};
Physical Surface("Semicon_out", 58) = {6};
Physical Surface("Insul_out", 59) = {5};

Physical Curve("axis", 60) = {12};
Physical Curve("HV", 61) = {5, 2, 6, 16, 15, 18, 17, 4, 14};
Physical Curve("gnd", 62) = {8, 19, 24, 23, 22, 21};
Physical Curve("out", 63) = {9, 26, 25, 30, 29, 28, 32};

Mesh.MeshSizeMax = 3/1000;
Mesh 2;
Mesh.MshFileVersion = 2;
Save "joint_LIMES_scaled.msh";

