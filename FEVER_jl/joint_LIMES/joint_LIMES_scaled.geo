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

// Transfinite Curve {35,43} = 5 Using Progression 1;
// Transfinite Curve {14} = 250 Using Progression 1;

Physical Surface("Cond_in", 54) = {3};
Physical Surface("Cond_out", 53) = {1};
Physical Surface("Insul_in_L", 55) = {2};
Physical Surface("Insul_in_R", 56) = {7};
Physical Surface("Semicon_in", 57) = {4};
Physical Surface("Semicon_out", 58) = {6};
Physical Surface("Insul_out", 59) = {5};
// 
Physical Curve("interf_cond_in", 60) = {5, 2, 7, 19, 18, 23, 22, 4, 16};
Physical Curve("interf_semicon_out", 61) = {9, 24, 33, 32, 31, 30, 42};
Physical Curve("axis", 62) = {14};
Physical Curve("interf_joint_out", 63) = {10, 35, 34, 45, 44, 43, 47};
//

Mesh.MeshSizeMax = 3/1000;
Mesh 2;
Mesh.MshFileVersion = 2;
Save "joint_LIMES_scaled.msh";//+

