SetFactory("OpenCASCADE");
Merge "GIUNTO CAD con superpian2.iges";

Mesh.MeshSizeMax = 2;

Transfinite Curve {35,43} = 5 Using Progression 1;
// Transfinite Curve {14} = 250 Using Progression 1;

Physical Surface("Cond_in", 54) = {3};
Physical Surface("Cond_out", 53) = {1};
Physical Surface("Insul_in_L", 55) = {2};
Physical Surface("Insul_in_R", 56) = {7};
Physical Surface("Insul_out", 57) = {4};
Physical Surface("Semicon", 58) = {6};
Physical Surface("Insul_out", 59) = {5};

Physical Curve("interf_cond_in", 60) = {5, 2, 1, 4, 16};
Physical Curve("interf_semicon_out", 61) = {9, 24, 33, 32, 31, 30};

Mesh 2;
Mesh.MshFileVersion = 2;
Save "joint_LIMES.msh";
