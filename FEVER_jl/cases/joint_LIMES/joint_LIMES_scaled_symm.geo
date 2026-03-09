SetFactory("OpenCASCADE");
Merge "GIUNTO CAD con superpian2.iges";

Dilate {{0, 0, 0}, {1/1000, 1/1000, 1/1000}} {
  Curve{:}; Surface{:};
}

Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Curve{:}; Surface{:}; 
}

Translate {-0.41121107898583415, +1.322018717455649960-0.6, 0} {
  Curve{:}; Surface{:}; 
}

Coherence;

// Rettangolo molto grande che rappresenta la metà superiore
Rectangle(1000) = {0, -0.6, 0, 0.6, 0.6};
// Tieni solo l'intersezione tra la geometria e il rettangolo superiore
BooleanIntersection{ Surface{:}; Delete; }{ Surface{1000}; Delete; }

Recursive Delete {
  Surface{8}; 
}

Coherence;

Physical Surface("Cond_in", 53) = {3};
Physical Surface("Cond_out", 54) = {1};
Physical Surface("Insul_in", 55) = {7};
Physical Surface("Semicon_in", 56) = {4};
Physical Surface("Semicon_out", 57) = {5};
Physical Surface("Insul_out", 58) = {5};
Physical Surface("Semicon_out", 57) += {6};

Physical Curve("axis", 59) = {39};
Physical Curve("out_T", 60) = {31, 30, 29, 28};
Physical Curve("gnd", 61) = {27, 28, 21, 22, 43};
Physical Curve("HV", 62) = {14, 4, 17, 18, 42};

Mesh.MeshSizeMax = 1.5/1000;
Mesh 2;
Mesh.MshFileVersion = 2;
Save "joint_LIMES_scaled_symm.msh";

