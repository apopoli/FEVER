// ----------------------------
// Parametri [m]
// ----------------------------
rc    = 20e-3;
r_sc1 = rc + 1e-3;
r_xlpe = r_sc1 + 30e-3;
r_sc2 = r_xlpe + 1e-3;
r_al  = r_sc2 + 2e-3;
r_ext = r_al + 4e-3;

Lz = 0.1;

// Mesh size (personalizza)
lc_core = 1.0e-3;
lc_thin = 0.5e-3;
lc_bulk = 2.0e-3;

// ----------------------------
// Punti: colonne a r = [0, rc, r_sc1, r_xlpe, r_sc2, r_al, r_ext]
// e righe a z = [0, Lz]
// ----------------------------
rlist[] = {0, rc, r_sc1, r_xlpe, r_sc2, r_al, r_ext};

p00[] = {}; // punti a z=0
p0L[] = {}; // punti a z=Lz

For i In {0:#rlist[]-1}
  r = rlist[i];
  p00[i] = newp; Point(p00[i]) = {r, 0, 0, 0.001};
  p0L[i] = newp; Point(p0L[i]) = {r, Lz, 0, 0.001};
EndFor

// ----------------------------
// Linee verticali (a r costante) e orizzontali (a z costante)
// ----------------------------
vline[] = {};
tline[] = {}; // top z=Lz
bline[] = {}; // bottom z=0

For i In {0:#rlist[]-1}
  vline[i] = newl; Line(vline[i]) = {p00[i], p0L[i]};
EndFor

For i In {0:#rlist[]-2}
  bline[i] = newl; Line(bline[i]) = {p00[i], p00[i+1]};
  tline[i] = newl; Line(tline[i]) = {p0L[i], p0L[i+1]};
EndFor

// ----------------------------
// Superfici: una per ogni strato (rettangoli tra due raggi)
// Ogni rettangolo: (r_i,z=0) -> (r_{i+1},z=0) -> (r_{i+1},z=Lz) -> (r_i,z=Lz)
// ----------------------------
surfaces[] = {};

For i In {0:#rlist[]-2}
  ll = newll;
  // attenzione all'orientazione: chiudiamo il loop in senso coerente
  Line Loop(ll) = { bline[i], vline[i+1], -tline[i], -vline[i] };
  s = news;
  Plane Surface(s) = {ll};
  surfaces[i] = s;
EndFor

// ----------------------------
// Physical Surfaces (materiali) - ordine:
// 0: cu, 1: sc_in, 2: xlpe, 3: sc_out, 4: al, 5: cover
// ----------------------------
Physical Surface("cu")     = {surfaces[0]};
Physical Surface("sc_in")  = {surfaces[1]};
Physical Surface("xlpe")   = {surfaces[2]};
Physical Surface("sc_out") = {surfaces[3]};
Physical Surface("al")     = {surfaces[4]};
Physical Surface("cover")  = {surfaces[5]};

// ----------------------------
// Physical Curves (BC)
// axis   : r=0  -> vline[0]
// outer  : r=r_ext -> vline[last]
// bottom : z=0 -> tutte bline[]
// top    : z=Lz -> tutte tline[]
// ----------------------------
Physical Curve("axis")  = {vline[0]};
Physical Curve("outer") = {vline[#vline[]-1]};

Physical Curve("bottom") = {bline[]};
Physical Curve("top")    = {tline[]};

// (opzionale) se vuoi anche i bordi interni per debug:
// Physical Curve("interfaces") = {vline[1], vline[2], vline[3], vline[4], vline[5]};

// ----------------------------
// Mesh options
// ----------------------------
Mesh.Algorithm = 6;
