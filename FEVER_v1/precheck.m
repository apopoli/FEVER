function [Area, gradL, refprop, Pdata, t] = precheck(p, t, ProbKind, ireg, materials, ifield, order)
 
% preprocessing per FEM
% Area = vettore di lunghezza nt = n° di elementi contenente le aree degli elementi 
% gradL = matrice [nt,2,3] contenente le matrici gradiente delle funzioni
% di forma per tutti gli elementi
% refprop = vettore di lunghezza np = né ni nodi, contenente il valore
% caratteristico locale della proprietà per ciascun nodo.
% elprop = vettore di lunghezza nt contenente il valore della proprietà su
% ciascun elemento.
% Pdata = struttura dati riguardante alcune info generali utili per la
% corsa:
%         Pdata.iprop = indice della proprietà in matlib
%         Pdata.invflag se = true -> la proprietà è il reciproco di quello che c'è in matlib
%         Pdata.RZflag se true -> problema assialsimmetrico
%         Pdata.cst = costante per la quale viene moltiplicato il termine
%         di sorgente
%         Pdata.order = ordine delle funzioni di forma
%         Pdata.Ksize = numero di righe nella matrice d'esemento
% t = matrice dei triangoli. preproc se necessario modifica l'ordine dei nodi 
% in modo che sia antiorario per tutti gli elementi

%   Detailed explanation goes here
nt = size(t,1);
np = size(p,1);
Area = zeros(nt,1);
elprop = zeros(nt,1);
refprop = zeros(np,1);
icnt = zeros(np,1);
gradL = zeros(nt,2,3);
KelSZ = [3 6];
switch ProbKind
    case 'magnetostatic'
        Pdata.iprop=2;  %indice della proprietà in matlib
        Pdata.invflag = true;   %true-> la proprietà è il reciproco di quello che c'è in matlib
        Pdata.RZflag = false; %true-> problema assialsimmetrico
        Pdata.cst = 1.2566370614359172e-6; % Vacuum permeability 
    case 'electrostatic'
        Pdata.iprop=1;
        Pdata.invflag = false;
        Pdata.RZflag = false;
        Pdata.cst = 1.129409067373019e+11;  % 1/Vacuum permittivity
    case 'electrostaticRZ'
        Pdata.iprop=1;
        Pdata.invflag = false;
        Pdata.RZflag = true;
        Pdata.cst = 1.129409067373019e+11;  % 1/Vacuum permittivity
    case 'thermal'
        Pdata.iprop=4;
        Pdata.invflag = false;
        Pdata.RZflag = false;
        Pdata.cst = 1.0;
    case 'thermalRZ'
        Pdata.iprop=4;
        Pdata.invflag = false;
        Pdata.RZflag = true;
        Pdata.cst = 1.0;
    case 'electrodynamicSS'
        Pdata.iprop=3;
        Pdata.invflag = false;
        Pdata.RZflag = false;
        Pdata.cst = 1.0;
    case 'electrodynamicSSRZ'
        Pdata.iprop=3;
        Pdata.invflag = false;
        Pdata.RZflag = true;
        Pdata.cst = 1.0;
end
Pdata.order = order;
Pdata.Ksize = KelSZ(order);

for iel = 1:nt
    pv = p(t(iel,:),:);
    PV = [ones(3,1), pv(1:3,:)];
    DetPV = det(PV);
    Area(iel) = 0.5 * DetPV;
    if Area(iel) < 0
        Area(iel) = abs(Area(iel));
        texc = t(iel,3);
        t(iel,3) = t(iel,2);
        t(iel,2) = texc;
        pv = p(t(iel,:),:);
        PV = [ones(3,1), pv(1:3,:)];
    end
    iPV = PV\eye(3); %inv(PV);
    GradL = iPV(2:3,:);
    gradL(iel,:,:) = GradL; 
    MatKind = materials(ireg(iel));
    Lg = ones(1,3)/3;
    varfield = 0;
    evalflg = true;
    derflag = false;
    [elprop(iel), ~, ~] = MatLib(MatKind, Pdata, 0, Lg, varfield, ifield(t(iel,:),:), evalflg, derflag);
                          
    refprop(t(iel,:))= refprop(t(iel,:)) + elprop(iel);
    icnt(t(iel,:)) = icnt(t(iel,:)) + 1;
end
refprop = refprop./icnt;
end

