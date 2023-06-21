function [Area, gradL, refprop, elprop, Pdata,t] = precheck_ter(p, t, ProbKind, ireg, materials, ifield)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nt = size(t,1);
np = size(p,1);
Area = zeros(nt,1);
elprop = zeros(nt,1);
refprop = zeros(np,1);
icnt = zeros(np,1);
gradL = zeros(nt,2,3);

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
    gradL(iel,:,:) = iPV(2:3,:); 
    MatKind = materials(ireg(iel));
    Lg = ones(1,3)/3;
    [elprop(iel), ~] = MatLib_ter(MatKind, Pdata.iprop, 0, Lg, ifield(t(iel,:),:), true, Pdata.invflag);

    refprop(t(iel,:))= refprop(t(iel,:)) + elprop(iel);
    icnt(t(iel,:)) = icnt(t(iel,:)) + 1;
end
refprop = refprop./icnt;
end

