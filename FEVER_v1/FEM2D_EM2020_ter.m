function [ phi, field, Kg] = FEM2D_EM2020( ProbKind, p, t, bedges, bcflag_e, bcval_e, ireg, ibe, materials, source, order, ifield )
%FEM2D_EN15 Finite element solution of div(prop(x,y) grad(phi)) = source(x,y)
% material(id_phys_reg) = MaterialKind assigned to the region id_phys_reg

eps0 = 8.854e-12;
mu0 = 4*pi*1.e-7;
np = size(p,1);
nt = size(t,1);
nbedges = size(bedges,1);
ndom = length(ireg);




% elprop = zeros(nt,1);

src = zeros(nt,1);
IntSource = zeros(ndom,1);
grdL = zeros(2,3);
% Field_e = zeros(nt,2);

BCflag_p = zeros(np,1); %initialize B.C. flag at the points
BCval_p = zeros(np,1);  %initialize B.C. value at the points

Kg = zeros(np); %spalloc(np,np,ncf);    %initialize the coefficient matric Kg
RHSg = zeros(np,1);   %spalloc(np,1,np);    %initialize the right hand side array


[Area, gradL, refprop, elprop, Pdata,t] = precheck_ter(p, t, ProbKind, ireg, materials, ifield);
for iedge = 1 : nbedges
    if bcflag_e(iedge) ~= 3 %BC Neumann & Robin
        
        pv = p(bedges(iedge, :),:);     %load in pv the coordinates of the two
        %nodes of the edge (iedge)
        MatKind = materials(ibe(iedge));
        BCval = bcval_e(iedge,:);
        ifielde = ifield(bedges(iedge,:),:);
        [RHS1D, Kel1D] = OneDEl21_ter(pv, MatKind, Pdata, order, bcflag_e(iedge), BCval, ifielde);
        RHSg(bedges(iedge, :)) = RHSg(bedges(iedge, :)) +  RHS1D(:);        
        Kg(bedges(iedge, :), bedges(iedge, :)) = Kg(bedges(iedge, :), bedges(iedge, :)) +  Kel1D(:, :);
               
        BCflag_p(bedges(iedge, :)) = bcflag_e(iedge);
        BCval_p(bedges(iedge, :)) = BCval_p(bedges(iedge, :)) + RHS1D;
    end
end

for iedge = 1:nbedges   % ciclo per C.C. di Dirichlet
    if bcflag_e(iedge) == 3
        for iv =1:2
            if BCflag_p(bedges(iedge,iv)) == 1 || BCflag_p(bedges(iedge,iv)) == 2
                BCflag_p(bedges(iedge,iv)) = 3;
                BCval_p(bedges(iedge,iv)) = bcval_e(iedge, 1);
            else
                BCflag_p(bedges(iedge,iv)) = 3;
                BCval_p(bedges(iedge,iv)) = BCval_p(bedges(iedge,iv)) + 0.5*bcval_e(iedge, 1);
            end
        end
        if order == 2
            BCflag_p(bedges(iedge,3)) = 3;
            BCval_p(bedges(iedge,3)) = bcval_e(iedge);
        end
    end
end





for iel=1:nt
    pv = p(t(iel,:),:);     %load in pv the coordinates of the three
    %nodes of the triangle (iel)
    MatKind = materials(ireg(iel));
    
    grdL(:,:) = gradL(iel,:,:);
    ifieldt = ifield(t(iel,:),:);
    [Kel, RHSel] = ElTri21_ter( pv, grdL, Area(iel), MatKind, Pdata, ireg(iel), source, order, ifieldt ); %find the element matrix
   
    %and the element array (and other useful stuff).
    
    Kg(t(iel,:), t(iel,:)) = Kg(t(iel,:), t(iel,:)) + Kel(:,:); %assembla
    %the element matrix Kel in the global matrix K
    RHSg(t(iel,:)) = RHSg(t(iel,:)) +  RHSel(:); %assemble
    %the element rhs RHSel in the global rhs RHSg
end

Kg = Kg./refprop;
RHSg = RHSg./refprop;

%apply the boundary conditions
for ip = 1:np
    switch BCflag_p(ip)
        case(1) %Neumann BC
%             RHSg(ip) = RHSg(ip)/refprop(ip);
        case(3)
            Kg(ip,:) = 0.0;
            Kg(ip, ip) = 1.;
            RHSg(ip) = BCval_p(ip);
    end
end

phi = Kg\RHSg;

gradphi_e = zeros(nt,2);
gradphiprop = zeros(np,2);
grdphi = zeros(np,2);
icount = zeros(np,1);

% energy = zeros(nt,1);

%  phi = zeros(np,1);
%  phi(:) = phis;

%calculate gradient
%  gradphi_e(iel,:) = field gradient components on a generic element iel.
%  this is obtained multipling the shape function gradient matrix gradN(iel,:,:)
%  of the element iel by the nodal values of the phi function.

%  gradphi(ip,:) = field gradient components on a generic element iel. This is
%  obtained as an average of the field gradient components of all the
%  elements having the node ip as vertex

for iel = 1 : nt
    
    grdL(:,:) = gradL(iel,:,:);
%     energy(iel) = CompEnergy(phi(t(iel,:),:), prop_e(iel), Area(iel), grdL, order);
    
    gradphi_e(iel,:) = grdL * phi(t(iel,1:3));
    IntSource(ireg(iel)) = IntSource(ireg(iel)) + src(iel)*Area(iel);
    for ip = 1 : 3
 
        grdphi(t(iel,ip),:) = grdphi(t(iel,ip),:) + gradphi_e(iel,:);
        gradphiprop(t(iel,ip),:) = gradphiprop(t(iel,ip),:) + gradphi_e(iel,:)*elprop(iel);
        icount(t(iel,ip)) = icount(t(iel,ip)) + 1;
    end
end



grdphi(:,1) = grdphi(:,1)./icount;
grdphi(:,2) = grdphi(:,2)./icount;

gradphiprop(:,1) = gradphiprop(:,1)./icount;
gradphiprop(:,2) = gradphiprop(:,2)./icount;


field = zeros(np,3,5);
switch ProbKind
    case {'electrostatic', 'electrostaticRZ'}
        field(:,1,1) = -grdphi(:,1);
        field(:,2,1) = -grdphi(:,2);
        field(:,1,2) = -gradphiprop(:,1)*eps0;
        field(:,2,2) = -gradphiprop(:,2)*eps0;
%         energy = energy * eps0; 
        
    case 'magnetostatic'
        field(:,1,3) = grdphi(:,2);
        field(:,2,3) = -grdphi(:,1);
        field(:,1,4) = gradphiprop(:,2)/mu0;
        field(:,2,4) = -gradphiprop(:,1)/mu0;
%         energy = energy/mu0; 
    case {'electrodynamicSS', 'electrodynamicSSRZ'}
        field(:,1,1) = -grdphi(:,1);
        field(:,2,1) = -grdphi(:,2);
        field(:,1,5) = -gradphiprop(:,1);
        field(:,2,5) = -gradphiprop(:,2);
end

%calculate the total field energy

% gradmod2 = (sum(gradphi_e.^2,2));
% W = 0.5 * prop_e .* gradmod2/mu0;
% energy = Area' * W;

% Energy = sum(energy);

