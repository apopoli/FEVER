function [ phi, field, IntSource, iter, err ] = FEM2DNLv1( TipoProblema, p, t, bedges, bcflag_e, bcval_e, ireg, iregbe, materiali, sorgente, frequenza, ifield, order)
%FEM2D_EN15 Finite element solution of div(prop(x,y) grad(phi)) = source(x,y)
% material(id_phys_reg) = MaterialKind assigned to the region id_phys_reg

eps0 = 8.8541878128e-12;
mu0 = 4*pi*1.e-7;
np = size(p,1);
nt = size(t,1);
nbedges = size(bedges,1);
ndom = length(unique(ireg));

prop_e = zeros(nt,1);
sigma = zeros(nt,1);
src = zeros(nt,1);  %sorgente sull'elemento
IntSource = zeros(ndom,1);
omega = 2*pi*frequenza;

BCflag_p = zeros(np,1); %initialize B.C. flag at the points
BCval_p = zeros(np,1);  %initialize B.C. value at the points

[Area, GradL, refprop, Pdata, t] = precheck(p, t, TipoProblema, ireg, materiali, ifield, order);

RHS0 = zeros(np,1);      %initialize the right hand side array
% K0 = zeros(np);

for iedge = 1 : nbedges
    if bcflag_e(iedge) == 1 %BC Neumann & Robin
        
        pv = p(bedges(iedge, :),:);     %load in pv the coordinates of the two
                                        %nodes of the edge (iedge)
                                        
          %Solo BC di Neumann omogenee
          RHS1D = [0; 0];
          
%         MatKind = materials(ibe(iedge));
%         BCval = bcval_e(iedge,:);
%         ifielde = ifield(bedges(iedge,:),:);
%         [RHS1D, Kel1D] = OneDEl21(pv, MatKind, Pdata, order, bcflag_e(iedge), BCval, ifielde);
%         RHS0(bedges(iedge, :)) = RHS0(bedges(iedge, :)) +  RHS1D(:);        
%         K0(bedges(iedge, :), bedges(iedge, :)) = K0(bedges(iedge, :), bedges(iedge, :)) +  Kel1D(:, :);
        
        BCflag_p(bedges(iedge, :)) = bcflag_e(iedge);
        BCval_p(bedges(iedge, :)) = BCval_p(bedges(iedge, :)) + RHS1D(:,:);
    end
end

for iedge = 1:nbedges   % ciclo per C.C. di Dirichlet
    if bcflag_e(iedge) == 2
        for iv =1:2
            if BCflag_p(bedges(iedge,iv)) == 1
                BCflag_p(bedges(iedge,iv)) = 2;
                BCval_p(bedges(iedge,iv)) = bcval_e(iedge, 1);
            else
                BCflag_p(bedges(iedge,iv)) = 2;
                BCval_p(bedges(iedge,iv)) = BCval_p(bedges(iedge,iv)) + 0.5*bcval_e(iedge, 1);
            end
        end
        if order == 2
            BCflag_p(bedges(iedge,3)) = 2;
            BCval_p(bedges(iedge,3)) = bcval_e(iedge);
        end
    end
end


for iel=1:nt
    pv = p(t(iel,:),:);

%   Calcolo del termine noto

    [RHSel] = ElTriRHS( pv, Area(iel), Pdata, ireg(iel), sorgente, ifield(t(iel,:),:));
    RHS0(t(iel,:)) = RHS0(t(iel,:)) + RHSel(:);   
end

% phi0 = K0\RHS0; %soluzione di primo tentativo valutata come soluzione con proprietà a campo nullo

phi0 = zeros(np,1);
tol = 1e-8;
itermax = 30;
% p, t, Area, ireg, GradL, materiali, RHS0, BCflag_p, BCval_p, omega);
% p, t, Area, ireg, gradN, materiali, RHSg, BCflag_p, BCval_p, ~)
[phi, iter, err]= NR_finale1(@NRFunFEM, phi0 , tol, itermax, p, t, Pdata, Area, ireg, GradL, materiali, ifield, RHS0, BCflag_p, BCval_p, omega);

gradphi_e = zeros(nt,2);
src_p = zeros(np,1);
gradphiprop = zeros(np,2);
grdphi = zeros(np,2);
icount = zeros(np,1);
icounts = zeros(np,1);
grdn = zeros(2,3);

%  phi = zeros(np,1);
%  phi(:) = phis;

%calculate gradient
%  gradphi_e(iel,:) = field gradient components on a generic element iel.
%  this is obtained multipling the shape function gradient matrix gradN(iel,:,:)
%  of the element iel by the nodal values of the phi function.

%  gradphi(ip,:) = field gradient components on a generic element iel. This is
%  obtained as an average of the field gradient components of all the
%  elements having the node ip as vertex

Lg = [1 1 1]/3;
evalflg = true;
derflag = false;
xyg = 0;
for iel = 1 : nt
    
    grdn(:,:) = GradL(iel,:,:);    
    gradphi_e(iel,:) = grdn * phi(t(iel,:));
    ifieldg = ifield(t(iel,:),:);
    TipoMateriale = materiali(ireg(iel));
    varfield = norm(gradphi_e(iel,:));
    [prop_e(iel), ~, ~] = MatLib(TipoMateriale, Pdata, xyg, Lg, varfield, ifieldg, evalflg, derflag);

    if TipoProblema == "QmagnetostaticoSin"
        ipsave = Pdata.iprop;
        Pdata.iprop = 3;
        [sigma(iel), ~, ~] = MatLib(TipoMateriale, Pdata, xyg, Lg, varfield, ifield, evalflg, derflag);
        Pdata.iprop = ipsave;
    end
    
    src(iel) = sorgente(Lg, xyg,  ifield, ireg(iel), Pdata.order) * Area(iel);
    %ATTENZIONE! Prima questa parte c'era, ma forse non è essenziale e deve
    %essere corretta
   % IntSource(ireg(iel)) = IntSource(ireg(iel)) + src(iel);
    for ip = 1 : 3
        if TipoProblema == "QmagnetostaticoSin"
            indcd_curr_density = - 1i*omega*sigma(iel)*phi(t(iel,ip));
            src_p(t(iel,ip)) = src_p(t(iel,ip)) + indcd_curr_density;
            IntSource(ireg(iel)) = IntSource(ireg(iel)) + indcd_curr_density*Area(iel)/3;
        end
        src_p(t(iel,ip)) = src_p(t(iel,ip)) + src(iel);
        if src_p(t(iel,ip)) ~= 0
             icounts(t(iel,ip)) = icounts(t(iel,ip)) + 1;
        end
        grdphi(t(iel,ip),:) = grdphi(t(iel,ip),:) + gradphi_e(iel,:);
        gradphiprop(t(iel,ip),:) = gradphiprop(t(iel,ip),:) + gradphi_e(iel,:)*prop_e(iel);
        icount(t(iel,ip)) = icount(t(iel,ip)) + 1;
    end
end

grdphi(:,1) = grdphi(:,1)./icount;
grdphi(:,2) = grdphi(:,2)./icount;

gradphiprop(:,1) = gradphiprop(:,1)./icount;
gradphiprop(:,2) = gradphiprop(:,2)./icount;

ips = find(src_p ~= 0);
src_p(ips) = src_p(ips)./icounts(ips);

field = zeros(np,3,5);
switch TipoProblema
    case 'elettrostatico'
        field(:,1,1) = -grdphi(:,1);
        field(:,2,1) = -grdphi(:,2);
        field(:,1,2) = -gradphiprop(:,1)*eps0;
        field(:,2,2) = -gradphiprop(:,2)*eps0;
    case 'magnetostatico'
        field(:,1,3) = grdphi(:,2);
        field(:,2,3) = -grdphi(:,1);
        field(:,1,4) = gradphiprop(:,2)/mu0;
        field(:,2,4) = -gradphiprop(:,1)/mu0;
        field(:,3,5) = src_p;
    case 'QmagnetostaticoSin'
        field(:,1,3) = grdphi(:,2);
        field(:,2,3) = -grdphi(:,1);
        field(:,1,4) = gradphiprop(:,2)/mu0;
        field(:,2,4) = -gradphiprop(:,1)/mu0;
        field(:,3,5) = src_p;
    case 'electrodynamicSSRZ'
        field(:,1,1) = -grdphi(:,1);
        field(:,2,1) = -grdphi(:,2);
        field(:,1,5) = -gradphiprop(:,1);
        field(:,2,5) = -gradphiprop(:,2);        
end



end

