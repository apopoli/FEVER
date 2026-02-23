function [F, Jac] = NRFunFEM(phi, flagJac, p, t, Pdata, Area, ireg, gradN, materiali, ifield, RHSg, BCflag_p, BCval_p, ~) %ultimo input omega che al momento non è utilizzato
Kg = zeros(size(p,1));
Jac = zeros(size(p,1));
np = size(p,1);
nt = size(t,1);
for iel=1:nt

    phiv = phi(t(iel,:));
    pv = p(t(iel,:),:);
    gradL(:,:) = gradN(iel,:,:);
    TipoMateriale = materiali(ireg(iel));
%     if flagJac
% %                       ElTri_nl( pv, phiv, Grad_Lel, Area, MatKind, Pdata, ifield, Jacflag)
%         [
% %
% %         [Kel, Jel1] = JacEl1(Area(iel), gradL, TipoMateriale, TipoProblema, omega, phiv);
%         
%         Jac(t(iel,:),t(iel,:)) = Jac(t(iel,:),t(iel,:)) + Jel(:,:);
%     else
%         ifield = 0;
%         [ Kel, Jel ] = ElTri_nl( pv, phip, grdL, Area, MatKind, Pdata, ireg, src, order, ifield );
% %        
% %         [Kel] = ElTri(Area(iel), gradL, TipoMateriale, TipoProblema, omega, ifield);
%     end

ifieldv = ifield(t(iel,:),:);

[Kel, Jel ] = ElTri_nl( pv, phiv, gradL, Area(iel), TipoMateriale, Pdata, ifieldv,  flagJac);
%assemblaggio
Jac(t(iel,:),t(iel,:)) = Jac(t(iel,:),t(iel,:)) + Jel(:,:);
Kg(t(iel,:),t(iel,:)) = Kg(t(iel,:),t(iel,:)) + Kel(:,:);
% iel

end

%Condizioni al contorno
for ip = 1:np
    switch BCflag_p(ip)
        case(1) %Neumann BC
            RHSg(ip) = RHSg(ip) + BCval_p(ip);
        case(2) %Dirichlet 
            Kg(ip,:) = 0.0;
            Kg(ip, ip) = 1.;
            if flagJac
                Jac(ip,:) = 0.0;
                Jac(ip, ip) = 1.;
            end
            
            RHSg(ip) = BCval_p(ip);
    end
end

F = Kg * phi-RHSg;

end

