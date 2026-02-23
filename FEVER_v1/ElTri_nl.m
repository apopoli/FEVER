function [ Kel, Jel ] = ElTri_nl( pv, phiv, Grad_Lel, Area, MatKind, Pdata, ifield, flagJac)
%ELTRI calculation of Kel and Jel for a triangular element
%Jel = element Jacobian matrix (it is assumed that non linearity resides
%only in the material property defining the coefficient matrix and not in
%the source term)

un_terzo = 0.3333333333333333;
order = Pdata.order;
% [Mprop, Mdeg, dprop_dfield] = MatLib(MatKind, iprop, xyg, Lg, GrdN, varfield, ifield, evalflg, invflag, derflag)
% retrieving the polynomia order of the Mprop function
Lg = [un_terzo un_terzo un_terzo];     %dummy arg
GrdN = 0;       %dummy arg
varfield =0;    %dummy arg
% varfield = phip; 
xyg = sum(pv,2)*un_terzo;
evalflg =  false; 
derflag = false; 
               
[~, Mdeg, ~] = MatLib(MatKind, Pdata, xyg, Lg, varfield, ifield, evalflg, derflag);  

if Mdeg < 0
    Mdeg = order-1;
end
if Pdata.RZflag
    deg = Mdeg + 2*(order - 1) + 1;
else
    deg = Mdeg + 2*(order - 1);
end
gradl2 = Grad_Lel'*Grad_Lel;

% [int] = TriGauss_Int(deg, f, pv, S, varargin) varargin : MatKind, Pdata, Grad_Lel, gradl2, phiv,  ifield

KJel = TriGauss_Int(deg, @KJelG, pv, Area, MatKind, Pdata, Grad_Lel, gradl2, phiv,  ifield, flagJac);
sizeKel = Pdata.Ksize;
Kel = KJel(1:sizeKel,:);
if flagJac
    Jel = KJel(sizeKel + 1:end,:);
    
    Jel = Jel + Kel;
else
    Jel = Kel;
end
end

