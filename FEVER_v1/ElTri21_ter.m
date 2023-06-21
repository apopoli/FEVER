function [ Kel, RHSel ] = ElTri21_ter( pv, Grad_Lel, Area, MatKind, Pdata, ireg, src, order, ifield )
%ELTRI calculation of Kel and RHSel for a triangular element



[~, Mdeg] = MatLib_ter(MatKind, Pdata.iprop, [0 0], 0, 0, false, false);  %ask for the polynomial degree ofthe function 
if Mdeg < 0
    Mdeg = abs(Mdeg)*order;
end
if Pdata.RZflag
    deg = Mdeg + order + 1;
else
    deg = Mdeg + order;
end
gradl2 = Grad_Lel'*Grad_Lel;
Kel = TriGauss_Int_ter(deg, @KelG_ter, pv(1:3,:), Area, MatKind, Pdata, gradl2, ifield, order);

[~, tdeg] = src(0, 0, 0, ireg, order);
if Pdata.RZflag
    deg = (1+ tdeg) * order + 1;
else
    deg = (1+ tdeg) * order;
end

RHSel = Pdata.cst * TriGauss_Int_ter(deg, @RHSelG_ter, pv(1:3,:), Area, src, ireg, Pdata, ifield, order);


end

