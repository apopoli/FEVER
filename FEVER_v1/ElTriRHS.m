function [ RHSel ] = ElTriRHS( pv, Area, Pdata, ireg, src, ifield )
%Evaluate the element right hand side

[~, tdeg] = src(0, 0, 0, ireg, Pdata.order);
if Pdata.RZflag
    deg = (1 + tdeg) * Pdata.order + 1;
else
    deg = (1 + tdeg) * Pdata.order;
end

RHSel = Pdata.cst * TriGauss_Int(deg, @RHSelG, pv, Area, src, ireg, Pdata, ifield);


end

