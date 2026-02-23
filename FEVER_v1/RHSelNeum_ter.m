function [val] = RHSelNeum_ter(Lcsi, pv, MatKind, order, Pdata, ifield)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   pflag = problem flag 
% RZflag = false --> problema piano
% RZflag = true --> problema RZ

[N] = NL1D_ter(Lcsi, order);
evalflg = true;
pg = Lcsi * pv;
[Mprop, ~] = MatLib_ter(MatKind, Pdata.iprop, pg, Lcsi, ifield, evalflg, Pdata.invflag);

if Pdata.RZflag 
    rg = pg(2);
    val = N' * rg * Mprop;   
else
    val = N' * Mprop;
end


end

