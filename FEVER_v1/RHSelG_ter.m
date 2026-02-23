function [RHSeval] = RHSelG_ter(Lg, xyg, t, ireg, Pdata, ifield, order)


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[N] = NL_ter(Lg, order);
[sourceg, ~]  = t(Lg, xyg,  ifield, ireg, order);
if Pdata.RZflag    
    RHSeval = N * xyg(2) * sourceg; 
else
    RHSeval = N * sourceg;
end 
end

