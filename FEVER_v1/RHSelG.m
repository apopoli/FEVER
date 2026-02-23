function [RHSeval] = RHSelG(Lg, xyg, src, ireg, Pdata, ifield)


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[N] = NL(Lg, Pdata.order);
[sourceg, ~]  = src(Lg, xyg,  ifield, ireg, Pdata.order);
if Pdata.RZflag    
    RHSeval = N * xyg(2) * sourceg; 
else
    RHSeval = N * sourceg;
end 
end

