function [val] = RHSelRob_ter(Lcsi, p, order, Pdata)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   pflag = problem flag 
% RZflag = false --> problema piano
% RZflag = true --> problema RZ

[N] = NL1D_ter(Lcsi, order);


if Pdata.RZflag 
    rg = Lcsi * p(:,2);
    val = N' * rg;   
else
    val = N';
end


end

