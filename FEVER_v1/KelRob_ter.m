function [val] = KelRob(Lcsi, p, order, Pdata)
%UNTITLED Summary of this function goes here
%   pflag = problem flag 
% RZflag = false --> problema piano
% RZflag = true --> problema RZ

[N] = NL1D_ter(Lcsi, order);

if Pdata.RZflag 
    rg = Lcsi * p(:,2);
    val = rg * (N' * N);
else
   val = (N' * N);
end

end

