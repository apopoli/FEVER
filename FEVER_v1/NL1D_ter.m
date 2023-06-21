function [N] = NL1D_ter(L, order)
%1D shape functions 
switch order
    case 1
      N = L;
    case 2
      N = [(L(:,1).*(2*L(:,1)-1)),(L(:,2).*(2*L(:,2)-1)),(4*L(:,1).*L(:,2))];             
             
end
end

