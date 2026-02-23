function [N] = NL(L, order)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
switch order
    case 1
      N = [L(1); L(2); L(3)];
    case 2
      N = [L(1)*(2*L(1)-1);
           L(2)*(2*L(2)-1); 
           L(3)*(2*L(3)-1); 
           4*L(1)*L(2); 
           4*L(2)*L(3); 
           4*L(3)*L(1)];  
end
end

