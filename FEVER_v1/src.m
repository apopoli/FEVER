function [srcval, srcdeg] = src(L, ~, ~, ireg, ~)
% src(Lg, xyg,  ifield, ireg, order)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
switch ireg
    case 1
        srcval = 0;%1e6;%
        srcdeg = 0;
    case 2
        srcval = 0;%1e6;%
        srcdeg = 0;
        %         = 1/(1+xy(2))^2;
    case 3
        srcval = 0;%1e6;%
        srcdeg = 0;
        %         = 1/(1+xy(2))^2;
       otherwise
        srcval = 0;
        srcdeg = 0;
end

