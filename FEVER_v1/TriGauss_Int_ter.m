function [int] = TriGauss_Int_ter(deg, f, pv, S, varargin)
%----------------------------------------------------------------------------
% Function TriGauss_Int compute the integral of the function f on a triangle  
% of area S using a Gauss quadrature rule of order p
% for the Gaussian quadrature for . %
% %
% On Input: 
% deg  = order of the Gauss Quadrature 
% f  = function to be integrated (must be given as a funcion of the 
%       barycentric coordinates &/or the cartesian coordinates)
% dmf = dimension of f. f is a coloumn array with dmf entries
% pv = vertex x & y coordinates
% S  = Area of the triangle on which the integral is calculated
% varargin = opional arguments (used by f) 
% %
% Output:
% Int = result of the integration
%
%   Author: A. Cristofolini 09/02/2017
%----------------------------------------------------------------------------

[ng, Lg, W] = TriGauss_P_W_ter(deg);    %load coordinates and weight

xyg = Lg * pv;       %x and y coordinated of the Gauss points
% fg = zeros([ng, dmf]);

int = 0;
for ig=1: ng
    int = int + W(ig)*f(Lg(ig,:), xyg(ig,:), varargin{:});
end

int = int*S;
end
    




