function Y = mmil_gauss2(x1,x2,beta)
%function Y = mmil_gauss2(x1,x2,beta)
%
% Purpose: evaluate 2D Gaussian function of the form
%   Y = a*exp(-(p*(x1-u1).^2 + 2*q*(x1-u1).*(x2-u2) + r*(x2-u2).^2))
%   with p = cos(t)^2/(2*s1^2) + sin(t)^2/(2*s2^2)
%        q = -sin(2*t)/(4*s1^2) + sin(2*t)/(4*s2^2)
%        r = sin(t)^2/(2*s1^2) + cos(t)^2/(2*s2^2)
%
% Required Input:
%   x1: vector of values for independent variable 1
%   x2: vector of values for independent variable 2
%   beta: vector of Gaussian fit parameters [a,u1,s1,u2,s2,t]
%
% Output:
%   Y: output matrix with size [nx1,nx2]
%
% Created:  08/13/12 by Don Hagler
% Last Mod: 08/14/12 by Don Hagler
%

Y = [];
if ~mmil_check_nargs(nargin,3), return; end;

x1 = mmil_rowvec(x1);
x2 = mmil_rowvec(x2);

% create grid of x1 and x2 values
[X2,X1] = meshgrid(x2,x1);

Y = mmil_gauss2_eval(X1,X2,beta);

return;

