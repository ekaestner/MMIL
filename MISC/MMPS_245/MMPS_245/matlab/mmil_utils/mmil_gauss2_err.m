function err = mmil_gauss2_err(x1,x2,Y,beta,W)
%function err = mmil_gauss2_err(x1,x2,Y,beta,W)
%
% Purpose: calculate error of 2D Gaussian fit of the form
%   Y = a*exp(-(p*(x1-u1).^2 + 2*q*(x1-u1).*(x2-u2) + r*(x2-u2).^2))
%   with p = cos(t)^2/(2*s1^2) + sin(t)^2/(2*s2^2)
%        q = -sin(2*t)/(4*s1^2) + sin(2*t)/(4*s2^2)
%        r = sin(t)^2/(2*s1^2) + cos(t)^2/(2*s2^2)
%
% Required Input:
%   x1: vector of values for independent variable 1
%   x2: vector of values for independent variable 2
%   Y: dependent variable matrix with size [nx1,nx2]
%   beta: vector of Gaussian fit parameters [a,u1,s1,u2,s2,t]
%
% Optional Input:
%   W: matrix of weights applied to (Y-Yfit)
%
% Output:
%   err: sum of squared error
%
% Created:  08/13/12 by Don Hagler
% Last Mod: 08/13/12 by Don Hagler
%

err = [];
if ~mmil_check_nargs(nargin,4), return; end;
if ~exist('W','var') || isempty(W), W = ones(numel(x1),numel(x2)); end;

Yfit = mmil_gauss2(x1,x2,beta);

err = sum(mmil_rowvec((W.*(Y - Yfit)).^2));

return;

