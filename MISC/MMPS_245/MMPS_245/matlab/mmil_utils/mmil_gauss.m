function y = mmil_gauss(x,beta)
%function y = mmil_gauss(x,beta)
%
% Purpose: 1D Gaussian function
%   of the form y = a * exp(-(x-b)^2 / (2*c^2)))
%
% Required Input:
%   x: independent variable
%   beta: vector of Gaussian fit parameters [a,b,c]
%
% Output:
%   y: dependent variable
%
% Created:  08/12/12 by Don Hagler
% Last Mod: 08/12/12 by Don Hagler
%


y = [];
if ~mmil_check_nargs(nargin,2), return; end;

y = beta(1)*exp(-(x-beta(2)).^2/(2*beta(3)^2));

return;

