function [beta,yfit,err] = mmil_gauss_lfit(x,y,thresh)
%function [beta,yfit,err] = mmil_gauss_lfit(x,y,thresh)
%
% Purpose: fit a 1D Gaussian function using least squares
%   uses log transform of y to use linear polynomial fit
%   of the form y = a * exp( -(x-b)^2 / (2*c^2)
%
% Required Input:
%   x: independent variable
%   y: dependent variable
%
% Optional Input:
%   thresh: relative threshold applied to y
%     to exclude small values from fit
%     {default = 0.1}
%
% Output:
%   beta: vector of Gaussian fit parameters [a,b,c]
%   yfit: Gaussian fit of y
%   err: sqrt(sum((y-yfit).^2))
%
% Created:  06/21/12 by Vijay Venkatraman
% Last Mod: 08/12/12 by Don Hagler
%

% informed by:
% http://terpconnect.umd.edu/~toh/spectrum/CurveFitting.html

beta = []; yfit = []; err = [];
if ~exist('thresh','var') || isempty(thresh), thresh = 0.1; end;

x = mmil_rowvec(x);
y = mmil_rowvec(y);

x0 = x;
y0 = y;

if thresh
  thresh = thresh*max(y);
  x = x(y>=thresh);
  y = y(y>=thresh);
end;

y=log(y);
coef=polyfit(x,y,2);
a=coef(3);
b=coef(2);
c=coef(1);
beta(1) = exp(a-c*(b/(2*c))^2); % A (magnitude)
beta(2) = -b/(2*c); % B (mean)
beta(3) = 1/(sqrt(2)*sqrt(-c)); % C (variance)

yfit = mmil_gauss(x0,beta);
err = mmil_gauss_err(x0,y0,beta);

return;

