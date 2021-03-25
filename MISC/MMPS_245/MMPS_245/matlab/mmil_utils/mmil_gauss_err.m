function err = mmil_gauss_err(x,y,beta,w)
%function err = mmil_gauss_err(x,y,beta,w)
%
% Purpose: calculate error of 1D Gaussian fit
%   of the form y = a * exp( -(x-b)^2 / (2*c^2)
%
% Required Input:
%   x: independent variable
%   y: dependent variable
%   beta: vector of Gaussian fit parameters [a,b,c]
%
% Optional Input:
%   w: vector of weights applied to (y-yfit)
%
% Output:
%   err: sum of squared error
%
% Created:  08/12/12 by Don Hagler
% Last Mod: 08/13/12 by Don Hagler
%

err = [];
if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('w','var') || isempty(w), w = ones(size(x)); end;

x = mmil_rowvec(x);
y = mmil_rowvec(y);
w = mmil_rowvec(w);

yfit = mmil_gauss(x,beta);

err = sum(mmil_rowvec((w.*(y - yfit)).^2));

return;
