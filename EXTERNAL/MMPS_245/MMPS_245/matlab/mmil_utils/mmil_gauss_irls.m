function [beta,yfit,w,err] = mmil_gauss_irls(x,y,varargin)
%function [beta,yfit,w,err] = mmil_gauss_irls(x,y,[options])
%
% Purpose: uses iteratively reweighted least squares
%   to fit a 1D Gaussian function
%   of the form y = a * exp( -(x-b)^2 / (2*c^2)
%
% Required Input:
%   x: independent variable
%   y: dependent variable
%
% Optional Input:
%   'lb': vector of lower bounds for 3 fit parameters
%     {default = []}
%   'ub': vector of upper bounds for 3 fit parameters
%     {default = []}
%   'maxiter': maximum number of iterations
%     {default = 100}
%   'errfact': scaling factor used to calculate weights from error
%     {default = 2}
%   'tol': fit tolerance
%     {default = 1e-7}
%   'verbose': [0|1] display status messages
%     {default = 0}
%
% Output:
%   beta: vector of Gaussian fit parameters [a,b,c]
%   w: vector of weights for each x-y pair
%
% Created:  08/12/12 by Don Hagler
% Last Mod: 08/13/12 by Don Hagler
%

beta = []; yfit = []; w = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{,...
  'lb',[],[],...
  'ub',[],[],...
  'maxiter',100,[10,1000],...
  'errfact',2,[0.1,10],...
  'tol',1e-7,[1e-20,1e-1],...
  'verbose',false,[false true],...
});

x = mmil_rowvec(x);
y = mmil_rowvec(y);

tiny_s = 1e-6 * max(abs(y));
if tiny_s==0, tiny_s = 1; end
ind_nonzero = find(y>tiny_s);
ind_zero = find(y<=tiny_s);

[beta,yfit] = mmil_gauss_nlfit(x,y,'lb',parms.lb,'ub',parms.ub);
ylast = yfit;

for i=1:parms.maxiter
  if parms.verbose
    fprintf('%s: reweight iter %d\n',mfilename,i);
  end;
  % calculate error
  err = abs(y - yfit);
  % calculate robust estimate of standard deviation using
  %   median absolute deviation
  s = max(median(err(ind_nonzero))/0.6745,tiny_s);
  % calculate weights from normalized error using
  %   Tukey's bisquare weight function
  w = err/(s*parms.errfact);
  w = sqrt((abs(w<1)) .* (1-w.^2).^2);
  w(ind_zero) = 0;
  % Gaussian fit with weighting factors
  [beta,yfit] = mmil_gauss_nlfit(x,y,,'w',w,'lb',parms.lb,'ub',parms.ub);
  % check for convergence
  d = mean(abs(yfit - ylast));
  if d < parms.tol
    if parms.verbose
      fprintf('%s: mean difference in fittted values is less than tolerance\n',...
        mfilename);
    end;
    break;
  end;
  ylast = yfit;
end;
if i>=parms.maxiter
  fprintf('%s: WARNING: irls did not converge within %d iterations\n',...
    mfilename,parms.maxiter);
elseif parms.verbose
  fprintf('%s: converged in %d iterations\n',mfilename,i);
end;

return;

