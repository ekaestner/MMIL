function [beta,Yfit,W,err] = mmil_gauss2_irls(x1,x2,Y,varargin)
%function [beta,Yfit,W,err] = mmil_gauss2_irls(x1,x2,Y,[options])
%
% Purpose: uses iteratively reweighted least squares
%   to fit Gaussian function of the form
%   Y = a*exp(-(p*(x1-u1).^2 + 2*q*(x1-u1).*(x2-u2) + r*(x2-u2).^2))
%   with p = cos(t)^2/(2*s1^2) + sin(t)^2/(2*s2^2)
%        q = -sin(2*t)/(4*s1^2) + sin(2*t)/(4*s2^2)
%        r = sin(t)^2/(2*s1^2) + cos(t)^2/(2*s2^2)
%
% Required Input:
%   x1: vector of values for independent variable 1
%   x2: vector of values for independent variable 2
%   Y: dependent variable matrix with size [nx1,nx2]
%
% Optional Input:
%   'lb': vector of lower bounds for 6 fit parameters
%     {default = []}
%   'ub': vector of upper bounds for 6 fit parameters
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
%   beta: vector of Gaussian fit parameters [a,u1,s1,u2,s2,t]
%   W: vector of weights for each x-y pair
%
% Created:  08/13/12 by Don Hagler
% Last Mod: 08/13/12 by Don Hagler
%

beta = []; Yfit = []; w = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{,...
  'lb',[],[],...
  'ub',[],[],...
  'maxiter',100,[10,1000],...
  'errfact',2,[0.1,10],...
  'tol',1e-7,[1e-20,1e-1],...
  'verbose',false,[false true],...
});

x1 = mmil_rowvec(x1);
x2 = mmil_rowvec(x2);

tiny_s = 1e-6 * max(abs(Y(:)));
if tiny_s==0, tiny_s = 1; end
ind_nonzero = find(Y>tiny_s);
ind_zero = find(Y<=tiny_s);

[beta,Yfit] = mmil_gauss2_nlfit(x1,x2,Y,'lb',parms.lb,'ub',parms.ub);
Ylast = Yfit;

for i=1:parms.maxiter
  if parms.verbose
    fprintf('%s: reweight iter %d\n',mfilename,i);
  end;
  % calculate error
  E = abs(Y - Yfit);
  % calculate robust estimate of standard deviation using
  %   median absolute deviation
  s = max(median(E(ind_nonzero))/0.6745,tiny_s);
  % calculate weights from normalized error using
  %   Tukey's bisquare weight function
  W = E/(s*parms.errfact);
  W = sqrt((abs(W<1)) .* (1-W.^2).^2);
  W(ind_zero) = 0;
  % Gaussian fit with weighting factors
  [beta,Yfit] = mmil_gauss2_nlfit(x1,x2,Y,'W',W,'lb',parms.lb,'ub',parms.ub);
  % check for convergence
  d = mean(abs(Yfit - Ylast));
  if d < parms.tol
    if parms.verbose
      fprintf('%s: mean difference in fittted values is less than tolerance\n',...
        mfilename);
    end;
    break;
  end;
  Ylast = Yfit;
end;
if i>=parms.maxiter
  fprintf('%s: WARNING: irls did not converge within %d iterations\n',...
    mfilename,parms.maxiter);
elseif parms.verbose
  fprintf('%s: converged in %d iterations\n',mfilename,i);
end;

return;

