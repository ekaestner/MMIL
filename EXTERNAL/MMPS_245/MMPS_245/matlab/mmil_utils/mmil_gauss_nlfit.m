function [beta,yfit,err] = mmil_gauss_nlfit(x,y,varargin)
%function [beta,yfit,err] = mmil_gauss_nlfit(x,y,[options])
%
% Purpose: use constrained nonlinear optimization to
%   to fit Gaussian function
%   of the form y = a * exp( -(x-b)^2 / (2*c^2)
%
% Required Input:
%   x: independent variable
%   y: dependent variable
%
% Optional Input:
%   'w': vector of weights applied to y-yfit
%     If empty, will weight all points equally
%     {default = []}
%   'lb': vector of lower bounds for 3 fit parameters
%     {default = []}
%   'ub': vector of upper bounds for 3 fit parameters
%     {default = []}
%   'maxiter': maximum number of iterations
%     {default = 500}
%   'tolfun': stop fitting when error changes less than this
%     {default = 1e-5}
%   'tolx': stop fitting when parameters change less than this
%     {default = 1e-6}
%   'thresh': relative threshold for initial linear fit
%     {default = 0.05}
%   'verbose': [0|1] display status messages
%     {default = 0}
%
% Output:
%   beta: vector of Gaussian fit parameters [a,b,c]
%   yfit: Gaussian fit of y
%   err: sum((y-yfit).^2)
%
% Created:  06/21/12 by Vijay Venkatraman
% Last Mod: 08/14/12 by Don Hagler
%

beta = []; yfit = []; err = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{,...
  'w',[],[],...
  'lb',[],[],...
  'ub',[],[],...
  'maxiter',500,[10,1e5],...
  'tolfun',1e-5,[],...
  'tolx',1e-6,[],...
  'thresh',0.05,[0,1],...
  'verbose',false,[false true],...
});

if isempty(parms.w)
  parms.w=ones(size(x));
end;
if isempty(parms.lb)
  parms.lb = [0.1*max(y),min(x),eps];
end;
if isempty(parms.ub)
  parms.ub = [1.1*max(y),max(x),max(x)];
end;

% initial parameters
beta_init = mmil_gauss_lfit(x,y,parms.thresh);

% fit options
fit_options = optimset('Display','off','Algorithm','active-set' ,...
  'MaxFunEvals',Inf,'MaxIter',parms.maxiter,...
  'TolFun',parms.tolfun,'TolX',parms.tolx);

[beta,err,exit_flag,output] = fmincon(@(p) mmil_gauss_err(x,y,p,parms.w),...
  beta_init,[],[],[],[],parms.lb,parms.ub,[],fit_options);

if any(isnan(beta))
  if parms.verbose
    fprintf('%s: WARNING: nonlinear fit failed:\n%s\n',...
      mfilename,output.message);
    fprintf('%s:   returning linear fit parameters\n',mfilename);
  end;
  beta = beta_init;
end;

yfit = mmil_gauss(x,beta);

return;

