function [beta,Yfit,err] = mmil_gauss2_nlfit(x1,x2,Y,varargin)
%function [beta,Yfit,err] = mmil_gauss2_nlfit(x1,x2,Y,[options])
%
% Purpose: use constrained nonlinear optimization to
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
%   'beta_init': initial estimates
%     If empty, will set based on 1D Gaussian fits for each independent variable
%     {default = []}
%   'W': matrix of weights applied to (Y-Yfit)
%     If empty, will weight all points equally
%     {default = []}
%   'lb': vector of lower bounds for 6 fit parameters
%     {default = []}
%   'ub': vector of upper bounds for 6 fit parameters
%     {default = []}
%   'maxiter': maximum number of iterations
%     {default = 500}
%   'tolfun': stop fitting when error changes less than this
%     {default = 1e-5}
%   'tolx': stop fitting when parameters change less than this
%     {default = 1e-6}
%   'thresh': relative threshold for initial 1D linear fits
%     {default = 0.05}
%   'verbose': [0|1] display status messages
%     {default = 0}
%
% Output:
%   beta: vector of Gaussian fit parameters [a,u1,s1,u2,s2,t]
%   Yfit: Gaussian fit of y
%   err: sum((Y-Yfit).^2)
%
% Created:  08/13/12 by Don Hagler
% Last Mod: 08/14/12 by Don Hagler
%

beta = []; Yfit = []; err = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{,...
  'beta_init',[],[],...
  'W',[],[],...
  'lb',[],[],...
  'ub',[],[],...
  'maxiter',500,[10,1e5],...
  'tolfun',1e-5,[],...
  'tolx',1e-6,[],...
  'thresh',0.05,[0,1],...
  'verbose',false,[false true],...
...
  'nlfit_tags',{'maxiter','tolfun','tolx','thresh','verbose'},[],...
});

if isempty(parms.W)
  parms.W = ones(numel(x1),numel(x2));
end;
if isempty(parms.lb)
  parms.lb = [0.1*max(Y(:)),min(x1),1e-3*max(x1),min(x2),1e-3*max(x2),-pi/4];
end;
if isempty(parms.ub)
  parms.ub = [1.1*max(Y(:)),max(x1),0.5*max(x1),max(x2),0.5*max(x2),pi/4];
end;

if isempty(parms.beta_init)
  % use 1D Gaussian fits to get initial parameters
  [y1,ind1] = max(Y,[],2);
  [y2,ind2] = max(Y,[],1);
  w1_ind = sub2ind(size(Y),1:length(y1),ind1')';
  w2_ind = sub2ind(size(Y),ind2,1:length(y2))';
  w1 = parms.W(w1_ind);
  w2 = parms.W(w2_ind);
  args = mmil_parms2args(parms,parms.nlfit_tags);
  beta1 = mmil_gauss_nlfit(x1,y1,'w',w1,args{:});
  beta2 = mmil_gauss_nlfit(x2,y2,'w',w2,args{:});
  beta_init = [mean([beta1(1),beta2(1)]),beta1(2:3),beta2(2:3),0];
else
  beta_init = parms.beta_init;
end;

% fit options
if parms.verbose
  dispopt = 'iter';
else
  dispopt = 'off';
end;
fit_options = optimset('Display',dispopt,'Algorithm','active-set' ,...
  'MaxFunEvals',Inf,'MaxIter',parms.maxiter,...
  'TolFun',parms.tolfun,'TolX',parms.tolx);

[beta,err,exit_flag,output] = fmincon(@(p) mmil_gauss2_err(x1,x2,Y,p,parms.W),...
  beta_init,[],[],[],[],parms.lb,parms.ub,[],fit_options);

if any(isnan(beta))
  fprintf('%s: WARNING: returning initial fit parameters',mfilename);
  if parms.verbose
    fprintf(':\n%s\n',output.message);
  else
    fprintf('\n');
  end;
  beta = beta_init;
end;

Yfit = mmil_gauss2(x1,x2,beta);

return;

