function results = mmil_glm_calc(X,D,varargin)
%function results = mmil_glm_calc(X,D,[options])
%
% Purpose: calculate general linear model fit to data
%
% Usage:
%  results = mmil_glm_calc(X,D,'key1', value1,...);
%
% Required Input:
%  X: design matrix   size = [nsubs,nparams] (e.g. nsubs = # of subjects x # of conditions)
%  D: data matrix     size = [nsubs,nvals]   (e.g. nvals = # of values to be separately tested)
%
% Optional Input:
%  'singularity_flag': [0|1]  eliminate non-identifiable params from X
%     to remove potential singularities
%    {default = 0}
%  'contrast_vectors': cell array of contrast vectors
%    number of elements for each vector must match nparams
%    {default = []}
%  'contrast_names': cell array of contrast names
%    number of elements must match contrast_vectors
%    {default = []}
%  'nconds': number of conditions
%    {default = 1}
%   
% Output:
%   results: structure containing GLM betas, contrasts, etc.
%
% Created:  10/21/11 by Don Hagler
% Last Mod: 05/04/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: add to header once issues are resolved
%  'weights': weighting vector with size = [nsubs,1]
%     for weighted least squares regression
%     {default = []}
%  'reweight_flag' [0|1]: perform iteratively reweighted least squares (IRLS)
%    {default = 0}
%  'reweight_fact': reweighting factor (multiple of median absolute deviation)
%    {default = 4.685}
%  'reweight_maxiter': maximum number of reweighting iterations
%    {default = 100}
%  'reweight_tol': tolerance value for change in weights (stop iterations)
%    {default = 10^-7}
%  'reweight_leverage_flag': [0|1] for IRLS, adjust residuals using leverage
%    {default = 1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: calculate rank of X:
  % from statrobustfit:
  %[Q,R,perm] = qr(X,0);
  %tol = abs(R(1)) * max(n,p) * eps(class(R));
  %xrank = sum(abs(diag(R)) > tol);

%% todo: calculate robust sigma
%% todo: adjust p values for weighted least squares?

%% todo: r values are incorrect for multiple regression
%% todo: repeated measures design (off-diagonals in beta_cov?)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin, 1), return; end;

parms = check_input(X,D,varargin);
results = init_results(parms);

if ~parms.reweight_flag
  results = calc_results(X,D,results,parms);
else
  results = calc_results_IRLS(X,D,results,parms);
end;  

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(X,D,options)
  parms = mmil_args2parms(options,{...
    'weights',[],[],...
    'reweight_flag',false,[false true],...
    'reweight_fact',4.685,[0.1,100],...
    'reweight_maxiter',100,[1,1000],...
    'reweight_tol',1e-7,[1e-100,1e-2],...
    'reweight_leverage_flag',true,[false true],...
    'singularity_flag',false,[false true],...
    'contrast_vectors',[],[],...
    'contrast_names',[],[],...
    'nconds',1,[],...
  });

  % get info about design matrix
  parms.nsubs = size(X,1);
  parms.nparams = size(X,2);
  %% todo: calculate rank of design matrix? (see statrobustfit)
  % check for correspondence with data matrix
  if size(D,1) ~= parms.nsubs
    error('size of D does not match X');
  end;
  parms.nvals = size(D,2);
  % check weights vector
  if ~isempty(parms.weights)
    if numel(parms.weights)~=parms.nsubs
      error('number of elements in weights vector (%d) does not match data (%d)\n',...
        numel(parms.weights),parms.nsubs);
    end;
    parms.weights = reshape(parms.weights,[parms.nsubs,1]);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results.nsubs = parms.nsubs;
  results.nparams = parms.nparams;
  results.nvals = parms.nvals;
  results.weights = parms.weights;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_results(X,D,results,parms)
  % calculate inverse and estimate parameters
  results = fit_data(X,D,results,parms);
  % calculate summary statistics
  if isempty(results.weights)
    results = calc_stats(D,results,parms);
  else
    results = calc_stats(D0,results,parms);
  end;
  % calculate statistics for contrasts
  results = calc_contrasts(results,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = fit_data(X,D,results,parms)
  results.X = X;
  % apply weights to X and D
  if ~isempty(results.weights)
    [X,D] = apply_weights(X,D,results.weights);
  end;
  % calculate inverse matrix
  results.Xi = pinv(X);
  % eliminate non-identifiable params to remove potential singularities
  if parms.singularity_flag
    %% todo: see statrobustfit for different approach (calculate rank)
    XiX = results.Xi*X; % Crosstalk matrix
    xix = diag(XiX);
    results.Xi(find(xix<0.9),:) = 0; 
  end;
  % calculate parameter estimates
  results.betas = results.Xi*D;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_results_IRLS(X,D,results,parms)
  % calculate adjustment factor based on leverage
  if parms.reweight_leverage_flag
    results.leverage_fact = calc_leverage_fact(X);
  else
    results.leverage_fact = 1;
  end;
  D0 = D;
  results.X0 = X;
  if isempty(results.weights)
    results.weights = ones(results.nsubs,1);
  end;
  results.madsigma = std(D,0,1);
  % iteratively fit data recalculate weights until they converge
  for i=1:parms.reweight_maxiter
    % calculate inverse and estimate parameters
    results = fit_data(X,D,results,parms);
    % calculate weights from residuals
    [w,s] = calc_weights(X,D,results,parms);
    % check for convergence (difference in weights is less than tolerance)
    wd = results.weights - w;
    wd = mean(abs(wd(:)))/max(abs(results.weights(:)));
    if wd<parms.reweight_tol
%      fprintf('%s: difference in estimates is less than tolerance\n',mfilename);
      break;
    end;
    results.weights = w;
    results.madsigma = s;
  end;
  if i>=parms.reweight_maxiter
    fprintf('%s: WARNING: reweighted least squares did not converge within %d iterations\n',...
      mfilename,parms.reweight_maxiter);
  else
    fprintf('%s: reweighted least squares converged after %d iterations\n',...
      mfilename,i);
  end;
  % calculate summary statistics
  results = calc_stats(D,results,parms);
  % calculate statistics for contrasts
  results = calc_contrasts(results,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,D] = apply_weights(X,D,w)
  afun = @(a,b) a.*b;
  sw = sqrt(w);
  X = bsxfun(afun,X,sw);
  D = bsxfun(afun,D,sw);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function levrage_fact = calc_leverage_fact(X)
  % from Barret and Gray, 1997
  H = X*pinv(X'*X)*X';
  % from matlab statrobustfit
  h = min(.9999, sum(H.*H,2));
  levrage_fact = 1 ./ sqrt(1-h);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w,s] = calc_weights(X,D,results,parms)
  % calculate residual error
  r = D - X*results.betas;
  % collapse across values (e.g. vertices, voxels, ROIs, or measures)
  if results.nvals>1
    r = mean(abs(r),1);
    r = r - min(r);
  end;
  % adjust for "leverage"
  r = r .* results.leverage_fact;
  % set minimum value for standard deviation estimate
  tiny_s = 1e-6 * max(std(D,0,1));
  if tiny_s==0, tiny_s = 1; end
  % calculate median absolute deviation
  %% todo: use rank instead of nparams?
  s = max(madsigma(r,results.nparams),tiny_s);
  % normalize errors
  w = r/(s*parms.reweight_fact);
  % apply Tukey's bisquare weight function
  w = (abs(w)<1) .* (1 - w.^2).^2;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = madsigma(r,p)
  rs = sort(abs(r));
  s = median(rs(max(1,p):end)) / 0.6745;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_stats(D,results,parms)
  Dh = results.X*results.betas;
  results.ssq_data = sum(D.^2,1);
  results.var_data = var(D,0,1);
  results.ssq_fit = sum(Dh.^2,1);
  results.var_fit = var(Dh,0,1);
  results.ssq_err = sum((Dh-D).^2,1);
  results.var_err = var(Dh-D,0,1);
  results.exp_var = results.var_fit ./ results.var_data;
  if isempty(results.weights)
    results.nsubs_incl = results.nsubs;
  else
    results.nsubs_incl = sum(results.weights>0);
  end;
  results.dof = max(1,results.nsubs_incl - results.nparams);
  if results.nsubs_incl == results.nparams
    fprintf('%s: WARNING: same number of parameters as data points\n',mfilename);
  elseif results.nparams > results.nsubs_incl
    fprintf('%s: WARNING: more parameters than data points\n',mfilename);
  end;
  %% todo: robust sigma (see statrobustsigma private function in stats toolbox)
  results.sigma = sqrt(results.ssq_err/results.dof);
  results.beta_cov = results.Xi*results.Xi';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_contrasts(results,parms)
  results.contrasts = [];
  for i=1:length(parms.contrast_vectors)
    contrast_vec = parms.contrast_vectors{i};
    results.contrasts(i).name = parms.contrast_names{i};
    results.contrasts(i).vec = contrast_vec;
    sigma2_contrast = contrast_vec*results.beta_cov*contrast_vec';  
    sigma_contrast = sqrt(sigma2_contrast);
    stderr = sigma_contrast*results.sigma;
    beta_contrast = contrast_vec*results.betas;
    x = results.X * contrast_vec';
    sx = std(x);
    % if multiple conditions within subject, do not count subject multiple times
    %   divide by number of conditions, unless we are making contrast between conditions
    tmp_cont_vec = abs(contrast_vec);
    ind_nonzero = find(tmp_cont_vec);
    tmp_nconds = parms.nconds;
    tmp_nparams = results.nparams;
    if any(ismember(ind_nonzero,1+[1:parms.nconds])) % contrast of one or more conditions
      tmp_nconds = 1;
    else % not a condition regressor
      % don't count the extra conditions as free parameters if we are collapsing across
      %  shouldn't be penalized for averaging
      tmp_nparams = results.nparams - (parms.nconds-1);
    end;
    results.contrasts(i).n = results.nsubs_incl/tmp_nconds;
    results.contrasts(i).dof = max(1,results.contrasts(i).n-tmp_nparams);
    if tmp_nparams > results.contrasts(i).n
      fprintf('%s: WARNING: more parameters than data points for %s\n',...
        mfilename,results.contrasts(i).name);
    end;
    results.contrasts(i).mean = contrast_vec*results.betas;
    results.contrasts(i).stderr = stderr;
    results.contrasts(i).stdv = stderr * sqrt(results.contrasts(i).n);
    results.contrasts(i).mean(find(abs(results.contrasts(i).mean)<10*eps)) = 0;
  %  results.contrasts(i).ssq = diag(beta_contrast' * pinv(sigma2_contrast) * ...
  %                                  beta_contrast)';
    results.contrasts(i).ssq = (beta_contrast.^2)/sigma2_contrast;
    results.contrasts(i).r = results.contrasts(i).mean .* sx ./ ...
                             sqrt(results.var_data); % stdv of data
  %                           results.contrasts(i).stdv; % stdv of data
    results.contrasts(i).Fstat = results.contrasts(i).ssq ./ ...
                                 (results.ssq_err/results.dof);
    results.contrasts(i).tstat = results.contrasts(i).mean./(stderr+eps);
    % NOTE: tcdf is non-monotonic for some extreme values of -t -- fixed in R2010a(7.3)
    results.contrasts(i).pval = ...
      (results.contrasts(i).tstat<0).* (2*tcdf(results.contrasts(i).tstat,...
                                                results.contrasts(i).dof)) +...
      (results.contrasts(i).tstat>0).* (2*tcdf(-results.contrasts(i).tstat,...
                                                results.contrasts(i).dof));
    results.contrasts(i).log10_pval = ...
      (results.contrasts(i).tstat<0).* (log10(2*tcdf(results.contrasts(i).tstat,...
                                                      results.contrasts(i).dof))) +...
      (results.contrasts(i).tstat>0).*(-log10(2*tcdf(-results.contrasts(i).tstat,...
                                                      results.contrasts(i).dof)));
  end;
return;
