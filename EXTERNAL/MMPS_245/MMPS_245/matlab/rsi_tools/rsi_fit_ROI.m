function results = rsi_fit_ROI(vals,qmat,varargin)
%function results = rsi_fit_ROI(vals,qmat,[options])
%
% Purpose: perform nonlinear fit on multi-shell diffusion MRI data for ROI
%   with multi-compartment difusion models
% 
% Required input:
%   vals: 2D volume with size = [nvox,ndirs]
%   qmat: matrix of diffusion direction vectors with size = [ndirs,3]
%
% Optional parameters:
%   'bvals': vector of b values
%     one for all, or one for each diffusion direction
%     {default = 1000}
%   'nob0_flag': [0|1] whether to exclude b=0 images from fits
%     {default = 0}
%   'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit for DT and MFI fitting
%     {default = 1}
%   'mean_flag': [0|1] fit to data averaged across voxels
%     {default = 0}
%   'fit_ROI_flag': [0|1] perform nonlinear fitting for parameters
%     such as ADC_long, ADC_free, etc.
%     {default = 1}
%   'verbose': [0|1] display frequent status updates
%     {default = 1}
%
% Optional parameters for ROI trimming:
%   'trim_flag': [0|1] exclude outlier voxels from ROI based on
%     FA, MD, and T2 calculated from tensor fit
%     {default = 0}
%   'trim_fact': scaling factor applied to median absolute deviations
%     {default = 4.7}
%   'trim_thresh': threshold applied to weighting factors for trimming
%     {default = 0}
%
% Optional parameters for ROI fitting:
%   'fit_ROI_parms': cell array of parameters to be nonlinearly optimized
%     must be one or more of the default values
%     {default = {'ADC_long','ADC_hindered','ADC_free'}}
%   'ADC_long_range':
%     {default = [0.5e-3,1.5e-3]}
%   'ADC_hindered_range'
%     {default = [0.5e-3,2e-3]}
%   'ADC_free_range'
%     {default = [2.5e-3,4e-3]}
%
% Optional parameters for multi-scale FOD plus isotropic fitting:
%   'lambda': regularization constant
%     {default = 0.1}
%   'ADC_long': longitudinal ADC
%     {default = 1e-3}
%   'ADC_trans_min': minimum transverse ADC
%     {default = 0}
%   'ADC_trans_max': maximum transverse ADC
%     {default = 0}
%   'num_ADC_trans': number of transverse ADC size scales
%     {default = 1}
%   'SH_order': spherical harmonic order -- must be even
%     {default = 4}
%   'iso_restricted_flag': [0|1] model isotropic diffusion of restricted water
%     {default = 0}
%   'iso_hindered_flag': [0|1] model isotropic diffusion of hindered water
%     {default = 1}
%   'iso_free_flag': [0|1] model isotropic diffusion of free water
%     {default = 1}
%   'ADC_hindered': ADC of isotropic hindered water (e.g. edema)
%     {default = 1.2e-3}
%   'ADC_free': apparent diffusion coefficient (ADC) of
%               isotropic free water (e.g. CSF)
%     {default = 3.3e-3}
%   'ADC_iso_min': minimum isotropic ADC
%     {default = 0}
%   'ADC_iso_max': maximum isotropic ADC
%     {default = 3e-3}
%   'num_ADC_iso': number of isotropic ADC size scales
%     {default = 0}
%   'ADC_iso_vals': vector of isotropic ADC values
%     if empty, will be set according to
%       ADC_iso_min, ADC_iso_max, and num_ADC_iso
%     {default = []}
%   'norm_flag': normalize data to b=0 image
%     {default = 0}
%
% Output:
%   results: struct containing parameter estimates
%
% Created:  03/03/14 by Don Hagler
% Last Mod: 07/31/15 by Don Hagler
%

%% todo: run in batches if nvox is large

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];

% check input arguments
if ~mmil_check_nargs(nargin, 2), return; end;
parms = check_input(vals,qmat,varargin);

% perform tensor fit for each voxel
if parms.trim_flag
  if parms.verbose
    fprintf('%s: trimming ROI outliers based on tensor fit for %d voxels...\n',...
      mfilename,parms.nvox);
  end;
  parms = trim_ROI(parms);
end;

% average across voxels
if parms.mean_flag
  if parms.verbose
    fprintf('%s: average across voxels %d voxels...\n',...
      mfilename,parms.nvox);
  end;
  parms = average_ROI(parms);
end;

% initialize results
results = init_results(parms);

% nonlinear fit for ROI parameters
if parms.fit_ROI_flag
  if parms.verbose
    fprintf('%s: fitting ROI parameters...\n',mfilename);
  end;
  results = fit_ROI(parms,results);
end;

% fit MFI
if parms.verbose
  fprintf('%s: fitting data for %d voxels with multi-scale FOD plus iso model...\n',...
    mfilename,parms.nvox);
end;
results = fit_MFI(parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(vals,qmat,options)
  parms = mmil_args2parms(options,{...
    'vals',vals,[],...
    'qmat',qmat,[],...
  ...
    'bvals',1000,[],...
    'nob0_flag',false,[false true],...
    'nonlin_flag',true,[false true],...
    'mean_flag',false,[false true],...
    'fit_ROI_flag',true,[false true],...
    'verbose',true,[false true],...
  ... % ROI trimming
    'trim_flag',false,[false true],...
    'trim_fact',4.7,[0,Inf],...
    'trim_thresh',0,[0,1],...
  ... % ROI fitting
    'fit_ROI_parms',{'ADC_long','ADC_hindered','ADC_free',},[],...
    'ADC_long_range',[0.5e-3,1.5e-3],[],...
    'ADC_hindered_range',[0.5e-3,2e-3],[],...
    'ADC_free_range',[2.5e-3,4e-3],[],...
  ... % MFI fitting
    'lambda',0.1,[],...
    'ADC_long',1e-3,[],...
    'ADC_trans_min',0,[],...
    'ADC_trans_max',0,[],...
    'num_ADC_trans',1,[],...
    'iso_free_flag',true,[false true],...
    'iso_hindered_flag',true,[false true],...
    'iso_restricted_flag',false,[false true],...
    'ADC_hindered',1.2e-3,[],...
    'ADC_free',3.3e-3,[],...
    'ADC_iso_min',0,[],...
    'ADC_iso_max',3e-3,[],...
    'num_ADC_iso',0,[],...
    'ADC_iso_vals',[],[],...
    'SH_order',4,[],...
    'norm_flag',false,[false true],...
    'nlfit_method','fmincon',{'lsqnonlin','fmincon'},...
    'nlfit_display','none',[],...
    'nlfit_F0_bounds',[0,100],[-Inf,Inf],...
    'nlfit_FX_bounds',[-15,15],[-Inf,Inf],...
    'nlfit_FX_flag',false,[false true],...
  ... % parameter names to be passed on
    'DT_tags',{'bvals','nob0_flag','nonlin_flag'},[],...
    'MFI_tags',{'volmask','bvals','nob0_flag',...
                'lambda','ADC_long','ADC_trans_min','ADC_trans_max',...
                'num_ADC_trans','ADC_iso_min','ADC_iso_max','num_ADC_iso',...
                'ADC_iso_vals','iso_free_flag','iso_hindered_flag',...
                'iso_restricted_flag','ADC_free','ADC_hindered',...
                'SH_order','norm_flag','nonlin_flag',...
                'nlfit_method','nlfit_display',...
                'nlfit_F0_bounds','nlfit_FX_bounds',...
                'nlfit_FX_flag'},[],...
  });
  if length(size(parms.vals))>2
    error('input vals must be 2D matrix');
  end;
  parms.nvox = size(parms.vals,1);
  parms.ndirs = size(parms.vals,2);

  % initialize fit ROI parameters
  parms.num_ROI_parms = length(parms.fit_ROI_parms);
  if parms.fit_ROI_flag
    fit_ROI_parms = [];
    if ~parms.iso_hindered_flag
      parms.fit_ROI_parms = setdiff(parms.fit_ROI_parms,{'ADC_hindered'});
    end;
    if ~parms.iso_free_flag
      parms.fit_ROI_parms = setdiff(parms.fit_ROI_parms,{'ADC_free'});
    end;
    if ~parms.num_ADC_trans
      parms.fit_ROI_parms = setdiff(parms.fit_ROI_parms,{'ADC_long'});
    end;
    parms.num_ROI_parms = length(parms.fit_ROI_parms);
    for i=1:parms.num_ROI_parms
      beta_name = parms.fit_ROI_parms{i};
      beta_range = parms.([beta_name '_range']);
      if beta_range(1)>beta_range(2)
        error('upper bound is less than lower bound for %s',beta_name);
      end;
      if beta_range(1)~=beta_range(2)
        fit_ROI_parms{end+1} = beta_name;
      end;
    end;
    parms.fit_ROI_parms = fit_ROI_parms;
    parms.num_ROI_parms = length(parms.fit_ROI_parms);
    parms.fit_ROI_flag = (parms.num_ROI_parms~=0);
  end;
  if parms.num_ROI_parms
    switch parms.nlfit_method
      case 'lsqnonlin'
        opts = {'Display',parms.nlfit_display,...
                'MaxIter',100,...
                'TolFun',1e-10};
        try % for R2009b
          parms.nlfit_options = optimset(opts{:},...
            'Algorithm','levenberg-marquardt');
        catch % for R2007a
          parms.nlfit_options = optimset(opts{:},...
            'LargeScale','off');
        end;
      case 'fmincon'
        opts = {'Display',parms.nlfit_display,...
                'MaxFunEvals',Inf,...
                'MaxIter',100,...
                'TolFun',1e-6,...
                'TolX',1e-7};
        try % for R2009b
          parms.nlfit_options = optimset(opts{:},...
            'Algorithm','active-set');
        catch % for R2007a
          parms.nlfit_options = optimset(opts{:},...
            'LargeScale','off');
        end;
    end;
  end;
return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)  
  if isfield(parms,'ind_excl')
    results.vals_orig = parms.vals_orig;
    results.nvox_orig = parms.nvox_orig;
    results.ind_excl = parms.ind_excl;
    results.ind_keep = parms.ind_keep;
  end;
  results.vals = parms.vals;
  results.nvox = parms.nvox;
  results.ndirs = parms.ndirs;
  results.qmat = parms.qmat;
  results.bvals = parms.bvals;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = trim_ROI(parms)
  % tensor fit
  args = mmil_parms2args(parms,parms.DT_tags);
  results = rsi_fit_ROI_DT(parms.vals,parms.qmat,args{:});
  % calculate weighting factors based on difference from median
  w_MD = calc_weights(results.MD.vals,parms.trim_fact);
  w_FA = calc_weights(results.FA.vals,parms.trim_fact);
  w_T2 = calc_weights(results.b0.vals,parms.trim_fact);
  w = power(w_MD .* w_FA .* w_T2,1/3);
  % exclude outliers
  parms.ind_excl = find(w<=parms.trim_thresh);
  parms.ind_keep = find(w>parms.trim_thresh);
  parms.vals_orig = parms.vals;
  parms.nvox_orig = parms.nvox;
  parms.vals = parms.vals_orig(parms.ind_keep,:);
  parms.nvox = size(parms.vals,1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = average_ROI(parms)
  parms.vals = mean(parms.vals,1);
  parms.nvox = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = calc_weights(x,df)
  mu = median(x);
  nu = df*median(abs(x - mu))/0.6745; % scaled median absolute deviation
  w = max(0,1-((x-mu)/nu).^2); % Tukey's bisquare function
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = fit_MFI(parms,results)
  if ~exist('results','var'), results = []; end;
  args = mmil_parms2args(parms,parms.MFI_tags);
  tmp_results = rsi_fit_ROI_MFI(parms.vals,parms.qmat,args{:});
  tmp_fieldnames = fieldnames(tmp_results);
  for i=1:length(tmp_fieldnames)
    results.(tmp_fieldnames{i}) = tmp_results.(tmp_fieldnames{i});
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,results] = fit_ROI(parms,results)
  % initialize parameter estimates and lower and upper bounds
  [beta_init,lb,ub] = init_ROI_parms(parms);
  % nonlinear fit
  switch parms.nlfit_method
    case 'lsqnonlin'
      beta = lsqnonlin(@(beta) costfunc(beta,parms),...
                      beta_init,lb,ub,parms.nlfit_options);
    case 'fmincon'
      beta = fmincon(@(beta) sum(costfunc(beta,parms).^2),...
                      beta_init,[],[],[],[],lb,ub,[],parms.nlfit_options);
  end;
  % set parms from beta
  parms = set_ROI_parms(beta,parms);
  for i=1:parms.num_ROI_parms
    beta_name = parms.fit_ROI_parms{i};
    results.(beta_name) = parms.(beta_name);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta,lb,ub] = init_ROI_parms(parms)
  beta = []; lb = []; ub = [];
  for i=1:parms.num_ROI_parms
    beta_name = parms.fit_ROI_parms{i};
    beta_val = parms.(beta_name);
    beta_range = parms.([beta_name '_range']);
    beta_val = min(max(beta_val,beta_range(1)),beta_range(2));
    beta(i) = beta_val;
    lb(i) = beta_range(1);
    ub(i) = beta_range(2);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = costfunc(beta,parms)
  cost = [];
  % set parms from beta
  parms = set_ROI_parms(beta,parms);
  % fit data with multi-scale FOD plus iso model
  results = fit_MFI(parms);
  % residual error vector for entire ROI
  cost = mmil_rowvec(results.err);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_ROI_parms(beta,parms)
  for i=1:parms.num_ROI_parms
    parms.(parms.fit_ROI_parms{i}) = beta(i);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

