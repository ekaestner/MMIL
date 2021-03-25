function results = rsi_fit_ROI_MI(vals,qmat,varargin)
%function results = rsi_fit_ROI_MI(vals,qmat,[options])
%
% Purpose: perform multi-scale isotropic fit
%   on diffusion MRI data for ROI
%
% Required Parameters:
%   vals: 2D volume with size = [nvox,nf]
%   qmat: matrix of diffusion direction vectors with size = [nf,3]
%
% Optional Parameters:
%   'bvals': vector of b values
%     one for all, or one for each diffusion direction
%     {default = 1000}
%   'nob0_flag': [0|1] whether to exclude b=0 images from fits
%     {default = 0}
%   'ADC_iso_min': minimum isotropic ADC
%     {default = 0}
%   'ADC_iso_max': maximum isotropic ADC
%     {default = 3e-3}
%   'num_ADC_iso': number of isotropic ADC size scales
%     {default = 10}
%   'ADC_iso_vals': vector of isotropic ADC values
%     if empty, will be set according to
%       ADC_iso_min, ADC_iso_max, and num_ADC_iso
%     {default = []}
%   'norm_flag': normalize data to b=0 image
%     {default = 0}
%   'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit
%     {default = 1}
%
% Output:
%   results: struct containing parameter estimates and other results
%
% Created:  02/25/14 by Don Hagler
% Last Mod: 03/04/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];

% check input arguments
if ~mmil_check_nargs(nargin, 2), return; end;
parms = check_input(vals,qmat,varargin);

% initialize results
results = init_results(parms);

% perform multi-scale isotropic diffusion fit for each voxel
results = fit_MI(parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(vals,qmat,options)
  parms = mmil_args2parms(options,{...
    'vals',vals,[],...
    'qmat',qmat,[],...
  ...
    'bvals',1000,[],...
    'nob0_flag',false,[false true],...
    'lambda',0.1,[],...
    'ADC_iso_min',0,[],...
    'ADC_iso_max',3e-3,[],...
    'num_ADC_iso',0,[],...
    'ADC_iso_vals',[],[],...
    'norm_flag',false,[false true],...
    'nonlin_flag',true,[false true],...
    'nlfit_method','fmincon',{'lsqnonlin','fmincon'},...
    'nlfit_display','none',[],...
    'nlfit_bounds',[0,100],[-Inf,Inf],...
  ...
    'MI_tags',{'volmask','bvals','nob0_flag',...
               'lambda','ADC_iso_min','ADC_iso_max','num_ADC_iso',...
               'ADC_iso_vals','norm_flag','nonlin_flag',...
               'nlfit_method','nlfit_display','nlfit_bounds'},[],...
  });
  if length(size(parms.vals))>2
    error('input vals must be 2D matrix');
  end;
  parms.nvox = size(parms.vals,1);
  parms.nf = size(parms.vals,2);
  parms.vol = reshape(parms.vals,[parms.nvox,1,1,parms.nf]);
  parms.volmask = ones(size(parms.vals,1),1,1);
return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results.vals = parms.vals;
  results.qmat = parms.qmat;
  results.bvals = parms.bvals;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = fit_MI(parms,results)
  % calculate multi-shell isotropic fit
  args = mmil_parms2args(parms,parms.MI_tags);
  MIfit = rsi_fit_iso(parms.vol,parms.qmat,args{:});
  % calculate stats for iso measures
  results.I = calc_stats(MIfit.volMI);
  % calculate fit and residual error
  MI = reshape(MIfit.volMI,[parms.nvox,MIfit.nb]);
  results.fit = (MIfit.A*MI')';
  results.err = results.vals - results.fit;
  results.norm_err = mean(sqrt(sum(results.err.^2,1)) ./ sqrt(sum(results.vals.^2,1)));
  % save relevant parameters
  results.ADC_iso_vals = MIfit.ADC_iso_vals;
  % save fit and measures
  results.MIfit = MIfit;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = calc_stats(vals)
  vals = reshape(vals,[size(vals,1),size(vals,4)]);
  output = [];
  output.vals = vals;
  output.mean = mean(vals);
  output.std = std(vals);
  output.med = median(vals);
return;

