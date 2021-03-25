function results = rsi_fit_ROI_MFI(vals,qmat,varargin)
%function results = rsi_fit_ROI_MFI(vals,qmat,[options])
%
% Purpose: perform multi-scale FOD plus isotropic fit
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
%   'lambda': regularization constant
%     {default = 0.1}
%   'ADC_long': longitudinal ADC
%     {default = 1e-3}
%   'ADC_trans_min': minimum transverse ADC
%     {default = 0}
%   'ADC_trans_max': maximum transverse ADC
%     {default = 0.9e-3}
%   'num_ADC_trans': number of transverse ADC size scales
%     {default = 5}
%   'SH_order': spherical harmonic order -- must be even
%     {default = 4}
%   'iso_restricted_flag': [0|1] model isotropic diffusion of restricted water
%     {default = 1}
%   'iso_hindered_flag': [0|1] model isotropic diffusion of hindered water
%     {default = 1}
%   'iso_free_flag': [0|1] model isotropic diffusion of free water
%     {default = 1}
%   'ADC_hindered': ADC of isotropic hindered water (e.g. edema)
%     {default = 1.5e-3}
%   'ADC_free': apparent diffusion coefficient (ADC) of
%               isotropic free water (e.g. CSF)
%     {default = 3e-3}
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
%   'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit
%     {default = 1}
%
% Output:
%   results: struct containing parameter estimates and other results
%
% Created:  02/25/14 by Don Hagler
% Last Mod: 07/31/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];

% check input arguments
if ~mmil_check_nargs(nargin, 2), return; end;
parms = check_input(vals,qmat,varargin);

% initialize results
results = init_results(parms);

% perform multi-scale FOD plus isotropic fit for each voxel
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
    'lambda',0.1,[],...
    'ADC_free',3e-3,[],...
    'ADC_hindered',1.5e-3,[],...
    'ADC_long',1e-3,[],...
    'ADC_trans_min',0,[],...
    'ADC_trans_max',0.9e-3,[],...
    'num_ADC_trans',5,[],...
    'iso_free_flag',true,[false true],...
    'iso_hindered_flag',true,[false true],...
    'iso_restricted_flag',true,[false true],...
    'ADC_iso_min',0,[],...
    'ADC_iso_max',3e-3,[],...
    'num_ADC_iso',0,[],...
    'ADC_iso_vals',[],[],...
    'SH_order',4,[],...
    'norm_flag',false,[false true],...
    'nonlin_flag',true,[false true],...
    'nlfit_method','fmincon',{'lsqnonlin','fmincon'},...
    'nlfit_display','none',[],...
    'nlfit_F0_bounds',[0,100],[-Inf,Inf],...
    'nlfit_FX_bounds',[-15,15],[-Inf,Inf],...
    'nlfit_FX_flag',false,[false true],...
  ...
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

function results = fit_MFI(parms,results)
  % calculate multi-shell FOD fit
  args = mmil_parms2args(parms,parms.MFI_tags);
  MFfit = rsi_fit_MFOD(parms.vol,parms.qmat,args{:});
  % calculate MF fit measures
  MFmeas = rsi_calc_MFmeas(MFfit,'mask_flag',0);
  % calculate stats for MF measures
  results.F0 = calc_stats(MFmeas.volF0);
  results.F2 = calc_stats(MFmeas.volF2);
  results.F4 = calc_stats(MFmeas.volF4);
  results.I  = calc_stats(MFmeas.volI);
  results.Ir = calc_stats(MFmeas.volIr);
  results.Ih = calc_stats(MFmeas.volIh);
  results.If = calc_stats(MFmeas.volIf);
  results.AU = calc_stats(MFmeas.volAU);
  % calculate fit and residual error
  MF = reshape(MFfit.volMF,[parms.nvox,MFfit.nb]);
  results.fit = (MFfit.A*MF')';
  results.err = results.vals - results.fit;
  results.norm_err = mean(sqrt(sum(results.err.^2,1)) ./ sqrt(sum(results.vals.^2,1)));
  % save relevant parameters
  results.ADC_trans_vals = MFfit.ADC_trans_vals;
  results.ADC_iso_vals = MFfit.ADC_iso_vals;
  % save fit and measures
  results.MFfit = MFfit;
  results.MFmeas = MFmeas;
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

