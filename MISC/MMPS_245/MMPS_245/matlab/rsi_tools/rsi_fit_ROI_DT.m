function results = rsi_fit_ROI_DT(vals,qmat,varargin)
%function results = rsi_fit_ROI_DT(vals,qmat,[options])
%
% Purpose: perform tensor fit
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
%   'nonlin_flag': [0|1] use nonlinear fitting for multi-voxel fits
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

% perform tensor fit for each voxel
results = fit_DT(parms,results);

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
  ...
    'DT_tags',{'volmask','bvals','nob0_flag','nonlin_flag'},[],...
  });
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

function results = fit_DT(parms,results)
  % calculate tensor fit
  args = mmil_parms2args(parms,parms.DT_tags);
  DTfit = dti_fit_tensor(parms.vol,parms.qmat,args{:});
  % calculate tensor fit measures
  DTmeas = dti_calc_DTmeas(DTfit,0);
  % calculate stats for DT measures
  results.FA = calc_stats(DTmeas.volFA);
  results.LD = calc_stats(DTmeas.volE(:,:,:,1));
  results.TD = calc_stats(mean(DTmeas.volE(:,:,:,2:3),4));
  results.MD = calc_stats(mean(DTmeas.volE,4));
  results.b0 = calc_stats(DTmeas.volb0);
  % calculate fit and residual error
  T = reshape(DTfit.volDT,[parms.nvox,size(DTfit.Q,2)]);
  results.fit = exp(DTfit.Q*T')';
  results.err = results.vals - results.fit;
  results.norm_err = mean(sqrt(sum(results.err.^2,1)) ./ sqrt(sum(results.vals.^2,1)));
  % save fit and measures
  results.DTfit = DTfit;
  results.DTmeas = DTmeas;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

