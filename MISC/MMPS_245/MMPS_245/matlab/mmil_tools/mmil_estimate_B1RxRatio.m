function vol_ratio = mmil_estimate_B1RxRatio(volA,volB,varargin)
%function vol_ratio = mmil_estimate_B1RxRatio(volA,volB,[options])
%
% Purpose: estimate smoothly varying intensity variations between two volumes
%   assumes images are well-registered
%
% Required Input:
%   volA: input volume A
%   volB: input volume B
%
% Optional Input:
%   'thresh': threshold for noise, multiplied by estimated noise level
%     {default = 4}
%   'smooth_sigma': sigma of smoothing kernel before calculating ratio
%     {default = 32}
%   'verbose': [0|1] display status messages
%     {default = 0}
%
% Output:
%  vol_ratio: volume of ratio between volA and volB (after smoothing)
%
% Created:  01/08/15 by Don Hagler
% Last Mod: 01/20/15 by Don Hagler
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_ratio = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'thresh',4,[],...
  'smooth_sigma',32,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% smooth images with padding at edges to avoid wrap-around
volA = mmil_smooth3d_pad(volA,parms.smooth_sigma);
volB = mmil_smooth3d_pad(volB,parms.smooth_sigma);

% identify non-noise voxels
noiseA = mmil_calc_noise_and_signal_levels(volA);
noiseB = mmil_calc_noise_and_signal_levels(volB);
ind_mask = find(volA>parms.thresh*noiseA & volB>parms.thresh*noiseB);
ind_noise = find(volA<=parms.thresh*noiseA | volB<=parms.thresh*noiseB);

% calculate ratio of means
mu_A = mean(volA(ind_mask));
mu_B = mean(volB(ind_mask));
mu_ratio = mu_A/mu_B;

% calculate ratio for each voxel
vol_ratio = sqrt(abs((volA.^2-noiseA^2)./(volB.^2-noiseB^2)));
vol_ratio = max(0.5*mu_ratio,min(2*mu_ratio,vol_ratio));
vol_ratio(ind_noise) = 1;

