function fs_erode_mask(fname_in,fname_out,varargin)
%function fs_erode_mask(fname_in,fname_out,[options])
%
% Purpose: Creating eroded ROIs to avoid boundaries
%    and prevent partial voluming in ROI analysis
% Description: input file is binarized, smoothed
%    and then thresholded, creating an ROI with an outer boundary some
%    distance inside the "true" boundary
%
% Required Input:
%   fname_in: full path to input mask volume
%   fname_out: full path to output eroded mask volume
%
% Optional Parameters:
%   'nvoxels': number of voxels to erode (iteratively smooth and threshold)
%     {default = 1}
%   'sigma': smoothing kernel sigma (voxels)
%     {default = 1}
%   'thresh_init': threshold value for binarizing input image
%     {default = 1e-5}
%   'thresh': threshold value for binarizing smoothed masks
%     {default = 0.99}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': [0|1] toggle overwrite of existing output file
%     {default = 0}
%
% Created:  04/01/09 by Don Hagler
% Prev Mod: 08/15/12 by Don Hagler
% Last Mod: 07/13/17 by Don Hagler
%

if ~mmil_check_nargs(nargin, 2), return; end;
parms = mmil_args2parms(varargin,{...
  'nvoxels',1,[1:100],...
  'sigma',1,[1e-2,10],...
  'thresh_init',1e-5,[0,Inf],...
  'thresh',0.99,[0.1,0.999],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
});

if ~exist(fname_in,'file')
  error('input file %s not found',fname_in);
end;
if exist(fname_out,'file') & ~parms.forceflag, return; end;

% load input segmentation volume
if parms.verbose
  fprintf('%s: loading input mask volume %s...\n',mfilename,fname_in);
end;
[vol,M,mr_parms,volsz] = fs_load_mgh(fname_in);

% binarize
vol = 1.0 * (vol>1-parms.thresh_init);
if parms.verbose
  fprintf('%s: eroding mask...\n',mfilename);
end;
for j=1:parms.nvoxels
  vol = mmil_smooth3d(vol,parms.sigma,parms.sigma,parms.sigma);
  vol = 1.0*(vol>=parms.thresh);
end;  

if parms.verbose
  fprintf('%s: saving output mask volume %s...\n',mfilename,fname_out);
end;
fs_save_mgh(vol,fname_out,M,mr_parms);

