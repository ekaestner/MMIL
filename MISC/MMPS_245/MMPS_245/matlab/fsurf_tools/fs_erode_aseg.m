function fs_erode_aseg(fname_in,fname_out,varargin)
%function fs_erode_aseg(fname_in,fname_out,[options])
%
% Purpose: Creating eroded segmentation ROIs to avoid boundaries
%    and prevent partial voluming in ROI analysis
% Description: Each ROI in the input segmentation file is smoothed
%    and then thresholded, creating an ROI with an outer boundary some
%    distance inside the "true" boundary
%
% Required Input:
%   fname_in: full path to input segmentation volume (e.g. aseg.mgz)
%   fname_out: full path to output segmentation volume (e.g. aseg_eroded.mgz)
%
% Optional Parameters:
%   'nvoxels': number of voxels to erode (iteratively smooth and threshold)
%     {default = 1}
%   'sigma': smoothing kernel sigma (voxels)
%     {default = 1}
%   'thresh': threshold value for binarizing smoothed masks
%     {default = 0.99}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': [0|1] toggle overwrite of existing output file
%     {default = 0}
%
% Created:  05/14/07 by Don Hagler
% Last Mod: 03/01/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'nvoxels',1,[1:100],...
  'sigma',1,[1e-2,10],...
  'thresh',0.99,[0.1,0.999],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
});

results = [];

if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;
if exist(fname_out,'file') & ~parms.forceflag, return; end;

% load input segmentation volume
if parms.verbose
  fprintf('%s: loading input segmentation volume %s...\n',mfilename,fname_in);
end;
[vol_aseg,M,mr_parms,volsz] = fs_load_mgh(fname_in);

% find unique roi numbers in vol_aseg
aseg_roicodes = unique(vol_aseg(:));
aseg_roicodes = aseg_roicodes(find(aseg_roicodes>0));
nrois = length(aseg_roicodes);

if parms.verbose
  fprintf('%s: eroding %d ROIs...\n',mfilename,nrois);
end;
vol_aseg2 = zeros(volsz);
for i=1:nrois
  tmp_vol = zeros(volsz);
  tmp_vol(vol_aseg==aseg_roicodes(i)) = 1.0;
  for j=1:parms.nvoxels
    tmp_vol = mmil_smooth3d(tmp_vol,parms.sigma,parms.sigma,parms.sigma);
    tmp_vol = 1.0*(tmp_vol>=parms.thresh);
  end;  
  vol_aseg2(tmp_vol>0) = aseg_roicodes(i);
end;

if parms.verbose
  fprintf('%s: saving output segmentation volume %s...\n',mfilename,fname_out);
end;
fs_save_mgh(vol_aseg2,fname_out,M,mr_parms);

