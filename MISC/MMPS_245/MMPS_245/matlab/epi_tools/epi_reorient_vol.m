function epi_reorient_vol(fname_in,fname_out,out_orient,forceflag)
%function epi_reorient_vol(fname_in,fname_out,out_orient,forceflag)
%
% Purpose: reorient volume to desired orientation
%
% Required:
%   fname_in:  full path name of input mgh/mgz file
%   fname_out: full path name of output mgh/mgz file
%
% Optional:
%   out_orient: desired orientation of output file
%      e.g. 'RAS', 'LPI', 'PRI', etc.
%     {default: 'LPS'}
%   forceflag: [0|1] whether to overwrite existing output
%     {default: 0}
%
% NOTE: when flipping dimensions, center coordinates in vox2ras matrix
%  will be adjusted for 1 voxel shift, but output volume will not be changed
%  other than permutation and flipping dimensions
%
% Created:  02/14/10 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('out_orient','var') || isempty(out_orient), out_orient='LPS'; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag=0; end;

if exist(fname_out,'file') && ~forceflag
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load input volume
[vol,M,mrp,volsz] = fs_load_mgh(fname_in);

% reorient volume
[vol_out,M_out] = fs_reorient(vol,M,out_orient);

% save output volume
fs_save_mgh(vol_out,fname_out,M_out,mrp);

