function mmil_normvol_asegroi(fname_in,fname_out,fname_aseg,roilist,forceflag)
%function mmil_normvol_asegroi(fname_in,fname_out,fname_aseg,roilist,forceflag)
%
% Purpose: normalize a volume by the average value inside voxels specified
%  by aseg ROI numbers in a vector 'roilist'
%
% Required Input:
%   fname_in: file name of volume to be normalized
%   fname_out: output file name
%   fname_aseg: file name of freesurfer aseg volume
%     (full paths all)
%   roilist: vector of aseg ROI numbers
%
% Optional Parameters:
%   forceflag: whether to overwrite existing fname_out
%
% Created:  03/09/09 by Matt Erhart
% Last Mod: 10/27/10 by Don Hagler
%

if ~mmil_check_nargs(nargin,4), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

if exist(fname_out) && ~forceflag, return; end;

if ~exist(fname_in), error('fname_in %s does not exist',fname_in); end;
if ~exist(fname_aseg), error('fname_aseg %s does not exist',fname_aseg); end;

% load fname_in and fname_aseg
[vol,M,tmp,volsz] = fs_load_mgh(fname_in);
[vol_aseg,M_aseg,tmp,volsz_aseg] = fs_load_mgh(fname_aseg);

% check that vol sizes and coordinates match
if (any(volsz ~= volsz_aseg)) || any(M(:)~=M_aseg(:))
  error('mismatch between volume in fname_in and fname_aseg');
end;

% find voxels in vol_aseg with values matching the values specified in roilist? (a vector of integers)
k = find((ismember(vol_aseg, roilist)));

% normalize vol by average of roi voxels
vol = vol/mean(vol(k));

% save normalized volume to fname_out
fs_save_mgh(vol,fname_out,M);


