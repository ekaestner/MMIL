function [vol,M]= fs_reorient(vol,M,out_orient)
%function [vol_out,M_out]= fs_reorient(vol,M,out_orient)
%
% Purpose: reorient volume to desired orientation
%
% Required:
%   vol: 4D volume matrix
%   M  : 4x4 vox2ras transform
%
% Optional:
%   out_orient: desired orientation of output file
%      e.g. 'RAS', 'LPI', 'PRI', etc.
%     {default: 'LPS'}
%
% Output:
%   vol_out: reoriented 4D output volume 
%     M_out: reoriented 4x4 vox2ras transform
%
% NOTE: when flipping dimensions, output volume will not be shifted by 1 voxel
%   and center coordinates in vox2ras matrix will not be adjusted
%
% Created:  02/24/12 by Vijay Venkatraman
% Last Mod: 05/21/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('out_orient','var') || isempty(out_orient), out_orient='LPS'; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for valid out_orient
out_orient = upper(out_orient);
if any(~ismember(out_orient(:),{'R','L','A','P','S','I'})) ||...
   length(unique(out_orient(:)))~=3 ||...
   isempty(regexp(out_orient,'[LR]')) ||...
   isempty(regexp(out_orient,'[AP]')) ||...
   isempty(regexp(out_orient,'[SI]'))
  error('invalid output orientation %s',out_orient);
end;

% get orientation for input volume
in_orient = fs_read_orient([],M);

if strcmp(in_orient,out_orient), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine how to permute and flip dims
[permvec,flipvec] = fs_compare_orient(in_orient,out_orient);

% permute volume
vol = permute(vol,[permvec,4]);

% convert to LPH scanner coordinates
M = M_RAS_TO_LPH*M;

% reorder vox2ras matrix
M = M(:,[permvec,4]);

% calculate center of volume
volsz= size(vol);
in_ctr_vox = (volsz(permvec))/2+1;
in_ctr = M(1:3,:)*[in_ctr_vox 1]';

% flip dimensions as necessary
out_ctr_vox = (volsz(permvec))/2+1;
for j=1:3
  if flipvec(j)<0
    vol = flipdim(vol,j);
    M(:,j) = -1*M(:,j);
    out_ctr_vox(j) = out_ctr_vox(j) - 1;
  end;
end;

% adjust for difference in center of coordinates
out_ctr = M(1:3,:)*[out_ctr_vox 1]';
M(1:3,4) = M(1:3,4) - (out_ctr-in_ctr);

% convert to RAS brain coordinates
M = M_LPH_TO_RAS*M;

