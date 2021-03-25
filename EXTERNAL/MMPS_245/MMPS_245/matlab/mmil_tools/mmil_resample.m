function [vol_res,M_res] = mmil_resample(vol,M,varargin)
%function [vol_res,M_res] = mmil_resample(vol,M,[options])
%
% Purpose: resamples a multi-frame volume with specified
%   registration matrix and output vox2ras matrix
%
% Required Input:
%   vol: 4D volume matrix
%   M: 4x4 vox2ras matrix
%
% Optional Parameters:
%   'nvoxels': vector of [nx,ny,nz] (in numbers of voxels)
%     if empty, keep input nvoxels
%     {default = []}
%   'resolution': vector of [resx,resy,resz] (in mm)
%     if empty, keep input resolution
%     {default = []}
%   'orient': output slice orientation
%      e.g. 'RAS', 'LPI', 'PRI', etc.
%      if not specified, keep input orientation
%     {default = []}
%   'M_ref': 4x4 vox2ras matrix for reference volume
%     if not specified, use M constructed from orient or M from fname
%     {default = []}
%   'nvox_ref': number of voxels [nx,ny,nz] for reference volume
%     ignored if M_ref is empty
%     if not specified, will use number of voxels for input volume
%     used to calculate offset for output M
%     {default = []}
%   'M_reg': registration matrix in LPH (scanner) coordinates applied to volume
%     {default = eye(4)}
%   'deoblique_flag': [0|1] resample oblique slices to on-axis
%     {default = 1}
%   'smooth': isotropic smoothing kernel sigma (in voxels)
%     {default = 0}
%   'interpm': [0|1|2|3|4|5] interpolation method
%      0:nearest  1:linear  2:cubic  3:Key's spline  4:cubic spline  5: hamming sinc
%     {default = 2}
%   'bclamp' : [0|1] set negative values to zero
%     {default = 1}
%
% Output:
%   vol_res: resampled volume
%   M_res: vox2ras matrix for resampled volume
%
% Created:   05/21/12 by Don Hagler
% Last Mod:  07/29/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
vol_res = []; M_res = [];
parms = mmil_args2parms(varargin, { ...
  'nvoxels',[],[],...
  'resolution',[],[],...
  'orient',[],[],...
  'M_ref',[],[],...
  'nvox_ref',[],[],...
  'M_reg',eye(4),[],...
  'deoblique_flag',true,[false,true],...
  'smooth',0,[0,100],...
  'interpm',2,[0:5],...
  'bclamp',true,[false,true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,ny,nz,nf] = size(vol);
nvoxels = [nx,ny,nz];
if isempty(parms.nvoxels)
  parms.nvoxels = nvoxels;
elseif length(parms.nvoxels)>3
  parms.nvoxels = parms.nvoxels(1:3);
end;
if isempty(parms.resolution)
  parms.resolution = fs_voxsize(M);
end;

if isempty(parms.nvox_ref) || isempty(parms.M_ref)
  parms.nvox_ref = parms.nvoxels;
end;
if length(parms.nvox_ref)>3
  parms.nvox_ref = parms.nvox_ref(1:3);
elseif length(parms.nvox_ref)~=3
  error('nvox_ref has %d elements (must have three)',length(parms.nvox_ref));
end;

if isempty(parms.M_ref)
  parms.M_ref = M;
end;

if ~isempty(parms.orient)
  M_res = mmil_construct_M('orient',parms.orient,'scale',parms.resolution);
else
  [orient,oblique_flag] = fs_read_orient([],parms.M_ref);
  if oblique_flag && ~parms.deoblique_flag
    M_res = parms.M_ref;
    voxel_sizes = fs_voxsize(M_res);
    if any(parms.resolution~=voxel_sizes)
      % adjust voxel dimensions
      for i=1:3
        M_res(1:3,i) = M_res(1:3,i) * parms.resolution(i)/voxel_sizes(i);
      end;
    end;
  else
    M_res = mmil_construct_M('orient',orient,'scale',parms.resolution);
  end;
end;
% calculate offset
M_res(1:3,4) = 0; % set to zero so that next line gives edge coordinates
M_res(1:3,4) = parms.M_ref(1:3,:)*[parms.nvox_ref/2+1 1]' - ...
               M_res(1:3,:)*[parms.nvoxels/2+1 1]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample each frame
vol_res = zeros([parms.nvoxels,nf],'single');
for f=1:nf
  tmp_vol = single(squeeze(vol(:,:,:,f)));
  tmp_vol = mmil_resample_vol(tmp_vol,M,...
    'M_ref',M_res,'nvox_ref',parms.nvoxels,...
    'interpm',parms.interpm,'bclamp',parms.bclamp,...
    'M_reg',parms.M_reg);
  if parms.smooth>1
    tmp_vol = mmil_smooth3d(tmp_vol,parms.smooth,parms.smooth,parms.smooth);
  end;
  vol_res(:,:,:,f) = single(tmp_vol);
end;

