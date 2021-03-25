function fs_maskvol(fname_in,fname_mask,fname_out,varargin)
%function fs_maskvol(fname_in,fname_mask,fname_out,varargin)
%
% Purpose: apply a mask volume to another volume
%
% Usage: fs_maskvol(fname_in,fname_mask,fname_out,'key1', value1,...);
%
% Required Input:
%  fname_in: full path of input volume file (mgh format)
%  fname_mask: full path of mask volume file (mgh format)
%  fname_out: full path of output volume file (mgh format)
%
% Optional Input:
%  'thresh': threshold applied to mask
%    {default: 0.5}
%  'forceflag': [0|1] whether to overwrite fname_mask if it already exists
%    {default: 0}
%
% NOTE: input volume and mask must have same dimensions and be
%   registered voxel by voxel (ignores vox2ras matrices)
%
% Created:  09/08/09 by Don Hagler
% Last Mod: 09/08/09 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,3)) return; end;
parms = mmil_args2parms(varargin, { ...
  'thresh',0.5,[0,Inf],...
  'forceflag',false,[false true],...
});

if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;
if ~exist(fname_mask,'file')
  error('file %s not found',fname_mask);
end;

if ~exist(fname_out,'file') || parms.forceflag
  [vol_in,M_in,mp,volsz_in] = fs_load_mgh(fname_in);
  [vol_mask,M_mask,mp,volsz_mask] = fs_load_mgh(fname_mask);
  if any(volsz_in(1:3)~=volsz_mask(1:3))
    error('dimensions of vol and mask do not match')
  end;
  vol_mask(vol_mask<parms.thresh)=0;
  vol_mask(vol_mask>=parms.thresh)=1;
  vol_out = zeros(volsz_in);
  for f=1:size(vol_in,4)
    vol_out(:,:,:,f) = vol_in(:,:,:,f).*vol_mask;
  end;
  fs_save_mgh(vol_out,fname_out,M_in);
end;

