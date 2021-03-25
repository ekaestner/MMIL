function [sf,r] = mmil_normvol_linfit(fname_in,fname_out,fname_ref,varargin)
%function [sf,r] = mmil_normvol_linfit(fname_in,fname_out,fname_ref,[options])
%
% Purpose: normalize a volume by a constant scaling factor
%   derived from a linear fit with another input volume
%
% Required Input:
%   fname_in: file name of volume to be normalized
%   fname_out: output file name
%   fname_ref: file name of reference volume
%
% Optional Parameters ('key',value pairs):
%   'fname_mask': file name of brain mask volume
%     if not supplied, will create one using mmil_quick_brainmask
%     {default = []}
%   'fname_sf': output file name for scaling factor (sf) and correlation (r)
%     {default = []}
%   'scalefact': additional scaling factor applied to output volume
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   sf: scaling factor applied to input volume (combined with input scalefact)
%   r: correlation between input and ref values for brain voxels
%
% Created:  10/27/12 by Don Hagler
% Last Mod: 03/24/14 by Don Hagler
%

sf = nan; r = nan;
if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'fname_sf',[],[],...
  'scalefact',1,[],...
  'forceflag',false,[false true],...
...
  'mask_thresh',0.5,[0,1],...
  'vol_thresh',1e-6,[0,1e6],...
  'ref_thresh',1e-6,[0,1e6],...
});

skip_flag = 1;
if ~exist(fname_out,'file'), skip_flag = 0; end;
if ~isempty(parms.fname_sf) && ~exist(parms.fname_sf,'file'), skip_flag = 0; end;
if skip_flag && ~parms.forceflag, return; end;

if ~exist(fname_in), error('fname_in %s does not exist'); end;
if ~exist(fname_ref), error('fname_ref %s does not exist'); end;

% load fname_in and fname_aseg
[vol,M,tmp,volsz] = fs_load_mgh(fname_in);
[vol_ref,M_ref,tmp,volsz_ref] = fs_load_mgh(fname_ref);

%check that vol sizes and coordinates match
if (any(volsz ~= volsz_ref)) || any(M(:)~=M_ref(:))
  error('mismatch between volume in fname_in and fname_ref');
end;

% load or create brain mask
if ~isempty(parms.fname_mask)
  if ~exist(fname_mask), error('fname_mask %s does not exist'); end;
  [vol_mask,M_mask,tmp,volsz_mask] = fs_load_mgh(parms.fname_mask);
  if (any(volsz ~= volsz_mask)) || any(M(:)~=M_mask(:))
    error('mismatch between volume in fname_in and fname_mask');
  end;
else
  vol_ctx = ctx_mgh2ctx(vol,M);
  vol_mask = mmil_quick_brainmask(vol_ctx);
  vol_mask = vol_mask.imgs;
end;

% get non-zero voxels inside brain mask
ind = find(vol_mask >= parms.mask_thresh &...
           vol >=  parms.vol_thresh &...
           vol_ref >= parms.ref_thresh);
vec = vol(ind);
vec_ref = vol_ref(ind);

% linear fit to get sf
sf = parms.scalefact*(vec\vec_ref);

% correlation as measure of goodness of fit
r = corr(vec,vec_ref);

% normalize vol
vol = sf*vol;

% save normalized volume to fname_out
fs_save_mgh(vol,fname_out,M);

% save sf and r
if ~isempty(parms.fname_sf)
  save(parms.fname_sf,'sf','r');
end;

