function MMIL_Long_Mask_DTI(fname_DTI,fname_T1_mask,fname_DTI_mask,forceflag)
%function MMIL_Long_Mask_DTI(fname_DTI,fname_T1_mask,fname_DTI_mask,forceflag)
%
% Required Input:
%  fname_DTI: full path file name of DTI image
%  fname_T1_mask: full path file name of brainmask created from structural
%  fname_DTI_mask: full path file name ot output DTI mask image
%
% Optional Input:
%   forceflag: overwrite if it already exists
%
% Created:  10/15/09 by Don Hagler
% Last Mod: 10/15/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,3)) return; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag = 0; end;

thresh = 1e-4;
smooth = 5;

if ~exist(fname_T1_mask,'file')
  error('file %s not found',fname_T1_mask);
end
if ~exist(fname_DTI,'file')
  error('file %s not found',fname_DTI);
end

if ~exist(fname_DTI_mask,'file') | forceflag
  vol_T1_mask = ctx_load_mgh(fname_T1_mask);
  [vol_DTI,mr_parms] = ctx_load_mgh(fname_DTI);
  vol_DTI_mask = vol_T1_mask;
  vol_DTI_mask.imgs(find(vol_DTI.imgs < thresh)) = 0;
  if smooth>0
    vol_DTI_mask = mmil_dilate_mask(vol_DTI_mask,...
      'smooth1',0,'smooth2',0,'smooth3',smooth);
  end;
  ctx_save_mgh(vol_DTI_mask,fname_DTI_mask,mr_parms);
end;

