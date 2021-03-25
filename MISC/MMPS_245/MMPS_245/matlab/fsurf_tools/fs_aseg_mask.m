function [vol_mask,M] = fs_aseg_mask(fname_aseg,varargin)
%function [vol_mask,M] = fs_aseg_mask(fname_aseg,[options])
%
% Description: A mask is created from a vector of aseg codes
%
% Required Input:
%   fname_aseg: full path to input segmentation volume (e.g. aseg.mgz)
%
% Optional Input:
%   'fname_mask': full or relative path to output segmentation volume
%     if empty, will not save file
%     {default = []}
%   'aseg_codes':  FreeSurfer aseg code numbers to include in mask
%     If empty, include any non-zero code (creates a brain mask)
%     {default: []}
%   'exclude_flag': create reverse mask so voxels with aseg_codes
%     are assigned 0, all others 1
%     {default = 0}
%   'forceflag': [0|1] toggle overwrite of existing output file
%     {default = 0}
%
% Output:
%   vol_mask: 3D mask volume
%   M: vox2ras matrix
%
% Created:  08/24/09 by Don Hagler
% Last Mod: 10/15/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

vol_mask = []; M = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_aseg',fname_aseg,[],...
  'fname_mask',[],[],...
  'aseg_codes',[],[],...
  'exclude_flag',false,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input and output files
if ~exist(fname_aseg,'file')
  error('file %s not found',parms.fname_aseg);
end;
if ~isempty(parms.fname_mask) && ...
   exist(parms.fname_mask,'file') && ~parms.forceflag
  if nargout>0
    [vol_mask,M] = fs_load_mgh(parms.fname_mask);
  end;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load input segmentation volume
[vol_aseg,M,mr_parms,volsz] = fs_load_mgh(parms.fname_aseg);
if isempty(parms.aseg_codes)
  parms.aseg_codes = unique(vol_aseg(vol_aseg>0));
end;

% create mask from set of aseg codes
vol_mask = zeros(volsz);
if parms.exclude_flag
  % do not include non-brain
  parms.aseg_codes = [0,parms.aseg_codes];
  ind = find(~ismember(vol_aseg,parms.aseg_codes));
else
  ind = find(ismember(vol_aseg,parms.aseg_codes));
end;
vol_mask(ind) = 1.0;

% save output mask volume
if ~isempty(parms.fname_mask)
  fs_save_mgh(vol_mask,parms.fname_mask,M,mr_parms);
end;

