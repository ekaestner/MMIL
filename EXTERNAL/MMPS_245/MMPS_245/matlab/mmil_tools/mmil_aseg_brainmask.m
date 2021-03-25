function mmil_aseg_brainmask(fspath,varargin)
%function mmil_aseg_brainmask(fspath,[options])
%
% Purpose: create a filled brainmask from a FreeSurfer aseg volume
%
% Usage: mmil_aseg_brainmask(fspath,'key1', value1,...);
%
% Required Input:
%  fspath: full path of freesurfer recon
%
% Optional Input:
%  'fname_mask': output file name for brain mask
%    {default = 'brainmask.aseg.mgz'}
%  'brain_flag': [0|1] whether output should be a 
%                masked brain (1) or just a mask (0)
%    {default = 1}
%  'fill_flag': [0|1] whether to do smooth and threshold steps to fill holes
%    If 0, smooth and thresh parameters are set to 0
%    {default = 1}
%  'edit_flag': [0|1] whether to use aseg.mgz (0) or aseg_edit.mgz (1)
%    {default = 0}
%  'smooth1': smoothing sigma (voxels) for initial fill step
%    {default = 20}
%  'thresh1': threshold applied to mask after first smoothing step
%    {default = 0.5}
%  'smooth2': smoothing sigma (voxels) for second dilation step
%    {default = 40}
%  'thresh2': threshold applied to mask after second smoothing step
%    {default = 0.2}
%  'smooth3': smoothing sigma (voxels) for third dilation step
%    {default = 10}
%  'forceflag': [0|1] whether to overwrite fname_mask if it already exists
%    {default = 0}
%
% See Also:
%   mmil_dilate_mask: to create mask from directly specified aseg volume
%   mmil_masksurf: to create surface mesh from mask
%   mmil_dct_brainmask: to create mask from image volume using dct morph to atlas
%
% Created:  02/02/07 by Don Hagler
% Last Mod: 02/07/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_mask','brainmask.aseg.mgz',[],...
  'brain_flag',true,[false true],...
  'fill_flag',true,[false true],...
  'edit_flag',false,[false true],...
  'smooth1',20,[0,100],...
  'thresh1',0.5,[0,1],...
  'smooth2',40,[0,100],...
  'thresh2',0.2,[0,1],...
  'smooth3',10,[0,100],...
  'forceflag',false,[false true],...
});

if isempty(parms.fname_mask)
  error('fname_mask is empty');
end;
if exist(parms.fname_mask,'file') && ~parms.forceflag
  return;
end;

if parms.edit_flag
  fname_aseg = sprintf('%s/mri/aseg_edit.mgz',fspath);
else
  fname_aseg = sprintf('%s/mri/aseg.mgz',fspath);
end;
if ~exist(fname_aseg,'file')
  error('file %s not found',fname_aseg);
end;

if parms.brain_flag
  fname_T1 =  sprintf('%s/mri/nu.mgz',fspath);
  if ~exist(fname_T1,'file')
    error('file %s not found',fname_T1);
  end;
end;

if ~parms.fill_flag
  parms.smooth1 = 0;
  parms.thresh1 = 0;
  parms.smooth2 = 0;
  parms.thresh2 = 0;
  parms.smooth3 = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_mask,'file') || parms.forceflag
  volmask = mmil_dilate_mask([],...
    'fname_in',fname_aseg,...
    'smooth1',parms.smooth1,'thresh1',parms.thresh1,...
    'smooth2',parms.smooth2,'thresh2',parms.thresh2,...
    'smooth3',parms.smooth3,...
    'forceflag',parms.forceflag);
  if parms.brain_flag
    % create masked T1 image
    vol_T1 = ctx_load_mgh(fname_T1);
    volmask.imgs = vol_T1.imgs.*volmask.imgs;
  end;
  ctx_save_mgh(volmask,parms.fname_mask);
end;

