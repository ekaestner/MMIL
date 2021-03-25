function mmil_warp_wmmask_to_atlas(FSContainerPath,varargin);
%function mmil_warp_wmmask_to_atlas(FSContainerPath,[options]);
%
% Required Parameters:
%  FSContainerPath: full path of FreeSurfer recon container
%
% Optional Parameters:
%  'outdir': output directory relative to FreeSurfer container dir
%    unless full path
%    {default = 'atlas'}
%  'forceflag': [0|1] overwrite existing files
%    {default = 0}
%
% Created:  10/21/14 by Don Hagler
% Last Mod: 10/21/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir','atlas',[],...
  'forceflag',false,[false true],...
...
  'wm_codes',[2,41],[],...
  'fname_T1','nu.mgz',[],...
  'extra_files',{'orig.mgz','wm.mgz','brainmask.mgz'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that Freesurfer recon exists
if ~exist(FSContainerPath,'dir')
  error('Freesurfer recon path %s not found',FSContainerPath);
end;

% create output directory
if mmil_isrelative(parms.outdir)
  outdir = sprintf('%s/%s',FSContainerPath,parms.outdir);
else
  outdir = parms.outdir;
end;
mmil_mkdir(outdir);

% check cortical ribbon file
fname_ribbon = sprintf('%s/mri/ribbon.mgz',FSContainerPath);
if ~exist(fname_ribbon,'file')
  error('file %s not found',fname_ribbon);
end;

% check T1-weighted structural image
fname_T1 = sprintf('%s/mri/%s',FSContainerPath,parms.fname_T1);
if ~exist(fname_T1,'file')
  error('file %s not found',fname_T1);
end;

% check extra files
for i=1:length(parms.extra_files)
  fname_extra = sprintf('%s/mri/%s',FSContainerPath,parms.extra_files{i});
  if ~exist(fname_extra,'file')
    error('file %s not found',fname_extra);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create white matter mask
fname_mask = sprintf('%s/wmmask.mgz',outdir);
fs_aseg_mask(fname_ribbon,'fname_mask',fname_mask,...
  'aseg_codes',parms.wm_codes,'forceflag',parms.forceflag);

% register intensity normalized T1-weighted MRI to atlas
%   and warp volume to atlas
mmil_warp_to_atlas(fname_T1,'fname_T1',fname_T1,...
  'outdir',outdir,'forceflag',parms.forceflag);

% warp wm mask to atlas
mmil_warp_to_atlas(fname_mask,'fname_T1',fname_T1,...
  'outdir',outdir,'forceflag',parms.forceflag);

% warp extra images to atlas
for i=1:length(parms.extra_files)
  fname_extra = sprintf('%s/mri/%s',FSContainerPath,parms.extra_files{i});
  mmil_warp_to_atlas(fname_extra,'fname_T1',fname_T1,...
  'outdir',outdir,'forceflag',parms.forceflag);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

