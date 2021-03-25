function dti_create_asegmask(fspath,varargin)
%function dti_create_asegmask(fspath,[options])
%
% Required Input:
%   fspath: full path of freesurfer recon
%
% Optional Parameters:
%   'outdir': output directory
%     {default = pwd}
%   'M_T1_to_DTI': registration matrix between DTI data and high-res T1 data
%     if not empty, will be used to resample aseg to DTI resolution
%     {default = []}
%   'fname_DTI': full path name of DTI data file
%     if supplied, volsz_DTI and M_DTI are ignored
%     {default = []}
%   'volsz_DTI': dimensions of DTI data
%     {default = [128 128 128]}
%   'M_DTI': vox2ras matrix of DTI data
%     {default = [-2,0,0,130;0,-2,0,130;0,0,-2,130;0,0,0,1]}
%   'smooth1': smoothing sigma (voxels) for initial fill step
%     {default = 20}
%   'thresh1': threshold applied to mask after first smoothing step
%     {default = 0.5}
%   'smooth2': smoothing sigma (voxels) for second dilation step
%     {default = 40}
%   'thresh2': threshold applied to mask after second smoothing step
%     {default = 0.2}
%   'smooth3': smoothing sigma (voxels) for third dilation step
%     {default = 10}
%   'nii_flag': [0|1] convert brain mask to nii format
%     {default = 0}
%   'nii_out_orient': output orientation for nii file (e.g. LAS, LPS, etc.)
%     If empty, use native
%     {default = []}
%   'xcg_flag': [0|1] create mask excluding CSF and cortical gray matter
%     {default = 1}
%   'xcg_codes': aseg ROI codes for CSF and cortical gray matter
%     {default = [0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63]}
%   'res_outfix': output string added to images resampled to DTI space
%     {default = 'resDTI'}
%   'interpm' - interpolation method
%      0 = nearest neighbor, 1 = linear, 2 = cubic
%      3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%     {default = 2}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Created:  10/15/12 by Don Hagler
% Last Mod: 01/15/13 Don Hagler
%

% based on DTI_MMIL_Create_Aseg_Mask, created 08/10/12 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: fname_aseg and fname_T1 parameters, pass to mmil_aseg_brainmask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters and check for problems
parms = check_input(fspath,varargin);

% create brain mask from aseg, resample to DTI, convert to nii
create_brain_mask(parms);

% create xcg mask from aseg, resample to DTI, convert to nii
if parms.xcg_flag
  create_xcg_mask(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fspath,options)
  parms_filter = {...
    'fspath',fspath,[],...
  ...
    'outdir',pwd,[],...
    'M_T1_to_DTI',[],[],...
    'fname_DTI',[],[],...
    'volsz_DTI',[128 128 128],[],...
    'M_DTI',[-2,0,0,130;0,-2,0,130;0,0,-2,130;0,0,0,1],[],...
    'smooth1',20,[0,100],...
    'thresh1',0.5,[0,1],...
    'smooth2',40,[0,100],...
    'thresh2',0.2,[0,1],...
    'smooth3',10,[0,100],...
    'nii_flag',false,[false true],...
    'nii_out_orient',[],[],...
    'xcg_flag',false,[false true],...
    'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
    'res_outfix','resDTI',[],...
    'interpm',2,[0:5],...
    'forceflag',false,[false true],...
  ...
    'aseg_name','aseg',[],...
    'xcg_suffix','xcg',[],...
  ...
    'mask_tags',{'fname_mask','brain_flag','fill_flag',...
                 'edit_flag','smooth1','thresh1','smooth2',...
                 'thresh2','smooth3','forceflag'},[],...
    'resamp_tags',{'M_reg','M_ref','nvox_ref','interpm'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % check that input exists
  if ~exist(parms.fspath,'dir')
    error('FreeSurfer recon dir %s not found',parms.fspath);
  end;
  % check that aseg exists
  parms.fname_aseg = sprintf('%s/mri/%s.mgz',parms.fspath,parms.aseg_name);
  if ~exist(parms.fname_aseg,'file')
    error('file %s not found',parms.fname_aseg);
  end
  if ~isempty(parms.M_T1_to_DTI)
    % get volsz and M from DTI file
    if ~isempty(parms.fname_DTI)
      [parms.M_DTI,parms.volsz_DTI] = mmil_load_mgh_info(parms.fname_DTI,...
        parms.forceflag,parms.outdir);
    end;
    % set parameters for resampling volumes
    parms.M_reg = inv(parms.M_T1_to_DTI);
    parms.M_ref = parms.M_DTI;
    parms.nvox_ref = parms.volsz_DTI;
  end;
  % create output directory if necesary
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_brain_mask(parms)
  % generate brain mask
  fname_brainmask = sprintf('%s/%s_brainmask.mgz',parms.outdir,parms.aseg_name);
  parms.fname_mask = fname_brainmask;
  args = mmil_parms2args(parms,parms.mask_tags);
  mmil_aseg_brainmask(parms.fspath,args{:});
  % resample brain mask
  if ~isempty(parms.M_T1_to_DTI)
    fname_brainmask_res = sprintf('%s/%s_brainmask_%s.mgz',...
      parms.outdir,parms.aseg_name,parms.res_outfix);
    if ~exist(fname_brainmask_res,'file') || parms.forceflag
      [vol,M] = fs_load_mgh(fname_brainmask);
      args = mmil_parms2args(parms,parms.resamp_tags);
      [vol,M] = mmil_resample_vol(vol,M,args{:});
      fs_save_mgh(vol,fname_brainmask_res,M);
    end
    fname_brainmask = fname_brainmask_res;
  end;
  % convert to nii format
  if parms.nii_flag
    fname_nii = regexprep(fname_brainmask,'mgz$','nii');
    fs_mri_convert(fname_brainmask_res,fname_nii,...
      'out_orient',parms.nii_out_orient,'forceflag',parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_xcg_mask(parms)
  % generate mask excluding CSF and gray matter
  fname_xcg = sprintf('%s/%s_%s_mask.mgz',...
    parms.outdir,parms.aseg_name,parms.xcg_suffix);
  fs_aseg_mask(parms.fname_aseg,'fname_mask',fname_xcg,...
    'exclude_flag',1,'aseg_codes',parms.xcg_codes,...
    'forceflag',parms.forceflag);
  % resample to DTI resolution
  if ~isempty(parms.M_T1_to_DTI)
    fname_xcg_res = sprintf('%s/%s_%s_mask_%s.mgz',...
      parms.outdir,parms.aseg_name,parms.res_outfix,parms.xcg_suffix);
    if ~exist(fname_xcg_res,'file') || parms.forceflag
      [vol,M] = fs_load_mgh(fname_xcg);
      args = mmil_parms2args(parms,parms.resamp_tags);
      [vol,M] = mmil_resample_vol(vol,M,args{:});
      vol = 1.0 * (vol == 1);
      fs_save_mgh(vol,fname_xcg_res,M);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
