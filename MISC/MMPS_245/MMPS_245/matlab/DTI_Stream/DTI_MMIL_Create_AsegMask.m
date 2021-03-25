function errcode = DTI_MMIL_Create_AsegMask(ContainerPath,FSContainerPath,varargin)
%function errcode = DTI_MMIL_Create_AsegMask(ContainerPath,FSContainerPath,varargin)
%
% Usage:
%  errcode = DTI_MMIL_Create_AsegMasks(ContainerPath,FSContainerPath,'key1', value1,...);
%
% Required Parameters:
%   ContainerPath: full path of directory containing processed
%    diffusion data
%   FSContainerPath: full path of directory containing freesurfer recon
%
% Optional Parameters:
%   'outdir': output directory
%     {default = ContainerPath}
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
%     {deafult = 0}
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
% Optional Parameters that determine input diffusion data:
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     e.g. 'corr_resDTI' or 'corr_regT1'
%     {default = []}
%   'revflag': [0|1|2] specify whether to use sournon-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%
% Created:  08/10/12 by Don Hagler
% Last Mod: 01/15/13 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
errcode = 0;

parms = mmil_args2parms(varargin,{...
  'outdir',ContainerPath,[],...
  'smooth1',20,[0,100],...
  'thresh1',0.5,[0,1],...
  'smooth2',40,[0,100],...
  'thresh2',0.2,[0,1],...
  'smooth3',10,[0,100],...
  'nii_flag',false,[false true],...
  'nii_out_orient',[],[],...
  'xcg_flag',false,[false true],...
  'xcg_suffix','xcg',[],...
  'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
  'res_outfix','resDTI',[],...
  'interpm',2,[0:5],...
  'forceflag',false,[false true],...
  'infix',[],[],...
  'revflag',0,[0,1,2],...
...
  'aseg_name','aseg',[],...
  'reg_tags',{'infix','revflag'},[],...
  'mask_tags',{'outdir','M_T1_to_DTI','fname_DTI','volsz_DTI','M_DTI',...
               'smooth1','thresh1','smooth2','thresh2','smooth3',...
               'nii_flag','nii_out_orient','xcg_flag','xcg_codes',...
               'res_outfix','interpm','forceflag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that processed data container exists
if ~exist(ContainerPath,'dir')
  fprintf('%s: ERROR: %s not found\n',mfilename,ContainerPath);
  errcode = 1; return;
end;

% load registration info
args = mmil_parms2args(parms,parms.reg_tags);
[RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,args{:});
if errcode
  fprintf('%s: ERROR: *_regT1.mat file not found\n',mfilename);
  errcode = 1; return;
end;
parms.volsz_DTI = RegInfo.volsz_T2;
parms.M_DTI = RegInfo.M_T2;
parms.M_T1_to_DTI = RegInfo.M_T1_to_T2;

% check that Freesurfer recon exists
if ~exist(FSContainerPath,'dir')
  fprintf('%s: ERROR: %s not found\n',mfilename,FSContainerPath);
  errcode = 1; return;
end

% check that aseg exists
fname_aseg = sprintf('%s/mri/%s.mgz',FSContainerPath,parms.aseg_name);
if ~exist(fname_aseg,'file')
  fprintf('%s: ERROR: %s not found\n',mfilename,fname_aseg);
  errcode = 1; return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_args2parms(parms,parms.mask_tags);
dti_create_asegmask(FSContainerPath,args{:});

