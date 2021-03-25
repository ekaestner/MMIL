function DTI_MMIL_CSD_Tracto_Exam(ContainerPath,varargin)
%function DTI_MMIL_CSD_Tracto_Exam(ContainerPath,varargin)
%
% Purpose:
%   Runs CSD tractography from Alexander Leeman's ExploreDTI
%
% Usage:
%  DTI_MMIL_CSD_Tracto_Exam(ContainerPath,'key1', value1,...);
%
% Required Input Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%                  (mgh/mgz format)
%
% Optional Parameters specifying which data to load:
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr','corr_resDTI','corr_regT1'
%     {default = 'corr_regT1'}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'min_bval': minimum b value a scan must have to be included
%     {default = 1}
%   'flex_flag': [0|1] DTI_flex scans included
%     {default = 0}
%
% Optional Parameters for tracking:
%   'see_point_sampling': Seed point sampling in x, y, and z direction (in mm)
%     {default = [2 2 2]}
%   'step_size': Step size (in mm)
%     {default = 1}
%   'FOD_thresh': FOD threshold
%     {default = 0.1}
%   'angle threshold': angle deviation threshold
%     {default = 45}
%   'fiber_length_range': min and max fiber length (in mm)
%     {default = [50 500]}
%   'max_order': highest order spherical harmonic (even integers)
%     {default = 4}
%   'FA_thresh': minimal FA to determine response function
%     {default = 0.8}
%
% Optional Parameters for output
%  'outdir': output directory (full path or relative to ContainerPath)
%    {default: 'CSD_Tracto'}
%  'forceflag': [0|1] whether to run conversions even if output files exist
%    {default: 0}
%
% Created:  08/28/09 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'infix',[],[],...
  'revflag',0,[0,1,2],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
...
  'seed_point_sampling',[2 2 2],[],...
  'step_size',1,[],...
  'FOD_thresh',0.1,[],...
  'angle_thresh',45,[],...
  'fiber_length_range',[50 500],[],...
  'max_order',4,[2:2:100],...
  'FA_thresh',0.8,[0 1],...
...
  'outdir',ContainerPath,[],...
  'outstem','CSD',[],...
  'orient',[],[],...
  'forceflag',false,[false true],...
...
  'min_nb0',1,[],...
  'min_ndirs',6,[],...
...
  'data_tags',{'snums','infix','revflag',...
               'min_nb0','min_ndirs','min_bval','flex_flag'},[],...
  'CSD_tags',{'outdir','outstem','orient','forceflag',...
              'seed_point_sampling','step_size','FOD_thresh',...
              'angle_thresh','fiber_length_range','max_order',...
              'FA_thresh'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_CSD = sprintf('%s/%s_data_Tracts_CSD.mat',parms.outdir,parms.outstem);
if ~exist(fname_CSD,'file') || parms.forceflag
  % load data
  args = mmil_parms2args(parms,parms.data_tags);
  [vol,M,qmat,bvals]=DTI_MMIL_Load_Data(ContainerPath,args{:});
  % run CSD tractography  
  args = mmil_parms2args(parms,parms.CSD_tags);
  dti_CSD_tracto(vol,M,qmat,bvals,args{:});
end;

