function [M_T1_to_DTI,RegInfo,errcode] = DTI_MMIL_Register_to_T1(ContainerPath,varargin)
%function [M_T1_to_DTI,RegInfo,errcode] = ...
%  DTI_MMIL_Register_to_T1(ContainerPath,[options])
%
% Usage:
%  DTI_MMIL_Register_to_T1(ContainerPath,'key1', value1,...);
%
% Required Parameters:
%   ContainerPath: full path of DTIPROC directory containing
%     processed diffusion data (mgz format)
%
% Optional Parameters specifying DTI image:
%   'fname_DTI': full path name of 4d diffusion data or b=0 volume
%     If empty, will choose reference DTI file based on the following
%     parameters using DTI_MMIL_Get_ScanInfo's reg_ref_for or reg_ref_rev
%     {default = []}
%   'snums': list of scan numbers to use
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr'
%     {default = []}
%   'min_ndirs': require at least this many diffusion directions
%     for processing
%     {default = 6}
%   'min_bval': minimum b value a scan must have to be processed
%     {default = 1}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': completely ignore scans with fewer than this
%     number of b=0 images
%     {default = 1}
%   'revflag': [0|1|2] specify whether to process non-rev or rev data
%     0: register forward phase-encode polarity data
%     1: register reverse phase-encode polarity data
%     2: register either forward or reverse data,
%        depending on which is more common
%     {default = 2}
%
% Optional Parameters specifying T1 image:
%   'T1ContainerPath': ContainerPath for T1-weighted image volume
%     if empty, will use ContainerPath
%     {default = []}
%   'fname_T1': file name of T1-weighted image volume
%     if empty, will look for MPR_res.mgz or hiFA_res.mgz
%     (preference depends on T1type)
%     {default = []}
%   'T1type': which T1 series to reg the DTI to
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%   'FSContainerPath': full path of directory containing freesurfer recon
%     will use mri/nu.mgz for fname_T1
%     Required if bbregflag=1, otherwise ignored if fname_T1 is supplied
%     {default = []}
%
% Other Optional Parameters:
%   'fname_reg': output file name for mat file with RegInfo struct
%     If empty, will generate based on fname_DTI
%     {default = []}
%   'bbegflag': [0|1] use FreeSurfer's bbregister (requires FSContainerPath)
%     {default = 0}
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%   'forceflag': run calculations even if output files exist
%     {default = 0}
%
% Created:  01/01/07 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output
errcode = 0;
M_T1_to_DTI = [];
RegInfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_DTI',[],[],...
  'snums',[],[],...
  'infix',[],[],...
  'min_nb0',1,[],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
  'revflag',2,[0,1,2],...
...
  'fname_T1',[],[],...
  'T1ContainerPath',[],[],...
  'T1type',2,[0,1,2,3],... 
  'FSContainerPath',[],[],...
  'fname_reg',[],[],...
  'bbregflag',false,[false true],...
  'forceflag',false,[false true],...
... % for jpdf
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],... 
  'cleanup_flag',true,[false true],...
... % for creating brain mask if does not exist
  'mask_smooth',15,[0,Inf],...
...
  'fnamestem','DTI',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

parms.T2_type = ['T2_' parms.fnamestem];

if isempty(parms.T1ContainerPath), parms.T1ContainerPath = ContainerPath; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose DTI file
ref_snum = [];
if isempty(parms.fname_DTI)
  % load ContainerInfo, get DTI scan info, determine valid scans, reference scans
  tags = {'min_nb0','min_ndirs','min_bval','flex_flag','revflag'};
  args = mmil_parms2args(parms,tags);
  [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
  if errcode ~= 0, return; end;
  parms.fname_DTI = sprintf('%s/%s%d',...
    ContainerPath,parms.fnamestem,SessInfo.regT1_ref);
  if SessInfo.revflag
    parms.fname_DTI = [parms.fname_DTI '_rev'];
  end;
  if ~isempty(parms.infix)
    parms.fname_DTI = [parms.fname_DTI '_' parms.infix];
  end;
  parms.fname_DTI = [parms.fname_DTI parms.ext];
end;
if ~exist(parms.fname_DTI,'file')
  fprintf('%s: ERROR: DTI file %s not found\n',mfilename,parms.fname_DTI);
  errcode = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: using %s scan %d to register to T1...\n',...
  mfilename,parms.fnamestem,SessInfo.regT1_ref);
tags = {'fname_T1','T1ContainerPath','T1type',...
  'FSContainerPath','bbregflag',...
  'atlasdir', 'atlasname','T2_type','cleanup_flag',...
  'mask_smooth','ext','forceflag'};
args = mmil_parms2args(parms,tags);
[M_T1_to_DTI,RegInfo,errcode] = ...
  MMIL_Register_T2_to_T1(parms.fname_DTI,args{:});

