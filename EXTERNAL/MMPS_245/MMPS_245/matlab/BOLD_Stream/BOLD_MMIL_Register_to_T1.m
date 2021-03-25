function [M_T1_to_BOLD,RegInfo,errcode] = BOLD_MMIL_Register_to_T1(ContainerPath,varargin)
%function [M_T1_to_BOLD,RegInfo,errcode] = BOLD_MMIL_Register_to_T1(ContainerPath,[options])
%
% Purpose: Register BOLD scans to T1-weighted volume and save register.dat files
%
% Usage:
%  BOLD_MMIL_Register_to_T1(ContainerPath,'key1', value1,...);
%
% Required Input
%  ContainerPath: full path of BOLDROC Container Directory with BOLD scans (mgh format)
%
% Optional Input:
%   'snums' - vector of scan numbers to process (if empty, do all)
%     {default = []}
%   'refsnum' - reference scan number to register to T1 volume
%     if empty, will register reference as determined by BOLD_MMIL_Get_ScanInfo
%      {default = []}
%   'infix' - BOLD file name infix
%      e.g. '', 'corr', 'corr_resBOLD', 'corr_regT1'
%     {default = []}
%   'T1ContainerPath': ContainerPath for T1-weighted image volume
%     if empty, will use ContainerPath
%     {default = []}
%   'fname_T1': file name of T1 volume
%     if empty, will look for MPR_res.mgh or hiFA_res.mgh
%     (preference depends on T1type)
%     {default = []}
%   'T1type': which T1 series to reg the DTI to
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%   'FSContainerPath': full path of directory containing freesurfer recon
%     will use mri/nu.mgz for fname_T1
%     (ignored if fname_T1 is also supplied)
%     {default = []}
%   'bbregflag': [0|1] use FreeSurfer's bbregister (requires FSContainerPath)
%     {default = 0}
%   'forceflag': run calculations even if output files exist
%     {default = 0}
%
% Created:  02/06/07 by Don Hagler
% Rcnt Mod: 04/12/12 by Don Hagler
% Last Mod: 09/10/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output
errcode = 0;
M_T1_to_BOLD = [];
RegInfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'refsnum',[],[],...
  'infix',[],[],...
...
  'fname_T1',[],[],...
  'T1ContainerPath',[],[],...
  'T1type',2,[0,1,2,3],... 
  'FSContainerPath',[],[],...
  'bbregflag',false,[false true],...
  'forceflag',false,[false true],...
... % for jpdf
  'cleanup_flag',true,[false true],...
... % for creating brain mask if does not exist
  'mask_smooth',15,[0,Inf],...
...
  'fnamestem','BOLD',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

parms.T2_type = ['T2_' parms.fnamestem];

if isempty(parms.T1ContainerPath), parms.T1ContainerPath = ContainerPath; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
  'snums',parms.snums,'fnamestem',parms.fnamestem);
if errcode || isempty(ScanInfo), return; end;

% choose ref scan num
if isempty(parms.refsnum) % use SessInfo ref
  parms.refsnum = SessInfo.regT1_ref;
end;

% check for bad scan nums
if ~ismember(parms.refsnum,SessInfo.snums_valid)
  fprintf('%s: ERROR: bad BOLD scan num (%d)\n',mfilename,parms.refsnum);
  errcode = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose reference BOLD scan
fname_BOLD = sprintf('%s/%s',ContainerPath,ScanInfo(parms.refsnum).fstem);
if ~isempty(parms.infix), fname_BOLD = [fname_BOLD '_' parms.infix]; end;
fname_BOLD = [fname_BOLD parms.ext];
if ~exist(fname_BOLD,'file')
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname_BOLD);
  errcode = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: using %s scan %d to register to T1...\n',...
  mfilename,parms.fnamestem,parms.refsnum);
tags = {'fname_T1','T1ContainerPath','T1type',...
  'FSContainerPath','bbregflag','ext',...
  'T2_type','forceflag','cleanup_flag','mask_smooth'};
args = mmil_parms2args(parms,tags);
[M_T1_to_BOLD,RegInfo,errcode] = ...
  MMIL_Register_T2_to_T1(fname_BOLD,args{:});

