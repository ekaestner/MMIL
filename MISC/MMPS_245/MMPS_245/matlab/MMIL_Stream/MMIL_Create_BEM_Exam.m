function MMIL_Create_BEM_Exam(ContainerDir,RootDirs,varargin)
%function MMIL_Create_BEM_Exam(ContainerDir,RootDirs,[options])
%
% Purpose: create BEM surfaces and copy to fsurf and fsico containers
%
% Required Parameters:
%  ContainerDir: Input MRIPROC Container Dir
%  RootDirs: a struct which must contain the fields fsurf and fsico
%    specifying the container root dirs
%
% Optional Parameters:
%  'BEMtype': source of BEM surfaces (e.g. 'PD','PD_NFT','T1','aseg')
%    for 'aseg' type, 'T1' is used for outer skull and outer scalp
%      and dilated FreeSurfer aseg is used for inner skull
%    for 'PD_NFT' type, surfaces from PD are used in combination with
%      FreeSurfer aseg and NFT mesh generation to create 4-shell BEM
%    {default = 'T1'}
%  'ico': icosahedral order number (1-7)
%    Will copy BEM surfaces to fsico directory if found
%    {default = 4}
%  'make_forceflag': [0|1] whether to overwrite existing output
%    in proc container
%    {default = 0}
%  'copy_forceflag': [0|1] whether to overwrite existing output
%    in fsurf and fsico containers
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output
%    overrides both 'make_forceflag' and 'copy_forceflag'
%    {default = 0}
%
% Created:  02/16/11 by Don Hagler
% Last Mod: 06/16/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'BEMtype','T1',{'PD','PD_NFT','T1','aseg'},...
  'ico',4,[1:7],...
  'make_forceflag',false,[false true],...
  'copy_forceflag',false,[false true],...
  'forceflag',false,[false true],...
...
  'proc_typestr','MRIPROC',[],...
  'fsurf_typestr','FSURF',[],...
  'fsico_typestr','FSICO',[],...  
});

if parms.forceflag
  parms.make_forceflag = 1;
  parms.copy_forceflag = 1;
end;

if parms.ico~=7
  parms.fsico_typestr = sprintf('%s%d',parms.fsico_typestr,parms.ico);
end;

if ~isfield(RootDirs,'fsico') | isempty(RootDirs.fsico)
  RootDirs.fsico = RootDirs.fsurf;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if proc container exists
ContainerPath = [RootDirs.proc '/' ContainerDir];
if isempty(ContainerDir) || ~exist(ContainerPath,'dir')
  fprintf('%s: ERROR: Processed Container %s not found\n',...
    mfilename,ContainerPath);
  return;
end;

% check if fsurf container exists
fsurf_ContainerDir = regexprep(ContainerDir,...
  parms.proc_typestr,parms.fsurf_typestr);
fsurf_ContainerPath = [RootDirs.fsurf '/' fsurf_ContainerDir];
if isempty(fsurf_ContainerDir) || ~exist(fsurf_ContainerPath,'dir')
  fprintf('%s: WARNING: FreeSurfer Container %s not found\n',...
    mfilename,fsurf_ContainerPath);
  fsurf_ContainerPath = [];
end;

% check if fsico container exists
fsico_ContainerDir = regexprep(ContainerDir,...
  parms.proc_typestr,parms.fsico_typestr);
fsico_ContainerPath = [RootDirs.fsico '/' fsico_ContainerDir];
if isempty(fsico_ContainerDir) || ~exist(fsico_ContainerPath,'dir')
  fprintf('%s: WARNING: FreeSurfer Ico Container %s not found\n',...
    mfilename,fsurf_ContainerPath);
  fsico_ContainerPath = [];
end;

MMIL_Make_BEM_Surfs(ContainerPath,fsurf_ContainerPath,parms.make_forceflag);

if ~isempty(fsurf_ContainerPath)
  MMIL_Copy_BEM_Surfs(ContainerPath,fsurf_ContainerPath,...
    parms.BEMtype,parms.copy_forceflag);
end;

if ~isempty(fsico_ContainerPath)
  MMIL_Copy_BEM_Surfs(ContainerPath,fsico_ContainerPath,...
    parms.BEMtype,parms.copy_forceflag);
end;

