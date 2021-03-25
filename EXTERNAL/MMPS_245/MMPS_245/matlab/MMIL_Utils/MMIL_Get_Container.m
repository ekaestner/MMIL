function [ContainerPath,ContainerDir,ContainerRootDir] = ...
  MMIL_Get_Container(RootDirs,VisitID,ContainerType)
%function [ContainerPath,ContainerDir,ContainerRootDir] = ...
% MMIL_Get_Container(RootDirs,VisitID,[ContainerType])
%
% Required:
%   RootDirs:
%    a struct which may contain the following fields:
%         home, batch, orig, raw,
%         proc, proc_dti, proc_bold,
%         fsurf, fsico,
%         orig_meg, raw_meg, proc_meg,
%         orig_pet, raw_pet, proc_pet
%    these specify the locations of batchdirs and data
%  VisitID:
%    string specifying Visit ID (name of orig data dir)
%      may be a MEG_VisitID or PET_VisitID
%
% Optional:
%   ContainerType:
%     allowed: orig, raw, proc, proc_dti, proc_bold,
%              fsurf, fsico, fsico4,
%              orig_meg, raw_meg, proc_meg,
%              orig_pet, raw_pet, proc_pet
%     {default = 'proc'}
%
% Note: If multiple containers found with same VisitID,
%   will use first one found
%
% Created:  03/31/09 by Don Hagler
% Last Mod: 12/05/16 by Don Hagler
%

ContainerPath = [];
ContainerDir = [];
ContainerRootDir = [];

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('ContainerType','var') || isempty(ContainerType)
  ContainerType = 'proc';
end;

% container types that have 8 digit date plus time attached as unique SessID
SessID_types = {'raw' 'proc' 'proc_dti' 'proc_bold'...
  'fsurf' 'fsico' 'fsico1' 'fsico2' 'fsico3' 'fsico4'...
  'fsico5' 'fsico6' 'raw_pet' 'proc_pet'};

if ~isempty(regexp(ContainerType,'fsico'))
  DirType = 'fsico';
else
  DirType = ContainerType;
end;
if ~isfield(RootDirs,DirType)
  error('no field %s in RootDirs',DirType);
end;

ContainerRootDir = getfield(RootDirs,DirType);
switch ContainerType
  case 'orig'
    ContainerStem = [];
  case 'raw'
    ContainerStem = 'MRIRAW_';
  case 'proc'
    ContainerStem = 'MRIPROC_';
  case 'proc_dti'
    ContainerStem = 'DTIPROC_';
  case 'proc_bold'
    ContainerStem = 'BOLDPROC_';
  case 'fsurf'
    ContainerStem = 'FSURF_';
  case 'long'
    ContainerStem = 'LONG_';
  case 'fsico'
    ContainerStem = 'FSICO_';
  case 'fsico1'
    ContainerStem = 'FSICO1_';
  case 'fsico2'
    ContainerStem = 'FSICO2_';
  case 'fsico3'
    ContainerStem = 'FSICO3_';
  case 'fsico4'
    ContainerStem = 'FSICO4_';
  case 'fsico5'
    ContainerStem = 'FSICO5_';
  case 'fsico6'
    ContainerStem = 'FSICO6_';
  case 'raw_asl'
    ContainerStem = [];
  case 'proc_asl'
    ContainerStem = 'ASLPROC_';
  case 'orig_meg'
    ContainerStem = [];
  case 'raw_meg'
    ContainerStem = 'MEGRAW_';
  case 'proc_meg'
    ContainerStem = 'MEGPROC_';
  case 'orig_pet'
    ContainerStem = [];
  case 'raw_pet'
    ContainerStem = 'PETRAW_';
  case 'proc_pet'
    ContainerStem = 'PETPROC_';
  otherwise
    fprintf('%s: ERROR: invalid ContainerType: %s\n',...
      mfilename,ContainerType);
    return;
end;

% find directories with VisitID
dirlist = dir(sprintf('%s/%s%s*',ContainerRootDir,ContainerStem,VisitID));
if isempty(dirlist), return; end;
dirlist = {dirlist.name};
if ismember(ContainerType,SessID_types)
  % make sure VisitID is immediately followed by "_", 8 digit date, and "."
  reg_pat = sprintf('%s%s_\\d{8}\\..+',ContainerStem,VisitID);
else
  % make sure ContainerStem + VisitID is entire ContainerName
  reg_pat = sprintf('%s%s$',ContainerStem,VisitID);
end;
% replace problem characters
reg_pat = regexprep(reg_pat,'\^','\\\^');
% check the format matches
n = find(~cellfun('isempty',regexp(dirlist,reg_pat)));
if isempty(n), return; end;
dirlist = dirlist(n);
if length(dirlist)>1
  fprintf('%s: WARNING: multiple Containers with VisitID %s\n',...
    mfilename,VisitID);
end;
ContainerDir = dirlist{1};
ContainerPath = [ContainerRootDir '/' ContainerDir];

if isempty(ContainerDir) || ~exist(ContainerPath,'dir')
  ContainerDir = [];
  ContainerPath = [];
end;

