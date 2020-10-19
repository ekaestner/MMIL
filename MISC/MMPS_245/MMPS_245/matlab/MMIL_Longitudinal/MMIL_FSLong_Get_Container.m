function [ContainerPath,ContainerDir,ContainerRootDir] = ...
  MMIL_FSLong_Get_Container(RootDirs,SubjID,ContainerType)
%function [ContainerPath,ContainerDir,ContainerRootDir] = ...
% MMIL_FSLong_Get_Container(RootDirs,SubjID,[ContainerType])
%
% Required:
%   RootDirs:
%    a struct which may contain the following fields:
%         fsurf, fslong
%  SubjID:
%    string specifying Subject ID 

% Optional:
%   ContainerType:
%     allowed: fsurf, fslong
%     {default = 'fslong'}
%
% Created:  05/05/16 by Sean Hatton
% Last Mod: 05/12/16 by Sean Hatton
%

ContainerPath = [];
ContainerDir = [];
ContainerRootDir = [];

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('ContainerType','var') || isempty(ContainerType)
  ContainerType = 'fslong';
end;

SessID_types = {'fslong'};

if ~isfield(RootDirs,ContainerType)
  error('no field %s in RootDirs',ContainerType);
end;

ContainerRootDir = getfield(RootDirs,ContainerType);
switch ContainerType
  case 'orig'
    ContainerStem = [];
  case 'fslong'
    ContainerStem = 'FSLONG_';
  otherwise
    fprintf('%s: ERROR: not a FreeSurfer ContainerType: %s\n',...
      mfilename,ContainerType);
    return;
end;

% find directories with SubjID
dirlist = dir(sprintf('%s/%s%s*',ContainerRootDir,ContainerStem,SubjID));
if isempty(dirlist), return; end;
dirlist = {dirlist.name};
reg_pat = sprintf('%s%s_base_v*',ContainerStem,SubjID);
% replace problem characters
reg_pat = regexprep(reg_pat,'\^','\\\^');
% check the format matches
n = find(~cellfun('isempty',regexp(dirlist,reg_pat)));
if isempty(n), return; end;
dirlist = dirlist(n);
ContainerDir = dirlist{1};
ContainerPath = [ContainerRootDir '/' ContainerDir];

if ~exist(ContainerPath,'dir')
  ContainerDir = [];
  ContainerPath = [];
end;

