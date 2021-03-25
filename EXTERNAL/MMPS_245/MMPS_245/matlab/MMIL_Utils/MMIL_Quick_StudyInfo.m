function [StudyInfo,RootDirs,ProjInfo]=MMIL_Quick_StudyInfo(ProjID,varargin)
%function [StudyInfo,RootDirs,ProjInfo]=MMIL_Quick_StudyInfo(ProjID,[options])
%
%  Usage: [StudyInfo,RootDirs,ProjInfo]=...
%               MMIL_Quick_StudyInfo(ProjID,'key',value,...)
%
% Required Input:
%   ProjID: project ID string
%
% Optional Parameters
%  'VisitIDs': cell array of VisitIDs to be included in StudyInfo
%    if not supplied, include all VisitIDs
%    {default = []}
%  'SubjIDs': cell array of SubjIDs to be included in StudyInfo
%    if not supplied, include all SubjIDs
%    {default = []}
%  'ignore_VisitInfo_flag': [0|1] ignore VisitInfo csv file
%     when creating StudyInfo (get VisitIDs from containers in RootDirs.orig)
%    {default = 0}
%  'RootDirs': struct containing locations of root data dirs
%     (may be empty if ProjID is not)
%    {default: []}
%  'user': user name of account with StudyInfo (i.e. ProjInfo, VisitInfo)
%    If not supplied, will use $USER environment variable
%    {default = []}
%
% Output:
%   StudyInfo: struct array with info for each "study" (scan session)
%   RootDirs: struct containing locations of data directories
%   ProjInfo: struct containing parameters for this ProjID
%
% Created:  04/03/16 by Don Hagler
% Prev Mod: 04/03/16 by Don Hagler
% Last Mod: 09/03/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'VisitIDs',[],[],...
  'SubjIDs',[],[],...
  'RootDirs',[],[],...
  'ignore_VisitInfo_flag',false,[false true],...
  'user',[],[],...
...
  'numvec_tags',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StudyInfo = [];
RootDirs = [];
ProjInfo = [];

if isempty(ProjID)
  error('ProjID is empty');
end;
if iscell(ProjID) & length(ProjID)>1
  error('multiple ProjID''s not allowed');
elseif iscell(ProjID)
  ProjID = ProjID{1};
end;

% get parameter values for this project as well as root directories
if ~isempty(parms.RootDirs)
  RootDirs = parms.RootDirs;
else
  [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID,parms.numvec_tags,parms.user);
end;

% load visit info file if exists
fname_visit = sprintf('%s/ProjInfo/%s/%s_VisitInfo.csv',...
  RootDirs.home,ProjID,ProjID);
if exist(fname_visit,'file') && ~parms.ignore_VisitInfo_flag
  fprintf('%s: creating StudyInfo from %s...\n',mfilename,fname_visit);
  StudyInfo = MMIL_Read_StudyInfo(fname_visit);
end

% get VisitIDs from orig dir
if isempty(StudyInfo)
  % get VisitIDs from orig dir
  dirlist = dir(RootDirs.orig);
  % get rid of . and .. in dirlist
  n = cellfun('isempty',regexp({dirlist.name},'^\.'));
  dirlist = dirlist(n);
  orig_dirs = {dirlist.name};
  % replace any '.' in VisitID with '_'
  %  (causes problems when inserted in job names and for ADNI)
  VisitIDs = regexprep(orig_dirs,'\.','_');
  StudyInfo = struct('orig',orig_dirs,'VisitID',VisitIDs);
end;

% restrict to specified VisitIDs
StudyInfo = restrict_info(StudyInfo,parms.VisitIDs,'VisitID');

% restrict to specified SubjIDs
StudyInfo = restrict_info(StudyInfo,parms.SubjIDs,'SubjID');

% check STRUCT_VisitID
StudyInfo = check_info(StudyInfo);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StudyInfo = restrict_info(StudyInfo,IDs,IDname)
  if ~isempty(IDs) &&...
     (isempty(StudyInfo) || isfield(StudyInfo,IDname))
    if isempty(StudyInfo)
      ind_add = [1:length(IDs)];
    else
      all_IDs = {StudyInfo.(IDname)};
      ind_keep = find(ismember(all_IDs,IDs));
      [tmp,ind_add] = setdiff(IDs,all_IDs);
      if ~isempty(ind_keep)
        % keep StudyInfo for matching IDs
        StudyInfo = StudyInfo(ind_keep);
      else
        StudyInfo = [];
      end;
      if ~isempty(ind_add)
        fprintf('%s: WARNING: %s not found in VisitInfo: %s\n',...
          mfilename,IDname,sprintf('%s ',IDs{ind_add}));
      end;  
    end;
    % add missing IDs to StudyInfo
    for i=1:length(ind_add)
      StudyInfo(end+1).(IDname) = IDs{ind_add(i)};
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StudyInfo = check_info(StudyInfo)
  for i=1:length(StudyInfo)
    if ~isfield(StudyInfo,'STRUCT_VisitID') ||...
        isempty(StudyInfo(i).STRUCT_VisitID)
      StudyInfo(i).STRUCT_VisitID = StudyInfo(i).VisitID;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

