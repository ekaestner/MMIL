function [ProjInfo,StudyInfo,RootDirs]=MMIL_Quick_Check_ProjID(ProjID,varargin)
%function [ProjInfo,StudyInfo,RootDirs]=MMIL_Quick_Check_ProjID(ProjID,[options])
%
% Purpose: Check validity of ProjID, StudyInfo, and RootDirs
%   If StudyInfo and RootDirs are supplied, ProjID is ignored
%   If not, ProjID is used to look for ProjInfo and VisitInfo files
%     (e.g. '/home/user/ProjInfo/MMIL_ProjInfo.csv'
%           '/home/user/ProjInfo/ProjID/ProjID_VisitInfo.csv' )
%
% Required Input:
%   ProjID: Project ID string
%     (may be empty if RootDirs is supplied)
%
% Optional Input:
%  'RootDirs': struct containing locations of root data dirs
%     (may be empty if ProjID is not)
%    {default: []}
%  'StudyInfo': struct array containing info for each subject
%     (may be empty if ProjID or RootDirs is supplied)
%    {default: []}
%  'VisitIDs': cell array of VisitIDs to be included in StudyInfo
%    if not supplied, include all VisitIDs
%    {default = []}
%  'SubjIDs': cell array of SubjIDs to be included in StudyInfo
%    if not supplied, include all SubjIDs
%    {default = []}
%  'user': user name of account with StudyInfo (i.e. ProjInfo, VisitInfo)
%    If not supplied, will use $USER environment variable
%    {default = []}
%  'ignore_VisitInfo_flag': [0|1] ignore VisitInfo csv file
%     when creating StudyInfo (get VisitIDs from containers in RootDirs.orig)
%    {default = 0}
%
% Output:
%   ProjInfo: struct containing parameters for this ProjID
%   StudyInfo: struct array with info for each "study" (scan session)
%   RootDirs: struct containing locations of data directories
%
% Created:  11/22/16 by Don Hagler
% Last Mod: 11/22/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'VisitIDs',[],[],...
  'SubjIDs',[],[],...
  'user',[],[],...
  'ignore_VisitInfo_flag',false,[false true],...
...
  'numvec_tags',[],[],...
...
  'info_tags',{'VisitIDs','SubjIDs','RootDirs','ignore_VisitInfo_flag','user','numvec_tags'},[],...
});

ProjInfo = [];
StudyInfo = [];
RootDirs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(ProjID)
  if isempty(parms.StudyInfo)
    args = mmil_parms2args(parms,parms.info_tags);
    [StudyInfo,RootDirs,ProjInfo] = MMIL_Quick_StudyInfo(ProjID,args{:});
  else
    try
      [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID,parms.numvec_tags,parms.user);
    catch me
      fprintf('%s: WARNING: %s\n',mfilename,me.message);
    end;
  end;
end;

if isempty(RootDirs), error('RootDirs is empty'); end;
if isempty(StudyInfo), error('StudyInfo is empty'); end;

