function [ProjInfo,StudyInfo,RootDirs]=MMIL_Check_ProjID(ProjID,varargin)
%function [ProjInfo,StudyInfo,RootDirs]=MMIL_Check_ProjID(ProjID,[options])
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
%
% Optional Parameters for checking RootDirs and StudyInfo:
%  'ignore_VisitInfo_flag': [0|1] ignore VisitInfo csv file
%     when creating StudyInfo (get VisitIDs from containers in RootDirs.orig)
%    {default = 0}
%  'checkflag': [0|1] check each visit for required containers, etc.
%    and get session ID; otherwise, only get VisitID from orig dirs
%    {default = 1}
%  'qcflag': only include subjects with StudyInfo.QC = 1
%    {default = 0}
%  'modality': limit returned StudyInfo to sessions of this modality only
%    allowed values: 'MRI', 'MEG', 'PET'
%    {default = []}
%  'required_containers': cell array of container names required to exist
%      or QC flag set to 0
%    e.g. 'raw', 'proc', 'proc_dti', 'proc_bold',
%         'fsurf', 'fsico', 'meg_raw', 'meg_proc', etc.
%    if speficied, all entries will be automatically added to 'required_rootdirs'
%    {default = []}
%  'required_rootdirs': cell array of root directories required to exist
%    e.g. 'raw', 'proc', 'proc_dti', 'proc_bold',
%         'fsurf', 'fsico', 'meg_raw', 'meg_proc', etc.
%    {default = []}
%  'ico': icosahedral order number (e.g 1-7)
%    used to specify name of fsico containers 
%    {default = 7}
%
% Optional Parameters for manual QC (determines whether overal QC=1):
%   'QC_raw' : [0|1] use good raw QC if it exists for this project
%     {default = 0}
%   'QC_recon' : [0|1] use recon QC if it exists for this project
%     {default = 0}
%   'QC_dv' : [0|1] use dv QC if it exists for this project
%     {default = 0}
%   'QC_DTI' : [0|1] use DTI QC if it exists for this project
%     {default = 0}
%   'QC_BOLD' : [0|1] use BOLD QC if it exists for this project
%     {default = 0}
%   'QC_PET' : [0|1] use PET QC if it exists for this project
%     {default = 0}
%
% Output:
%   ProjInfo: struct containing parameters for this ProjID
%   StudyInfo: struct array with info for each "study" (scan session)
%   RootDirs: struct containing locations of data directories
%
% Created:  03/12/11 by Don Hagler
% Last Mod: 04/03/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'VisitIDs',[],[],...
  'SubjIDs',[],[],...
  'user',[],[],...
...
  'ignore_VisitInfo_flag',false,[false true],...
  'checkflag',true,[false true],...
  'qcflag',false,[false true],...
  'modality',[],{'MRI','MEG','PET'},...
  'required_containers',[],[],...
  'required_rootdirs',[],[],...
  'ico',7,[1:7],...
  'numvec_tags',[],[],...
...
  'QC_raw',false,[false true],...
  'QC_recon',false,[false true],...
  'QC_dv',false,[false true],...
  'QC_DTI',false,[false true],...
  'QC_BOLD',false,[false true],...
  'QC_PET',false,[false true],...
});

ProjInfo = [];
StudyInfo = [];
RootDirs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(ProjID)
  if isempty(parms.StudyInfo)
    args = MMIL_Args(parms,'MMIL_Get_StudyInfo');
    [StudyInfo,RootDirs,ProjInfo] = MMIL_Get_StudyInfo(ProjID,args{:});
  else
    try
      [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID,parms.numvec_tags,parms.user);
    catch me
      fprintf('%s: WARNING: %s\n',mfilename,me.message);
    end;
  end;
end;

if ~isempty(parms.StudyInfo) || ~isempty(parms.RootDirs)
  if ~isempty(parms.StudyInfo), StudyInfo = parms.StudyInfo; end;
  if ~isempty(parms.RootDirs), RootDirs = parms.RootDirs; end;
  args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
  [StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs,args{:}); 
end;

if isempty(RootDirs), error('RootDirs is empty'); end;
if isempty(StudyInfo), error('StudyInfo is empty'); end;

