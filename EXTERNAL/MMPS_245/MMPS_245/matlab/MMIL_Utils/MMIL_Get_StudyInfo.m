function [StudyInfo,RootDirs,ProjInfo]=MMIL_Get_StudyInfo(ProjID,varargin)
%function [StudyInfo,RootDirs,ProjInfo]=MMIL_Get_StudyInfo(ProjID,[options])
%
%  Usage: [StudyInfo,RootDirs,ProjInfo]=...
%               MMIL_Get_StudyInfo(ProjID,'key',value,...)
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
%         'fsurf', 'fsico', 'raw_meg', 'proc_meg', etc.
%    if speficied, all entries will be automatically added to 'required_rootdirs'
%    {default = []}
%  'required_rootdirs': cell array of root directories required to exist
%    e.g. 'raw', 'proc', 'proc_dti', 'proc_bold',
%         'fsurf', 'fsico', 'raw_meg', 'proc_meg', etc.
%    {default = []}
%  'ico': icosahedral order number (e.g 1-7)
%    used to specify names of fsico containers
%    {default = 7}
%  'numvec_tags': cell array of column headers specifying which columns
%     should be treated as numeric vectors
%    {default = []}
%  'user': user name of account with StudyInfo (i.e. ProjInfo, VisitInfo)
%    If not supplied, will use $USER environment variable
%    {default = []}
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
%   StudyInfo: struct array with info for each "study" (scan session)
%   RootDirs: struct containing locations of data directories
%   ProjInfo: struct containing parameters for this ProjID
%
% Created:  03/18/11 by Don Hagler
% Last Mod: 04/03/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'VisitIDs',[],[],...
  'SubjIDs',[],[],...
  'ignore_VisitInfo_flag',false,[false true],...
  'checkflag',true,[false true],...
  'qcflag',false,[false true],...
  'modality',[],{'MRI','MEG','PET'},...
  'required_containers',[],[],...
  'required_rootdirs',[],[],...
  'ico',7,[1:7],...
  'numvec_tags',[],[],...
  'user',[],[],...
...
  'QC_raw',false,[false true],...
  'QC_recon',false,[false true],...
  'QC_dv',false,[false true],...
  'QC_DTI',false,[false true],...
  'QC_BOLD',false,[false true],...
  'QC_PET',false,[false true],...
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
[ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID,parms.numvec_tags,parms.user);

% load visit info file if exists
fname_visit = sprintf('%s/ProjInfo/%s/%s_VisitInfo.csv',...
  RootDirs.home,ProjID,ProjID);
if exist(fname_visit,'file') && ~parms.ignore_VisitInfo_flag
  fprintf('%s: creating StudyInfo from %s...\n',mfilename,fname_visit);
  StudyInfo = MMIL_Read_StudyInfo(fname_visit);
end

% find containers for each visit, or create studyinfo based on raw containers
args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs,args{:});

% if empty, MRI project has not been unpacked, or no orig data
if isempty(StudyInfo), return; end;

% get manual QC information if available
if ~isempty(ProjInfo) &&...
    (parms.QC_raw || parms.QC_recon || parms.QC_dv ||...
     parms.QC_DTI  || parms.QC_BOLD || parms.QC_PET)
  args = MMIL_Args(parms,'MMIL_Get_QCInfo');
  StudyInfo = MMIL_Get_QCInfo(StudyInfo,RootDirs,ProjInfo,args{:});
  if parms.qcflag
    StudyInfo = StudyInfo([StudyInfo.QC]==1);
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
