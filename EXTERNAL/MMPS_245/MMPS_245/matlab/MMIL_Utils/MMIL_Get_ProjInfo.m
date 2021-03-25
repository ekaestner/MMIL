function [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID,numvec_tags,user)
%function [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID,[numvec_tags],[user])
%
% Purpose: Load information from ProjInfo spreadsheet file
%   named /home/<user>/ProjInfo/MMIL_ProjInfo.csv
%
% Optional Input:
%   ProjID: project identifier string
%
% Optional Input:
%   numvec_tags: cell array of column headers specifying which columns
%     should be treated as numeric vectors
%   user: user name of account with StudyInfo (i.e. ProjInfo, VisitInfo)
%    If not supplied, will use $USER environment variable
%    {default = []}
%
% Output:
%   ProjInfo: struct containing parameter settings for ProjID
%   RootDirs: struct containing locations of data directories
%
% Created:   03/18/11 by Don Hagler
% Last Mod:  03/12/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if isempty(ProjID), error('ProjID is empty'); end;
if ~ischar(ProjID), error('ProjID must be a string.'); end
if ~exist('numvec_tags','var'), numvec_tags = []; end;
if ~exist('user','var') || isempty(user), user = getenv('USER'); end;

fname = ['/home/' user '/ProjInfo/MMIL_ProjInfo.csv'];
ProjInfo = MMIL_Read_ProjInfo(fname,numvec_tags);

if ~isempty(ProjID)
  ProjIDs = {ProjInfo.ProjID};
  [ProjIDs,ind] = intersect(ProjIDs,ProjID);
  if isempty(ind)
    error('specified ProjID %s not found in ProjInfo',ProjID);
  else
    ProjInfo = ProjInfo(ind);
  end
end;

RootDirs = MMIL_RootDirs_from_ProjInfo(ProjInfo);

if ~strcmp(user,getenv('USER'))
  RootDirs.home = sprintf('/home/%s',user);
end;
