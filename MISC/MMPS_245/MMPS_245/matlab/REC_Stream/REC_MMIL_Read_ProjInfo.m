function ProjInfo = REC_MMIL_Read_ProjInfo(ProjID)
%function ProjInfo = REC_MMIL_Read_ProjInfo(ProjID)
%
% Optional Input:
%   ProjID: project identifier string
%     If not supplied, will retrun struct array with all projects
%
% Created:   03/03/09 by Alain Koyama
% Last Mod:  03/23/11 by Don Hagler
%

if ~exist('ProjID','var'), ProjID = ''; end;

if ~ischar(ProjID),
   error('ProjID must be a string.');
end

fname = ['/home/' getenv('USER') '/ProjInfo/MMIL_ProjInfo.csv'];
if ~exist(fname,'file') % temporary fix until we complete transition away from REC
  fname = ['/home/' getenv('USER') '/ProjInfo/REC_MMIL_ProjInfo.csv'];
end;
ProjInfo = MMIL_Read_ProjInfo(fname);

if ~isempty(ProjID)
  ProjIDs = {ProjInfo.ProjID};
  [ProjIDs,ind] = intersect(ProjIDs,ProjID);
  if isempty(ind)
    error('specified ProjID %s not found in ProjInfo',ProjID);
  else
    ProjInfo = ProjInfo(ind);
  end
end;

