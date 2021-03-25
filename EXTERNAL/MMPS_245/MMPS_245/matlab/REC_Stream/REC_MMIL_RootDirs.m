function [RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID)
%function [RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID)
%
% Required Parameters:
%   ProjID: project ID string
%
% Output:
%   RootDirs: struct containing locations of data directories
%   ProjInfo: struct containing parameters for this ProjID
%
% Created:   06/04/09 by Don Hagler
% Last Mod:  09/17/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
if ~ischar(ProjID)
   error('ProjID must be a project ID string');
end

ProjInfo = REC_MMIL_Read_ProjInfo(ProjID);
RootDirs = MMIL_RootDirs_from_ProjInfo(ProjInfo);

