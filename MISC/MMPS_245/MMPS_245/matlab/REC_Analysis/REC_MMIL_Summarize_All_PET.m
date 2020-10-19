function REC_MMIL_Summarize_All_PET(ProjID,normflag)
%function REC_MMIL_Summarize_All_PET(ProjID,normflag)
%
% Required Input:
%   ProjID: project ID string
%
% Optional Input:
%   normflag: normalize to pons (default = 1)
%     (0 = absolute values, 1 = norm to pons, 2 = both)
%
% Last Mod:  01/28/11 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
if isempty(ProjID), error('empty ProjID'); end;
if ~exist('normflag','var'), normflag = 1; end

REC_MMIL_Summarize_PET_aseg(ProjID,normflag);
REC_MMIL_Summarize_PET_aparc(ProjID,normflag);
REC_MMIL_Summarize_CortThick_forPET(ProjID);
REC_MMIL_Summarize_Volume_forPET(ProjID);
