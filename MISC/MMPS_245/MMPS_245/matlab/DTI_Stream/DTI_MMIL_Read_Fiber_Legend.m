function FiberLegend = DTI_MMIL_Read_Fiber_Legend(fname)
%function FiberLegend = DTI_MMIL_Read_Fiber_Legend([fname])
%
% Created:  03/22/10 by Don Hagler
% Last Mod: 01/14/13 by Don Hagler
%

if ~exist('fname','var') || isempty(fname)
  dir_parms = getenv('MMPS_PARMS');
  fname = [dir_parms '/DTI_Fiber/DTI_Fiber_Legend.csv'];
end;

FiberLegend = MMIL_Read_StudyInfo(fname);

