function MMIL_Long_Mask_T1(fname_aseg,fname_mask,forceflag)
%function MMIL_Long_Mask_T1(fname_aseg,fname_mask,forceflag)
%
% Created:  10/15/09 by Don Hagler
% Last Mod: 10/15/09 by Don Hagler
%

smooth1 = 2;
thresh1 = 0.2;
smooth2 = 10;
thresh2 = 0.1;
smooth3 = 5;
thresh3 = 0.01;

if (~mmil_check_nargs(nargin,2)) return; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag = 0; end;

mmil_dilate_mask([],...
  'fname_in',fname_aseg,'fname_out',fname_mask,...
  'smooth1',smooth1,'thresh1',thresh1,...
  'smooth2',smooth2,'thresh2',thresh2,...
  'smooth3',smooth3,'thresh3',thresh3,...
  'forceflag',forceflag);

