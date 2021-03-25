function mismatch_flag = mmil_check_series_nimgs(serinfo,i_series)
%function mismatch_flag = mmil_check_series_nimgs(serinfo,i_series)
%
% Required Input:
%   serinfo: struct containing information for i_series
%   i_series: series index
%
% Output:
%   mismatch_flag: [0|1] whether there is a mismatch between
%     serinfo.ImagesInAcquisition and serinfo.nimgs
%
% Created:  09/12/12 by Don Hagler
% Last Mod: 09/12/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

mismatch_flag = 0;

if ~isfield(serinfo,'ImagesInAcquisition') || serinfo.ImagesInAcquisition==0
  fprintf('%s: WARNING: ImagesInAcquisition field missing for series %d\n',...
    mfilename,i_series);
elseif serinfo.ImagesInAcquisition~=serinfo.nimgs
  fprintf('%s: WARNING: ImagesInAcquisition (%d) does not match number of images (%d) for series %d... skipping\n',...
    mfilename,serinfo.ImagesInAcquisition,serinfo.nimgs,i_series);
  mismatch_flag = 1;
end;

return;

