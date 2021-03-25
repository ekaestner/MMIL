function SeriesInfo = mmil_slicetiming_serinfo(SeriesInfo)
%function SeriesInfo = mmil_slicetiming_serinfo(SeriesInfo)
%
% Purpose: modify SeriesInfo by setting slice_tpattern and slice_order
%   for series with 'BOLD' in the SeriesType
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series including SeriesType
%
% Output:
%   SeriesInfo: struct array with slice timing information for BOLD series
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 11/14/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  if strncmp(SeriesInfo(s).SeriesType,'BOLD',4)
    if mmil_Philips(SeriesInfo(s))
      [SeriesInfo(s).slice_tpattern,SeriesInfo(s).slice_order] = ...
        mmil_slicetiming_Philips(SeriesInfo(s));
    elseif mmil_GE(SeriesInfo(s))
      [SeriesInfo(s).slice_tpattern,SeriesInfo(s).slice_order] = ...
        mmil_slicetiming_GE(SeriesInfo(s));
    elseif mmil_Siemens(SeriesInfo(s))
      [SeriesInfo(s).slice_tpattern,SeriesInfo(s).slice_order] = ...
        mmil_slicetiming_Siemens(SeriesInfo(s));
    end;
  end;
end;

