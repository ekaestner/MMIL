function SeriesInfo = mmil_classify_by_code(SeriesInfo)
%function SeriesInfo = mmil_classify_by_rules(SeriesInfo)
%
% Purpose: classify series of dicom files based on codified rules
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series returned from mmil_classify_dicoms
%
% Output:
%   SeriesInfo: will contain updated SeriesType field plus
%     additional information set in code
%
% Note: This function may be replaced with a customized version
%   for specialized classification
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 12/03/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if 0

  if mmil_Philips(SeriesInfo)
    SeriesInfo = mmil_classify_Philips(SeriesInfo);
  elseif mmil_GE(SeriesInfo)
    SeriesInfo = mmil_classify_GE(SeriesInfo);
  elseif mmil_Siemens(SeriesInfo)
    SeriesInfo = mmil_classify_Siemens(SeriesInfo);
  end;

end;


