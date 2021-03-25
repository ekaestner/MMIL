function SeriesInfo = mmil_sort_serinfo(SeriesInfo)
%function SeriesInfo = mmil_sort_serinfo(SeriesInfo)
%
% Purpose: sort series by StudyDate, StudyTime, SeriesDate, and SeriesTime
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series returned from mmil_serinfo
%
% Output:
%   SeriesInfo: series information in chronological order of acquisition
%
% Created:  11/13/12 by Don Hagler
% Prev Mod: 02/08/13 by Cooper Roddey
% Last Mod: 12/04/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% sort series by StudyDate, StudyTime, SeriesDate, and SeriesTime
table = {};
for s=1:length(SeriesInfo)
  if isempty(SeriesInfo(s).StudyDate) || ...
     (isnumeric(SeriesInfo(s).StudyDate) &&...
       SeriesInfo(s).StudyDate==0)
    table{s,1} = '0';
  else
    table{s,1} = SeriesInfo(s).StudyDate;
  end;
  if isempty(SeriesInfo(s).StudyTime) ||...
      (isnumeric(SeriesInfo(s).StudyTime) &&...
       SeriesInfo(s).StudyTime==0)
    table{s,2} = '0';
  else
    table{s,2} = SeriesInfo(s).StudyTime;
  end;
  if isempty(SeriesInfo(s).SeriesDate) ||...
      (isnumeric(SeriesInfo(s).SeriesDate) &&...
       SeriesInfo(s).SeriesDate==0)
    table{s,3} = '0';
  else
    table{s,3} = SeriesInfo(s).SeriesDate;
  end;
  if isempty(SeriesInfo(s).SeriesTime) ||...
      (isnumeric(SeriesInfo(s).SeriesTime) &&...
       SeriesInfo(s).SeriesTime==0)
    table{s,4} = '0';
  else
    table{s,4} = SeriesInfo(s).SeriesTime;
  end;
end
[sortedtable,rowindx] = sortrows(table,1:size(table,2));
SeriesInfo = SeriesInfo(rowindx);


