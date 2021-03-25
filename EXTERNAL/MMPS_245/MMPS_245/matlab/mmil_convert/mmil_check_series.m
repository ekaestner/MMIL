function [serinfo,valid_flag,missing_flag] = ...
  mmil_check_series(SeriesInfo,i_series,ContainerInPath,valid_SeriesTypes)
%function [serinfo,valid_flag,missing_flag] = ...
%  mmil_check_series(SeriesInfo,i_series,ContainerInPath,valid_SeriesTypes)
%
% Required Input:
%   SeriesInfo: SeriesInfo struct array from input Container
%   i_series: series index
%   ContainerInPath: full path of input Container
%   valid_SeriesTypes: cell array of series types to be processed
%
% Output:
%   serinfo: struct containing information for i_series
%     this is identical to SeriesInfo(i_series) with the addition of fnames
%   valid_flag: [0|1] whether series is valid for processing
%   missing_flag: [0|1] whether series is missing first or last file
%
% Created:  09/12/12 by Don Hagler
% Last Mod: 02/19/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,4), return; end;

scaninfo = []; valid_flag = 1; missing_flag = 0;

serinfo = SeriesInfo(i_series);
if ~ismember(serinfo.SeriesType,valid_SeriesTypes) || serinfo.ignore
  valid_flag = 0;
  return;
end;

SeriesDescription = mmil_getfield(serinfo,'SeriesDescription','?'); 

if ~exist(serinfo.FileNames{1},'file') || ~exist(serinfo.FileNames{end},'file')
  missing_flag = 1;
  fprintf('%s: skipping %s series %d (%s) with missing files\n',...
    mfilename,serinfo.SeriesType,i_series,SeriesDescription);
  return;
end;

return;

