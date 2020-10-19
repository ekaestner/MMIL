function SeriesInfo = mmil_modify_serinfo_rfrxcoiltype(SeriesInfo)
%function SeriesInfo = mmil_modify_serinfo_rfrxcoiltype(SeriesInfo)
%
% Purpose: determine type of receive coil from dicom info
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Created:  11/15/12 by Don Hagler
% Last Mod: 11/15/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

PA_types = {'8HRBRAIN','SENSE-HEAD','HighResolutionHead',...
            'HeadMatrix','8_Channel_Head','NeckMatrix','32Ch_Head'};
HC_types = {'HEAD','CP_HeadArray','CP_SpineArray'};
BODY_tpyes = {'Body'};

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  ReceiveCoilName = ...
    upper(mmil_getfield(SeriesInfo(s),'ReceiveCoilName','UNKNOWN'));

  ind = find(~cellfun(@isempty,regexpi(ReceiveCoilName,PA_types)));
  if ~isempty(ind)
    SeriesInfo(s).rfrxcoiltype = 'PA';
    continue;
  end;

  ind = find(~cellfun(@isempty,regexpi(ReceiveCoilName,HC_types)));
  if ~isempty(ind)
    SeriesInfo(s).rfrxcoiltype = 'HC';
    continue;
  end;

  ind = find(~cellfun(@isempty,regexpi(ReceiveCoilName,BODY_tpyes)));
  if ~isempty(ind)
    SeriesInfo(s).rfrxcoiltype = 'BODY';
    continue;
  end;

  SeriesInfo(s).rfrxcoiltype = 'UNKNOWN';
  fprintf('%s: ReceiveCoil type UNKNOWN for series %d\n',...
    mfilename,SeriesInfo(s).SeriesNumber);
end

