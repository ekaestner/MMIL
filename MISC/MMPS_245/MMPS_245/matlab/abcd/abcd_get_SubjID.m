function [SubjID,pGUID,EventName,SessionType] = abcd_get_SubjID(fstem,series_info)
%function [SubjID,pGUID,EventName,SessionType] = abcd_get_SubjID(fstem,series_info)
%
% Purpose: extract SubjID from fstem or series_info
%
% Input:
%   fstem: string from json/tgz file stem
%   series_info: struct containing metadata from dicom header
%
% Output:
%   SubjID: ABCD Subject ID
%     underscores, blanks, and minus signs will be replaced with "^"
%     "NDAR_" prefix will be removed
%   pGUID: globally unique ID
%   EventName: name of data acquisition "event"
%     e.g. 'baseline_year_1_arm_1'
%   SessionType: type of data acquisition session
%     e.g. 'SessionA1'
%
% Created:  07/28/16 by Don Hagler
% Last Mod: 12/25/16 by Don Hagler
%

SubjID = []; pGUID = []; EventName = []; SessionType = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('series_info','var'), series_info = []; end;

% make sure fstem is not a full file name
[tmp,fstem] = fileparts(fstem);

% try to get pGUID from fstem
[pGUID,EventName,SessionType] = abcd_get_pGUID(fstem);

if isempty(pGUID) && ~isempty(series_info)
  if isfield(series_info,'PatientName') && ~isempty(series_info.PatientName)
    if isstruct(series_info.PatientName) &&...
      isfield(series_info.PatientName,'FamilyName') &&...
      ~isempty(series_info.PatientName.FamilyName)
      pGUID = series_info.PatientName.FamilyName;
    elseif ischar(series_info.PatientName)
      pGUID = series_info.PatientName;
    end;
  elseif isfield(series_info,'PatientFamilyName') && ~isempty(series_info.PatientFamilyName)
    pGUID = series_info.PatientFamilyName;
  end;
  if isempty(pGUID)
    pGUID = series_info.PatientID;
  end;
  if strcmp(pGUID(1:4),'NDAR') && length(pGUID)>16
    pGUID = pGUID(1:16);
  end;  
end;

SubjID = pGUID;

if isempty(SubjID), return; end;

% remove 'NDAR' prefix if present
SubjID = regexprep(SubjID,'NDAR','');

% remove study date if present
k = regexp(SubjID,'(?<ID>\w+)_\d{8,10}','names');
if ~isempty(k), SubjID = k.ID; end;

% remove underscores, dots, spaces, minus signs, or plus signs
SubjID = regexprep(SubjID,'[_\.\s-+]','');

