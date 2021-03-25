function SeriesInfo = mmil_serinfo_Philips(SeriesInfo)
%function SeriesInfo = mmil_serinfo_Philips(SeriesInfo)
%
% Purpose: gather information from dicom header for manufacturer Philips
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Output:
%   SeriesInfo: struct containing additional information
%     extracted from Philips-specific private tags
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 05/01/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  SeriesInfo(s).readoutFOV = mmil_getfield(SeriesInfo(s).info,...
    'ReconstructionDiameter');
  if isfield(SeriesInfo(s).info,'PercentPhaseFieldOfView')
    SeriesInfo(s).phaseFOV = ...
      SeriesInfo(s).readoutFOV*SeriesInfo(s).info.PercentPhaseFieldOfView/100;
  end;
  if isfield(SeriesInfo(s).info,'Private_2001_1018')
    SeriesInfo(s).nslices = SeriesInfo(s).info.Private_2001_1018(1);
  end;
  tmp = mmil_getfield(SeriesInfo(s).info,'Private_2005_140f',[]);
  if ~isa(tmp,'struct')
    SeriesInfo(s).SequenceName = ...
      upper(char(mmil_getfield(SeriesInfo(s).info,'Private_2001_1020')));
  elseif ~isempty(tmp) && isa(tmp,'struct') &&...
          isfield(tmp,'Item_1') && ~isempty(tmp.Item_1)
    SeriesInfo(s).SequenceName = ...
      upper(char(mmil_getfield(tmp.Item_1,'PulseSequenceName')));
  end;
end;

