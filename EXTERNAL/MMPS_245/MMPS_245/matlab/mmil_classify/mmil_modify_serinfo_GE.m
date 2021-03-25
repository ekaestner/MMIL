function SeriesInfo = mmil_serinfo_GE(SeriesInfo)
%function SeriesInfo = mmil_serinfo_GE(SeriesInfo)
%
% Purpose: gather information from dicom headers for manufacturer GE
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Output:
%   SeriesInfo: struct containing additional information
%     extracted from GE-specific private tags
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  SeriesInfo(s).SequenceLabel = ...
    mmil_rowvec(char(mmil_getfield(SeriesInfo(s).info,'Private_0019_109c')));
  SeriesInfo(s).SequenceName = ...
    mmil_rowvec(char(mmil_getfield(SeriesInfo(s).info,'Private_0019_109e')));
  if isfield(SeriesInfo(s).info,'ImagesInAcquisition')
    SeriesInfo(s).ImagesInAcquisition = SeriesInfo(s).info.ImagesInAcquisition;
  elseif isfield(SeriesInfo(s).info,'Private_0025_1007')
    SeriesInfo(s).ImagesInAcquisition = SeriesInfo(s).info.Private_0025_1007;
  else
    fprintf('%s: WARNING: DICOM tag for # of images in series not found\n',...
      mfilename);
  end
  SeriesInfo(s).readoutFOV = ...
    mmil_getfield(SeriesInfo(s).info,'ReconstructionDiameter');
  if isfield(SeriesInfo(s).info,'PercentPhaseFieldOfView')
    SeriesInfo(s).phaseFOV = ...
      SeriesInfo(s).readoutFOV*SeriesInfo(s).info.PercentPhaseFieldOfView/100;
  end;
  SeriesInfo(s).nslices = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0021_104f',0);
  SeriesInfo(s).ndiffdirs = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0019_10e0',0);
  SeriesInfo(s).pepolar = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0019_10b3',-1);
  SeriesInfo(s).tensor_fnum = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0019_10b6');
  SeriesInfo(s).EchoSpacing = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0043_102c',-1);
  SeriesInfo(s).AcquisitionRows = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0027_1060');
  SeriesInfo(s).AcquisitionColumns = ...
    mmil_getfield(SeriesInfo(s).info,'Private_0027_1061');
  SeriesInfo(s).magnitude_flag = ...
    ~isempty(findstr('ORIGINAL\PRIMARY\M\ND',SeriesInfo(s).ImageType)) | ...
    ~isempty(findstr('ORIGINAL\PRIMARY\OTHER',SeriesInfo(s).ImageType));
end;

