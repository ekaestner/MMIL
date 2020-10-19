function SeriesInfo = mmil_init_serinfo(SeriesInfo)
%function SeriesInfo = mmil_init_serinfo(SeriesInfo)
%
% Purpose: intiialize SeriesInfo by gathering information from dicom header
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series returned from mmil_unpack_dicoms
%
% Output:
%   SeriesInfo: struct array with information extracted from dicom headers
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

default_values = [];
default_values.ImageType = 'UNDEFINED';
default_values.Modality = 'UNDEFINED';
default_values.StudyInstanceUID = [];
default_values.StudyDate = 0;
default_values.StudyTime = 0;
default_values.SeriesInstanceUID = [];
default_values.SeriesNumber = 0;
default_values.SeriesDate = 0;
default_values.SeriesTime = 0;
default_values.SeriesType = 'UNKNOWN';
default_values.SeriesDescription = '';
default_values.SequenceName = '';
default_values.SequenceLabel = '';
default_values.SequenceVariant = '';
default_values.Manufacturer = 'UNDEFINED';
default_values.ManufacturersModelName = 'UNDEFINED';
default_values.DeviceSerialNumber = '';
default_values.MagneticFieldStrength = '';
default_values.PatientID = '';
default_values.PatientName = '';
default_values.PatientBirthDate = '';
default_values.PatientSex = '';
default_values.PatientAge = '';  
default_values.ScanningSequence = '';
default_values.ScanOptions = '';
default_values.ImagesInAcquisition = 0;
default_values.AcquisitionRows = [];
default_values.AcquisitionColumns = [];
default_values.Rows = 0;
default_values.Columns = 0;
default_values.nslices = 0;
default_values.nc = [];
default_values.nr = [];
default_values.PixelSpacing = [];
default_values.XPixelSpacing = [];
default_values.YPixelSpacing = [];
default_values.SliceThickness = 1;
default_values.NumberOfPhaseEncodingSteps = 0;
default_values.InPlanePhaseEncodingDirection = 'UNKNOWN';
default_values.InPlaneRot = [];
default_values.readoutFOV = [];
default_values.phaseFOV = [];
default_values.FlipAngle = NaN;
default_values.RepetitionTime = NaN;
default_values.InversionTime = NaN;
default_values.EchoSpacing = NaN;
default_values.MRAcquisitionType = 'UNDEFINED';
default_values.ReceiveCoilName = 'UNDEFINED';
default_values.PixelBandwidth = [];
default_values.NumberOfAverages = 1;
default_values.diffdirs = [];
default_values.ndiffdirs = 0;
default_values.nb0 = 0;
default_values.bval = 0;
default_values.pepolar = -1;
default_values.tensor_fnum = -1;
default_values.SliceLocation = [];
default_values.gradwarpinfo = struct;
default_values.rfrxcoiltype = '';
default_values.slice_tpattern = [];
default_values.slice_order = [];
default_values.magnitude_flag = [];
default_values.pe_rev_flag = []; % for Siemens DTI and BOLD to determine pepolar
default_values.pe_for_flag = []; % for Siemens DTI and BOLD to determine pepolar
default_values.pe_extra_val = []; % for Siemens DTI and BOLD to determine pepolar

tags = fieldnames(default_values);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:length(SeriesInfo)
  SeriesInfo(s).nimgs = length(SeriesInfo(s).FileNames);
  fname = SeriesInfo(s).FileNames{1};
  if ~exist(fname,'file')
    fprintf('%s: WARNING: file %s not found, skipping series %d...\n',...
      mfilename,fname,s);
    SeriesInfo(s).SeriesType = 'MISSING';
    SeriesInfo(s).ignore = 1;
    continue;
  end;
  SeriesInfo(s).info = dicominfo(fname);
  SeriesInfo(s).ignore = 0;
  if isempty(SeriesInfo(s).info)
    SeriesInfo(s).ignore = 1;
  else
    % set values in SeriesInfo(s) using SeriesInfo(s).info or default_values
    for t=1:length(tags)
      SeriesInfo(s).(tags{t}) = ...
        mmil_getfield(SeriesInfo(s).info,tags{t},default_values.(tags{t}));
    end;
  end;
  % identify invalid dicom file series
  if ~SeriesInfo(s).ignore
    SeriesInfo(s).ignore = check_series(SeriesInfo(s));
  end;
  if SeriesInfo(s).ignore, SeriesInfo(s).SeriesType = 'JUNK'; end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ignore = check_series(serinfo)
  ignore = serinfo.ignore;
  fname = serinfo.FileNames{1};
  % exclude other invalid series
  if ~isempty(findstr('DERIVED',serinfo.ImageType))
    fprintf('%s: ignoring DERIVED image %s\n',...
      mfilename,fname);
    ignore = 1;
  elseif isempty(findstr('MR',serinfo.Modality))
    fprintf('%s: ignoring file %s (Modality %s != MR)\n',...
      mfilename,fname,serinfo.Modality);
    ignore = 1;
  elseif isempty(serinfo.InstanceNumber)
    fprintf('%s: ignoring file %s with missing InstanceNumber\n',...
      mfilename,fname);
    ignore = 1;
  elseif isfield(serinfo,'errmsg') && ~isempty(serinfo.errmsg),
    ignore = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

