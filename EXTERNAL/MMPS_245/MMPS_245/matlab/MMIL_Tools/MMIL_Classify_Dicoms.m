function errcode = MMIL_Classify_Dicoms(ContainerPath,RejectSeries,classify_file,forceflag)
%function errcode = MMIL_Classify_Dicoms(ContainerPath,RejectSeries,classify_file,forceflag)
%
% Required Input:
%   ContainerPath: full path name of MRIRAW container
%
% Optional Parameters:
%   RejectSeries: vector containing indeces of manually rejected series
%     {default = []}
%   classify_file: csv file containing classification rules
%     see mmil_read_classify_rules and mmil_classify_by_rules
%     {default = []}
%   forceflag: [0|1] overwrite existing classification
%     {default = 0}
%
% Created:  12/19/08 by Don Hagler
% Last Mod: 02/11/14 by Don Hagler
%

errcode = 0;

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('RejectSeries','var'), RejectSeries = []; end
if ~exist('classify_file','var'), classify_file = ''; end
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

% fields to be added to ContainerInfo from SeriesInfo(1)
container_tags = {...
  'StudyInstanceUID'...
  'Manufacturer'...
  'ManufacturersModelName'...
  'DeviceSerialNumber'...
  'MagneticFieldStrength'...
  'PatientID'...
  'PatientName'...
  'PatientBirthDate'...
  'PatientSex'...
  'PatientAge'...
};

fprintf('%s: classifying dicoms for %s...\n',mfilename,ContainerPath);

[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);

if isfield(ContainerInfo,'ClassificationCompleted')
  fprintf('%s: Container %s already classified on %s\n',...
    mfilename,ContainerPath,datestr(ContainerInfo.ClassificationCompleted));
end;
if ~isfield(ContainerInfo,'ClassificationCompleted') || forceflag
  % read csv file containing classification rules
  if ~isempty(classify_file)
    ContainerInfo.classify_rules = mmil_read_classify_rules(classify_file);
  else
    ContainerInfo.classify_rules = [];
  end;

  % get information from headers and determine series types
  SeriesInfo = ...
    mmil_classify_dicoms(ContainerInfo.SeriesInfo,ContainerInfo.classify_rules);

  % copy values of several fields to ContainerInfo from first non-ignored series
  for s=1:length(SeriesInfo)
    if ~SeriesInfo(s).ignore, break; end;
  end;
  for t=1:length(container_tags)
    ContainerInfo.(container_tags{t}) = SeriesInfo(s).(container_tags{t});
  end;

  % handle any manually rejected series
  if ~isempty(RejectSeries),
    if ~isnumeric(RejectSeries) |...
       ~isempty(setdiff(RejectSeries,[1:length(SeriesInfo)]))
      fprintf('%s: WARNING: invalid entries in RejectSeries. Skipping series-reject step.\n',...
        mfilename);
    else
      for i=1:length(RejectSeries)
        SeriesInfo(RejectSeries(i)).SeriesType = [SeriesInfo(RejectSeries(i)).SeriesType '_REJECT'];
      end;
      fprintf('%s: %d series manually rejected: SeriesIndex = [ ',mfilename,length(RejectSeries));
      fprintf('%d ',RejectSeries); fprintf(']\n');
    end
  end

  % display results
  fprintf('%s: Manufacturer = %s, ManufacturersModelName = %s, DeviceSerialNumber = %s, MagneticFieldStrength = %f\n',...
    mfilename,ContainerInfo.Manufacturer, ContainerInfo.ManufacturersModelName,...
    ContainerInfo.DeviceSerialNumber, ContainerInfo.MagneticFieldStrength);

  for s=1:length(SeriesInfo)
    fprintf('Series %d (%s): %d images\n',s,...
      SeriesInfo(s).SeriesDirPath,SeriesInfo(s).nimgs);

    unwarpflag  = mmil_getfield(SeriesInfo(s).gradwarpinfo,'unwarpflag',-1);
    gwtype = mmil_getfield(SeriesInfo(s).gradwarpinfo,'gwtype',-1);
    ambiguousgwtype = mmil_getfield(SeriesInfo(s).gradwarpinfo,'ambiguousgwtype',0);

    fprintf('  ----> SeriesDescription=%s, SequenceName=%s, SequenceVariant=%s, SequenceLabel=%s, ScanOptions=%s, rfrxcoiltype=%s, unwarpflag=%d, gwtype=%d, ambiguousgwtype=%d, FlipAngle=%d, SeriesType=%s\n',...
      SeriesInfo(s).SeriesDescription,...
      SeriesInfo(s).SequenceName,...
      SeriesInfo(s).SequenceVariant,...
      SeriesInfo(s).SequenceLabel,...
      SeriesInfo(s).ScanOptions,...
      SeriesInfo(s).rfrxcoiltype,...
      unwarpflag,gwtype,ambiguousgwtype,...
      SeriesInfo(s).FlipAngle,...
      SeriesInfo(s).SeriesType);
  end

  % save some info to SeriesInfo.csv file
  fname = sprintf('%s/SeriesInfo.csv',ContainerPath);
  mmil_summarize_serinfo(SeriesInfo,fname,1);

  ContainerInfo.SeriesInfo = SeriesInfo;
  ContainerInfo.ClassificationCompleted = datestr(now);
  errcode = MMIL_Save_ContainerInfo(ContainerPath,ContainerInfo);
end;

