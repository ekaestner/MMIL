function N = ReadPETdcmfiles(fnames)
%function N = ReadPETdcmfiles(fnames)
%
% Created:  11/02/06 by Anders Dale
% Last Mod: 11/15/12 by Don Hagler
%

wtmp = warning('query','all');

N = {};
for i = 1:length(fnames)
  fname = fnames{i};
  warning('off','all');
  info = dicominfo(fname);
  warning(wtmp);
  Modality = info.Modality;
  ImageType = info.ImageType;
  if strcmp(Modality,'PT') & ~isempty(findstr('PRIMARY',ImageType)) % Check if PET image & PRIMARY
    if isempty(N)
      NumberOfSlices = double(info.NumberOfSlices);
      Rows = info.Rows;
      Columns = info.Columns;
      if isfield(info,'NumberOfTimeSlices')
          NumberOfTimeSlices = double(info.NumberOfTimeSlices);
          if NumberOfTimeSlices == 0
              NumberOfTimeSlices = floor(length(fnames)/NumberOfSlices);
              fprintf('%s: WARNING: DICOM field NumberOfTimeSlices = 0, assuming NumberOfTimeSlices = %d\n',...
                mfilename,NumberOfTimeSlices);
          elseif (NumberOfTimeSlices*NumberOfSlices > length(fnames))
              fprintf('%s: WARNING: NumberOfSlices = %d * NumberOfTimeSlices = %d > length(fnames) = %d\n',...
                mfilename,NumberOfSlices,NumberOfTimeSlices,length(fnames));
              NumberOfSlices = floor(length(fnames)/NumberOfTimeSlices);
              fprintf('    assuming NumberOfSlices = %d\n',NumberOfSlices);
          end
      else
          NumberOfTimeSlices = floor(length(fnames)/NumberOfSlices);
          fprintf('%s: WARNING: DICOM field NumberOfTimeSlices missing, assuming NumberOfTimeSlices = %d\n',...
            NumberOfTimeSlices);
      end
      PatientID = info.PatientID;
      StudyDate = info.StudyDate;
      StudyTime = info.StudyTime;
      SeriesDate = info.SeriesDate;
      SeriesTime = info.SeriesTime;
      SeriesType = info.SeriesType;
      PixelSpacing = info.PixelSpacing;
      SliceThickness = info.SliceThickness;
%      if isfield(info,'FieldOfViewDimensions')
%        FieldOfViewDimensions = info.FieldOfViewDimensions;
%      else
%        fprintf('WARNING: DICOM field FieldOfViewDimensions missing\n');
%      end
      SeriesNumber = info.SeriesNumber;
      SeriesDescription = info.SeriesDescription;
      SeriesInstanceUID = info.SeriesInstanceUID;
      for j = 1:NumberOfTimeSlices
        N{j}.vol = zeros(Rows,Columns,NumberOfSlices);
      end
%      SliceLocations = zeros(NumberOfTimeSlices,NumberOfSlices);
      ImagePositions = zeros(NumberOfTimeSlices,NumberOfSlices,3);
      ImageOrientations = zeros(NumberOfTimeSlices,NumberOfSlices,6);
    end
    if strcmp(info.SeriesInstanceUID,SeriesInstanceUID) % Make sure images are from same series
      InstanceNumber = double(info.InstanceNumber);
      timeslicenum = 1+floor((InstanceNumber-1)/NumberOfSlices);
      slicenum = 1+mod(InstanceNumber-1,NumberOfSlices);
%      fprintf('imnum=%d: InstanceNumber=%d, timeslicenum=%d, slicenum=%d\n',i,InstanceNumber,timeslicenum,slicenum);
%      SliceLocations(timeslicenum,slicenum) = info.SliceLocation;
      ImagePositions(timeslicenum,slicenum,:) = info.ImagePositionPatient;
      ImageOrientations(timeslicenum,slicenum,:) = info.ImageOrientationPatient;
      warning('off','all');
      im = dicomread(fname);
      warning(wtmp);
      RescaleOffset = mmil_getfield(info,'RescaleOffset',0.0);
      RescaleSlope = mmil_getfield(info,'RescaleSlope',1.0);
      im = RescaleSlope*im+RescaleOffset;
      N{timeslicenum}.vol(:,:,slicenum) = im;
      if isfield(info,'ActualFrameDuration')
        framedur{timeslicenum} = double(info.ActualFrameDuration)/1000;
      else
        framedur{timeslicenum} = 1.0;
        if i==1
          fprintf('%s: WARNING: DICOM field ActualFrameDuration missing, assuming ActualFrameDuration = %0.1f\n',...
            mfilename,framedur{timeslicenum});
        end;
      end;
    end
  else
    fprintf('%s: WARNING: invalid PET DICOM file %s: Modality = %s, ImageType = %s\n',...
      mfilename,fname,Modality,ImageType);
  end
end

if length(N)>0
  slsp = squeeze(ImagePositions(1,2,:)-ImagePositions(1,1,:));
  M = eye(4);
  M(1:3,1:3) = [squeeze(ImageOrientations(1,1,4:6))*PixelSpacing(2) squeeze(ImageOrientations(1,1,1:3))*PixelSpacing(1) slsp];
  %M(1:3,4) = squeeze(ImagePositions(1,1,:)) - M(1:3,:)*[1 1 1 1]';
  M(1:3,4) = -M(1:3,:)*[size(N{1}.vol)/2 1]'; % Set origin of coord system to center of scanned volume
  M = M_LPH_TO_RAS*M;

  StudyDateNum = datenum(sprintf('%s %s',SeriesDate,SeriesTime),'yyyymmdd HHMMSS');
  for i = 1:NumberOfTimeSlices
    N{i}.M = M;
    N{i}.filetype = 'DICOM';
    N{i}.PatientID = PatientID;
    N{i}.StudyDateNum = StudyDateNum;
    N{i}.StudyDateStr = datestr(N{i}.StudyDateNum);
    N{i}.SeriesNumber = SeriesNumber;
    N{i}.SeriesDescription = SeriesDescription;
    N{i}.timing.tspace = framedur{i};
  end
else
% fprintf('WARNING: No PET volumes read\n');
end

