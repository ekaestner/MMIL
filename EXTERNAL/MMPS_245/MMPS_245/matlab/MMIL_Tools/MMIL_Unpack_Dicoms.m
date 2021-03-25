function [ContainerOutDir,errcode] = MMIL_Unpack_Dicoms(...
  SourceRootDir,SourceDir,ContainerRootDir,VisitID,linkflag,forceflag)
% function [ContainerOutDir,errcode] = MMIL_Unpack_Dicoms(...
%   SourceRootDir,SourceDir,ContainerRootDir,VisitID,linkflag,[forceflag])
%
% Purpose: Unpack DICOM directory into MMIL data container
%
% Required Input:
%   SourceRootDir: root directory for original dicom files
%   SourceDir: subdirectory within OrigRootDir containing dicom files
%   ContainerRootDir: root directory for output raw dicom files
%   VisitID: unique identifier string for scan session
%
% Optional Input:
%   linkflag: whether to create symbolic links (otherwise make copies)
%     {default = 1}
%   forceflag: whether to overwrite existing output
%     {default = 0}
%
% Created:  12/18/08 by Don Hagler
% Last Mod: 11/14/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,4), return; end;
errcode = 0;
if ~exist('linkflag','var'), linkflag = true; end
if ~exist('forceflag','var'), forceflag = false; end
batch_limit = 250; % larger value unpacks faster, too high results in error
                   %   because of unix command string length limit

max_attempts = 10;
license_pause = 600;

ContainerInfo = [];
StudyInstanceUIDlist = {};
SeriesInstanceUIDlist = {};
SeriesInfo = [];
ContainerOutDir = [];
ContainerOutPath = [];

indir = sprintf('%s/%s',SourceRootDir,SourceDir);

% get info from first valid dicom file
d = recursive_dir(char(indir));
if isempty(d)
  fprintf('%s: ERROR: no files found in %s\n',mfilename,indir);
  errcode = 1;
  return;
end;
for dindx = 1:length(d)
  fname = char(d(dindx));
  if ~mmil_isdicomfile(fname), continue; end;
  for i=1:max_attempts
    try
      tmpinfo = dicominfo(fname); 
    catch ERR
      tmpinfo = [];
      if ~isempty(regexp(ERR.message,'License Manager Error'))
        if i==max_attempts
          rethrow(ERR);
        end;
        fprintf('%s: WARNING: Matlab License is unavailable, pausing for %d seconds...',...
          mfilename,license_pause);
        pause(license_pause);
      else
        break;
      end;
    end
    if ~isempty(tmpinfo), break; end;
  end;
  if ~isempty(tmpinfo)
    if ~isfield(tmpinfo,'Modality'), tmpinfo.Modality = ''; end;
    if ~isempty(findstr('DERIVED',mmil_getfield(tmpinfo,'ImageType',[])))
      continue;
    end;
    if isempty(findstr('MR',mmil_getfield(tmpinfo,'Modality','')))
      continue;
    end;
    ContainerUID = '1'; % Hardcode for now, need system call to obtain unique container ID
    ContainerInfo.ContainerType = 'MRIRAW';
    ContainerInfo.ContainerUID = ContainerUID;
    ContainerInfo.ContainerCreationDate = datestr(now);
    ContainerInfo.VisitID = VisitID;
    ContainerInfo.StudyDate = tmpinfo.StudyDate;
    ContainerInfo.StudyTime = tmpinfo.StudyTime;
    ContainerOutDir = sprintf('MRIRAW_%s_%s.%s_%s',...
      VisitID,ContainerInfo.StudyDate,ContainerInfo.StudyTime,ContainerUID);
    ContainerOutPath = sprintf('%s/%s',ContainerRootDir,ContainerOutDir);
    ContainerInfo.SourceRootDir = SourceRootDir;
    ContainerInfo.SourceDir = SourceDir;
    break;
  end
end

if isempty(ContainerOutPath)
  fprintf('%s: ERROR: no valid dicom files found in %s\n',mfilename,indir);
  errcode = 1;
  return;
end;  

fname = sprintf('%s/ContainerInfo.mat',ContainerOutPath);
fname_dir = dir(fname);
if exist(fname,'file') && fname_dir.bytes == 0
  fprintf('%s: WARNING: removing empty ContainerInfo file %s\n',mfilename,fname);
  delete(fname);
end
if ~exist(fname,'file') || forceflag
  [SeriesInfo,errcode,msg] = mmil_unpack_dicoms(indir,...
    'outdir',ContainerOutPath,'batch_limit',batch_limit,...
    'linkflag',linkflag,'forceflag',forceflag);
  if errcode
    fprintf('%s: WARNING: dicom unpacking failed:\n%s\n',...
      mfilename,msg);
    return;
  end;
  ContainerInfo.StudyInstanceUIDlist = unique({SeriesInfo.StudyInstanceUID});
  ContainerInfo.SeriesInstanceUIDlist = unique({SeriesInfo.SeriesInstanceUID});
  ContainerInfo.SeriesInfo = SeriesInfo;
  MMIL_Save_ContainerInfo(ContainerOutPath,ContainerInfo);
end;

