function [ContainerOutDir,errcode] = MMIL_Convert_NonDicom_MRIRAW_to_MRIPROC(...
  ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag)
%function [ContainerOutDir,errcode] = MMIL_Convert_NonDicom_MRIRAW_to_MRIPROC(...
% ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag)
%
% Required Input:
%  ContainerRootInDir
%  ContainerInDir
%  ContainerRootOutDir
%
% Optional Input:
%   forceflag: [0|1] whether to overwrite existing output
%
% supported scan types:
%   MPR, FLASHhi, FLASHlo
%
% Created:  05/28/09 by Alain Koyama
% Last Mod: 02/12/13 by Don Hagler
%

errcode = 0;
ContainerOutDir = [];

if ~exist('forceflag','var'), forceflag = false; end;

fprintf('%s(''%s'',''%s'',''%s'',%d)\n',mfilename,ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag);

ContainerInPath = sprintf('%s/%s',ContainerRootInDir,ContainerInDir);
try
  fname = sprintf('%s/ContainerInfo.mat',ContainerInPath);
  load(fname);
catch
  fprintf('%s: ERROR: cannot read file %s\n',mfilename,fname);
  errcode = 1;
  return
end
ContainerInInfo = ContainerInfo;
clear ContainerInfo;

ContainerOutInfo.SourceDir = ContainerInDir;
ContainerOutUID = '1'; % Hardcode for now, need system call to obtain unique container ID
ContainerOutInfo.ContainerType = 'MRIPROC';
ContainerOutInfo.ContainerUID = ContainerOutUID;
ContainerOutInfo.VisitID = ContainerInInfo.VisitID;
ContainerOutInfo.StudyDate = ContainerInInfo.StudyDate;
ContainerOutInfo.StudyTime = ContainerInInfo.StudyTime;
ContainerOutInfo.ContainerCreationDate = datestr(now);
ContainerOutInfo.MagneticFieldStrength = ContainerInInfo.MagneticFieldStrength;

ContainerOutDir = sprintf('MRIPROC_%s_%s.%s_%s',ContainerOutInfo.VisitID,...
  ContainerOutInfo.StudyDate,ContainerOutInfo.StudyTime,ContainerOutUID);
ContainerOutPath = sprintf('%s/%s',ContainerRootOutDir,ContainerOutDir);

if ~exist(ContainerOutPath,'file'), mkdir(ContainerOutPath); end

MPR_cntr = 0;
MEDIChi_cntr = 0;
MEDIClo_cntr = 0;
XetaT2_cntr = 0;
FLASHhi_cntr = 0;
FLASHlo_cntr = 0;
B1HC_cntr = 0;
B1BC_cntr = 0;
GEB1CAL_cntr = 0;

SeriesInfo = ContainerInInfo.SeriesInfo;

ScanInfo = [];

for s=1:length(SeriesInfo)
  serinfo = SeriesInfo(s);

  indir = sprintf('%s/%s',...
    ContainerInPath,serinfo.SeriesDirPath);
  fnames = dir(indir);
  fnames = fnames(cellfun('isempty', regexp({fnames.name},'^\.')));
  if isempty(fnames), continue; end; % if no files in this subdir

  n = regexp({fnames.name},'BRIK$'); n = find(~cellfun('isempty', n)==1);
  if ~isempty(n) % if BRIK so load .BRIK not .HEAD into mri_convert
    fnames = {fnames(n).name};
  else % sequentially named series of non-DICOM's
    fnames = {fnames.name};
  end
  cmd = sprintf('mri_convert %s/%s %s',...
    indir,fnames{1},fname);
  [status, result] = system(cmd);
  if status
    fprintf('%s: ERROR: %s\n',mfilename,result);
    return;
  end

  switch serinfo.SeriesType
    case 'MPR'
      MPR_cntr = MPR_cntr+1;
      fname = sprintf('%s/MPR%d.mgz',ContainerOutPath,MPR_cntr);
      ScanInfo.MPR(MPR_cntr) = struct(...
        'SeriesIndex',s,...
        'TE',1,'gradwarpinfo',serinfo.gradwarpinfo); 
    case 'FLASHhi'
      FLASHhi_cntr = FLASHhi_cntr+1;
      fname = sprintf('%s/FLASHhi%d.mgz',ContainerOutPath,FLASHhi_cntr);
      ScanInfo.FLASHhi(FLASHhi_cntr) = struct(...
        'SeriesIndex',s,...
        'TE',1,'gradwarpinfo',serinfo.gradwarpinfo); 
    case 'FLASHlo'
      FLASHlo_cntr = FLASHlo_cntr+1;
      fname = sprintf('%s/FLASHlo%d.mgz',ContainerOutPath,FLASHlo_cntr);
      ScanInfo.FLASHlo(FLASHlo_cntr) = struct(...
        'SeriesIndex',s,...
        'TE',1,'gradwarpinfo',serinfo.gradwarpinfo); 
  end
end


ContainerOutInfo.MPR_cntr = MPR_cntr;
ContainerOutInfo.FLASHhi_cntr = FLASHhi_cntr;
ContainerOutInfo.FLASHlo_cntr = FLASHlo_cntr;
ContainerOutInfo.XetaT2_cntr = XetaT2_cntr;
ContainerOutInfo.MEDIChi_cntr = MEDIChi_cntr;
ContainerOutInfo.MEDIClo_cntr = MEDIClo_cntr;
ContainerOutInfo.B1HC_cntr = B1HC_cntr;
ContainerOutInfo.B1BC_cntr = B1BC_cntr;
ContainerOutInfo.GEB1CAL_cntr = GEB1CAL_cntr;
ContainerOutInfo.ScanInfo = ScanInfo;
ContainerOutInfo.ProcessingCompleted = now;

errcode = MMIL_Save_ContainerInfo(ContainerOutPath,ContainerOutInfo);

% save SeriesIndex to ScanInfo.csv file
fname = sprintf('%s/ScanInfo.csv',ContainerOutPath);
fid = fopen(fname,'wt');
if fid<0
  fprintf('%s: WARNING: failed to open %s for writing\n',mfilename,fname);
  errcode = 1;
  return;
end;
fprintf(fid,'"ScanType","ScanIndex","SeriesIndex"\n');
ScanTypes = fieldnames(ScanInfo);
for t=1:length(ScanTypes)
  type = ScanTypes{t};
  tmpInfo = ScanInfo.(type);
  for s=1:length(tmpInfo)
    fprintf(fid,'"%s",%d,%d\n',type,s,tmpInfo(s).SeriesIndex);
  end;
end;
fclose(fid);

