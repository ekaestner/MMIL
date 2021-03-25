function MMIL_Import_PET_Exam(ContainerPath,RootOutDir,VisitID,forceflag)
%function MMIL_Import_PET_Exam(ContainerPath,RootOutDir,VisitID,forceflag)
%
% Early Mod: 04/06/10 by Don Hagler
% Last Mod:  03/24/10 by Don Hagler
%

%% todo: does each series need to have its own Container? (DJH)
%% todo: change so can process not just .dcm files if needed (AK)

fprintf('MMIL_Import_PET_Exam(''%s'',''%s'',''%s'',%d)\n',...
  ContainerPath,RootOutDir,VisitID,forceflag);

if ~exist('forceflag','var'), forceflag = 0; end;

statfname = sprintf('%s/.unpackstatus',ContainerPath);
if exist(statfname,'file') & ~forceflag
    fprintf('%s: WARNING: Study already unpacked\n',mfilename);
    return;
end

dirlist = recursive_dir_regexp(ContainerPath,'\.dcm$');
for i=1:length(dirlist)
  dirlist{i} = fileparts(dirlist{i});
end
dirlist = unique(dirlist);

for h=1:length(dirlist)
  indir = dirlist{h};
  dlist = dir(sprintf('%s/*.dcm',indir));
  if ~isempty(dlist)
    for i = 1:length(dlist)
      fnames{i} = sprintf('%s/%s',indir,dlist(i).name);
    end
    vols = ReadPETdcmfiles(fnames);
  else
    fprintf('%s: WARNING: No PET files found in dir %s\n',mfilename,indir);
    return;
  end
  if isempty(vols)
    fprintf('%s: WARNING: no PET data read in dir %s\n',mfilename,indir);
    return;
  end

  ContainerUID = '1'; % Hardcode for now, need system call to obtain unique container ID
  ContainerInfo.ContainerType = 'PETRAW';
  ContainerInfo.ContainerUID = ContainerUID;
  ContainerInfo.ContainerCreationDate = datestr(now);
  ContainerInfo.VisitID = VisitID;
  if ~isfield(vols{1},'PatientID')
    fprintf('%s: WARNING: PatientID field missing for VisitID %s (possible anonymization problem)\n',...
      mfilename,VisitID);
  elseif ~strcmp(VisitID,vols{1}.PatientID)
    fprintf('%s: WARNING: VisitID = %s does not match PatientID = %s (possible anonymization problem)\n',...
      mfilename,VisitID,vols{1}.PatientID);
  end
  StudyDateNum = vols{1}.StudyDateNum;
  ContainerInfo.StudyDateNum = StudyDateNum;
  ContainerInfo.StudyDate = datestr(StudyDateNum,'yyyymmdd'); % Convert to DICOM date format
  ContainerInfo.StudyTime = datestr(StudyDateNum,'HHMMSS'); % Convert to DICOM time format
  ContainerInfo.filetype = vols{1}.filetype;
  ContainerInfo.SeriesNumber = vols{1}.SeriesNumber;
  ContainerInfo.SeriesDescription = vols{1}.SeriesDescription;
  ContainerInfo.SourceDir = indir;
  ContainerInfo.PIBflag = ~isempty(regexp(indir,'PIB'));
  ContainerOutDir = sprintf('PETRAW_%s_%s.%s_%s',...
    VisitID,ContainerInfo.StudyDate,ContainerInfo.StudyTime,ContainerUID);
  ContainerOutPath = sprintf('%s/%s',RootOutDir,ContainerOutDir);
  if ~exist(ContainerOutPath,'dir'), mkdir(ContainerOutPath); end
  if exist('ContainerInfo','var') & exist('ContainerOutPath','var')
    sp_idx = regexp(ContainerInfo.SeriesDescription,' ');
    SeriesName = ContainerInfo.SeriesDescription(1:sp_idx-1);
%% todo: use MMIL_Save_ContainerInfo
    save(sprintf('%s/ContainerInfo_%s.mat',ContainerOutPath,SeriesName),'ContainerInfo');
%    xml_save(sprintf('%s/ContainerInfo_%s.xml',ContainerOutPath,SeriesName),ContainerInfo);
  end

  tmp_indir = regexprep(indir,'\s','\\ '); % replace spaces
  tmp_indir = regexprep(tmp_indir,'(','\\('); % replace parenthesis
  tmp_indir = regexprep(tmp_indir,')','\\)'); % replace parenthesis
  cmd = sprintf('cp -r %s/* %s\n',tmp_indir,ContainerOutPath);

  fprintf('%s: cmd = %s\n',mfilename,cmd);
  [s,r] = unix(cmd);
  if s, fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r); end;

  volsfname = sprintf('%s/vols_%s.mat',ContainerOutPath,SeriesName);
  save(volsfname,'vols');
  fid_stat = fopen(statfname,'w');
  fprintf(fid_stat,'unpacked by %s at %f (%s) to\n%s',mfilename,now,datestr(now),ContainerOutPath);
  fclose(fid_stat);
end

%% todo: use MMIL_Stamp_Container
