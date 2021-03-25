function REC_MMIL_Process_PET_Exams(ProjID,forceflag)
%function REC_MMIL_Process_PET_Exams(ProjID,[forceflag])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Input:
%   forceflag: overwrite existing output
%     {default = 0}
%
% Created:  09/24/08 by Alain Koyama
% Last Mod: 10/05/12 by Don Hagler
%

%% todo: rename as MMIL_Process_PET_Exams
%% todo: use MMIL_Check_ProjID
%% tood: use StudyInfo instead of creating dirlist here

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;
if isempty(ProjID)
  error('Empty ProjID');
end

[StudyInfo,RootDirs,ProjInfo]=REC_MMIL_Get_StudyInfo(ProjID)

batchdir = sprintf('%s/%s_%s',RootDirs.batch,ProjID,mfilename);

PET_ContainerDirs = {};
FS_ContainerDirs = {};
PET_ProcContainerDirs = {};

logdir = sprintf('%s/MetaData/%s',RootDirs.home,ProjID);
if ~exist(logdir,'dir')
  mkdir(logdir);
end
logfile = sprintf('%s/MetaData/%s/%s_RegisterPETtoMRI_Exams.log',...
  RootDirs.home,ProjID,ProjID);
flog = fopen(logfile,'wt');
if flog==-1
  fprintf('%s: ERROR: unable to create log file\n',mfilename);
  return;
end;

FS_dirlist = dir(sprintf('%s/FSURF*',RootDirs.fsurf));
PET_dirlist = dir(sprintf('%s/PETRAW*',RootDirs.raw_pet));

if isempty(RootDirs.raw_pet) | isempty(RootDirs.proc_pet)
  fprintf('%s: Study %s does not appear to have PET scans\n',mfilename,ProjID);
  continue;
end

for i = 1:length(PET_dirlist)
  PET_ContainerDir = char(PET_dirlist(i).name);
  %% todo: use regexprep, not (11:end)
  PET_ProcContainerDir = sprintf('%s/PETPROC_%s',...
    RootDirs.proc_pet,PET_ContainerDir(11:end));
  PET_files = dir(sprintf('%s/%s/vols*',RootDirs.raw_pet,PET_ContainerDir));
  n = regexp(char(PET_ContainerDir),'PETRAW_(?<SubjID>[^]+)_\d{8}','names');
  SubjID = n.SubjID;
  n = regexp(char(PET_ContainerDir),'(?<StudyDate>\d{8}(?=\.))','names');
  StudyDate = n.StudyDate;
  numvols = 0;
  volnames = {};
  petregstems = {};
  for j=1:length(PET_files)
    numvols = numvols + 1;
    volnames{end+1} = PET_files(j).name(6:end-4);
  end

  for j=1:length(volnames)
    petregstems{end+1} = sprintf('PET_reg_%s',volnames{j});
  end

  if ~numvols
    fprintf('%s: WARNING: missing vols.mat... skipping %s_%s\n',...
      mfilename,SubjID,StudyDate);
    fprintf(flog,'WARNING: missing vols.mat... skipping %s_%s\n',...
      SubjID,StudyDate);
    continue;  
  end

  PET_ContainerPath = sprintf('%s/%s',RootDirs.raw_pet,PET_ContainerDir);
  if isempty(PET_ContainerPath)
    fprintf('%s: WARNING: missing PETRAW Container... skipping %s_%s\n',...
      mfilename,SubjID,StudyDate);
    fprintf(flog,'WARNING: missing PETRAW Container... skipping %s_%s\n',...
      SubjID,StudyDate);
    continue;
  end;

  if isempty(PET_ProcContainerDir)
    mkdir(PET_ProcContainerDir);
  end;

  studymatch = {};
  for j=1:length(FS_dirlist) %find matching studies and get earliest one as baseline
    n = regexp(char(FS_dirlist(j).name),'FSURF_(?<SubjID>[^]+)_\d{8}','names');
    SubjIDFS = n.SubjID;
    if strcmp(SubjIDFS,SubjID)
      studymatch{end+1} = FS_dirlist(j).name;
    end
  end
  if ~isempty(studymatch)
    studymatch = sort(studymatch);
    FS_ContainerPath = sprintf('%s/%s',RootDirs.fsurf,studymatch{1});
  else
    fprintf('%s: WARNING: no FS recon for baseline MRI... skipping %s\n',...
      mfilename,SubjID);
    fprintf(flog,'WARNING: no FS recon for baseline MRI... skipping %s\n',...
      SubjID);
    continue;
  end

  fname = sprintf('%s/mri/nu.mgz',FS_ContainerPath);
  if ~exist(fname,'file')
    fprintf('%s: WARNING: missing %s... skipping %s_%s\n',...
      mfilename,fname,SubjID,StudyDate);
    fprintf(flog,'WARNING: missing %s... skipping %s_%s\n',...
      fname,SubjID,StudyDate);
    continue;
  end;

  for j=1:length(petregstems)
    fname = sprintf('%s/%s.mgh',PET_ProcContainerDir,petregstems{j});
    if exist(fname,'file')
      regfile = sprintf('%s/%s.txt',PET_ProcContainerDir,petregstems{j});
      if exist(regfile,'file')
        fid = fopen(regfile,'rt');
        if fid==-1
          fprintf('%s: ERROR: unable to open PET registration file %s\n',...
            mfilename,regfile);
          return;
        end;
        tline = fgetl(fid);
        fclose(fid);
        k = findstr(tline,'/');
        tmp_FS_ContainerPath = tline(k(1):k(end-1)-1);
        if strcmp(tmp_FS_ContainerPath,FS_ContainerPath) & ~forceflag
           % if already registered to correct FS Container, do not generate script below
          continue;
        elseif ~strcmp(tmp_FS_ContainerPath,FS_ContainerPath)
          fprintf('%s: WARNING: %s_%s registered to wrong FS Container\n',...
            mfilename,SubjID,StudyDate);
          continue;
        end
      end
    end
    FS_ContainerDirs{i} = FS_ContainerPath;
    PET_ContainerDirs{i} = PET_ContainerPath;
    PET_ProcContainerDirs{i} = PET_ProcContainerDir;
  end
end
fclose(flog);

scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'file')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf(1,'cmd = %s',cmd);
  unix(cmd);
else
  mkdir(batchdir);
end
fid = fopen(scriptlistfname,'w'); fclose(fid);
jobnum = 0;

for i = 1:length(PET_ContainerDirs)
  PET_ContainerDir = PET_ContainerDirs{i};
  FS_ContainerDir = FS_ContainerDirs{i};
  PET_ProcContainerDir = PET_ProcContainerDirs{i};
  if ~isempty(FS_ContainerDir) & ~isempty(PET_ContainerDir)
    jobnum = jobnum+1;
    jobID = sprintf('job_%03d',jobnum);
    jobfname = sprintf('%s/%s.m',batchdir,jobID);
    fid = fopen(jobfname,'w');
    fprintf(fid,'MMIL_Process_PET_Exam(''%s'',''%s'',''%s'',%d);\n',...
      PET_ContainerDir,PET_ProcContainerDir,FS_ContainerDir,forceflag);
    fprintf(fid,'exit\n');
    fclose(fid);
    fid = fopen(scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
  end
end

fprintf('%%% Now login to mmilcluster and run this:\n',mfilename);
fprintf('    qmatjobs %s_%s\n',mfilename,ProjID);
