function REC_MMIL_Import_PET_Exams(ProjID,forceflag)
% function REC_MMIL_Import_PET_Exams(ProjID,[forceflag])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Input:
%   forceflag: overwrite existing output
%
% Created:  09/24/08 by Alain Koyama
% Last Mod: 10/05/12 by Don Hagler
%

%% todo: rename as MMIL_Import_PET_Exams

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;
if isempty(ProjID)
   error('Empty ProjID');
end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
homedir = RootDirs.home;

batchdir = sprintf('%s/batchdirs/%s_%s',homedir,mfilename,ProjID);
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'file')
    cmd = sprintf('rm -rf %s/*\n',batchdir);
    fprintf('cmd = %s',cmd);
    system(cmd);
else
    mkdir(batchdir);
end
fid = fopen(scriptlistfname,'w'); fclose(fid);

jobnum = 0;
SourceRootDir = RootDirs.orig_pet;
ContainerRootDir = RootDirs.raw_pet;
if isempty(SourceRootDir) | isempty(ContainerRootDir)
   fprintf('Study %s does not appear to have PET scans\n',ProjInfo(p).ProjID);
   continue;
end
dirlist = dir(SourceRootDir);

for i=1:length(dirlist)
  SourceDir = dirlist(i).name;
  n = regexp(upper(SourceDir),'^[^\.]\w+');
  if dirlist(i).isdir & ~isempty(n)
    jobnum = jobnum+1;
    SubjID = SourceDir;
    jobID = sprintf('job_%03d',jobnum);
    jobfname = sprintf('%s/%s.m',batchdir,jobID);
    fid = fopen(jobfname,'w');
    fprintf(fid,'MMIL_Import_PET_Exam(''%s'',''%s'',''%s'',%d);\n',...
        sprintf('%s/%s',SourceRootDir,SourceDir),ContainerRootDir,SubjID,forceflag);
    fprintf(fid,'exit\n');
    fclose(fid);
    fid = fopen(scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
  end
end

fprintf('%%%% Now login to mmilcluster and run this:\n');
fprintf('    qmatjobs %s_%s\n',mfilename,ProjID);

