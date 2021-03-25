function REC_MMIL_Analyze_PET_Exams(ProjID,normflag,forceflag)
%function REC_MMIL_Analyze_PET_Exams(ProjID,normflag,forceflag)
%
% Required Parameters:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Parameters:
%   normflag - [0|1|2] whether to normalize to pons
%    (0 = absolute values, 1 = norm to pons, 2 = both)
%    {default = 1}
%   forceflag: [0|1] overwrite existing files
%    {default = 0}
%
% Created:  11/26/08 by Alain Koyama
% Last Mod: 10/05/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if isempty(ProjID)
  error('empty ProjID');
end;
if ~exist('normflag','var'), normflag = 1; end
if ~exist('forceflag','var'), forceflag = false; end;

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);

batchname = sprintf('%s_%s',ProjID,mfilename);
batchdir = sprintf('%s/%s',RootDirs.batch,batchname);
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

ContainerRootDir = RootDirs.proc_pet;
dirlist = dir(sprintf('%s/PETPROC*',ContainerRootDir));
if isempty(dirlist)
  fprintf('%s: Project %s has no PET Processed directory\n',mfilename,ProjID);
  continue;
end
for i = 1:length(dirlist)
  jobnum = jobnum+1;
  ContainerDir = char(dirlist(i).name);
  jobID = sprintf('job_%03d',jobnum);
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  fid = fopen(jobfname,'w');
  fprintf(fid,'REC_MMIL_Analyze_PET_Exam(''%s'',''%s'',%d,%d);\n',...
    ContainerRootDir,ContainerDir,normflag,forceflag);
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

fprintf('%%%% Now login to mmilcluster and run this:\n',mfilename);
fprintf('    qmatjobs %s\n',batchname);

