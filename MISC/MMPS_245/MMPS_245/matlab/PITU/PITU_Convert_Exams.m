function PITU_Convert_Exams(forceflag)

if ~exist('forceflag','var')
  forceflag = false;
end

ContainerRootDir = '/space/pogo1/4/pub/tmp/dnr_atle';

dirlist = dir(sprintf('%s/SUBJECT*',ContainerRootDir));
batchdir = '/home/dale/batchdirs/PITU_Convert_Exams';
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'file')
  cmd = sprintf('rm -f %s/*\n',batchdir);
  fprintf(1,'cmd = %s',cmd);
  system(cmd);
else
  mkdir(batchdir);
end
fid = fopen(scriptlistfname,'w'); fclose(fid);
jobnum = 0;
for i = 1:length(dirlist)
  jobnum = jobnum+1;
  ContainerDir = char(dirlist(i).name);
  jobID = sprintf('job_%03d',jobnum);
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  fid = fopen(jobfname,'w');
  fprintf(fid,'PITU_Convert_Exam(''%s/%s'',%d);\n',ContainerRootDir,ContainerDir,forceflag);
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

