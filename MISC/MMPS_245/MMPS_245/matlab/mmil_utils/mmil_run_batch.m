function mmil_run_batch_jobs(batchdir)
%function mmil_run_batch_jobs(batchdir)
%
% Purpose: run batch jobs serially in single instance
%  of matlab
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/07/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if mmil_isrelative(batchdir)
  batchdir = ['/home/' getenv('USER') '/batchdirs/' batchdir];
end;

fname_scriptlist = [batchdir '/scriptlist.txt'];
if ~exist(fname_scriptlist,'file')
  error('file %s not found',fname_scriptlist);
end;

scriptlist = mmil_readtext(fname_scriptlist);

olddir = pwd;
cd(batchdir);

for i=1:length(scriptlist)
  job = scriptlist{i};
  orig_job = job;

  job_fname = [batchdir '/' job '.m'];
  if ~exist(job_fname,'file')
    fprintf('%s: WARNING: job script %s not found\n',...
      mfilename,job_fname);
    continue;
  end;
  
  % read script
  job_cmds = [];
  fid = fopen(job_fname);
  while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    job_cmds{end+1} = tline;
  end
  fclose(fid);

  % exclude final "exit"
  ncmds = length(job_cmds);
  if ~isempty(regexp(job_cmds{ncmds},'exit'))
    exit_flag = 1;
    tmp_job_fname = [batchdir '/' job '_tmp.m'];
    fid = fopen(tmp_job_fname,'wt');
    if fid<0, error('failed to create temporary file %s',tmp_job_fname); end;
    for i=1:ncmds-1
      fprintf(fid,'%s\n',job_cmds{i});
    end;
    fclose(fid);
    job = [job '_tmp'];
  else
    exit_flag = 0;
  end;

  err_flag = 0;
  fprintf('#####################################################\n');
  fprintf('#running job %s...\n\n',orig_job);
  try
    run(job);
  catch
    fprintf('%s: WARNING: job %s ended with error:\n%s\n',...
      mfilename,orig_job,lasterr);
    err_flag = 1;
  end;
  fprintf('\n#finished job %s\n',orig_job);

  % remove tmp job
  if exit_flag && ~err_flag
    delete(tmp_job_fname);
  end;
end;

cd(olddir);

