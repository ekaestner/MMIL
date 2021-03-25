function abcd_pc_all(indir,varargin)
%function abcd_pc_all(indir,[options])
%
% required input:
%   indir: input directory

% optional input:
%   'outdir': output directory
%     {default = indir}
%   'batchname': name of batch directory
%     {default = 'pc'}
%   'forceflag': [0|1] force recheck even if completed previously
%     {default = 0}
%
% Created:  08/15/16 by Jose Antolin
% Prev Mod: 09/29/16 by Jose Antolin
% Last Mod: 08/07/17 by Don Hagler
%

% NOTE: based on QC_PC_paramsbatch, last mod 08/15/16, created by Jose Antolin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'indir',indir,[],...
  'outdir',indir,[],...
  'batchname','pc',[],...
  'forceflag',false,[false true],...
});

% create output batch directory
root_batchdir = sprintf('%s/batchdirs',getenv('HOME'));
mmil_mkdir(root_batchdir);
batchdir = sprintf('%s/batchdirs/%s',getenv('HOME'),parms.batchname);
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: cmd %s failed:\n%s',mfilename,cmd,result);
  end;
end;

% get all directories in indir
fprintf('%s: finding directories in %s...\n',mfilename,parms.indir);
dlist = dir(parms.indir);
dlist = dlist([dlist.isdir] & ~strncmpi('.', {dlist.name}, 1)); 
ndirs = length(dlist);
if ~ndirs
  fprintf('%s: ERROR: no input dirs in %s!\n',mfilename,parms.indir);
  return;
end;

% create scriptlist
mmil_mkdir(batchdir)
fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

% create jobs
j=1;
for i=1:ndirs
  VisitID = dlist(i).name;
  st_indir = sprintf('%s/%s',parms.indir,VisitID);
  st_outdir = sprintf('%s/%s',parms.outdir,VisitID);
  fname_check = sprintf('%s/.finished',st_outdir);
  if exist(fname_check,'file') && ~parms.forceflag
    fprintf('%s: skipping %s (already completed)...\n',mfilename,VisitID);
    continue;
  end;
  % create job
  jstem = VisitID;
  if length(jstem)>30, jstem = jstem(1:30); end;
  jstem = regexprep(jstem,'[\^\-]','_');
  fprintf('%s: creating job for %s...\n',mfilename,jstem);
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  mmil_write_script(jobfname,'abcd_pc',{st_indir,st_outdir},[],[]);
  % add to list
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if j==1
  fprintf('%s: WARNING: no jobs created for %s\n',mfilename,parms.indir);
else
  fprintf('\n%%%% Now login to a cluster and run this:\n');
  fprintf('    qmatjobs %s\n',parms.batchname);
  fprintf('\n%%%% Or run this:\n');
  fprintf('    bmatjobs %s\n',parms.batchname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

