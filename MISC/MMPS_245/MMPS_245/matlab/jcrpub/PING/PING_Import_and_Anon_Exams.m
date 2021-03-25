function PING_Import_and_Anon_Exams(subjid_pattern,setup_file,forceflag)

%% import from rsync to orig for PING

%addpath('~cooper/MMPS_mat/jcrpub/PING');

if ~exist('subjid_pattern','var') || isempty(subjid_pattern) || ~ischar(subjid_pattern),
   error('SUBJID_PATTERN arg must be a valid string.');
end

%forceflag = 1;

STAT_OK = 0;
STAT_ERROR = 1;
STAT_PRINT = 2;

anon_fname = 'anon.status';

% load dir paths
run(setup_file);

if 0
rsyncdir = '/space/md10/10/data/MMILDB/PING/rsync';
anondir = '/space/md2/1/data/MMILDB/PING/orig_anon';
%anondir = '/space/md17/8/data/MMILDB/PING/tmp_anon';
unzip_topdir = '/home/mmilping/anonymize/scratch';
logdir = '/home/mmilping/anonymize/logs';
exclude_dir = '/home/mmilping/anonymize/exclude';
end

%fprintf(1,'SubjID pattern = %s\n',subjid_pattern);
if isempty(strfind(subjid_pattern,'*'))
   subjectdirs = {subjid_pattern};
else
   subjectdirs = dir(sprintf('%s/%s',rsyncdir, subjid_pattern));
   subjectdirs = {subjectdirs.name};
end
if isempty(subjectdirs)
   fprintf(1,'subjectdirs is empty.\n');
end

dstr = datestr(now,'yyyymmdd');
[status,ret_str] = unix(sprintf('mktemp -p %s PING_anon_%s.log.XXXX',...
                   logdir,dstr));
if status,
   fprintf(1,'%s: Unable to create logfile in %s',mfilename,logdir);
end
anonlog = deblank(ret_str);
unix(sprintf('chmod a+r %s',anonlog));
%fprintf(1,'anonlog = %s\n',anonlog);

for s = 1:length(subjectdirs)
    subjid = subjectdirs{s};
    subjpath = sprintf('%s/%s',rsyncdir,subjid);
    visit_dirs = dir(subjpath);
    for ii=1:length(visit_dirs)
       visdir = visit_dirs(ii).name;
       if isempty(regexp(visdir,'\d\d\d')), continue; end
       vispath = sprintf('%s/%s',subjpath,visdir);
       attempt_dirs = dir(vispath);
       for jj=1:length(attempt_dirs)
            attdir = attempt_dirs(jj).name;
            if isempty(regexp(attdir,'\d\d')), continue; end
            attpath = sprintf('%s/%s',vispath,attdir);
            outdir = sprintf('%s_%s_%s',subjid,visdir,attdir);
            proc_status(STAT_PRINT,sprintf('Processing %s...',outdir),anonlog);

            % Check if this attempt is in our exclude list
            if exist(sprintf('%s/%s',exclude_dir,outdir),'file'), 
               proc_status(STAT_PRINT,sprintf('%s is in exclude list...skipping.',outdir),anonlog);
               continue;  
            end

            anondirpath = sprintf('%s/%s',anondir,outdir);
            if exist(anondirpath,'dir')
               if forceflag,
                  proc_status(STAT_PRINT,sprintf('Anon dir %s exists...reprocessing now.',...
                         outdir),anonlog);
                  [status,result] = unix(sprintf('rm -rf %s',anondirpath));
                  proc_status(status,deblank(result),anonlog);
                  if status, continue; end
               else
                  proc_status(STAT_PRINT,sprintf('Anon dir %s exists...skipping.',...
                         outdir),anonlog);
                  continue;   % Already anonymized
               end
            end

            % Try unzipping the zipfile to the scratch dir and copy ScanInfo.txt also
            [unzippath,status,result] = PING_unzip_exam(attpath,unzip_topdir,outdir);

            if isempty(unzippath) || status,
               errmsg = sprintf('%s: unzip halted: %s',outdir,result);
               proc_status(STAT_ERROR,errmsg,anonlog);
               continue;
            end

            Anon = PING_anonymize(unzippath);
            anon_statfile = sprintf('%s/%s',unzippath,anon_fname);

            if Anon.AnonFileCount == 0,
               errmsg = sprintf('%s: anon halted: No DICOMs anonymized.',outdir);
               if isempty(Anon.ErrMsg), Anon.ErrMsg = errmsg;
               else Anon.ErrMsg = sprintf('%s | %s.',Anon.ErrMsg,errmsg);
               end
            end
    
            if ~isempty(Anon.ErrMsg)
               proc_status(STAT_ERROR,Anon.ErrMsg,anonlog);
               PING_Save_Anon_Status(anon_statfile,Anon);
               continue;
            end

            Anon.AnonOK = 1;
            errmsg = PING_Save_Anon_Status(anon_statfile,Anon);
            if ~isempty(errmsg)
               proc_status(STAT_ERROR,errmsg,anonlog);
               continue;
            end

            % cp unzippath dir to anondir
            [status,result] = unix(sprintf('cp -rp %s %s',unzippath,anondirpath));
            proc_status(status,deblank(result),anonlog);
            if status, continue; end

            % Give a+rx access to anondirpath + subdirs, and set go-w permission
            cmd = sprintf('find %s -type d -exec chmod u+rx,g+rx,g-ws,o+rx,o-w ''{}'' \\;',anondirpath);
            [status,result] = unix(cmd);
            proc_status(status,deblank(result),anonlog);
            if status, continue; end

            % Allow owner to write to top-level anon dir (hopefully not necessary in future)
            [status,result] = unix(sprintf('chmod u+w %s',anondirpath));
            [status,result] = unix(cmd);
            proc_status(status,deblank(result),anonlog);
            if status, continue; end

            % Set g- and o- wxs permission on all files in anondirpath
            cmd = sprintf('find %s -type f -exec chmod u+r,u-x,g+r,g-wxs,o+r,o-wxs ''{}'' \\;',anondirpath);
            [status,result] = unix(cmd);
            proc_status(status,deblank(result),anonlog);
            if status, continue; end

            % Set perms for anon.status file
            cmd = sprintf('chmod u+r,u-x,g+r,g-wxs,o+r,o-wxs %s/%s',anondirpath,anon_fname);
            [status,result] = unix(cmd);
            proc_status(status,deblank(result),anonlog);
            if status, continue; end

            % rm unzippath
            cd(unzip_topdir);
            [status,result] = unix(sprintf('rm -rf %s',unzippath));
            proc_status(status,deblank(result),anonlog);
            if status, continue; end

            proc_status(STAT_PRINT,sprintf('Proc completed for %s.',outdir),anonlog);
       end
   end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function proc_status(status,result,logfile)

 if status,
   result = strrep(result,sprintf('\n'),'  ');
   fprintf(1,'%s\n',result);
   unix(sprintf('echo "%s" >> %s',result,logfile));
 end

