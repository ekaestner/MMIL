function abcd_import_all(varargin)
%function abcd_import_all([options])
%
% Optional Input: ('key', value pairs)
%   'batchname': name of batch directory in /home/mmilrec/batchdirs
%     {default = 'abcd_import_all'}
%   'batch_append_flag': [0|1] add jobs to existing batchdir
%     otherwise delete any existing jobs
%     {default = 0}
%   'batch_msg_flag': [0|1] display message after creating jobs
%     instructing user to submit jobs to cluster
%     {default = 1}
%   'indir': input directory containing json
%     and tar/tgz files from ABCD FIONA
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD/incoming'}
%   'unpackdir': output directory to be filled
%     with contents of tar/tgz files
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD/unpack'}
%   'outdir': output directory that will contain unpacked data
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD/orig'}
%   'fname_sites': spreadsheet containing site names and IDs
%     {default = '/home/mmilrec/ProjInfo/ABCD/ABCD_Sites.csv'}
%   'phantom_flag': [0|1|2] whether to import phantom scans
%     0: exclude phantom scans (pGUID matches ABCDPhantomXXX)
%     1: allow phantom scans
%     2: only phantom scans
%     {default = 0}
%   'humanphantom_flag': [0|1|2] whether to import human phantom scans
%     0: exclude human phantom scans (pGUID contains PhantomTravelingHuman or PhantomLocalHuman))
%     1: allow human phantom scans
%     2: only human phantom scans
%     {default = 0}
%   'testQA_flag': [0|1|2] whether to import test QA scans
%     0: exclude test QA scans (pGUID contains TestQA)
%     1: allow test QA scans
%     2: only test QA scans
%     {default = 0}
%   'PMU_flag':  [0|1|2] whether to import test PMU files
%     0: exclude PMU files (SeriesDescription contains PMU)
%     1: allow PMU files
%     2: only PMU files
%     {default = 0}
%   'require_subjid_flag': [0|1] whether to skip series without a SubjID in fstem
%     {default = 1}
%   'mb_recon_flag': [-1|0|1|2] whether to perform multi-band recon
%     -1: unpack only, this is for reconned data
%     0: no multi-band recon
%     1: unpack and perform multi-band recon if required
%     2: do not unpack unless multi-band recon is required
%     {default = 1}
%   'mb_tmpdir': temporary directory that will contain untar'd P-files
%     {default = '/scratch'}
%   'mb_dicoms_flag': [0|1] whether to write dicoms after multi-band recon
%     0: produce mat files only
%     1: produce mat files and dicoms
%     2: write dicoms only (do not produce jobs if recon has not been run yet)
%     {default = 1}
%   'mb_cleanup_flag': [0|1|2] whether to cleanup temporary mb output
%     0: perform no cleanup
%     1: remove P-file when complete
%     2: remove mb_tmpdir when complete
%     {default = 2}
%   'separate_series_flag': [0|1] create separate output directories
%     for each series
%     {default = 0}
%   'check_dates_flag': [0|1] compare date of unpack to tgz in incoming
%     force import if unpack is older than tgz
%     {default = 0}
%   'check_md5sum_flag': [0|1] compare contents of md5sum files in incoming
%     force import if there is a mismatch with file in unpack
%     {default = 0}
%   'check_md5sum_server_flag': [0|1] compare contents of md5sum and md5sum_server files in incoming
%     skip import if there is a mismatch or one or more does not exist
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%   'force_err_flag': [0|1] overwrite existing output if previous errors
%     {default = 0}
%   'force_unpack_err_flag': [0|1] overwrite existing output if previous unpack errors
%     {default = 0}
%   'force_recon_err_flag': [0|1] overwrite existing output if previous recon errors
%     {default = 0}
%   'force_dicoms_err_flag': [0|1] overwrite existing output if previous dicoms errors
%     for multi-band recon only
%     {default = 0}
%   'force_rename_err_flag': [0|1] overwrite existing output if previous rename errors
%     {default = 0}
%   'logfile': name of file created to log messages
%     {default = []}
%
% Output:
%   datadir: data directory containing imported dicoms
%
% Created:  07/17/16 by Don Hagler
% Prev Mod: 08/25/17 by Don Hagler
% Last Mod: 10/18/17 by Feng Xue
%

%% todo: get json file info from csv
%%       identify tgz files without json files, extract for their json info
%%       then make jobs, without loading each json file individually in loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);

% create output batch directory
parms = init_batchdir(parms);

% create logfile
parms = create_logfile(parms);

% find tgz files
parms = find_tgz_files(parms);
if isempty(parms.fnames)
  close_logfile(parms);
  return;
end;

% create batchdir
parms = create_batchdir(parms);

% create jobs
parms = create_jobs(parms);

% display message prompting to submit jobs to cluster
display_batch_msg(parms);

% close log file
close_logfile(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'batchname','abcd_import_all',[],...
    'batch_append_flag',false,[false true],...
    'batch_msg_flag',true,[false true],...
    'indir','/space/syn05/1/data/MMILDB/DAL_ABCD/incoming',[],...
    'unpackdir','/space/syn05/1/data/MMILDB/DAL_ABCD/unpack',[],...
    'outdir','/space/syn05/1/data/MMILDB/DAL_ABCD/orig',[],...
    'fname_sites','/home/mmilrec/ProjInfo/ABCD/ABCD_Sites.csv',[],...
    'phantom_flag',0,[0:2],...
    'humanphantom_flag',0,[0:2],...
    'testQA_flag',0,[0:2],...
    'PMU_flag',0,[0:2],...
    'require_subjid_flag',true,[false true],...
    'mb_recon_flag',1,[-1:2],...
    'mb_tmpdir','/scratch',[],...
    'mb_dicoms_flag',1,[0:2],...
    'mb_cleanup_flag',2,[0:2],...
    'separate_series_flag',false,[false true],...
    'check_dates_flag',false,[false true],...
    'check_md5sum_flag',false,[false true],...
    'check_md5sum_server_flag',false,[false true],...
    'forceflag',false,[false true],...
    'force_err_flag',false,[false true],...
    'force_unpack_err_flag',false,[false true],...
    'force_recon_err_flag',false,[false true],...
    'force_dicoms_err_flag',false,[false true],...
    'force_rename_err_flag',false,[false true],...
    'logfile',[],[],...
    ...
    'md5sum_ext','.md5sum',[],...
    'md5sum_server_ext','.md5sum_server',[],...
    'phantom_pat_list',{'ABCDPhantom'},[],...
    'humanphantom_pat_list',{'PhantomTravelingHuman','PhantomLocalHuman'},[],...
    'testQA_pat_list',{'TestQA'},[],...
    'fsize_thresh',1,[],...
    ...
    'tags',{'indir','unpackdir','outdir','fname_sites',...
            'mb_recon_flag','mb_tmpdir','mb_dicoms_flag','mb_cleanup_flag',...
            'separate_series_flag','forceflag','force_unpack_err_flag',...
            'force_recon_err_flag','force_dicoms_err_flag',...
            'force_rename_err_flag'},[],...
  });

  if parms.phantom_flag==2 && parms.humanphantom_flag==2
    error('phantom_flag==2 and humanphantom_flag==2 are mutually exclusive');
  end;
  if parms.phantom_flag==2 && parms.testQA_flag==2
    error('phantom_flag==2 and humanphantom_flag==2 are mutually exclusive');
  end;
  if parms.humanphantom_flag==2 && parms.testQA_flag==2
    error('phantom_flag==2 and humanphantom_flag==2 are mutually exclusive');
  end;
  if parms.force_err_flag
    parms.force_unpack_err_flag = 1;
    parms.force_recon_err_flag = 1;
    parms.force_dicoms_err_flag = 1;
    parms.force_rename_err_flag = 1;
  end;
  parms.force_err_flag = (parms.force_unpack_err_flag ||...
                          parms.force_recon_err_flag ||...
                          parms.force_dicoms_err_flag ||...
                          parms.force_rename_err_flag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = init_batchdir(parms)
  root_batchdir = sprintf('%s/batchdirs',getenv('HOME'));
  mmil_mkdir(root_batchdir);
  parms.batchdir = sprintf('%s/batchdirs/%s',getenv('HOME'),parms.batchname);
  parms.scriptlistfname = sprintf('%s/scriptlist.txt',parms.batchdir);
  if exist(parms.batchdir,'dir') && ~parms.batch_append_flag
    cmd = sprintf('rm -rf %s\n',parms.batchdir);
    fprintf('cmd = %s',cmd);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: WARNING: cmd %s failed:\n%s',mfilename,cmd,result);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_batchdir(parms)
  mmil_mkdir(parms.batchdir)
  if parms.batch_append_flag
    fid = fopen(parms.scriptlistfname,'a');
    if fid==-1
      error('failed to open scriptlist file %s for appending\n',parms.scriptlistfname);
    end;
    fclose(fid);
    % read file, get last job number
    try
      joblist = mmil_readtext(parms.scriptlistfname);
    catch
      joblist = [];
    end;
    if isempty(joblist)
      parms.j = 1;
    else
      jobID = joblist{end};
      n = regexp(jobID,'job_(?<jobnum>\d+)','names');
      if isempty(n)
        error('failed to read job number from %s',parms.scriptlistfname);
      end;
      parms.j = str2num(n.jobnum) + 1;
    end;
  else
    fid = fopen(parms.scriptlistfname,'w');
    if fid==-1
      error('failed to open scriptlist file %s for writing\n',parms.scriptlistfname);
    end;
    fclose(fid);
    parms.j=1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_logfile(parms)
  if ~isempty(parms.logfile)
    logdir = fileparts(parms.logfile);
    mmil_mkdir(logdir);
    parms.flog = fopen(parms.logfile,'a');
    if parms.flog==-1
      error('failed to open logfile %s for writing\n',parms.logfile);
    end;
  else
    parms.flog = -1;
  end;
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('%s: starting abcd import for %s at %s...\n',...
    mfilename,parms.indir,datestr(now));
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  if parms.flog>0
    fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(parms.flog,'%s: starting abcd import for %s at %s...\n',...
      mfilename,parms.indir,datestr(now));
    fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = find_tgz_files(parms)
  parms.fnames = [];
  
  % check for tgz files
  dlist = dir(sprintf('%s/*.tgz',parms.indir));
  if isempty(dlist)
    fprintf('%s: ERROR: no tgz files in %s!\n',mfilename,parms.indir);
    if parms.flog>0
      fprintf(parms.flog,'%s: ERROR: no tgz files in %s!\n',...
        mfilename,parms.indir);
    end;
    return;
  end;
  fnames = {dlist.name};

  % identify fnames that have phantom pattern (or not)
  if ismember(parms.phantom_flag,[0,2])
    ind_match = find_matches(fnames,parms.phantom_pat_list,~parms.phantom_flag);
    fnames = fnames(ind_match);
  end;

  % identify fnames that have human phantom pattern (or not)
  if ismember(parms.humanphantom_flag,[0,2])
    ind_match = find_matches(fnames,parms.humanphantom_pat_list,~parms.humanphantom_flag);
    fnames = fnames(ind_match);
  end;

  % identify fnames that have testQA pattern (or not)
  if ismember(parms.testQA_flag,[0,2])
    ind_match = find_matches(fnames,parms.testQA_pat_list,~parms.testQA_flag);
    fnames = fnames(ind_match);
  end;

  % update list of files
  parms.fnames = fnames;
  
  if isempty(parms.fnames)
    fprintf('%s: ERROR: no tgz files meeting criteria in %s!\n',...
      mfilename,parms.indir);
    if parms.flog>0
      fprintf(parms.flog,'%s: ERROR: no tgz files meeting criteria in %s!\n',...
        mfilename,parms.indir);
    end;
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_match = find_matches(fnames,pat_list,nonmatch_flag)
  ind_match = [];
  for i=1:length(pat_list)
    pat = pat_list{i};
    ind_tmp = find(~cellfun(@isempty,regexp(fnames,pat)));
    if ~isempty(ind_tmp)
      ind_match = union(ind_match,ind_tmp);
    end;
  end;
  if nonmatch_flag
    ind_match = setdiff([1:length(fnames)],ind_match);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ismatch = check_subj_pat(SubjID,pat_list)
  ismatch = 0;
  for i=1:length(pat_list)
    pat = pat_list{i};
    if ~isempty(regexp(SubjID,pat))
      ismatch = 1;
      return;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function display_batch_msg(parms)
  if parms.j==1
    fprintf('%s: WARNING: no jobs created for %s\n',mfilename,parms.indir);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no jobs created for %s\n',...
        mfilename,parms.indir);
    end;
  elseif parms.batch_msg_flag
    fprintf('\n%%%% Now login to a cluster and run this:\n');
    fprintf('    qmatjobs %s\n',parms.batchname);
    fprintf('\n%%%% Or run this:\n');
    fprintf('    bmatjobs %s\n',parms.batchname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: more subfunctions
function parms = create_jobs(parms)
  for i=1:length(parms.fnames)
    fname_tgz = sprintf('%s/%s',parms.indir,parms.fnames{i});

    % check that tgz file does not contain k-space data
    fstem_tgz = check_tgz_file(fname_tgz,parms);
    if isempty(fstem_tgz), continue; end;

    % check json file, extract from tgz if necessary
    fprintf('%s: checking json file for %s...\n',mfilename,fname_tgz)
    if parms.flog>0
      fprintf(parms.flog,'%s: checking json file for %s...\n',...
        mfilename,fname_tgz);
    end;
    tic;
    [jinfo,fstem_json,fname_tar,fname_json] = check_json_file(fname_tgz,parms);
    t = toc;
    fprintf('%0.2f sec\n',t);
    if parms.flog>0
      fprintf(parms.flog,'%0.2f sec\n',t);
    end;
    if isempty(jinfo), continue; end;    
    
    % set directory to contain output of untarring
    unpackdir = sprintf('%s/%s',parms.unpackdir,fstem_json);

    % skip if this is a json file for P file tar
    if jinfo.pfile_flag
      fprintf('%s: skipping import for %s because is a P file json file\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because is a P file json file\n',...
          mfilename,fstem_json);
      end;
      continue;
    end;

    % NOTE: commented this out to allow dicom writing if recon is complete
    %       but P-file tgz is no longer present
    % check for corresponding tar or tgz file
%    if isempty(fname_tar) || ~exist(fname_tar,'file')
%      fprintf('%s: skipping import for %s because missing tar/tgz file\n',...
%        mfilename,fstem_json);
%      if parms.flog>0
%        fprintf(parms.flog,'%s: skipping import for %s because missing tar/tgz file\n',...
%          mfilename,fstem_json);
%      end;
%      continue;
%    end;

    % skip if this json file does not have SeriesDescription (e.g. protocol compliance check json)
    if isempty(mmil_getfield(jinfo,'SeriesDescription'))
      fprintf('%s: skipping import for %s because json file does not have SeriesDescription\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because json file does not have SeriesDescription\n',...
          mfilename,fstem_json);
      end;
      continue;
    end;

    % skip if this is a json file for ScreenSave
    isScreenSave = ~isempty(regexp(jinfo.ImageType,'SCREEN SAVE'));
    if isScreenSave
      fprintf('%s: skipping import for %s because ImageType is SCREEN SAVE\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because ImageType is SCREEN SAVE\n',...
          mfilename,fstem_json);
      end;
      continue;
    end;

    % skip if this is a json file for PMU
    isPMU = ~isempty(regexp(jinfo.SeriesDescription,'PMU'));
    if parms.PMU_flag==0 && isPMU
      fprintf('%s: skipping import for %s because is PMU and PMU_flag = 0\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because is PMU and PMU_flag = 0\n',...
          mfilename,fstem_json);
      end;
      continue;
    elseif parms.PMU_flag==2 && ~isPMU
      fprintf('%s: skipping import for %s because is not PMU and PMU_flag = 2\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because is not PMU and PMU_flag = 2\n',...
          mfilename,fstem_json);
      end;
      continue;
    end;
    
    if parms.check_md5sum_server_flag
      % check whether transfer was complete
      %   by comparing md5sum and md5sum_server files
      if jinfo.recon_required_flag
        fname = fname_tar;
      else
        fname = fname_json;
      end;
      if ~isempty(fname)
        % find corresponding md5sum file
        [fpath,fstem,fext] = fileparts(fname);
        fname_md5sum_in = sprintf('%s/%s%s',fpath,fstem,parms.md5sum_ext);
        if ~exist(fname_md5sum_in)
          fprintf('%s: WARNING: md5sum file %s is missing\n',...
            mfilename,fname_md5sum_in);
          if parms.flog>0
            fprintf(parms.flog,'%s: WARNING: md5sum file %s is missing\n',...
              mfilename,fname_md5sum_in);
          end;
        else
          fname_md5sum_server = sprintf('%s/%s%s',...
            fpath,fstem,parms.md5sum_server_ext);
          if ~exist(fname_md5sum_server)
            fprintf('%s: WARNING: skipping import for %s because md5sum_server file %s is missing\n',...
              mfilename,fstem_json,fname_md5sum_server);
            if parms.flog>0
              fprintf(parms.flog,'%s: WARNING: skipping import for %s because md5sum_server file %s is missing\n',...
                mfilename,fstem_json,fname_md5sum_server);
            end;
            continue;
          end;
          mstrin = read_md5sum_file(fname_md5sum_in,parms);
          mstrsv = read_md5sum_file(fname_md5sum_server,parms);
          if ~isempty(mstrin) && ~isempty(mstrsv) && ~strcmp(mstrin,mstrsv)
            fprintf('%s: WARNING: skipping import for %s because md5sum %s and md5sum_server file %s do not match\n',...
              mfilename,fstem_json,fname_md5sum_in,fname_md5sum_server);
            if parms.flog>0
              fprintf(parms.flog,'%s: WARNING: skipping import for %s because md5sum %s and md5sum_server file %s do not match\n',...
                mfilename,fstem_json,fname_md5sum_in,fname_md5sum_server);
            end;
            continue;
          end;
        end;
      end;
    end;
    
    % get SubjID from fstem_json
    SubjID = abcd_get_SubjID(fstem_json);
    if isempty(SubjID)
      if parms.require_subjid_flag
        fprintf('%s: skipping import for %s because of missing SubjID\n',...
          mfilename,fstem_json);
        if parms.flog>0
          fprintf(parms.flog,'%s: skipping import for %s because of missing SubjID\n',...
            mfilename,fstem_json);
        end;
        continue;
      else
        SubjID = sprintf('S%d',i);
      end;
    end;

    % determine whether mb recon is required
    if parms.mb_recon_flag==2 && ~jinfo.recon_required_flag
      fprintf('%s: skipping import for %s because recon not required and mb_recon_flag = 2\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because recon not required and mb_recon_flag = 2\n',...
          mfilename,fstem_json);
      end;
      continue;
    elseif parms.mb_recon_flag==0 && jinfo.recon_required_flag
      fprintf('%s: skipping import for %s because recon required and mb_recon_flag = 0\n',...
        mfilename,fstem_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: skipping import for %s because recon required and mb_recon_flag = 0\n',...
          mfilename,fstem_json);
      end;
      continue;
    end;

    % skip if only writing dicoms and recon not complete
    if parms.mb_dicoms_flag==2 && jinfo.recon_required_flag
      fname_out = sprintf('%s/%s.reconned',unpackdir,fstem_json);
      if ~exist(fname_out,'file')
        fprintf('%s: skipping import for %s because recon is required, is not complete, and mb_dicoms_flag = 2\n',...
          mfilename,fstem_json);
        if parms.flog>0
          fprintf(parms.flog,'%s: skipping import for %s because recon is required, is not complete, and mb_dicoms_flag = 2\n',...
            mfilename,fstem_json);
        end;
        continue;
      end;
    else
      % NOTE: added this to allow dicom writing if recon is complete
      %       but P-file tgz is no longer present
      % check for corresponding tar or tgz file
      if isempty(fname_tar) || ~exist(fname_tar,'file')
        fprintf('%s: skipping import for %s because missing tar/tgz file\n',...
          mfilename,fstem_json);
        if parms.flog>0
          fprintf(parms.flog,'%s: skipping import for %s because missing tar/tgz file\n',...
            mfilename,fstem_json);
        end;
        continue;
      end;
    end;

    % check for previous errors
    if parms.force_err_flag
      err_type = [];
      if parms.force_unpack_err_flag
        dlist_unpack = dir(sprintf('%s/%s.unpack_error',unpackdir,fstem_json));
        err_type = 'upack';
      else
        dlist_unpack = [];
      end;
      if parms.force_recon_err_flag
        dlist_recon = dir(sprintf('%s/%s.recon_error',unpackdir,fstem_json));
        err_type = 'recon';
      else
        dlist_recon = [];
      end;
      if parms.force_rename_err_flag
        dlist_rename = dir(sprintf('%s/%s.rename_error',unpackdir,fstem_json));
        err_type = 'rename';
      else
        dlist_rename = [];
      end;
      if ~isempty(dlist_unpack) || ~isempty(dlist_recon) || ~isempty(dlist_rename)
        fprintf('%s: NOTE: forcing import for %s despite previous %s errors\n',...
          mfilename,fstem_json,err_type);
        % delete unpack dir
        fprintf('%s: deleting %s ...\n',...
          mfilename,unpackdir);
        cmd = sprintf('rm -r %s',unpackdir);
        [s,r] = unix(cmd);
        if s, error('cmd %s failed:\n%s',cmd,r); end;
        if parms.flog>0
          fprintf(parms.flog,'%s: NOTE: forcing import for %s despite previous %s errors\n',...
            mfilename,fstem_json,err_type);
          fprintf(parms.flog,'%s: deleting %s ...\n',...
            mfilename,unpackdir);
        end;
      end;
    end;

    % check whether file was imported previously
    dlist = cat(1,dir(sprintf('%s/%s.unpack*',unpackdir,fstem_json)),...
                  dir(sprintf('%s/%s.recon*',unpackdir,fstem_json)));
    if ~isempty(dlist)
      % compare dates of check files
      if parms.check_dates_flag
        for d=1:length(dlist)
          fname_check = sprintf('%s/%s',unpackdir,dlist(d).name);
          fs = dir(fname_check);
          ts = dir(fname_tar);
          if ts.datenum>fs.datenum
            [fdir,fstem,fext] = fileparts(fname_check);
            fprintf('%s: NOTE: forcing import for %s because %s file is older than tar file (%s vs %s)\n',...
              mfilename,fstem_json,fext,fs.date,ts.date);
            % delete unpack dir
            fprintf('%s: deleting %s ...\n',...
              mfilename,unpackdir);
            cmd = sprintf('rm -r %s',unpackdir);
            [s,r] = unix(cmd);
            if s, error('cmd %s failed:\n%s',cmd,r); end;
            if parms.flog>0
              fprintf(parms.flog,'%s: NOTE: forcing import for %s because %s file is older than tar file (%s vs %s)\n',...
                mfilename,fstem_json,fext,fs.date,ts.date);
              fprintf(parms.flog,'%s: deleting %s ...\n',...
                mfilename,unpackdir);
            end;
            break;
          end;
        end;
      end;

      % check md5sum file
      if parms.check_md5sum_flag
        if jinfo.recon_required_flag
          fname = fname_tar;
        else
          fname = fname_json;
        end;
        % find corresponding md5sum file
        [fpath,fstem,fext] = fileparts(fname);
        fname_md5sum_in = sprintf('%s/%s%s',fpath,fstem,parms.md5sum_ext);
        if ~exist(fname_md5sum_in)
          fprintf('%s: WARNING: md5sum file %s is missing\n',...
            mfilename,fname_md5sum_in);
          if parms.flog>0
            fprintf(parms.flog,'%s: WARNING: md5sum file %s is missing\n',...
              mfilename,fname_md5sum_in);
          end;
        else
          % replace spaces with underscores
          fstem = regexprep(fstem,' ','_');
          fname_md5sum_out = sprintf('%s/%s%s',unpackdir,fstem,parms.md5sum_ext);
          if ~exist(fname_md5sum_out)
            cmd = sprintf('cp -p %s %s',...
              regexprep(fname_md5sum_in,' ','\\ '),fname_md5sum_out);
            fprintf('%s: copying md5sum file %s to %s...\n',...
              mfilename,fname_md5sum_in,fname_md5sum_out);
            if parms.flog>0
              fprintf(parms.flog,'%s: copying md5sum file %s to %s...\n',...
                mfilename,fname_md5sum_in,fname_md5sum_out);
            end;
            [s,r] = unix(cmd);
            if s
              error('failed to copy md5sum file:\n%s\n%s\n',cmd,r);
            end;
          end;
          % compare contents of files
          mstrin = read_md5sum_file(fname_md5sum_in,parms);
          mstrout = read_md5sum_file(fname_md5sum_out,parms);
          if ~isempty(mstrin) && ~isempty(mstrout) && ~strcmp(mstrin,mstrout)
            fprintf('%s: NOTE: forcing import for %s because md5sum files %s and %s do not match\n',...
              mfilename,fstem_json,fname_md5sum_in,fname_md5sum_out);
            % delete unpack dir
            fprintf('%s: deleting %s ...\n',...
              mfilename,unpackdir);
            cmd = sprintf('rm -r %s',unpackdir);
            [s,r] = unix(cmd);
            if s, error('cmd %s failed:\n%s',cmd,r); end;
            if parms.flog>0
              fprintf(parms.flog,'%s: NOTE: forcing import for %s because md5sum files do not match\n',...
                mfilename,fstem_json);
              fprintf(parms.flog,'%s: deleting %s ...\n',...
                mfilename,unpackdir);
            end;
          else
            fprintf('%s: NOTE: md5sum files matched for %s\n',...
              mfilename,fstem_json);
            if parms.flog>0
              fprintf(parms.flog,'%s: md5sum files matched for %s...\n',...
                mfilename,fstem_json);
            end;
          end;
        end;
      end;

      % check for errors
      dlist = dir(sprintf('%s/%s.*_error',unpackdir,fstem_json));
      if ~isempty(dlist)
        fprintf('%s: WARNING: skipping import for %s because of previous errors\n',...
          mfilename,fstem_json);
        if parms.flog>0
          fprintf(parms.flog,'%s: WARNING: skipping import for %s because of previous errors\n',...
            mfilename,fstem_json);
        end;
        continue;
      end;

      % check for unpacked file again
      dlist = dir(sprintf('%s/%s.unpack*',unpackdir,fstem_json));
      if ~isempty(dlist)
        if parms.mb_recon_flag && ~parms.mb_dicoms_flag && jinfo.recon_required_flag
          dlist = dir(sprintf('%s/%s.recon*',unpackdir,fstem_json));
          if ~isempty(dlist) && ~parms.forceflag
            fprintf('%s: skipping import for %s because recon complete\n',...
              mfilename,fstem_json);
            if parms.flog>0
              fprintf(parms.flog,'%s: skipping import for %s because recon complete\n',...
                mfilename,fstem_json);
            end;
            continue;
          end;
        else
          dlist = dir(sprintf('%s/%s.rename*',unpackdir,fstem_json));
          % check if already renamed
          if ~isempty(dlist) && ~parms.forceflag
            fprintf('%s: skipping import for %s because rename complete\n',...
              mfilename,fstem_json);
            if parms.flog>0
              fprintf(parms.flog,'%s: skipping import for %s because rename already complete\n',...
                mfilename,fstem_json);
            end;
            continue;
          end;
        end;
      end;
    end;

    % get other info
    StudyDate = jinfo.StudyDate;
    SeriesNumber = jinfo.SeriesNumber;

    % create job
    jstem = sprintf('%s_%s_se%s',SubjID,StudyDate,SeriesNumber);
    jstem = regexprep(jstem,'\^','_');
    fprintf('%s: creating job for %s...\n',mfilename,jstem);
    if parms.flog>0
      fprintf(parms.flog,'%s: creating job for %s...\n',mfilename,jstem);
    end;
    jobID = sprintf('job_%03d_%s',parms.j,jstem);
    parms.j = parms.j + 1;
    jobfname = sprintf('%s/%s.m',parms.batchdir,jobID);
    mmil_write_script(jobfname,'abcd_import',...
      fname_json,parms.tags,parms);
    % add to list
    fid = fopen(parms.scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem_tgz = check_tgz_file(fname_tgz,parms)
  [indir,fstem_tgz] = fileparts(fname_tgz);
  % check if filestem contains SUID because these contain k-space data
  if ~isempty(regexp(fstem_tgz,'SUID'))
    fprintf('%s: WARNING: skipping import for %s because file stem contains ''SUID''\n',...
      mfilename,fstem_tgz);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: skipping import for %s because file stem contains ''SUID''\n',...
        mfilename,fstem_tgz);
    end;
    fstem_tgz = [];
    return;
  end;
  % check tgz file size and exclude if too big (contains k-space data)
  d = dir(fname_tgz);
  fsize = d.bytes/1e9; % in GB
  if fsize>parms.fsize_thresh
    fprintf('%s: WARNING: skipping import for %s because file size is %0.1f GB\n',...
      mfilename,fstem_tgz,fsize);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: skipping import for %s because file size is %0.1f GB\n',...
        mfilename,fstem_tgz,fsize);
    end;
    fstem_tgz = [];
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jinfo,fstem_json,fname_tar,fname_json] = check_json_file(fname_tgz,parms)
  jinfo = []; fstem_json = []; fname_tar = [];

  [indir,fstem_json] = fileparts(fname_tgz);
  fname_json = sprintf('%s/%s.json',parms.indir,fstem_json);

  % check if json file exists and is not empty
  if ~exist(fname_json,'file')
    fprintf('%s: NOTE: json file missing for %s\n',...
      mfilename,fstem_json);
    if parms.flog>0
      fprintf(parms.flog,'%s: NOTE: json file missing for %s\n',...
        mfilename,fstem_json);
    end;    
    ext_flag = 1;
  else
    k = dir(fname_json);
    if k.bytes==0
      fprintf('%s: NOTE: json file %s is empty\n',mfilename,fname_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: NOTE: json file %s is empty\n',mfilename,fname_json);
      end;
      ext_flag = 1;
    else
      ext_flag = 0;
    end;     
  end;
  
  % extract json file from tgz into unpack dir
  if ext_flag
    unpackdir = sprintf('%s/%s',parms.unpackdir,fstem_json);
    mmil_mkdir(unpackdir);
    
    fname_json = sprintf('%s/%s.json',unpackdir,fstem_json);
    if ~exist(fname_json,'file') || parms.forceflag
      fprintf('%s: extracting json file from %s to %s...\n',...
        mfilename,fname_tgz,fname_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: extracting json file from %s...\n',...
          mfilename,fstem_json);
      end;
      errcode = extract_json(fname_tgz,fname_json,parms);
      if errcode
        fprintf('%s: skipping %s because failed to extract json...\n',...
          mfilename,fstem_json);
        if parms.flog>0
          fprintf(parms.flog,'%s: skipping %s because of missing json...\n',...
            mfilename,fstem_json);
        end;
        fname_json = [];
        return; 
      end;
    else
      fprintf('%s: using extracted json file %s...\n',...
        mfilename,fname_json);
      if parms.flog>0
        fprintf(parms.flog,'%s: using extracted json file %s...\n',...
          mfilename,fname_json);
      end;
    end;
  end;

  % check json file, find corresponding tar file
  me = [];
  try
    [jinfo,fstem_json,fname_tar,fname_pfile_json] = ...
      abcd_check_json(fname_json,parms.indir);
  catch me
    jinfo = [];
  end;

  % warn about bad json file
  if isempty(jinfo)
    fprintf('%s: WARNING: skipping import for %s because of error reading json file\n',...
      mfilename,fstem_json);
    if ~isempty(me), fprintf(' %s\n',me.message); end;
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: skipping import for %s because of error reading json file\n',...
        mfilename,fstem_json);
      if ~isempty(me), fprintf(parms.flog,' %s\n',me.message); end;
    end;
  elseif jinfo.invalid_flag
    fprintf('%s: WARNING: skipping import for %s because of invalid json file\n',...
      mfilename,fstem_json);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: skipping import for %s because of invalid json file\n',...
        mfilename,fstem_json);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = extract_json(fname_tgz,fname_json,parms)
  errcode = 0;
  fprintf('%s: checking for json in %s...\n',mfilename,fname_tgz);
  if parms.flog>0
    fprintf(parms.flog,'%s: checking for json in %s...\n',mfilename,fname_tgz);
  end;
  cmd = sprintf('echo `tar -ztf %s | grep json`',fname_tgz);
  [errcode,r] = unix(cmd);
  if errcode
    fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
    end;
    return;
  end;
  if isempty(r)
    fprintf('%s: WARNING: json file not found inside %s\n',mfilename,fname_tgz);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: json file not found inside %s\n',mfilename,fname_tgz);
    end;
    errcode = 1;
    return;
  end;
  fprintf('%s: extracting json from %s...\n',mfilename,fname_tgz);
  if parms.flog>0
    fprintf(parms.flog,'%s: extracting json from %s...\n',mfilename,fname_tgz);
  end;
  fname_tmp = strtrim(r);
  cmd = sprintf('tar -zxf ''%s'' ''%s'' -O > ''%s''',fname_tgz,fname_tmp,fname_json);
  [errcode,r] = unix(cmd);
  if errcode
    fprintf('%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
    if parms.flog>0
      fprintf(parms.flog,'%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mstr = read_md5sum_file(fname,parms)
  mstr = [];
  try
    mstr = textread(fname,'%s');
    mstr = mstr{1};
  catch
    fprintf('%s: WARNING: md5sum file %s could not be opened\n',...
      mfilename,fname);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: md5sum file %s could not be opened\n',...
        mfilename,fname);
    end;
    return;
  end;
  if isempty(mstr)
    fprintf('%s: WARNING: md5sum file %s is empty\n',...
      mfilename,fname);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: md5sum file %s is empty\n',...
        mfilename,fname);
    end;
    return;
  end;
  k = regexp(mstr,'(?<md5sum>\w+)[\W\s]','names');
  if ~isempty(k)
    mstr = k(1).md5sum;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function close_logfile(parms)
  if parms.flog>0
    fclose(parms.flog);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

