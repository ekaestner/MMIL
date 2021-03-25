function ABCD_Compile_Exams(ProjID,varargin)
%function ABCD_Compile_Exams(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional parameters:
%   'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%   'batchname': name of output batchdir
%     {default = 'ABCD_Compile_Exams'}
%   'source_user': user ID for input data and QC/PC info
%     {default: 'mmilrec14'}
%   'source_ProjID': project ID string for input data
%     {default = 'DAL_ABCD_QC'}
%   'rootdir': root directory containing project directories with orig dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%   'max_exams': maximum number of exams to compile
%     {default = Inf}
%   'time_limit': number of days after which a visit will be compiled
%     {default = 30}
%   'linkflag': [0|1] create symbolic links to directories instead of copying
%     {default = 0}
%   'follow_links_flag': [0|1] when copying directories, turn links into copies
%     ignored if linkflag = 1
%     {default = 0}
%   'ndar_flag': [0|1] only include visits with pGUIDs that start with 'NDAR_'
%     {default = 1}
%   'compliant_flag': [0|1] only include series that are ABCD protocol compliant
%     {default = 0}
%   'complete_flag': [0|1] only include series that are complete
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%
% Created:  11/20/16 by Don Hagler
% Last Mod: 05/31/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: document instem_info
%% todo: make QC filter optional

%% todo: check that a visit is complete based on series types
%%       and pre-set rules about what is a complete visit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(ProjID,varargin);

% create output batch directory
parms = create_batchdir(parms);

% create scriptlist
parms = create_scriptlist(parms);

% create log file and display starting message
parms = create_log(parms);

% initialize series info
series_info = get_series_info(parms);

notes_fields = get_notes_fields(series_info);

% loop over subjects, creating jobs to compile exams
create_jobs(series_info,parms,notes_fields);

% close log file
if parms.flog>0, fclose(parms.flog); end;

% check available disk space
MMIL_Check_Usage(parms.rootdir);

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms = mmil_args2parms(options,{...
    'ProjID',ProjID,[],...
  ...
    'batchrootdir',[],[],...
    'batchname','ABCD_Compile_Exams',[],...
  ...
    'source_user','mmilrec14',[],...
    'source_ProjID','DAL_ABCD_QC',[],...
    'rootdir','/space/syn05/1/data/MMILDB',[],...
    'max_exams',Inf,[1,Inf],...
    'time_limit',30,[1,Inf],...
    'linkflag',false,[false true],...
    'follow_links_flag',false,[false true],...
    'ndar_flag',true,[false true],...
    'compliant_flag',false,[false true],...
    'complete_flag',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'instem_info','merged_pcqcinfo',[],...
    'root_outdir',[],[],...
    'fname_info',[],[],...
    'logfile',[],[],...
    'SubjID_pattern','[GSPX]\d{3}_(?<SubjID>[^_]+)_\d{8}',[],...
  ...
    'tags',{'ProjID','root_outdir',...
            'linkflag','follow_links_flag','forceflag'},[],...
  });
  if ~isempty(ProjID)
    parms.batchname = [ProjID '_' parms.batchname];
  end;
  if isempty(parms.root_outdir)
    parms.root_outdir = sprintf('%s/%s',parms.rootdir,ProjID);
  end;
  if isempty(parms.fname_info)
    parms.fname_info = sprintf('/home/%s/MetaData/%s/%s_%s.csv',...
      parms.source_user,parms.source_ProjID,...
      parms.source_ProjID,parms.instem_info);
  end;
  if ~exist(parms.fname_info,'file')
    error('info file %s not found',parms.fname_info);
  end;
  if isempty(parms.logfile)
    parms.logfile = sprintf('%s/logs/%s_compile_exams.log',...
      getenv('HOME'),ProjID);
  end;
  % attempt to create output directory
  %   (opportunity to error if lacking permission)
  mmil_mkdir(parms.root_outdir);
  % set output directory for abcd_load_csv
  parms.cachedir = sprintf('/home/%s/MetaData/%s/cache',...
      getenv('USER'),parms.ProjID);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_batchdir(parms)
  if isempty(parms.batchrootdir)
    parms.batchrootdir = sprintf('%s/batchdirs',getenv('HOME'));
  end;
  mmil_mkdir(parms.batchrootdir);
  parms.batchdir = sprintf('%s/%s',parms.batchrootdir,parms.batchname);
  if exist(parms.batchdir,'dir')
    cmd = sprintf('rm -rf %s\n',parms.batchdir);
    fprintf('cmd = %s',cmd);
    [s,r] = unix(cmd);
    if s
      fprintf('%s: WARNING: cmd %s failed:\n%s',mfilename,cmd,r);
    end;
  end;
  mmil_mkdir(parms.batchdir)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_scriptlist(parms)
  parms.scriptlistfname = sprintf('%s/scriptlist.txt',parms.batchdir);
  fid = fopen(parms.scriptlistfname,'w');
  if fid==-1
    error('failed to open scriptlist file %s for writing',parms.scriptlistfname);
  end;
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_log(parms)
  % create logfile
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
  fprintf('%s: starting abcd compile visits at %s...\n',...
    mfilename,datestr(now));
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  if parms.flog>0
    fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(parms.flog,'%s: starting abcd compile visits at %s...\n',...
      mfilename,datestr(now));
    fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = get_series_info(parms)
  fprintf('%s: reading %s...\n',mfilename,parms.fname_info);
  if parms.flog>0
    fprintf(parms.flog,'%s: reading %s...\n',mfilename,parms.fname_info);
  end;
  % read PC and QC info
%  series_info = mmil_csv2struct(parms.fname_info);
  %series_info = abcd_load_csv(parms.fname_info,parms.cachedir);
  series_info = abcd_load_csv(parms.fname_info);
  % exclude empty pGUID
  i_valid = find(~cellfun(@isempty,{series_info.pGUID}) |...
                 ~cellfun(@isnumeric,{series_info.pGUID}));
  series_info = series_info(i_valid);
  % replace 'invalid' pGUID with SubjID inferred from VisitID
  for i=1:length(series_info)
    pGUID = series_info(i).pGUID;
    if strcmp(pGUID,'invalid')
      VisitID = series_info(i).VisitID;
      n = regexp(VisitID,parms.SubjID_pattern,'names');
      if ~isempty(n)
        series_info(i).pGUID = n.SubjID;
      else
        fprintf('%s: WARNING: VisitID %s does not match SubjID_pattern\n',...
          mfilename,VisitID);
      end;
    end;
  end;
  % exclude invalid pGUID
  i_valid = find(cellfun(@isempty,regexp({series_info.pGUID},'invalid')));
  series_info = series_info(i_valid);
  % exclude non-NDAR pGUID
  if parms.ndar_flag
    i_valid = find(~cellfun(@isempty,regexp({series_info.pGUID},'^NDAR')));
    series_info = series_info(i_valid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_jobs(series_info,parms,notes_fields)
  RootDirs = [];
  RootDirs.orig = [parms.root_outdir '/orig'];
  RootDirs.raw = [parms.root_outdir '/raw'];

  % loop over unique subject IDs
  pGUIDs = unique({series_info.pGUID});
  nexams = 0;
  j = 1;
  switchoff = 0;
  for i=1:length(pGUIDs)
    pGUID = pGUIDs{i};
    %if ~strcmp(pGUID,'NDAR_INV07RPB2TU') && ~switchoff
    %  switchoff = 1;
    %  continue;
    %end;
    % identify all series for this subject
    subj_info = select_subj_series(pGUID,series_info,parms);
    if isempty(subj_info), continue; end;
    % loop over events
    EventNames = {subj_info.EventName};
    uniq_EventNames = unique(EventNames);
    for k=1:length(uniq_EventNames)
      if nexams >= parms.max_exams
        fprintf('%s: NOTE: created jobs for %d exams (max_exams = %d)\n',...
          mfilename,nexams,parms.max_exams);
        if parms.flog>0
          fprintf(parms.flog,'%s: NOTE: created jobs for %d exams (max_exams = %d)\n',...
            mfilename,nexams,parms.max_exams);
        end;
        return;
      end;
      EventName = uniq_EventNames{k};
      i_event = find(strcmp(EventNames,EventName));
      event_info = subj_info(i_event);
      % calculate days since each date of acquisition
      time_delay = calc_time_delay(event_info);

      %% todo: check if all necessary series are included
      %%       or at least one of T1, T2, DTI and at least one task fMRI
      
      % copy link to dicoms from source orig dirs to dest orig dir
      if time_delay > parms.time_limit
        % set compiled VisitID
        VisitID = set_VisitID(event_info,parms);
        % check if output already exists
        if ~parms.forceflag
          OrigContPath = MMIL_Get_Container(RootDirs,VisitID,'orig');
          RawContPath = MMIL_Get_Container(RootDirs,VisitID,'raw');
          if ~isempty(OrigContPath) && ~isempty(RawContPath)
            fname_orig_check = sprintf('%s/contents.mat',OrigContPath);
            fname_raw_check = sprintf('%s/contents.mat',RawContPath);
            if exist(fname_orig_check,'file') &&...
               exist(fname_raw_check,'file')
              fprintf('%s: NOTE: skipping VisitID %s because already compiled\n',...
                mfilename,VisitID);
              if parms.flog>0
                fprintf(parms.flog,'%s: NOTE: skipping VisitID %s because already compiled\n',...
                  mfilename,VisitID);
              end;
              continue;
            end;
          end;
        end;
        %%%%%%%%%%%%%%%%%%%%%
        %Filter T1/T2 here instead
        event_info = clean_subj_struct_series(event_info,notes_fields,'T1');
        event_info = clean_subj_struct_series(event_info,notes_fields,'T2');
        if isempty(event_info), continue; end;
        VisitIDs = {event_info.VisitID};
        SeUIDs = {event_info.SeriesInstanceUID};
        root_indirs = regexprep({event_info.fname_pc_json},'/pc/.*','');
        % create job
        jstem = VisitID;
        if length(jstem)>30, jstem = jstem(1:30); end;
        jstem = regexprep(jstem,'[\^\-]','_');
        fprintf('%s: creating job for %s...\n',mfilename,jstem);
        if parms.flog>0
          fprintf(parms.flog,'%s: creating job for %s...\n',mfilename,jstem);
        end;
        jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
        jobfname = sprintf('%s/%s.m',parms.batchdir,jobID);
        mmil_write_script(jobfname,'ABCD_Compile_Exam',{SeUIDs,VisitIDs,VisitID,root_indirs},...
                          parms.tags,parms);
        % add to list
        fid = fopen(parms.scriptlistfname,'a');
        fprintf(fid,'%s\n',jobID);
        fclose(fid);
        nexams = nexams + 1;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subj_info = select_subj_series(pGUID,series_info,parms)
  subj_info = [];
  pGUIDs = {series_info.pGUID};
  i_subj = find(strcmp(pGUIDs,pGUID));
  subj_info = series_info(i_subj);
  PC = {subj_info.ABCD_Compliant};
  CC = {subj_info.Completed};
  QC = {subj_info.QC};
  % skip subject if series have not been checked yet
  if parms.compliant_flag
    i_reviewed = find(~cellfun(@isempty,PC) & ~cellfun(@isempty,QC));
  else
    i_reviewed = find(~cellfun(@isempty,QC));
  end;
  if isempty(i_reviewed)
    fprintf('%s: WARNING: no reviewed series for pGUID %s\n',mfilename,pGUID);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no reviewed series for pGUID %s\n',mfilename,pGUID);
    end;
    subj_info = [];
    return;
  end;
  % skip subject if no good series
  i_good = find([QC{i_reviewed}]);
  if parms.compliant_flag
    i_good = intersect(i_good,find(~cellfun(@isempty,regexp(PC(i_reviewed),'Yes'))));
  end;
  if parms.complete_flag
    i_good = intersect(i_good,find([CC{i_reviewed}]));
  end;
  if isempty(i_good)
    fprintf('%s: WARNING: no good series for pGUID %s\n',mfilename,pGUID);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no good series for pGUID %s\n',mfilename,pGUID);
    end;
    subj_info = [];
    return;
  end;
  subj_info = subj_info(i_reviewed(i_good));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subj_info = clean_subj_struct_series(subj_info,notes_fields,SType);
  SeriesType = {subj_info.SeriesType};
  idx_norm = strcmp(SeriesType,sprintf('%s_NORM',SType));
  idx_raw = strcmp(SeriesType,SType);
  len_norm = length(find(idx_norm));
  len_raw = length(find(idx_raw));
  switch len_norm
    case 0
      if len_raw < 2
        subj_info = subj_info;
      else
        subj_info = remove_junk_series_by_notes(subj_info,idx_raw,notes_fields);
      end;
    case 1
      subj_info = subj_info(~idx_raw);
    otherwise
      subj_info = subj_info(~idx_raw);
      SeriesType = {subj_info.SeriesType};
      idx_norm = strcmp(SeriesType,sprintf('%s_NORM',SType));
      subj_info = remove_junk_series_by_notes(subj_info,idx_norm,notes_fields);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subj_info = remove_junk_series_by_notes(subj_info,idx,notes_fields)
  subj_info_tmp = subj_info(idx);
  hasnotes=[];
  for i=1:length(notes_fields)
    hasnotes=[hasnotes; cellfun(@isempty, {subj_info_tmp.(cell2mat(notes_fields(i)))})];
  end
  hasnotes = sum(hasnotes);
  least_hasnotes = min(hasnotes);
  least_hasnotes_idx = find(hasnotes==least_hasnotes);
  least_hasnotes_idx = least_hasnotes_idx(end);
  realid = find(idx);
  junk_realidx = setdiff(realid,realid(least_hasnotes_idx));
  good_realidx = realid(least_hasnotes_idx);
  idx = ~idx;
  idx(good_realidx) = 1;
  subj_info = subj_info(find(idx==1));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function time_delay = calc_time_delay(event_info)
  time_delay = Inf;
  StudyDates = {event_info.StudyDate};
  StudyDates = cellfun(@num2str,StudyDates,'UniformOutput',false);
  dt = floor(now - datenum(StudyDates,'yyyymmdd'));
  time_delay = min(dt);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitID = set_VisitID(event_info,parms)
  VisitID = event_info(1).VisitID;
  n = regexp(VisitID,'(?<Site>[^_]+)_(?<Subj>[^_]+)_(?<Date>[^_]+)_','names');
  if isempty(n)
    fprintf('%s: WARNING: series VisitID %s does not match expected pattern\n',...
      mfilename,VisitID);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: series VisitID %s does not match expected pattern\n',...
        mfilename,VisitID);
    end;
    return;
  end;
  VisitID = sprintf('%s_%s_%s',n.Site,n.Subj,n.Date);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function notes_fields = get_notes_fields(series_info)
  fieldnames_list = fieldnames(series_info);
  idx = cellfun(@isempty, regexp(fieldnames_list,'notes_'));
  notes_fields = fieldnames_list(~idx);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
