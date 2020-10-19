function MMIL_Process_MEG_Exams(ProjID,varargin)
% function MMIL_Process_MEG_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    Use this option to limit processing to specific subjects
%    Should have a field 'MEG_VisitID' that indicates the name of the
%       input data directory (i.e. in orig_meg data directory)
%    If MEG_VisitID field is missing, will use SubjID field instead
%    May also specify subject-specific parameter values
%     but note that these will override ProjInfo and command line input
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: raw_meg, proc_meg
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'batchname': name of ougtput batchdir
%    {default = 'MMIL_Process_MEG_Exams'}
%  'procstep': [0|1] at which stage to start processing
%    0 = start from orig_meg, 1 = start from raw_meg
%    {default = 0}
%  'newflag': [0|1] whether to create jobs only for new data
%    ignored if forceflag=1
%    {default = 0}
%  'qcflag': [0|1] whether to exclude visits with StudyInfo.QC = 0
%    {default = 1}
%  'forceflag': [0|1] whether to overwrite existing output
%    If 1, will override RAW_forceflag and PROC_forceflag
%    {default = 0}
%
% Optional Import Parameters:
%  'RAW_rawflag': [0|1|2] link (1) or copy (2) data from orig to raw
%    ignored if RAW_sssflag = 1 or RAW_moveflag = 1
%    {default = 0}
%  'RAW_sssflag': [0|1] whether to run maxfilter with sss
%    {default = 0}
%  'RAW_moveflag': [0|1] whether to run maxfilter with maxmove
%    {default = 0}
%  'RAW_format': output data format (short, long, or float)
%     {default = 'float'}
%  'RAW_forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Optional Processing Parameters:
%  'PROC_prefix': prefix of all output files
%    {default: 'proc'}
%  'PROC_valid_event_codes': vector of event codes to be averaged
%    {default: []} (if empty, treat all event codes as valid)
%  'PROC_max_num_trials': maximum number of trials per condition
%    {default: Inf} (infinite)
%  'PROC_prestim_dur': duration of prestimulus period (msec)
%     {default: 100}
%  'PROC_poststim_dur': duration of poststimulus period (msec)
%     {default: 400}
%  'PROC_stim_delay': duration of stimulus onset delay after trigger (msec)
%     {default: 0}
%  'PROC_reject_mag': automatic rejection threshold for magnetometer channels (fT)
%     if 0, rejection based on magnetometers is disabled
%     {default: 10000}
%  'PROC_reject_grad': auto-rejection threshold for gradiometer channels (fT/cm)
%     if 0, rejection based on gradiometers is disabled
%     {default: 6000}
%  'PROC_reject_eeg': auto-rejection threshold for EEG channels (uV)
%     if 0, rejection based on eeg is disabled
%     {default: 0}
%  'PROC_reject_eog': auto-rejection threshold for EOG channel (uV)
%     if 0, rejection based on eog is disabled
%     {default: 200}
%  'PROC_bandpass_flag': [1|0] Toggle bandpass fft filter before averaging
%     {default: 0}
%  'PROC_bandpass_low_cf': low center frequency (high-pass filter) (Hz)
%     {default: 0.2}
%  'PROC_bandpass_low_tb': low transition band (Hz)
%     {default: 0.4}
%  'PROC_bandpass_high_cf': high center frequency (low-pass filter) (Hz)
%     {default: 100}
%  'PROC_bandpass_high_tb': high transition band (Hz)
%     {default: 10}
%  'PROC_notch_flag': [1|0] Toggle notch fft filter before averaging
%     {default: 0}
%  'PROC_notch_cf': notch center frequency (notch filter) (Hz)
%     {default: 60}
%  'PROC_notch_tb': notch transition band (Hz)
%     {default: 4}
%  'PROC_dsfact': downsampling factor -- must be an integer
%     { default: 1 (no downsampling) }
%  'PROC_detrend_flag': [1|0] Toggle detrending of single trials
%     before averaging
%     {default: 1}
%  'PROC_baseline_flag': [1|0] Toggle baseline subtraction of single trials
%     before averaging
%     {default: 1}
%  'PROC_baseline_start': start time of baseline period (msec)
%     { default: -Inf } (start at beginning of prestimulus period)
%  'PROC_baseline_end': end time of baseline period (msec)
%     { default: 0 } (end at trigger onset)
%  'PROC_ncov_ex_evnts': vector of event codes that should not
%     be used in calculating the noise covariance matrix
%     { default: [] }
%  'PROC_badchans': vector of bad channel indices -- will be set to zero
%     { default: [] }
%  'PROC_badchanfile': name of text file containing bad channel labels
%    {default: 'badchans.txt'}
%  'PROC_post_combcond_flag': [0|1] post-averaging combination of conditions
%     {default = 0}
%  'PROC_post_combcond_append_flag': [0|1] for combination of conditions
%     add new conditions to existing, otherwise, only include new conditions
%     {default = 1}
%  'PROC_fname_combcond_info': name of csv (comma-separated-value) file
%     containing 'neweventcodes' and 'combinations'
%     May be full path or relative to /home/user/ProjInfo/ProjID or orig MEG container
%     {default = []}
%  'PROC_forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  02/16/11 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'procstep',0,[0:1],...
  'newflag',false,[false true],...
  'qcflag',false,[false true],...
  'forceflag',false,[false true],...
  'batchname','MMIL_Process_MEG_Exams',[],...
...
  'RAW_rawflag',0,[0:2],...
  'RAW_sssflag',false,[false true],...
  'RAW_moveflag',false,[false true],...
  'RAW_format','float',{'short','long','float'},...
  'RAW_st',6,[0.1 100],...
  'RAW_corr',0.95,[0 1],...
  'RAW_forceflag',false,[false true],...
...
  'PROC_prefix','proc',[],...
  'PROC_browseraw',0,[0 Inf],...
  'PROC_saveepochs_flag',false,[false true],...
  'PROC_write_fif_flag',true,[false true],...
  'PROC_fifname_evcode_flag',true,[false true],...
  'PROC_trigchan','STI101',[],...
  'PROC_events_fnames',[],[],...
  'PROC_evcode_offset',[],[],...
  'PROC_valid_event_codes',[],[],...
  'PROC_event_recode_rules',[],[],...
  'PROC_code_excl',[],[],...
  'PROC_time_excl_pre',0,[],...
  'PROC_time_excl_post',0,[],...
  'PROC_trig_minduration',5,[],...
  'PROC_max_num_trials',Inf,[1 Inf],...
  'PROC_stim_delay',0,[],...
  'PROC_prestim_dur',100,[],...
  'PROC_poststim_dur',400,[],...
  'PROC_reject_mag',10000,[],...
  'PROC_reject_grad',6000,[],...
  'PROC_reject_eeg',0,[],...
  'PROC_reject_eog',200,[],...
  'PROC_bandpass_flag',false,[false true],...
  'PROC_bandpass_low_cf',0.2,[],...
  'PROC_bandpass_low_tb',0.4,[],...
  'PROC_bandpass_high_cf',100,[],...
  'PROC_bandpass_high_tb',10,[],...
  'PROC_notch_flag',false,[false true],...
  'PROC_notch_cf',60,[],...
  'PROC_notch_tb',4,[],...
  'PROC_dsfact',1,[],...
  'PROC_detrend_flag',true,[false true],...
  'PROC_baseline_flag',true,[false true],...
  'PROC_baseline_start',-Inf,[-Inf,Inf],...
  'PROC_baseline_end',0,[-Inf,Inf],...
  'PROC_ncov_ex_evnts',[],[],...
  'PROC_badchans',[],[],...
  'PROC_badchanfile','badchans.txt',[],...
  'PROC_readtrans_flag',true,[false true],...
  'PROC_post_subnull_flag',false,[false true],...
  'PROC_null_event',[],[],...
  'PROC_post_combcond_flag',false,[false true],...
  'PROC_post_combcond_append_flag',false,[false true],...
  'PROC_fname_combcond_info',[],[],...
  'PROC_post_bandpass_flag',false,[false true],...
  'PROC_post_notch_flag',false,[false true],...
  'PROC_post_dsfact',1,[],...
  'PROC_post_detrend_flag',false,[false true],...
  'PROC_post_baseline_flag',false,[false true],...
  'PROC_post_stim_delay',0,[],...
  'PROC_post_badchans',[],[],...
  'PROC_post_badchanfile',[],[],...
  'PROC_post_rm_badchans_flag',false,[false true],...
  'PROC_ICA_auto_flag',false,[false true],...
  'PROC_ICA_manual_flag',false,[false true],...
  'PROC_ICA_ref_chan','EOG061',[],...
  'PROC_ICA_chantype','all',{'all', 'mag' 'grad1' 'grad2' 'eeg', 'other', 'grad', 'meg'},...
  'PROC_ICA_maxsteps',20,[],...
  'PROC_ICA_ntrial',5,[],...
  'PROC_ICA_ncomponents',80,[],...
  'PROC_ICA_rescale_flag',true,[false true],...
  'PROC_ICA_sorttrials',false,[false true],...
  'PROC_saveperevent',false,[false true],...
  'PROC_datatype','single',{'single','double'},...
  'PROC_forceflag',false,[false true],...
...
  'FSContainerPath',[],[],...
...
  'required_rootdirs',{'orig_meg','raw_meg','proc_meg'},[],...
};
parms = mmil_args2parms(varargin,parms_filter);

% excl_tags are fields that should not be passed to MMIL_Process_MEG_Exam
excl_tags = {'procstep','newflag','qcflag','forceflag',...
  'RootDirs','StudyInfo','batchname','required_rootdirs'};
tags = setdiff(fieldnames(parms),excl_tags);

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

if parms.forceflag
  parms.newflag = 0;
  parms.RAW_forceflag = 1;
  parms.PROC_forceflag = 1;
end;

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

switch parms.procstep
  case 0
    RootDir = RootDirs.orig_meg;
    regpat = '(?<MEG_VisitID>^[^\.].+)';
  case 1
    RootDir = RootDirs.raw_meg;
    regpat = '^MEGRAW_(?<MEG_VisitID>^[^\.].+)';
end;
dirlist = dir(sprintf('%s/*',RootDir));

j = 1;
for i=1:length(dirlist)
  ContainerDir = char(dirlist(i).name);
  n = regexp(ContainerDir,regpat,'names','once');
  if isempty(n) || ~dirlist(i).isdir
    continue;
  end;
  if parms.procstep==0
    % replace any '.' in MEG_VisitID with '_'
    %  (causes problems when inserted in job names and for ADNI)
    MEG_VisitID = regexprep(ContainerDir,'\.','_');
    % skip if this is a raw_meg or proc_meg container in orig RootDir (e.g. shared root dir)
    if ~isempty(regexp(MEG_VisitID,'MEGRAW')) ||...
       ~isempty(regexp(MEG_VisitID,'MEGPROC'))
      continue;
    end;    
    StudyDate = [];
  else
    MEG_VisitID = n.MEG_VisitID;
    StudyDate = n.StudyDate;
  end;
  ContainerPath = [RootDir '/' ContainerDir];

  tmp_parms = parms;
  if ~isempty(StudyInfo)
    % check that MEG_VisitID is in StudyInfo  
    ind = find(strcmp(MEG_VisitID,{StudyInfo.MEG_VisitID}));
    if isempty(ind)
      fprintf('%s: WARNING: MEG_VisitID %s not found in StudyInfo... skipping\n',...
        mfilename,MEG_VisitID);
      continue;
    end;
    if length(ind)>1
      fprintf('%s: WARNING: MEG_VisitID %s is found %d times in StudyInfo... using first entry only\n',...
        mfilename,MEG_VisitID,length(ind));
      ind = ind(1);
    end;

    % replace values in parms with (non-empty) values from StudyInfo
    for t=1:length(tags)
      tmp_val = mmil_getfield(StudyInfo(ind),tags{t},[]);
      if ~isempty(tmp_val), tmp_parms.(tags{t}) = tmp_val; end;
    end;

    % FSContainerPath is not used for MEG processing, but is added
    %   to ContainerInfo to enable MEG_MMIL_Reg2MRI
    if isfield(RootDirs,'fsurf') && ~isempty(StudyInfo(ind).fsurf)
      tmp_parms.FSContainerPath = sprintf('%s/%s',...
        RootDirs.fsurf,StudyInfo(ind).fsurf);
    else
      tmp_parms.FSContainerPath = [];
    end;
  end;

  % PROC_fname_combcond_info file could be full path, or relative to ProjInfo or ContainerPath
  if tmp_parms.PROC_post_combcond_flag && ~isempty(tmp_parms.PROC_fname_combcond_info)
    if mmil_isrelative(tmp_parms.PROC_fname_combcond_info)
      tmp_fname = [ContainerPath '/' tmp_parms.PROC_fname_combcond_info];
      if exist(tmp_fname,'file')
        tmp_parms.PROC_fname_combcond_info = tmp_fname;
      else
        tmp_fname = [RootDirs.home '/ProjInfo/' ProjID '/'...
          tmp_parms.PROC_fname_combcond_info];
        if exist(tmp_fname,'file')
          tmp_parms.PROC_fname_combcond_info = tmp_fname;
        else
          fprintf('%s: WARNING: missing PROC_fname_combcond_info file %s for %s\n',...
            mfilename,tmp_parms.PROC_fname_combcond_info,MEG_VisitID);
          tmp_parms.PROC_post_combcond_flag = 0;
        end;
      end;
    end;
  end;

  % skip a session if output already exists
  if parms.newflag && parms.procstep<2
    [tmp,ContainerOut] = ...
      MMIL_Get_Container(RootDirs,MEG_VisitID,'proc_meg');
    if ~isempty(ContainerOut)
      fprintf('%s: WARNING: %s already exists... skipping\n',...
        mfilename,ContainerOut);
      continue;
    end;
  end;

  % create script
  jstem = regexprep(MEG_VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  matfname = [batchdir '/' jobID '.mat'];
  save(matfname,'RootDirs');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_Process_MEG_Exam(''%s'',...\n',ContainerPath);
  fprintf(fid,'  ''RootDirs'',RootDirs,...\n');
  mmil_write_tags(fid,tags,tmp_parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

