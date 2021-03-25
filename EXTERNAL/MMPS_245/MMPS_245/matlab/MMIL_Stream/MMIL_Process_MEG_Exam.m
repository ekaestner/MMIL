function MMIL_Process_MEG_Exam(ContainerPath,varargin)
%function MMIL_Process_MEG_Exam(ContainerPath,[options])
%
% Required Input:
%  ContainerPath:  full path of input Container (orig, raw, or proc MEG)
%
% Optional Parameters
%  'FSContainerPath': full path FreeSurfer recon
%     If supplied, will add to ContainerInfo, enabling MEG_MMIL_Reg2MRI
%     {default = []}
%  'RootDirs':
%    a struct which must contain the following fields:
%         raw_meg, proc_meg
%    and may contain the following fields:
%         orig, raw, proc
%         fsurf, fsico, fsclean, long
%    these specify the locations of data
%    If empty, will place output directories in root dir of ContainerPath
%     {default = []}
%  'ProjID': project name
%     {default = []}
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
%  'PROC_prefix' - prefix of all output files
%    {default: 'proc'}
%  'PROC_valid_event_codes' - vector of event codes to be averaged
%    {default: []} (if empty, treat all event codes as valid)
%  'PROC_max_num_trials' - maximum number of trials per condition
%    {default: Inf} (infinite)
%  'PROC_prestim_dur'  - duration of prestimulus period (msec)
%     {default: 100}
%  'PROC_poststim_dur' - duration of poststimulus period (msec)
%     {default: 400}
%  'PROC_stim_delay' - duration of stimulus onset delay after trigger (msec)
%     {default: 0}
%  'PROC_reject_mag'  - automatic rejection threshold for magnetometer channels (fT)
%     if 0, rejection based on magnetometers is disabled
%     {default: 10000}
%  'PROC_reject_grad' - auto-rejection threshold for gradiometer channels (fT/cm)
%     if 0, rejection based on gradiometers is disabled
%     {default: 6000}
%  'PROC_reject_eeg' - auto-rejection threshold for EEG channels (uV)
%     if 0, rejection based on eeg is disabled
%     {default: 0}
%  'PROC_reject_eog' - auto-rejection threshold for EOG channel (uV)
%     if 0, rejection based on eog is disabled
%     {default: 200}
%  'PROC_bandpass_flag' - [1|0] Toggle bandpass fft filter before averaging
%     {default: 0}
%  'PROC_bandpass_low_cf'  - low center frequency (high-pass filter) (Hz)
%     {default: 0.2}
%  'PROC_bandpass_low_tb'  - low transition band (Hz)
%     {default: 0.4}
%  'PROC_bandpass_high_cf' - high center frequency (low-pass filter) (Hz)
%     {default: 100}
%  'PROC_bandpass_high_tb' - high transition band (Hz)
%     {default: 10}
%  'PROC_notch_flag' - [1|0] Toggle notch fft filter before averaging
%     {default: 0}
%  'PROC_notch_cf'  - notch center frequency (notch filter) (Hz)
%     {default: 60}
%  'PROC_notch_tb'  - notch transition band (Hz)
%     {default: 4}
%  'PROC_dsfact'  - downsampling factor -- must be an integer
%     { default: 1 (no downsampling) }
%  'PROC_detrend_flag' - [1|0] Toggle detrending of single trials
%     before averaging
%     {default: 1}
%  'PROC_baseline_flag' - [1|0] Toggle baseline subtraction of single trials
%     before averaging
%     {default: 1}
%  'PROC_baseline_start' - start time of baseline period (msec)
%     { default: -Inf } (start at beginning of prestimulus period)
%  'PROC_baseline_end'   - end time of baseline period (msec)
%     { default: 0 } (end at trigger onset)
%  'PROC_ncov_ex_evnts' - vector of event codes that should not
%     be used in calculating the noise covariance matrix
%     { default: [] }
%  'PROC_badchans'    - vector of bad channel indices -- will be set to zero
%     { default: [] }
%  'PROC_badchanfile' - name of text file containing bad channel labels
%    {default: 'badchanfile'}
%  'PROC_post_combcond_flag': [0|1] post-averaging combination of conditions
%     {default = 0}
%  'PROC_post_combcond_append_flag': [0|1] for combination of conditions
%     add new conditions to existing, otherwise, only include new conditions
%     {default = 1}
%  'PROC_fname_combcond_info': full path name of csv (comma-separated-value) file
%     containing 'neweventcodes' and 'combinations'
%     {default = []}
%  'PROC_forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% For a complete list of processing options, see ts_process_fif_data
%
% Created:  02/17/11 by Don Hagler
% Last Mod: 08/01/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'FSContainerPath',[],[],...
  'RootDir',pwd,[],...
  'RootDirs',[],[],...
  'ProjID',[],[],...
... % import parameters
  'RAW_rawflag',0,[0:2],...
  'RAW_sssflag',false,[false true],...
  'RAW_moveflag',false,[false true],...
  'RAW_format','float',{'short','long','float'},...
  'RAW_st',6,[0.1 100],...
  'RAW_corr',0.95,[0 1],...
  'RAW_forceflag',false,[false true],...
... % processing parameters
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
});

required_rootdirs = {'orig_meg','raw_meg','proc_meg'};
if isempty(parms.RootDirs)
  parms.RootDirs = MMIL_Set_Common_RootDirs(parms.RootDir,required_rootdirs);
end;
RootDirs = MMIL_Check_RootDirs(parms.RootDirs,required_rootdirs);

[ContainerRootDir,ContainerDir,tmp_ext] = fileparts(ContainerPath);
ContainerDir = [ContainerDir tmp_ext];
if isempty(ContainerRootDir)
  ContainerRootDir = pwd;
  ContainerPath = [ContainerRootDir '/' ContainerDir];
end;

% determine the current stage of processing
procstep = 0;
VisitID = ContainerDir;
ContainerType = 'orig_meg';
n = regexp(ContainerDir,'^(?<ContainerType>[^_]+)_(?<VisitID>\w+)','names');
if ~isempty(n)
  switch n.ContainerType
    case 'MEGRAW'
      procstep = 1;
      ContainerType = 'raw_meg';
      VisitID = n.VisitID;
    case 'MEGPROC'
      procstep = 2;
      ContainerType = 'proc_meg';
      VisitID = n.VisitID;
  end;
end;

% check for badchanfile in ContainerPath
if ~isempty(parms.PROC_badchanfile) && mmil_isrelative(parms.PROC_badchanfile)
  rel_badchanfile = parms.PROC_badchanfile;
  parms.PROC_badchanfile = [ContainerPath '/' rel_badchanfile];
  if exist(parms.PROC_badchanfile,'file')
    fprintf('%s: found badchanfile %s\n',mfilename,parms.PROC_badchanfile);
  else
    parms.PROC_badchanfile = rel_badchanfile;
  end;
end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if procstep==0 % start from orig_meg
  args = MMIL_Args(parms,'MMIL_Import_MEG');
  [ContainerPath,errcode] = MMIL_Import_MEG(ContainerPath,RootDirs,args{:});
  if errcode ~= 0, return; end;  
  procstep = 1;

  % add ProjID and MMPSVER to ContainerInfo
  [errcode,warncode] = MMIL_Stamp_Container(ContainerPath,parms.ProjID);
  if errcode ~= 0, return; end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if procstep==1 % start from raw_meg
  args = MMIL_Args(parms,'MMIL_Process_MEG');
  [ContainerPath,errcode] = MMIL_Process_MEG(ContainerPath,RootDirs,args{:});
  if errcode ~= 0, return; end;  
  procstep = 2;

  % add ProjID and MMPSVER to ContainerInfo
  [errcode,warncode] = MMIL_Stamp_Container(ContainerPath,parms.ProjID);
  if errcode ~= 0, return; end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

