function [ContainerOutPath,errcode] = MMIL_Process_MEG(ContainerPath,RootDirs,varargin)
%function [ContainerOutPath,errcode] = MMIL_Process_MEG(ContainerPath,RootDirs,varargin)
%
% Required Input:
%  ContainerPath: full path of directory containing imported, raw MEG data
%  RootDirs:
%    a struct which must contain the following fields:
%         orig_meg, raw_meg, proc_meg
%    these specify the locations of data
%   orig_meg: full path root dir for orig MEG Containers
%   raw_meg: full path root dir for raw MEG Containers
%   proc_meg: full path root dir for processed MEG Containers
%
% Optional Input (partial list):
%  'FSContainerPath': full path FreeSurfer recon
%    If supplied, will add to ContainerInfo, enabling MEG_MMIL_Reg2MRI
%    {default = []}
%  'prefix': prefix of all output files
%    {default = 'proc'}
%  'valid_event_codes': vector of event codes to be averaged
%    {default = []} (if empty, treat all event codes as valid)
%  'max_num_trials': maximum number of trials per condition
%    {default = Inf} (infinite)
%  'stim_delay': stimulus delay (msec) -- use to compensate for lag between
%    trigger and actual stimulus presentation
%    {default = 0}
%  'prestim_dur': duration of prestimulus period (msec)
%    {default = 100}
%  'poststim_dur': duration of poststimulus period (msec)
%    {default = 400}
%  'stim_delay': duration of stimulus onset delay after trigger (msec)
%    {default = 0}
%  'reject_mag': automatic rejection threshold for magnetometer channels (fT)
%    if 0, rejection based on magnetometers is disabled
%    {default = 10000}
%  'reject_grad': auto-rejection threshold for gradiometer channels (fT/cm)
%    if 0, rejection based on gradiometers is disabled
%    {default = 6000}
%  'reject_eeg': auto-rejection threshold for EEG channels (uV)
%    if 0, rejection based on eeg is disabled
%    {default = 0}
%  'reject_eog': auto-rejection threshold for EOG channel (uV)
%    if 0, rejection based on eog is disabled
%    {default = 200}
%  'bandpass_flag': [1|0] Toggle bandpass fft filter before averaging
%    {default = 0}
%  'bandpass_low_cf': low center frequency (high-pass filter) (Hz)
%    {default = 0.2}
%  'bandpass_low_tb': low transition band (Hz)
%    {default = 0.4}
%  'bandpass_high_cf': high center frequency (low-pass filter) (Hz)
%    {default = 100}
%  'bandpass_high_tb': high transition band (Hz)
%    {default = 10}
%  'notch_flag': [1|0] Toggle notch fft filter before averaging
%    {default = 0}
%  'notch_cf': notch center frequency (notch filter) (Hz)
%    {default = 60}
%  'notch_tb': notch transition band (Hz)
%    {default = 4}
%  'dsfact': downsampling factor -- must be an integer
%    {default = 1 (no downsampling)}
%  'detrend_flag': [1|0] Toggle detrending of single trials
%     before averaging
%    {default = 1}
%  'baseline_flag': [1|0] Toggle baseline subtraction of single trials
%     before averaging
%    {default = 1}
%  'baseline_start': start time of baseline period (msec)
%    {default = -Inf} (start at beginning of prestimulus period)
%  'baseline_end'   - end time of baseline period (msec)
%    {default = 0} (end at trigger onset)
%  'ncov_ex_evnts': vector of event codes that should not
%     be used in calculating the noise covariance matrix
%    {default = []}
%  'badchans'    - vector of bad channel indices -- will be set to zero
%    {default = []}
%  'badchanfile': name of text file containing bad channel labels
%    full path or relative to ContainerPath
%    if badchanfile not found, all channels will be included
%    {default = 'badchans.txt'}
%  'post_combcond_flag': [0|1] post-averaging combination of conditions
%     {default = 0}
%  'post_combcond_append_flag': [0|1] for combination of conditions
%     add new conditions to existing, otherwise, only include new conditions
%     {default = 1}
%  'fname_combcond_info': full path name of csv (comma-separated-value) file
%     containing 'neweventcodes' and 'combinations'
%     {default = []}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% For a complete list of processing options, see ts_process_fif_data
%
% Output:
%   ContainerOutPath: full path of output MEGPROC Container
%   errcode: 0 if no errors, 1 if errors
%
% Created:  02/17/11 by Don Hagler
% Last Mod: 02/19/14 by Don Hagler
%

ContainerOutDir = [];
ContainerOutPath = [];
errcode = 0;

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'FSContainerPath',[],[],...
  'prefix','proc',[],...
  'browseraw',0,[0 Inf],...
  'saveepochs_flag',false,[false true],...
  'write_fif_flag',true,[false true],...
  'fifname_evcode_flag',true,[false true],...
  'trigchan','STI101',[],...
  'events_fnames',[],[],...
  'evcode_offset',[],[],...
  'valid_event_codes',[],[],...
  'event_recode_rules',[],[],...
  'code_excl',[],[],...
  'time_excl_pre',0,[],...
  'time_excl_post',0,[],...
  'trig_minduration',5,[],...
  'max_num_trials',Inf,[1 Inf],...
  'stim_delay',0,[],...
  'prestim_dur',100,[],...
  'poststim_dur',400,[],...
  'reject_mag',10000,[],...
  'reject_grad',6000,[],...
  'reject_eeg',0,[],...
  'reject_eog',200,[],...
  'bandpass_flag',false,[false true],...
  'bandpass_low_cf',0.2,[],...
  'bandpass_low_tb',0.4,[],...
  'bandpass_high_cf',100,[],...
  'bandpass_high_tb',10,[],...
  'notch_flag',false,[false true],...
  'notch_cf',60,[],...
  'notch_tb',4,[],...
  'dsfact',1,[],...
  'detrend_flag',true,[false true],...
  'baseline_flag',true,[false true],...
  'baseline_start',-Inf,[-Inf,Inf],...
  'baseline_end',0,[-Inf,Inf],...
  'ncov_ex_evnts',[],[],...
  'badchans',[],[],...
  'badchanfile','badchans.txt',[],...
  'readtrans_flag',true,[false true],...
  'post_subnull_flag',false,[false true],...
  'null_event',[],[],...
  'post_combcond_flag',false,[false true],...
  'post_combcond_append_flag',false,[false true],...
  'fname_combcond_info',[],[],...
  'post_bandpass_flag',false,[false true],...
  'post_notch_flag',false,[false true],...
  'post_dsfact',1,[],...
  'post_detrend_flag',false,[false true],...
  'post_baseline_flag',false,[false true],...
  'post_stim_delay',0,[],...
  'post_badchans',[],[],...
  'post_badchanfile',[],[],...
  'post_rm_badchans_flag',false,[false true],...
  'ICA_auto_flag',false,[false true],...
  'ICA_manual_flag',false,[false true],...
  'ICA_ref_chan','EOG061',[],...
  'ICA_chantype','all',{'all', 'mag' 'grad1' 'grad2' 'eeg', 'other', 'grad', 'meg'},...
  'ICA_maxsteps',20,[],...
  'ICA_ntrial',5,[],...
  'ICA_ncomponents',80,[],...
  'ICA_rescale_flag',true,[false true],...
  'ICA_sorttrials',false,[false true],...
  'forceflag',false,[false true],...
	'saveperevent',false,[false true],...
  'datatype','single',{'single','double'},...
});

% excl_tags are fields that should not be passed to ts_process_fif_data
excl_tags = {'FSContainerPath'};

fprintf('%s: processing %s...\n',mfilename,ContainerPath);

[tmp,ContainerDir,text] = fileparts(ContainerPath);
if ~isempty(text)
  ContainerDir = [ContainerDir text];
end;

[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);

% check input files
input_data_files = ContainerInfo.output_data_files;
nfiles = length(input_data_files);
if ~nfiles
  fprintf('%s: ERROR: no input fif files in %s\n',...
    mfilename,ContainerPath);
  errcode = 1;
  return;
end;
for f=1:length(nfiles)
  if ~exist(input_data_files{f},'file')
    fprintf('%s: ERROR: input fif files %s not found\n',...
      mfilename,input_data_files{f});
    errcode = 1;
    return;
  end;
end;

% create output directory
ContainerOutDir = regexprep(ContainerDir,'MEGRAW_','MEGPROC_');
ContainerOutPath = [RootDirs.proc_meg '/' ContainerOutDir];
[succ,msg] = mkdir(ContainerOutPath);

% check for badchanfile in MEGPROC container
if ~isempty(parms.badchanfile)
  if mmil_isrelative(parms.badchanfile)
    parms.badchanfile = [ContainerOutPath '/' parms.badchanfile];
  else
    % if badchanfile is somewhere else, copy to MEGPROC container
    [tpath,tstem,text] = fileparts(parms.badchanfile);
    if ~strcmp(tpath,ContainerOutPath)
      proc_badchanfile = [ContainerOutPath '/' tstem text];
      if ~exist(proc_badchanfile,'file')
        cmd = sprintf('cp %s %s',parms.badchanfile,proc_badchanfile);
        [s,r] = unix(cmd);
        if s
          error('cmd %s failed:\n%s',cmd,r);
        end;
      end;
      parms.badchanfile = proc_badchanfile;
    end;
  end;
  if ~exist(parms.badchanfile,'file')
    fprintf('%s: did not find badchanfile %s\n',mfilename,parms.badchanfile);
    parms.badchanfile = [];
  else
    fprintf('%s: found badchanfile %s\n',mfilename,parms.badchanfile);
  end;
end;  

% create ContainerInfo.mat
ContainerOutInfo = [];
ContainerOutInfo.ContainerType = 'MEGPROC';
ContainerOutInfo.VisitID = ContainerInfo.VisitID;
ContainerOutInfo.input_data_files = input_data_files;
ContainerOutInfo.FSContainerPath = parms.FSContainerPath;
errcode = MMIL_Save_ContainerInfo(ContainerOutPath,ContainerOutInfo);

% run ts_process_fif_data for all files
parms.rootoutdir = ContainerOutPath;
tags = setdiff(fieldnames(parms),excl_tags);
args = mmil_parms2args(parms,tags);
ts_process_fif_data(input_data_files,args{:});

