% auto_run_process_data

% automatically generated on 13-Feb-2008

%%%%%%%%%%%%%%%%%%%%%%%%%
% set variables
rawdatadir = '/home/halgdev/incoming/MEG_UCSD/DJH_EG01_070519';
datafile={};datafile{end+1} = '/home/halgdev/incoming/MEG_UCSD/DJH_EG01_070519/eg01_070519_onlineavg.fif';
datafile{end+1} = '/home/halgdev/incoming/MEG_UCSD/DJH_EG01_070519/eg01_070519_raw.fif';

badchanfile = [];
% specify a text file listing the bad channels
valid_event_codes = [];
% specify event codes to be averaged
%   (empty means include all events found in files)
event_recode_rules = [];
% specify rules to change event condition values
%   (empty means include all events found in files)
code_excl = [];
% event code around which other events are excluded
%   (empty means do not exclude events like this)
time_excl_pre = 500;
% duration before code_excl that other events are excluded
time_excl_post = 500;
% duration after code_excl that other events are excluded
trig_minduration = 5;
% minimum duration for data to appear on a trigger channel
% before the trigger is considered to be "on"
stim_delay = 32.5;
% delay between trigger and stimulus onset 
prestim_dur = 100;
% pre-stimulus duration (msec)
poststim_dur = 400;
% post-stimulus duration (msec)
bandpass_low_cf = 0.2;
% bandpass filter low cut-off frequency (Hz)
bandpass_low_tb = 0.4;
% bandpass filter low transition band (Hz)
bandpass_high_cf = 50;
% bandpass filter high cut-off frequency (Hz)
bandpass_high_tb = 10;
% bandpass filter high transition band (Hz)
dsfact = 4;
% downsampling factor
baseline_start = -80;
% start of baseline period (msec) relative to stim onset
baseline_end = -5;
% end of baseline period (msec) relative to stim onset
null_event = [];
% event code used for "null" subtraction
post_dsfact = 1;
% downsampling factor applied after averaging

bandpass_flag = 1;
% whether (0 or 1) to apply bandpass filtering before averaging
detrend_flag = 1;
% whether (0 or 1) to detrend before averaging
baseline_flag = 1;
% whether (0 or 1) to baseline subtract before averaging
post_subnull_flag = 0;
% whether (0 or 1) to do null subtraction after averaging
post_bandpass_flag = 0;
% whether (0 or 1) to apply bandpass filtering after averaging
post_detrend_flag = 0;
% whether (0 or 1) to detrend after averaging
post_baseline_flag = 0;
% whether (0 or 1) to baseline subtract after averaging

saveepochs_flag = 0;
% whether (0 or 1) to save trial-by-trial data

%%%%%%%%%%%%%%%%%%%%%%%%%
% run ts_process_fif_data
ts_process_fif_data(...
  datafile,...
  'browseraw',2,...
  'badchanfile',badchanfile,...
  'valid_event_codes',valid_event_codes,...
  'event_recode_rules',event_recode_rules,...
  'trig_minduration',trig_minduration,...
  'stim_delay',stim_delay,...
  'prestim_dur',prestim_dur,...
  'poststim_dur',poststim_dur,...
  'bandpass_low_cf',bandpass_low_cf,...
  'bandpass_low_tb',bandpass_low_tb,...
  'bandpass_high_cf',bandpass_high_cf,...
  'bandpass_high_tb',bandpass_high_tb,...
  'dsfact',dsfact,...
  'baseline_start',baseline_start,...
  'baseline_end',baseline_end,...
  'noise_start',baseline_start,...
  'noise_end',baseline_end,...
  'null_event',null_event,...
  'post_dsfact',post_dsfact,...
  'code_excl',code_excl,...
  'time_excl_pre',time_excl_pre,...
  'time_excl_post',time_excl_post,...
  'bandpass_flag',bandpass_flag,...
  'detrend_flag',detrend_flag,...
  'baseline_flag',baseline_flag,...
  'post_subnull_flag',post_subnull_flag,...
  'post_bandpass_flag',post_bandpass_flag,...
  'post_detrend_flag',post_detrend_flag,...
  'post_baseline_flag',post_baseline_flag,...
  'saveepochs_flag',saveepochs_flag);
