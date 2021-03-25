function filenames = ts_process_fif_data(datafile,varargin)
%function [filenames] = ts_process_fif_data(datafile,[options])
%
% Purpose: process neuromag MEG/EEG data
%
% Usage:
%  ts_process_fif_data(datafile,'key1', value1,...);
%  e.g. ts_process_fif_data('/usr/data/test.fif','valid_event_codes',[1:6]);
%
% Required Input:
%  datafile: full or relative path name of Neuromag fif raw data file
%    For multiple files that are to be treated as a continuous acquisition,
%    input a cell array of file names e.g.:
%      datafile{1}='/usr/data/test.fif';
%      datafile{2}='/usr/data/test-1.fif';
%      datafile{3}='/usr/data/test-2.fif';
%      ts_process_fif_data(datafile,'valid_event_codes',[1:6]);
%
% Optional parameters:
%  'browseraw': [0,1,...] launch browser for one of the raw data files
%     after loading and editing events
%     If 0, no browsing, if 1, browse datafile{1}, if 2, browse datafile{2},...
%     Processing will stop after browser is launched
%     To restart, set this flag to 0 and rerun this function
%     (allows opportunity for manual rejection and event marking)
%     {default = 0 }
%  'saveepochs_flag': [1|0] instead of calculating averages, extract
%     epochs (single trial time courses)
%     {default = 0 }
%  'prefix': prefix of all output files
%    {default = 'proc'}
%  'rootoutdir': root output directory (several subdirectories will be created)
%    {default = pwd} (current working directory)
%  'write_fif_flag': [1|0] whether to output average sensor waveforms
%     as fif files (displayable by Neuromag's xplotter)
%     {default = 1}
%  'fifname_evcode_flag': whether to append output fif file names with
%     event codes (1) or with condition numbers (0)
%     {default = 1}
%  'events_fnames': name(s) of txt file(s) containing events
%    (see ts_import_events)
%    -- use 'events_fnames', followed by the filename
%       (in single quotes if bare text, not if a variable name)
%       (or cell array for multiple datafiles)
%       to import events into matlab from tab-delimited
%    -- or let this program read the events from the raw data's trigger channel
%    {default = []}
%  'trigchan': trigger channel label
%    {default = 'STI101'}
%  'evcode_offset': value subtracted from event codes
%    use this, for example, to correct for an always-on button
%    {default = []}
%  'event_recode_rules': cell array of strings defining how to 
%    change one event's condition to a new value.  See
%    ts_recode_events for more details.
%  'valid_event_codes': vector of event codes to be averaged
%    event codes are checked against this vector after recoding (if done)
%    so the valid event code is the newly recoded one
%    {default = []} (if empty, treat all event codes as valid)
%  'code_excl': exclude events surrounding this event code(s) (can be a vector)
%    use this, for example, to exclude events surrounding a button press
%    {default = []}
%  'time_excl_pre': time (msec) before code_excl to exclude other events
%    {default = 0}
%  'time_excl_post': time (msec) after code_excl to exclude other events
%    {default = 0}
%  'trig_minduration': minimum # of samples for a trigger code to be ON
%    before it is considered ACTUALLY ON and not just noise.
%    {default = []} (if empty, treat all event codes as valid)
%  'max_num_trials': maximum number of trials per condition
%    {default = Inf} (infinite)
%  'stim_delay': stimulus delay (msec) -- use to compensate for lag between
%    trigger and actual stimulus presentation
%    {default = 0}
%  'prestim_dur' : duration of prestimulus period (msec)
%     {default = 100}
%  'poststim_dur': duration of poststimulus period (msec)
%     {default = 400}
%  'stim_delay': duration of stimulus onset delay after trigger (msec)
%     {default = 0}
%  'reject_mag' : automatic rejection threshold for magnetometer channels (fT)
%     if 0, rejection based on magnetometers is disabled
%     {default = 10000}
%  'reject_grad': auto-rejection threshold for gradiometer channels (fT/cm)
%     if 0, rejection based on gradiometers is disabled
%     {default = 6000}
%  'reject_eeg': auto-rejection threshold for EEG channels (uV)
%     if 0, rejection based on eeg is disabled
%     {default = 0}
%  'reject_eog': auto-rejection threshold for EOG channel (uV)
%     if 0, rejection based on eog is disabled
%     {default = 200}
%  'bandpass_flag': [1|0] bandpass fft filter before averaging
%     {default = 0}
%  'bandpass_low_cf' : low center frequency (high-pass filter) (Hz)
%     {default = 0.2}
%  'bandpass_low_tb' : low transition band (Hz)
%     {default = 0.4}
%  'bandpass_high_cf': high center frequency (low-pass filter) (Hz)
%     {default = 100}
%  'bandpass_high_tb': high transition band (Hz)
%     {default = 10}
%  'notch_flag': [0|1] notch fft filter before averaging
%     {default = 0}
%  'notch_cf' : notch center frequency (notch filter) (Hz)
%     {default = 60}
%  'notch_tb' : notch transition band (Hz)
%     {default = 4}
%  'dsfact' : downsampling factor -- must be an integer
%     data is downsampled to lower sampling frequency
%     e.g. original sampling freq = 1000, dsfact = 4,
%         resulting sampling freq = 250
%     {default = 1 (no downsampling) }
%  'detrend_flag': [0|1] detrending of single trials before averaging
%     {default = 1}
%  'baseline_flag': [0|1] baseline subtraction of single trials before averaging
%     {default = 1}
%  'baseline_start': start time of baseline period (msec)
%     relative to trigger onset; negative times occur before trigger
%     {default = -Inf} (start at beginning of prestimulus period)
%  'baseline_end'  : end time of baseline period (msec)
%     relative to trigger onset; negative times occur before trigger
%     {default = 0} (end at trigger onset)
%  'ncov_ex_evnts': vector of event codes that should not
%     be used in calculating the noise covariance matrix
%     use this, for example, to exclude events with pre-stimulus activity
%     that is not noise -- for example, a button press
%     {default = []}
%  'badchans'   : vector of bad channel indices -- will be set to zero
%     {default = []}
%  'badchanfile': name of text file containing bad channel labels
%    {default = []}
%  'readtrans_flag' : [0|1] read device2head transform from fif file
%     if not found in fif file, will cause core dump
%     {default = 1}
%  'post_subnull_flag': [0|1] post-averaging null subtraction
%     i.e. subtract a specified "null" condition from all other conditions
%     {default = 0}
%  'null_event': event code of "null" event that is optionally subtracted
%     from all other conditions
%     if no event with with code is found, null subtraction will be skipped
%     {default = []}
%  'post_combcond_flag': [0|1] post-averaging combination of conditions
%     {default = 0}
%  'post_combcond_append_flag': [0|1] for combination of conditions
%     add new conditions to existing, otherwise, only include new conditions
%     {default = 1}
%  'fname_combcond_info': full path name of csv (comma-separated-value) file
%     containing 'neweventcodes' and 'combinations'
%     if not specified, post_combcond_flag will be set to 0
%     see ts_load_combcond_info and ts_combine_conditions
%     {default = []}
%  'post_bandpass_flag': [0|1] bandpass fft filter after averaging
%     {default = 0}
%  'post_notch_flag': [0|1] notch fft filter after averaging
%     {default = 0}
%  'post_dsfact' : downsampling factor applied after averaging
%     {default = 1 (no downsampling) }
%  'post_detrend_flag': [0|1] detrending of averaged waveforms
%     {default = 0}
%  'post_baseline_flag': [0|1] baseline subtraction of averaged
%     waveforms
%     {default = 0}
%  'post_stim_delay': stimulus delay (msec) -- correct for stimulus delay
%    after averaging
%    {default = 0}
%  'post_badchans': vector of bad channel indices
%     set to zero during post-averaging processing
%     {default = []}
%  'post_badchanfile': name of text file containing bad channel labels
%     set to zero during post-averaging processing
%    {default = []}
%  'post_rm_badchans_flag': [0|1] whether to completely remove badchans
%     from data matrices during post-processing -- otherwise just set to zero
%     {default = 0}
%  'ICA_auto_flag': [0|1] whether to use automatic ICA (independent
%     components analysis) for artifact (e.g. blink, EKG) removal
%     If 1, epochs are extracted and ICA is performed on epoch data
%     automatically using EOG/EKG/STIM reference channel
%     {default = 0}
%  'ICA_manual_flag': [0|1] whether to use manual ICA (independent
%     components analysis) for artifact (e.g. blink, EKG) removal
%     If 1, epochs are extracted and ICA is performed on epoch data
%     The user then selects ICs to remove and trials are averaged
%       if saveepochs_flag = 0
%     {default = 0}
%   ICA_ref_chan: channel name of EOG/EKG/STIM reference channel for auto ICA
%     {default = 'EOG061'}
%  'ICA_chantype': cell array of channel types to process (each set will be
%     processed individually)
%     {default = 'all'}
%  'ICA_maxsteps': maximum number of steps for ICA
%     {default = 20}
%  'ICA_ntrials': number of trials to display on screen for IC selection
%      only applies to 'activations' plot
%     {default = 5}
%  'ICA_ncomponents': the number of components to display on screen for IC
%      selection.  Note: use [] to display all
%     {default = 80}
%  'ICA_rescale_flag': [0|1] whether to perform rescaling on the ICA
%      processed data
%     {default = 1}
%  'ICA_sorttrials': [0|1] whether to sort the trials as well as components
%      before viewing
%     {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output files
%     {default = 0}
%  'saveperevent': [0|1] whether to save epoch_data with all events or one
%			epoch_data structure per event
%     {default = 0}
%
% Output:
%   filenames: cell array of output files created
%
%  Note: This program is designed so that the output of each processing step
%      is stored as a mat file (in the matfiles dir that will be automatically
%      created).
%    If you rerun this program with forceflag=0 (the default), it will skip any
%      steps that have already been run (if the mat file exists).
%    If you rerun this program with forceflag=1, it will repeat all steps
%      and overwrite any existing output files.
%    If you change parameters and try to rerun, they will be ignored unless you
%     use forceflag=1 or delete the proc_parms.mat file
%    If you want to rerun the averaging with different parameters, but do not
%      want to reread the event codes from the raw data file, you can delete the
%      proc_parms.mat and proc_avg_data*.mat files.
%    If you want to redo the post-processing but keep the original average,
%      delete proc_avg_data_post.mat and rerun this program.
%
% created:  04/26/06 by Don Hagler
% last mod: 11/04/14 by Don Hagler
%

%% todo: use MNE matlab toolbox instead of fiff access
%%  for reading events, data, etc.

%% todo: save parms in avg_data?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenames = {};
if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = { ...
  'rootoutdir',pwd,[],...
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
  'badchanfile',[],[],...
  'readtrans_flag',true,[false true],...
  'post_subnull_flag',false,[false true],...
  'null_event',[],[],...
  'post_combcond_flag',false,[false true],...
  'post_combcond_append_flag',true,[false true],...
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
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

% check input parameters
parms = ts_proc_check_input(datafile,varargin,parms_filter);

% load "events" from fif file or tab-delimited text file
ts_proc_load_events(parms);

% edit events by excluding invalid codes, recoding revents, etc.
ts_proc_edit_events(parms);

% optionally browse raw data
if parms.browseraw>0
  ts_proc_browse_raw(parms);
  fprintf('%s: to resume processing,\n change browseraw parameter to 0 and rerun\n',...
    mfilename);
  return;
end;

% extract epochs and calculate averages (if save_epochs_flag = 0)
[filenames,datamatfile] = ts_proc_avg_data(parms,filenames);

% combine averages (or epochs) from multiple fif files
[filenames,datamatfile] = ts_proc_combine_avg(parms,filenames,datamatfile);

% use ICA to remove blink and EKG artifacts
if parms.ICA_flag
  [filenames,datamatfile] = ts_proc_ICA(parms,filenames,datamatfile);
end;

% subtract null event from each condition
if parms.post_subnull_flag
  [filenames,datamatfile] = ts_proc_subnull(parms,filenames,datamatfile);
end;

% combine conditions (i.e. averaging or subtraction)
if parms.post_combcond_flag
  [filenames,datamatfile] = ts_proc_combcond(parms,filenames,datamatfile);
end;

% postprocessing: filter, downsample, detrend, baseline
if parms.post_dsfact>1 | parms.post_stim_delay~=0 |...
     parms.post_bandpass_flag | parms.post_notch_flag |...
     parms.post_baseline_flag | parms.post_detrend_flag |...
    ~isempty(parms.post_badchans) | ~isempty(parms.post_badchanfile)
  [filenames,datamatfile] = ts_proc_post(parms,filenames,datamatfile);
end;

% write to fif
if parms.write_fif_flag
  ts_proc_write_fif(parms,datamatfile);
end;

fprintf('%s: finished.\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = ts_proc_check_input(datafile,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  matdir = [parms.rootoutdir '/matfiles'];
  mmil_mkdir(matdir);
  if iscell(datafile)
    parms.num_data_files=length(datafile);
  else
    datafile={datafile};
    parms.num_data_files=1;
  end;
  if(~parms.num_data_files)
    error('no data files specified');
  else
    fprintf('%s: %d data file(s) to process\n',mfilename,parms.num_data_files);
  end
  parms.datafile = datafile;
  % check events files exist
  if ~isempty(parms.events_fnames)
    if ~iscell(parms.events_fnames)
      parms.events_fnames={parms.events_fnames};
    end;
    % length of events_fnames must match numdatafiles
    if length(parms.events_fnames) ~= parms.num_data_files
      error('# of events files (%d) does not match # datafiles (%d)',...
        length(parms.events_fnames),parms.num_data_files);
    end;
    for i=1:length(parms.events_fnames)
      if ~exist(parms.events_fnames{i},'file')
        error('events file %s not found',...
          parms.events_fnames{i});
      end;
    end;
  end;
  parms.ICA_flag = (parms.ICA_auto_flag || parms.ICA_manual_flag);
  if parms.post_subnull_flag    
    if isempty(parms.null_event)
      fprintf('%s: WARNING: null event not specified... setting post_subnull_flag = 0\n',...
        mfilename);
      parms.post_subnull_flag = 0;
    elseif length(parms.null_event)>1
      error('null event must be a single event code, not a vector');
    end;
  end;
  if parms.post_combcond_flag
    if isempty(parms.fname_combcond_info)
      fprintf('%s: WARNING: fname_combcond_info not specified... setting post_subnull_flag = 0\n',...
        mfilename);
      parms.post_combcond_flag = 0;
    else
      [parms.neweventcodes,parms.combinations,parms.combcond_info] = ...
        ts_load_combcond_info(parms.fname_combcond_info);
    end;
  end;
  if parms.write_fif_flag
    mmil_mkdir([parms.rootoutdir '/fifs']);
  end
  matfile = sprintf('%s/%s_parms.mat',matdir,parms.prefix);
  if ~exist(matfile,'file') || parms.forceflag
    save(matfile,'parms');
  else
    fprintf('%s: NOTICE: not overwriting parameter file %s\n',...
      mfilename,matfile);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ts_proc_load_events(parms)
  fprintf('%s: loading events...\n',mfilename);
  for f=1:parms.num_data_files
    fname = parms.datafile{f};
    matfile=sprintf('%s/matfiles/%s_events_%d.mat',...
      parms.rootoutdir,parms.prefix,f);
    evnts = []; hdr = [];
    if ~isempty(parms.events_fnames)
      % load hdr from events mat file read from raw data
      matfile=sprintf('%s/matfiles/%s_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
      if ~exist(matfile,'file') || parms.forceflag
        % read header
        fprintf('%s: reading header...\n',mfilename);
        tic;
        hdr = ts_read_fif_header(fname,1);
        toc;
        if isempty(hdr)
          error('hdr structure is empty');
        end;
      else
        load(matfile);
      end;
      % load events from tab-delimited text file
      fprintf('%s: importing events...\n',mfilename);
      matfile=sprintf('%s/matfiles/%s_imported_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
      evnts=ts_import_events(parms.events_fnames{f});
      save(matfile, 'evnts','hdr');
    else
      if ~exist(matfile,'file') || parms.forceflag
        if(~exist(fname,'file'))
          error('data file %s not found',fname);
        end
        fprintf('%s: reading event codes from %s...\n',mfilename,fname);
        tic;
        [hdr,evnts] = ...
          ts_read_fif_events(fname,parms.trigchan,parms.evcode_offset, parms.trig_minduration);
        toc;
        if isempty(evnts)
%          error('evnts structure is empty -- check trigger channel');
          fprintf('%s: WARNING: evnts structure is empty -- check trigger channel\n',...
            mfilename);
        end;
        matfile=sprintf('%s/matfiles/%s_events_%d.mat',...
          parms.rootoutdir,parms.prefix,f);
        save(matfile, 'evnts','hdr');
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ts_proc_edit_events(parms)
  for f=1:parms.num_data_files
    if ~isempty(parms.events_fnames)
      matfile=sprintf('%s/matfiles/%s_edited_imported_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
    else
      matfile=sprintf('%s/matfiles/%s_edited_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
    end;
    if ~exist(matfile,'file')  || parms.forceflag
      if ~isempty(parms.events_fnames)
        matfile=sprintf('%s/matfiles/%s_imported_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
      else
        matfile=sprintf('%s/matfiles/%s_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
      end;
      if(~exist(matfile,'file'))
        error('mat file %s not found',matfile);
      end
      evnts = []; hdr = [];
      load(matfile);
      fprintf('%s: editing event codes...\n',mfilename);
      if ~isempty(evnts)
        evnts = ts_edit_events(evnts, ...
                               hdr.sfreq,...
                               parms.valid_event_codes,...
                               parms.code_excl, ...
                               parms.time_excl_pre, ...
                               parms.time_excl_post, ...
                               parms.event_recode_rules);
        if isempty(evnts)
%          error('edited evnts structure is empty');
          fprintf('%s: WARNING: edited evnts structure is empty\n',...
            mfilename);
        end;
      end;
      if isempty(evnts)
        % create dummy event
        evnts(1).type = 'trigger';
        evnts(1).duration = 1;
        evnts(1).condition = -1;
        evnts(1).latency = 1 + round(parms.prestim_dur * hdr.sfreq / 1000);
      end;
      if ~isempty(parms.events_fnames)
        matfile=sprintf('%s/matfiles/%s_edited_imported_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
      else
        matfile=sprintf('%s/matfiles/%s_edited_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
      end;
      save(matfile, 'evnts','hdr');
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ts_proc_browse_raw(parms)
  f = parms.browseraw;
  if f > parms.num_data_files
    error('bad value for browseraw: %d (must be less than or equal to number of data files (%d)',...
      f,parms.num_data_files);
  end;
  fname = parms.datafile{f};
  fprintf('%s: browsing data file %s...\n',mfilename,fname);
  % get existing events
  if ~isempty(parms.events_fnames)
    % events imported from csv file
    matfile=sprintf('%s/matfiles/%s_edited_imported_events_%d.mat',...
      parms.rootoutdir,parms.prefix,f);
  else
    % manually edited events
    matfile=sprintf('%s/matfiles/%s_manual_events_%d.mat',...
      parms.rootoutdir,parms.prefix,f);
    if ~exist(matfile,'file')
      % automatically edited events (possibly recoded)
      matfile=sprintf('%s/matfiles/%s_edited_events_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
    end;
  end;
  if(~exist(matfile,'file'))
    error('mat file %s not found',matfile);
  end
  load(matfile);
  matfile=sprintf('%s/matfiles/%s_manual_events_%d.mat',...
    parms.rootoutdir,parms.prefix,f);
  if parms.reject_mag || parms.reject_grad ||...
     parms.reject_eeg || parms.reject_eog
    autoreject = 1;
  else
    autoreject = 0;
  end;
  ts_browseraw(fname,'events',evnts,'hdr',hdr,...
    'autoreject',autoreject,...
    'reject_mag',parms.reject_mag,...
    'reject_grad',parms.reject_grad,...
    'reject_eeg',parms.reject_eeg,...
    'reject_eog',parms.reject_eog,...
    'filter',parms.bandpass_flag,...
    'dsfact',parms.dsfact,...
    'filter_low_cf',parms.bandpass_low_cf,...
    'filter_low_tb',parms.bandpass_low_tb,...
    'filter_high_cf',parms.bandpass_high_cf,...
    'filter_high_tb',parms.bandpass_high_tb,...
    'filter_notch',parms.notch_flag,...
    'filter_notch_cf',parms.notch_cf,...
    'filter_notch_tb',parms.notch_tb,...
    'prestim_dur',parms.prestim_dur,...
    'poststim_dur',parms.poststim_dur,...
    'badchans',parms.badchans,...
    'badchanfile',parms.badchanfile,...
    'outfile',matfile);
  % if user saves events, will write to matfile
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,datamatfile] = ts_proc_avg_data(parms,filenames)
  datamatfile = [];
  % average data
  fprintf('%s: processing data...\n',mfilename);
  for f=1:parms.num_data_files
    fname = parms.datafile{f};
    if parms.saveepochs_flag || parms.ICA_flag
      % if using ICA artifact removal, save_epochs no matter what
      matfile=sprintf('%s/matfiles/%s_epoch_data_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
    else
      matfile=sprintf('%s/matfiles/%s_avg_data_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
    end;
    if exist(matfile,'file')  && ~parms.forceflag
      % only run averager if output does not already exist (unless forceflag)
      fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
    else
      if ~isempty(parms.events_fnames)
        matfile=sprintf('%s/matfiles/%s_edited_imported_events_%d.mat',...
          parms.rootoutdir,parms.prefix,f);
      else
        matfile=sprintf('%s/matfiles/%s_manual_events_%d.mat',...
          parms.rootoutdir,parms.prefix,f);
        if ~exist(matfile,'file')
          matfile=sprintf('%s/matfiles/%s_edited_events_%d.mat',...
            parms.rootoutdir,parms.prefix,f);
        end;
      end;
      if(~exist(matfile,'file'))
        error('mat file %s not found',matfile);
      end
      evnts = []; hdr = [];
      load(matfile);
      if ~isempty(evnts)
        if(~exist(fname,'file'))
          error('data file %s not found',fname);
        end
        tic;
        if parms.saveepochs_flag || parms.ICA_flag
          fprintf('%s: extracting epochs...\n',mfilename);
        else
          fprintf('%s: averaging data...\n',mfilename);
        end;
        tags = {'prestim_dur','poststim_dur','baseline_flag',...
          'baseline_start','baseline_end','detrend_flag',...
          'stim_delay','ncov_ex_evnts','badchans','badchanfile',...
          'dsfact','bandpass_flag','bandpass_low_cf','bandpass_low_tb',...
          'bandpass_high_cf','bandpass_high_tb','notch_flag',...
          'notch_cf','notch_tb','reject_mag','reject_grad',...
          'reject_eeg','reject_eog','readtrans_flag','max_num_trials','datatype'};
        args = mmil_parms2args(parms,tags);
        avg_data = ts_avg_fif_data(fname,args{:},'evnts',evnts,'hdr',hdr,...
          'save_epochs_flag',(parms.saveepochs_flag || parms.ICA_flag));
        toc;
      else
        avg_data = [];
      end;
      if parms.saveepochs_flag || parms.ICA_flag
        epoch_data = avg_data;
			  if parms.saveperevent
          epoch_data.opt = parms;
				  tmpfile = ts_proc_save_per_event(epoch_data,parms.rootoutdir,parms.prefix,'');
                  filenames = {filenames{:} tmpfile};
			  else
	        matfile=sprintf('%s/matfiles/%s_epoch_data_%d.mat',...
  	      parms.rootoutdir,parms.prefix,f);
      	  save(matfile, 'epoch_data');
          filenames = {filenames{:} matfile};
			  end
      else
        matfile=sprintf('%s/matfiles/%s_avg_data_%d.mat',...
        parms.rootoutdir,parms.prefix,f);
        save(matfile, 'avg_data');
        filenames = {filenames{:} matfile};
      end;
    end;
    datamatfile=matfile;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,datamatfile] = ts_proc_combine_avg(parms,filenames,datamatfile)
  % combine averages (or epochs)
  % prevent empty valid_event_codes
  if parms.saveperevent &&...
      isempty(parms.valid_event_codes) && exist('evnts','var')
    parms.valid_event_codes = unique([evnts.condition]);
  end;
  if(parms.num_data_files>1) || (length(parms.valid_event_codes)>1 && parms.saveperevent)
    if parms.saveepochs_flag || parms.ICA_flag
      matfile=sprintf('%s/matfiles/%s_epoch_data.mat',...
        parms.rootoutdir,parms.prefix);
    else
      matfile=sprintf('%s/matfiles/%s_avg_data.mat',...
        parms.rootoutdir,parms.prefix);
    end;
    if exist(matfile,'file')  && ~parms.forceflag
      fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
    else
      if parms.saveepochs_flag || parms.ICA_flag
        fprintf('%s: combining epochs from multiple files...\n',mfilename);
        if parms.saveperevent   
          for f=1:length(parms.valid_event_codes)
            matfile=sprintf('%s/matfiles/%s_epoch_data_cond%d.mat',...
              parms.rootoutdir,parms.prefix,parms.valid_event_codes(f));
            if(~exist(matfile,'file'))
              error('mat file %s not found',matfile);
            end
            load(matfile);
            epoch_data_arr(f) = epoch_data;
          end
        else
          for f=1:parms.num_data_files
            matfile=sprintf('%s/matfiles/%s_epoch_data_%d.mat',...
              parms.rootoutdir,parms.prefix,f);
            if(~exist(matfile,'file'))
              error('mat file %s not found',matfile);
            end;
            load(matfile);
            epoch_data_arr(f) = epoch_data;
          end;
        end;
        epoch_data = ts_combine_data(epoch_data_arr);
        matfile=sprintf('%s/matfiles/%s_epoch_data.mat',...
          parms.rootoutdir,parms.prefix);
        save(matfile, 'epoch_data');
        filenames = {filenames{:} matfile};
      else
        fprintf('%s: combining averages from multiple files...\n',mfilename);
        for f=1:parms.num_data_files
          matfile=sprintf('%s/matfiles/%s_avg_data_%d.mat',...
            parms.rootoutdir,parms.prefix,f);
          if(~exist(matfile,'file'))
            error('mat file %s not found',matfile);
          end
          load(matfile);
          avg_data_arr(f) = avg_data;
        end;
        avg_data = ts_combine_data(avg_data_arr);
        matfile=sprintf('%s/matfiles/%s_avg_data.mat',...
          parms.rootoutdir,parms.prefix);
        save(matfile, 'avg_data');
        filenames = {filenames{:} matfile};
      end;
    end;
    datamatfile=matfile;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,datamatfile] = ts_proc_ICA(parms,filenames,datamatfile)
  % ICA artifact removal
  matfile=sprintf('%s/matfiles/%s_epoch_data_ICA.mat',...
    parms.rootoutdir,parms.prefix);
  if exist(matfile,'file')  && ~parms.forceflag
    fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
  else
    if ~exist(datamatfile,'file')
      error('mat file %s not found',datamatfile);
    else
      load(datamatfile);
    end;
    % perform automatic ICA rejection for blinks
    if parms.ICA_auto_flag
      epoch_data = ts_autoICA(epoch_data,...
        'ICA_ref_chan',parms.ICA_ref_chan,...
        'chantype',parms.ICA_chantype,...
        'rescale',parms.ICA_rescale_flag);
    end
    save(matfile,'epoch_data');
    filenames = {filenames{:} matfile};

    % perform manual ICA rejection for blinks and/or EKG
    if parms.ICA_manual_flag
      epoch_data = ts_manualICA(epoch_data,...
        'maxsteps',parms.ICA_maxsteps,...
        'ntrial',parms.ICA_ntrial,...
        'ncomponents',parms.ICA_ncomponents,...
        'chantype',parms.ICA_chantype,...
        'rescale',parms.ICA_rescale_flag,...
        'sorttrials',parms.ICA_sorttrials,...
        'allconditions',1);
    end
    if parms.saveperevent
      epoch_data.opt = parms;
      tmpfile = ts_proc_save_per_event(epoch_data,parms.rootoutdir,parms.prefix,'_ICA');
      filenames = {filenames{:} tmpfile{:}};
    else
      save(matfile,'epoch_data');
      filenames = {filenames{:} matfile};
    end
  end;
  datamatfile = matfile;

  % convert epoch to avg if saveepochs_flag=0
  if ~parms.saveepochs_flag
    datamatfile = sprintf('%s/matfiles/%s_epoch_data_ICA.mat',...
      parms.rootoutdir,parms.prefix);
    matfile=sprintf('%s/matfiles/%s_avg_data_ICA.mat',...
      parms.rootoutdir,parms.prefix);
    if exist(matfile,'file')  && ~parms.forceflag
      fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
    else
      if(~exist(datamatfile,'file'))
        error('mat file %s not found',datamatfile);
      else
        load(datamatfile);
      end;
      avg_data = ts_epoch2avg(epoch_data);
      save(matfile,'avg_data');
      filenames = {filenames{:} matfile};
    end;
    datamatfile = matfile;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,datamatfile] = ts_proc_subnull(parms,filenames,datamatfile)
  %% todo: null subtraction of epochs
  if parms.saveepochs_flag
    fprintf('%s: WARNING: null subtraction for epochs not currently supported\n',...
      mfilename);
    return;
  end;
  matfile=sprintf('%s/matfiles/%s_avg_data_subnull.mat',...
    parms.rootoutdir,parms.prefix);
  if exist(matfile,'file')  && ~parms.forceflag
    fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
  else
    if ~exist(datamatfile,'file')
      error('mat file %s not found',datamatfile);
    else
      load(datamatfile);
      evcodes = cell2mat({avg_data.averages.event_code});
      null_cond = find(evcodes==parms.null_event);
      if isempty(null_cond)
        fprintf('%s: WARNING: null event %d not found... skipping null subtraction\n',...
          mfilename,parms.null_event);
        matfile = datamatfile;
      else
        fprintf('%s: subtracting null event from averages...\n',mfilename);
        nconds=length(avg_data.averages);
        avg_data_subnull = avg_data;
        for j=1:nconds
          avg_data_subnull.averages(j).data = avg_data.averages(j).data - ...
            avg_data.averages(null_cond).data;
        end
        avg_data = avg_data_subnull;
        save(matfile,'avg_data');
        filenames = {filenames{:} matfile};
      end;
    end
  end;
  datamatfile=matfile;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,datamatfile] = ts_proc_combcond(parms,filenames,datamatfile)
  if parms.saveepochs_flag
    fprintf('%s: WARNING: combination of conditions for epochs not supported\n',...
      mfilename);
    return;
  end;
  matfile=sprintf('%s/matfiles/%s_avg_data_combcond.mat',...
    parms.rootoutdir,parms.prefix);
  if exist(matfile,'file')  && ~parms.forceflag
    fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
  else
    if ~exist(datamatfile,'file')
      error('mat file %s not found',datamatfile);
    else
      load(datamatfile);
      fprintf('%s: combining averages of multiple conditions...\n',mfilename);
      avg_data = ts_combine_conditions(avg_data,...
        'combinations',parms.combinations,...
        'neweventcodes',parms.neweventcodes);
      if ~parms.post_combcond_append_flag
        evcodes = cell2mat({avg_data.averages.event_code});
        ind = find(ismember(evcodes,parms.neweventcodes));
        avg_data.averages = avg_data.averages(ind);
      end;
      save(matfile,'avg_data');
      filenames = {filenames{:} matfile};
    end
  end;
  datamatfile=matfile;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,datamatfile] = ts_proc_post(parms,filenames,datamatfile)
  %% todo: postprocessing of epochs
  if parms.saveepochs_flag
    fprintf('%s: WARNING: postprocessing for epochs not currently supported\n',...
      mfilename);
    return;
  end;
  matfile=sprintf('%s/matfiles/%s_avg_data_post.mat',...
    parms.rootoutdir,parms.prefix);
  if exist(matfile,'file')  && ~parms.forceflag
    fprintf('%s: NOTICE: not overwriting output %s\n', mfilename, matfile);
  else
    if ~exist(datamatfile,'file')
      error('mat file %s not found',datamatfile);
    else
      load(datamatfile);
      avg_data = ts_postprocess_avg(avg_data,...
        'baseline_flag',parms.post_baseline_flag,...
        'baseline_start',parms.baseline_start,...
        'baseline_end',parms.baseline_end,...
        'detrend_flag',parms.post_detrend_flag,...
        'stim_delay',parms.post_stim_delay,...
        'badchans',parms.post_badchans,...
        'badchanfile',parms.post_badchanfile,...
        'rm_badchans_flag',parms.post_rm_badchans_flag,...
        'dsfact',parms.post_dsfact,...
        'bandpass_flag',parms.post_bandpass_flag,...
        'bandpass_low_cf',parms.bandpass_low_cf,...
        'bandpass_low_tb',parms.bandpass_low_tb,...
        'bandpass_high_cf',parms.bandpass_high_cf,...
        'bandpass_high_tb',parms.bandpass_high_tb,...
        'notch_flag',parms.post_notch_flag,...
        'notch_cf',parms.notch_cf,...
        'notch_tb',parms.notch_tb,...
        'datatype',parms.datatype);
      save(matfile,'avg_data');
      filenames = {filenames{:} matfile};
    end;
  end;
  datamatfile=matfile;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ts_proc_write_fif(parms,datamatfile)
  if parms.saveepochs_flag
    fprintf('%s: WARNING: writing to fif for epochs not supported\n',...
      mfilename);
    return;
  end;
  matfile=datamatfile;
  if(~exist(matfile,'file'))
    error('mat file %s not found',matfile);
  elseif (~exist(parms.datafile{1},'file'))
    error('fif file %s not found',parms.datafile{1});
  else
    load(matfile);
    outstem=sprintf('%s/fifs/%s_avg',parms.rootoutdir,parms.prefix);
    if ~isempty(avg_data)
      ts_avg2fif(avg_data,parms.datafile{1},outstem,...
        parms.fifname_evcode_flag,parms.forceflag);
    end;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filenames = ts_proc_save_per_event(data,rootoutdir,prefix,tag)
  filenames = {};
  epoch_hdr = rmfield(data,'epochs');
  evcodes = [data.epochs.event_code];
  for c = 1:length(evcodes)
	  epoch_data = epoch_hdr;
	  epoch_data.epochs = data.epochs(c);
	  matfile = sprintf('%s/matfiles/%s_epoch_data%s_cond%d.mat',...
      rootoutdir,prefix,tag,evcodes(c));
    epoch_data.opt.filename{1} = matfile;
	  save(matfile,'epoch_data');
      filenames = {filenames{:} matfile};
	  clear epoch_data;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

