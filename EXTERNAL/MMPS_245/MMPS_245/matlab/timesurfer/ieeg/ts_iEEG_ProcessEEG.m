function filename = ts_iEEG_ProcessEEG (varargin)
% function varargout = ts_iEEG_ProcessEEG (varargin)

% Use: ts_iEEG_ProcessEEG ('eegpath',eegpath,...);
%
% Used to process .eeg NeuroScan data and save TimeSurfer avg_data and
% epoch_data .mat files.  Note that each epoch is saved in its own matfile
% due to memory limitation when dealing with large number of trials.
% Preprocessing options are only applicable to the average data.
%
% Required Input:
%
%   <None> - User will be prompted for required inputs.
%
% Optional Input:
%
%   eegpath       - directory containing eegfile
%   eegfile       - .eeg NeuroScan formated file(s) or .vhdr BrainVision
%                    file
%   savefile      - path/name to use as basis to save avg_data and
%                   epoch_data
%   badchanfile   - path/name of text file containing bad channels
%   channamefile  - path/name of text file containing names of channels
%   saveavgs      - 'yes' or 'no' {default = 'yes'}
%   saveepochs    - 'yes','no','events','conditions' {default = 'no' }
%   oldeventcodes - array or cell of arrays (multi file) of old event codes
%   neweventcodes - array or cell of arrays of new event codes
%   rejectfile    - name of mat file containing information on previous
%                   artifact rejection
%   reject        - 'yes','no','before_ica','after_ica','both'
%                    manual rejection processing {default = 'no'}
%                   Note: Ignored if provided a rejectfile.
%   reject_method - method option for rejectvisual: 'summary', 'channel', 'trial'
%                   {default = 'summary'}
%   reject_metric - metric option for rejectvisual: 'var', 'min', 'max',
%                    'absmax', 'rang' {default = 'var'}
%   auto_ica      - 'yes' or 'no' to perform automatic ica rejection
%                    using specified reference channel.
%   refchan       - name(s) of channels to use as reference in auto_ica
%   manual_ica    - perform manual ica rejection
%   ica_rescale   - rescale the ica data after rejection
%   ica_allconds  - combine all conditions prior to ICA
%   ica_sorttrials- sort trials by variance prior to display
%   stim_delay    - Adjusts the time vector so that the time point
%                   specified by stim_delay (in sec) is set to time point 0
%                   and adjusts the rest of the time vector accordingly.
%                   This is done prior to any processing so any time
%                   specifications should be based on new timing.
%   noise_start   - start of noise period in sec (default - -.08sec)
%   noise_end     - end of noise period in sec (defualt - 0 sec)
%   timelimits    - start and end time of trials when loading a BrainVision
%                   file (default = [-2 2])
%   combinations  - a list of combinations to produce.  Each combination is
%                   supplied as a string.  The numbers within the string are
%                   treated as event codes. Valid examples are:
%                   To combine events 1,2 & 3:   '1+2+3'
%                    To subtract event 10 from 2: '2-10' (Note: Only two events at a time
%                                                         is valid for subtractions)
%   comboeventcodes - a list of event codes for the new conditions
%   combo_calc     - 'weighted' or 'avg' - specifies whether to perform a weighted
%                   average or a straight average when combining conditions using the '+'
%                    operation
%   recode_rules   - string or cell array of strings, of recoding rules
%     Each rule specifies a sequence of events and the new event code which is
%       assigned to the first event in the seqeuence
%     Examples:
%       '2->5=101' : event 2 followed directly by event 5 gets recoded as event 101
%       '[2,3,4]->5=101' : events 2, 3, or 4, followed by event 5
%       '[2]->[5,10]=101' : event 2 followed by event 5 or 10
%       '2->200:1500[5]=101' : event 2 followed by event 5 within 200 and 1500 msec
%     '->' separates events in the sequence
%     '<-' same as about but instead of the first events following it
%          specified if those events come before the others
%     '=' precedes the new event code
%     ':' is used to separate start and end times (in msec) of time window
%     Lists of event codes are placed in square brackets [], separated by commas
%     Square brackets are also required when an event code follows a time window
%   additioncombos    - specified just like combinations above but this will
%                       ADD averages not create weighted averages
%   additioneventcodes- the list of events codes for the additioncomboes
%   
%
%   FieldTrip preprocessing options available:
%     
%     padding, lpfilter, hpfilter, bpfilter, lnfilter, dftfilter, 
%     medianfilter, lpfreq, hpfreq, bpfreq, lnfreq, dftfreq, lpfiltord,
%     hpfiltord, bpfiltord, lnfiltord, lpfilttype, hpfilttype, bpfilttype,
%     lbfiltdir, hpfiltdir, bpfiltdir, medianfiltord, detrend, blc,
%     blcwindow, hilbert, rectify, precision
%
% Created:  09/20/07 by Rajan Patel
% Rcnt Mod: 05/15/08 by Rajan Patel
% Last Mod: 09/15/12 by Don Hagler
%
% See also: preprocessing, ts_iEEG_eeg2epoch, ts_autoICA, ts_manualICA,
% ts_combine_conditions

%% Check Options
%  THIS IS THE OLD WAY.  THIS SHOULD BE UPDATED TO THE FORMAT
%  OF OPT = MMIL_ARGS2PARMS(... WITH ALL THE PARAMETERS LISTED
%
% todo: set badchannels to zero instead of nan

try
    options = varargin;
    for index = 1:length(options)
        if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
    end;
    if ~isempty( varargin ), opt=struct(options{:});
    else opt = []; end;
catch
    fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
    return;
end;

if isfield(opt,'precision') && ~isempty(opt.precision)
  precision = opt.precision;
else
  precision = 'single';
end

% set up condition key for later use if present
if isfield(opt,'conditionkey') && ~isempty(opt.conditionkey)
  conditionkey = opt.conditionkey;
  if ~exist(conditionkey,'file')
    warning('Cannot find condition key: %s.\ncond_key.mat will not be made.\n',conditionkey);
    conditionkey = [];
  end
  fprintf('Making condition key from %s.\n',conditionkey);
  ts_iEEG_makecondkey(conditionkey);
end

filename = {};
if isfield(opt,'filename')
  opt.savefile = opt.filename{1};
end
if isfield(opt,'datapath')
  opt.eegpath = opt.datapath;
end
if isfield(opt,'datafile')
  if ~iscell(opt.datafile), opt.datafile={opt.datafile}; end
  for i = 1:length(opt.datafile)
    [x1 x2 x3 x4] = fileparts(opt.datafile{i});
    if isfield(opt,'eegfile') && ~iscell(opt.eegfile)
      opt.eegfile = {opt.eegfile};
    end
    opt.eegfile{i} = [x2 x3];
  end
end

% set up input file

if isfield(opt,'eegpath')
  eeg_path = opt.eegpath;
else
  eeg_path = [];
end

if isfield(opt,'eegfile')
    eeg_file = opt.eegfile;
else
    [eeg_file,eeg_path,index]=uigetfile({'*.mat'},'Select data file','MultiSelect','on');
end
if ~iscell(eeg_file), eeg_file = {eeg_file}; end;

if isfield(opt,'rootoutdir') && isfield(opt,'outpath')
  [f_path, f_name, ext] = fileparts(eeg_file{1});
  f_path = fullfile(opt.rootoutdir,opt.outpath);
elseif isfield(opt,'rootoutdir')
  [f_path, f_name, ext] = fileparts(eeg_file{1});
  f_path = opt.rootoutdir;
elseif isfield(opt,'savefile')
   [f_path, f_name, ext] = fileparts(opt.savefile);
else
   [f_path, f_name, ext] = fileparts(eeg_file{1});
end

if ~exist(f_path,'dir')
    fprintf('Creating directory: %s\n',f_path);
    mkdir(f_path);
end

if isfield(opt,'badchanfile') && ~isempty(opt.badchanfile)
    pa.badchanfile=opt.badchanfile;
    if ~exist(pa.badchanfile,'file')
      fprintf('Unable to locate bad channel file: %s\n',pa.badchanfile);
      pa.badchanfile = [];
    else
      fid = fopen(pa.badchanfile);
      tline = fgetl(fid);
      if ~ischar(tline)
        fprintf('Bad channel file is empty: %s\n',pa.badchanfile);
        pa.badchanfile = [];
      end
      fclose(fid);
    end
else
    pa.badchanfile = [];
end

if isfield(opt,'channamefile') && ~isempty(opt.channamefile)
  pa.channamefile = opt.channamefile;
  if ~exist(pa.channamefile,'file')
    fprintf('Unable to locate channel name file: %s\n',pa.channamefile);
    pa.channamefile = [];
  else
    fid = fopen(pa.channamefile);
    tline = fgetl(fid);
    if ~ischar(tline)
      fprintf('Channel name file is empty: %s\n',pa.channamefile);
      pa.channamefile = [];
    end
    fclose(fid);
  end
else
  pa.channamefile = [];
end

save_avgs = 1;
if isfield(opt,'saveavgs')
 if strcmpi(opt.saveavgs,'no'), save_avgs = 0; end
end

save_epochs = 0;
if isfield(opt,'saveepochs')
 if strcmpi(opt.saveepochs,'yes') || strcmpi(opt.saveepochs,'events')
   save_epochs = 1; 
   epoch_tag   = 'event';
 elseif strcmpi(opt.saveepochs,'conditions')
   save_epochs = 1;
   epoch_tag   = 'cond';
 end
end

if isfield(opt,'oldeventcodes'), 
  oldeventcodes = opt.oldeventcodes;
  if ~iscell(oldeventcodes), oldeventcodes={oldeventcodes}; end
else
 for i = 1:length(eeg_file)
  oldeventcodes{i} = []; 
 end
end
    
if isfield(opt,'neweventcodes'), 
  neweventcodes = opt.neweventcodes;
  if ~iscell(neweventcodes), neweventcodes={neweventcodes}; end
else
 for i = 1:length(eeg_file)
  neweventcodes{i} = []; 
 end
end

if length(oldeventcodes)~=length(eeg_file) || length(neweventcodes)~=length(eeg_file)
  fprintf('The number of lists of old or new event codes do not match number of eeg_files.\n');
  return;
end

if isfield(opt,'stim_delay')
  stim_delay = opt.stim_delay;
else
  stim_delay = [];
end

if isfield(opt,'rejectfile')
  reject_file = opt.rejectfile;
  if ~isempty(reject_file) && ~exist(reject_file,'file')
    fprintf('Could not locate rejection file: %s.\n',reject_file);
    reject_file = [];
  end
else
  reject_file = [];
end

if isfield(opt,'timelimits')
    timelimits = opt.timelimits;
else
    timelimits = [];
end

reject = [];
prepro = [];
prepro.feedback = 'no';
if isfield(opt,'padding'),      prepro.padding       = opt.padding;end
if isfield(opt,'lpfilter'),     prepro.lpfilter      = opt.lpfilter; end
if isfield(opt,'hpfilter'),     prepro.hpfilter      = opt.hpfilter; end
if isfield(opt,'bpfilter'),     prepro.bpfilter      = opt.bpfilter; end
if isfield(opt,'lnfilter'),     prepro.lnfilter      = opt.lnfilter; end
if isfield(opt,'dftfilter'),    prepro.dftfilter     = opt.dftfilter; end
if isfield(opt,'medianfilter'), prepro.medianfilter  = opt.medianfilter; end
if isfield(opt,'lpfreq'),       prepro.lpfreq        = opt.lpfreq;
                                reject.high_cf       = opt.lpfreq; end                            
if isfield(opt,'hpfreq'),       prepro.hpfreq        = opt.hpfreq; 
                                reject.low_cf        = opt.hpfreq; end
if isfield(opt,'bpfreq'),       prepro.bpfreq          = opt.bpfreq;
                                reject.bandpass_low_cf = opt.bpfreq(1);
                                reject.bandpass_high_cf = opt.bpfreq(2); end
if isfield(opt,'lnfreq'),       prepro.lnfreq        = opt.lnfreq; 
                                reject.linefreq      = opt.lnfreq;  end
if isfield(opt,'dftfreq'),      prepro.dftfreq       = opt.dftfreq; end
if isfield(opt,'lpfiltord'),    prepro.lpfiltord     = opt.lpfiltord; end
if isfield(opt,'hpfiltord'),    prepro.hpfiltord     = opt.hpfiltord; end
if isfield(opt,'bpfiltord'),    prepro.bpfiltord     = opt.bpfiltord; end
if isfield(opt,'lnfiltord'),    prepro.lnfiltord     = opt.lnfiltord; end
if isfield(opt,'lpfilttype'),   prepro.lpfilttype    = opt.lpfilttype; end
if isfield(opt,'hpfilttype'),   prepro.hpfilttype    = opt.hpfilttype; end
if isfield(opt,'bpfilttype'),   prepro.bpfilttype    = opt.bpfilttype; end
if isfield(opt,'lpfiltdir'),    prepro.lpfiltdir     = opt.lpfiltdir; end
if isfield(opt,'hpfiltdir'),    prepro.hpfiltdir     = opt.hpfiltdir; end
if isfield(opt,'bpfiltdir'),    prepro.bpfiltdir     = opt.bpfiltdir; end
if isfield(opt,'medianfiltord'),prepro.medianfiltord = opt.medianfiltord; end
if isfield(opt,'detrend'),      prepro.detrend       = opt.detrend; end
if isfield(opt,'blc'),          prepro.blc           = opt.blc; end
if isfield(opt,'blcwindow'),    prepro.blcwindow     = opt.blcwindow; 
                                reject.baseline_start = opt.blcwindow(1)*1000;
                                reject.baseline_end   = opt.blcwindow(2)*1000; end
if isfield(opt,'hilbert'),      prepro.hilbert       = opt.hilbert; end
if isfield(opt,'rectify'),      prepro.rectify       = opt.rectify; end
if isfield(opt,'precision'),    prepro.precision     = opt.precision; end
if isfield(opt,'polyremoval'),  prepro.polyremoval   = opt.polyremoval; end
if isfield(opt,'polyorder'),    prepro.polyorder     = opt.polyorder; end

% ONLY ADD THE PROCESSING OPTIONS TO THE AVG_DATA NAME AS THE AVG_DATA SET
% IS THE ONLY TO GET FILTERS APPLIED

a_name = f_name;

if isfield(prepro,'lpfilter')
    if strcmpi(prepro.lpfilter,'yes')
        a_name = [a_name '.lp' strrep(num2str(prepro.lpfreq),'.','')];
%        reject.lowpass = 1;
    end
end

if isfield(prepro,'hpfilter')
    if strcmpi(prepro.hpfilter,'yes')
        a_name = [a_name '.hp' strrep(num2str(prepro.hpfreq),'.','')];
%        reject.highpass = 1;
    end
end

if isfield(prepro,'bpfilter')
  if strcmpi(prepro.bpfilter,'yes')
    a_name = [a_name '.bp' strrep(num2str(prepro.bpfreq(1)),'.','') '-' strrep(num2str(prepro.bpfreq(2)),'.','')];
%    reject.bandpass = 1;
  end
end

if isfield(prepro,'lnfilter')
    if strcmpi(prepro.lnfilter,'yes')
        a_name = [a_name '.ln' strrep(num2str(prepro.lnfreq),'.','')];
%        reject.linefilt = 1;
    end
end

if isfield(prepro,'blc')
    if strcmpi(prepro.blc,'yes')
        a_name = [a_name '.bc'];
%        reject.baseline_sub = 1;
    end
end

if isfield(prepro,'detrend')
    if strcmpi(prepro.detrend,'yes')
        a_name = [a_name '.dt'];
        reject.detrend_events = 1;
    end
end

auto_ica = 0;
manual_ica = 0;
if (isfield(opt,'auto_ica') && strcmpi(opt.auto_ica,'yes')) || ...
   (isfield(opt,'manual_ica') && strcmpi(opt.manual_ica,'yes'))
    a_name = [a_name '.ica'];
end

if (isfield(opt,'auto_ica') && strcmpi(opt.auto_ica,'yes'))
    auto_ica = 1;
    if isfield(opt,'refchan')
        refchan = opt.refchan;
    else
        error('%s: No ''refchan'' option specified for auto ica.',mfilename);
    end
end

if (isfield(opt,'manual_ica') && strcmpi(opt.manual_ica,'yes'))
    manual_ica = 1;
end

if isfield(opt,'ica_rescale') && strcmpi(opt.ica_rescale,'yes')
  ica_rescale = 1;
else
  ica_rescale = 0;
end

if isfield(opt,'ica_sorttrials') && strcmpi(opt.ica_sorttrials,'no')
    ica_sorttrials = 0;
else
    ica_sorttrials = 1;
end

if isfield(opt,'ica_allconds') && strcmpi(opt.ica_allconds,'no')
    ica_allconds = 0;
else
    ica_allconds = 1;
end

reject_flag = 0; 
if isfield(opt, 'reject')
   if strcmpi(opt.reject,'yes') || ...
      strcmpi(opt.reject,'before_ica') || ...
      strcmpi(opt.reject,'after_ica') ||...
      strcmpi(opt.reject,'both')
       a_name = [a_name '.rej'];
       reject.feedback = 'no';
       reject_flag = 1;
       reject.chantype = 'eeg';
       if isfield(opt,'reject_metric'), reject.metric = opt.reject_metric;
       else reject.metric = 'var'; end
       if isfield(opt,'reject_method'), reject.method = opt.reject_method;
       else reject.method = 'summary'; end    
   end
end

if isfield(opt,'noise_start') && ~isempty(opt.noise_start)
  noise_start = opt.noise_start;
else
  noise_start = -0.080;
end

if isfield(opt,'noise_end') && ~isempty(opt.noise_end)
  noise_end = opt.noise_end;
else
  noise_end = 0;
end

if isfield(opt,'recode_rules')
    recode_rules = opt.recode_rules;
else
    recode_rules = [];
end

if isfield(opt,'combinations') && ~isempty(opt.combinations)
    combinations = opt.combinations;
    if isfield(opt,'combo_calc') && ~isempty(opt.combo_calc)
        calc = opt.combo_calc;
    else
        calc = 'weighted';
    end
    if isfield(opt,'comboeventcodes') && ~isempty(opt.comboeventcodes)
        comboeventcodes = opt.comboeventcodes;
    else
        comboeventcodes = [];
    end
else
    combinations = [];
end

if isfield(opt,'additioncombos') && ~isempty(opt.additioncombos)
    additioncombos = opt.additioncombos;
    if isfield(opt,'additioneventcodes') && ~isempty(opt.additioneventcodes)
        additioneventcodes = opt.additioneventcodes;
    else
        additioneventcodes = [];
    end
else
    additioncombos = [];
end

%% Load Data / Sort events across files

% Load and Convert Data into TimeSurfer

for i = 1:length(eeg_file)
  fprintf('Loading EEG file: %s\n',fullfile(eeg_path,eeg_file{i}));
  % CONVERT THE INPUT DATA TO TIMESURFER FORMAT
	[jk1 jk2 filetype jk3] = fileparts(eeg_file{i});
	if strcmp(filetype,'.mat')
%     if strcmp(who('-file',fullfile(eeg_path,eeg_file{i}),'epoch_data'),'epoch_data')
      epoch_data_array(i) = getfield(load(fullfile(eeg_path,eeg_file{i}),'epoch_data'),'epoch_data');
%       if length(tmp_epoch.epochs) > 1
%       else        
%       end
%     else
%       if ~isfield(opt,'mat_dim_order'), opt.mat_dim_order = ''; end
%       if ~isfield(opt,'mat_sfreq'), 		opt.mat_sfreq = ''; 		end
%       if ~isfield(opt,'mat_tlims'), 		opt.mat_tlims = ''; 		end
%       epoch_data_array(i)  = ts_iEEG_mat2epochs(fullfile(eeg_path,eeg_file{i}),opt.mat_dim_order,opt.mat_sfreq,opt.mat_tlims);
%     end
	else
  	epoch_data_array(i)  = ts_iEEG_eeg2epoch (fullfile(eeg_path,eeg_file{i}),...
                                           'badchanfile',pa.badchanfile,'channamefile',pa.channamefile,...
                                           'noise_start',noise_start,'noise_end',noise_end,...
                                           'recode_rules',recode_rules,...
                                           'timelimits',timelimits);                                        
	end
  if i == 1                                                                                                       
    sensor_info = epoch_data_array(i).sensor_info;                                            
  else                                                                                                          
    % CHECK CONSISTENCY OF CHANNELS ACROSS FILES
    if ~isempty(setdiff({sensor_info.label},{epoch_data_array(i).sensor_info.label}))                                         
      error('ERROR: Channels between files are not consistent!!!\n');
    end
  end
  % SETUP THE CHANNEL INFORMATION FOR WHEN PREPROCESSING WITH FIELDTRIP
  pa.channels = 1:length(epoch_data_array(i).sensor_info);           
  pa.chantype = epoch_data_array(i).sensor_info(1).typestring;
  % CHANGE EVENT CODES IF NECESSARY
  for j=1:length(epoch_data_array(i).epochs)         
    if ~strcmp(class(epoch_data_array(i).epochs(j).data),precision)
      if strcmpi(precision,'double')
        epoch_data_array(i).epochs(j).data = double(epoch_data_array(i).epochs(j).data);
      else
        epoch_data_array(i).epochs(j).data = single(epoch_data_array(i).epochs(j).data);
      end
    end
    if ~isempty(neweventcodes{i})                                                                                  
      fprintf('Changing event code %s',num2str(epoch_data_array(i).epochs(j).event_code));
      if find(oldeventcodes{i} == epoch_data_array(i).epochs(j).event_code)
       epoch_data_array(i).epochs(j).event_code  = neweventcodes{i}(find(oldeventcodes{i} == epoch_data_array(i).epochs(j).event_code));
      end
      fprintf(' to %s for data from %s.\n',num2str(epoch_data_array(i).epochs(j).event_code),eeg_file{i});      
    end    
  end
end

if length(epoch_data_array) > 1
  fprintf('Combining across eeg files...\n');
  all_epoch_data = ts_combine_data(epoch_data_array);
else
  all_epoch_data = epoch_data_array(1);
end
clear epoch_data_array

[badtrials{1:length(all_epoch_data.epochs)}] = deal([]);

%% Stimulus Delay
%  SIMPLY MODIFIES THE TIME VECTOR ACROSS THE CONDITIONS SO THAT THE 0 TIME
%  POINT NOW CORRESPONDS TO THE STIMULUS RATHER THAN TRIGGER ONSET

if ~isempty(stim_delay)
  fprintf('Introducing stimulus delay of %s seconds.\n',num2str(stim_delay))
  sampling_rate = all_epoch_data.sfreq;
  time_steps    = 1/sampling_rate;
  for j = 1:length(all_epoch_data.epochs)
    orig_time     = all_epoch_data.epochs(j).time;
    prestim_samp  = find(orig_time <= stim_delay,1,'last')-1;
    poststim_samp = length(orig_time) - (prestim_samp + 1);
    time_min      = -(prestim_samp*time_steps);
    time_max      =  (poststim_samp*time_steps);
    all_epoch_data.epochs(j).time = [];
    all_epoch_data.epochs(j).time = [time_min:time_steps:time_max];
    fprintf('Condition %2s, prestimulus period: %s sec / poststimulus period: %s sec.\n',num2str(j),num2str(time_min),num2str(time_max));
  end
  clear orig_time prestim_samp poststim_samp time_min time_max
end
clear sampling_rate time_steps


%% Remove Data as indicated by Reject File
%  WHEN REPROCESSING WITH A KNOWN REJECTION CRITERIA

if ~isempty(reject_file)
  load(reject_file);
  reject_flag = 1;
  [all_epoch_data.sensor_info(reject_data.badchans).badchan] = deal(1);
  badchannels = find([all_epoch_data.sensor_info.badchan]==1);
  for i = 1:length(all_epoch_data.epochs)
    old_data = all_epoch_data.epochs(i).data;
    all_epoch_data.epochs(i).data = [];
    goodtrials = setdiff(1:size(old_data,3),reject_data.badtrials{i});
    badtrials{i} = reject_data.badtrials{i};
    all_epoch_data.epochs(i).data = old_data(:,:,goodtrials);
    all_epoch_data.epochs(i).num_trials = size(all_epoch_data.epochs(i).data,3);
    all_epoch_data.epochs(i).num_rejects.manual = all_epoch_data.epochs(i).num_rejects.manual + ...
                                                  length(badtrials{i});
    if ~isempty(badchannels), all_epoch_data.epochs(i).data(badchannels,:,:) = nan; end
    clear old_data goodtrials
  end
end

%% Manual Rejection - will only happen if no reject file - BEFORE ICA
if isempty(reject_file) && reject_flag && ...
   (strcmpi(opt.reject,'yes') || strcmpi(opt.reject,'before_ica') || strcmpi(opt.reject,'both'))   % Manual Rejection
  rej_args =  mmil_parms2args(reject);
  [all_epoch_data,badchannels,badtrials] = ts_vis_reject(all_epoch_data,rej_args{:});
  % Save badchannel and trials as matlab file
  reject_data = [];
  reject_data.badchans = badchannels;
  reject_data.badtrials= [];
  for i = 1:length(all_epoch_data.epochs)
    reject_data.badtrials{i} = badtrials{i};
  end
elseif reject_flag
    % make a blank rejection info structure for use after ica
  reject_data.badchans = [];
  for i = 1:length(all_epoch_data.epochs)
    reject_data.badtrials{i} = [];
  end
end% if reject_flag

%% ICA Rejection

if auto_ica
  all_epoch_data = ts_autoICA(all_epoch_data,'ICA_ref_chan',refchan,'rescale',ica_rescale,'notch',true);
end

if manual_ica
  all_epoch_data = ts_manualICA(all_epoch_data,'ncomponents',[],'sorttrials',ica_sorttrials,'allconditions',ica_allconds,'rescale',ica_rescale,'notch',true);
end

%% Manual Rejection - will only happen if no reject file - AFTER ICA

if isempty(reject_file) && reject_flag && ...
   (strcmpi(opt.reject,'after_ica') || strcmpi(opt.reject,'both'))   % Manual Rejection
  rej_args =  mmil_parms2args(reject);
  [all_epoch_data,badchannels,badtrials] = ts_vis_reject(all_epoch_data,rej_args{:});
  % Save badchannel and trials as matlab file
  reject_data.badchans = unique(reject_data.badchans,badchannels);
  for i = 1:length(all_epoch_data.epochs)
    reject_data.badtrials{i} = unique(reject_data.badtrials{i}, badtrials{i});
  end
  [rej_path, reject_file] = fileparts(f_name);
  reject_file = [reject_file '.reject_data.mat'];
  reject_file = fullfile(f_path,rej_path,reject_file);
  save(reject_file,'reject_data');
end % if reject_flag

%% Save out the rejection info into a mat file

if reject_flag && ~exist(reject_file,'file')
    [rej_path, reject_file] = fileparts(f_name);
    reject_file = [reject_file '.reject_data.mat'];
    reject_file = fullfile(f_path,rej_path,reject_file);
    save(reject_file,'reject_data');
end

%% Save Rejection info to text file

badchannels = find([all_epoch_data.sensor_info.badchan] == 1);

[rej_path reject_file] = fileparts(f_name);
if reject_flag
  reject_file = [reject_file '.rej.rejection_info.txt'];
else
  reject_file = [reject_file '.rejection_info.txt'];
end
reject_file = fullfile (f_path,rej_path,reject_file);

fprintf('Saving rejections info: %s\n',reject_file);

fid = fopen(reject_file,'w+');
fprintf(fid,'Rejected Channels\n\n');
if ~isempty(badchannels)
  for i = 1:length(badchannels)
    fprintf(fid,'%s\n',sensor_info(badchannels(i)).label);
  end
else
  fprintf(fid,'No bad channels.\n');
end
fprintf(fid,'\n');
fprintf(fid,'Trial Info\n\n');
fprintf(fid,'         \t Good\t    Rejected     \t Bad\n');
fprintf(fid,'Condition\tTrials\tEEG File\tManual\tTrials\n\n');

for j = 1:length(all_epoch_data.epochs)
 fprintf(fid,'%-9s\t%-6s\t%-8s\t%-6s\t%s\n',num2str(j),num2str(all_epoch_data.epochs(j).num_trials),...
                                                   num2str(all_epoch_data.epochs(j).num_rejects.eeg),...
                                                   num2str(all_epoch_data.epochs(j).num_rejects.manual),...
                                                   num2str(badtrials{j}));
end

fclose(fid);

%% Save TimeSurfer Epoch Data
% keep track of original number of conditions for averaging
original_conds  = length(all_epoch_data.epochs);
if save_epochs      % Save each epoch in its own mat file due to file size limits
  if ~isempty(combinations)
    all_epoch_data         = ts_combine_conditions(all_epoch_data,...
                                                   'combinations',combinations,...
                                                   'calc',calc,...
                                                   'reference','events',...
                                                   'neweventcodes',comboeventcodes);
  end
  epoch_data.num_sensors = all_epoch_data.num_sensors;
  epoch_data.sensor_info = all_epoch_data.sensor_info;
  epoch_data.coor_trans  = all_epoch_data.coor_trans;
  epoch_data.sfreq       = all_epoch_data.sfreq;
  epoch_data.noise       = all_epoch_data.noise;
  try epoch_data.parms = opt; end
  for j = 1:length(all_epoch_data.epochs)                                                              % Split off individual epochs for saving
    epoch_data.epochs(1) = all_epoch_data.epochs(j);
    fprintf ('Saving event %s epoch: ',num2str(all_epoch_data.epochs(j).event_code));
    e_name = f_name;
    if manual_ica || auto_ica,    e_name = [e_name '.ica']; end   
    if reject_flag, e_name = [e_name '.rej']; end
    try
      if strcmpi(epoch_tag,'cond')
          e_name = [e_name '.epoch.cond' num2str(j) '.mat'];
      end
      if strcmpi(epoch_tag,'event')
          e_name = [e_name '.epoch.event' num2str(all_epoch_data.epochs(j).event_code) '.mat'];
      end
      fprintf('%s\n',e_name);
      try epoch_data.parms.filename{1} = fullfile(f_path,e_name); end
      save(fullfile(f_path,e_name),'epoch_data');
      filename = {filename{:},fullfile(f_path,e_name)};
    catch
      error('ERROR: Unable to save epoch_data: %s\n',e_name);
    end;
    epoch_data = rmfield(epoch_data,'epochs');                                                            % Clear this epoch
  end
  clear epoch_data
end

%% Compute Averages and Save

if save_avgs                 % Average each condition.
  % only process and average the original conditions
  for j = 1:original_conds
    fprintf ('Preprocessing and averaging event: %s\n',num2str(all_epoch_data.epochs(j).event_code));  
    % convert to fieldtrip first
    all_epoch_data.epochs(j).data = double(all_epoch_data.epochs(j).data);
    ft_epochs = ts_epoch2fieldtrip(all_epoch_data,'condition',j,'dimord','chan_time','chantype',pa.chantype);
    prepro.keeptrials         = 'no';
    % use FieldTrip preprocessing to apply filters, ect.
    warning off all
    fieldtrip_data{j}= timelockanalysis(prepro, ft_epochs);
    warning on all
    clear ft_epochs
  end
  % convert back to TimeSurfer
  avg_data = ts_fieldtrip2data (fieldtrip_data,'averages',all_epoch_data);              
  % Create the combinations for the averages
  if ~isempty(combinations)  
    avg_data = ts_combine_conditions(avg_data,...
                                     'combinations',combinations,...
                                     'calc',calc,...
                                     'reference','events',...
                                     'neweventcodes',comboeventcodes);
  end
  if ~isempty(additioncombos)
    avg_data = ts_iEEG_addavg(avg_data,...
                              'combinations',additioncombos,...
                              'neweventcodes',additioneventcodes);
  end

  try                                                                                                      % Saving data
    a_name = [a_name '.avg.mat'];
    fprintf('Saving avg_data to: %s\n',a_name);
    avg_data.parms = opt; 
    avg_data.parms.filename{1} = fullfile(f_path,a_name);
    save(fullfile(f_path,a_name),'avg_data');
    filename = {filename{:},fullfile(f_path,a_name)};
  catch
    fprintf('ERROR: Unable to save avg_data.\n');
  end
end

% if nargout == 1
%   varargout{1} = filename;
% end

