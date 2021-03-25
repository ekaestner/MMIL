function process_iEEG_data (setupfile,varargin)
%
% Usage: process_iEEG_data('setup.m');
%
% Script for running through the iEEG processing stream. Loads the
% designated setup file and runs throught the processing tream
%
% Required Input:
%
%  setupfile - a string containing name of the setup.m file
%  containing the settings to run the stream.  The file should be a MATLAB
%  M-file created in a similar format to that of the sample_setup.m file as
%  shown below
%
% %%% Experiment Settings
% 
% subjects_dir   = '';       % the directory for the experiment
% conditionkey   = '';       % name of the condition key .csv file
% rootoutdir     = '';       % the directory for saving data and figures
% 
% %% What to do?
% 
% processEEGfile = 'no';     % process the input data
% TFAnalysis     = 'no';     % time frequency analysis
% TFCombos       = 'no';     % combine time frequency conditions
% statistics     = 'no';     % perform statistical analysis on avg data
% TFstatistics   = 'no';     % time-frequency statistics - this is not yet a working option
% plotAverages   = 'no';     % plot averages
% plotTF         = 'no';     % plot time frequency data
% 
% %% Data Handling 
% 
% eegfiles       = '';      % input file(s) - placed in a 'eeg' directory of the subjects_dir
% timelimits     = [0 0];   % time limits for the trials when inputing BrainVision data
% prefix         = '';      % prefix to use for new files
% badchanfile    = [];      % name of the bad channel text file
% channamefile   = [];      % name of new channel name text file
% rejectfile     = [];      % name of manual rejection .mat file
% oldeventcodes  = [];      % event condes in the input data set
% neweventcodes  = [];      % what to change the old event codes to
% saveavgs       = 'no';    % save out averaged data
% saveepochs     = 'no';    % save out epoch data (data with trials)
% 
% %% Preprocessing Options
% 
% stim_delay       = 0;     % delay to the stimulus
% lpfilter         = 'no';  % low pass filter
% hpfilter         = 'no';  % high pass filter
% lnfilter         = 'no';  % line filter (notch filter)
% bpfilter         = 'no';  % band pass filter
% lpfreq           = 0;     % frequency of low pass filter in Hz
% hpfreq           = 0;     % frequency of high pass filter in Hz
% lnfreq           = 0;     % frequency of line filter in Hz
% bpfreq           = [0 0]; % the high and low pass settings for the band pass filter
% blc              = 'no';  % baseline correct
% blcwindow        = [0 0]; % time period of baseline (in seconds)
% detrend          = 'no';  % detrend the data
% visual_reject    = 'no';  % perform visual rejection: 'before_ica' or 'after_ica' or 'yes'
% auto_ica         = 'no';  % perform automatic ica rejection
% refchan          = '';    % reference channel for auto ica
% manual_ica       = 'no';  % perofrm manual ica rejection
% ica_rescale      = 'no';  % rescale the data after removing bad components
% ica_allconds     = 'yes'; % combine all conditions prior to component analysis
% ica_sorttrials   = 'yes'; % sort trials from high to low variance before displaying components
% combinations     = {};    % list of combination rules - creates weghted avgs
% comboeventcodes  = [];    % list of event codes to give the new conditions
% additioncombos   = {};    % list of combination rules as above but adds two averages
% additioneventcodes = [];  % list of event codes to give the additionconditions
% recode_rules     = {};    % recoding rules (create new conditions based on order of events)
% 
% %% Frequency Analysis
% 
% tf.events          = [];  % the events to perform frequency analysis on
% tf.baseline_events = [];  % the events to use for the baseline - this should be only the original events (not the combinations or recode events)
% tf.fa_method       =  ''; % the FieldTrip method to use for frequency analysis
% tf.foi             = [];  % list of frequencies (in Hz) to process
% tf.width           = [];  % width of the wavelets for each of the frequencies
% tf.toi             = [];  % the time over which to perform the analysis
% tf.reject          = 'no'; % perform rejection in the time-frequency domain
% tf.bl_threshold    = 10;   % the threshold of the baseline rejection
% tf.trial_threshold = 20;   % the threshold to rejecting trials
% tf.chan_threshold  = 90;   % the threshold to rejecting channels
% tf.reject_exclude  = 'bychannel'; % how to exclude trials: 'bychannel','inall','groupchannels'
% tf.baseline_file   = ''; % the name of a previously calculated baseline data set
% 
% %% TF Combos
% 
% tf_combos.combinations  = {}; % list of rules for creating time-frequency combinations
% tf_combos.neweventcodes = []; % the new event codes for the combinations
% 
% %% Statistics
% 
% stat.events           = {[]};  % list of statistical tests you want to perform ([] means test between the events within the brackets)
% stat.method           = '';    % the FieldTrip statistical method
% stat.statistic        = '';    % type "help timelockstatistics" for more information
% stat.correctm         = '';    % of the rest of these options
% stat.alpha            = [];    % consult the manual (or wiki) for some basic setups
% stat.tail             = 0;
% stat.ivar             = 1;
% stat.uvar             = [];
% stat.wvar             = [];
% stat.clusterstatistic = '';
% stat.clusterthreshold = '';
% stat.clusteralpha     = [];
% stat.clustercrtival   = [];
% stat.clustertail      = [];
% stat.avgoverchan      = [];
% stat.avgovertime      = [];
% stat.numrandomization = 0;
% stat.minnbchan        = [];
% stat.neighbours       = [];
% stat.design           = [];   % for the most part the design is automatically made
% 
% %% TF Statistics - not yet working
% 
% tf_stat.events           = {[]};  % list of statistical tests you want to perform ([] means test between the events within the brackets)
% tf_stat.method           = '';    % the FieldTrip statistical method
% tf_stat.statistic        = '';    % type "help freqstatistics" for more information
% tf_stat.correctm         = '';    % of the rest of these options
% tf_stat.alpha            = [];    % consult the manual (or wiki) for some basic setups
% tf_stat.tail             = 0;
% tf_stat.ivar             = 1;
% tf_stat.uvar             = [];
% tf_stat.wvar             = [];
% tf_stat.clusterstatistic = '';
% tf_stat.clusterthreshold = '';
% tf_stat.clusteralpha     = [];
% tf_stat.clustercrtival   = [];
% tf_stat.clustertail      = [];
% tf_stat.avgoverchan      = [];
% tf_stat.avgovertime      = [];
% tf_stat.numrandomization = 0;
% tf_stat.minnbchan        = [];
% tf_stat.neighbours       = [];
% tf_stat.design           = [];
% 
% %% Common Plotting Options
% 
% clearscreen   = 'no';      % do you want to clear each of the figures off the screen
% savefigures   = 'no';      % what type of file to save the figures as: 'jpg','eps',ect..
% resolution    = 300;       % resolution of the saved image
% showaxis      = 'yes';     % show the axis on the plot
% channels      = {};        % limit the plotting to this list of channels 
% showbadchans  = 'no';      % show the data of bad channels
% layoutfile    = [];        % use this layout file (.asc or .lay) for the figures
% notes         = '';        % add these notes to the plots
% whichfigs     = [];        % plot only these pages
% 
% %% Plot Waveforms
%  
% plotavg.events     = {[]};  % which events to plot (those within [] are plotted together)
% plotavg.xlim       = [];    % time limits of the plot
% plotavg.ylim       = [];    % amplitude limits of the plot
% plotavg.statistics = {''};  % list of statistic files to use for each of the above set of plots
% plotavg.alpha      = [];    % the p-value to consider significant
% 
% %% Plot TF
% 
% plotwaves.events      = [];  % which events to plot
% plotwaves.xlim        = [];  % time limits of the plot
% plotwaves.ylim        = [];  % frequency limits of the plot
% plotwaves.zlim        = [];  % the color limits of the plot
% plotwaves.freqcorr    = 'no'; % when plotting power correct for drop off at higher frequencies ignored if any of the next two are set to yes
% plotwaves.plotzscores = 'no'; % plot z-scores
% plotwaves.baseline    = 'no'; % what type of baseline plot ('relative','relchange','absolute') ignored if plotting z-scores
% plotwaves.timedsfact  = 1; % downsample across time by this factor (1 = no downsampling)
% plotwaves.freqdsfact  = 1; % downsample across frequencies by this factor (1 = no downsampling)
% plotwaves.statistics  = {''};  % list of statistic files to use for each of the above plots
%
% Created:  10/24/07 by Rajan Patel
% Last Mod: 09/13/08 by Jason Sherfey
% Last Mod: 09/15/12 by Don Hagler
%

% See also: sample_setup.m

% TO DO:
%
% - Double check the ability to run TF statistics
% - Check the passing of appropriate TF stats file names to TF Plotting
%   function
% - Laminar processing
% 
% Revision 01-Aug-2008 by Jason Sherfey
% Added mkdir for rootoutdir
%
% Revision 05-Aug-2008 by Jason Sherfey
% Added rootindir option: tf.freq_path  = fullfile(rootindir,'freq_analysis');

%% Check input and Load the Setup File (run it as a script)

% if ~mmil_check_nargs(nargin, 1, 2), return; end;

if ~exist(setupfile,'file')
  fprintf('ERROR: Specified setup file not found: %s\n\n',setupfile);
  return;
end

[path,settingsfile,ext] = fileparts(setupfile);      % get rid of .m extension
fprintf('Setup file: %s\n',[path settingsfile ext]);
if strcmp(ext,'.xls')
	ts_read_excel_setup(setupfile,1);
	settingsfile = ['auto_' settingsfile];
end

eval([path settingsfile]);                               % load settings

%% Check for Processing Flags and Subject Settings

if ~exist('subjects_dir','var'),     subjects_dir   = pwd;          end  
if ~exist('rootindir','var'),        rootindir      = subjects_dir; end
if ~exist('rootoutdir','var'),       rootoutdir     = rootindir;    end 
if ~exist('conditionkey','var'),     conditionkey   = [];           end
if ~exist('processEEGfile','var'),   processEEGfile = 'no';         end
%if ~exist('combinations','var'),     combinations   = 'no';         end
if ~exist('combinations','var'),     combinations   = [];         end 		% changed on 11-Nov-2008 by Jason Sherfey
if ~exist('TFAnalysis','var'),       TFAnalysis     = 'no';         end
if ~exist('TFCombos','var'),         TFCombos       = 'no';         end
if ~exist('statistics','var'),       statistics     = 'no';         end
if ~exist('TFstatistics','var'),      TFstatistics   = 'no';        end
if ~exist('plotAverages','var'),     plotAverages   = 'no';         end
if ~exist('plotTF','var'),           plotTF         = 'no';         end
if ~exist('prefix','var'),           prefix         = 'iEEGdata';   end
if ~exist('badchanfile','var'),      badchanfile    = [];           end
if ~exist('channamefile','var'),     channamefile   = [];           end
if ~exist('rejectfile','var'),       rejectfile     = [];           end


%% Initiate arguments that will be passed to the processing functions

eeg_args     = {};
stat_args    = {};
combo_args   = {};
tf_args      = {};
tf_combo_args= {};
plotavg_args = {};
plottf_args  = {};

%% Verify Files / Data Handling - make sure everything exists
%  ADD the necessary options to the argument list by using the subfunction
%  addargs

if ~exist(subjects_dir,'dir'),
  error('Subject''s directory "%s" cannot be found.\n',subjects_dir);
end

if ~exist(rootindir,'dir'),
  error('Input directory "%s" cannot be found.\n',rootindir);
end

if ~exist(rootoutdir,'dir'),
    fprintf('Making output directory: %s.\n',rootoutdir);
    unix(['mkdir -p ' rootoutdir]);
end

if ~isempty(conditionkey)
  conditionkey = fullfile(subjects_dir,conditionkey);
  if ~exist(conditionkey,'file')
    error('Cannot find condition key: %s.\n',conditionkey);
    conditionkey = [];
  end
  fprintf('Making condition key from %s.\n',conditionkey);
  ts_iEEG_makecondkey(conditionkey);
end

if ~isempty(badchanfile)
  badchanfile = fullfile(subjects_dir,badchanfile);
  if ~exist(badchanfile,'file')
    error('Cannot find bad channel file: %s.\n',badchanfile);
    badchanfile = [];
  end
end

if ~isempty(channamefile)
  channamefile = fullfile(subjects_dir,channamefile);
  if ~exist(channamefile,'file')
    error('Cannot find channel name file: %s.\n',channamefile);
    channamefile = [];
  end
end

if ~isempty(rejectfile)
  rejectfile = fullfile(rootindir,'matfiles',rejectfile);
  if ~exist(rejectfile,'file')
    error('Cannot find rejection file: %s.\n',rejectfile);
  end
end

if exist('saveavgs','var') && (strcmpi(saveavgs,'yes') || strcmpi(saveavgs,'no'))
  eeg_args = addargs(eeg_args,'saveavgs',saveavgs);
else
  eeg_args = addargs(eeg_args,'saveavgs','yes');
end
if exist('saveepochs','var') && (strcmpi(saveepochs,'yes') || strcmpi(saveepochs,'no') || ...
                                 strcmpi(saveepochs,'events') || strcmpi(saveepochs,'conditions'))
  eeg_args = addargs(eeg_args,'saveepochs',saveepochs);
  %% THIS DETERMINES HOW EPOCHS ARE SAVED OUT
  %% THE END OF THE EPOCH FILE NAMES WILL HAVE
  %% .event#.epoch.mat OR .cond#.epoch.mat
  %% THE CONDITION NUMBER DESIGNATION IS FOR OLD VERSIONS REALLY
  switch saveepochs
    case {'yes','no','events'}
      epoch_tag = 'event';
    case {'conditions'}
      epoch_tag = 'cond';
  end
else
  eeg_args = addargs(eeg_args,'saveepochs','yes');
  epoch_tag = 'event';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse Pre-processing options

% HERE WE GO THROUGH EACH OF THE PREPROCESSING OPTIONS AND WHERE NECESSARY
% ADD IT TO THE LIST OF ARGUMENTS BEING PASSED TO THE PREPROCESSING
% FUNCTION AS WELL AS THE STATS FUNCITON.
%
% THE FILE NAME IS ALSO DETERMINED HERE.  THIS FILE NAME MUST
% CORRESPOND TO THE FILE NAME THAT IS CREATED IN ts_iEEG_ProcessEEG.  
%
% TO DO: MAKE THE FILE NAME HERE AND THEN JUST USE THE FILENAME CREATED AND
% SEND THAT TO THE PROCESSING FUNCTION AS AN OUTPUT FILE NAME OR PREFIX
% SO THAT IT CAN BE CHANGED EASILY IN HERE RATHER THAN FIGURING IT OUT IN
% EACH FUNCTION AS WELL AS IN THIS INITIAL FUNCTION.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nil prefix] = fileparts(prefix);
% THIS IS THE STARTING FILE NAME
avg_name   = prefix;
epoch_name = prefix;

if ~exist('lpfilter','var'),      lpfilter = 'no';      end;
if ~exist('hpfilter','var'),      hpfilter = 'no';      end;
if ~exist('lnfilter','var'),      lnfilter = 'no';      end;
if ~exist('bpfilter','var'),      bpfilter = 'no';      end;
if ~exist('blc','var'),           blc      = 'no';      end;
if ~exist('detrend','var'),       detrend  = 'no';      end;
if ~exist('visual_reject','var'), visual_reject = 'no'; end;
if ~exist('auto_ica','var'),      auto_ica = 'no';      end;
if ~exist('manual_ica','var'),    manual_ica = 'no';    end;
if ~exist('refchan','var'),       refchan  = [];        end;
if ~exist('ica_rescale','var'),   ica_rescale = 'no';   end;
if ~exist('ica_sorttrials','var'),ica_sorttrials = 'yes';end;
if ~exist('ica_allconds','var'),  ica_allconds = 'yes'; end;

if exist('stim_delay','var') && ~isempty(stim_delay)
    eeg_args = addargs(eeg_args,'stim_delay',stim_delay);
end

if exist('recode_rules','var') && ~isempty(recode_rules)
    eeg_args = addargs(eeg_args,'recode_rules',recode_rules);
end

if strcmpi(lpfilter,'yes')
 if ~exist('lpfreq','var') || isempty(lpfreq), 
   fprintf('Setting low pass frequency to 100 Hz.\n');
   lpfreq = 100; 
 end  
 eeg_args  = addargs(eeg_args,'lpfilter','yes');
 eeg_args  = addargs(eeg_args,'lpfreq',lpfreq);
 stat_args = addargs(stat_args,'lpfilter','yes');
 stat_args = addargs(stat_args,'lpfreq',lpfreq);
 avg_name = [avg_name '.lp' strrep(num2str(lpfreq),'.','')];
end

if strcmpi(hpfilter,'yes')
 if ~exist('hpfreq','var') || isempty(hpfreq), 
   fprinf('Setting high pass frequency to 0 Hz.\n');
   hpfreq = 0; 
 end
 eeg_args  = addargs(eeg_args,'hpfilter','yes');
 eeg_args  = addargs(eeg_args,'hpfreq',hpfreq);
 stat_args = addargs(stat_args,'hpfilter','yes');
 stat_args = addargs(stat_args,'hpfreq',hpfreq); 
 avg_name= [avg_name '.hp' strrep(num2str(hpfreq),'.','')];
end

if strcmpi(bpfilter,'yes')
  if ~exist('bpfreq','var') || isempty(bpfreq), 
    fprintf('Setting band pass filter high pass to 0 Hz and low pass to 100 Hz.\n');
    bpfreq = [0 100]; 
  end
 eeg_args  = addargs(eeg_args,'bpfilter','yes');
 eeg_args  = addargs(eeg_args,'bpfreq',bpfreq);
 stat_args = addargs(stat_args,'bpfilter','yes');
 stat_args = addargs(stat_args,'bpfreq',bpfreq);     
 avg_name = [avg_name '.bp' strrep(num2str(bpfreq(1)),'.','') '-' strrep(num2str(bpfreq(2)),'.','')];
end

if strcmpi(lnfilter,'yes')
  if ~exist('lnfreq','var') || isempty(lnfreq)
    fprintf('Setting line frequency to 60 Hz.\n');
    lnfreq = 60; 
  end
 eeg_args  = addargs(eeg_args,'lnfilter','yes');
 eeg_args  = addargs(eeg_args,'lnfreq',lnfreq);
 stat_args = addargs(stat_args,'lnfilter','yes');
 stat_args = addargs(stat_args,'lnfreq',lnfreq);   
 avg_name = [avg_name '.ln' strrep(num2str(lnfreq),'.','')];
end

if strcmpi(blc,'yes')
 eeg_args  = addargs(eeg_args,'blc','yes');
 stat_args = addargs(stat_args,'blc','yes');  
 avg_name = [avg_name '.bc'];
end

if ~exist('blcwindow','var') || isempty(blcwindow) || length(blcwindow) ~= 2
    fprintf('Setting baseline correction window to -0.080 sec - 0 sec.\n');
    blcwindow = [-.08 0];
end
eeg_args  = addargs(eeg_args,'blcwindow',blcwindow);
stat_args = addargs(stat_args,'blcwindow',blcwindow);

if exist('noise_start','var') && ~isempty(noise_start)
  eeg_args = addargs(eeg_args,'noise_start',noise_start);
else
  eeg_args = addargs(eeg_args,'noise_start',blcwindow(1));
end

if exist('noise_end','var') && ~isempty(noise_end)
  eeg_args = addargs(eeg_args,'noise_end',noise_end);
else
  eeg_args = addargs(eeg_args,'noise_end',blcwindow(2));
end

if exist('timelimits','var') && ~isempty(timelimits)
   eeg_args = addargs(eeg_args,'timelimits',timelimits);
end

if strcmpi(detrend,'yes')
 eeg_args = addargs(eeg_args,'detrend','yes');
 stat_args = addargs(stat_args,'detrend','yes');
 avg_name = [avg_name '.dt'];
end

if strcmpi(auto_ica,'yes') || strcmpi(manual_ica,'yes')
  eeg_args = addargs(eeg_args,'auto_ica',auto_ica);
  eeg_args = addargs(eeg_args,'manual_ica',manual_ica);
  eeg_args = addargs(eeg_args,'refchan',refchan);
  eeg_args = addargs(eeg_args,'ica_rescale',ica_rescale);
  eeg_args = addargs(eeg_args,'ica_allconds',ica_allconds);
  eeg_args = addargs(eeg_args,'ica_sorttrials',ica_sorttrials);
  avg_name = [avg_name '.ica'];
  epoch_name = [epoch_name '.ica'];
end

if strcmpi(visual_reject,'yes') ||...
   strcmpi(visual_reject,'before_ica') ||...     
   strcmpi(visual_reject,'after_ica') ||...
   strcmpi(visual_reject,'both') ||...
   ~isempty(rejectfile)
 eeg_args = addargs(eeg_args,'reject',visual_reject);
 avg_name = [avg_name '.rej'];
 epoch_name = [epoch_name '.rej'];
end

if exist('combinations','var') && ~isempty(combinations)
    eeg_args = addargs(eeg_args,'combinations',combinations);
    if exist('comboeventcodes','var')
        eeg_args = addargs(eeg_args,'comboeventcodes',comboeventcodes);
    end
    if exist('combo_calc','var')
        eeg_args = addargs(eeg_args,'combo_calc',combo_calc);
    end       
end

if exist('additioncombos','var') && ~isempty(additioncombos)
    eeg_args = addargs(eeg_args,'additioncombos',additioncombos);
    if exist('additioneventcodes','var')
        eeg_args = addargs(eeg_args,'additioneventcodes',additioneventcodes);
    end
end

epoch_name = [epoch_name '.epoch'];

%% Individualize Settings and Run
%  MAKE SURE THE ARGUMENTS AGREE WITH THE INPUT PARAMATERS OF EACH FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process EEG

if strcmpi(processEEGfile,'yes')
  eeg_args = addargs(eeg_args,'eegpath', fullfile(subjects_dir,'eeg'));
  eeg_args = addargs(eeg_args,'savefile', fullfile(rootoutdir,'matfiles',prefix));
  eeg_args = addargs(eeg_args,'badchanfile', badchanfile);
  eeg_args = addargs(eeg_args,'channamefile', channamefile);
  eeg_args = addargs(eeg_args,'rejectfile', rejectfile);
  if exist('eegfiles','var') && ~isempty(eegfiles)
    eeg_args = addargs(eeg_args, 'eegfile', eegfiles);
  else
    error('Must specify name of eeg file.');
  end
  if exist('oldeventcodes','var') && ~isempty(oldeventcodes)
    eeg_args = addargs(eeg_args, 'oldeventcodes', oldeventcodes);
  end
  if exist('neweventcodes','var') && ~isempty(neweventcodes)
    eeg_args = addargs(eeg_args, 'neweventcodes', neweventcodes);
  end

	if ~exist('mat_dim_order','var'), mat_dim_order = ''; end
	if ~exist('mat_sfreq','var'),			mat_sfreq  = ''; 		end
	if ~exist('mat_tlims','var'), 		mat_tlims = ''; 		end
	eeg_args = addargs(eeg_args,'mat_dim_order',mat_dim_order);
	eeg_args = addargs(eeg_args,'mat_sfreq',mat_sfreq);
	eeg_args = addargs(eeg_args,'mat_tlims',mat_tlims);
	
  ts_iEEG_ProcessEEG(eeg_args{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combinations  
%  THIS IS THE OLD WAY... LEFT IN TO WORK WITH OLD SETUP FILES

if strcmpi(combinations,'yes')
  combo_args = addargs(combo_args,'f_path',fullfile(rootoutdir,'matfiles'));
  if exist('saveavgs','var') && strcmpi(saveavgs,'no')
  else
     combo_args = addargs(combo_args,'avg_name',[avg_name '.avg.mat']);
  end
  if exist('saveepochs','var') && strcmpi(saveepochs,'no')
  elseif exist('saveepochs','var') && (strcmpi(saveepochs,'conditions'))
     combo_args = addargs(combo_args,'epoch_name',[epoch_name '.cond']);
  else   
     combo_args = addargs(combo_args,'epoch_name',[epoch_name '.event']);
  end
  if exist('combos','var')
    if isfield(combos,'combinations') && ~isempty(combos.combinations)
      combo_args = addargs(combo_args,'combinations',combos.combinations);
    else
      error('Must provide a list of combinations to produce.');
    end
    if isfield(combos,'neweventcodes') && ~isempty(combos.neweventcodes)
      combo_args = addargs(combo_args,'neweventcodes',combos.neweventcodes);
    end
    if ~isfield(combos,'calc') || isempty(combos.calc)
      combos.calc = 'weighted';
    end
    combo_args = addargs(combo_args,'calc',combos.calc);      
    if isfield(combos,'reference')
      combo_args = addargs(combo_args,'reference',combos.reference);
    end
  else
    error('No options for creating combinations was provided.\n Please refer to the sample setup file.');
  end
  ts_iEEG_combos(combo_args{:});
end 

% OLD - UPDATE THE FILE NAME TO THE COMBINATION FILE NAME FOR PLOTTING IF A
% COMBINATION AVG_DATA FILE EXISTS

if exist(fullfile(rootoutdir,'matfiles',sprintf('%s.combo.avg.mat',avg_name)),'file');
  avg_name = [avg_name '.combo'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TF Analysis

if strcmpi(TFAnalysis,'yes')   
  if exist('calcwaves','var') % backward compatability
    tf = calcwaves; 
    tf.conditions = tf.file_tags; 
    tf = rmfield(tf,'file_tags');
  end
  epoch_files = {};
	
	% check for varargin input of tf.width & tf.gwidth
	if nargin>1 
		if exist('varargin') && ~isfield(tf,'width') && ~isfield(tf,'gwidth')
			tf_in = mmil_args2parms(varargin,...
                       {...
                         'width',[],[],...
												 'gwidth',[],[]...
                       },...
                       false);
			tf.width 	= tf_in.width;
			tf.gwidth = tf_in.gwidth;
		end											 
	end

  if isfield(tf,'events') && ~isempty(tf.events)
      for i = 1:length(tf.events)
          epoch_files{i} = sprintf('%s.event%s.mat',epoch_name,num2str(tf.events(i)));
      end
      %% Add any baseline files if it's not included in the events for
      %% processing
      if isfield(tf,'baseline_events') && ~isempty(tf.baseline_events)
          add_bl_files = setdiff(tf.baseline_events,tf.events);
          for i = 1:length(add_bl_files)
              epoch_files{end+1} = sprintf('%s.event%s.mat',epoch_name,num2str(tf.baseline_events(i)));
          end
      end
      tf = rmfield(tf,'events');
  elseif isfield(tf,'conditions') && ~isempty(tf.conditions)
    for i = 1:length(tf.conditions)
      epoch_files{i} = sprintf('%s.cond%s.mat',epoch_name,num2str(tf.conditions(i)));
    end
    tf = rmfield(tf,'conditions');
  else
    error('Must provide either eventcodes or condition numbers for frequency analysis.');    
  end
 
  if isfield(tf,'baseline_file') && ~isempty(tf.baseline_file)
		tf.baseline_file = sprintf('matfiles/freq_analysis/%s',tf.baseline_file);
		if exist(fullfile(rootindir,tf.baseline_file),'file')
			tf.baseline_file = fullfile(rootindir,tf.baseline_file);
		elseif exist(fullfile(rootoutdir,tf.baseline_file),'file')
			tf.baseline_file = fullfile(rootoutdir,tf.baseline_file);
		elseif exist(fullfile(subjects_dir,tf.baseline_file),'file')
			tf.baseline_file = fullfile(subjects_dir,tf.baseline_file);
		end
  end
  if ~isfield(tf,'baseline_start') || isempty(tf.baseline_start)
    tf.baseline_start = blcwindow(1);
  end
  if ~isfield(tf,'baseline_end') || isempty(tf.baseline_end)
    tf.baseline_end   = blcwindow(end);
  end
  tf.f_path     = fullfile(rootindir,'matfiles');
  tf.freq_path  = fullfile(rootoutdir,'matfiles','freq_analysis');
  tf.f_name     = epoch_files;
  tf_args_opt   = mmil_parms2args(tf);
	ts_timefreq_analysis(tf_args_opt{:},tf_args{:});
%  ts_iEEG_FreqAnalysis(tf_args_opt{:},tf_args{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Echo Names

avg_name   = [avg_name '.avg.mat'];
epoch_name = [epoch_name '.' epoch_tag];
tf_name    = epoch_name(1:strfind(epoch_name,'epoch')-2);
try 
 tf_name    = [tf_name ...
              '.foi' strrep(num2str(tf.foi(1)),'.','') '-' strrep(num2str(tf.foi(end)),'.','') ...
              '.toi' strrep(num2str(tf.toi(1)),'.','') '-' strrep(num2str(tf.toi(end)),'.','')];
end

if exist('tf','var') && isfield(tf,'reject') && strcmpi(tf.reject,'yes')
    if isfield(tf,'exclude')
        switch tf.exclude
            case 'inall'
                tfrejtag = '.tfrejall';
            case 'groupchannels'
                tfrejtag = '.tfrejgroups';
            case 'bychannel'
                tfrejtag = '.tfrejbychan';
        end
    else
        tfrejtag = '.tfrejbychan';
    end
    if isfield(tf,'trial_threshold')
        tf_name = [tf_name tfrejtag num2str(tf.trial_threshold)];
    else
        tf_name = [tf_name tfrejtag '50'];
    end
    if isfield(tf,'bl_threshold')
      tf_name = [tf_name '.blthresh' num2str(tf.bl_threshold)];
    else
      tf_name = [tf_name '.blthresh10'];
    end
	  if isfield(tf,'chan_threshold')
			tf_name = [tf_name '.chanthresh' num2str(tf.chan_threshold)];
		else
			tf_name = [tf_name '.chanthresh80'];
		end
end
bl_name    = [tf_name '.baseline.mat'];
tf_name    = [tf_name '.timefreq.' epoch_tag];

fprintf('\nFiles(s) will have these names if generated:\n\n');
fprintf('avg_data file: %s\n',avg_name);
fprintf('epoch_data file: %s\n',epoch_name);
fprintf('timefreq_data file: %s\n',tf_name);
fprintf('baseline_data file: %s\n\n',bl_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TimeFrequency Combinations

if strcmpi(TFCombos,'yes')
  if ~exist('tf_combos','var')
      error('Must provide TFCombos variables as tf_combos structure.');
  end
  tf_combo_args = addargs(tf_combo_args,'f_path',fullfile(rootoutdir,'matfiles','freq_analysis'));
  if exist('saveepochs','var') && (strcmpi(saveepochs,'conditions'))
     tf_combo_args = addargs(tf_combo_args,'timefreq_name',[tf_name '.cond']);
  else   
     tf_combo_args = addargs(tf_combo_args,'timefreq_name',[tf_name '.event']);
  end
  if exist('tf_combos','var')
    if isfield(tf_combos,'combinations') && ~isempty(tf_combos.combinations)
      tf_combo_args = addargs(tf_combo_args,'combinations',tf_combos.combinations);
    else
      error('Must provide a list of combinations to produce.');
    end
    if isfield(tf_combos,'neweventcodes') && ~isempty(tf_combos.neweventcodes)
      tf_combo_args = addargs(tf_combo_args,'neweventcodes',tf_combos.neweventcodes);
    end
    if ~isfield(tf_combos,'calc') || isempty(tf_combos.calc)
      tf_combos.calc = 'weighted';
    end
    tf_combo_args = addargs(tf_combo_args,'calc',tf_combos.calc);      
    if isfield(tf_combos,'reference')
      tf_combo_args = addargs(tf_combo_args,'reference',tf_combos.reference);
    end
  else
    error('No options for creating combinations was provided.\n Please refer to the sample setup file.');
  end
  ts_iEEG_combos(tf_combo_args{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics

if strcmpi(statistics,'yes')
  if ~exist('stat','var'), error('Must provided stat options to run statistics.'); end
  if isfield(stat,'file_tags')  % backward compatability
    stat.conditions = stat.file_tags;
    stat = rmfield(stat,'file_tags');
  end
  stat.f_path   = fullfile(rootindir,'matfiles');
  stat.stat_path = fullfile(rootoutdir,'matfiles','stats');
  stat.f_name   = {};
  if isfield(stat,'events') && ~isempty(stat.events)
    statevents = stat.events;
    stat.events= [];
  elseif isfield(stat,'conditions') && ~isempty(stat.conditions)  
    statconditions = stat.conditions;
    stat.conditions= [];
  else
    error('Must provide either eventcodes or condition numbers for statistical testing.');
  end
  if exist('statevents','var')
    if ~iscell(statevents), statevents = {statevents}; end
    for j = 1:length(statevents)
			if isempty(statevents{j}), continue; end;
      for i = 1:length(statevents{j})
				
        stat.f_name{i} = sprintf('%s%s.mat',epoch_name,num2str(statevents{j}(i)));
      end
      stat.events   = statevents{j};
      stat_args_opt = mmil_parms2args(stat);
			ts_stat_analysis(stat_args_opt{:},stat_args{:});
%      ts_iEEG_StatAnalysis(stat_args_opt{:},stat_args{:});
      stat.f_name = {};
      stat.events = [];
    end
  else
    if ~iscell(statconditions), statsconditions = {statconditions}; end
    for j = 1:length(statconditions)
			if isempty(statconditions{j}), continue; end;
      for i = 1:length(statconditions{j})
        stat.f_name{i} = sprintf('%s%s.mat',epoch_name,num2str(statconditions{j}(i)));
      end
      stat.conditions = statconditions{j};
      stat_args_opt = mmil_parms2args(stat);
			ts_stat_analysis(stat_args_opt{:},stat_args{:});
%      ts_iEEG_StatAnalysis(stat_args_opt{:},stat_args{:});
      stat.f_name = {};
      stat.conditions = [];
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TF Stats

if strcmpi(TFstatistics,'yes')
  %%% REMOVE ERROR STATEMENT WHEN READY!!!!
  warning('Time-Frequency Stats is still undergoing testing.');
  if ~exist('tf_stat','var'), error('Must provided tf_stat options to run statistics.'); end
  tf_stat.f_path   = fullfile(rootindir,'matfiles','freq_analysis');
  tf_stat.f_name   = {};
  tf_stat.stat_path = fullfile(rootoutdir,'matfiles','tf_stats');
  if isfield(tf_stat,'events') && ~isempty(tf_stat.events)
    tfstatevents = tf_stat.events;
    tf_stat.events= [];
  elseif isfield(tf_stat,'conditions') && ~isempty(tf_stat.conditions)  
    tfstatconditions = tf_stat.conditions;
    tf_stat.conditions= [];
  else
    error('Must provide either eventcodes or condition numbers for timefrequency statistical testing.');
  end
  if exist('tfstatevents','var')
    if ~iscell(tfstatevents), tfstatevents = {tfstatevents}; end
    for j = 1:length(tfstatevents)
			if isempty(tfstatevents{j}), continue; end;
      for i = 1:length(tfstatevents{j})
        tf_stat.f_name{i} = sprintf('%s%s.mat',tf_name,num2str(tfstatevents{j}(i)));
      end
      tf_stat.events   = tfstatevents{j};
      tf_stat_args_opt = mmil_parms2args(tf_stat);
			ts_stat_analysis(tf_stat_args_opt{:},stat_args{:});	
%      ts_iEEG_StatAnalysis(tf_stat_args_opt{:},stat_args{:});
      tf_stat.f_name = {};
      tf_stat.events = [];
    end
  else
    if ~iscell(tfstatconditions), tfstatconditions = {tfstatconditions}; end
    for j = 1:length(tfstatconditions)
			if isempty(tfstatconditions{j}), continue; end;
      for i = 1:length(tfstatconditions{j})
        tf_stat.f_name{i} = sprintf('%s%s.mat',tf_name,num2str(tfstatconditions{j}(i)));
      end
      tf_stat.conditions =tfstatconditions{j};
      tf_stat_args_opt = mmil_parms2args(tf_stat);
			ts_stat_analysis(tf_stat_args_opt{:},stat_args{:});			
%      ts_iEEG_StatAnalysis(tf_stat_args_opt{:},stat_args{:});
      tf_stat.f_name = {};
      tf_stat.conditions = [];
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('clearscreen','var')
  plotwaves.cls = clearscreen;
  plotavg.cls   = clearscreen;
end

if exist('savefigures','var')
  plotwaves.savefigs = savefigures;
  plotavg.savefigs   = savefigures;
end

if exist('showaxis','var')
  plotwaves.showaxis = showaxis;
  plotavg.showaxis   = showaxis;
end

if exist('channels','var') && ~isempty(channels)
  plotwaves.channels = channels;
  plotavg.channels   = channels;
end

if exist(fullfile(subjects_dir,'cond_key.mat'),'file')
  plotwaves.condkey = fullfile(subjects_dir,'cond_key.mat');
  plotavg.condkey   = fullfile(subjects_dir,'cond_key.mat');
end

if exist('layoutfile','var') && ~isempty(layoutfile)
  plotwaves.layoutfile = fullfile(subjects_dir,layoutfile);
  plotavg.layoutfile   = fullfile(subjects_dir,layoutfile);
end

if exist('notes','var') && ~isempty(notes)
    plotwaves.notes  = notes;
    plotavg.notes    = notes;
end

if exist('whichfigs','var') && ~isempty(whichfigs)
    plotwaves.whichfigs = whichfigs;
    plotavg.whichfigs = whichfigs;
end

if exist('showbadchans','var') && ~isempty(showbadchans)
    plotwaves.showbadchans = showbadchans;
    plotavg.showbadchans   = showbadchans;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot Averages

if strcmpi(plotAverages,'yes')
 plotavg.f_path = fullfile(rootindir,'matfiles');
 plotavg.f_name = avg_name;
 plotavg.badchanfile = badchanfile;
 plotavg.figpath = fullfile(rootoutdir,'images');
 if isfield(plotavg,'events')
   eventstoplot = plotavg.events;
   plotavg = rmfield(plotavg,'events');
 elseif isfield(plotavg,'conditions')
   conditionstoplot = plotavg.conditions;
   plotavg = rmfield(plotavg,'conditions');
 end
 if isfield(plotavg,'statistics') && ~isempty(plotavg.statistics)
   if ~iscell(plotavg.statistics), plotavg.statistics = {plotavg.statistics}; end
   if exist('eventstoplot','var')
     numberofplots = length(eventstoplot);
   elseif exist('conditionstoplot','var')
     numberofplots = length(conditionstoplot);
   end
   listofstatfiles = plotavg.statistics;
   if length(listofstatfiles) ~= numberofplots
     fprintf('NOTICE: The number of statistics files is not the same as the number of plots to create.\n');
     clear listofstatfiles
   else
     for i = 1:length(listofstatfiles)
       listofstatfiles{i} = fullfile(rootindir,'matfiles','stats',listofstatfiles{i});
     end
   end
   plotavg = rmfield(plotavg,'statistics');
 end
 if exist('eventstoplot','var')
   if ~iscell(eventstoplot), eventstoplot = {eventstoplot}; end
   for i = 1:length(eventstoplot)
     plotavg.eventvals = eventstoplot{i};
     if exist('listofstatfiles','var'), 
       plotavg.statistics = listofstatfiles{i};
     end
     plotavg_args = mmil_parms2args(plotavg);
     ts_iEEG_PlotAvg (plotavg_args{:});
     clear plotavg_args;
   end
 elseif exist('conditionstoplot','var')
   if ~iscell(conditionstoplot), conditionstoplot = {conditionstoplot}; end
   for i = 1:length(conditionstoplot)
     plotavg.conditions = conditionstoplot{i};
     if exist('listofstatfiles','var')
       plotavg.statistics = listofstatfiles{i}; 
     end
     plotavg_args = mmil_parms2args(plotavg);
     ts_iEEG_PlotAvg (plotavg_args{:});
     clear plotavg_args;
   end
 else
   fprintf('Nothing specified to plot.\n');
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot TF

if strcmpi(plotTF,'yes')
  if ~exist('plotwaves','var'), error('Must provide options to plot wavelets.'); end
  plotwaves.f_path = fullfile(rootindir,'matfiles','freq_analysis');
  plotwaves.badchanfile = badchanfile;
  plotwaves.figpath = fullfile(rootoutdir,'images','freq_analysis');
  if ~isfield(plotwaves,'baseline_start') || isempty(plotwaves.baseline_start)
      plotwaves.baseline_start = blcwindow(1);
  end
  if ~isfield(plotwaves,'baseline_end') || isempty(plotwaves.baseline_end)
      plotwaves.baseline_end = blcwindow(end);
  end
  plotwaves.baseline_file = bl_name;
 if isfield(plotwaves,'statistics') && ~isempty(plotwaves.statistics)
   if ~iscell(plotwaves.statistics), plotwaves.statistics = {plotwaves.statistics}; end
   listoftfstatfiles = plotwaves.statistics;
   for i = 1:length(listofstatfiles)
     listoftfstatfiles{i} = fullfile(rootindir,'matfiles','tf_stats',listoftfstatfiles{i});
   end
   plotwaves = rmfield(plotwaves,'statistics');
 end
  if isfield(plotwaves,'events') && ~isempty(plotwaves.events)
    for i = 1:length(plotwaves.events)
     plotwaves.f_name = sprintf('%s%s.mat',tf_name,num2str(plotwaves.events(i)));
     try
        plotwaves.statistics = listoftfstatfiles{i};
     end
     plottf_args = mmil_parms2args(plotwaves);
     ts_iEEG_PlotWavelets(plottf_args{:});
    end
  elseif isfield(plotwaves,'conditions') && ~isempty(plotwaves.conditions)
    for i = 1:length(plotwaves.conditions)
     plotwaves.f_name = sprintf('%s%s.mat',tf_name,num2str(plotwaves.conditions(i)));
     try
         plotwaves.statistics = listofstatfiles{i};
     end
     plottf_args = mmil_parms2args(plotwaves);
     ts_iEEG_PlotWavelets(plottf_args{:});
    end
  else
    fprintf('No event codes given to plot.\n');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SUB FUNCTIONS

function new_arg_list = addargs (old_arg_list,option,value)
% old_arg_list - the argument list is a cell array with 'option' value
% option       - the string name of the option
% value        - the value to associate with the option
  new_arg_list = old_arg_list;
  clear old_arg_list
  new_arg_list{end+1} = option;
  new_arg_list{end+1} = value;
