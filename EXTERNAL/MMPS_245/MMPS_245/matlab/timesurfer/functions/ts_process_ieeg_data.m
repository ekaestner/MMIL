function [filenames,alldata] = ts_process_ieeg_data(datafile,varargin)
% Outputs:
%   alldata: epoch_data or avg_data
% Inputs:
%   datafile: cell array of files containing Neuroscan (*.eeg),
%             epoch_data (*.mat) data
%   Parameters:
%     pre_xxx: FT or TS preprocessing parameters to run BEFORE artifact rejection
%     xxx: includes FT or TS preprocessing parameters to run after artifact rejection
%     save_epochs_raw_flag: save epochs after loading *.eeg file
%     saveepochs_flag: save epochs after artifact rejection
%     saveepochs_post_flag: save epochs after the preprocessing following
%                           artifact rejection
%     saveaverages_flag:    save averages after everything
global filenames
parms = mmil_args2parms(varargin,...
						{'prefix','proc',[],...
             'rootoutdir',pwd,[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...             
						 'precision','double',{'single','double'},...
             'conditionkey','',[],...
             'conditionfile','',[],...
             'datafile',[],[],...
             'filename',[],[],...
             'badchanfile',[],[],...
             'channamefile',[],[],...
             'rejectfile',[],[],...
             'channels',[],[],...
             'chantype','all',[],...
             'saveepochs_flag',1,{0,1},...
             'saveepochs_raw_flag',0,{0,1},...
             'saveepochs_post_flag',0,{0,1},...
             'saveaverages_flag',1,{0,1},...
             'oldeventcodes',[],[],...
             'neweventcodes',[],[],...
             'stim_delay',[],[],...
             'timelimits',[],[],...
             'ICA_auto_flag',0,{0,1},...
             'ICA_manual_flag',0,{0,1},...
             'ICA_ref_chan','EOG061',[],...
             'ICA_chantype','all',[],...
             'ICA_maxsteps',20,[],...
             'ICA_ntrials',5,[],...
             'ICA_ncomponents',80,[],...
             'ICA_rescale_flag',1,{0,1},...
             'ICA_sorttrials',0,{0,1},...
             'visualreject',0,{0,1},...
             'reject_method','summary',[],...
             'reject_metric','var',[],...
             'method',[],[],...
             'metric',[],[],...
             'noise_start',-.08,[],...
             'noise_end',0,[],...
             'recode_rules',[],[],...
             'combinations',[],[],...
             'epoch_combinations',[],[],...
             'epoch_neweventcodes',[],[],...
             'average_combinations',[],[],...
             'average_neweventcodes',[],[],...
             'combo_calc','weighted',[],...
             'comboeventcodes',[],[],...
             'additioncombos',[],[],...
             'additioneventcodes',[],[],...
             'write_avg_flag',0,{0,1},...
             'bandpass_flag',false,[false true],...
             'bandpass_low_cf',0.2,[],...
             'bandpass_low_tb',0.4,[],...
             'bandpass_high_cf',50,[],...
             'bandpass_high_tb',10,[],...
             'dsfact',1,[],...
             'detrend_flag',false,[false true],...
             'baseline_flag',true,[false true],...
             'baseline_start',-Inf,[-Inf,Inf],...
             'baseline_end',0,[-Inf,Inf],...     
             'keeptrials',[],[],...
             'combinations',[],[],...
             'comboeventcodes',[],[],...
             'calc','weighted',{'weighted','avg','sum'},...
             'cfg',[],[],...
             'feedback','none',{'non','gui','dial','textbar','text','textcr','textnl','no','none','yes'},...
             'lpfilter','no',{'yes','no'},...
             'hpfilter','no',{'yes','no'},...
             'bpfilter','no',{'yes','no'},...
             'bsfilter','no',{'yes','no'},...
             'lnfilter','no',{'yes','no'},...
             'dftfilter','no',{'yes','no'},...
             'medianfilter','no',{'yes','no'},...
             'lpfreq',[],[],...
             'hpfreq',[],[],...
             'bpfreq',[],[],...
             'lnfreq',50,[],...
             'dftfreq',[50 100 150],[],...              
             'lpfiltord',6,[],...
             'hpfiltord',6,[],...
             'bpfiltord',4,[],...
             'bsfiltord',4,[],...
             'lnfiltord',4,[],...
             'lpfilttype','but',{'but','fir'},...
             'hpfilttype','but',{'but','fir'},...
             'bpfilttype','but',{'but','fir'},...
             'bsfilttype','but',{'but','fir'},...
             'lpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'hpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'bpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'bsfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'medianfiltord',9,[],...
             'blc','no',{'yes', 'no'},...
             'blcwindow',[-inf 0],[],...
             'detrend','no',{'yes','no'},...     
             'polyremoval','no',{'yes','no'},...
             'latency',[],[],...
             'polyorder',2,[],...
             'hilbert','no',{'no','abs','complex','real','imag','absreal','absimag','angle'},...
             'rectify','no',{'yes','no'},...
             'precision','double',{'single','double'},...
             'reref','no',{'yes','no'},...
             'refchannel',[],[],...
             'implicitref',[],[],...
             'montage','no',[],...          
             'fieldtrip_flag',0,{0,1},...
             'threshold_flag',0,{0,1},...
             'events',[],[],...
             'overwrite',0,{0,1},...
             'coordinatefile',[],[],...
             'coordfield','coords',[],...
             'keepbadchans',0,{0,1},...
						},false);
          
parms     = backcompatible(parms,varargin{:});      
preparms = get_preparms(varargin{:});
% postparms = get_postparms(varargin{:});

filenames = {};
rejflag = parms.ICA_auto_flag || parms.ICA_manual_flag || parms.threshold_flag || (parms.visualreject && isempty(parms.rejectfile));
epoflag = parms.saveepochs_flag || parms.saveepochs_raw_flag || parms.saveepochs_post_flag;
preflag = parms.saveaverages_flag || (epoflag && ischar(parms.keeptrials) && strcmp(parms.keeptrials,'yes'));
%preflag = parms.saveaverages_flag || ~parms.fieldtrip_flag || (epoflag && ischar(parms.keeptrials) && strcmp(parms.keeptrials,'yes'));

if ~iscell(datafile) && exist(datafile,'file'), datafile = {datafile}; end
if ~isempty(parms.conditionkey) && exist(parms.conditionkey,'file')
  ts_makecondkey(parms.conditionkey);
end

if exist(parms.badchanfile,'file')
  fid=fopen(parms.badchanfile,'rt');
  j = 1;
  while (~feof(fid))
   badchans{j} = fgetl(fid);
   j=j+1;
  end
  fclose(fid);
else
  badchans = [];
end
if length(badchans)==1 && isequal(badchans{1},-1)
  badchans = [];
end

%% Load and Convert Data into TimeSurfer format
for f = 1:length(datafile)
  thisfile = datafile{f};
  if parms.verbose
    mmil_logstr(parms,'Loading data file: %s\n',thisfile);
%     fprintf('Loading data file: %s\n',thisfile);
  end
	[jk1 jk2 filetype] = fileparts(thisfile);
    switch filetype
        case '.mat'
%     if strcmp(who('-file',thisfile,'epoch_data'),'epoch_data')
          if f == 1
            data_arr(f) = getfield(load(thisfile,'epoch_data'),'epoch_data');
          else
            data_arr(f) = orderfields(getfield(load(thisfile,'epoch_data'),'epoch_data'),data_arr(1));
          end
          if ~isempty(badchans)
            [jnk idx] = match_str(badchans,{data_arr(f).sensor_info.label});
            [data_arr(f).sensor_info(idx).badchan] = deal(1);
          end
        case {'.eeg','.vhdr','.set'}
          data_arr(f)  = ts_iEEG_eeg2epoch (thisfile,...
                                           'badchanfile',parms.badchanfile,'channamefile',parms.channamefile,...
                                           'noise_start',parms.noise_start,'noise_end',parms.noise_end,...
                                           'recode_rules',parms.recode_rules,...
                                           'timelimits',parms.timelimits);  
        case '.cnt'
            data_arr(f) = ts_loadcnt(thisfile);
        case '.avg'
          tmp         = ts_matrix2avg(loadavg(thisfile));
          tmp.epochs  = rmfield(tmp.averages,'std_dev');
          data_arr(f) = rmfield(tmp,'averages');
          clear tmp
    end
  clear jk1 jk2 filetype
  % Added by JSS
  if (any([data_arr(f).sensor_info.badchan]) || ~isempty(parms.channels) || ~strcmpi(parms.chantype,'all')) && parms.keepbadchans == 0
    data_arr(f) = orderfields(ts_data_selection(data_arr(f),'removebadchans',1,'channels',parms.channels,'chantype',parms.chantype),data_arr(f));
  end
    
  if f == 1                                         
    sensor_info = data_arr(f).sensor_info;
  else                                                                                                          
    % CHECK CONSISTENCY OF CHANNELS ACROSS FILES
    if ~isempty(setdiff({sensor_info.label},{data_arr(f).sensor_info.label}))                                         
      error('ERROR: Channels between files are not consistent!!!\n');
    end
  end
  % SETUP THE CHANNEL INFORMATION FOR WHEN PREPROCESSING WITH FIELDTRIP
  parms.channels = 1:length(data_arr(f).sensor_info);
  parms.chantype = data_arr(f).sensor_info(1).typestring;
  % CHANGE EVENT CODES IF NECESSARY
  for j=1:length(data_arr(f).epochs)   
    if ~strcmp(class(data_arr(f).epochs(j).data),parms.precision)
      if strcmpi(parms.precision,'double')
        data_arr(f).epochs(j).data = double(data_arr(f).epochs(j).data);
      else
        data_arr(f).epochs(j).data = single(data_arr(f).epochs(j).data);
      end
    end
    if ~isempty(parms.neweventcodes) && ~(iscell(parms.neweventcodes) && isempty(parms.neweventcodes{f}))
      mmil_logstr(parms,'Changing event code %s',num2str(data_arr(f).epochs(j).event_code));
%       fprintf('Changing event code %s',num2str(data_arr(f).epochs(j).event_code));
      if find(parms.oldeventcodes{f} == data_arr(f).epochs(j).event_code)
        data_arr(f).epochs(j).event_code  = parms.neweventcodes{f}(find(parms.oldeventcodes{f} == data_arr(f).epochs(j).event_code));
      end
      mmil_logstr(parms,' to %s for data from %s.\n',num2str(data_arr(f).epochs(j).event_code),thisfile);      
%       fprintf(' to %s for data from %s.\n',num2str(data_arr(f).epochs(j).event_code),thisfile);
    end    
  end
end

% combine across raw data files
if length(data_arr) > 1
  % this double for loop is a temporary addition to deal with a anomalous
  % data set
  for f = 1:length(data_arr)
      if isfield(data_arr(f).epochs,'name')
        data_arr(f).epochs = rmfield(data_arr(f).epochs,'name');
      end
  end
  mmil_logstr(parms,'%s: combining data across files...\n',mfilename);
%   fprintf('%s: combining data across files...\n',mfilename);
  alldata = ts_combine_data(data_arr);
else
  alldata = data_arr(1);
end
clear data_arr
if ~isempty(parms.events)
  alldata = ts_data_selection(alldata,'events',parms.events);
end
[badtrials{1:length(alldata.epochs)}] = deal([]);

%% Stimulus delay (Rajan's code copied from ts_iEEG_ProcessEEG)
%  SIMPLY MODIFIES THE TIME VECTOR ACROSS THE CONDITIONS SO THAT THE 0 TIME
%  POINT NOW CORRESPONDS TO THE STIMULUS RATHER THAN TRIGGER ONSET
if ~isempty(parms.stim_delay)
  mmil_logstr(parms,'Introducing stimulus delay of %s seconds.\n',num2str(parms.stim_delay))
%   fprintf('Introducing stimulus delay of %s seconds.\n',num2str(parms.stim_delay))
  sampling_rate = alldata.sfreq;
  time_steps    = 1/sampling_rate;
  for j = 1:length(alldata.epochs)
    orig_time     = alldata.epochs(j).time;
    prestim_samp  = find(orig_time <= parms.stim_delay,1,'last')-1;
    poststim_samp = length(orig_time) - (prestim_samp + 1);
    time_min      = -(prestim_samp*time_steps);
    time_max      =  (poststim_samp*time_steps);
    alldata.epochs(j).time = [];
    alldata.epochs(j).time = [time_min:time_steps:time_max];
    mmil_logstr(parms,'Condition %2s, prestimulus period: %s sec / poststimulus period: %s sec.\n',num2str(j),num2str(time_min),num2str(time_max));
%     fprintf('Condition %2s, prestimulus period: %s sec / poststimulus period: %s sec.\n',num2str(j),num2str(time_min),num2str(time_max));
  end
  clear orig_time prestim_samp poststim_samp time_min time_max sampling_rate time_steps
end

%% Reject file
if ~isempty(parms.rejectfile) && exist(parms.rejectfile) % && ~parms.visualreject
  mmil_logstr(parms,'loading reject file: %s\n',parms.rejectfile);
%   fprintf('loading reject file: %s\n',parms.rejectfile);
    load(parms.rejectfile);
    parms.visualreject = 0;
  goodchans = sort(setdiff(1:alldata.num_sensors,reject_data.badchans));
  mmil_logstr(parms,'%s: removing %g channels\n',mfilename,length(reject_data.badchans))
%   fprintf('%s: removing %g channels\n',mfilename,length(reject_data.badchans))
  alldata.sensor_info = alldata.sensor_info(goodchans);
  alldata.num_sensors = length(goodchans);
  for k = 1:length(reject_data.badtrials)
    if isfield(reject_data,'event_code')
      if ~ismember(reject_data.event_code(k),[alldata.epochs.event_code])
        continue;
      end
      kk = find(reject_data.event_code(k)==[alldata.epochs.event_code]);
      goodtrials = sort(setdiff(1:alldata.epochs(kk).num_trials,reject_data.badtrials{k}));
      alldata.epochs(kk).data        = alldata.epochs(kk).data(goodchans,:,goodtrials);
      alldata.epochs(kk).num_rejects.manual = alldata.epochs(kk).num_rejects.manual +...
             (alldata.epochs(kk).num_trials - length(goodtrials));
      alldata.epochs(kk).num_trials  = length(goodtrials);
      mmil_logstr(parms,'%s: removing %g trials from event %g\n',mfilename,length(reject_data.badtrials{k}),reject_data.event_code(k));      
%       fprintf('%s: removing %g trials from event %g\n',mfilename,length(reject_data.badtrials{k}),reject_data.event_code(k));      
    else
      goodtrials = sort(setdiff(1:alldata.epochs(k).num_trials,reject_data.badtrials{k}));
      alldata.epochs(k).data        = alldata.epochs(k).data(goodchans,:,goodtrials);
      alldata.epochs(k).num_rejects.manual = alldata.epochs(k).num_rejects.manual +...
             (alldata.epochs(k).num_trials - length(goodtrials));
      alldata.epochs(k).num_trials  = length(goodtrials);
    end
  end
end

%% Preprocessing (filter, baseline correction, detrending, downsampling)
if ~isempty(preparms)
  args = mmil_parms2args(rmfield(parms,{'saveaverages_flag','saveepochs_flag','combinations'}));
  alldata = ts_preproc(alldata,args{:});
  clear args;
end

%% Artifact Rejection (visual & ICA)
if rejflag
  args = mmil_parms2args(rmfield(parms,'rejectfile'));
  [alldata reject_data] = ts_reject (alldata,args{:});
  clear args;
else 
  reject_data = [];
end

%% Combine epochs (concatenation)
if ~isempty(parms.epoch_combinations) && epoflag
  if isempty(parms.epoch_neweventcodes) && ~isempty(parms.comboeventcodes)
    parms.epoch_neweventcodes = parms.comboeventcodes;
  end
  if length(parms.epoch_combinations) == length(parms.epoch_neweventcodes)  
    alldata = ts_combine_conditions(alldata,'combinations',parms.epoch_combinations,'neweventcodes',parms.epoch_neweventcodes,'calc',parms.combo_calc);
  else
    mmil_logstr(parms,'%s: warning: skipping combinations. # epoch combos does not equal # codes.\n',mfilename);
%     fprintf('%s: warning: skipping combinations. # epoch combos does not equal # codes.\n',mfilename);
  end
end

%% Add electrode coordinates
if ~isempty(parms.coordinatefile)
  if exist(parms.coordinatefile,'file')
    alldata = ts_read_coords(parms.coordinatefile,alldata,parms.coordfield);
  else
    mmil_logstr(parms,'%s: coordinate file not found: %s\n',mfilename,parms.coordinatefile);
%     fprintf('%s: coordinate file not found: %s\n',mfilename,parms.coordinatefile);
  end
end

if parms.saveepochs_flag
  save_epochs('');
end

%% Postprocessing (filter, baseline correction, detrending, downsampling)
% if ~isempty(postparms)
%   args = mmil_parms2args(rmfield(postparms,{'saveaverages_flag','saveepochs_flag','combinations'}));
%   alldata = ts_preproc(alldata,args{:});
%   clear args;
% end
% if parms.saveepochs_post_flag
%   save_epochs('.post');
% end
args = mmil_parms2args(rmfield(parms,{'saveaverages_flag','saveepochs_flag','combinations'}));
if preflag
  alldata = ts_preproc(alldata,args{:});
end
clear args;
if parms.saveepochs_post_flag
  save_epochs('.post');
end

%% Data selection (conditions, channels, trials)
% alldata = ts_data_selection(alldata,varargin{:});
alldata = ts_data_selection(alldata,'events',parms.events);
if parms.saveepochs_raw_flag
  save_epochs('.raw');
end

%% Save results
if ~isempty(reject_data)
  save_reject;
end
% Save averages
if parms.saveaverages_flag
  avg_data = ts_trials2avg(alldata,'stdev','data');
  if nargout < 2, clear alldata; end
  %% Combine averages (addition, subtraction, ...)
  if ~isempty(parms.average_combinations)
    if isempty(parms.average_neweventcodes) && ~isempty(parms.comboeventcodes)
      parms.average_neweventcodes = parms.comboeventcodes;
    end
    if length(parms.average_combinations) == length(parms.average_neweventcodes)
      avg_data = ts_combine_conditions(avg_data,'combinations',parms.average_combinations,'neweventcodes',parms.average_neweventcodes,'calc',parms.combo_calc);
    else
      mmil_logstr(parms,'%s: warning: skipping combinations. # average combos does not equal # codes.\n',mfilename);
%       fprintf('%s: warning: skipping combinations. # average combos does not equal # codes.\n',mfilename);
    end
  end
    if isempty(parms.filename)
      parms.filename = sprintf('%s/%s_averages_cond%g.mat',parms.rootoutdir,parms.prefix,c);
    else
      if ~iscell(parms.filename), parms.filename = {parms.filename}; end
      [fpath fname] = fileparts(parms.filename{1});
      filename = sprintf('%s/%s.mat',fpath,fname);
      if ~isempty(strfind(filename,'epoch'))
        filename = strrep(filename,'epoch','avg');
      elseif isempty(strfind(filename,'avg'))
        filename = strrep(filename,'.mat','.avg.mat');
      end
    end
    filenames         = {filenames{:} filename};
    avg_data.parms    = parms;
    avg_data.parms.filename = {filename};   
    mmil_logstr(parms,'%s: saving file: %s\n',mfilename,filename)    
%     fprintf('%s: saving file: %s\n',mfilename,filename)    
    check_path(filename);
    if ~exist(avg_data.parms.filename{1},'file') || parms.overwrite
      save(avg_data.parms.filename{1},'avg_data');
    else
      mmil_logstr(parms,'%s: not overwriting %s\n',mfilename,avg_data.parms.filename{1});
%       fprintf('%s: not overwriting %s\n',mfilename,avg_data.parms.filename{1});
    end
    clear avg_data filename
end

% optionally save Neuroscan averages (*.avg)
if parms.write_avg_flag 
  % TODO: add this functionality
end
  function save_reject
  if isempty(parms.comboeventcodes)
    if iscell(parms.average_neweventcodes)
      parms.comboeventcodes = [parms.comboeventcodes parms.average_neweventcodes{:}];
    else
      parms.comboeventcodes = [parms.comboeventcodes parms.average_neweventcodes];
    end
    if iscell(parms.epoch_neweventcodes)
      parms.comboeventcodes = [parms.comboeventcodes parms.epoch_neweventcodes{:}];
    else
      parms.comboeventcodes = [parms.comboeventcodes parms.epoch_neweventcodes];
    end  
    parms.comboeventcodes = unique(parms.comboeventcodes);
  end  
  % Save reject data
  if ~isempty(reject_data)
    % save info to mat file
    if ~isempty(filenames)
      [fpath fname] = fileparts(filenames{1});
      fname = fname(1:find(fname == '.',1,'last'));
      fname = [fname 'reject_data.mat'];
    else
      fpath = parms.rootoutdir;
      fname = sprintf('%s_reject_data.mat',parms.prefix);
    end
    check_path(fullfile(fpath,fname));
    save(fullfile(fpath,fname),'reject_data');
    filenames = {filenames{:} fullfile(fpath,fname)};

    % save info to text file (copied from Rajan's code in ts_iEEG_ProcessEEG)
    [datatype datafield dataparam] = ts_object_info(alldata);
    fname = strrep(fname,'.mat','.txt');
    fid = fopen(fullfile(fpath,fname),'w+');
    badchans = find([alldata.sensor_info.badchan]);  
    fprintf(fid,'Rejected Channels\n\n');
    if ~isempty(badchans)
      for j = 1:length(badchans)
        fprintf(fid,'%s\n',alldata.sensor_info(badchans(j)).label);
      end
    else
      fprintf(fid,'No bad channels.\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Trial Info\n\n');
    fprintf(fid,'         \t Good\t    Rejected     \t Bad\n');
    fprintf(fid,'Condition\tTrials\tEEG File\tManual\tTrials\n\n');
    for j = 1:length(setdiff([alldata.epochs.event_code],parms.comboeventcodes))
      try
        if ~ismember(alldata.(datafield)(j).event_code,[reject_data.event_code])
          continue;
        end
      end
      if j > length(reject_data.badtrials), break; end
     fprintf(fid,'%-9s\t%-6s\t%-8s\t%-6s\t%s\n',num2str(j),num2str(alldata.(datafield)(j).num_trials),...
                                                       num2str(alldata.(datafield)(j).num_rejects.eeg),...
                                                       num2str(alldata.(datafield)(j).num_rejects.manual),...
                                                       num2str(reject_data.badtrials{j}));
    end
    fclose(fid);  
  end
  end

  function save_epochs(tag)
    % Save epochs
    if parms.saveepochs_flag
      for c = 1:length(alldata.epochs)
        epoch_data = rmfield(alldata,'epochs');
        epoch_data.epochs(1) = alldata.epochs(c);
        if isempty(parms.filename)
          filename = sprintf('%s/%s_epochs%s_cond%02i.mat',parms.rootoutdir,parms.prefix,tag,epoch_data.epochs.event_code);
        else
          if ~iscell(parms.filename), parms.filename = {parms.filename}; end
          [fpath fname] = fileparts(parms.filename{1});
          filename = sprintf('%s/%s%s_event%02i.mat',fpath,fname,tag,epoch_data.epochs.event_code);
        end
        filenames         = {filenames{:} filename};
        epoch_data.parms  = parms;
        epoch_data.parms.filename = {filename};
        mmil_logstr(parms,'%s: saving file %g of %g: %s\n',mfilename,c,length(alldata.epochs),filename)        
%         fprintf('%s: saving file %g of %g: %s\n',mfilename,c,length(alldata.epochs),filename)        
        check_path(filename);
        if ~exist(epoch_data.parms.filename{1},'file') || parms.overwrite
          save(epoch_data.parms.filename{1},'epoch_data');
        else
          mmil_logstr(parms,'%s: not overwriting %s\n',mfilename,epoch_data.parms.filename{1});          
%           fprintf('%s: not overwriting %s\n',mfilename,epoch_data.parms.filename{1});
        end
        clear epoch_data filename
      end
    end
    clear c
  end
end

%% Subfunctions
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
     'auto_ica','',{'yes','no'},...
     'refchan','',[],...
     'manual_ica','',{'yes','no'},...
     'ica_rescale','',{'yes','no'},...
     'ica_sorttrials','',{'yes','no'},...
     'reject','',{'yes','after_ica','before_ica','both','no'},...     
     'reject_method','',[],...
     'reject_metric','',[],...
     'combinations','',[],...
     'saveepochs','',{'yes','no'},...
     'saveavgs','',{'yes','no'},...
     },false);
   
if ~isempty(opt.auto_ica) && strcmp(opt.auto_ica,'yes')
  parms.ICA_auto_flag = 1;
end 
if ~isempty(opt.manual_ica) && strcmp(opt.manual_ica,'yes')
  parms.ICA_manual_flag = 1;
end
if ~isempty(opt.ica_rescale) && strcmp(opt.ica_rescale,'no')
  parms.ICA_rescale_flag = 0;
end
if ~isempty(opt.ica_sorttrials) && strcmp(opt.ica_sorttrials,'yes')
  parms.ICA_sorttrials = 1;
end
if ~isempty(opt.refchan) && ischar(opt.refchan)
  parms.ICA_ref_chan = opt.refchan;
end
if ~isempty(opt.reject) && ~strcmp(opt.reject,'no')
  parms.visualreject = 1;
end
if ~isempty(opt.reject_method) && ischar(opt.reject_method)
  parms.method = opt.reject_method;
end
if ~isempty(opt.reject_metric) && ischar(opt.reject_metric)
  parms.metric = opt.reject_metric;
end
if ~isempty(opt.combinations)
  parms.average_combinations = opt.combinations;
end
if ~isempty(opt.saveepochs) && ischar(opt.saveepochs) && strcmp(opt.saveepochs,'no')
  parms.saveepochs_flag = 0;
end
if ~isempty(opt.saveavgs) && ischar(opt.saveavgs) && strcmp(opt.saveavgs,'no')
  parms.saveaverages_flag = 0;
end
end

function parms = get_preparms(varargin)
parms = [];
for i = 1:length(varargin)
  if ischar(varargin{i}) && mod(i,2) && ~isempty(findstr(varargin{i},'pre_'))
    fld = strrep(varargin{i},'pre_','');
    try parms.(fld) = varargin{i+1}; end
  end
end
end

% function parms = get_postparms(varargin)
% parms = [];
% for i = 1:length(varargin)
%   if ischar(varargin{i}) && mod(i,2) && ~isempty(findstr(varargin{i},'post_'))
%     fld = strrep(varargin{i},'post_','');
%     try parms.(fld) = varargin{i+1}; end
%   end
% end
% end

function check_path(filename)
 pathstr = fileparts(filename);
if ~exist(pathstr,'dir')
  mmil_logstr(parms,'making output directory: %s\n',pathstr);
%   fprintf('making output directory: %s\n',pathstr);
  unix(['mkdir -p ' pathstr]);
end      
end
