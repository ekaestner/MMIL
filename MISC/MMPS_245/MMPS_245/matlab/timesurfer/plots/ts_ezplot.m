function zlims = ts_ezplot(varargin)
% Generate multiplots and topoplots of time and time-frequency data + stats
%
% ts_ezplot(avg_data)      ==> multiplot overlay & topoplot sequences for all events
% ts_ezplot(epoch_data)    ==> multiplot overlay & topoplot sequences for all events
% ts_ezplot(timefreq_data) ==> multiplot overlay & topoplot sequences for all events
% ts_ezplot(stat_data)     ==> plot the mask
%
% EXAMPLES:
% Note: you can find these and more at neuroimaging.pbworks.com
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Parameters & data to use for all examples
% rootoutdir        = pwd;
% datafile          =
% '/space/mdkm1/9/kmdev/projects/Cog/MEG/nosss/P19j/P19j_SL_nosss.fif';
% trigchan          = 'STI101';                  % amplitude-coded or binary trigger channels
% chanlabels        = {'MEG0113','MEG0122','MEG0132'}; % labels of channels to process (set to [] for all); ex. {'MEG0113','MEG0122','MEG0132'}
% events            = [3 4];                     % valid event codes
% foi               = 'MEGlong';                 % Hz, frequencies of interest; ex. alpha band (8-12Hz)
% sf                = [];                        % Hz, spectral resolution
% conditionkey      = [];
% statfile          = [];
% rejectfile        = [];
% layoutfile        = '/home/jsherfey/svn/dev/onestream/vv_planar_2.lay'; % sensor locations (used for plotting)
%   % for MEG data:
%   %   /home/jsherfey/svn/dev/onestream/vv_planar.lay    <= grads, labels with whitespace (ex. MEG 0113)
%   %   /home/jsherfey/svn/dev/onestream/vv_planar_2.lay  <= grads, labels without whitespace (ex. MEG0113)
%   %   /home/jsherfey/svn/dev/onestream/vv_eeg.lay
%   %   /home/jsherfey/svn/dev/onestream/vv_all.lay
%   % for iEEG data: layoutfile = []
%  
% % Load & epoch data using amplitude code for TF analysis
% epoch_data = ts_process_data(datafile,'trigchan',trigchan,'chanlabels',chanlabels,'events',events);
% % Calc & return average power
% timefreq_data = ts_freqanalysis(epoch_data,'foi',foi,'sf',sf,'trials_flag',0,'complex_flag',0,'save_flag',0);
%   % NOTE: while timefreq_data contains average data, it also contains a list
%   % of the files with trial power that will be used to load all the TF data.
% % Preprocess epochs before calculating ERP or stats
% epoch_data = ts_preproc(epoch_data,'bandpass_flag',1,'bpfreq',[.1 40],'blc',1,'blcwindow',[-.2 0],'dsfact',2);
% % Average waveforms over trials
% avg_data = ts_trials2avg(epoch_data);
% % Limit stat analysis to select channels
% epoch_data = ts_data_selection(epoch_data,'chanlabels',chanlabels);
% % Calc nonparametric stats + temporal clustering
% stat_data = ts_statistics_wrapper(epoch_data,'events',[3 4]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Examples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1. Plotting epoch_data and avg_data
% %     - multiplot ERP waveforms
% %     - multiplot overlay of epoch waveforms
% %     - single ERP topoplot
% %     - series of ERP topoplots
% % 2. Plotting timefreq_data
% %     - multiplot TF spectral power
% %     - multiplot TF spectral power (z-score)
% %     - multiplot TF spectral power (frequency-normalized)
% %     - single TF spectral power topoplot
% %     - series of TF spectral power topoplots
% % 3. Customizing plots and other options (demonstrated with ERPs)
% %     - plotting ERPs + statistics
% %     - plotting ERPs + condition labels (legend)
% %     - plotting ERPs with different layouts
% %     - customizing figure properties
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% 1. Plotting epoch_data and avg_data
% % multiplot ERP waveforms
% ts_ezplot(avg_data);
% ts_ezplot(epoch_data);
% 
% % multiplot overlay of epoch waveforms
% ts_ezplot(epoch_data,'trials_flag',1);
% 
% % single ERP topoplot
% ts_ezplot(avg_data,'topoplot',1);
% 
% % series of ERP topoplots
% ts_ezplot(avg_data,'topoplot',1,'nrows',5,'ncols',5);
% 
% %% 2. Plotting timefreq_data
% % multiplot TF spectral power
% ts_ezplot(timefreq_data);
% 
% % multiplot TF spectral power (z-score)
% ts_ezplot(timefreq_data,'blc',1,'baselinetype','zscore','blcwindow',[-.2 0]);
% 
% % multiplot TF spectral power (frequency-normalized)
% ts_ezplot(timefreq_data,'freqcorr',1,'freqcorr_alpha',2);
% 
% % single TF spectral power topoplot (12-30Hz)
% ts_ezplot(timefreq_data,'foilim',[12 30],'topoplot',1);
% 
% % series of TF spectral power topoplots (12-30Hz)
% ts_ezplot(timefreq_data,'foilim',[12 30],'topoplot',1,'nrows',5,'ncols',5);
% 
% %% 3. Customizing plots and other options (demonstrated with ERPs)
% % plotting ERP multiplot + statistics
% ts_ezplot(avg_data,'stat_data',stat_data,'events',events);  % stat_data from ts_statistics or ts_statistics_wrapper
% 
% % If the stat_data structure is saved in statfile:
%   % ts_ezplot(avg_data,'statfile',statfile,'events',events);
% 
% % plotting series of ERP topoplots + statistics
% ts_ezplot(avg_data,'topoplot',1,'nrows',5,'ncols',5,'stat_data',stat_data);
% 
% % plotting ERPs + condition labels in the legend
% ts_ezplot(avg_data,'conditionkey',conditionkey);
% ts_ezplot(avg_data,'cond_labels',{'cond1','cond2'},'events',events);
% 
% % plotting ERPs with different layouts
% ts_ezplot(avg_data,'layout',layoutfile); % meaningful layout
% ts_ezplot(avg_data,'layout','ordered');  % grid of subplots
% 
% % saving figures
% ts_ezplot(avg_data,'save',1,'format','jpg' ,'close',0);  % jpeg
% ts_ezplot(avg_data,'save',1,'format','depsc','close',1);  % eps
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL NOTES
%
% condition labels:
%   cond_key        matlab structure with condition #s, event #s, & labels
%   cond_key_file   name of mat-file containing cond_key
%   conditionfile   name of csv file with cond_key info
%   cond_labels     cell array of condition labels
%     note: if cond_labels is given then the label order must correspond to
%     a sorted list of the unique event codes in the total event list
%
% baseline for zscore:
%   baseline_data   matlab structure
%   baselinefile    mat-file with baseline_data
%   
% stats for masking:
%   stat_data       matlab structure
%   statfile        mat-file with stat_data
%   mask            mask from stat_data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Jason Sherfey
% Modified: SRD 3/14/12 - to allow for F-tests and other >2 condition tests.
%                       - to change default graph colors to avoid yellow
            
 

scanflag = 0;  % do not search for data
if mod(nargin,2) && isstruct(varargin{1})
  % odd # of parameters AND 1st parm is data struct
  data = varargin{1};
  if nargin > 1
    varargin = varargin(2:end);
  else
    varargin = {};
  end
  data = ts_checkdata_header(data,varargin{:});
elseif mod(nargin,2)
  % odd # of parameters AND 1st parm is not data struct
  error('invalid specification of input parameters.');
else
  % even # of parameters (assume key/value pairs & search for data)
  scanflag = 1;
end
  
parms = mmil_args2parms(varargin,...
						{'datapath',pwd,[],...
             'rootoutdir',pwd,[],...
						 'datafile',[],[],...
						 'datatype',[],[],...
						 'datafield',[],[],...
						 'dataparam',[],[],...
             'zparam',[],[],...
             'object',[],[],...         
             'funspecs',[],[],...
             'badchanfile',[],[],...
             'rejectfile',[],[],...
             'baseline','no',[],...
             'baselinefile',[],[],...
             'baseline_data',[],[],...
             'baseline_start',-inf,[],...
             'baseline_end',0,[],...
             'statfile',[],[],...
             'stat_data',[],[],...
             'mask',[],[],...
             'maskstyle','transparency',[],...
             'baseline_mask',[],[],...
             'alpha',.05,[],...
             'conditionfile',[],[],...
             'conditionkey',[],[],...
             'cond_key_file',[],[],...          
             'condkey',[],[],...
             'cond_key',[],[],...
						 'cond_labels',[],[],...             
             'sensor_info',[],[],...
             'opt',[],[],...             
						 'cond',[],[],...             
             'condition',[],[],...
             'conditions',[],[],...
             'event',[],[],...
             'events',[],[],...
             'channel',[],[],...
             'channels',[],[],...
             'chantype','all', [],...
             'chanlabel',[],[],...
             'badchans',[],[],...
             'trial',[],[],...
             'trials',[],[],...
             'toilim',[],[],...
             'toi',[],[],...
             'foilim',[],[],...
             'foi',[],[],...       
						 'zlim','maxmin',[],...
						 'xlim','maxmin',[],...
						 'ylim','maxmin',[],... 
             'trials_flag',0,{1,0},...           
             'showerror_flag',0,{1,0},...
             'plotAvg',0,{1,0},...
             'plotwaveforms',0,{1,0},...
             'ezmultiplot',1,{1,0},...
						 'multiplot',1,{1,0},...
             'topoplot',0,{1,0},...
             'toprows',1,[],...
             'topcols',1,[],...
             'matchany',0,{1,0},...
             'matchall',0,{1,0},...
						 'singleplot',0,{1,0},...
						 'zscore',0,[],...             
             'baselinetype','zscore',[],...
             'overlay',1,{1,0},...             
             'close',0,{1,0},...             
						 'save',0,{1,0},...
             'overwrite',0,{1,0},...
						 'format','jpg',[],...
             'resolution',150,[],...
						 'outpath',[],[],... 
             'prefix',[],[],...
						 'layout',[],[],...
             'layoutfile',[],[],...
						 'nr',8,[],...
						 'nc',8,[],...
						 'showlabels','no',[],...             
             'showbadchans','yes',[],...
             'removebadchans',1,[],...
             'keepbadtrials',0,{1,0},...
						 'interactive','no',[],...
						 'stim','yes',[],...
						 'colorbar','no',[],...
						 'title','',[],...
             'comment','detailed',{'auto','xlim','detailed','details','none'},...
             'commentpos','leftbottom',[],...
						 'fontsize',8,[],...
             'axisfontsize',4,[],...
						 'skip_multievent_files',0,[],...     
             'newfig',1,[],...
             'get_zlims',0,[],...
             'preproc_flag',0,[],...
             'freqcorr',0,[],...
             'freqcorr_alpha',2,[],...
             'graphcolor','brgkmcywbrgkmcywbrgkmcywbrgkmcywbrgkmcywbrgkmcyw',[],...  % was: 'brgkywrgbkywrgbkywrgbkywbrgkywrgbkywrgbkywrgbkyw' SRD 3/14/12
             'statcolor','b',[],...
             'clip','yes',[],...
             'toffset',0,[],...
             'style',[],[],...
             'electrodes',[],[],...
             'box',[],[],...
             'highlight','off',[],...
             'hlmarker',[],[],...
             'hlcolor',[],[],...
             'hcolor',[],[],...
             'hlmarkersize',[],[],...
             'hllinewidth',[],[],...
             'markersize',20,[],...
             'linestyle','-',[],...
             'linewidth',.5,[],...
             'interpolation',[],[],...
             'interplimits',[],[],...
             'shading',[],[],...
             'contournum',[],[],...
             'colormap',[],[],...
             'avgovertime','no',[],...
             'axis','tight',[],...
             'axes',[],[],...
             'blcwindow',[-inf 0],[],...
             'autoscale',0,[],...
             'transparency',.1,[],...
             'footnote',[],[],...
             'figname',date,[],...
             'background','w',[],...
             'refchan',[],[],...
             'padfactor',.05,[],...
             'vline',[],[],...
             'hline',[],[],...    
             'vlinewidth',.2,[],...
             'hlinewidth',.2,[],...
             'vlinecolor','k',[],...
             'hlinecolor','k',[],...             
             'zerolinetimes',0,[],...
             'zerolinewidth',.25,[],...
             'zerolinecolor','k',[],...
             'guimode',0,[],...
             'ticklength',.1,[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...                
						},false);
parms = backcompatible(parms,varargin{:});

if scanflag
  scan_files(parms); 
  return;
end

if parms.singleplot || parms.topoplot, parms.multiplot = 0; parms.ezmultiplot = 0; end

% get metadata for the input data structure
[datatype,datafield,dataparam] = ts_object_info(data,varargin{:});
if isempty(data.(datafield))
  mmil_logstr(parms,'Structure does not contain data to plot.');
  if parms.verbose, fprintf('Structure does not contain data to plot.\n'); end
  return; 
end
if ischar(parms.dataparam) && ismember(parms.dataparam,dataparam)
  dataparam = parms.dataparam;
else
  dataparam = dataparam{1};
end
parms.datatype  = datatype;
parms.datafield = datafield;
parms.dataparam = dataparam;

% check time limits
if isempty(parms.toilim) || isequal(parms.toilim,'all')
  parms.toilim  = [data.(datafield)(1).time(1) data.(datafield)(1).time(end)];
end
if isequal(parms.xlim,'maxmin') || isempty(parms.xlim)
  parms.xlim    = parms.toilim;
end
if isequal(parms.zscore,1) || isequal(parms.zscore,'yes')
  if ~isempty(parms.blcwindow) && isnumeric(parms.blcwindow) && parms.blcwindow(1)<parms.toilim(1)
    parms.toilim(1) = max(parms.blcwindow(1),data.(datafield)(1).time(1));
  end
  if ~isempty(parms.blcwindow) && isnumeric(parms.blcwindow) && parms.blcwindow(2)>parms.toilim(2)
    parms.toilim(2) = min(parms.blcwindow(2),data.(datafield)(1).time(end));
  end
  if ~isempty(parms.xlim) && isnumeric(parms.xlim) && parms.xlim(1)<parms.toilim(1)
    parms.toilim(1) = max(parms.xlim(1),data.(datafield)(1).time(1));
  end
  if ~isempty(parms.xlim) && isnumeric(parms.xlim) && parms.xlim(2)>parms.toilim(2)
    parms.toilim(2) = min(parms.xlim(2),data.(datafield)(1).time(end));
  end  
end
if any(ismember('toilim',varargin(1:2:end)))
  ind = find(cellfun(@(x)(isequal(x,'toilim')),varargin));
  varargin{ind+1} = parms.toilim;
end

% select data to plot
args = mmil_parms2args(parms);
[data,jnk,err] = ts_data_selection(data,args{:});  
if err, return; end

if any(strcmp(datatype,{'timefreq_data','timefreq_avg','tfaverage','avg_timefreq','avg_power'}))
  powflag       = 1;
  parms.overlay = 0;
  dataparam = 'power';
  if isempty(findstr(parms.prefix,'power'))
    parms.prefix = sprintf('%s_power',parms.prefix); 
  end
  parms.newfig = 1;
else
  powflag = 0;
end
if isfield(data.(datafield),'cmplx') && ~isfield(data.(datafield),'power')
  for i = 1:length(data.(datafield))
    data.(datafield)(i).power = 2*abs(data.(datafield)(i).cmplx).^2;
  end
    data.(datafield) = rmfield(data.(datafield),'cmplx');
end
if any(strcmp(datatype,{'epoch_data'})) && parms.preproc_flag
  data = ts_preproc(data,varargin{:});
end
maskflag = ~(isempty(parms.stat_data) && isempty(parms.statfile) && isempty(parms.mask));
parms.maskflag = maskflag;
ieegflag = (any(strcmpi(parms.chantype,'ieeg')) || ismember('eeg',unique({data.sensor_info.typestring})) ...
           && isempty(intersect({'grad1','grad2','mag','meg','other'},{data.sensor_info.typestring})) );%MJE_20110630

% set up event and condition lists
if isempty(parms.events)  && ~isempty(parms.event), parms.events = parms.event; end
if ~isempty(parms.events) && ~iscell(parms.events), parms.events = {parms.events}; end
if iscell(parms.events),     evlist = parms.events;     parms.events     = unique([parms.events{:}]);     end
if iscell(parms.condition),  cdlist = parms.condition;  parms.condition  = unique([parms.condition{:}]);  end 
if iscell(parms.conditions), cdlist = parms.conditions; parms.conditions = unique([parms.conditions{:}]); end 
if ~exist('evlist','var'), evlist = {[data.(datafield).event_code]}; end
if ~exist('cdlist','var') || (length(evlist) ~= length(cdlist))
  for k = 1:length(evlist)
    [jnk loc] = ismember(evlist{k},[data.(datafield).event_code]);
    cdlist{k} = loc(loc~=0);
  end
end
if length(cdlist)==1 && ~parms.overlay
  cdlist = num2cell([cdlist{:}]); 
  evlist = num2cell([evlist{:}]); 
end

% remove index specification already used for data selection
parms.condition  = [];
parms.conditions = [];
parms.channels   = [];

% convert data to averages if trials are present
if ~parms.trials_flag && (strcmp(datatype,'epoch_data') || ...
   (powflag && ndims(data.(datafield)(1).(dataparam)) > 3))
  mmil_logstr(parms,'averaging data over epochs'); 
  if powflag, datafield = 'timefreq'; else datafield = 'averages'; end
  parms.datafield = datafield;
  args = mmil_parms2args(parms);
  data = ts_trials2avg(data,args{:});
  parms.datafield = datafield;
  clear args
end
% now: data is either [ERP] = chan x time, or [power] = chan x time x freq

if parms.plotAvg && ~powflag
  if isempty(parms.datafile)
    error('''plotAvg'' requires the full filename of an epoch_data or avg_data structure.\nNote: pass the data structure using ''multiplot'' or the filename using ''plotAvg.''');
  end
  cfg  = ezplot2PlotAvg(parms,ieegflag,evlist);
  args = mmil_parms2args(cfg);  
  ts_iEEG_PlotAvg(args{:});
  return;
elseif parms.plotwaveforms && ~powflag
  args = mmil_parms2args(parms);
%   ts_iEEG_Plot_Waveforms(data,args{:});
  ts_iEEG_Plot_Waveforms(data,'baseline',parms.baseline,'ylim',parms.ylim);
  return;
end

% split channels by type
if ieegflag
  chantype = {'ieeg'};
elseif strcmp(parms.chantype,'meg')
  chantype = intersect({'grad1','grad2','mag'},unique({data.sensor_info.typestring}));
elseif strcmp(parms.chantype,'all')
  chantype = unique({data.sensor_info.typestring});
elseif strcmp(parms.chantype,'grad')
  if parms.topoplot% || powflag
    chantype = {'grad1','grad2'};
  else
    chantype = {'grad'};
  end
else
  chantype = {parms.chantype};
end
conds_orig  = parms.cond_labels;
parms_orig  = parms;
% sensor_orig = data.sensor_info;
zlims    = [0 0];
clearlay = 0;
%% main loop
% loop over channel types
for g = 1:length(chantype)
  
  % loop over condition vectors in condition list
  for k = 1:length(cdlist)
    % set parameters to plot channels in this group for events in this condition vector
    parms.chantype = chantype{g};              % this group
    conds          = cdlist{k};                % this condition vector
    if ~isempty(conds_orig) && (length(conds_orig) == length(unique([cdlist{:}])))
      parms.cond_labels = conds_orig(sort(conds)); % = conds_orig(conds);
    end
    parms.events   = evlist{k};
    parms          = load_cond_labels(parms);  % get 'cond_labels' => cell array of condition labels
%     title          = parms.title;          
    condstr        = sprintf('\n');
    comment        = '';%sprintf('\n%s',parms.comment);
    % todo: add chantype & event_codes to the comment if cond_labels is empty
    
    % data for this condition vector
    dat = rmfield(data,datafield);
    dat.(datafield) = data.(datafield)(conds);
    dat = ts_data_selection(dat,'chantype',parms.chantype,'removebadchans',parms.removebadchans,'keepbadtrials',parms.keepbadtrials); % added 3/26/2010
    sensor_orig = dat.sensor_info;
    % zscore 
    if parms.zscore
      dat     = ts_zscore(dat,'toilim',parms.toilim,'blcwindow',parms.blcwindow,'baselinefile',parms.baselinefile,...
                'baseline_data',parms.baseline_data,'baselinetype',parms.baselinetype,'zparam',parms.zparam,'removebadchans',parms.removebadchans,'verbose',parms.verbose);
      comment = sprintf(' \n%s%s',strrep(parms.baselinetype,'absolute','blc'),comment);
    end

    % add typestring prefix if label contains digits only
    if any(cellfun('isempty',regexp({dat.sensor_info.label},'\D+')))
      ix = find(cellfun('isempty',regexp({dat.sensor_info.label},'\D+')));
      for j = 1:length(ix)
        dat.sensor_info(ix(j)).label = [dat.sensor_info(ix(j)).typestring dat.sensor_info(ix(j)).label];
      end
    end
    
    % iEEG layout
    if ieegflag
      if ~isempty(parms.outpath)
        outpath = parms.outpath;
        if ~strcmp(outpath(1),'/'), outpath = fullfile(parms.rootoutdir,outpath); end
      else
        outpath = parms.datapath;
      end      
      if isempty(parms.layout)
        try
          parms.layout = write_layout_ieeg(outpath,dat.sensor_info,parms.nr,parms.nc,parms.padfactor);
          npages   = length(parms.layout);
          clearlay = 1;
        catch
          npages = 1;
          parms.layout = 'ordered';
        end  
      elseif ischar(parms.layout)
        % do nothing
      elseif ~isstruct(parms.layout) && ~iscell(parms.layout)
        try
          parms.layout = write_layout_ieeg(outpath,dat.sensor_info,parms.nr,parms.nc,parms.padfactor,parms.layout);
          npages   = length(parms.layout);
          clearlay = 1;
        catch
          npages = 1;
          parms.layout = 'ordered';
        end  
      elseif ~iscell(parms.layout)
        npages   = 1;
        clearlay = 0;
      end
      clear outpath
    else
      label = {dat.sensor_info.label};
      ix = find(cellfun('isempty',regexp(label,'MEG\s+\d+')));
      for j = 1:length(ix)
%         dat.sensor_info(ix(j)).label = strrep(dat.sensor_info(ix(j)).label,'MEG','MEG ');
      end
    end
    if ~isempty(parms.layout) && iscell(parms.layout)
      npages = length(parms.layout);
    else
      npages = 1;
    end
    if ieegflag && ~ismember('ieeg',unique({dat.sensor_info.typestring}))
      parms.chantype = 'eeg';
    end    
%%    
    for c = 1:length(conds)
      if parms.freqcorr ~= 0 && powflag
        mmil_logstr(parms,'correcting for 1/f^%g scaling of power spectrum',parms.freqcorr_alpha);
        freqcorr = permute(dat.(datafield)(c).frequencies.^(parms.freqcorr_alpha),[1 3 2]); % 1 x 1 x freqs
        freqcorr = repmat(freqcorr,[dat.num_sensors length(dat.(datafield)(c).time) 1]);
        dat.(datafield)(c).power = dat.(datafield)(c).power .* freqcorr;  
        clear freqcorr
      end
    end
    if parms.ezmultiplot
      % prepare TS data for plotting
%       for c = 1:length(conds)
%         dat.(datafield)(c).time = dat.(datafield)(c).time + parms.toffset;
%         if c == 1 && isnumeric(parms.xlim) && ~isempty(parms.xlim)
%           parms.xlim = parms.xlim + parms.toffset;
%         end
%       end
      layout  = parms.layout;
      parms.evcode  = [dat.(datafield).event_code];
      if maskflag
        sensor_temp     = dat.sensor_info;
        dat.sensor_info = sensor_orig;
        [dat,parms] = apply_mask(dat,parms,[dat.(datafield).event_code]);
        [jnk cix] = intersect({sensor_orig.label},{dat.sensor_info.label});
        dat.sensor_info = sensor_temp(sort(cix));
        clear sensor_temp
        comment = sprintf('%s\nSig: p < %g\n',comment,parms.alpha);
        mask    = dat.(datafield)(c).mask;        
%         mask    = dat.(datafield)(conds(1)).mask;
      else
        mask = [];
      end
      for c = 1:length(conds)
        dat.(datafield)(c).time = dat.(datafield)(c).time + parms.toffset;
        if c == 1 && isnumeric(parms.xlim) && ~isempty(parms.xlim)
          parms.xlim = parms.xlim + parms.toffset;
        end
      end      
      if ~strcmpi(parms.comment,'none')
        parms.comment = sprintf('%s\n%s\n%s',condstr,parms.chantype,comment);      
      end
        
    else
      % convert TS to FT & prepare for plotting
      for c = 1:length(conds)
        % whitening (correction for 1/f^2 scaling of power spectrum)
%         if parms.freqcorr ~= 0 && powflag
%           freqcorr = permute(dat.timefreq(c).frequencies.^(parms.freqcorr_alpha),[1 3 2]); % 1 x 1 x freqs
%           freqcorr = repmat(freqcorr,[dat.num_sensors length(dat.timefreq(c).time) 1]);
%           dat.timefreq(c).power = dat.timefreq(c).power .* freqcorr;  
%           clear freqcorr
%         end

        evcode(c) = dat.(datafield)(c).event_code;
        % apply statistical mask
        if maskflag
          sensor_temp     = dat.sensor_info;
          dat.sensor_info = sensor_orig;          
          if parms.topoplot || parms.overlay
            [dat,parms] = apply_mask(dat,parms,[dat.(datafield).event_code]);
          else
            [dat,parms] = apply_mask(dat,parms,evcode(c));
          end
          try comment = sprintf('%s\nSig: p < %g\n',comment,parms.alpha); end
          [jnk jnk2 cix] = intersect({sensor_orig.label},{dat.sensor_info.label});
          dat.sensor_info = sensor_temp(sort(cix));
          clear sensor_temp          
        end

        % convert to fieldtrip
        ftdata{c} = ts_data2fieldtrip(dat,'condition',c,'chantype',parms.chantype);%,'channels',parms.channels);    
        evcode(c) = dat.(datafield)(c).event_code;

        % add condition labels to title string
        if ~isempty(parms.cond_labels) && length(parms.cond_labels) == length([cdlist{:}])
          condstr = [condstr sprintf('%s\nNum Trials: %g\n',parms.cond_labels{parms.events==evcode(c)},dat.(datafield)(c).num_trials)];
        end
      end
      parms.evcode = evcode;
      if isempty(condstr), condstr = sprintf('events %s\n',num2str(evcode)); end      
    end
%%    
    % loop over pages    
    for p = 1:npages
      parms.page = p;
      
      % plot TS data
      if parms.ezmultiplot
        try parms.layout = layout{p}; catch parms.layout = layout; end
        args = mmil_parms2args(rmfield(parms,{'stat_data','statfile','mask','toilim','channels','rejectfile','badchans','badchanfile'}));
        [jnk,err] = ts_ezmultiplot(dat,'mask',mask,args{:});      
        if ~err
          if parms.save,  parms.toi=dat.(datafield)(1).time; save_plots(parms,'multiplot'); end
        else
          mmil_logstr(parms,'an ezmultiplot error occurred while processing page %g',p);
        end
        if parms.close, close; end        
        continue;
      end
      
      % plot FT data
      % set up fieldtrip configuration
      cfg = [];
      cfg = ftconfig(parms);
      if iscell(parms.layout)% && length(parms.layout) > 1
        cfg.layout = parms.layout{p};
      elseif ischar(parms.layout) || isstruct(parms.layout)
        cfg.layout = parms.layout;
      end
      if strcmpi(parms.comment,'detailed') || strcmpi(parms.comment,'details')
%         if ieegflag
%           cfg.comment = sprintf('%s\n%s\n%s',parms.title,condstr,'intracranial EEG',comment);
%         else
          cfg.comment = sprintf('%s\n%s\n%s',parms.title,condstr,parms.chantype,comment);
%         end
        try cfg = rmfield(cfg,'title'); end
      elseif strcmpi(parms.comment,'xlim') || strcmpi(parms.comment,'auto')
        cfg.comment = parms.comment;
      end
      
      if parms.newfig && ~parms.get_zlims
        screensize = get(0,'ScreenSize');
        figure('Name',parms.figname,'NumberTitle','on','Color',parms.background,'Position',[1 1 .9*screensize(3) .9*screensize(4)]);
%         figure; 
      end
      % call plotting function
      if parms.topoplot
        try 
          cfg.get_zlims = parms.get_zlims;
          z = plot_topoplot(cfg,ftdata,parms,powflag,ieegflag);
          try
            if abs(z(1)) > abs(zlims(1))
              zlims = z;
            end
          end
        catch
          mmil_logstr(parms,'topoplot failed!');
          close; %continue;
        end
      end
      if parms.multiplot
          plot_multiplot(cfg,ftdata,parms,powflag);
      end          
      if parms.singleplot
        try
          plot_singleplot(cfg,ftdata,parms,powflag);
        catch
          mmil_logstr(parms,'singleplot failed!');          
          close;
        end
      end        
      if parms.close
        close;% all force; 
      end      
    end
    clear ftdata evcode cfg
    parms = parms_orig;  
  end   % end loop over condition vectors in one condition list
  if clearlay
    parms.layout = [];
  end
%   parms = parms_orig;
end     % end loop over channel types

%% subfunctions

function zlims = plot_topoplot(cfg,ftdata,parms,powflag,ieegflag)
zlims = [];
% add ability to handle iEEG grid (hard code grid based on labels Gxx)
%     topoplotER_iEEG(cfg,ftdata);
%     ts_MEGtopoplot();     

if ieegflag
  topoplotER_iEEG(cfg,ftdata{1});
  if parms.save, parms.toi=ftdata{1}.time; save_plots(parms,'topoplot'); end
  return;
end

if parms.toffset ~= 0
  for i = 1:length(ftdata)
    ftdata{i}.time = ftdata{i}.time + parms.toffset;
  end
  if isnumeric(parms.xlim) && ~isempty(parms.xlim)
    parms.xlim = parms.xlim + parms.toffset;
  end
  
end

nr = parms.toprows;
nc = parms.topcols;
N  = nr*nc;
comment = cfg.comment;
cfg = rmfield(cfg,'colorbar');

if isfield(cfg,'electrodes')
  if ~isempty(cfg.electrodes)
    cfg = rmfield(cfg,'showlabels');
  else
    cfg = rmfield(cfg,'electrodes'); 
  end
end
% toi
if N > 1
  t   = ftdata{1}.time;
  if isnumeric(parms.xlim) && ~isempty(parms.xlim)
    t0 = parms.xlim(1);
    tf = parms.xlim(end);
  else
    t0 = t(1);
    tf = t(end);
  end
  if strcmp(parms.avgovertime,'yes')
    dt  = (tf-t0) / N;
    toi = t0 + dt*[0:N];    
  else
    dt  = (tf-t0) / (N-1);
    toi = t0 + dt*[0:N-1];
  end
else
  cfg.xlim = [ftdata{1}.time(1) ftdata{1}.time(end)];
end

% foi
if powflag % [ftdata.powspctrm] = chan x freq x time
  cfg.ylim = [ftdata{1}.freq(1) ftdata{1}.freq(end)];
  zparam = 'powspctrm';
else
  zparam = 'avg';
end

% zlim
if strcmp(parms.zlim,'maxmin') || strcmp(parms.zlim,'absmax') || strcmp(parms.zlim,'symmetric')
  zmin = min(min(min(ftdata{1}.(zparam))));
  zmax = max(max(max(ftdata{1}.(zparam))));
  if strcmp(parms.zlim,'symmetric') || strcmp(parms.zlim,'absmax') 
    zmax = max(abs(zmin),abs(zmax));
    zmin = -zmax;
  end
  cfg.zlim = [zmin zmax];
  if cfg.get_zlims
      zlims = cfg.zlim;
      mmil_logstr(parms,'%s zlimits = [%g %g]',parms.zlim,zlims(1),zlims(2));
    return;
  end
end

% difference condition
if length(ftdata) == 2
  ftdata{1}.(zparam) = ftdata{1}.(zparam) - ftdata{2}.(zparam);
end

% plot series of topoplots
for k=1:N
  if N > 1
    subplot(nr,nc,k)
    if strcmp(parms.avgovertime,'yes')
      cfg.xlim = [toi(k) toi(k+1)];
    else
      cfg.xlim = [toi(k) toi(k)];
    end
  end
%   cfg.comment = sprintf('%0s\ntime=[%.3g %.3g]',comment,cfg.xlim(1),cfg.xlim(2));
  if k == N && strcmp(cfg.comment,'xlim'), cfg.comment = 'auto'; end
  if parms.autoscale~=0
    t = ftdata{1}.time;
    if powflag
      zmin = min(min(min(ftdata{1}.(zparam)(:,:,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))))));
      zmax = max(max(max(ftdata{1}.(zparam)(:,:,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))))));      
    else
      zmin = min(min(min(ftdata{1}.(zparam)(:,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))))));
      zmax = max(max(max(ftdata{1}.(zparam)(:,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))))));
    end
    cfg.zlim = [zmin zmax];    
  end
  if powflag
%     cfg.comment = sprintf('%0s\nfreq=[%.3g %.3g]',cfg.comment,cfg.ylim(1),cfg.ylim(2));
%     ts_topoplotTFR(cfg,ftdata{:});
    if isfield(parms,'maskparameter') 
      t   = ftdata{1}.time;
      if ~(ischar(parms.mask) && ~strcmp(parms.mask,'off'))
        cfg.mask = squeeze(ftdata{1}.(cfg.maskparameter)(:,1,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))));
        if N==1 || strcmp(parms.avgovertime,'yes')
          cfg.mask = any(any(cfg.mask,3),2);
        end        
      end
      if ischar(parms.highlight) && strcmp(parms.highlight,'on')
        mdata = squeeze(ftdata{1}.(cfg.maskparameter)(:,1,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))));
        cfg.highlight = find(any(any(mdata,3),2));
      end
    end
    try 
      if strcmp(parms.mask,'off'), cfg = rmfield(cfg,'mask'); end; 
    end
    topoplotTFR(cfg,ftdata{1});
  else
    if isfield(parms,'maskparameter') 
      t   = ftdata{1}.time;
      if ~(ischar(parms.mask) && strcmp(parms.mask,'off'))
        cfg.mask = squeeze(ftdata{1}.(cfg.maskparameter)(:,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))));
        if N==1 || strcmp(parms.avgovertime,'yes')
          cfg.mask = any(cfg.mask,2);
        end
      end
      if ischar(parms.highlight) && strcmp(parms.highlight,'on')
        mdata = squeeze(ftdata{1}.(cfg.maskparameter)(:,nearest(t,cfg.xlim(1)):nearest(t,cfg.xlim(2))));
        cfg.highlight = find(any(mdata,2));
      end      
    end
    try 
      if strcmp(parms.mask,'off'), cfg = rmfield(cfg,'mask'); end; 
    end
    topoplotER(cfg,ftdata{1});
  end
%   if k==1 && strcmpi(parms.colorbar,'yes'), colorbar; end;
  if isfield(cfg,'mask'), cfg = rmfield(cfg,'mask'); end
end
% add title
if ~isempty(parms.title)
  parms.title = strrep(parms.title,'_','\_');
  annotation('textbox','Position',[0 .9 1 .1],'VerticalAlignment','middle','HorizontalAlignment','center',...
             'Color','k','FontSize',10,'FontWeight','bold','FitHeightToText','on','LineStyle','none','String',parms.title);
end
% add colorbar
if strcmp(parms.colorbar,'yes')
  if N > 1, subplot(nr,nc,N); end
  colorbar('location','EastOutside'); 
end
% save figure
if parms.save
  if ischar(parms.mask) && strcmp(parms.mask,'off'), parms.maskflag = 0; end
  parms.toi=ftdata{1}.time; 
  save_plots(parms,'topoplot'); 
end
zlims = cfg.zlim;

function plot_multiplot(cfg,ftdata,parms,powflag)
  if powflag
    % [ftdata.powspctrm] = chan x freq x time
    multiplotTFR(cfg,ftdata{:});
  else
    % [ftdata.avg] = chan x time
    multiplotER(cfg,ftdata{:});
  end
  if parms.save, parms.toi=ftdata{1}.time; save_plots(parms,'multiplot'); end
     
function plot_singleplot(cfg,ftdata,parms,powflag)
  if powflag
    % [ftdata.powspctrm] = chan x freq x time
    singleplotTFR(cfg,ftdata{:});
    xlabel('time (s)'); ylabel('frequency (Hz)')
  else
    % [ftdata.avg] = chan x time
    singleplotER(cfg,ftdata{:});
  end
  if parms.save, parms.toi=ftdata{1}.time; save_plots(parms,'singleplot'); end  
  
function [data,parms] = apply_mask(data,parms,evcode) 
baseline_mask = []; 
% load mask
if ~isempty(parms.mask) && isnumeric(parms.mask)
  mask = parms.mask;
  c = find([data.(parms.datafield).event_code]==evcode);
  data.(parms.datafield)(c).(parms.dataparam) = data.(parms.datafield)(c).(parms.dataparam) .* mask;
  return;
end
if ~isempty(parms.statfile) 
  mmil_logstr(parms,'loading statistics file: %s',parms.statfile);
  load(parms.statfile);
elseif ~isempty(parms.stat_data) 
  stat_data = parms.stat_data; 
else
  return;
end

% 3/14/12 SRD & ML, changed if from '== 2' to '>= 2' so that F-Tests etc. will work.
if length(evcode) >= 2 
  if all(ismember(evcode,[stat_data.stats.event_code]))
    mask = stat_data.stats(1).prob <= parms.alpha;
%     baseline_mask{1} = stat_data.stats([stat_data.stats.event_code]==evcode(1)).prob <= parms.alpha;
%     baseline_mask{2} = stat_data.stats([stat_data.stats.event_code]==evcode(2)).prob <= parms.alpha;
    for i=1:length(evcode) baseline_mask{i} = stat_data.stats([stat_data.stats.event_code]==evcode(i)).prob <= parms.alpha; end  % replaces above, SRD & ML 3/14/12
  else 
    return;
  end
elseif length(evcode) == 1 && ismember(evcode,[stat_data.stats.event_code]) 
  mask = stat_data.stats([stat_data.stats.event_code]==evcode).prob <= parms.alpha;
elseif length(evcode) == 1 
  mask = stat_data.stats(1).prob <= parms.alpha;
else
  return;
end
% if issubfield(stat_data.stats(1),'cfg.alpha')
%   parms.alpha = stat_data.stats(1).cfg.alpha;
% end
parms.maskparameter = 'mask';

% restrict to common channels
[sel1 sel2] = match_str({data.sensor_info.label},{stat_data.sensor_info.label});
[jnk rm1]  = intersect(sel1,find([data.sensor_info.badchan]));
[jnk rm2]  = intersect(sel2,find([stat_data.sensor_info.badchan]));
sel1([rm1 rm2])  = [];
sel2([rm1 rm2])  = [];
data.sensor_info = data.sensor_info(sel1);
data.num_sensors = length(sel1);
stat_data.sensor_info = stat_data.sensor_info(sel2);
stat_data.num_sensors = length(sel2);
% [sel1 sel2] = match_str({data.sensor_info.label},{stat_data.sensor_info.label});

% match dimensions b/w data & mask for all relevant conditions
c = find([data.(parms.datafield).event_code]==evcode);
for i = 1:length(c)
  cc = c(i);
  % truncate times
  xmin = max(min(data.(parms.datafield)(cc).time),min(stat_data.stats(1).time));
  xmax = min(max(data.(parms.datafield)(cc).time),max(stat_data.stats(1).time));
%   xix1 = data.(parms.datafield)(cc).time >= xmin & data.(parms.datafield)(cc).time <= xmax;
%   xix2 = stat_data.stats(1).time         >= xmin & stat_data.stats(1).time        <= xmax;
  xix1 = nearest(data.(parms.datafield)(cc).time,xmin):nearest(data.(parms.datafield)(cc).time,xmax);
  xix2 = nearest(stat_data.stats(1).time,xmin):nearest(stat_data.stats(1).time,xmax);
  if length(xix1) < length(xix2), 
    xix2 = xix2(1:end-(length(xix2)-length(xix1)));
  end
  data.(parms.datafield)(cc).(parms.dataparam) = data.(parms.datafield)(cc).(parms.dataparam)(sel1,xix1,:,:);
  data.(parms.datafield)(cc).time              = data.(parms.datafield)(cc).time(xix1);
  data.(parms.datafield)(cc).mask              = mask(sel2,xix2,:,:);
  if ~isempty(baseline_mask) && iscell(baseline_mask) && length(c)==length(baseline_mask)
    baseline_mask{i} = baseline_mask{i}(sel2,xix2,:,:);
  end    
%   if ~isempty(baseline_mask) && iscell(baseline_mask)
%     for j = 1:length(baseline_mask)
%       baseline_mask{j} = baseline_mask{j}(sel2,xix2,:,:);
%     end
%   end    
  if isfield(data.(parms.datafield),'stdev')
    data.(parms.datafield)(cc).stdev = data.(parms.datafield)(cc).stdev(sel1,xix1,:,:);
  end
end
try x1 = parms.xlim(1); catch x1 = -inf; end
try x2 = parms.xlim(2); catch x2 = inf;  end
if strcmp(parms.xlim,'maxmin')
    parms.xlim = [max([cellfun(@min,{data.(parms.datafield)(c).time}) -inf]) ...
                  min([cellfun(@max,{data.(parms.datafield)(c).time}),inf])];
    
else
    parms.xlim = [max([cellfun(@min,{data.(parms.datafield)(c).time}) x1]) ...
                  min([cellfun(@max,{data.(parms.datafield)(c).time}),x2])];    
end

% account for frequencies
if isfield(data.(parms.datafield),'frequencies')
  ymin = max([cellfun(@min,{data.(parms.datafield).frequencies}),...
              cellfun(@min,{stat_data.stats.frequencies})]);
  ymax = min([cellfun(@max,{data.(parms.datafield).frequencies}),...
              cellfun(@max,{stat_data.stats.frequencies})]);
  data = ts_data_selection(data,'foilim',[ymin ymax]);
  yidc = nearest(stat_data.stats(1).frequencies,ymin):nearest(stat_data.stats(1).frequencies,ymax);
% 	yidc = stat_data.stats(1).frequencies >= ymin & stat_data.stats(1).frequencies <= ymax;
%   mask = mask(:,:,yidc,:);
  mask = mask(sel2,xix2,yidc,:);
  nf = size(mask,3);
  if parms.matchany
    k = any(mask,3);
    mask(:,:,1:nf) = repmat(k,[1 1 nf]);
  elseif parms.matchall
    k = all(mask,3);
    mask(:,:,1:nf) = repmat(k,[1 1 nf]);    
  end
else
%   data = ts_data_selection(data,'toilim',toilim);
end

if ndims(mask) == 3 && strcmp(parms.maskstyle,'contour')
  % TF mask (convert to contour)
  nt   = size(mask,2);
  nf   = size(mask,3);
  lt0 = (mask(:,1:nt-2,2:nf-1)==0);
  rt0 = (mask(:,3:nt,2:nf-1)==0);
  tp0 = (mask(:,2:nt-1,3:nf)==0);
  bt0 = (mask(:,2:nt-1,1:nf-2)==0);
  MASK = mask(:,2:nt-1,2:nf-1).*((lt0 + rt0 + tp0 + bt0) > 0);
  MASK = cat(2,mask(:,1,2:nf-1),MASK,mask(:,nt,2:nf-1)); % chan x time x freq
  MASK = cat(3,mask(:,:,1),MASK,mask(:,:,nf)); % chan x time x freq
  mask = MASK;
  clear MASK;
end

if ~isfield(data.(parms.datafield),'mask')
  for i = 1:length(c)
    data.(parms.datafield)(c(i)).mask = mask;
  end
end
if isempty(parms.baseline_mask)
  parms.baseline_mask = baseline_mask;
end
% apply mask
% data.(parms.datafield)(c).(parms.dataparam) = data.(parms.datafield)(c).(parms.dataparam) .* mask;

function parms = load_cond_labels(parms)
% % load a cell array of condition labels
if ~isempty(parms.conditionkey) && isempty(parms.conditionfile)
  parms.conditionfile = parms.conditionkey;
end
if ~isempty(parms.condkey) && isempty(parms.cond_key)
  parms.cond_key_file = parms.condkey;
end
if ~isempty(parms.conditionfile) && exist(parms.conditionfile,'file')
    ts_makecondkey(parms.conditionfile);
    [pathstr,name] = fileparts(parms.conditionfile);
    load(fullfile(pathstr,'cond_key.mat'));
    idx = cellfun(@(x)find([cond_key.event_code]==x),num2cell(parms.events));
    parms.cond_labels = [cond_key.name];
    parms.cond_labels = parms.cond_labels(idx);
elseif ~isempty(parms.cond_key)
    idx = cellfun(@(x)find([parms.cond_key.event_code]==x),num2cell(parms.events));
    parms.cond_labels = [parms.cond_key.name];  
    parms.cond_labels = parms.cond_labels(idx);
elseif ~isempty(parms.cond_key_file) && exist(parms.cond_key_file,'file')
    load(parms.cond_key_file);
    idx = cellfun(@(x)find([cond_key.event_code]==x),num2cell(parms.events));
    parms.cond_labels = [cond_key.name];
    parms.cond_labels = parms.cond_labels(idx);
end 
if isfield(parms,'cond_labels') && ~isempty(parms.cond_labels)
  for i = 1:length(parms.cond_labels)
    parms.cond_labels{i} = strrep(parms.cond_labels{i},'_','\_');
  end
end
return;

function cfg = ftconfig(parms)
pn = {'showlabels','xlim','ylim','zlim','interactive',...
      'colorbar','stim','fontsize','layout','maskparameter','mask',...
      'style','electrodes','title','box','highlight','hlmarker','hcolor','hlcolor','hlmarkersize','hllinewidth',...
      'interpolation','interplimits','shading','contournum','colormap',...
      'commentpos'};
for p = 1:length(pn)
  try
    if isfield(parms,pn{p}) && ~isempty(parms.(pn{p}))
      cfg.(pn{p}) = parms.(pn{p});
    end
  end
end

function scan_files(parms)
% get data from datafiles
% find matfiles
if isempty(parms.datafile) || (~iscellstr(parms.datafile) && ~exist(parms.datafile,'file'))
  datafile = what(parms.datapath);
  datafile = datafile.mat;
else
  datafile = parms.datafile;
end
if ~iscell(datafile), datafile = {datafile}; end
% load data
for f = 1:length(datafile)
  parms.datafile = datafile{f};
  mmil_logstr(parms,'Scanning file %i of %i',f,length(datafile));
  if parms.verbose, fprintf('Scanning file %i of %i\n',f,length(datafile)); end
  data  = load(parms.datafile);
  itype = fieldnames(data);
  itype = itype{1};
  data  = data.(itype);
  % check that data is a valid timesurfer structure
  try
    args = mmil_parms2args(parms);
    [datatype,datafield,dataparam] = ts_object_info(data,args{:});
    if ~strcmp(datatype,itype), continue; end
  catch
    continue;
  end
  args = mmil_parms2args(parms);
  ts_ezplot(data,args{:});
  clear data;
end
return;

function save_plots(parms,plottype)
if isempty(parms.outpath), parms.outpath = fullfile(parms.rootoutdir,'images'); end
if ~strcmp(parms.outpath(1),'/'), parms.outpath = fullfile(parms.rootoutdir,parms.outpath); end
if ~exist(parms.outpath,'file'), mkdir(parms.outpath); end
if ischar(parms.resolution)
  res = sprintf('-r%s',parms.resolution); % dpi
else
  res = sprintf('-r%g',parms.resolution); % dpi
end
if iscell(parms.format) && strcmp(parms.format{1}(1),'-')
  printcmd     = parms.format(1:end-1);
  tmp          = parms.format;
  parms.format = tmp{end};
elseif iscell(parms.format)
  printcmd     = {['-' parms.format{1}]}; 
  tmp          = parms.format; 
  parms.format = tmp{2};
elseif any(strcmp(parms.format,'jpg')),  printcmd = {'-djpeg',res};
elseif any(strcmp(parms.format,'eps')),  printcmd = {'-depsc','-tiff',res};  
elseif any(strcmp(parms.format,'eps2')), printcmd = {'-depsc2','-tiff',res}; parms.format = 'eps'; 
elseif any(strcmp(parms.format,'depsc')),printcmd = {'-depsc',res}; parms.format = 'eps'; 
elseif any(strcmp(parms.format,'tif')),  printcmd = {'-dtiff'}; parms.format = 'tiff';
else
  printcmd  = {['-' parms.format]};
end
if     ~isempty(parms.prefix),  prefix = parms.prefix;
elseif ~isempty(parms.datafile), [jnk prefix] = fileparts(parms.datafile{1});
else   prefix = parms.datafield;
end

evstring = '';
if isfield(parms,'cond_labels') && iscell(parms.cond_labels)
  try
    for c = 1:length(parms.evcode)
      cstr = parms.cond_labels{parms.events==parms.evcode(c)};
      cstr = strrep(cstr,' ','_');
      cstr = strrep(cstr,'\_','_');
      cstr = strrep(cstr,'\\','');
      evstring = sprintf('%sevent%g_%s_',evstring,parms.evcode(c),cstr);
    end
  end
end
if isempty(evstring), evstring = sprintf('events%s_',strrep(num2str(parms.evcode),'  ','_')); end
if     ~isempty(parms.foi),    foistring = sprintf('_foi%g-%g',parms.foi(1),parms.foi(end));
elseif ~isempty(parms.foilim), foistring = sprintf('_foi%g-%g',parms.foilim(1),parms.foilim(end));
else                           foistring = '';
end
if parms.maskflag, maskstring = '_masked';
else               maskstring = '';
end
if isnumeric(parms.xlim)
  if parms.xlim(1)   > parms.toi(1),   parms.toi(1)   = parms.xlim(1);   end
  if parms.xlim(end) < parms.toi(end), parms.toi(end) = parms.xlim(end); end
end
if exist(parms.rejectfile,'file')
  prefix  = [prefix '_rej'];
end
if isnumeric(parms.zscore) && parms.zscore~=0
  prefix  = [prefix '_' parms.baselinetype];
end
if isnumeric(parms.freqcorr) && parms.freqcorr ~=0
  prefix  = [prefix '_freq-corrected'];
end
if parms.autoscale
  prefix  = [prefix '_autoscale'];
elseif isnumeric(parms.zlim)
  prefix = sprintf('%s_zlim%g-%g',prefix,parms.zlim(1),parms.zlim(end));
elseif ischar(parms.zlim)
  prefix = [prefix '_zlim-' parms.zlim];
end
outfile = sprintf('%s/%s_%s_%s_%stoi%.2g-%.2g%s%s_page%i.%s',parms.outpath,prefix,plottype,parms.chantype,evstring,parms.toi(1),parms.toi(end),foistring,maskstring,parms.page,parms.format);
[tmppath tmpname] = fileparts(outfile);
if length(tmpname) > 150
  mmil_logstr(parms,'Warning: filename exceeds maximum length. Shortening to 150 characters.');
%   fprintf('Warning: filename exceeds maximum length. Shortening to 150 characters.\n');
  outfile = fullfile(tmppath,[tmpname(1:148) '_' num2str(parms.page) '.' parms.format]);
end
if ~isempty(parms.footnote)
  if strcmp(parms.footnote,'filename')
%     annotation('textbox','Position',[0 0 1 .05],'VerticalAlignment','bottom','HorizontalAlignment','left',...
%              'Color','k','FontSize',4,'FitHeightToText','on','LineStyle','none','String',strrep(outfile,'_','\_'));
    annotation('textbox','Position',[.5 0 .5 .05],'VerticalAlignment','bottom','HorizontalAlignment','right',...
             'Color','k','FontSize',6,'LineStyle','none','String',strrep(outfile,'_','\_'));           
  else
    annotation('textbox','Position',[.5 0 .5 .05],'VerticalAlignment','bottom','HorizontalAlignment','right',...
             'Color','k','FontSize',6,'FitHeightToText','on','LineStyle','none','String',strrep(parms.footnote,'_','\_'));
  end
end

if exist(outfile,'file') && ~parms.overwrite
  mmil_logstr(parms,'not saving figure.  File already exists: %s',outfile);
  return;
end

%% questionable block
% orient portrait
% set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation' ,'portrait',...
        'PaperType'        ,'usletter',...
        'PaperUnits'       ,'inches'  ,...
        'PaperPositionMode','auto');  
ppsize = get(gcf,'PaperSize');
width  = ppsize(2) - 2;
height = ppsize(1) - 2;
left   = (ppsize(1) - width ) / 2;
bottom = (ppsize(2) - height) / 2;
set(gcf,'PaperPosition',[left bottom width height]);

if parms.guimode
  pos = get(gcf,'Position');
  if (pos(3)/pos(4)) > (11/8.5)
    pos(3) = (11/8.5)*pos(4);
    set(gcf,'position',pos);
  end
  set(gcf,'PaperPositionMode','manual');  % auto (Ned)
  set(gcf,'InvertHardcopy','off');
  printcmd = {printcmd{:},'-loose'};
end
if ~isempty(maskstring) && strcmp(parms.format,'eps')
  printcmd = {printcmd{:},'-painters'};
end
try 
  mmil_logstr(parms,'Saving figure: %s',outfile);
  print(gcf,printcmd{:},outfile);
catch
  mmil_logstr(parms,'failed to save figure: %s',outfile);
end

function cfg = ezplot2PlotAvg(parms,ieegflag,evlist)
cfg = [];
if ~isempty(parms.layout)
    cfg.layoutfile = parms.layout;
elseif isempty(parms.layout) && ~ieegflag
%  cfg.layoutfile = 'vv_planar.lay';
    cfg.layoutfile = 'vv_planar_nospace.lay';
else
  cfg.layoutfile = [];
end
[cfg.f_path name ext] = fileparts(parms.datafile{1});
cfg.f_name = [name ext];
if iscell(evlist)
  cfg.eventvals = cell2mat(evlist);
else
  cfg.eventvals = evlist;
end
cfg.figpath = parms.outpath;
if parms.save
  cfg.savefigs = parms.format;
else
  cfg.savefigs = 'no';
end
if ~ieegflag, cfg.title = parms.title; end
if ~isempty(parms.conditionkey)
  cfg.condkey = parms.conditionkey;
elseif ~isempty(parms.conditionfile)
  cfg.condkey = parms.conditionfile;
end
if ~isempty(parms.statfile)
  cfg.statistics = parms.statfile;
end

function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
        'blc',[],[],...
        'overlay',[],[],...
        'clip',[],[],...
        'showlabels',[],[],...
        'axes',[],[],...
        'zerolines',[],[],...
        'autoscale',[],[],...
        'newfig',[],[],...
        },false);

if ~isempty(opt.blc)
  if (ischar(opt.blc)    && strcmp(opt.blc,'yes')) ||...
     (isnumeric(opt.blc) && opt.blc==1)
   parms.zscore = 1;
  end
end
if ischar(opt.overlay)
  if     strcmp(opt.overlay,'yes'), parms.overlay = 1;
  elseif strcmp(opt.overlay,'no'),  parms.overlay = 0; end
end
if ischar(opt.newfig)
  if     strcmp(opt.newfig,'yes'), parms.newfig = 1;
  elseif strcmp(opt.newfig,'no'),  parms.newfig = 0; end
end
if ischar(opt.autoscale)
  if     strcmp(opt.autoscale,'yes'), parms.autoscale = 1;
  elseif strcmp(opt.autoscale,'no'),  parms.autoscale = 0; end
end
if isnumeric(opt.clip) && ~isempty(opt.clip)
  if     opt.clip == 1,  parms.clip = 'yes';
  elseif opt.clip == 0,  parms.clip = 'no'; end
end
if isnumeric(opt.showlabels) && ~isempty(opt.showlabels) && opt.showlabels == 1
  parms.showlabels = 'yes';
elseif ischar(opt.showlabels) && strcmp(opt.showlabels,'yes')
  parms.showlabels = 'yes';
end
if isnumeric(opt.axes) && ~isempty(opt.axes)
  if     opt.axes == 1,  parms.axes = 'yes';
  elseif opt.axes == 0,  parms.axes = 'no'; end
end
if isnumeric(opt.zerolines) && ~isempty(opt.zerolines)
  if opt.zerolines == 1
    parms.zerolines = 'yes';
  elseif opt.zerolines == 0
    parms.zerolines = 'no';
  end
elseif ischar(opt.zerolines)
  parms.zerolines = opt.zerolines;
end
