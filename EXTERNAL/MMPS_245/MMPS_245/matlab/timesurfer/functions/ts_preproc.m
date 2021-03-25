function data = ts_preproc(data,varargin)
% Purpose: proprocess epoch data and return epoch or average data
% 
% Usage: data = ts_preproc(data,'key1',val1, ... )
%
% Example: epoch_data = (epoch_data,...
%                             'bandpass_flag'   ,1 ,...
%                             'bandpass_low_cf' ,.5,...
%                             'bandpass_high_cf',40,...
%                             'bandpass_low_tb' ,.2,...
%                             'bandpass_high_tb',5 ,...);
%
% Inputs: 
%     data: timesurfer data struct of epoch data 
%
% Outputs:
%     data: timesurfer data struct of epoch data
%       -OR-
%     data: timesurfer data struct of average data
%
% Parameters: default : options :  description 
%     saveepochs_flag: 0 : 0,1 : whether to save epoch data 
%     saveaverages_flag: 0 : 0,1 : whether to save average data 
%     returnepochs_flag: 1 : 0,1 : whether to return epoch data 
%     returnaverages_flag: 0 : 0,1 : whether to return average data 
%     prefix: 'preproc' :: specifies a prefix to use for all output file names
%     rootoutdir: pwd :: All output files will be saved in the rootoutdir directory (default is working directory)
%     verbose: 1 : 1,0 : whether to produce fulll output
%     logfile: [] :: filename for log file (stdout will be appended)
%     logfid: 1 :: file FID for log file
%     timesurfer_flag: 1 : 0,1 : whether to use timesurfer preprocessing 
%     fieldtrip_flag: 0 : 0,1 : whether to use fieldtrip preprocessing 
%        Timesurfer Paramers:
%     calc_ncov_flag: 0 : 0,1 : whether to calculate/recalculate noise
%                               covariance matrix after epoching/re-epoching
%     bandpass_flag: 0 : 0,1 : whether to apply a bandpass filter 
%     bandpass_low_cf: 0.2 :: low cutoff frequency (high-pass filter) (Hz)
%     bandpass_low_tb: 0 :: low cutoff transition band (Hz)
%     bandpass_high_cf: 30 :: high cutoff frequency (low-pass filter) (Hz)
%     bandpass_high_tb: 0 :: high cutoff transition band (Hz)
%     bandpass_detrend_flag: 1 : 0,1 : whether to subtract a linear fit before filtering
%     notch_flag: 0 : 0,1 : whether to apply a notch filter
%     notch_cf: [] :: cutoff frequency (notch filter) (Hz)
%     notch_tb: 5 :: notch filter transition band (Hz)
%     lowpass_flag: 0 : 0,1 : whether to apply a lowpass filter  
%     lowpass_cf: [] :: high cutoff frequency (low-pass filter) (Hz)
%     lowpass_tb: 0 :: high cutoff transition band (Hz)
%     highpass_flag: 0 : 0,1 : whether to apply a highpass filter 
%     highpass_cf: [] :: low cutoff frequency (high-pass filter) (Hz)
%     highpass_tb: 0 :: low cutoff transition band (Hz)
%     dsfact: 1 : downsampling factor -- must be an integer
%                 data is downsampled to lower sampling frequency
%                 e.g. original sampling freq = 1000, dsfact = 4,
%                 resulting sampling freq = 250
%     detrend_flag: 0 : 0,1 : whether to detrend single trials before averaging
%     baseline_flag: 0 : 0,1 : whether to subtract a linear fit before filtering
%     blcwindow: [] :: [begin end] in seconds, the default is the complete trial
%     combinations: [] :: this is the list of combinations to produce 
%                         provided as a cell array.  Each combination is
%                         supplied as a string.  The numbers within the string are
%                         treated either as event codes (if using the reference =
%                         'events') or as condition numbers (if using the reference =
%                         'conditions'), for more information see ts_combine_conditions
%     comboeventcodes: [] :: a list of new event codes to assign each of the new
%                              combinations, if none is supplied it will default creating new event
%                              codes starting with the largest current event code in the data
%                              structure conditions.  It is recommended to specify your own.
%     calc: 'weighted' : 'weighted','avg','sum' : specify whether to perform
%                                                 a weighted average, straight average,  
%                                                 or addition conditions with '+' operation
%     cfg: [] :  a structure specifying the configuration for MC clustering
%     feedback: 'none' : 'non','gui','dial','textbar','text','textcr','textnl','none','no','yes' :
%        Fieldtrip Parameters:
%     lpfilter: 'no' : 'no','yes' :  lowpass filter
%     hpfilter: 'no' : 'no','yes' : highpass filter
%     bpfilter: 'no' : 'no','yes' :  bandpass filter
%     bsfilter: 'no' : 'no','yes' : bandstop filter
%     lnfilter: 'no' : 'no','yes' : line noise removal using notch filter
%     dftfilter: 'no' : 'no','yes' : line noise removal using discrete fourier transform
%     medianfilter: 'no' : 'no','yes' : jump preserving median filter
%     lpfreq: 30 :: lowpass frequency in Hz
%     hpfreq: [] :: highpass frequency in Hz
%     bpfreq: [] :: bandpass frequency range, specified as [low high] in Hz
%     bsfreq: [] :: bandstop frequency range, specified as [low high] in Hz 
%     lnfreq: 60 :: line noise frequency in Hz
%     dftfreq: [60 120 180] : line noise frequencies for DFT filter in Hz
%     lpfiltord: 6 :: lowpass  filter order
%     hpfiltord: 6 :: highpass filter order
%     bpfiltord: 4 :: bandpass filter order
%     bsfiltord: 4 :: bandstop filter order
%     lnfiltord: 4 :: line noise notch filter order
%     lpfilttype: 'but' : 'but','fir' : lowpass digital filter type
%     hpfilttype: 'but' : 'but','fir' : highpass digital filter type
%     bpfilttype: 'but' : 'but','fir' : bandpass digital filter type
%     bsfilttype: 'but' : 'but','fir' : bandstop digital filter type
%     lpfiltdir: 'twopass' : 'onepass','onepass-reverse',twopass' : lowpass filter direction
%     hpfiltdir: 'twopass' : 'onepass','onepass-reverse',twopass' : highpass filter direction
%     bpfiltdir: 'twopass' : 'onepass','onepass-reverse',twopass' : bandpass filter direction
%     bsfiltdir: 'twopass' : 'onepass','onepass-reverse',twopass' : bandstop filter direction
%     medianfiltord: 9 :: length of median filter
%     blc: 'no': 'no','yes' : whether to subtract a linear fit before filtering
%     blcwindow: [] :: [begin end] in seconds, the default is the complete trial
%     detrend: 'no' : 'no','yes' : whether to detrend single trials before averaging
%     polyremoval: 'no' : 'no','yes' : this is done on the complete trial
%     polyorder: 2 :: polynome order
%     hilbert: 'no': 'np','abs','complex','real','imag','absreal','absimag','angle' :
%     rectify: 'no' : 'no','yes' : 
%     precision: 'double' : 'single','double' : whether to save data with
%     single precision instead of double
%        EEG only preprocessing options:
%     montage : 'no' : 'no','yes' :
%     reref: 'no': 'no',yes' : 
%     refchannel: [] :: cell-array with new EEG reference channel(s)
%     implicitref: [] : 'channel_label': add the implicit EEG reference as zeros
%     montage : 'no',montage_stucture  : montage structure 
%
% Created:  04/10/09 by Jason Sherfey
% Last Mod: 08/11/14 by Don Hagler
%

parms = mmil_args2parms(varargin,...
						{'timesurfer_flag',1,{0,1},...
             'events',[],[],...
             'fieldtrip_flag',0,{0,1},...
             'calc_ncov_flag',false,[false true],...
             'bandpass_flag',false,[false true],...
             'bandpass_detrend_flag',true,[false true],...
             'bandpass_baseline_flag',false,[false true],...
             'bandpass_low_cf',0.2,[],...
             'bandpass_low_tb',0,[],...
             'bandpass_high_cf',30,[],...
             'bandpass_high_tb',0,[],...
             'notch_flag',false,[false true],...
             'notch_cf',[],[],...
             'notch_tb',5,[],...
             'lowpass_flag',false,[false true],...
             'lowpass_cf',[],[],...
             'lowpass_tb',0,[],...
             'highpass_flag',false,[false true],...
             'highpass_cf',[],[],...
             'highpass_tb',0,[],...
             'dsfact',1,[],...
             'detrend_flag',false,[false true],...
             'baseline_flag',false,[false true],...
             'baseline_start',-Inf,[-Inf,Inf],...
             'baseline_end',Inf,[-Inf,Inf],...     
             'combinations',[],[],...
             'comboeventcodes',[],[],...
             'neweventcodes',[],[],...
             'calc','weighted',{'weighted','avg','sum'},...
             'cfg',[],[],...
             'feedback','none',{'non','gui','dial','textbar','text','textcr','textnl','none','no','yes'},...
             'lpfilter','no',{'yes','no'},...
             'hpfilter','no',{'yes','no'},...
             'bpfilter','no',{'yes','no'},...
             'bsfilter','no',{'yes','no'},...
             'lnfilter','no',{'yes','no'},...
             'dftfilter','no',{'yes','no'},...
             'medianfilter','no',{'yes','no'},...
             'lpfreq',30,[],...
             'hpfreq',[],[],...
             'bpfreq',[],[],...
             'lnfreq',60,[],...
             'dftfreq',[60 120 180],[],...              
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
             'blc','no',[],...
             'blcwindow',[],[],...
             'detrend','no',{'yes','no'},...     
             'polyremoval','no',{'yes','no'},...
             'polyorder',2,[],...
             'hilbert','no',{'no','abs','complex','real','imag','absreal','absimag','angle'},...
             'rectify','no',{'yes','no'},...
             'precision',[],{'single','double'},...
             'reref','no',{'yes','no'},...
             'refchannel',[],[],...
             'implicitref',[],[],...
             'montage','no',[],...
             'saveepochs_flag',0,{0,1},...
             'saveaverages_flag',0,{0,1},...
             'returnepochs_flag',1,{0,1},...
             'returnaverages_flag',0,{0,1},...
             'saveavgs','no',{'yes','no'},...
             'saveepochs','no',{'yes','no'},...
             'returnavgs','no',{'yes','no'},...
             'returnepochs','no',{'yes','no'},...
             'rootoutdir',pwd,[],...
             'prefix','preproc',[],...
             'filename',[],[],...
             'overwrite',0,{0,1},...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...                
						 },false);

data  = ts_checkdata_header(data,'events',parms.events);

% store original data structure parms
if isfield(data,'parms'), parms.previous = data.parms; end

[datatype,datafield,dataparam] = ts_object_info(data,varargin{:});
% data = ts_data_selection(data,varargin{:});

if isempty(parms.comboeventcodes), parms.comboeventcodes = parms.neweventcodes; end
if isempty(parms.precision),       parms.precision = class(data.(datafield)(1).data); end
if strcmp(parms.saveavgs,'yes'),   parms.saveaverages_flag = 1; end
if strcmp(parms.saveepochs,'yes'), parms.saveepochs_flag = 1;   end
if strcmp(parms.returnavgs,'yes'),   parms.returnaverages_flag = 1; end
if strcmp(parms.returnepochs,'yes'), parms.returnepochs_flag = 1;   end

avgflag = parms.saveaverages_flag || parms.returnaverages_flag;
epoflag = parms.saveepochs_flag || parms.returnepochs_flag;

mmil_logstr(parms,'preprocessing %g conditions, %g channels',length(data.(datafield)),data.num_sensors);

%% TimeSurfer preprocessing    
if parms.timesurfer_flag && ~parms.fieldtrip_flag
  % convert fieldtrip parameters to timesurfer parameters
  if parms.highpass_flag
      parms.hpfilter = 'yes';
  end
  if parms.baseline_flag || strcmp(parms.blc,'yes')
    if strcmp(parms.blc,'yes'), 
      parms.baseline_flag = 1;
    end
    if isnumeric(parms.blcwindow) && ~isempty(parms.blcwindow)
      parms.baseline_start = parms.blcwindow(1);
      parms.baseline_end   = parms.blcwindow(end);
    elseif ischar(parms.blcwindow) && strcmp(parms.blcwindow,'all')
      parms.baseline_start = -Inf;
      parms.baseline_end   =  Inf;
    else
      parms.baseline_start = parms.baseline_start / 1000;
      parms.baseline_end   = parms.baseline_end / 1000;
    end
  end
  if parms.lowpass_flag || strcmp(parms.lpfilter,'yes')
    if strcmp(parms.lpfilter,'yes')
      parms.lowpass_flag = 1;
    end
    if isempty(parms.lowpass_cf)
      parms.lowpass_cf = parms.lpfreq;
    end
  end
  if parms.bandpass_flag || strcmp(parms.bpfilter,'yes')
    if strcmp(parms.bpfilter,'yes')
      parms.bandpass_flag = 1;
    end
    if ~isempty(parms.bpfreq)
      parms.bandpass_low_cf  = parms.bpfreq(1);
      parms.bandpass_high_cf = parms.bpfreq(2);
    end    
  end
  if parms.notch_flag || strcmp(parms.lnfilter,'yes')
    if strcmp(parms.lnfilter,'yes')
      parms.notch_flag = 1;
    end
    if isempty(parms.notch_cf)
      parms.notch_cf = parms.lnfreq;
    end
  end  
  if strcmp(parms.detrend,'yes')
    parms.detrend_flag = 1;
  end
  
  % do anything?
  proc_flag = parms.baseline_flag || (parms.dsfact~=1) || parms.lowpass_flag || ...
    parms.bandpass_flag || strcmpi(parms.hpfilter,'yes') || parms.notch_flag || ...
    parms.detrend_flag || parms.baseline_flag || avgflag;
  
  % loop over conditions
  for c = 1:length(data.(datafield))
    if ~proc_flag || (~isempty(parms.events) && ~ismember(data.(datafield)(c).event_code,parms.events))
      continue;
    end
    mmil_logstr(parms,'preprocessing %g of %g',c,length(data.(datafield)));
%     if ischar(parms.blcwindow) && strcmp(parms.blcwindow,'all')
%       parms.baseline_start = data.(datafield)(c).time(1);
%       parms.baseline_end   = data.(datafield)(c).time(end);
%     end
    % get dimensions for preallocation
    sx = size(data.(datafield)(c).data,1);
    sy = size(data.(datafield)(c).data,2);
    sz = size(data.(datafield)(c).data,3);
    if parms.dsfact ~= 1
      data.(datafield)(c).time = downsample(data.(datafield)(c).time,parms.dsfact);
      sy = length(data.(datafield)(c).time);
    end
    % preallocate space for the preprocessed data
    tmp = zeros(sx,sy,sz);
    % convert to single precision
    if ~strcmp(parms.precision,'double')
      tmp = cast(tmp,parms.precision);
    end    
    % loop over trials
    for trl = 1:sz
      % extract one trial
      dat = squeeze(data.(datafield)(c).data(:,:,trl));
      if trl == sz
        data.(datafield)(c).data = [];
      end
      % force to double precision
      if ~isa(dat,'double')
        dat = double(dat);
      end
      % lowpass filter
      if parms.lowpass_flag
        dat = ts_freq_filt(dat',data.sfreq,parms.lowpass_cf,parms.lowpass_tb,'lowpass')';
      end      
      % bandpass filter
      if parms.bandpass_flag
        if parms.bandpass_baseline_flag
          dat = blc(dat);
        end        
        szthresh = 100; % arbitrary size threshold
        if numel(dat)*8*2/10E6 > szthresh
          % loop over channels to handle large matrices
          if parms.bandpass_detrend_flag
            for ch = 1:size(dat,1)
              dat(ch,:) = detrend(dat(ch,:)')';
            end
          end
          for ch = 1:size(dat,1)
            dat(ch,:) = ts_freq_filt( dat(ch,:)',data.sfreq,[parms.bandpass_low_cf,parms.bandpass_high_cf],...
                                [parms.bandpass_low_tb,parms.bandpass_high_tb],'bandpass' )';
          end
        else
          if parms.bandpass_detrend_flag
            dat = detrend(dat')';
          end
          dat = ts_freq_filt( dat',data.sfreq,[parms.bandpass_low_cf,parms.bandpass_high_cf],...
                              [parms.bandpass_low_tb,parms.bandpass_high_tb],'bandpass' )';        
        end
      end
      % highpass filter
      if strcmpi(parms.hpfilter,'yes')
        dat = ts_freq_filt(dat',data.sfreq,parms.hpfreq,parms.highpass_tb,'highpass')';
      end
      % notch filter
      if parms.notch_flag
        for fc = 1:length(parms.notch_cf)
          dat = ts_freq_filt(dat',data.sfreq,parms.notch_cf(fc),parms.notch_tb,'notch')';
        end
      end
      % convert to single precision
      if ~strcmp(parms.precision,'double')
        dat = cast(dat,parms.precision);
      end      
      % downsample
      if parms.dsfact ~= 1
        dat = downsample(dat',parms.dsfact)'; 
      end
      % linear trend removal
      if parms.detrend_flag
        dat = detrend(dat')';
      end
      % baseline correction
      if parms.baseline_flag
        bid = [nearest(data.(datafield)(c).time,parms.baseline_start),...
               nearest(data.(datafield)(c).time,parms.baseline_end)];
        bmu = mean(dat(:,bid(1):bid(2)),2);
        dat = dat - bmu*ones(1,size(dat,2));
        clear bmu bid
      end
      tmp(:,:,trl) = dat;
      clear dat
    end
%     % convert to single precision
%     if ~strcmp(parms.precision,'double')
%       tmp = cast(tmp,parms.precision);
%     end
    % remove data from structure prior to copy to prevent out of memory error
    data.(datafield)(c).data = [];
    data.(datafield)(c).data = tmp;
    clear tmp
  end
  data.sfreq = data.sfreq / parms.dsfact;
  if avgflag
    data = ts_trials2avg(data);
  end
end

%% FieldTrip preprocessing
if parms.fieldtrip_flag
  n = 1;
  for c = 1:length(data.(datafield))
    if ~isempty(parms.events) && ~ismember(data.(datafield)(c).event_code,parms.events)
      continue;
    end
    mmil_logstr(parms,'preprocessing %g of %g',c,length(data.(datafield)));
    if ~isempty(parms.cfg)
      cfg = parms.cfg;
    else
      cfg = parms;
    end
    % force to double precision
    if ~isa(data.(datafield)(c).data,'double')
      data.(datafield)(c).data = double(data.(datafield)(c).data);
    end
    tmp = ts_data2fieldtrip(data,'condition',c,'dimord','chan_time');
    if avgflag
      cfg.keeptrials = 'no';
    else
      cfg.keeptrials = 'yes';
    end
    ftout{n} = timelockanalysis(cfg,tmp);
    clear tmp;
    if ~strcmp(parms.precision,'double')
      tmp = cast(tmp,parms.precision);
    end
    n = n + 1;
  end
  data = rmsubfield(data,'epochs.data');
  if avgflag
    data = ts_fieldtrip2data(ftout,'averages',data);
  else
    data = ts_fieldtrip2data(ftout,'epochs',data);
  end
  clear ftout
end

% calculate noise covariance matrix
if parms.calc_ncov_flag
  data = ts_calc_ncov(data,varargin{:}); 
end

% calculate composite conditions
if ~isempty(parms.combinations)
  data = ts_combine_conditions(data,'combinations',parms.combinations,'neweventcodes',parms.comboeventcodes,'calc',parms.calc);
end

if isfield(parms,'previous'), data.parms = parms; end

if parms.saveepochs_flag && ~avgflag
  for c = 1:length(data.epochs)
    epoch_data = rmfield(data,'epochs');
    epoch_data.epochs(1) = data.epochs(c);
    if isempty(parms.filename) || ~ischar(parms.filename)
      parms.filename = sprintf('%s/%s_epochs_cond%g.mat',parms.rootoutdir,parms.prefix,c);
    end
    if exist(parms.filename) && ~parms.overwrite
      mmil_logstr(parms,'not overwriting %s',parms.filename);
    else
      mmil_logstr(parms,'saving epoch_data to %s',parms.filename);
      save(parms.filename,'epoch_data');
    end
    clear epoch_data;
  end
end
if parms.saveaverages_flag
%   data = ts_trials2avg(data);
  for c = 1:length(data.averages)
    avg_data = rmfield(data,'averages');
    avg_data.averages(1) = data.averages(c);
    if isempty(parms.filename) || ~ischar(parms.filename)
      parms.filename = sprintf('%s/%s_averages_cond%g.mat',parms.rootoutdir,parms.prefix,c);
    end
    if exist(parms.filename) && ~parms.overwrite
      mmil_logstr(parms,'not overwriting %s',parms.filename);
    else
      mmil_logstr(parms,'saving avg_data to %s',parms.filename);
      save(parms.filename,'avg_data');
    end
    clear avg_data;
  end  
end
