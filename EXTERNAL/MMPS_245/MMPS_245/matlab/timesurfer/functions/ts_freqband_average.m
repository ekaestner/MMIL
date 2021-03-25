function tf_wave_data = ts_freqband_average(varargin)
% function outdata = ts_freqband_average(varargin)
% Purpose: average TF power over frequencies
% 
% Inputs (one of the following):
% tf_cont_data  (3D, num_trials = 1; chan x time x freq)
% tf_avg_data   (3D, num_trials > 1; chan x time x freq)
% tf_epoch_data (4D, num_trials > 1; chan x time x freq x trial)
%
% Output:
% tf_wave_data (3D, num_trials >= 1; chan x time x trial)
%
% common_baseline_flag - whether to use one baseline for all trials from
% all conditions or a different baseline for all trials for each condition.

% check inputs & load data if necessary
% Last Mod: 12.09.11 SRD: Skipped tf data doubling in header check, but
%                    compensated in pow calc.  Used looping to save Memory
%                    in power calculation.  More memory saving is possible.

data = [];
if mod(nargin,2)      % odd number of inputs
  data = varargin{1};  %%%%% Bad Move - freq_data large. Call with filename (not the data) to avoid this!  (Penalty: reload of freq data.)
  if nargin > 1
    varargin = varargin(2:end) ;
  else
    varargin = {};
  end
  if ~(issubfield(data,'timefreq.cmplx') && ndims(data.timefreq(1).cmplx)==4)
    try     % load TF trial spectra
      data = ts_load_timefreq_data(data,varargin{:});
    catch   % load average TF spectra
%       warning('Failed to load TF trial spectra!');
      data = ts_checkdata_header(data,varargin{:});
    end
  end
end
% => [tfr] = chan x time x freq x trial  ||  chan x time x freq

% Set up parameters
parms = mmil_args2parms(varargin,...
    {'conditions'   , [],[],...
     'events'       , [],[],...
     'chanlabels'    , [],[],...
     'chantype',      'all', {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
     'foi',[],[],...
     'foilim',[],[],...
     'exc_freqs',[],[],...
     'blc',[],[],...
     'baselinetype','zscore',{'absolute','relative','relchange','zscore'},...
     'baselinefile',[],[],...
     'common_baseline_flag',0,{0,1},...
     'blcwindow',[-inf inf],[],...
     'baseline_data',[],[],...
     'freqcorr',0,{0,1},...
     'freqcorr_alpha',2,[],...
     'toilim',[],[],...
     'rejectfile',[],[],...
     'reject_data',[],[],...
     'verbose',1,{0,1},...
     'logfile',      [],[],...
     'logfid',       [1],[], ...        
    },...
    false);
%        'frequency', 'all',[],...
%      'findex',[],[],...

  
% Backwards compatibility
parms = backcompatible(parms,varargin{:});     

% force baseline correction if common_baseline_flag is set to 1
if parms.common_baseline_flag && ~isequal(parms.blc,'yes')
  parms.blc = 'yes';
end
  
% Select data to process
data  = ts_checkdata_header(data,'precision','skip','events',parms.events,'verbose',0);  %%%%% double --> skip, SRD, 12.6.14  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data  = ts_data_selection(data,'chanlabels',parms.chanlabels,'chantype',parms.chantype,...
  'events',parms.events,'foi',parms.foi,'foilim',parms.foilim); 

% Get info on data
[datatype datafield dataparam] = ts_object_info(data);
if isempty(data) || isempty(data.(datafield)) % was any data selected?
  error('No data selected.');
end
T     = data.(datafield)(1).time;
F     = data.(datafield)(1).frequencies;
Fix   = find(~ismember(F,parms.exc_freqs));
ncond          = length(data.(datafield));
nchan          = data.num_sensors;
ntime          = length(T);
nfreq          = length(Fix);
ntrials        = [data.(datafield).num_trials];
parms.channels = 1:data.num_sensors;
bl_index       = find(T>=parms.blcwindow(1)&T<=parms.blcwindow(2));

% Flags
blc_flag = isequal(parms.blc,'yes')      || isequal(parms.blc,1);
fc_flag  = isequal(parms.freqcorr,'yes') || isequal(parms.freqcorr,1);
use_cmplx_flag = isfield(data.(datafield),'cmplx') && ~isfield(data.(datafield),'power');
 
% initialize output data structure
rmfields = {dataparam{:} 'frequencies'};
tf_wave_data = rmfield(data,datafield); 
tf_wave_data.epochs             = rmfield(data.(datafield),rmfields);
[tf_wave_data.epochs.time] = deal(T);

% remove bad trials but keep bad channels
data = ts_data_selection(data,'reject_data',parms.reject_data,'rejectfile',parms.rejectfile,'keepbadchans_flag',1);

% calculate common baseline if desired
if parms.common_baseline_flag
  if parms.verbose
    mmil_logstr(parms,'Calculating one mean and standard deviation across all trials and conditions for each channel and frequency.'); 
  end
  % concatenate power from all trials in all conditions
 
  for c = 1:ncond
    if use_cmplx_flag
      pow{c} = double(data.(datafield)(c).cmplx(:,bl_index,Fix,:));
      pow{c} = abs(pow{c}).^2;
    else
      pow{c} = data.(datafield)(c).power(:,bl_index,Fix,:);
    end
  end
  % concatenate trials from all conditions
  pow = cat(4,pow{:}); % [chan x time x freq x all trials]
  for ch = 1:nchan
    % calculate mean and stdev for each channel and frequency
    for freq = 1:length(Fix)
      tmp          = pow(1,:,freq,:);
      common_mu(ch,1,freq) = mean(tmp(:));
      common_sd(ch,1,freq) = std(tmp(:));
    end
  end
  clear pow tmp c ch
end

% preallocate result matrix
for c = 1:ncond, tf_wave_data.epochs(c).data = zeros(nchan,ntime,ntrials(c),'single'); end  %%%%% Wasteful iterative memory allocation

% memuse = whos('tf_wave_data'); message = sprintf('Memory used for tf_wave_data reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%% Debug tool

% Freq average
% loop over conditions
freqcorr = (permute(F(Fix).^(parms.freqcorr_alpha),[1 3 2])); % 1 x 1 x freqs   SRD: Create once in case used, moved out of freqcorr loop.  %%%%%%%%% double permute %%%%%%%%%%%%%%%%%%%%%%
for c = 1:ncond
  % TODO: eliminate loop over channels
  for ch = 1:nchan
    if use_cmplx_flag
      pow = double(data.(datafield)(c).cmplx(ch,:,Fix,:)); % This was already doubled?
      pow = abs(pow).^2;
    else
      pow = double(data.(datafield)(c).power(ch,:,Fix,:)); %%%%% Doubled precision here to make up for changed ts_checkdata_header parm, SRD, 12.6.14
    end
    
    % frequency-correction  %%%%% SRD: Modified to loop over pow rather than create an equivalent sized freqcorr matrix
    if fc_flag
        ntri = data.(datafield)(c).num_trials ;
        ntim = size(data.(datafield)(c).time,2) ;
        for ntr=1:ntri %%%%%
            for nti=1:ntim  %%%%%
                for nfrq=1:nfreq
                     pow(1,nti,nfrq,ntr) =pow(1,nti,nfrq,ntr) .* freqcorr(1,1,nfrq); %%%%% End of mod. SRD, 12.9.11 
                end
            end
        end
    end
    if blc_flag
      if parms.common_baseline_flag
        mu = common_mu(ch,:,:);
        sd = common_sd(ch,:,:);
      else
        for freq = 1:length(Fix)
          tmp          = pow(1,bl_index,freq,:);
          mu(1,1,freq) = mean(tmp(:));
          sd(1,1,freq) = std(tmp(:));
        end
      end
      switch parms.baselinetype %%%%% SRD: Should do the same this section, but not tried yet.
        case 'relative'
          shift = 0;
          scale = repmat(mu,[1 ntime 1 ntrials(c)]);
        case 'relchange'
          shift = repmat(mu,[1 ntime 1 ntrials(c)]);
          scale = repmat(mu,[1 ntime 1 ntrials(c)]);
        case 'absolute'
          shift = repmat(mu,[1 ntime 1 ntrials(c)]);
          scale = 1;
        case 'zscore'
          shift = repmat(mu,[1 ntime 1 ntrials(c)]);
          scale = repmat(sd,[1 ntime 1 ntrials(c)]);
        otherwise
          shift = 0;
          scale = 1;
      end
      % shift & scale the data
      pow = (pow - shift) ./ scale;
      
% memuse = whos('shift'); message = sprintf('Memory used for shift reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%%
% memuse = whos('scale'); message = sprintf('Memory used for scale reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%%
% memuse = whos('pow'); message = sprintf('Memory used for pow reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%%

      clear mu sd tmp
    end
    tf_wave_data.epochs(c).data(ch,:,:) = squeeze(mean(single(pow),3));
    
% memuse = whos('tf_wave_data'); message = sprintf('Memory used for tf_wave_data reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%%

    clear tmp mu sd shift scale
  end
end
clear freqcorr  %%%%% SRD: OK to delete now.

% remove bad channels
if exist(parms.rejectfile,'file')
  res = load(parms.rejectfile,'reject_data');
  parms.reject_data = res.reject_data;
  clear res
end
if ~isempty(parms.reject_data)
  parms.reject_data.badtrials = [];
  tf_wave_data = ts_data_selection(tf_wave_data,'reject_data',parms.reject_data);
  
% memuse = whos('tf_wave_data'); message = sprintf('Memory used for tf_wave_data reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%%

end

% whos('data','tf_wave_data','pow','shift','scale'); %%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
        'event_codes',[],[],...
        'blc',[],[],...
        },false);
  
if isempty(parms.events) && ~isempty(opt.event_codes)
  parms.events = opt.event_codes;
end
if ~isempty(opt.blc)
  if isequal(opt.blc,false)
    parms.blc = 'no';
  elseif isequal(opt.blc,true)
    parms.blc = 'yes';
  end
end
