function [data,varargout] = ts_visual_reject (data,varargin)
 
% Usage: [epoch_data] = ts_vis_reject (epoch_data,'option',value....
%        [epoch_data,badchannels] = ts_vis_reject (epoch_data,'option',value....
%        [epoch_data,badchannels,badtrials] = ts_vis_reject (epoch_data,'option',value....
%
% Perform rejection processing on epoch data using fieldtrip's rejectvisual
% function.  The output data will contain NaN's for any bad channels across
% all conditions.  The trials selected as bad will be removed from the
% output data set.  It would be necessary to recalculate the noise matrix
% in the calling function based on the noise time period selected by the
% user.
%
% Required input:
%
%  epoch_data - TimeSurfer epoch data structure.
%
% Optional input:
% 
%  chantype- which channels do you want to look at 
%            'mag','grad1','grad2','eeg','grad', 'all'
%            default = 'all' (all channels but 'other')
%  method  - describes how the data should be shown
%            'summary' - show a single number for each channel and trial
%            (default)
%            'channel' - show the data per channel, all trials at once
%            'trial'   - show the data per trial, all channels at once
%  metric  - metric that should be computed in summar mode for each channel
%            in each trial
%            'var'    - variance within each channel (default)
%            'min'    - minimum value in each channel
%            'max'    - maximum value each channel
%            'absmax' - maximum absolute value in each channel
%            'range'  - range from min to max in each channel
%  alim    - determines amplitude scaling for the channel and trial
%            display, if empty then the scaling is automatic (default = [])
%            can also be one of the following:
%            'max'    - use maximum value over all conds
%            'min'    - use amplitude of the minimum value over all conds
%            'absmax' - use the maximum of the absolute value over all conds
%  ylim    - same as alim
%  latency - [begin end] in seconds or 'maxperlength' (default), 'minperlength', 
%            'prestim' ,'poststim'
%  logfile -
%  logid   -
%  
%  Option Preprocessing Options - These preprocessing steps will be applied
%  to data before viewing but will not be applied to the output data.
%
%  bandpass        - 0/1 turn bandpass filter on or off {default = 0}
%  bandpass_low_cf - the value for the high-pass filter {default = 0}
%  bandpass_high_cf- the value for the low-ass filter {default = 100}
%  detrend_events  - 0/1 to detrend the data {default = 0}
%  baseline_sub    - 0/1 to perform baseline correction {default = 0}
%  baseline_start  - start of baesline in ms {default = -80}
%  baseline_end    - end of baseline in ms {default = -5}
%  highpass        - 0/1 to turn on high pass filter {default = 0}
%  low_cf          - value of high pass filter {default = 0}
%  lowpass         - 0/1 to turn on low pass filter {default = 0}
%  high_cf         - value of low pass filter {default = 100}
%  linefilt        - 0/1 to turn on line filter {default = 0}
%  linefreq        - freq of line noise {default = 60}
%
% Required Output: 
%
%   epoch_data  - epoch_data TimeSurfer data structure
%                 trials selected as bad will be removed from data
%                 channels selected as bad will be replaced by NaN's across
%                 all conditions
%
% Optional Output:
%
%   badchannels - index of bad channels
%   badtrials   - cell array the length of the no. of conditions each with a
%                 list of indices for the trials that were removed for a given condition
%                 (indices are in ref. to original data_in)
%
% Created:       11/14/2007 Rajan Patel
% Last Modified: 05/15/2007 Rajan Patel
%
% See also: rejectvisual 


%% Verify Inputs

if (~mmil_check_nargs(nargin, 1))
    return;
end

opt = mmil_args2parms(varargin,...
                      {'chantype','all',{'all','mag','grad1','grad2','eeg','grad'},...  
                       'method','summary',{'summary','channel','trial'},...        % reject visual options
                       'metric','var',{'var','min','max','absmax','range'},...
                       'alim',[],[],...
                       'ylim',[],[],...
                       'keepchannel','no',{'no','yes','nan'},...
                       'feedback','no',[],...
                       'latency','maxperlength',[],...
                       'bandpass',0,[],...          % preprocessing options
                       'bandpass_low_cf',0,[],...
                       'bandpass_high_cf',100,[],...
                       'detrend_events',0,[],...
                       'baseline_sub',0,[],...
                       'baseline_start',-80,[],...
                       'baseline_end',-5,[],...
                       'highpass',0,[],...
                       'lowpass',0,[],...
                       'low_cf',[],[],...
                       'high_cf',100,[],...
                       'linefilt',0,[],...
                       'linefreq',60,[],...
                       'channel','all',[],...
                       'logfile',      [],       [],...
                       'logfid',       1,        [] ...
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
                       'lnfreq',[],[],...
                       'dftfreq',[],[],...              
                       'lpfiltord',[],[],...
                       'hpfiltord',[],[],...
                       'bpfiltord',[],[],...
                       'bsfiltord',[],[],...
                       'lnfiltord',[],[],...
                       'lpfilttype','but',{'but','fir'},...
                       'hpfilttype','but',{'but','fir'},...
                       'bpfilttype','but',{'but','fir'},...
                       'bsfilttype','but',{'but','fir'},...
                       'lpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'hpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'bpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'bsfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'medianfiltord',[],[],...
                       'blc','no',{'yes', 'no'},...
                       'blcwindow',[-inf 0],[],...
                       'detrend','no',{'yes','no'},...     
                       'polyremoval','no',{'yes','no'},...
                       'polyorder',2,[],...
                       'hilbert','no',{'no','abs','complex','real','imag','absreal','absimag','angle'},...
                       'rectify','no',{'yes','no'},...         
                       'verbose',1,{0,1},...
                       'logfile',      [],[],...
                       'logfid',       [1],[], ...                          
                      },...
                      false);
   
errors = ts_checkdata(data);   
if (length(errors)~=0),     mmil_error(opt, 'Errors supplied data: %s.', sprintf('\t%s\n', errors{:})); end
if ~isfield(data,'epochs'), mmil_error(opt,'Input data must be epoch_data.\n');                         end    
error(nargoutchk(1,3,nargout,'string'));

if ~isempty(opt.ylim), opt.alim = opt.ylim; end
  
cfg = rmfield(opt,{'chantype','bandpass','bandpass_low_cf','bandpass_high_cf',...
  'detrend_events','baseline_sub','baseline_start','baseline_end','highpass','lowpass',...
  'low_cf','high_cf','linefilt','linefreq'});

% processing info
if opt.bandpass       == 1, cfg.bpfilter = 'yes'; cfg.bpfreq    = [opt.bandpass_low_cf opt.bandpass_high_cf]; end
if opt.detrend_events == 1, cfg.detrend  = 'yes'; end
if opt.baseline_sub   == 1, cfg.blc      = 'yes'; cfg.blcwindow = [opt.baseline_start opt.baseline_end]/1000; end
if opt.highpass       == 1, cfg.hpfilter = 'yes'; cfg.hpfreq   = opt.low_cf;   end
if opt.lowpass        == 1, cfg.lpfilter = 'yes'; cfg.lpfreq   = opt.high_cf;  end
if opt.linefilt       == 1, cfg.lnfilter = 'yes'; cfg.lnfreq   = opt.linefreq; end

%% BEGIN JSS CODE
parms = opt;

% determine the initial selection of trials and channels
if strcmp(parms.method,'summary')
  ncond   = 1;
  alldata = data;
  data    = rmfield(alldata,'epochs');
  data.epochs = rmfield(alldata.epochs(1),'data');
  data.epochs.data = cat(3,alldata.epochs.data);
  data.epochs.num_trials = sum([alldata.epochs.num_trials]);
  if isfield(alldata.epochs,'trial_info')
    fld='number';        tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
    fld='latency';       tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
    fld='badtrial';      tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
    fld='event_code';    tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
    fld='duration';      tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
    fld='datafile';      tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
    fld='events_fnames'; tmp=arrayfun(@(x)x.trial_info.(fld),alldata.epochs,'uniformoutput',false); data.epochs.trial_info.(fld)=cat(2,tmp{:}); clear tmp
  end
else
  ncond   = length(data.epochs);
  if ischar(cfg.alim)
    tmp = cat(3,data.epochs.data);
    if strcmpi(cfg.alim,'max')
      cfg.alim = max(tmp(:));
    elseif strcmpi(cfg.alim,'min')
      cfg.alim = min(tmp(:));
    elseif strcmpi(cfg.alim,'absmax')
      cfg.alim = max(abs(tmp(:)));
    else
      cfg.alim = [];
    end
    clear tmp
  elseif length(opt.alim) > 1
    cfg.alim = max(cfg.alim(:)); 
  end
  
end
nchan     = data.num_sensors;
chsel     = [data.sensor_info.badchan] == 0;
badtrl    = {};

for k = 1:ncond
  mmil_logstr(opt,'%s: displaying condition %g of %g\n',mfilename,k,ncond);  
%   fprintf('%s: displaying condition %g of %g\n',mfilename,k,ncond);  
  ntrl   = data.epochs(k).num_trials;
  trlsel = true(1,ntrl);
  trl    = 1:ntrl;
  
  % convert TS to FT
  ftdata         = ts_data2fieldtrip(data,'condition',k,'dimord','chan_time','chantype',parms.chantype,'badchans',find(chsel == 0));
  [chix,ft_chix] = match_str({data.sensor_info.label},ftdata.label); % chix - TS indices to process channels
  
  % loop over rejectvisual
  reply   = 'n';
  while ~strcmpi(reply,'y')
    [ftdata, ft_chsel, ft_trlsel] = rejectvisual(cfg,ftdata); 
    reply = input('Are you happy with your rejection? Y/N: ','s');
%     close;
    % need to add code to update chansel/trlsel...
    trl   = trl(ft_trlsel == 1);
    chsel(chix(~ft_chsel)) = 0;
    [chix,ft_chix] = match_str({data.sensor_info.label},ftdata.label);
  end
  badtrl{k} = sort(setdiff(1:ntrl,trl));
  data.epochs(k).data       = data.epochs(k).data(:,:,trl);
  data.epochs(k).num_trials = length(trl);  
  data.epochs(k).num_rejects.manual = data.epochs(k).num_rejects.manual + length(badtrl{k});
end
badchan = find(chsel == 0);

if strcmp(parms.method,'summary')
  data = alldata;
  clear alldata
  offset = [data.epochs(1:end-1).num_trials 0];
  for k = 1:length(data.epochs)
    trix      = find(trl <= data.epochs(k).num_trials);
    tmptrl    = trl(trix);
    trl       = trl(length(trix)+1:end) - offset(k);
    data.epochs(k).data               = data.epochs(k).data(:,:,tmptrl);
    data.epochs(k).data(badchan,:,:)  = 0;    
    badtrl{k} = sort(setdiff(1:data.epochs(k).num_trials,tmptrl));    
    data.epochs(k).num_trials         = length(tmptrl);        
    data.epochs(k).num_rejects.manual = data.epochs(k).num_rejects.manual + length(badtrl{k});    
  end
  close(gcf)
else
  for k = 1:ncond
    data.epochs(k).data(badchan,:,:) = 0;
  end
end
[data.sensor_info(badchan).badchan] = deal(1);

if isfield(data.epochs,'trial_info') && ~isempty([badtrl{:}])
  for c = 1:length(data.epochs)
    data.epochs(c).trial_info.badtrial(badtrl{c}) = 1;
  end
end
    
if nargout > 1, varargout{1} = badchan; end
if nargout > 2, varargout{2} = badtrl;  end

% Remove bad trials from trial_info (later add keepbadtrials_flag)
% if keepbadtrials_flag
  for k = 1:length(data.epochs)
    % remove rejects from data
    idx = find(data.epochs(k).trial_info.badtrial == 1);
    if isempty(idx) || all(idx==0), continue; end
    data.epochs(k).trial_info.number(idx)        = [];
    data.epochs(k).trial_info.latency(idx)       = [];
    data.epochs(k).trial_info.badtrial(idx)      = [];
    data.epochs(k).trial_info.event_code(idx)    = [];
    data.epochs(k).trial_info.duration(idx)      = [];
    data.epochs(k).trial_info.datafile(idx)      = [];
    data.epochs(k).trial_info.events_fnames(idx) = [];
%     data.epochs(k).data(:,:,idx) = [];
%     data.epochs(k).num_trials = data.epochs(k).num_trials - length(idx);
  end
%   end