function timefreq_data = ts_timefreq_fieldtrip(epoch_data, varargin)
%function ts_timefreq( indata, varargin )
%
% Purpose:
%   Do time/frequency power analysis
%
% Required Param:
%   epoch_data: file from which to read raw data
%
% Optional Params:
%
%   foi:     frequencies of interest (frequencies at which to calculate
%            the analysis of your choice
%            {default: 1:1:60Hz}
%   toi:     epoch over which to analyze your data (in s)
%            {default: (epoch_data's epoch length)}
%   method:  algorithm for calculating data {wltconvol, tfr, mtmconvol,
%            wltconvol_coeff, tf_coeff, mtmconvol_coeff}.  See freqanalysis
%            for definition of the first 3; the second 3 use custom code to
%            extract the raw complex spectral coefficients (as opposed to
%            calculating power directly)
%            {default: wltconvol}
%   waveletwidth: controls size of curve to convolve data with
%            {default: 7}
%   downSampleRate:
%            {default: 1}
%   cfg:     raw fieldtrip parameters to additionally specify (for advanced
%            users)
%
%   conditions: conditions to process
%            {default: all}
%
%   channels: channels (by #) to process
%            {default: all}
%      OR
%   chantype: channel type to process
%            {default: all}
%
%   reject_flag: whether to show interactive tools for manual rejection
% 
%   verbose:  turning this option on makes execution slower
%
%
%  Created by:        Ben Cipollini  07/10/2007
%  Last Modified by:  A   Irimia     02/24/2009
%
% See also: freqanalysis, freqanalysis_wltconvol, freqanalysis_tfr,
% freqanalysis_mtmconvol

MIN_SRATE_TO_FOI_RATIO = 3; % minimum samplerate-to-max(foi)-ratio

if (~mmil_check_nargs(nargin, 1))
    return;
end;

% First, scrub and validate the params.
parms           = mmil_args2parms( varargin, ...
   {'foi',           1:20,        [],...
    'toi',           [],          [],...
    'method',        'wltconvol', [],...
    'cfg',           [],          [],...
    ...
    'event_codes',   {},          [],...
    'conditions',    [],          [],...
    'downSampleRate',1,           [],... %no downsampling by default
    'channels',      [],          [],...
    'chanType',      [],          [],...
    ...
    'reject_flag',   false,       sort([false true]),...
    ... %I/O args
    'logfile',       '',          [],...
    ... %user prefs
    'verbose',       false,        sort([false true]), ...
    'logfile',       [],          [],...
    'logfid',        1,           [],...
    'trials_flag',   0,           sort([false true]),...
    'verbose',       true,        sort([false true]),...
    'sf',            ones(1,20),  [],...
    'gwidth',        ones(1,20).*pi, []...
    },...
    false);

%%%%%%%%%%%%%%%%%%%%%% PARAMETER CHECKING AND DEFAULTS %%%%%%%%%%%%%%%%%%%
% Check input data type
if (~strcmp(ts_objecttype(epoch_data), 'epochs')), mmil_error(parms, 'Input object must be a valid epoch data structure.'); end;

% Check channel selections
if (~isempty(parms.chanType) && ~isempty(parms.channels)), 
    mmil_error(parms, 'Cannot specify BOTH chanType AND channels param values.');
    % Choose all channels by default
elseif (isempty(parms.chanType) && isempty(parms.channels))
    parms.channels = 1:epoch_data.num_sensors;
    % Convert all inputs directly to channel indices
elseif (isempty(parms.channels))
    parms.channels = find(strcmp({epoch_data.sensors.type}, parms.chanType));
end;

% Make event codes a cell array
if (~isempty(parms.event_codes) && ~isempty(parms.conditions))
    mmil_error(parms, 'Cannot specify BOTH event_codes AND conditions.');
    %specified none
elseif (isempty(parms.event_codes) && isempty(parms.conditions))
    parms.event_codes = { epoch_data.epochs.event_code };
    parms.conditions  = num2cell(1:length(epoch_data.epochs));
    % specified conditions
elseif (isempty(parms.event_codes))
    if (min(parms.conditions) < 1 || max(parms.conditions) > length(epoch_data.epochs))
        mmil_error(parms, 'Conditions are out of range; nConditions=%d', length(epoch_data.epochs));
    else
        parms.event_codes = { epoch_data.epochs(parms.conditions).event_code };
    end;
    % specified event_codes
else
    if (~isempty(setdiff([parms.event_codes{:}], [epoch_data.epochs.event_code])))
        mmil_error(parms, 'Event code doesn''t exist in epoch data: %d.', ...
            min(setdiff([parms.event_codes{:}], [epoch_data.epochs.event_code])));
    else
        [a,parms.conditions]=intersect([epoch_data.epochs.event_code], [parms.event_codes{:}]);
        parms.conditions    = num2cell(parms.conditions);
    end
end

% Make sure both conditions and event_codes are cell arrays
if (~iscell(parms.event_codes)), parms.event_codes = num2cell(parms.event_codes); end;
if (~iscell(parms.conditions)),  parms.conditions  = num2cell(parms.conditions);  end;

% Default: downsample to a min of 3x max foi
% Note: this code is currently disabled (default downsamplerate is not empty, but 1)
if isempty(parms.downSampleRate)
    parms.downSampleRate  = floor(epoch_data.sfreq/(MIN_SRATE_TO_FOI_RATIO*max(parms.foi)));
end;

% Verify that the sample frequency is high enough to determine the FOI
if (parms.downSampleRate==0 || ...
        floor(epoch_data.sfreq/(MIN_SRATE_TO_FOI_RATIO*max(parms.foi))) < parms.downSampleRate)
    fprintf(parms.logfid, '%s: ERROR: can''t compute frequencies up to %7.2f Hz with sfreq=%7.2f\n', ...
        mfilename, max(parms.foi), epoch_data.sfreq);
    fprintf('      for sfreq=%7.2f, max foi=%7.2f\n', epoch_data.sfreq, epoch_data.sfreq/MIN_SRATE_TO_FOI_RATIO);
    return;
end;

if parms.verbose,
    mmil_logstr(parms, 'Running t/f analysis, method=%s, on %d frequencies and %d times.',...
        parms.method, length(parms.foi), length(parms.toi));
end

% Force ts_timefreq to run ONLY with one event code at a time.
if (length(parms.event_codes) > 1)
    for i=1:length(parms.event_codes)
        newparms   = rmfield(parms, 'conditions');
        newparms.event_codes = parms.event_codes(i);
        newargs    = mmil_parms2args(newparms);
        if (~exist('timefreq', 'var'))
            timefreq        = ts_timefreq_fieldtrip(epoch_data, newargs{:});
        else
            timefreq(end+1) = ts_timefreq_fieldtrip(epoch_data, newargs{:});
        end;
    end;
    timefreq = ts_combine_data(timefreq);
    return;
end;

% at this point we're guaranteed exactly one event code,
% so let's just re-arrange the struct
parms.event_code = parms.event_codes{1};
parms.condition  = parms.conditions{1};
parms            = rmfield(parms, {'event_codes', 'conditions'});

% convert output object back to current data type
epoch_idx           = find([epoch_data.epochs.event_code] == parms.event_code);
cond_epoch          = epoch_data;
cond_epoch.epochs   = cond_epoch.epochs(epoch_idx);

% Finally... we can default toi
if (isempty(parms.toi)), parms.toi = cond_epoch.epochs.time; end;

% prepare single trial data
if parms.verbose, mmil_logstr(parms, 'Converting to Fieldtrip structure.'); end
ft_epochs  = ts_data2fieldtrip_08_11_07(epoch_data, 'dimord', 'chan_time', 'condition',  parms.condition);

cfg            = [];
cfg.keeptrials = 'yes';
cfg.channel    = {epoch_data.sensor_info(parms.channels).label};
cfg.removemean = 'no';
ft_tla         = timelockanalysis(cfg, ft_epochs);
ft_tla.trial   = single(ft_tla.trial);

% Visual rejection
if (parms.reject_flag)
    cfg             = [];
    cfg.channel     = num2cell(parms.channels);
    cfg.keepchannel = 'yes';
    cfg.method      = 'channel';
    rej_data        = ft_tla;
    rej_ok          = 2;
    while rej_ok == 2
        %rej_data = timelock2preprocessed(FTstructure_timelock);
        rej_data = rejectvisual(cfg,rej_data);
        rej_ok   = menu('Happy with artifact rejection?', 'Yes','No');
    end;
    ft_tla = rej_data;
end;

% time-frequency analysis (wavelet convolution) with freqanalysis_tfr
if parms.verbose, mmil_logstr(parms, 'Doing t/f analysis'); end
cfg            = parms.cfg;
cfg.channel    = 'all'; % we already limited channels above
cfg.keeptrials = 'no';
cfg.sf         = parms.sf;
cfg.gwidth     = parms.gwidth;

switch (parms.method)
    case {'wltconvol', 'tfr' 'mtmconvol'}
        cfg.method          = parms.method;
        cfg.output          = 'pow';
    case {'wltconvol_coeff', 'tfr_coeff'}
        cfg.method          = parms.method(1:length(parms.method)-length('_coeff'));
        cfg.output          = 'coeff';
    case {'mtmconvol_coeff'}
        cfg.method          = parms.method(1:length(parms.method)-length('_coeff'));
        cfg.output          = 'fourier';
    otherwise, mmil_error(parms, 'Unknown timefreq method: %s', parms.method);
end;

switch(cfg.method)
    case {'mtmconvol'}
        cfg.foi             = parms.foi;
        cfg.toi             = parms.toi;%(1):1/ft_tla.fsample:parms.toi(2);
        if (strcmp(cfg.output, 'pow'))
            cfg.t_ftimwin       = repmat( min(parms.waveletwidth./parms.foi), [length(parms.foi) 1]);
            %cfg.tapsmofrq       = min(cfg.foi/8);
            cfg.tapsmofrq       = 1./cfg.t_ftimwin;
            cfg.keeptrials      = 'no';
            cfg.keeptapers      = 'no';
        else
            cfg.t_ftimwin       = repmat( min(parms.waveletwidth./parms.foi), [length(parms.foi) 1]);
            %cfg.tapsmofrq       = min(cfg.foi/8);
            cfg.tapsmofrq       = 1./cfg.t_ftimwin;
            cfg.keeptrials      = 'yes';
            cfg.keeptapers      = 'yes';
        end;
        cfg.channelcmb          = {[], []};
    otherwise
        % set configuration for freqanalysis_tfr
        cfg.foi                 = parms.foi;
        cfg.toi                 = parms.toi;
end

% calculate power spectrum
ft_freq    = freqanalysis(cfg, ft_tla, parms.trials_flag);

% IF statement below removed by AI -- 02/24/2009 -- if trials_flag is off,
% if statement presence generates error
% if parms.trials_flag
    tmp           = ft_freq;
    timefreq_data = ts_fieldtrip2timefreq(tmp, cond_epoch); clear tmp;
    sz            = size(ft_freq.powspctrm);
    for k = 1:size(ft_freq.powspctrm,4)
        tmp = ft_freq;
        tmpfreq_data  = ts_fieldtrip2timefreq(tmp, cond_epoch); clear tmp;
        szz = min([size(timefreq_data.timefreq.power,1) size(tmpfreq_data.timefreq.power,1)]);
        if length(size(timefreq_data.timefreq.power)) > 4
            timefreq_data.timefreq.power(k,:,:,:) = single(tmpfreq_data.timefreq.power(1:szz,:,:));
            timefreq_data.timefreq.cmplx(k,:,:,:) = single(tmpfreq_data.timefreq.cmplx(1:szz,:,:));
        else
            timefreq_data.timefreq.power = single(tmpfreq_data.timefreq.power(1:szz,:,:));
            timefreq_data.timefreq.cmplx = single(tmpfreq_data.timefreq.cmplx(1:szz,:,:));            
        end
    end
% end
timefreq_data.timefreq.power = permute(timefreq_data.timefreq.power,[1 3 2]);
timefreq_data.timefreq.cmplx = permute(timefreq_data.timefreq.cmplx,[1 3 2]);
% Do post-processing rejections
if (parms.reject_flag), 
    if parms.verbose, mmil_logstr(parms, 'NOTICE: post-processed TFR rejections NYI.'); end
end
if parms.verbose, mmil_logstr(parms, 'Done.'); end
return
