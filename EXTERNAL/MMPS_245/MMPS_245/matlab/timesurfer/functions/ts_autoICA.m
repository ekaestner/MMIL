function epoch_data = ts_autoICA(epoch_data,varargin)
%function epoch_data = ts_autoICA (epoch_data,varargin)
%
% Usage:
%  epoch_data = ts_autoICA (epoch_data,'ICA_ref_chan',ICA_ref_chan,...
%
% Required Input:
%
%  epoch_data  - valid TimeSurfer epoch_data structure
%
% Optional Input:
%   ICA_ref_chan    - name of EOG/EKG/STIM reference channel(s)
%                     (default = 'EOG61')
%   event_codes
%       OR
%   conditions      - the event codes or conditions to inspect (default - all)
%   chantype        - cell array of channel types to process (each set will be
%                     processed individually)
%       OR
%   channels        - list of channels to run on 
%   notch           - flag to run 60 Hz notch filter
%   showcomponents  - flag to show components as rejecting
%
% Output:
%
%   epoch_data - epoch data without the selected components
%
%  If no output specfied script will save epoch_data to outfile.
%
% created:         03/18/08    by Andrei Irimia
% last modified:   05/29/08    by Rajan Patel

% 05/29/08 - Added notch filter and option to show components while
%            rejecting
% 04/14/08 - Allow multiple reference channels and handle correctly
%            Rescaling done as you go along

%% Check Inputs

if nargin < 1, help(mfilename); end

parms = mmil_args2parms(varargin,...
                        { ...
                         'ICA_ref_chan','EOG061',[],...
                         'chantype', 'all', {'all', 'mag' 'grad1' 'grad2' 'eeg', 'other', 'grad', 'meg'},...
                         'channels', [],[],...
                         'event_codes',[],[],...
                         'conditions',[],[],...
                         'showcomponents',false,sort([false true]),...
                         'notch',false,sort([false true]),...
                         'rescale',true,sort([false true]),...
                         'verbose',true,sort([false true]),...                         
                         'logfile',[],[],...
                         'logfid',1,[],...
                        },...
                        false);

epoch_data = ts_checkdata_header(epoch_data);
errors     = ts_checkdata(epoch_data);
if ~isempty(errors)
    mmil_error(parms,'Errors in provided data structure: %s.',errors);
end                   

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
    end;
end;

% Make sure both conditions and event_codes are cell arrays
if (~iscell(parms.event_codes)), parms.event_codes = num2cell(parms.event_codes); end;
if (~iscell(parms.conditions)),  parms.conditions  = num2cell(parms.conditions);  end;

if isempty(parms.channels)
    chans = [];
    if ~iscell(parms.chantype), parms.chantype = {parms.chantype}; end
    for i = 1:length(parms.chantype)
        switch parms.chantype{i}
            case {'mag' 'grad1' 'grad2' 'eeg', 'other'}
                chans = [chans find(strcmp(parms.chantype{i},{epoch_data.sensor_info.typestring}))];
            case {'grad'}
                chans = [chans find(strncmp(parms.chantype{i},{epoch_data.sensor_info.typestring},...
                    length(parms.chantype{i})))];
            case 'meg'
                [a,ch] = find(ismember({epoch_data.sensor_info.typestring}, ...
                    {'mag', 'grad1', 'grad2'}));
                chans  = [chans ch];
            case 'all'
                chans = [chans 1:epoch_data.num_sensors];
        end;
    end
    chans = unique(chans);
else
    chans = parms.channels;
end

%% Run ICA

% find bad channels
bchans  = find([epoch_data.sensor_info.badchan]);
chans   = setdiff(chans,bchans);

if isempty(chans)
    mmil_error(parms,'No good channels for auto ICA rejection.');
end

% find reference channel
if ~iscell(parms.ICA_ref_chan), parms.ICA_ref_chan = {parms.ICA_ref_chan}; end
ICA_ref_chan = [];
for r = 1:length(parms.ICA_ref_chan)
  found = 0;
  c = 1;
  while ~found && ~(c > length(epoch_data.sensor_info))
    if strcmp(epoch_data.sensor_info(c).label,parms.ICA_ref_chan{r})
        ICA_ref_chan = [ICA_ref_chan c];
        found        = 1;
    end
    c = c+1;
  end
  if ~found
      mmil_logstr(parms,'WARNING: Could not find the following reference channel: %s.',parms.ICA_ref_chan{r});
  end
end

if isempty(ICA_ref_chan)
    mmil_error(parms,'Could not find any of the reference channels specified for automatic ICA rejection.');
end

% run through the selected reference channels
for r = 1:length(ICA_ref_chan)
    mmil_logstr(parms,'Current reference channel to be used is %s.',...
        epoch_data.sensor_info(ICA_ref_chan(r)).label);
    proc_chans = setdiff(chans,ICA_ref_chan(r));
    tic
    for k = 1:length(parms.conditions)
        mmil_logstr(parms,'Performing auto-ICA rejection on %.0f channels for condition %.0f.',length(proc_chans),parms.conditions{k});
        [numchan,samples,numtrials] = size(epoch_data.epochs(parms.conditions{k}).data);
        ref_chan = reshape(epoch_data.epochs(parms.conditions{k}).data(ICA_ref_chan(r),:,:),1,samples*numtrials);
        if parms.notch
            ref_chan = permute(ref_chan,[2 1]);
            ref_chan = ts_freq_filt(ref_chan,epoch_data.sfreq,60,5,'notch');
            ref_chan = permute(ref_chan,[2 1]);
        end
        for c = 1:length(proc_chans)
            % mmil_logstr(parms,'Processing channel %s.',epoch_data.sensor_info(proc_chans(c)).label);
            data                  = reshape(epoch_data.epochs(parms.conditions{k}).data(proc_chans(c),:,:),1,samples*numtrials);
            if parms.notch
                data = permute(data,[2 1]);
                data = ts_freq_filt(data,epoch_data.sfreq,60,5,'notch');
                data = permute(data,[2 1]);
            end
            if max(data) ~= 0 && min(data) ~= 0
                data              = [data; ref_chan];
                [weights, sphere] = runica(data,'maxsteps',2,'verbose','off');
                if parms.showcomponents
                    close all
                    trls = 5;
                    IC = weights*sphere*data;
                    figure;
                    subplot(2,2,1);
                    plot(0:1/samples:trls-1/samples,ref_chan(1,1:trls*samples));
                    title(sprintf('Channel %s',epoch_data.sensor_info(ICA_ref_chan(r)).label));
                    subplot(2,2,2);
                    plot(0:1/samples:trls-1/samples,data(1,1:trls*samples));
                    title(sprintf('Channel %s',epoch_data.sensor_info(proc_chans(c)).label));
                    subplot(2,2,3);
                    plot(0:1/samples:trls-1/samples,IC(1,1:trls*samples));
                    title('IC 1');
                    subplot(2,2,4);
                    plot(0:1/samples:trls-1/samples,IC(2,1:trls*samples));
                    title('IC 2');  
                    pause(3);
                end
                rej_data          = icaproj(data,weights*sphere,2,mean(data')'); % leave component # 2
                if parms.rescale
                    epoch_data.epochs(parms.conditions{k}).data(proc_chans(c),:,:) = ...
                        ts_ICAfit_auto(epoch_data.epochs(parms.conditions{k}).data(proc_chans(c),:,:),...
                        reshape(rej_data(1,:),1,samples,numtrials));
                else
                    epoch_data.epochs(parms.conditions{k}).data(proc_chans(c),:,:) = reshape(rej_data(1,:),1,samples,numtrials);
                end
            end
        end
    end
    mmil_logstr(parms,'Time elapsed: %.0f seconds.',toc);
end
mmil_logstr('Done.\n');
    
