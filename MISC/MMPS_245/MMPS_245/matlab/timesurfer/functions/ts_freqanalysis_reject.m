function [data_out, badtrials] = ts_freqanalysis_reject (data_in,bl_in,varargin)
% function [data_out, badtrials] = ts_freqanalysis_reject (data_in,bl_in,varargin)
%
% Usage: [data_out, badtrials] = ts_freqanalysis_reject (data_in,'option',value...
%
% Rejects trials on a single channel of timefreq data.
%
% Rejection threshold determines z-score from baseline to reject at
%
% Required Inputs:
%
% timefreq_data - a valid timefreq data structure
% baseline_data - a valid timefreq data strucure with baseline 
%
% Optional Parameters
%
% trial_threshold - default = 40
% reject_exclude  - Not Yet Implemented here - only does it by channel
% 
% Outputs:
%
% timefreq_data - data averaged across trials after removing bad trials
% badtrials     - the list of bad baselines
% 
% Created on 04/187/08 by Rajan Patel
% Last Modified on 04/18/08 by Rajan Patel
%% Check Inputs

if nargin < 2, help(mfilename); end

error(nargoutchk(1,2,nargout,'string')); 

parms = mmil_args2parms(varargin,...
                        {...
                         'trial_threshold',40,[],...
                         'reject_exclude','bychannel',{'inall','groupchannels','bychannel'},...
                         'debug_mode',false,sort([false true]),...
                         'logfile',[],[],...
                         'logfid',1,[]...
                         },...
                         false);                    
                     
errors = ts_checkdata(data_in);
if ~isempty(errors)
    mmil_error(parms,'Errors in supplied timefreq data: %s.', sprintf('\t%s\n', errors{:}));
end; 
               
errors = ts_checkdata(bl_in);
if ~isempty(errors)
    mmil_error(parms,'Errors in supplied baseline data: %s.', sprintf('\t%s\n', errors{:}));
end; 

if ~isfield(data_in,'timefreq')
    mmil_error(parms,'Input data must be timefreq data.');
end

if ndims(data_in.timefreq.power) ~= 4
    mmil_error(parms,'Timefreq data should have trial information for rejection.');
end

if size(data_in.timefreq.power,1) ~= size(bl_in.timefreq.power,1)
    mmil_error(parms,'The channels of the two data sets must agree.');
end

if size(data_in.timefreq.power,3) ~= size(bl_in.timefreq.power,3)
    mmil_error(parms,'The number of frequencies do not match in the two data sets.');
end

%% Initialize data_out

data_out                = data_in;
data_out.timefreq.power = [];
data_out.timefreq.data  = [];

%% Perform rejection on a condition by condition and channel by channel basis

bl_end        = bl_in.timefreq.time(end);
bl_mean       = nanmean(bl_in.timefreq.power,2);     
bl_std        = nanstd (bl_in.timefreq.power,[],2);
% adjust threshold based on number of trials in baseline
% you are comparing the mean and standard deviation of multiple
% trials to one trial
warning off MATLAB:divideByZero
scales = 1./sqrt(1./bl_in.timefreq.num_trials);
warning on MATLAB:divideByZero
% make it good for all the channels if the number of trials are consistent
% across all the channels
if length(scales)==1, scales = repmat(scales,[1 length(data_in.sensor_info)]); end

for c=1:length(data_in.timefreq);
    % run on all time points including baseline
    time = 1:size(data_in.timefreq(c).power,2);    
    for ch = 1:length(data_in.sensor_info)
      if ~bl_in.sensor_info(ch).badchan
        badtrials_bychan{ch} = [];       
        % mmil_logstr(parms,'Threshold: %.2f.',parms.trial_threshold * scales(ch));       
        z_scores = (data_in.timefreq(c).power(ch,time,:,:) - ...
                    repmat(bl_mean(ch,:,:),[1 length(time) 1 size(data_in.timefreq(c).power,4)]))...
                    ./  ...
                    repmat(bl_std (ch,:,:),[1 length(time) 1 size(data_in.timefreq(c).power,4)]);
        for t = 1:data_in.timefreq(c).num_trials
            % find z_scores in this trial
            if find(z_scores(:,:,:,t) >= (parms.trial_threshold * scales(ch))), badtrials_bychan{ch}(end+1) = t; end
            if parms.debug_mode
                imagesc(squeeze(z_scores(:,:,:,t))',[0 parms.trial_threshold*scales(ch)]);
                pause
            end
        end
        % mmil_logstr(parms,'Channel %s - number of bad trials: %.0f',data_in.sensor_info(ch).label,length(badtrials_bychan{ch}));
        if size(data_in.timefreq(c).power,4) > length(badtrials_bychan{ch})
            data_out.timefreq(c).power(ch,:,:) = squeeze(mean(data_in.timefreq(c).power(ch,:,:,...
                setdiff(1:size(data_in.timefreq(c).power,4),...
                badtrials_bychan{ch})),4));
            data_out.timefreq(c).data (ch,:,:) = squeeze(mean(data_in.timefreq(c).data(ch,:,:,...
                setdiff(1:size(data_in.timefreq(c).power,4),...
                badtrials_bychan{ch})),4));          
        else %% if there are no trials left set to nans
            mmil_logstr(parms,'WARNING: No good trials left for this channel.');
            data_out.timefreq(c).power(ch,:,:) = nan(1,size(data_in.timefreq(c).power,2),size(data_in.timefreq(c).power,3));
            data_out.timefreq(c).data (ch,:,:) = nan(1,size(data_in.timefreq(c).data,2),size(data_in.timefreq(c).data,3));
        end
      else
        mmil_logstr(parms,'SKIPPING rejection: Channel %s, the baseline data for this channel is all bad.');
        data_out.timefreq(c).power(ch,:,:) = squeeze(mean(data_in.timefreq(c).power(ch,:,:,:),4));
        data_out.timefreq(c).data (ch,:,:) = squeeze(mean(data_in.timefreq(c).data(ch,:,:,:),4)); 
        badtrials_bychan{ch} = [];
      end
    end
    % Do not return a cell of length == 1
    if length(badtrials_bychan)==1, badtrials_percond{c} = badtrials_bychan{1}; end
end

if length(badtrials_percond) == 1, badtrials_percond = badtrials_percond{1}; end

badtrials = badtrials_percond;
