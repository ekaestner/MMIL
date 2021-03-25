function data = ts_calc_ncov(data,varargin)
% Purpose: input continuous or epoched data in TimeSurfer format along with
% event info and return epoched (or re-epoched) data.
%
% Usage: data = ts_calc_ncov(data,'key1',val1, ... )
%
% Example: epoch_data = ts_calc_ncov(epoch_data,'noisewindow',[-.2 0]);
%
% Inputs: 
%     data: timesurfer epoch_data struct
%
% Outputs:
%     data: timesurfer epoch_data struct with noise covariance matrix
%
% Parameters: default : options :  description 
%     noisewindow: [] :: noise covariance matrix calculation window, [begin end] in seconds,
%                        the default is the entire epoch ([-inf inf])
%     ncov_ex_evnts: [] :: vector of event codes that should not be used in calculating the noise covariance matrix
%                          e.g. to exclude events with pre-stimulus activity that is not noise -e.g. a button press
%     concat_ncov_flag: 0 : 0,1 : whether to concatenate noise before calculating noise covarience (1) or 
%                                 calculate noise covarience for each trial and then average (0)
%     max_noise_samples: 5000 :: max # of samples used for noise estimate
%                                Note: especially significant if concat_ncov_flag = 1
%     zparam : [] :: name of the data field to epoch (string; default
%                    depends on the data type)
%     keepbadchans_flag: 0 : 0,1 : whether to include channels previously marked as bad
%     keepbadtrials_flag: 0 : 0,1 : whether to include trials previously marked as bad
%
% Created:  28-Apr-2011, Jason Sherfey
% Last mod: 28-Apr-2011, Jason Sherfey

% consider: getting parameters from data.parms if not specified and parms
% is defined with the parameters specified.


% set default trig_minduration to 5?
parms = mmil_args2parms( varargin, ...
    {'noisewindow',[-inf inf],[],...
     'ncov_ex_evnts',[],[],...
     'max_noise_samples',5000,[],...
     'concat_ncov_flag',false,[false true],...
     'keepbadtrials_flag',false,[false true],...
     'keepbadchans_flag',false,[false true],...
     'zparam',[],[],...
    }, ...
    false );

[datatype,datafield,dataparam] = ts_object_info(data);
if isempty(parms.zparam) || ~ismember(parms.zparam,dataparam)
  parms.zparam = dataparam{1}; % name of data field to epoch
end

% add noise estimate
noise_events          = setdiff([data.(datafield).event_code],parms.ncov_ex_evnts);
noise                 = ts_data_selection(data,'events',noise_events,'toilim',parms.noisewindow,'verbose',0,'keepbadtrials_flag',parms.keepbadtrials_flag,'keepbadchans_flag',parms.keepbadchans_flag);
  % NOTE: this call to ts_data_selection always removes bad trials
noise                 = {noise.(datafield).(parms.zparam)};
data.noise            = [];
data.noise.num_trials = sum([data.(datafield).num_trials]);
if parms.concat_ncov_flag   % concatenate noise from all trials then calc cov
  noise = cellfun(@(x)reshape(x,[size(x,1) size(x,2)*size(x,3)]),noise,'uniformoutput',false);
  noise = cat(2,noise{:});  % chan x samp, contains noise samps from all trials and conditions
  if size(noise,2) > parms.max_noise_samples
    noise = noise(:,1:parms.max_noise_samples);
  end
  data.noise.num_samples = size(noise,2);
  data.noise.covar       = cov(noise');
else                        % average over single trial noise cov matrices
  noise = cat(3,noise{:});
  if size(noise,2) > parms.max_noise_samples
    noise = noise(:,1:parms.max_noise_samples);
  end
  [sx,sy,sz] = size(noise);
  if sz > 1
    noise = squeeze(mat2cell(noise,sx,sy,ones(1,sz)));
  else
    noise = squeeze(mat2cell(noise,sx,sy));
  end
%   noise                  = squeeze(mat2cell(noise,size(noise,1),size(noise,2),ones(1,size(noise,3))));
  data.noise.num_samples = sum(cellfun(@(x)size(x,2),noise));
  data.noise.covar       = cellfun(@(x)cov(x'),noise,'uniformoutput',false);
  data.noise.covar       = cat(3,data.noise.covar{:});
  data.noise.covar       = mean(data.noise.covar,3);
end       
  