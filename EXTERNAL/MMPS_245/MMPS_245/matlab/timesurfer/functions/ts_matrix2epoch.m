function epoch_data = ts_matrix2epoch(y,varargin)
% y - data matrix (1 condition)
% Optional:
%   time    - time vector in seconds
%   sfreq   - sampling frequency in Hz
%   sens    - TimeSurfer sensor_info structure
%   hdr     - TimeSurfer structure without datafield
%   tstart  - time of first point in the time vector
%     note: (sfreq & tstart) can determine the time vector
%   continuous - whether the data should be returned as continuous
%   mri2head
%   noise_covar
%
% Last modified on 25-Oct-2010 by JSS

parms = mmil_args2parms(varargin,...
						{'time',[],[],...
             'sfreq',100,[],...
						 'sens',[],[],...
             'hdr' ,[],[],...
             'tstart',0,[],...
             'continuous',0,[],...
             'device2head',eye(4),[],...
             'mri2head',[],[],...
             'noise_covar',[],[],...
             'datafield','epochs',[],...
             'trial_info',[],[],...
						},false);

t = parms.time;

% if [y] = 1 x samples x 1
if ndims(y) == 2 && any(size(y)==1)
  if size(y,1) > size(y,2)
    y = y';
  end
  nsen = 1;  
  ntrl = 1;
end

% if continuous
if ndims(y) == 2 && parms.continuous
  nsen = size(y,1);
  ntrl = 1;
elseif ndims(y) == 2 && ~any(size(y)==1)
  % assume 1 channel and multiple trials if not continuous, ndims(y) = 2 && ~any(size(y)==1)
  yy(1,:,:) = y;  clear y;
  y         = yy; clear yy;
end

% if [y] = 1 x samples x trials
if ndims(y) == 3 && any(size(y)==1)
  nsen = 1;
  ntrl = size(y,3);
end

% if [y] = channels x samples x trials
if ndims(y) == 3 && ~any(size(y)==1)
  nsen = size(y,1);
  ntrl = size(y,3);
end

% default time vector is shifted with pre-stimulus period = 1/3 data length
if isempty(t)
  N  = size(y,2);
  dt = 1 / parms.sfreq;
  t  = [0:N-1]*dt + parms.tstart;
%   t  = [1:N]*dt - N*dt/3;
else
  parms.sfreq = 1 / (t(2) - t(1));
end

nsmp = length(t);  

% set up epoch_data structure
% default channel has intracranial eeg parameters
epoch_data.num_sensors                = nsen;
if ~isempty(parms.sens) && isstruct(parms.sens) && length(parms.sens)==nsen
  epoch_data.sensor_info = parms.sens;
else
  for i = 1:nsen  
    epoch_data.sensor_info(i).label       = sprintf('chan%g',i);
    epoch_data.sensor_info(i).typestring  = 'eeg';
    epoch_data.sensor_info(i).type        = 1;
    epoch_data.sensor_info(i).kind        = 2;
    epoch_data.sensor_info(i).badchan     = 0;
    epoch_data.sensor_info(i).lognum      = i;
    epoch_data.sensor_info(i).loc         = eye(4);
  end
end
epoch_data.coor_trans.device2head     = parms.device2head;
epoch_data.coor_trans.mri2head        = parms.mri2head;
epoch_data.sfreq                      = parms.sfreq;
epoch_data.noise.num_trials           = ntrl;
epoch_data.noise.num_samples          = nsmp;
% if isempty(parms.noise_covar)
%   epoch_data.noise.covar              = cov(rand(nsmp*ntrl,nsen));
% else
  epoch_data.noise.covar              = parms.noise_covar;
% end
epoch_data.(parms.datafield).event_code          = 1;
epoch_data.(parms.datafield).num_trials          = ntrl;
epoch_data.(parms.datafield).num_rejects.mag     = 0;
epoch_data.(parms.datafield).num_rejects.grad    = 0;
epoch_data.(parms.datafield).num_rejects.eeg     = 0;
epoch_data.(parms.datafield).num_rejects.eog     = 0;
epoch_data.(parms.datafield).num_rejects.manual  = 0;
epoch_data.(parms.datafield).num_rejects.skip    = 0;
epoch_data.(parms.datafield).time                = t;
epoch_data.(parms.datafield).data                = y;
% epoch_data.epochs.label               = 'cond1';
if ~isempty(parms.trial_info)
  epoch_data.(parms.datafield).trial_info        = parms.trial_info;
else
  % TODO: add default trial_info
end
