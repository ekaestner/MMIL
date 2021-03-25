function avg_data = ts_matrix2avg(y,varargin)
% y - data matrix (1 condition)

parms = mmil_args2parms(varargin,...
						{'time',[],[],...
             'sfreq',1000,[],...
						 'sens',[],[],...
             'hdr' ,[],[],...
             'tstart',0,[],...
						},false);

t = parms.time;

% if [y] = 1 x samples
if ndims(y) == 2 && any(size(y)==1)
  if size(y,1) > size(y,2)
    y = y';
  end
end

% if [y] = channels x samples x trials
if ndims(y) == 3 && ~any(size(y)==1)
  error('%s: use ts_matrix2epoch() to convert to epoch_data\n',mfilename);
end

% default time vector is shifted with pre-stimulus period = 1/3 data length
if isempty(t)
  N  = size(y,2);
  dt = 1 / parms.sfreq;
  t  = [0:N-1]*dt + parms.tstart;
%   t  = [1:N]*dt - N*dt/3;
else
  parms.sfreq = round(1 / (t(2) - t(1)));
end

nsmp = length(t);  
nsen = size(y,1);
ntrl = 1;

% set up avg_data structure
% default channel has intracranial eeg parameters
avg_data.num_sensors                = nsen;
if ~isempty(parms.sens) && isstruct(parms.sens) && length(parms.sens)==nsen
  avg_data.sensor_info = parms.sens;
else
  for i = 1:nsen  
    avg_data.sensor_info(i).typestring  = 'eeg';
    avg_data.sensor_info(i).label       = sprintf('chan%g',i);
    avg_data.sensor_info(i).loc         = eye(4);
    avg_data.sensor_info(i).badchan     = 0;
    avg_data.sensor_info(i).type        = 1;
    avg_data.sensor_info(i).kind        = 2;
    avg_data.sensor_info(i).lognum      = i;
  end
end
avg_data.coor_trans.device2head     = eye(4);
avg_data.coor_trans.mri2head        = [];
avg_data.sfreq                      = parms.sfreq;
avg_data.noise.num_trials           = ntrl;
avg_data.noise.num_samples          = nsmp;
avg_data.noise.covar                = cov(rand(nsmp*ntrl,nsen));
avg_data.averages.event_code          = 1;
avg_data.averages.num_trials          = ntrl;
avg_data.averages.num_rejects.mag     = 0;
avg_data.averages.num_rejects.grad    = 0;
avg_data.averages.num_rejects.eeg     = 0;
avg_data.averages.num_rejects.eog     = 0;
avg_data.averages.num_rejects.manual  = 0;
avg_data.averages.num_rejects.skip    = 0;
avg_data.averages.time                = t;
avg_data.averages.data                = y;
avg_data.averages.std_dev             = zeros(size(y));
% avg_data.averages.label               = 'cond1';

