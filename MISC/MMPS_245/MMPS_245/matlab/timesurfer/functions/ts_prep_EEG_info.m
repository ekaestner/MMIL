function EEG_info = ts_prep_EEG_info(avg_data)
%function EEG_info = ts_prep_EEG_info(avg_data)
%
% Required Input:
%   avg_data: structure containing the following fields:
%      num_sensors    (int)
%      sfreq          (double)
%      sensor_info    (struct)
%        type         (string)
%        label        (string)
%        loc          (4x4 matrix of doubles)
%        badchan      (int: 0 or 1)
%      coor_trans     (struct)
%        device2head  (4x4 matrix of doubles) (or empty)
%        mri2head     (4x4 matrix of doubles) (or empty)
%      averages       (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples)
%        stdev        (matrix of doubles) (num_sensors x num samples)
%      noise          (struct)
%        num_trials   (int)
%        scale_fact   (double)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
% Note: minimally necessary information is "sensor_info"
%     field of avg_data structure, which is returned by ts_avg_fif_data
%
% Output:
%   EEG sensor information
%   Note: coordinates are returned in "head" space
%
%   EEG_info.sensor: 1 by nEEG structure. 
%       sensor(i).n the number of integration points in each channel
%       sensor(i).wei are the weight matrix for the all integration points in channel i; 
%       sensor(i).loc location of integration points;
%   EEG_info.intpnt: location of all integration points concatenated 
%
% Created:  07/31/06 by Don Hagler
% Last Mod: 10/11/13 by Don Hagler
%

EEG_info=[];
sensor=[];
intpnt=[];
labels={};

EEG_chans = find(strcmp('eeg',lower({avg_data.sensor_info.typestring})));
nEEG = length(EEG_chans);

coords = zeros(nEEG,3);
for c=1:nEEG
  k = EEG_chans(c);
  Trans=avg_data.sensor_info(k).loc;
  sensor(c).n=2;
  sensor(c).wei=diag([1 -1]); % subtract reference electrode
  sensor(c).loc=[Trans(1:3,4)';Trans(1:3,1)']; % 1st column of trans is ref elec
  intpnt=[intpnt;sensor(c).loc];
  labels{c} = avg_data.sensor_info(k).label;
end;

EEG_info.sensor = sensor;
EEG_info.intpnt = intpnt;
EEG_info.labels = labels;

return
