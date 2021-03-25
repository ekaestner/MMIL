function [avg_data_comb] = ts_comb_avg_data(avg_data_cellarr);
% [avg_data_comb] = ts_comb_avg_data(avg_data_cellarr);
%
% Usage:
%  [avg_data_comb] = ts_comb_avg_data({avg_data_1,avg_data_2,...});
%
% Required input:
%  avg_data_cellarr - cell array of avg_data structures
%                     averaged data structure (output from avg_fif_data.m)
%
% Output:
%   avg_data - structure containing the following fields:
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
%      noise          (struct)
%        num_trials   (int)
%        num_samples  (int)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%  created:       04/26/06   by Don Hagler
%  last modified: 01/04/07   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help(mfilename);
   return;
end

if ~iscell(avg_data_cellarr) | length(avg_data_cellarr)<2
  fprintf('%s: input should be a cell array of avg_data structures\n',mfilename);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keyboard

% make nconds equal to number of unique event codes across all avg_data's
unique_codes = [];
avg_data_raw = [];
k = 1;
num_avgs = length(avg_data_cellarr);
for a=1:num_avgs
  if ~isempty(avg_data_cellarr{a})
    codes = unique(cell2mat({avg_data_cellarr{a}.averages.event_code}));
    unique_codes = unique([unique_codes,codes]);
    avg_data_raw{k} = avg_data_cellarr{a};
    k = k+1;
  end;
end;

% initialize output
avg_data_comb = avg_data_raw{1};
data_size = [];
time = [];
nconds = length(unique_codes);
for j=1:nconds
  if length(avg_data_comb.averages)<j
    avg_data_comb.averages(j).data = [];
  elseif ~isempty(avg_data_comb.averages(j).data)
    if isempty(data_size)
      data_size = size(avg_data_comb.averages(j).data);
      time = avg_data_comb.averages(j).time;
    elseif any(size(avg_data_comb.averages(j).data)~=data_size)
      fprintf('%s: ERROR: data size is not consistent across files\n',mfilename);
      fprintf('%s:        Reprocess with raw files with different data sizes (channels x time)\n',mfilename);
      fprintf('%s:        separated into different directories\n',mfilename);
      return;
    end;
  end;
end;
for j=1:nconds
  avg_data_comb.averages(j).event_code = unique_codes(j);
  avg_data_comb.averages(j).num_trials = 0;
  avg_data_comb.averages(j).stdev = zeros(data_size);
  avg_data_comb.averages(j).data = zeros(data_size);
  avg_data_comb.averages(j).time = time;
  avg_data_comb.averages(j).num_rejects.mag = 0;
  avg_data_comb.averages(j).num_rejects.grad = 0;
  avg_data_comb.averages(j).num_rejects.eeg = 0;
  avg_data_comb.averages(j).num_rejects.eog = 0;
  avg_data_comb.averages(j).num_rejects.manual = 0;
  avg_data_comb.averages(j).num_rejects.skip = 0;
  avg_data_comb.averages(j).num_rejects.skip = 0;
end;
avg_data_comb.noise.num_trials = 0;
avg_data_comb.noise.num_samples = 0;
avg_data_comb.noise.covar = zeros(size(avg_data_comb.noise.covar));

num_avgs = length(avg_data_raw);
for a=1:num_avgs
  avg_data = avg_data_raw{a};
  event_codes = unique(cell2mat({avg_data.averages.event_code}));
  for j=1:nconds
    % is this condition found in avg_data?
    event_code = avg_data_comb.averages(j).event_code;
    k = find(event_code==event_codes);
    if isempty(k), continue; end;

    % num_trials
    n = avg_data.averages(k).num_trials;
    avg_data_comb.averages(j).num_trials = ...
      avg_data_comb.averages(j).num_trials + n;

    % num_rejects
    avg_data_comb.averages(j).num_rejects.mag = ...
      avg_data_comb.averages(j).num_rejects.mag + ...
      avg_data.averages(k).num_rejects.mag;
    avg_data_comb.averages(j).num_rejects.grad = ...
      avg_data_comb.averages(j).num_rejects.grad + ...
      avg_data.averages(k).num_rejects.grad;
    avg_data_comb.averages(j).num_rejects.eeg = ...
      avg_data_comb.averages(j).num_rejects.eeg + ...
      avg_data.averages(k).num_rejects.eeg;
    avg_data_comb.averages(j).num_rejects.eog = ...
      avg_data_comb.averages(j).num_rejects.eog + ...
      avg_data.averages(k).num_rejects.eog;
    avg_data_comb.averages(j).num_rejects.manual = ...
      avg_data_comb.averages(j).num_rejects.manual + ...
      avg_data.averages(k).num_rejects.manual;
    avg_data_comb.averages(j).num_rejects.skip = ...
      avg_data_comb.averages(j).num_rejects.skip + ...
      avg_data.averages(k).num_rejects.skip;

    % calculate sums of squares from stdev and average
    if ~isempty(avg_data_comb.averages(j).stdev) & ~isempty(avg_data_comb.averages(j).data)
      avg_data_comb.averages(j).stdev = ...
        avg_data_comb.averages(j).stdev + ...
        (n-1)*avg_data.averages(k).stdev.*avg_data.averages(k).stdev + ...
        n*avg_data.averages(k).data.*avg_data.averages(k).data;
    end;

    % calculate sums from average
    avg_data_comb.averages(j).data = ...
      avg_data_comb.averages(j).data + ...
      n*avg_data.averages(k).data;
  end;
  n = avg_data.noise.num_trials;
  avg_data_comb.noise.num_trials = avg_data_comb.noise.num_trials + n;
  avg_data_comb.noise.covar = avg_data_comb.noise.covar + n*avg_data.noise.covar;
  n = avg_data.noise.num_samples;
  avg_data_comb.noise.num_samples = avg_data_comb.noise.num_samples + n;
end;

% finish combining averages by dividing by total num_trials
nconds=length(avg_data_comb.averages);
for j=1:nconds
  N = avg_data_comb.averages(j).num_trials;
  if (N>0)
    % calculate average
    avg_data_comb.averages(j).data = ...
      avg_data_comb.averages(j).data / N;
  end;
  if (N>1)
    % calculate standard deviation
    avg_data_comb.averages(j).stdev = ...
      sqrt((avg_data_comb.averages(j).stdev - N*avg_data_comb.averages(j).data)/(N-1));
  end;
end;

% finish combining noise by dividing by total num_trials
avg_data_comb.noise.covar = avg_data_comb.noise.covar/avg_data_comb.noise.num_trials;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

