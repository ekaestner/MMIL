function [epoch_data_comb] = ts_comb_epoch_data(epoch_data_cellarr);
% [epoch_data_comb] = ts_comb_epoch_data(epoch_data_cellarr);
%
% Usage:
%  [epoch_data_comb] = ts_comb_epoch_data({epoch_data_1,epoch_data_2,...});
%
% Required input:
%  epoch_data_cellarr - cell array of epoch_data structures
%                       (output from avg_fif_data.m)
%
% Output:
%   epoch_data - structure containing the following fields:
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
%      epochs         (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples x num_trials)
%      noise          (struct)
%        num_trials   (int)
%        num_samples  (int)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%  created:       07/02/06   by Don Hagler
%  last modified: 01/04/07   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help(mfilename);
   return;
end

if ~iscell(epoch_data_cellarr) | length(epoch_data_cellarr)<2
  fprintf('%s: input should be a cell array of epoch_data structures\n',mfilename);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_inputs = length(epoch_data_cellarr);

% initialize output
epoch_data_comb = epoch_data_cellarr{1};
nconds = length(epoch_data_comb.epochs);
for j=1:nconds
  epoch_data_comb.epochs(j).num_trials = 0;
  epoch_data_comb.epochs(j).data = [];
  epoch_data_comb.epochs(j).num_rejects.mag = 0;
  epoch_data_comb.epochs(j).num_rejects.grad = 0;
  epoch_data_comb.epochs(j).num_rejects.eeg = 0;
  epoch_data_comb.epochs(j).num_rejects.eog = 0;
  epoch_data_comb.epochs(j).num_rejects.manual = 0;
  epoch_data_comb.epochs(j).num_rejects.skip = 0;
end;
epoch_data_comb.noise.num_trials = 0;
epoch_data_comb.noise.num_samples = 0;
epoch_data_comb.noise.covar = zeros(size(epoch_data_comb.noise.covar));

for a=1:num_inputs
  epoch_data = epoch_data_cellarr{a};
  nconds=min(length(epoch_data_comb.epochs),length(epoch_data.epochs));
  for j=1:nconds
    % num_rejects
    epoch_data_comb.epochs(j).num_rejects.mag = ...
      epoch_data_comb.epochs(j).num_rejects.mag + ...
      epoch_data.epochs(j).num_rejects.mag;
    epoch_data_comb.epochs(j).num_rejects.grad = ...
      epoch_data_comb.epochs(j).num_rejects.grad + ...
      epoch_data.epochs(j).num_rejects.grad;
    epoch_data_comb.epochs(j).num_rejects.eeg = ...
      epoch_data_comb.epochs(j).num_rejects.eeg + ...
      epoch_data.epochs(j).num_rejects.eeg;
    epoch_data_comb.epochs(j).num_rejects.eog = ...
      epoch_data_comb.epochs(j).num_rejects.eog + ...
      epoch_data.epochs(j).num_rejects.eog;
    epoch_data_comb.epochs(j).num_rejects.manual = ...
      epoch_data_comb.epochs(j).num_rejects.manual + ...
      epoch_data.epochs(j).num_rejects.manual;
    epoch_data_comb.epochs(j).num_rejects.skip = ...
      epoch_data_comb.epochs(j).num_rejects.skip + ...
      epoch_data.epochs(j).num_rejects.skip;

    % concatenate epochs
    n0 = epoch_data_comb.epochs(j).num_trials;
    n1 = epoch_data.epochs(j).num_trials;
    epoch_data_comb.epochs(j).data(:,:,n0+1:n0+n1) = ...
      epoch_data.epochs(j).data;
    epoch_data_comb.epochs(j).num_trials=n0+n1;
  end;
  n = epoch_data.noise.num_trials;
  epoch_data_comb.noise.num_trials = epoch_data_comb.noise.num_trials + n;
  epoch_data_comb.noise.covar = epoch_data_comb.noise.covar + n*epoch_data.noise.covar;
  n = epoch_data.noise.num_samples;
  epoch_data_comb.noise.num_samples = epoch_data_comb.noise.num_samples + n;
end;

% finish combining noise by dividing by total num_trials
epoch_data_comb.noise.covar = epoch_data_comb.noise.covar/epoch_data_comb.noise.num_trials;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

