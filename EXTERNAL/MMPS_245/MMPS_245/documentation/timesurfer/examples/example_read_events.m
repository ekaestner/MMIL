test_data = [repmat(0,[1 20]) repmat(2,[1 40]) repmat(10,[1 20])   repmat(8,[1 20]) repmat(0, [1 40]) ... %overlapping triggers
             repmat(3, [1 3]) repmat(0, [1 40]) ... %onset noise
             repmat(3, [1 20]) repmat(0, [1 3]) repmat(3, [1 20]) repmat(0, [1 40]) ... %offset noise
             repmat(20, [1 20])];
           
%figure; plot(test_data); title('test data');

trigger_diff = diff(test_data);

% check for "noise"--trigger changes that
% last less than some # of samples
change_idx = find(trigger_diff~=0);
onset_idx = find(trigger_diff >0);
change_lengths = diff(change_idx);

% Noise is where we see:
% - onset of some value
% - offset of SAME VALUE less than X samples later
% Otherwise, we just let it through; too complex
% to separate noise from fast async trigger value changes
noise_candidate_idx = change_idx(find(change_lengths < 5));
noise_onset_idx = [];
for i=1:length(noise_candidate_idx)
  trig_onset_idx = noise_candidate_idx(i);
  trig_code      = trigger_diff(trig_onset_idx);
  next_samples   = trigger_diff(trig_onset_idx:trig_onset_idx+5);
  trig_offset_idx= trig_onset_idx + find(next_samples == -trig_code) - 1;
  
  if (~isempty(trig_offset_idx) )
    trigger_diff(trig_onset_idx) = 0;
    trigger_diff(trig_offset_idx) = 0;
    noise_onset_idx(end+1) = trig_onset_idx;
  end;
end;

% Refresh onset index after noise removal
onset_idx = find(trigger_diff >0);
    
figure; hold on;
plot(test_data, 'k'); 
title('test data');
plot(onset_idx, test_data(onset_idx), 'b*');
plot(noise_onset_idx, test_data(noise_onset_idx), 'r*');
legend({'test data', 'trigger onsets', 'noise onsets'}, 'Location', 'NorthEast');
%offset_idx = find(trigger_diffs<0);

