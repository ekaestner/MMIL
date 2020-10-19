function [samp,cond] = extract_triggers(trigdata,varargin)
% Purpose: parse trigger channels and return trigger latencies and codes

parms = mmil_args2parms( varargin, ...
    { 'trigthresh',[],[],...
      'trig_minduration',5,[],...
      'trig_minseparation',.01,[],...
    }, ...
    false );
  
if trigdata.num_sensors == 1
  [samp,cond] = extract_amplitude_code(trigdata,parms);
else
  [samp,cond] = extract_binary_code(trigdata,parms);
end

function [samp,cond] = extract_amplitude_code(trigdata,parms)
  % based on ts_read_fif_events()
  X   = trigdata.epochs.data; % 1-D array
  dX  = diff(X);
  change_idx = find(dX~=0);
  change_idx = change_idx(change_idx+parms.trig_minduration<length(dX));
  change_len = diff(change_idx);
  onset_idx  = find(dX >0);
  onset_idx  = onset_idx(onset_idx+parms.trig_minduration<length(dX));
  noise_idx  = change_idx(change_len < parms.trig_minduration);
  noise_onset= [];
  for k = 1:length(noise_idx)
    trig_onset_idx = noise_idx(k);
    trig_code      = dX(trig_onset_idx);
    next_samples   = dX(trig_onset_idx+1:trig_onset_idx+parms.trig_minduration);
    trig_offset_idx= trig_onset_idx + find(next_samples == -trig_code) - 1;
    % noise found (remove from 'differences')
    if ~isempty(trig_offset_idx)
      dX(trig_onset_idx)  = 0;
      dX(trig_offset_idx) = 0;
      noise_onset(end+1)  = trig_onset_idx;
    end
  end
  onset_idx  = find(dX >0);   % refresh onset w/o noise
  trig_onset = onset_idx;% + 1; % adjust to account for the element lost by diff
  cond = dX(onset_idx);
  samp = trig_onset;      
  % try to correct for non-integer event code multiples
  if ~isequal(cond,round(cond))
    tmpcond = cond / min(cond);
    if isequal(tmpcond,round(tmpcond))
      cond = tmpcond;
    end
  end
% end

function [samp,cond] = extract_binary_code(trigdata,parms)
  % assume these channels carry a binary code for epoching
  Fs = trigdata.sfreq;
  % set threshold for trigger detection
  if isempty(parms.trigthresh) || ~isnumeric(parms.trigthresh)
    tmp              = trigdata.epochs.data(:);
    parms.trigthresh = mean(tmp) + std(tmp);
    clear tmp
  end
  % find threshold crossings
  t     = trigdata.epochs.time;
  nsamp = length(t);
  ntrigchan = trigdata.num_sensors;
  [crsind{1:ntrigchan}] = deal([]);
  for k = 1:ntrigchan
    tmpdat  = trigdata.epochs.data(k,:);
    ind     = crossing(tmpdat,[],parms.trigthresh);
    % only keep left edge of pulse
    sel     = (tmpdat(ind+1)-tmpdat(ind-1)) > 0;
    ind     = ind(sel);
%     % only keep crossings that occur > minduration since last crossing, if specified
%     if ~isempty(parms.trig_minduration)
%       sel = [1 find((t(ind(2:end))-t(ind(1:end-1))) >= parms.trig_minduration/Fs)];
%       ind = ind(sel);
%     end
    crsind{k} = ind;
    clear ind sel tmpdat
  end
  bincode = zeros(ntrigchan,nsamp);
  for k   = 1:length(crsind), bincode(k,crsind{k}) = 1; end
  detind  = any(bincode,1); % whether there's a detection in any trigchan
  % make sure detections for the same pulse are grouped before conversion
    delta   = round(parms.trig_minseparation*Fs/2);
    % get indices at which a detection occurs in any trigchan
    ind     = find(detind);
    % find detections with a second detection > delta steps before it
    tmp     = cellfun(@(x)any((ind>x-delta & (ind<x))),num2cell(ind));
    R_ind   = find(tmp);
    % find detections with a second detection < delta steps after it
    tmp     = cellfun(@(x)any((ind<x+delta & (ind>x))),num2cell(ind));
    L_ind   = find(tmp);
    % remove detections that are between two other detections in the same pulse
    Rsel    = ~ismember(R_ind,L_ind);
    Lsel    = ~ismember(L_ind,R_ind);
    R_ind   = R_ind(Rsel); clear Rsel
    L_ind   = L_ind(Lsel); clear Lsel
    % for each pair (L,R), set [ind(L):ind(R)] to 1 in each bincode row w/ a detection
    for k = 1:ntrigchan
      sel = cellfun(@(x,y)any(bincode(k,ind(x):ind(y))),num2cell(L_ind),num2cell(R_ind));
      sel = cellfun(@(x,y)ind(x):ind(y),num2cell(L_ind(sel)),num2cell(R_ind(sel)),'uniformoutput',false);
      sel = unique([sel{:}]);
      bincode(k,sel) = 1;
    end
  detind  = any(bincode,1); % whether there's a detection in any trigchan    
  % only keep the first pulse in contiguous detections
  samp    = find(detind(1:end-1)==0 & detind(2:end)==1) + 1;
  bincode = flipud(bincode(:,samp));
  bincode = mat2cell(bincode,ntrigchan,ones(1,length(samp)));
  % convert binary to decimal
  cond    = cellfun(@(x)polyval(x,2),bincode);
% end