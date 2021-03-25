function [data_comb] = ts_combine_data(varargin)
% [data_comb] = ts_comb_data(data_arr);
%
% Purpose: 
%   to concatenate the metadata and data of multiple
%   {average,epoch,timefreq} objects into a single object.
%   This combination keeps the fidelity of the original data.
%
% Usage:
%  [data_comb] = ts_comb_data([data_1,data_2,...]);
%
% Required input:
%  data_arr - array of data structures (average, epoch, timefreq, or other)
%                       (output from ts_avg_fif_data.m)
%
% Output:
%   data_comb - structure containing the following fields:
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
%      {average, epoch, or timefreq} (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples x {average, epoch, timefreq})
%        [frequencies (vector of doubles) (in hz)] (timefreq only)
%      noise          (struct)
%        num_trials   (int)
%        num_samples  (int)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%  created:       07/23/07   by Ben Cipollini
%  last modified: 04/16/08   by Rajan Patel 
 
% 04/16/08 - allowed num_trials to be an array = to number of sensors for
%            time frequency data (important after rejection)
% 04/03/08 - fixed combining across multiple timefreq data structures
% 09/15/08 - fixed combining across multiple timefreq data structures for power or data
% 					 matrices that vary in size between events
% 10/13/10 - added trial_info to the combine_common subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin, 1))
  return;
end;

try
  data_arr  = [varargin{:}];
catch
end;

if (~exist('data_arr','var') || ~isstruct(data_arr) | length(data_arr)<2)
  error('%s: Input should be an array of [average_data, epoch_data, or timefreq_data] structures',mfilename);
end;

% % for backward/forward compatibility
% if isfield(data_arr,'opt') || isfield(data_arr,'parms')
%   p = {};
%   try p = {data_arr.opt,p{:}};   end
%   try p = {data_arr.parms,p{:}}; end 
%   for i = 1:length(p), parm_orig(i) = p{i};   end
% end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine shared props
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%copy original structure & data
data_comb = rmfield(data_arr(1), intersect({'averages', 'epochs', 'timefreq'}, fields(data_arr(1))));

%combine noise
allnoise = [data_arr.noise];
data_comb.noise.num_trials  = sum([allnoise.num_trials]);
data_comb.noise.num_samples = sum([allnoise.num_samples]);
data_comb.noise.covar       = zeros(size(data_arr(1).noise.covar));

% weighted sum of noise covariance
for i=1:length(data_arr)
  weight                    = data_arr(i).noise.num_trials / data_comb.noise.num_trials;
  data_comb.noise.covar     = data_comb.noise.covar ...
                              + weight*data_arr(i).noise.covar;
end;

% added 09/15/08 by Jason Sherfey
% - to make "badchan" an array of badchan flags for all conditions if they differ
%allsensorinfo = [data_arr.sensor_info];
%[diff_badchan_flags,badchan_flags] = compare_badchan_flags( allsensorinfo );	
%if diff_badchan_flags, 
%	[data_comb.sensor_info.badchan] = deal(badchan_flags{:});
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isfield(data_arr, 'epochs'))

  %find set of event_codes
  allepochs         = [data_arr.epochs];
  
  % find a unique set of events, and order them
  epoch_event_codes = unique([allepochs.event_code]);

  for i=1:length(epoch_event_codes)
    
    % Find epochs across all objects with this event code
    epochs_to_include = find([allepochs.event_code]==epoch_event_codes(i));
    allepochs_cond    = allepochs(epochs_to_include);
    tot_num_trials    = sum([allepochs_cond.num_trials]);
    
    % Combine common data (rejects, num_trials, etc)
    cond_epoch        = combine_common( allepochs_cond );
    
    % concatenate epoch data
    cond_epoch.data   = cat(3, allepochs_cond.data);
  
    data_comb.epochs(i) = cond_epoch;
    
    clear('cond_epoch');
  end;
  
  clear('allepochs');
  clear('epoch_event_codes');
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isfield(data_arr, 'averages'))

  %find set of event_codes
  allaverages         = [data_arr.averages];
  average_event_codes = unique([allaverages.event_code]);

  for i=1:length(average_event_codes)
    
    % Find averages across all objects with this event code
    averages_to_include = find([allaverages.event_code]==average_event_codes(i));
    allaverages_cond    = allaverages(averages_to_include);
    tot_num_trials      = sum([allaverages_cond.num_trials]);
    
    % Combine common data (rejects, num_trials, etc)
    cond_average        = combine_common( allaverages_cond );
    
    % concatenate mean & stdev 
    cond_average.data   = zeros(size(allaverages(1).data));
    cond_average.stdev  = zeros(size(allaverages(1).stdev));
    
    for j=1:length(allaverages_cond)
      num_trials = allaverages_cond(j).num_trials;
      weight     = num_trials / tot_num_trials;

      % mean
      cond_average.data  = cond_average.data ...
                           + weight*allaverages_cond(j).data;
    
      % calculate sums of squares from stdev and average
      cond_average.stdev = cond_average.stdev ...
                           + (num_trials-1)*(allaverages_cond(j).stdev.^2) ...
                           + (num_trials)*(allaverages_cond(j).data.^2);
    end;
    
    % calculate standard deviation
    if (tot_num_trials>1)
      cond_average.stdev = ...
        sqrt((cond_average.stdev - tot_num_trials*cond_average.data)/(tot_num_trials-1));
    end;

  
    data_comb.averages(i)             = cond_average;
    
    clear('cond_average');
  end;
  
  clear('allepochs');
  clear('epoch_event_codes');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine timefreq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (isfield(data_arr, 'timefreq'))

  %find set of event_codes
  alltimefreq           = [data_arr.timefreq];
  timefreq_event_codes  = unique([alltimefreq.event_code]);

  for i=1:length(timefreq_event_codes)
     
    % Find timefreq across all objects with this event code
    alltimefreq_cond  = alltimefreq(find([alltimefreq.event_code]==timefreq_event_codes(i)));
    
    
    % Make sure all data have the same sets of time & frequency;
    % otherwise, we'd require extra metadata and an interpolation
    % algorithm to proceed.
    %
    % NOTE that this check ONLY needs to be done on
    % combinations on the same condition; e.g. different
    % conditions are allowed to have different frequencies or timepoints!
    for j=2:length(alltimefreq_cond)
      
      % Failed frequency consistnecy check
      if (~isempty(setdiff(alltimefreq_cond(1).frequencies, alltimefreq_cond(j).frequencies)))
        error('%s: Cannot combine timefreq objects with different frequencies.', mfilename);
        
      % Failed time consistency check
      elseif isfield(alltimefreq_cond,'power') && ((size(alltimefreq_cond(1).power,2) ~= size(alltimefreq_cond(j).power,2)))
        error('%s: Cannot combine timefreq objects with different #s of timepoints.', mfilename);
      elseif isfield(alltimefreq_cond,'cmplx') && ((size(alltimefreq_cond(1).cmplx,2) ~= size(alltimefreq_cond(j).cmplx,2)))
        error('%s: Cannot combine timefreq objects with different #s of timepoints.', mfilename);
      end;
    end;    
    % Combine common data (rejects, num_trials, etc)
    cond_timefreq                    = combine_common( alltimefreq_cond );
    cond_timefreq.frequencies        = alltimefreq_cond(1).frequencies;
    
    % Weighted sum of timefreq data,
    % which are complex vectors
    try cond_timefreq.data  = zeros(size(alltimefreq_cond(1).data)); end;
    try cond_timefreq.power = zeros(size(alltimefreq_cond(1).power)); end;
    try cond_timefreq.cmplx = zeros(size(alltimefreq_cond(1).cmplx)); end
%     try cond_timefreq.cmplx = zeros(size(alltimefreq_cond(1).cmplx)); end;    
    for j=1:length(alltimefreq_cond)
        % is there a 'data' field?
        if isfield(alltimefreq_cond(j),'data')
            datafield = 'data';
        elseif isfield(alltimefreq_cond(j),'power')
            datafield = 'power';
        elseif isfield(alltimefreq_cond(j),'cmplx')
            datafield = 'cmplx';            
        end
        
        weight = alltimefreq_cond(j).num_trials./cond_timefreq.num_trials;
		if ndims(alltimefreq_cond(j).(datafield)) == 3
      % no trials
	      if length(weight) == 1
	          weight = repmat(permute(weight,[2 1]),[size(alltimefreq_cond(j).(datafield),1)...
	                                                 size(alltimefreq_cond(j).(datafield),2)...
	                                                 size(alltimefreq_cond(j).(datafield),3)]);
	      else
	          weight = repmat(permute(weight,[2 1]),[1 ...
	                                                 size(alltimefreq_cond(j).(datafield),2)...
	                                                 size(alltimefreq_cond(j).(datafield),3)]);              
	      end
			elseif ndims(alltimefreq_cond(j).(datafield)) == 4
        % trials
	      if length(weight) == 1
	          weight = repmat(permute(weight,[2 1]),[size(alltimefreq_cond(j).(datafield),1)...
	                                                 size(alltimefreq_cond(j).(datafield),2)...
	                                                 size(alltimefreq_cond(j).(datafield),3)...
													 size(alltimefreq_cond(j).(datafield),4)]);
	      else
	          weight = repmat(permute(weight,[2 1]),[1 ...
	                                                 size(alltimefreq_cond(j).(datafield),2)...
	                                                 size(alltimefreq_cond(j).(datafield),3)...
													 size(alltimefreq_cond(j).(datafield),4)]);              
	      end			
		else
				error('%s: Time-frequency data must be 3-D or 4-D.',mfilename);
		end

      try cond_timefreq.data  = cond_timefreq.data ...
                                    + weight.*alltimefreq_cond(j).data; end;
      try cond_timefreq.power = cond_timefreq.power ...
                                    + weight.*alltimefreq_cond(j).power; end;
      try cond_timefreq.cmplx = cond_timefreq.cmplx ...
                                    + weight.*alltimefreq_cond(j).cmplx; end;    
                                
    end;
    data_comb.timefreq(i)  = cond_timefreq;
  end;
	
  clear('alltimefreq');
  clear('timefreq_event_codes');
end;

% if exist('parm_orig','var')
%   data_comb.parms = parm_orig;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%
  function obj   = combine_common( obj_arr )
  %
  %
    obj.event_code    = obj_arr(1).event_code;
    obj.num_trials    = obj_arr(1).num_trials;
    for i = 2:length(obj_arr)
      obj.num_trials  = obj.num_trials+obj_arr(i).num_trials;  
    end
    obj.num_rejects   = combine_rejects( [obj_arr.num_rejects] );
    obj.time          = obj_arr(1).time;

    if isfield(obj_arr,'trial_info')
      obj.trial_info  = combine_trial_info( [obj_arr.trial_info] );
    end
   
  %%%%%%%%%%%%%%%
  function rejects  = combine_rejects( num_rejects_arr )
  %
  %
    rejects.mag    = sum([num_rejects_arr.mag]);
    rejects.grad   = sum([num_rejects_arr.grad]);
    rejects.eeg    = sum([num_rejects_arr.eeg]);
    rejects.eog    = sum([num_rejects_arr.eog]);
    rejects.manual = sum([num_rejects_arr.manual]);
    rejects.skip   = sum([num_rejects_arr.skip]);

    %%%%%%%%%%%%%%%
	function [diff_badchan_flags,badchan_flags] = compare_badchan_flags( sens_arr )
	%
	%
		tf_labels = unique({sens_arr.label});		
		diff_badchan_flags = 0;
		for i=1:length(tf_labels)
			sens_idx 							= find(strcmp(tf_labels{i},{sens_arr.label}));
			badchan_flags{i} 			= [sens_arr(sens_idx).badchan];
			if length(unique(badchan_flags{i})) > 1, 
				diff_badchan_flags 	= 1; 
			end  	
    end  
    
    %%%%%%%%%%%%%%%
    function trl = combine_trial_info( trial_arr )
    trl   = trial_arr(1);
    for i = 2:length(trial_arr)
      trl.datafile      = {trl.datafile{:} trial_arr(i).datafile{:}};
      trl.events_fnames = {trl.events_fnames{:} trial_arr(i).events_fnames{:}};
      trl.latency       = [trl.latency trial_arr(i).latency];
      trl.event_code    = [trl.event_code trial_arr(i).event_code];
      trl.duration      = [trl.duration trial_arr(i).duration];
      trl.badtrial      = [trl.badtrial trial_arr(i).badtrial];
      trl.number        = [trl.number trial_arr(i).number];
      if isfield(trl,'epoch_num')
        trl.epoch_num   = [trl.epoch_num trial_arr(i).epoch_num];
      end
    end
      
%       files = trial_arr(i).datafile;
%       for j = 1:length(files)
%         if ismember(files{j},trl.datafile)
%           k = strmatch(files{j},trl.datafile);
%           trl.latency{k}      = [trl.latency{k} trial_arr(i).latency{j}];
%           trl.duration{k}     = [trl.duration{k} trial_arr(i).duration{j}];
%           trl.badtrial{k}     = [trl.badtrial{k} trial_arr(i).badtrial{j}];
%           trl.trial_number{k} = [trl.trial_number{k} trial_arr(i).trial_number{j}];
%         else
%           k = length(trl.datafile) + 1;
%           trl.datafile{k}       = trial_arr(i).datafile{j};
%           trl.events_fnames{k}  = trial_arr(i).events_fnames{j};
%           trl.latency{k}        = trial_arr(i).latency{j};
%           trl.duration{k}       = trial_arr(i).duration{j};
%           trl.badtrial{k}       = trial_arr(i).badtrial{j};
%           trl.trial_number{k}   = trial_arr(i).trial_number{j};          
%         end
%       end
%     end