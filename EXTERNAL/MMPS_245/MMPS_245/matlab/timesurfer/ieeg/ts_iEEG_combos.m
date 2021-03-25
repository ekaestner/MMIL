function ts_iEEG_combos (varargin)

% Usage: ts_iEEG_combos ('option',value,...
%
% This is a wrapper for the ts_comb_conditions functions for use in the
% iEEG stream.  
%
% For avg_data:
%
% It serves to produce weighted averages and subtraction
% conditions across multiple conditions for avg_data. Multiple files names
% may be specified if you want to create a avg_data
% file that contains the conditions from all the listed files.  If there
% are multiple files and no list of conditions or events given then a new
% file will be saved that combines the listed files.
%
% For epoch_data:
%
% The file name will be used to find each event/condtion file to create a 
% condition that has the trials from all the the files.  
% A new file will be created containing only that data.  
% You must provide the 'eventcode' option so a proper file name
% can be contructed.  
%
% Required input: 
%
%   f_path   - path to the files
%   avg_name - the name of the file containing avg_data
%     
%     AND/OR
%
%   epoch_name - the template name of the epoch_data file(s)
%   timefreq_name - the template name of the timefreq_data file(s)
%   neweventcodes - a list of event codes for the new conditions
%   combinations - a list of combinations to produce.  Each combination is
%                supplied as a string.  The numbers within the string are
%                treated either as event codes (if using the reference
%                'events') or as condition numbers (if using the reference
%                'conditions').  Valid examples are:
%    To combine events 1,2 & 3:   '1+2+3'
%    To subtract event 10 from 2: '2-10' (Note: Only two events at a time
%                                               is valid for subtractions)
%   
%    Only simple operations are available.  You cannot for instance specify
%    a combination that includes both combining and then subtracting.  There is
%    a workaround.  For instance...
%
%    If the data contained 10 events with event codes 1-10 and you wanted to
%    do the following: (1+2+3+4+5)-(6+7+8+9+10) the settings would be -
%      events       : {'1+2+3+4+5','6+7+8+9+10','11-12'}
%      neweventcodes: 11 12 13
%    So that you create the two averages across multiple events first and
%    then you subtract the two using the new event codes you supplied as
%    the reference.
%
% Optional Input:
%
%   calc - 'weighted' or 'avg' - specifies whether to perform a weighted
%   average or a straight average when combining conditions using the '+'
%   operation
%   reference - 'events' or 'conditions' - are the numbers in
%               'combinations' refering to event_codes or condition
%               numbers. (default = 'events');
%
%  Created:       10/26/2007 by Rajan Patel
%  Last Modified: 01/10/2007 by Rajan Patel
%
% See also: ts_iEEG_comb_conditions

%% Check Options
if ~mmil_check_nargs(nargin, 2), return; end;

opt = mmil_args2parms(varargin,...
                      {'f_path',       [],[],...
                       'avg_name'    , [],[],...
                       'epoch_name'  , [],[],...
                       'timefreq_name',[],[],...
                       'combinations', [],[],...
                       'neweventcodes',[],[],...
                       'reference'   , 'events',{'events','conditions'},...
                       'calc'      ,'weighted',{'weighted','avg'},...
                      },...
                      true);  
                    
run_field = opt.reference;

if ~exist(opt.f_path,'dir')
  error('Cannot find the given data directory: %s',opt.f_path);
end

add_calc = opt.calc;

if isempty(opt.combinations)
  error('There have been no combinations supplied.');
elseif length(opt.combinations) ~= length(opt.neweventcodes)
  fprintf('Number of Combinations: %s.\nNumber of new event codes supplied: %s.\n',num2str(length(opt.combinations)),num2str(length(opt.neweventcodes)));
  error('You must supply the same number of new event codes as combinations that will be created.');
end

error_list = {};
for i = 1:length(opt.combinations)
  if ~isempty(find(opt.combinations{i} ~= '+' & opt.combinations{i} ~= '-' & ~ismember(opt.combinations{i},['0':'9'])))
    error_list{end+1} = sprintf('The following equation contains invalid characters: %s.\n',opt.combinations{i});
  elseif find(opt.combinations{i}=='-')
    if (~all(find(opt.combinations{i}=='-') == find(~ismember(opt.combinations{i},'0':'9')))) || (length(find(opt.combinations{i}=='-'))~=1)
     error_list{end+1} = sprintf('You can only subtract between two conditions, the following is invalid: %s.\n',opt.combinations{i});
    end
  end
end
if ~isempty(error_list)
  error('%s',error_list{:});
end


%% Combine Averages First

if ~isempty(opt.avg_name)
  avg_file = fullfile (opt.f_path,opt.avg_name);
  if exist(avg_file,'file') && strcmp(who('-file',avg_file),'avg_data')
    fprintf('Loading: %s\n',avg_file);
    load(avg_file);
    numberofconditions = length(avg_data.averages);
    for i = 1:length(opt.combinations)                                     % go through each combo to be made
      fprintf('Calculating new average event %s: %s\n',num2str(opt.neweventcodes(i)),opt.combinations{i});
      curr_combo = opt.combinations{i};
      op_locs = find((~ismember(curr_combo,['0':'9']))==1);               % find out where the operation signs are
      inds = [];
      for j = 1:length(op_locs)                                           % extract the numbers of the events/conditions
        if j == 1
          inds(j) = str2num(curr_combo(1:op_locs(j)-1));
        else
          inds(j) = str2num(curr_combo(op_locs(j-1)+1:op_locs(j)-1));
        end
      end
      inds(end+1) = str2num(curr_combo(op_locs(j)+1:end));
      if curr_combo(op_locs(1)) == '-'                                    % set the operation to be performed
        curr_calc = 'sub';
      else
        curr_calc = add_calc;
      end
      avg_data = ts_calc_combine_conditions(avg_data,...
        run_field,inds,...
        'calc',curr_calc,...
        'eventcode',opt.neweventcodes(i));
    end
    fprintf('Saving data:\n');
    new_file = strrep(opt.avg_name,'avg.mat','combo.avg.mat');
    fprintf('%s\n',fullfile(opt.f_path,new_file));
    save(fullfile(opt.f_path,new_file),'avg_data');
  else
    fprintf('Cannot find average data file: %s.\n',avg_file);
  end
end

%% Combine Epochs

if ~isempty(opt.epoch_name)

  if strcmp(run_field,'events')
    if isempty(strfind(opt.epoch_name,'event'))                                                  % make sure the prefix for the file makes sense
      error('When supplying events the epoch name supplied should end with ''event''.');
    else
      opt.epoch_name = opt.epoch_name(1:strfind(opt.epoch_name,'event')+4);                      % remove anything extra in case given with a number
    end
  elseif strcmp(run_field,'conditions')
    if isempty(strfind(opt.epoch_name,'cond'))
      error('When supplying conditions the epoch name supplied should end with ''cond''.');
    else
      opt.epoch_name = opt.epoch_name(1:strfind(opt.epoch_name,'cond')+3);
    end
  end

  prefix = opt.epoch_name;

  for i=1:length(opt.combinations)
    curr_combo = opt.combinations{i};
    op_locs = find((~ismember(curr_combo,['0':'9']))==1);               % find out where the operation signs are
    inds = [];
    for j = 1:length(op_locs)                                           % extract the numbers of the events/conditions
      if j == 1
        inds(j) = str2num(curr_combo(1:op_locs(j)-1));
      else
        inds(j) = str2num(curr_combo(op_locs(j-1)+1:op_locs(j)-1));
      end
    end
    inds(end+1) = str2num(curr_combo(op_locs(j)+1:end));
    if curr_combo(op_locs(1)) == '-'                                    % set the operation to be performed
      curr_calc = 'sub';
    else
      curr_calc = add_calc;
    end
    if strcmp(curr_calc,'sub')                                          % skip combining epochs for subtractions
      fprintf('SKIPPING new event %s: Cannot subtract events for epoch(s): %s\n',num2str(opt.neweventcodes(i)),opt.combinations{i});
    else
      fprintf('Loading:\n');
      for f = 1:length(inds)
        curr_file = sprintf('%s%s.mat',prefix,num2str(inds(f)));
        curr_file = fullfile(opt.f_path,curr_file);
        fprintf('%s\n',curr_file)
        if ~exist(curr_file,'file')
          error('Cannot find the following epoch data file: %s.',curr_file);
        elseif  ~strcmp(who('-file',curr_file),'epoch_data')
          error('This file does not contain epoch_data: %s.',curr_file);
        end
        load(curr_file);
        data(f) = epoch_data;
        clear epoch_data
      end
      fprintf('Combining epochs...\n');
      data = ts_combine_data(data);
      fprintf('Calculating new epoch event %s: %s\n',num2str(opt.neweventcodes(i)),opt.combinations{i});
      conditions = 1:length(inds);
      data = ts_calc_combine_conditions(data,'conditions',conditions,'eventcode',opt.neweventcodes(i));
      % split off new condition
      epoch_data.num_sensors = data.num_sensors;
      epoch_data.sensor_info = data.sensor_info;
      epoch_data.coor_trans  = data.coor_trans;
      epoch_data.sfreq       = data.sfreq;
      epoch_data.noise       = data.noise;
      epoch_data.epochs.event_code = data.epochs(end).event_code;
      epoch_data.epochs.num_trials = data.epochs(end).num_trials;
      epoch_data.epochs.num_rejects= data.epochs(end).num_rejects;
      epoch_data.epochs.time       = data.epochs(end).time;
      epoch_data.epochs.data       = data.epochs(end).data;
      clear data;
      if strcmp(run_field,'conditions') && exist(numberofconditions,'var')
        numberofconditions = numberofconditions+1;
        new_file = sprintf('%s%s.mat',prefix,num2str(numberofconditions));
      else
        new_file = sprintf('%s%s.mat',prefix,num2str(opt.neweventcodes(i)));
      end
      new_file = fullfile(opt.f_path,new_file);
      fprintf('Saving data:\n%s\n',new_file);
      save(new_file,'epoch_data');
      clear epoch_data
    end
  end % for loop of combinations
end  % combine epochs

%% Combine timefreq_data

if ~isempty(opt.timefreq_name)

  if strcmp(run_field,'events')
    if isempty(strfind(opt.timefreq_name,'event'))                                                  % make sure the prefix for the file makes sense
      error('When supplying events the timefreq name supplied should end with ''event''.');
    else
      opt.timefreq_name = opt.timefreq_name(1:strfind(opt.timefreq_name,'event')+4);                      % remove anything extra in case given with a number
    end
  elseif strcmp(run_field,'conditions')
    if isempty(strfind(opt.timefreq_name,'cond'))
      error('When supplying conditions the timefreq name supplied should end with ''cond''.');
    else
      opt.timefreq_name = opt.timefreq_name(1:strfind(opt.timefreq_name,'cond')+3);
    end
  end

  prefix = opt.timefreq_name;

  for i=1:length(opt.combinations)
    curr_combo = opt.combinations{i};
    op_locs = find((~ismember(curr_combo,['0':'9']))==1);               % find out where the operation signs are
    inds = [];
    for j = 1:length(op_locs)                                           % extract the numbers of the events/conditions
      if j == 1
        inds(j) = str2num(curr_combo(1:op_locs(j)-1));
      else
        inds(j) = str2num(curr_combo(op_locs(j-1)+1:op_locs(j)-1));
      end
    end
    inds(end+1) = str2num(curr_combo(op_locs(j)+1:end));
    if curr_combo(op_locs(1)) == '-'                                    % set the operation to be performed
      curr_calc = 'sub';
    else
      curr_calc = add_calc;
    end
      fprintf('Loading:\n');
      for f = 1:length(inds)
        curr_file = sprintf('%s%s.mat',prefix,num2str(inds(f)));
        curr_file = fullfile(opt.f_path,curr_file);
        fprintf('%s\n',curr_file)
        if ~exist(curr_file,'file')
          error('Cannot find the following timefreq data file: %s.',curr_file);
        elseif  ~strcmp(who('-file',curr_file),'timefreq_data')
          error('This file does not contain timefreq_data: %s.',curr_file);
        end
        load(curr_file);
        data(f) = timefreq_data;
        clear timefreq_data
      end
      fprintf('Combining timefreq...\n');
      data = ts_combine_data(data);
      fprintf('Calculating new timefreq event %s: %s\n',num2str(opt.neweventcodes(i)),opt.combinations{i});
      conditions = 1:length(inds);
      data = ts_calc_combine_conditions(data,'conditions',conditions,'eventcode',opt.neweventcodes(i),'calc',curr_calc);
      % split off new condition
      timefreq_data.num_sensors = data.num_sensors;
      timefreq_data.sensor_info = data.sensor_info;
      timefreq_data.coor_trans  = data.coor_trans;
      timefreq_data.sfreq       = data.sfreq;
      timefreq_data.noise       = data.noise;
      timefreq_data.timefreq.event_code = data.timefreq(end).event_code;
      timefreq_data.timefreq.num_trials = data.timefreq(end).num_trials;
      timefreq_data.timefreq.num_rejects= data.timefreq(end).num_rejects;
      timefreq_data.timefreq.time       = data.timefreq(end).time;
      timefreq_data.timefreq.frequencies= data.timefreq(end).frequencies;
      timefreq_data.timefreq.data       = data.timefreq(end).data;
      timefreq_data.timefreq.power      = data.timefreq(end).power;
      clear data;
      if strcmp(run_field,'conditions') && exist(numberofconditions,'var')
        numberofconditions = numberofconditions+1;
        new_file = sprintf('%s%s.mat',prefix,num2str(numberofconditions));
      else
        new_file = sprintf('%s%s.mat',prefix,num2str(opt.neweventcodes(i)));
      end
      new_file = fullfile(opt.f_path,new_file);
      fprintf('Saving data:\n%s\n',new_file);
      save(new_file,'timefreq_data');
      clear timefreq_data
  end % for loop of combinations
end  % combine epochs