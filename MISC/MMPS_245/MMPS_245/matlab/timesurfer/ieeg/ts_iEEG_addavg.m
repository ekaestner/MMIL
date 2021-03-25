function out_data = ts_iEEG_addavg (in_data,varargin)
% Usage: out_data = ts_iEEG_addavg (in_data,'option',value....
%
% Adds two averages together and creates a new event and adds it to the avg
% data structure.  
%
% Required Input:
%
% in_data      - a valid TimeSurfer avg_data structure.
% combinations - this is the list of combinations to produce 
%                provided as a cell array.  Each combination is
%                supplied as a string.  The numbers within the string are
%                treated either as event codes (if using the reference =
%                'events') or as condition numbers (if using the reference =
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
%  
%  In this script the combinations will be additions not weighted averages.
%
% Optional Input:
%
% reference - 'events' or 'conditions' - should the numbers specified in
%   the combinations be treated as event codes or condition number (default
%   = 'events')
% neweventcodes - a list of new event codes to assign each of the new
%   combinations, if none is supplied it will default creating new event
%   codes starting with the largest current event code in the data
%   structure conditions.  It is recommended to specify your own.
%
% Created 05/29/2008 by Rajan Patel
% Last Modified 05/29/2008 by Rajan Patel

% TO DO: Calculate correct standard deviation
%% Check Inputs 

if nargin < 1, help(mfilename); end

opt = mmil_args2parms(varargin,...
                     {'combinations',[],[],...
                      'reference','events',{'events','conditions'},...
                      'neweventcodes',[],[],...
                      'verbose',true,sort([false true]),...
                      'logfile',[],[],...
                      'logfid',1,[],...
                      },...
                      false);                    

data_type = ts_objecttype(in_data);
      
switch data_type
    case 'average' 
         datafield = 'averages';
		case 'averages'
			datafield = 'averages';
%     case 'epoch'
%         data_field = 'epochs';
%     case 'timefreq'
%         data_field = 'timefreq';
    otherwise
     mmil_error(opt,'Input object must be a valid avg data structure.');
end                        
                  
if ~iscell(opt.combinations), opt.combinations = {opt.combinations}; end

if ~isempty(opt.neweventcodes) && (length(opt.neweventcodes) ~= length(opt.combinations))
  mmil_error(opt,'The number of new event codes and number of combinations do not match.');
elseif isempty(opt.neweventcodes)
  largesteventcode = max([in_data.(datafield).event_code]);
  opt.neweventcodes = (largesteventcode+1):(largesteventcode + length(opt.combinations));
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
  error('%s:\n%s',mfilename,error_list{:});
end

out_data = in_data;

%% Process Data

for i = 1:length(opt.combinations)                                     % go through each combo to be made
    mmil_logstr(opt,'Calculating new %s event %s: %s',data_type,num2str(opt.neweventcodes(i)),opt.combinations{i});
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
        curr_calc = 'add';
    end
    if strcmpi(opt.reference,'events')
        conditions = [];
        for j = 1:length(out_data.(datafield))                            % convert events to conditions
            if ismember(out_data.(datafield)(j).event_code,inds)
                conditions = cat(1,conditions,j);
            end
        end
        if length(conditions) ~= length(inds)                      % make sure found all events
            error('Invalid event codes given.');
        end
    elseif ~isempty(opt.conditions)
        conditions = inds;
    end
    if strcmp(curr_calc,'add')
        out_data.(datafield)(end+1).event_code       = opt.neweventcodes(i);
        out_data.(datafield)(end).time               = out_data.(datafield)(conditions(1)).time;
        out_data.(datafield)(end).stdev              = nan(size(out_data.(datafield)(conditions(1)).stdev));
        out_data.(datafield)(end).num_trials         = 0;
        out_data.(datafield)(end).data               = zeros(size(out_data.(datafield)(conditions(1)).data));
        out_data.(datafield)(end).num_rejects.mag    = 0;
        out_data.(datafield)(end).num_rejects.grad   = 0;
        out_data.(datafield)(end).num_rejects.eeg    = 0;
        out_data.(datafield)(end).num_rejects.eog    = 0;
        out_data.(datafield)(end).num_rejects.manual = 0;
        out_data.(datafield)(end).num_rejects.skip   = 0;
        for c = 1:length(conditions)
            out_data.(datafield)(end).num_trials         = out_data.(datafield)(end).num_trials + ...
                                                        out_data.(datafield)(conditions(c)).num_trials;
            out_data.(datafield)(end).data               = out_data.(datafield)(end).data + ...
                                                        out_data.(datafield)(conditions(c)).data;
            out_data.(datafield)(end).num_rejects.mag    = out_data.(datafield)(end).num_rejects.mag + ...
                                                        out_data.(datafield)(conditions(c)).num_rejects.mag;
            out_data.(datafield)(end).num_rejects.grad   = out_data.(datafield)(end).num_rejects.grad + ...
                                                        out_data.(datafield)(conditions(c)).num_rejects.grad;
            out_data.(datafield)(end).num_rejects.eeg    = out_data.(datafield)(end).num_rejects.eeg + ...
                                                        out_data.(datafield)(conditions(c)).num_rejects.eeg;
            out_data.(datafield)(end).num_rejects.eog    = out_data.(datafield)(end).num_rejects.eog + ...
                                                        out_data.(datafield)(conditions(c)).num_rejects.eog;
            out_data.(datafield)(end).num_rejects.manual = out_data.(datafield)(end).num_rejects.manual + ...
                                                        out_data.(datafield)(conditions(c)).num_rejects.manual;
            out_data.(datafield)(end).num_rejects.skip   = out_data.(datafield)(end).num_rejects.skip + ...
                                                        out_data.(datafield)(conditions(c)).num_rejects.skip;            
        end
    elseif strcmp(curr_calc,'sub')
        if length(conditions)~=2,
                error('%s: To calculate subtraction conditions only two condition numbers is appropriate.',mfilename);
        end      
        out_data.(datafield)(end+1).event_code       = opt.neweventcodes(i);
        out_data.(datafield)(end).time               = out_data.(datafield)(conditions(1)).time;
        out_data.(datafield)(end).stdev              = nan(size(out_data.(datafield)(conditions(1)).stdev));
        out_data.(datafield)(end).num_trials         = round(geomean(out_data.(datafield)(conditions(1)).num_trials + ...
                                                        out_data.(datafield)(conditions(2)).num_trials));
        out_data.(datafield)(end).data               = out_data.(datafield)(conditions(1)).data - ...
                                                        out_data.(datafield)(conditions(2)).data;
        out_data.(datafield)(end).num_rejects.mag    = out_data.(datafield)(conditions(1)).num_rejects.mag + ...
                                                        out_data.(datafield)(conditions(2)).num_rejects.mag;
        out_data.(datafield)(end).num_rejects.grad   = out_data.(datafield)(conditions(1)).num_rejects.grad + ...
                                                        out_data.(datafield)(conditions(2)).num_rejects.grad;
        out_data.(datafield)(end).num_rejects.eeg    = out_data.(datafield)(conditions(1)).num_rejects.eeg + ...
                                                        out_data.(datafield)(conditions(2)).num_rejects.eeg;
        out_data.(datafield)(end).num_rejects.eog    = out_data.(datafield)(conditions(1)).num_rejects.eog + ...
                                                        out_data.(datafield)(conditions(2)).num_rejects.eog;
        out_data.(datafield)(end).num_rejects.manual = out_data.(datafield)(conditions(1)).num_rejects.manual + ...
                                                        out_data.(datafield)(conditions(2)).num_rejects.manual;
        out_data.(datafield)(end).num_rejects.skip   = out_data.(datafield)(conditions(1)).num_rejects.skip + ...
                                                        out_data.(datafield)(conditions(2)).num_rejects.skip;
    else
    end
end  


  
  
