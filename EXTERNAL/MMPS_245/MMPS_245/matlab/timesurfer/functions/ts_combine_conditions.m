function out_data = ts_combine_conditions(in_data,varargin)
%function out_data = ts_combine_conditions(in_data,[options])
%
% Purpose: combine across various conditions using weighted
%  or unweighted averages; can also subtract two conditions.
%
% Usage: out_data = ts_combine_conditions(in_data,'key1','value1',...)
%
% Required input:
%   avg_data: valid TimeSurfer avg_data structure
%      OR
%   epoch_data: valid TimeSurfer epoch_data structure
%      OR
%   timefreq_data: valid TimeSurfer timefreq_data structure
%   combinations: cell array of combinations to produce 
%     Each combination is supplied as a string.
%       The numbers within the string are treated either as event codes
%       or condition numbers (if 'reference' = 'events' or 'conditions')
%     Valid examples are:
%       To combine events 1,2 & 3:   '1+2+3'
%       To subtract event 10 from 2: '2-10' (Note: Only two events at a time
%                                               is valid for subtractions)
%   
%    Only simple operations are available.  You cannot for instance specify
%    a combination that includes both combining and then subtracting.  There is
%    a workaround.  For instance...
%
%    If the data contained 10 events with event codes 1-10 and you wanted to
%    do the following: (1+2+3+4+5)-(6+7+8+9+10) the settings would be -
%      events       : {'1+2+3+4+5','6+7+8+9+10','11-12'}
%      neweventcodes: [11 12 13]
%    So that you create the two averages across multiple events first and
%    then you subtract the two using the new event codes you supplied as
%    the reference.
%
% Optional Input:
%   'calc': 'weighted', 'avg', or 'sum' - specify whether to perform
%     a weighted average, straight average, or addition when combining
%     conditions with '+' operation
%     {default = 'weighted'}
%   'reference': 'events' or 'conditions' - should the numbers specified in
%     the combinations be treated as event codes or condition number
%     {default = 'events'}
%   'neweventcodes': a list of new event codes to assign each of the new
%     combinations, if none is supplied it will default creating new event
%     codes starting with the largest current event code in the data
%     structure conditions.  It is recommended to specify your own.
%
% Created    01/09/08 by Rajan Patel
% Last Mod:  10/24/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if ~mmil_check_nargs(nargin, 2), return; end;
parms = mmil_args2parms(varargin,{...
  'combinations',[],[],...
  'calc','weighted',{'weighted','avg','sum'},...
  'reference','events',{'events','conditions'},...
  'neweventcodes',[],[],...
  'verbose',true,sort([false true]),...
  'logfile',[],[],...
  'logfid',1,[],...
},false);                    
   
in_data = ts_checkdata_header(in_data,'events',[]);               
data_type = ts_objecttype(in_data);

switch data_type
  case 'averages'
    data_field = 'averages';
	case 'average'
    data_field = 'average';
  case 'epoch'
    data_field = 'epochs';
	case 'epochs'
    data_field = 'epochs';
  case 'timefreq'
    data_field = 'timefreq';
  otherwise
    mmil_error(parms,'input object must be a valid epoch, avg or timefreq data structure.');
end
  
if isempty(parms.combinations)
  mmil_error(parms,'No set of combinations provided to produce.');
end

if ~iscell(parms.combinations), parms.combinations = {parms.combinations}; end

if ~isempty(parms.neweventcodes) && (length(parms.neweventcodes) ~= length(parms.combinations))
  mmil_error(parms,'The number of new event codes and number of combinations do not match.');
elseif isempty(parms.neweventcodes)
  largesteventcode = max([in_data.(data_field).event_code]);
  parms.neweventcodes = (largesteventcode+1):(largesteventcode + length(parms.combinations));
end

error_list = {};
for i = 1:length(parms.combinations)
  if ~isempty(find(parms.combinations{i} ~= '+' & parms.combinations{i} ~= '-' & ~ismember(parms.combinations{i},['0':'9'])))
    error_list{end+1} = sprintf('The following equation contains invalid characters: %s.\n',parms.combinations{i});
  elseif find(parms.combinations{i}=='-')
    if (~all(find(parms.combinations{i}=='-') == find(~ismember(parms.combinations{i},'0':'9')))) || (length(find(parms.combinations{i}=='-'))~=1)
     error_list{end+1} = sprintf('You can only subtract between two conditions, the following is invalid: %s.\n',parms.combinations{i});
    end
  end
end
if ~isempty(error_list)
  error('%s:\n%s',mfilename,error_list{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data

for i = 1:length(parms.combinations)                                     % go through each combo to be made
  mmil_logstr(parms,'Calculating new %s event %s: %s',data_type,num2str(parms.neweventcodes(i)),parms.combinations{i});
  curr_combo = parms.combinations{i};
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
    curr_calc = parms.calc;
  end
  if strcmp(curr_calc,'sub') && (strcmp(data_type,'epoch'))               % skip combining subtraction for epochs
    mmil_logstr(parms,'SKIPPING new event %s: Cannot subtract events using epoch data: %s',num2str(parms.neweventcodes(i)),parms.combinations{i});
  else
    in_data = ts_calc_combine_conditions(in_data,...
      parms.reference,inds,...
      'calc',curr_calc,...
      'eventcode',parms.neweventcodes(i));
  end
end

out_data = in_data;
