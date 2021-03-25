function [evnts,new_evnts]=ts_recode_events(evnts,recode_rules,sfreq)
%function [evnts,new_evnts]=ts_recode_events(evnts,recode_rules,[sfreq])
%
% Required Input:
%   evnts: structure containing events
%     events.type      = string
%     events.latency   = expressed in samples, first sample of file is 1
%     events.condition = numeric event code
%     events.duration  = duration of event, expressed in samples
%
%   recode_rules: string, or cell array of strings, of recoding rules
%     Each rule specifies a sequence of events and the new event code which is
%       assigned to the first event in the seqeuence
%     Examples:
%       '2->5=101' : event 2 followed directly by event 5 gets recoded as event 101
%       '[2,3,4]->5=101' : events 2, 3, or 4, followed by event 5
%       '[2]->[5,10]=101' : event 2 followed by event 5 or 10
%       '2->200:1500[5]=101' : event 2 followed by event 5 within 200 and 1500 msec
%     '->' separates events in the sequence
%     '=' precedes the new event code
%     ':' is used to separate start and end times (in msec) of time window
%     Lists of event codes are placed in square brackets [], separated by commas
%     Square brackets are also required when an event code follows a time window
%
% Optional Input:
%   sfreq: sampling frequency (Hz)
%     Used to convert time ranges to samples
%     {default = 1000}
%
% Output:
%   evnts: events structure containing all events with recodings (replacements)
%   new_evnts: events structure containing only recoded events
%   
% Feaures currently not supported:
%   "loose" sequence of events where an intervening event can be ignored
%   specifying a negative time window
%     (recoding an event if a particular event precedes it)
%
% Recoding is done independently for each rule, so a given event could be recoded
%   more than once (e.g. to combine conditions in different ways)
%
% See also: ts_read_fif_events
%
% Created:  06/25/07 by Ben Cipollini
% Last Mod: 03/23/11 by Don Hagler
%

if ~mmil_check_nargs(nargin, 2), return; end;

new_evnts = [];

if (~iscell(recode_rules)) recode_rules = {recode_rules}; end;
recode_rules_obj = ts_recode_events_parse_rules(evnts,recode_rules);
for i=1:length(recode_rules_obj)
  [evnts,tmp_new_evnts] = ts_recode_events_internal(evnts,...
    recode_rules_obj(i).event_sequence,sfreq);
  new_evnts = [new_evnts,tmp_new_evnts];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recode_rules_obj=ts_recode_events_parse_rules(evnts,recode_rules_str)
%
% recode_rules_str: must be cell array
%
  if (~iscell(recode_rules_str))
    error('Programming error: recode_rules_str must be a cell array.');
  end;
  recode_rules_obj = [];
  % syntax allows us to use mmil_splitstr A LOT
  for i=1:length(recode_rules_str)
    recode_rules_obj(i).event_sequence = [];
  
    % find match and assign portions of string.
    tmp = mmil_splitstr(recode_rules_str{i},'=');
    if (length(tmp) ~= 2)
      error('Parse error when trying to find ''='' for rule ''%s''',...
        recode_rules_str{i});
    end;
    match_str = tmp{1};
    assign_str = tmp{2};
    %find transition sequence
    match_rules = mmil_splitstr(match_str,'->');
    %% todo: allow loose sequences that allow for intervening events
    %%       specified by "," rather than "->"
    for j=1:length(match_rules)
      % rule checking is in order from lowest to highest priority,
      % so that reassignment of the 'method' ends with the assignment
      % of the highest priority method.
      cur_rule = [];
      % loop through all match rules; stop when we've found an applicable one
      tmp = {}; iRule = 1;
      while (isempty(tmp))
        switch (iRule)
          % find events
          case 1
            tmp=regexp(match_rules{j}, '^\[?([0-9,\s]+)\]?$', 'tokens');
            if (~isempty(tmp))
              cur_rule.conditions = double(mmil_cell2num(mmil_splitstr(tmp{1}{1},',')));
              if j==1
                cur_rule.mode = 'loose'; % does not have to be very next event
              else
                cur_rule.mode = 'strict';
              end;
              cur_rule.method = @match_by_sequence;
              cur_rule.epoch = [];
          end;
          % find epoch
          case 2
if 0 % this produces the strange behavior of requiring an intervening
     % event within the time window, but not caring which event
     %  according to Ben's scheme, that behavior should come from a rule like this:
     %    1->200:1500[*]=4 
     %  and 1->200:1500->[3]=4 should require no intervening events for that
     %    time window
            tmp=regexp(match_rules{j}, '^([0-9]*):([0-9]*)$', 'tokens');
            if (~isempty(tmp))
              cur_rule.conditions = [];
              cur_rule.mode = 'relative';
              cur_rule.method = @match_by_epoch;
              cur_rule.epoch = [str2double(tmp{1}{1}) str2double(tmp{1}{2})];
            end;
end;
          % find epoch AND events
          case 3
            tmp=regexp(match_rules{j},'^([0-9.]*):([0-9.]*)\[([0-9,\s]+)\]$', 'tokens');
            if (~isempty(tmp))
              cur_rule.conditions = mmil_cell2num(mmil_splitstr(tmp{1}{3:end},','));
              cur_rule.mode = 'relative';
              cur_rule.method = @match_by_epoch;
              cur_rule.epoch = [str2double(tmp{1}{1}) str2double(tmp{1}{2})];
            end;  
          otherwise
            tmp = 'failure!';
            error('failed to parse rule ''%s''.',recode_rules_str{i});
        end;
        iRule = iRule+1;
      end; % loop over match rules
      clear('tmp','iRule');
      % recode only first in sequence
      if j==1
        cur_rule.assign = str2num(assign_str);
      else
        cur_rule.assign = [];
      end;
      %% todo: how to recode only one event but allow it to be other than first?
      if (isempty(recode_rules_obj(i).event_sequence))
        recode_rules_obj(i).event_sequence = cur_rule;
      else
        recode_rules_obj(i).event_sequence(end+1) = cur_rule;
      end;
    end;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evnts,new_evnts]=ts_recode_events_internal(evnts,event_sequence,sfreq,cur_evt_idx)
  %	Rules define a procedure and parameters for FINDING an event sequence
  % This function is called recursively to check each event

  new_evnts = [];

  %	base case: no recode rules.  Just return as-is.
  if (isempty(event_sequence))
    return;
  elseif (~exist('cur_evt_idx','var'))
    cur_evt_idx = 1;
  end;

  % call the appropriate function to find the matches for the current rule.
  cur_rule = event_sequence(1);
  if (ischar(cur_rule.method))
    method = str2func(cur_rule.method);
  else
    method = cur_rule.method;
  end;
  matches = method(evnts,sfreq,event_sequence(1),cur_evt_idx);

  % if no matches to rule, return empty
  % but if this is the original call, return evnts struct unaltered
  if (isempty(matches))
    if (cur_evt_idx ~=1)
      evnts = [];
    end;
    return;
  end;

  %	recode if any full matches were found by "method"
  for i=1:length(matches)
    % do a recursive call starting just after the found event with the next rule
    %   this checks for matches to the entire sequence of events in one rule
    [tmp_evnts,tmp_new_evnts] = ...
      ts_recode_events_internal(evnts,event_sequence(2:end),sfreq,matches(i)+1);
    %	found a successful path; recode if necessary.
    if (~isempty(tmp_evnts))
      evnts	= tmp_evnts;
      if (isfield(cur_rule, 'assign') && ~isempty(cur_rule.assign))
        if (length(cur_rule.assign)==1)
          % assign every matched code to the same condition
          evnts(matches(i)).condition = cur_rule.assign;
        elseif (length(cur_rule.conditions)==length(cur_rule.assign))
          % assign a different code for different conditions
          %  (recode multiple conditions with a single rule)
          [m,loc] = ismember([evnts(matches(i)).condition],cur_rule.assign);
          evnts(matches(i)).condition = cur_rule.assign(loc);
        end;
        % copy new event
        if isempty(new_evnts)
          new_evnts = evnts(matches(i));
        else
          new_evnts(end+1) = evnts(matches(i));
        end;
      end;
    end;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epoch_event_idx	= match_by_epoch(evnts,sfreq,recode_rule,cur_evt_idx)
%
%	required:
%		latency	- in ms [min max]
%	optional:
%		conditions	- acceptable conditions codes (default: all)
%

%fprintf('%s: cur_evt_idx = %d, code = %d\n',...
%  mfilename,cur_evt_idx,evnts(cur_evt_idx).condition);

	if (~isfield(recode_rule,'epoch') || isempty(recode_rule.epoch))
		error('Specified match_in_time with no epoch set');
	elseif (~isfield(recode_rule,'mode') || isempty(recode_rule.mode))
    recode_rule.mode    = 'relative';
	end;

  switch recode_rule.mode
    case 'relative'
      %   When cur_evt_idx=1, relative is the same as absolute
      if (cur_evt_idx~=1)
        recode_rule.epoch   = recode_rule.epoch + (evnts(cur_evt_idx-1).latency-1)/(sfreq/1000);
      end;
    case 'absolute'
        %   Nothing to do
  end;
	%	grab all events within the epoch
	epoch_event_idx	= intersect(find( (([evnts(cur_evt_idx:end).latency]-1)/(sfreq/1000)) >= recode_rule.epoch(1) ), ...
                              find( (([evnts(cur_evt_idx:end).latency]-1)/(sfreq/1000)) <= recode_rule.epoch(2) ));

  % make event indices relative to cur_evt_idx into absolute indices
  %  (relative to 1)
  epoch_event_idx = epoch_event_idx + cur_evt_idx - 1;

  % limit the found events to a list
  if ~isfield(recode_rule,'conditions') || isempty(recode_rule.conditions)
    % allow any events
  elseif ~isempty(epoch_event_idx)
    % only look at next event
    epoch_event_idx = epoch_event_idx(1);
    if ~ismember(evnts(epoch_event_idx).condition,recode_rule.conditions)
      epoch_event_idx = [];
    end;
	end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evnt_idx	= match_by_sequence(evnts,sfreq,recode_rule,cur_evt_idx)
%
%	required:
%		conditions	- acceptable conditions codes (default: all)
%		mode - either of two values:
%     strict (e.g. the next event) 
%     loose (any number of events can intervene)
%

  if (~isfield(recode_rule,'conditions') || isempty(recode_rule.conditions))
		error('Specified match_in_space with no conditions set');
	elseif (~isfield(recode_rule,'mode') || isempty(recode_rule.mode))
    recode_rule.mode    = 'strict';
	end;

  % default 'found' indices to empty
  evnt_idx     = [];

	%	limit the found events to a list
  switch recode_rule.mode
    case 'strict'
      if (isempty(find(ismember([evnts(cur_evt_idx).condition], recode_rule.conditions))))
        evnt_idx     = [];
      else
        evnt_idx     = cur_evt_idx;
      end;
    case 'loose'
      evnt_idx	= find(ismember([evnts(cur_evt_idx:end).condition], recode_rule.conditions));
      evnt_idx    = evnt_idx + cur_evt_idx - 1;
    otherwise
      error('unknown recode mode: %s',recode_rule.mode);
  end;

