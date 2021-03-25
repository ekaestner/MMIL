function evnts=ts_edit_events(evnts,sfreq,valid_event_codes,code_excl,...
  time_excl_pre,time_excl_post,recode_rules)
%function edited_evnts=ts_edit_events(evnts,sfreq,[valid_event_codes],...
%  [code_excl],[time_excl_pre],[time_excl_post])
%
% Required Input:
%   evnts - structure containing events
%     events.type      = string
%     events.latency   = expressed in samples, first sample of file is 1
%     events.condition = numeric event code
%     events.duration  = duration of event, expressed in samples
%
% Optional Input:
%   sfreq - sampling frequency
%
%   valid_event_codes - vector containing valid event codes
%     events with all other codes will be excluded
%     if empty ([]), all events will be kept
%
%   code_excl - event code or vector of event codes
%     other events surrounding this event will be excluded
%
%   time_excl_pre  - exclusion time (msec) before code_excl
%   time_excl_post - exclusion time (msec) after code_excl
%
%   recode_rules - string, or cell array of strings, of recoding rules 
%                         see ts_recode_events for recode rules syntax.
%
% created:       04/27/06 by Don Hagler
% last modified: 12/17/07 by Ben Cipollini
% last modified: 02/14/08 by Don Hagler
%
% See also: ts_read_fif_events, ts_recode_events
%

if nargin < 2
   help(mfilename);
   return;
end

% initialize optional parameters
if ~exist('valid_event_codes')
  valid_event_codes = [];
end;
if ~exist('code_excl')
  code_excl = [];
end;
if ~exist('time_excl_pre')
  time_excl_pre = 0;
end;
if ~exist('time_excl_post')
  time_excl_post = 0;
end;
if ~exist('recode_rules')
  recode_rules = {};
end;
  
% todo: check for errors in parameters

if isempty(evnts)
  return;
end;

% separate evnts by type
reject_i = find(strcmp({evnts.type},'reject') |...
                strcmp({evnts.type},'skip'));
other_i = setdiff(1:length(evnts),reject_i);
reject_evnts = evnts(reject_i);
evnts = evnts(other_i);

% make sure all events have a condition
for j=1:length(evnts)
  if isempty(evnts(j).condition)
    evnts(j).condition=0;
  end;
end;

% exclude invalid event codes
if ~isempty(valid_event_codes)
  keep_i = find(ismember(cell2mat({evnts.condition}),valid_event_codes));
  evnts = evnts(keep_i);
end;

% recode some event codes to new ones
if ~isempty(recode_rules)
  evnts = ts_recode_events(evnts,recode_rules,sfreq);
end;

%% todo: this should be done by ts_recode_events
% exclude events surrounding code_excl (this will remove code_excl too!)
excl_events = evnts(find(ismember(cell2mat({evnts.condition}),code_excl)));
excl_pre_samp = time_excl_pre * sfreq/1000;
excl_post_samp = time_excl_post * sfreq/1000;
evnts_beg   = cell2mat({evnts.latency});
evnts_end   = evnts_beg + cell2mat({evnts.duration}) - 1;
exclude_i = [];
for j=1:length(excl_events)
  begsamp = excl_events(j).latency - excl_pre_samp;
  endsamp = excl_events(j).latency + excl_post_samp;
  exclude_i = [exclude_i ...
    find((evnts_beg >= begsamp & evnts_beg <= endsamp) | ...
         (evnts_end >= begsamp & evnts_end <= endsamp))];
end;
keep_i = setdiff(1:length(evnts),exclude_i);
evnts = evnts(keep_i);

evnts = [evnts reject_evnts];

return;
