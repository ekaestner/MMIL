function [FT_data] = ts_iEEG_avg2fieldtrip(avg_data,varargin);
% ts_iEEG_avg2fieldtrip - converts the avg_data TimeSurfer structure
%                         into a FieldTrip FT_data structure, which subsequently
%                         can be used for topoplotER and other FieldTrip tools
%
% Usage:
%  [FT_data] = ts_iEEG_avg2fieldtrip(avg_data,condition);
%
% Required input:
%  avg_data - average data structure
%
% Optional parameters:
%  condition - condition number (not event code) used to index avg_data.averages
%    {default: 1}
%  time0 - start time of data range (msec)
%    {default: first time point}
%  time1 - end time of data range (msec)
%    {default: last time point}
%  badchanfile - name of text file containing bad channel labels
%    {default: []}
%  chantype - channel type (e.g. 'mag','grad1','grad2','eeg', 'grad', or 'other')
%    {default: 'grad'} ('grad' returns both 'grad1' and 'grad2' channels)
%  channels - vector of channel indices - overrides chantype
%    {default: []}
%
% Output:
%   FT_data - FIELDTRIP structure
%
%
% based on FieldTrip's eeglab2fieldtrip by Robert Oostenveld
% See http://www.ru.nl/fcdonders/fieldtrip/
%
%  created:       04/20/06   by Don Hagler
%  last modified: 04/27/11   by Don Hagler
%

DEFAULT_BADCHANFILE = [];
DEFAULT_CONDITION = 1;
DEFAULT_CHANNELS = [];
DEFAULT_CHANTYPE = 'eeg';

if nargin < 1
  help(mfilename);
  return;
end;

FT_data = [];

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), opt=struct(options{:}); 
  else opt = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

try, opt.channels;      catch, opt.channels = DEFAULT_CHANNELS; end;
try, opt.chantype;      catch, opt.chantype = DEFAULT_CHANTYPE; end;
try, opt.badchanfile;   catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;         catch, opt.time0 = []; end;
try, opt.time1;         catch, opt.time1 = []; end;
try, opt.condition;     catch, opt.condition = DEFAULT_CONDITION; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'channels' 'chantype' 'badchanfile' 'time0' 'time1'...
         'condition'...
   }
   otherwise, error(['ts_avg_freqplot: unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
channels = opt.channels;
chantype = opt.chantype;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
condition = opt.condition;
clear opt optfields options;

if ~isempty(badchanfile)
  badchan_i = ts_read_txt_badchans(badchanfile,{avg_data.sensor_info.label});
else
  badchan_i = [];
end;

if length(condition)~=1 | ~isnumeric(condition) |...
   ~mmil_isint(condition) | condition<1
  fprintf('%s: condition must be a single integer >= 1\n',mfilename);
  return;
end;
numconds = length(avg_data.averages);
if condition > numconds
  fprintf('%s: condition must be <= number of conditions (%d)\n',...
    mfilename,numconds);
  return;
end;

if isempty(channels)
  if ~isstr(chantype)
    help(mfilename)
    fprintf('%s: chantype must be a string\n',mfilename);
    return;
  end;
end;

% set start and end samples
t0 = avg_data.averages(1).time(1)*1000;
t1 = avg_data.averages(1).time(end)*1000;
if isempty(time0)
  time0 = t0;
end;
if isempty(time1)
  time1 = t1;
end;
if time0<t0, time0=t0; end;
if time1>t1, time1=t1; end;
if time0>time1
  t0 = time0;
  t1 = time1;
  time1 = t0;
  time0 = t1;
end;
sfreq = avg_data.sfreq;
t_trigger = avg_data.averages(1).time(1)*1000;
t0 = round((time0 - t_trigger)*sfreq/1000);
t1 = round((time1 - t_trigger)*sfreq/1000);
if t0 <= 1, t0 = 1; end;
if t1 > length(avg_data.averages(1).time)
  t1 = length(avg_data.averages(1).time);
end;

if isempty(channels)
  switch chantype
  case {'mag' 'grad1' 'grad2' 'eeg' 'other'}
    chans = find(strcmp(chantype,{avg_data.sensor_info.typestring}));
  case {'grad'}
    chans = find(strncmp(chantype,{avg_data.sensor_info.typestring},...
      length(chantype)));
  otherwise
    help(mfilename)
    fprintf('%s: unsupported chantype (%s)\n',mfilename);
    return;
  end;
else
  chans = channels;
end;
chans = setdiff(chans,badchan_i);
if isempty(chans)
  fprintf('%s: no good channels selected\n',mfilename);
  return;
end;

FT_data.label = {avg_data.sensor_info(chans).label};
FT_data.fsample = avg_data.sfreq;

if isempty(channels)
  switch chantype
  case {'mag' 'grad1' 'grad2' 'grad' 'other'}
    FT_data.grad.pnt   = zeros(length(chans),3);
    FT_data.grad.label = cell(length(chans),1);
  case 'eeg'
    FT_data.elec.pnt   = zeros(length(chans),3);
    FT_data.elec.label = cell(length(chans),1);
  end;

  for c=1:length(chans)
    k=chans(c);
    loc = avg_data.sensor_info(k).loc;
    switch chantype
    case {'mag' 'grad1' 'grad2' 'grad' 'other'}
      % apply device2head transform
      if ~isempty(avg_data.coor_trans.device2head)
        T = avg_data.coor_trans.device2head;
        loc = T*loc;
      end;
      FT_data.grad.label{c} = avg_data.sensor_info(k).label;
      FT_data.grad.pnt(c,1:3) = loc(1:3,4);
    case 'eeg'
      % locations already in "head" space
      FT_data.elec.label{c} = avg_data.sensor_info(k).label;
      FT_data.elec.pnt(c,1:3) = loc(1:3,4);
    end;
  end;
end;

FT_data.avg  = avg_data.averages(condition).data(chans,t0:t1);
if ~isempty(avg_data.averages(condition).stdev), FT_data.var  = avg_data.averages(condition).stdev(chans,t0:t1); 
else FT_data.var  = zeros(chans,t0:t1); end;
% FT_data.var  = avg_data.averages(condition).stdev(chans,t0:t1);
FT_data.time = avg_data.averages(condition).time(t0:t1);
FT_data.dimord = 'chan_time';

try
  % get the full name of the function
  FT_data.cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  FT_data.cfg.version.name = st(i);
end

% add the version details of this function call to the configuration
FT_data.cfg.version.id   = '$Id: ts_avg2fieldtrip.m,v 0.0 2006/04/20 donh Exp $';

return
