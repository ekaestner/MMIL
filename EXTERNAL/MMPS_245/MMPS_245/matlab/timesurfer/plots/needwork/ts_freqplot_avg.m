function freq=ts_freqplot_avg(avg_data,varargin);
% ts_freqplot_avg - calculate FFT of average data and make color scale plot
%
% Usage:
%   freq=ts_freqplot_avg(avg_data,'key1', value1,...);
%
% Required input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  condition - condition number (not event code) used to index avg_data.averages
%    {default: 1}
%  time0 - start time of data range (msec)
%    {default: -100}
%  time1 - end time of data range (msec)
%    {default: 1000}
%  badchanfile - name of text file containing bad channel labels
%    {default: []}
%  channel - channel index
%    {default: 1}
%  scale_max - color scale max value
%    {default: min/max}
%
% Output: 
%   freq - field trip frequency 
%
%  created:       06/29/06   by Don Hagler
%  last modified: 07/31/06   by Don Hagler
%

%% todo: chantype option
%% todo: use multiplotTFR


freq = [];

DEFAULT_CHANNEL = 1;
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 100;
DEFAULT_TIME1 = 200;
DEFAULT_SCALE_MAX = 0;
DEFAULT_CONDITION = 1;

DEFAULT_LOWFREQ = 20;
DEFAULT_HIFREQ = 70;
DEFAULT_WAVELETWIDTH = 5;

if nargin < 1
  help(mfilename);
  return;
end;

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

try, opt.channel;       catch, opt.channel = DEFAULT_CHANNEL; end;
try, opt.badchanfile;   catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;         catch, opt.time0 = DEFAULT_TIME0; end;
try, opt.time1;         catch, opt.time1 = DEFAULT_TIME1; end;
try, opt.scale_max;     catch, opt.scale_max = DEFAULT_SCALE_MAX; end;
try, opt.condition;     catch, opt.condition = DEFAULT_CONDITION; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'channel' 'badchanfile' 'time0' 'time1'...
         'scale_max' 'condition'...
   },;
   otherwise, error([mfilename ': unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
channel = opt.channel;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
scale_max = opt.scale_max;
condition = opt.condition;
clear opt optfields options;

% convert to field trip format
FT_data = ts_avg2fieldtrip(avg_data,'condition',condition,'channels',channel,...
   'badchanfile',badchanfile,'time0',time0,'time1',time1);

if isempty(FT_data)
  fprintf('%s: error converting to fieldtrip - bad chan selected?\n',mfilename);
  return;
end;

% set configuration for freqanalysis_tfr
cfg = [];
cfg.method = 'tfr';
f1 = DEFAULT_LOWFREQ;
f2 = DEFAULT_HIFREQ;
cfg.foi = [f1:2:f2];
cfg.waveletwidth = DEFAULT_WAVELETWIDTH;

% calculate power spectrum
freq = freqanalysis(cfg,FT_data);
X = freq.time*1000;
Y = freq.freq;
C = squeeze(freq.powspctrm);

switch avg_data.sensor_info(channel).typestring
case {'grad1' 'grad2'}
  units = (10^-13)^2;
case 'mag'
  units = (10^-15)^2;
case 'eeg'
  units = (10^-6)^2;
end; 

%% todo: use multiplotTFR

if scale_max > 0
  clim = [0 scale_max]*units;
  imagesc(X,Y,C,clim)
else
  imagesc(X,Y,C)
end;
xlabel('time');
ylabel('freq');
axis xy;
colorbar;

return;

