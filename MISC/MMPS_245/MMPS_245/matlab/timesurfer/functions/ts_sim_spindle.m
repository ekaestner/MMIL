function spindle = ts_sim_spindle(ntpoints,freq,phas,t0,td,amp,sf)
%function spindle = ts_sim_spindle(ntpoints,freq,phas,t0,td,amp,sf)
%
% Required input:
%   ntpoints: number of time points
%   freq: frequency (Hz)
%   phas: phase (fraction of cycle)
%   t0: start time of hanning window (msec)
%   td: duration of hanning window (msec)
%   amp: amplitude
%   sf: sampling frequency (Hz)
%   
% Created:  08/30/13 by Don Hagler
% Last Mod: 09/11/13 by Don Hagler
%

spindle = [];

if ~exist('ntpoints','var') || isempty(ntpoints), ntpoints = 1000; end;
if ~exist('freq','var') || isempty(freq), freq = 13; end;
if ~exist('phas','var') || isempty(phas), phas = 0; end;
if ~exist('t0','var') || isempty(t0), t0 = 250; end;
if ~exist('td','var') || isempty(td), td = 500; end;
if ~exist('amp','var') || isempty(amp), amp = 1; end;
if ~exist('sf','var') || isempty(sf), sf = 1000; end;

y = zeros(1,ntpoints);
w = zeros(1,ntpoints);

t = [0:ntpoints-1]/sf;

t0 = t0/1000;
td = td/1000;
t1 = t0 + td;

[tmp,s0] = min(abs(t-t0));
[tmp,s1] = min(abs(t-t1));

if s1>s0
  y = amp*sin(2*pi*(t*freq + phas));
  m = s1-s0+1;
  w(s0:s1) = hanning(m);
  spindle = y.*w;
else
  spindle = y;
end;

