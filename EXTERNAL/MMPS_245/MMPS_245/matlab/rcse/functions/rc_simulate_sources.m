function waveforms = rc_simulate_sources(num_areas,sfreq,tmin,tmax,randseed);
%function waveforms = rc_simulate_sources([num_areas],[sfreq],[tmin],[tmax],[randseed]);
%
% Early Mod: 04/08/09 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

waveforms = [];

if ~exist('num_areas','var') || isempty(num_areas), num_areas = 4; end;
if ~exist('sfreq','var') || isempty(sfreq), sfreq = 1000; end;
if ~exist('tmin','var') || isempty(tmin), tmin = -0.1; end;
tmin = tmin*1000;
if ~exist('tmax','var') || isempty(tmax), tmax = 0.4; end;
tmax = tmax*1000;
if ~exist('randseed','var') || isempty(randseed), randseed = 0; end;

rand('state',randseed);

% these could be parameters
stim_delay = 32.5;
min_amp = 0.1;
max_amp = 0.2;
min_uamp = 0.2;
max_uamp = 0.5;
min_delay = 40;
max_delay = 130;
if num_areas>1
  step_delay = (max_delay-min_delay)/(num_areas-1);
else
  step_delay = 0;
end;
min_udelay = 10;
max_udelay = 40;
min_rise = 10;
max_rise = 30;
min_urise = 10;
max_urise = 30;
min_flat = 5;
max_flat = 15;
min_uflat = 5;
max_uflat = 20;
min_fall = 10;
max_fall = 30;
min_ufall = 30;
max_ufall = 80;

% set waveform properties for each area
amp = min_amp+rand(1,num_areas)*(max_amp-min_amp);
uamp = min_uamp+rand(1,num_areas)*(max_uamp-min_uamp);
if step_delay
  delay = [min_delay:step_delay:max_delay] - tmin + stim_delay;
else
  delay = min_delay - tmin + stim_delay;
end;
udelay = delay + min_udelay + rand(1,num_areas)*(max_udelay-min_udelay);
rise = min_rise+rand(1,num_areas)*(max_rise-min_rise);
urise = min_urise+rand(1,num_areas)*(max_urise-min_urise);
flat = min_flat+rand(1,num_areas)*(max_flat-min_flat);
uflat = min_uflat+rand(1,num_areas)*(max_uflat-min_uflat);
fall = min_fall+rand(1,num_areas)*(max_fall-min_fall);
ufall = min_ufall+rand(1,num_areas)*(max_ufall-min_ufall);

% convert from msec to samples
nsamps = round((tmax-tmin)*sfreq/1000)+1;
tstep = (tmax-tmin)/(nsamps-1);
delay = round(delay*sfreq/1000);
udelay = round(udelay*sfreq/1000);
rise = round(rise*sfreq/1000);
urise = round(urise*sfreq/1000);
flat = round(flat*sfreq/1000);
uflat = round(uflat*sfreq/1000);
fall = round(fall*sfreq/1000);
ufall = round(ufall*sfreq/1000);


waveforms = zeros(num_areas,nsamps);
for i=1:num_areas
  peak = rc_respfunc(amp(i),delay(i),rise(i),flat(i),fall(i),nsamps);
  upeak = rc_respfunc(uamp(i),udelay(i),urise(i),uflat(i),ufall(i),nsamps);
  waveforms(i,:) = peak(1:nsamps)-upeak(1:nsamps);
end;

