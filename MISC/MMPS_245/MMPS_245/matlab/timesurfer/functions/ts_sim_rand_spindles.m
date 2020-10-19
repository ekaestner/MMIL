function spindles = ts_sim_rand_spindles(varargin)
%function spindles = ts_sim_rand_spindles([options])
%
% Optional parameters:
%   'nchans': number of simulated channels
%     {default = 1}
%   'nepochs': number of simulated epochs
%     {default = 1}
%   'time_dur': duration of simulation time window (msec)
%     {default = 1000}
%   'sf': sampling frequency (Hz)
%     {default = 500}
%   'freq_mean': mean spindle frequency
%     {default = 13}
%   'freq_std': standard deviation of spindle frequency
%     {default = 1}
%   'phase_mean': mean spindle phase (cycles)
%     {default = 0}
%   'phase_std': st.dev. of spindle phase (cycles)
%     {default = 1}
%   'phase_gauss_flag': [0|1] use Gaussian distribution for spindle phases
%     otherwise use uniform distribution with range of phase_std
%     {default = 0}
%   'amp_mean': mean spindle amplitude (nA m)
%     {default = 1}
%   'amp_std': st.dev. of spindle amplitude (nA m)
%     {default = 0}
%   't0_min': minimum spindle start time (msec)
%     {default = 100}
%   't0_mean': mean spindle start time (msec)
%     {default = 250}
%   't0_std': st.dev. of spindle start time (msec)
%     {default = 50}
%   'td_min': minimum spindle duration (msec)
%     {default = 200}
%   'td_mean': mean spindle duration (msec)
%     {default = 500}
%   'td_std': st.dev. of spindle duration (msec)
%     {default = 50}
%   'source_noise': st.dev. of Gaussian noise added to each dipole (nA m)
%     {default = 0.1}
%
% Output:
%   spindles: matrix with size = [nchans,nepochs*ntpoints]
%     with ntpoints = round(time_dur*1000 / sf)
%
% Created:  04/23/14 by Don Hagler
% Last Mod: 04/23/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);

% initialize spindles matrix
spindles = zeros(parms.nchans,parms.ntpoints_total);

% loop over channels and epochs
for i=1:parms.nchans
  for j=1:parms.nepochs
    % simulate a spindle
    spindle = sim_spindle(parms);
    % insert spindle into spindles
    t0 = 1 + (j-1)*parms.ntpoints;
    t1 = t0 + parms.ntpoints - 1;
    spindles(i,t0:t1) = spindle;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'nchans',1,[1,1e4],...
    'nepochs',1,[1,1e4],...
    ... % time window
    'time_dur',1000,[1e2,1e10],...
    'sf',500,[1e2,1e4],...
    ... % spindle characteristics
    'freq_mean',13,[0,1e2],...
    'freq_std',1,[0,1e2],...
    'phase_mean',0,[0,1],...
    'phase_std',1,[0,1],... 
    'phase_gauss_flag',true,[false true],...
    'amp_mean',1,[-100,100],...
    'amp_std',0,[0,100],...
    't0_min',100,[-1e4,1e4],...
    't0_mean',250,[-1e4,1e4],...
    't0_std',50,[0,1e4],...
    'td_min',200,[-1e4,1e4],...
    'td_mean',500,[-1e4,1e4],...
    'td_std',50,[0,1e4],...
    'source_noise',0.1,[0,1e100],...
  });

  parms.ntpoints = round(parms.time_dur * parms.sf / 1000);
  parms.ntpoints_total = parms.nepochs*parms.ntpoints;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spindle = sim_spindle(parms)
  freq = parms.freq_mean + parms.freq_std*randn;
  if parms.phase_gauss_flag
    phas = parms.phase_mean + parms.phase_std*randn;
  else
    phas = parms.phase_mean + parms.phase_std*rand;
  end;
  amp = parms.amp_mean + parms.amp_std*randn;
  t0 = max(parms.t0_min,parms.t0_mean + parms.t0_std*randn);
  td = max(parms.td_min,parms.td_mean + parms.td_std*randn);
  spindle = ts_sim_spindle(parms.ntpoints,freq,phas,t0,td,amp,parms.sf);
  % add source noise
  spindle = spindle + parms.source_noise*randn(size(spindle));
return;

