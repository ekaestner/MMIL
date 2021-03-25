function wform = rc_model_wform(varargin)
%function wform = rc_model_wform([options])
%
% Optional Input:
%   'sfreq': sampling frequency
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 400}
%   'time': time vector (if supplied, sfreq, t0, and t1 are ignored)
%     {default = []}
%
% Optional Waveform Parameters:
%   'polarity': waveform polarity (1 or -1)
%     {default = -1}
%   'latency': peak latency (msec)
%     {default = 60}
%   'amplitude': vector of component amplitudes
%     {default = 10}
%   'rise_tc': rise time constant
%     {default = 2}
%   'fall_tc': fall time constant
%     {default = 4}
%
% Created:  05/10/11 by Don Hagler
% Last Mod: 05/16/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
...
  'polarity',-1,[-1,1,1],...
  'latency',60,[0,Inf],...
  'amplitude',10,[0,Inf],...  
  'rise_tc',2,[],...
  'fall_tc',4,[],...
...
  'wfit_func','sigmoid',{'sigmoid','gamma'},...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wform = [];

if isempty(parms.time)
  sd = 1000/parms.sfreq;
  parms.time = [parms.t0:sd:parms.t1];
end;
ntpoints = length(parms.time); 

switch parms.wfit_func
  case 'sigmoid'
    % model rise and fall with sigmoid functions
    rise = sigmoid((parms.time' - parms.latency)/parms.rise_tc);
    fall = 1-sigmoid((parms.time' - parms.latency)/parms.fall_tc);
    wform = rise .* fall;
  case 'gamma'
    wform = gamma(parms.time',parms.latency,parms.fall_tc);
end;

wform = wform*parms.polarity*parms.amplitude/max(wform);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wf = sigmoid(x)
  wf = 1./(1 + exp(-x));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wf = gamma(t,t0,tc)
  wf = zeros(size(t));
  ind = find(t>=t0);
  wf(ind) = (t(ind) - t0) .* exp(-(t(ind) - t0)/tc);
return;

