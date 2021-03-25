function results = rc_fit_wform_components(wform,varargin)
%function results = rc_fit_wform_components(wform,[options])
%
% Required Input:
%   wform: 1D waveform vector
%
% Optional Input:
%   'niters': number of iterations for random search
%     {default = 500}
%   'stepsize': size of random search step relative to range of bounds
%     {default = 0.05}
%   'rand_init_flag': [0|1] use random starting values for parameter estimates
%     otherwise use middle of bounds
%   'ncomponents': number of components to model
%     {default = 1}
%   'polarity': waveform component polarity (may be vector with ncomponents)
%     {default = -1}
%   'latency': initial estimate(s) for component latency (msec)
%     if ncomponents>1, should be vector with ncomponents
%     {default = []}
%   'amplitude': initial estimate(s) for component amplitude
%     if ncomponents>1, should be vector with ncomponents
%     {default = []}
%   'rise_tc':initial estimate(s) for exponential rise time constant (msec)
%     if ncomponents>1, should be vector with ncomponents
%     {default = []}
%   'fall_tc': initial estimate(s) for exponential fall time constant (msec)
%     if ncomponents>1, should be vector with ncomponents
%     {default = []}
%   'latency_bounds': vector of lower and upper bounds for latency parameter
%     May be matrix with size = [ncomponents,2]
%     {default = [40,120]}
%   'amplitude_bounds': lower and upper bounds for amplitude parameter
%     May be matrix with size = [ncomponents,2]
%     {default = [0.5,25]}
%   'rise_tc_bounds': component exponential rise time constant (msec)
%     May be matrix with size = [ncomponents,2]
%     {default = [1,20]}
%   'fall_tc_bounds': component exponential fall time constant (msec)
%     May be matrix with size = [ncomponents,2]
%     {default = [1,20]}
%   'sfreq': sampling frequency; used to calculate time vector
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 400}
%   'time': time vector (msec)
%     if supplied, sfreq, t0, and t1 are ignored
%     length of time vector must match length of input wform
%     {default = []}
%
% Created:  05/10/11 by Don Hagler
% Last Mod: 05/17/11 by Don Hagler
%

%% todo: allow restricted time range for fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'search_type','rand',{'rand','fmincon'},...
...
  'niters',500,[0,1e10],...
  'stepsize',0.05,[],...
  'rand_init_flag',true,[false true],...
...
  'ncomponents',1,[],...
  'polarity',-1,[],...
...
  'latency',[],[],...
  'amplitude',[],[],...
  'rise_tc',[],[],...
  'fall_tc',[],[],...
...
  'latency_bounds',[40,120],[],...
  'amplitude_bounds',[0.5,25],[],...
  'rise_tc_bounds',[1,20],[],...
  'fall_tc_bounds',[1,20],[],...
...
  'wfit_func','sigmoid',{'sigmoid','gamma'},...
  'single_tc_flag',false,[false true],...
...
  'delay_sf',2,[],...
...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
...
  'DiffMinChange',0.001,[],...
  'TolFun',1e-5,[],...
  'TolX',1e-6,[],...
};
results = [];
parms = check_input(wform,varargin,parms_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

betas = init_betas(parms);

switch parms.search_type
  case 'rand'
    [betas,min_err] = rand_search(betas,parms);
  case 'fmincon'
    [betas,min_err] = fmincon_search(betas,parms);
end;

results = compile_results(betas,min_err,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wform,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.wform = wform;
  parms.beta_names = {'latency','amplitude','rise_tc','fall_tc'};
  parms.nbeta = length(parms.beta_names);
  parms.wfsize = size(parms.wform);
  if length(parms.wfsize)>2 || all(parms.wfsize~=1)
    error('wform must be a vector');
  end;
  parms.ntpoints = length(parms.wform);
  if isempty(parms.time)
    sd = 1000/parms.sfreq;
    parms.time = [parms.t0:sd:parms.t1];
  end;
  if length(parms.time) ~= parms.ntpoints
    error('number of elements in time (%d) does not match wform (%d)',...
      length(parms.time),parms.ntpoints);
  end;
  for b=1:parms.nbeta
    val = parms.(parms.beta_names{b});
    if ~isempty(val) && numel(val)~=parms.ncomponents
      error('initial estimate for %s must be vector with %d values',...
        parms.beta_names{b},parms.ncomponents);
    end;

    bname = [parms.beta_names{b} '_bounds'];
    bounds = parms.(bname);
    bsz = size(bounds);
    if length(bsz)>2
      error('%s has size [%s]',bname,sprintf('%d ',bsz));
    end;
    if numel(bounds)==2
      bounds = ones(parms.ncomponents,1)*reshape(bounds,[1,2]);
    elseif any(bsz~=[parms.ncomponents,2])
      error('size of %s is [%s] but should be [%d 2]',...
        bname,sprintf('%d ',bsz),parms.ncomponents);
    end;
    parms.(bname) = bounds;
  end;

  if strcmp(parms.wfit_func,'gamma') || parms.single_tc_flag
    parms.beta_names = parms.beta_names([1,2,4]);
    parms.nbeta = length(parms.beta_names);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wform = model_wform(betas,parms)
  wform = single(0);
  for c=1:parms.ncomponents
    tmp_parms = [];
    tmp_parms.time = parms.time;
    tmp_parms.polarity = parms.polarity(c);
    for b=1:parms.nbeta
      tmp_parms.(parms.beta_names{b}) = betas(c,b);
    end;
    if parms.single_tc_flag
      tmp_parms.rise_tc = tmp_parms.fall_tc;
    end;
    tmp_parms.wfit_func = parms.wfit_func;

    args = mmil_parms2args(tmp_parms);
    wform = wform + rc_model_wform(args{:});

%    wform = wform + rc_model_wform(...
%      'time',parms.time,...
%      'polarity',parms.polarity(c),...
%      'latency',betas(c,1),...
%      'amplitude',betas(c,2),...
%      'rise_tc',betas(c,3),...
%      'fall_tc',betas(c,4));

  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function betas = init_betas(parms)
  betas = zeros(parms.ncomponents,parms.nbeta);
  for c=1:parms.ncomponents
    for b=1:parms.nbeta
      vals = parms.(parms.beta_names{b});
      bounds = parms.([parms.beta_names{b} '_bounds']);
      % require that latency bounds for c>1 are later
      if c>1 && strcmp(parms.beta_names{b},'latency')
        % set lower bound to be latency for earlier component
        %   plus offset time constant (fall_tc)
        lower_bound = betas(c-1,1) + parms.delay_sf*betas(c-1,end);
        bounds(c,1) = max(bounds(c,1),lower_bound);
      end;
      if ~isempty(vals)
        betas(c,b) = check_val(vals(c),bounds(c,:));
      elseif parms.rand_init_flag
        betas(c,b) = rand_range(bounds(c,:));
      else
        betas(c,b) = mean(bounds(c,:));
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function betas = rand_step_betas(betas,parms)
  for c=1:parms.ncomponents
    for b=1:parms.nbeta
      bounds = parms.([parms.beta_names{b} '_bounds']);
      % require that latency bounds for c>1 are later
      if c>1 && strcmp(parms.beta_names{b},'latency')
        % set lower bound to be latency for earlier component
        %   plus offset time constant (fall_tc)
        lower_bound = betas(c-1,1) + parms.delay_sf*betas(c-1,4);
        bounds(c,1) = max(bounds(c,1),lower_bound);
      end;
      betas(c,b) = rand_step(betas(c,b),bounds(c,:),parms.stepsize);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betas,min_err] = rand_search(betas,parms)
  min_err = 1e10;
  for i=1:parms.niters
    % randomly adjust parameter estimates
    tmp_betas = rand_step_betas(betas,parms);

    % calculate error
    err = calc_err(tmp_betas,parms);

    % update best fit
    if err < min_err
      min_err = err;
      betas = tmp_betas;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betas,min_err] = fmincon_search(betas,parms)
  min_err = [];

  try % for R2009b
    fit_options = optimset('Display','iter',...
      'Algorithm','active-set' ,...
      'MaxFunEvals',Inf,...
      'MaxIter',parms.niters,...
      'DiffMinChange',parms.DiffMinChange,...
      'TolFun',parms.TolFun,...
      'TolX',parms.TolX);
  catch % for R2007a
    fit_options = optimset('Display','iter',...
      'LargeScale','off',...
      'MaxFunEvals',Inf,...
      'MaxIter',parms.niters,...
      'DiffMinChange',parms.DiffMinChange,...
      'TolFun',parms.TolFun,...
      'TolX',parms.TolX);
  end;
  
  betas = reshape(betas,[parms.ncomponents*parms.nbeta,1]);
  [lbounds,ubounds] = reshape_bounds(parms);

  [betas,min_err,exitflag,output] = fmincon(@(b) calc_err(b,parms),betas,...
    [],[],[],[],lbounds,ubounds,[],fit_options);

  if exitflag~=1
    fprintf('%s: WARNING: fmincon returned with exitflag %d\n',...
      mfilename,exitflag);
  end;

  betas = reshape(betas,[parms.ncomponents,parms.nbeta]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = calc_err(betas,parms)
  % reshape betas if necessary
  if numel(betas)==length(betas)
    betas = reshape(betas,[parms.ncomponents,parms.nbeta]);
  end;
  % calculate waveform with new parameters
  wform = model_wform(betas,parms);
  % calculate error
  err_vec = (parms.wform - wform).^2;
  err_norm = err_vec/max(parms.wform.^2);
  err = 100*sum(err_norm)/length(err_norm);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lbounds,ubounds] = reshape_bounds(parms)
  total_nbeta = parms.ncomponents*parms.nbeta;
  lbounds = zeros(total_nbeta,1);
  ubounds = zeros(total_nbeta,1);
  j = 1;
  for b=1:parms.nbeta
    vals = parms.(parms.beta_names{b});
    bounds = parms.([parms.beta_names{b} '_bounds']);
    for c=1:parms.ncomponents
      lbounds(j) = bounds(c,1);
      ubounds(j) = bounds(c,2);
      j = j + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(betas,min_err,parms)
  results = [];
  results.beta_names = parms.beta_names;
  results.nbeta = parms.nbeta;
  results.ncomponents = parms.ncomponents;
  results.polarity = parms.polarity;
  for b=1:parms.nbeta
    results.(parms.beta_names{b}) = betas(:,b);
  end;
  results.time = parms.time;
  results.wform = parms.wform;
  results.wform_fit = model_wform(betas,parms);
  results.min_err = min_err;
  results.betas = betas;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = rand_range(bounds)
  val = bounds(1) + randn*range(bounds)/2;
  val = check_val(val,bounds);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = rand_step(val,bounds,stepsize)
  val = val + randn*range(bounds)*stepsize;
  val = check_val(val,bounds);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = check_val(val,bounds);
  val = max(val,bounds(1));
  val = min(val,bounds(2));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
