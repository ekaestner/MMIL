function regstruct = rc_landmark_reg(data,time,varargin)
%function regstruct = rc_landmark_reg(data,time,[options])
%
% Required Input:
%   data: time series data; with size [ntpoints,nconds]
%   time: time vector; with size [ntpoitns,1]
%
% Optional Parameters:
%   'padflag': [0|1|2] whether to zero-pad
%     0: no zero-padding
%     1: zero-padding at begining of time series
%     2: zero-padding at begining and end of time series
%     {default = 1}
%   'npad': number of zero-padded timepoints
%     {default = 50}
%   'ntgrid': number of evenly spaced nodes in time
%     {default = 25}
%   'extra_nodes_flag': [0|1] add extra nodes at large extrema
%     of average time course
%     {default = 1}
%   'spline_order': fda spline order
%     {default = 4}
%   'smooth_flag': [0|1] whether to do derivative constrained smoothing
%     {default = 1}
%   'smooth_order': derivative constrained smoothing order
%     {default = 2}
%   'smooth_lambda': regularization term for derivative constrained smoothing
%     {default = 1e-1}
%   'zero_landmark_flag': [0|1] include landmark at time = 0
%     {default = 0}
%   'delta': peak detection minimum deflection
%     {default = 0.1}
%   'peak_tmin': minimum time for landmark peak finding
%     {default = 50}
%   'peak_tmin': maximum time for landmark peak finding
%     {default = 160}
%   'peak_pol': peak polarity for landmark peak finding
%     {default = -1}
%
% Output:
%   regstruct: struct containing these fields:
%     ntpoints: number of time points
%     nconds: number of conditions
%     data: input data matrix; size = [ntpoints,nconds]
%     time: input time vector; size = [ntpoints,1]
%     data_sp: spline fit of data matrix; size = [ntpoints,nconds]
%     data_reg: registered data matrix; size = [ntpoints,nconds]
%     warpfun: time warping function; size = [ntpoints,1]
%
% Created:  03/24/16 by Don Hagler
% Last mod: 05/27/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regstruct = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'padflag',1,[0,1,2],...
  'npad',50,[],...
  'ntgrid',25,[],...
  'extra_nodes_flag',true,[false true],...
  'spline_order',4,[2:10],...
  'smooth_flag',1,[],...
  'smooth_order',2,[],...
  'smooth_lambda',1e-1,[],...
  'reg_lambda',1e-4,[],...
  'zero_landmark_flag',false,[false true],...
  'delta',0.1,[],...
  'peak_tmin',50,[],...
  'peak_tmax',160,[],...
  'peak_pol',-1,[-1,1],...
  ...
  'quiet_landmark_flag',true,[false true],...
  'monotone_reg_flag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure ndims of data is 2
regstruct.data = data;
if ndims(regstruct.data) ~= 2
  error('data must be a 2-dimensional matrix');
end;

% make sure time is a vector
regstruct.time = mmil_colvec(time);

% get dimensions of input data
[regstruct.ntpoints,regstruct.nconds] = size(data);

% check for mismatch between data and time
if regstruct.ntpoints ~= length(regstruct.time)
  error('number of columns in data matrix does not match number of elements in time vector');
end;

% check lengths of peak_tmin, peak_tmax, and peak_pol
parms.npeaks = length(parms.peak_tmin);
if length(parms.peak_tmax) ~= parms.npeaks
  error('mismatch between peak_tmin and peak_tmax');
end;
if length(parms.peak_pol) ~= parms.npeaks
  error('mismatch between peak_tmin and peak_pol');
end;
parms.range_tmin = min(parms.peak_tmin);
parms.range_tmax = max(parms.peak_tmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% zero pad data and time
if parms.padflag
  npad = parms.npad;
  ntpoints = regstruct.ntpoints + parms.padflag*npad;
  data = zeros(ntpoints,regstruct.nconds);
  data(npad+1:npad+regstruct.ntpoints,:) = regstruct.data;
  sfreq = 1000/(time(2) - time(1));
  t = 1000*[1:ntpoints]'/sfreq;
  t0 = time(1);
  t1 = t(npad+1);
  time = t - t1 + t0;
else
  ntpoints = regstruct.ntpoints;  
  time = regstruct.time;
  npad = 0;
end;

% calculate normalized average across conditions
data_avg = mean(data,2);
data_avg = data_avg/max(abs(data_avg));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define nodes for spline fitting
regstruct.bspline_time = time;
regstruct.bspline_npad = npad;
regstruct.bspline_padflag = parms.padflag;
regstruct.bspline_nodes = round(linspace(1,ntpoints,parms.ntgrid));
if parms.extra_nodes_flag
  [maxtab,mintab] = mmil_peakdet(data_avg,parms.delta);
  idx = union(maxtab(:,1),mintab(:,1));
  regstruct.bspline_nodes = union(regstruct.bspline_nodes,idx);
end;

% define basis functions for spline fitting
breaks = time(regstruct.bspline_nodes);
basis_range = [min(time) max(time)];
nbasis = length(breaks) + parms.spline_order - 2;
regstruct.bspline_obj =...
  create_bspline_basis(basis_range,nbasis,parms.spline_order,breaks);

% convert data to functional data object
data_fd = data2fd(time,data,regstruct.bspline_obj);

% apply derivative constrained smoothing
if parms.smooth_flag
  regstruct.smooth_order = parms.smooth_order;
  regstruct.smooth_obj = fdPar(data_fd,parms.smooth_order,parms.smooth_lambda);
  data_fd = smooth_fd(data_fd,regstruct.smooth_obj);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set number of landmarks
nlm = parms.npeaks;
if parms.zero_landmark_flag
  nlm = nlm + 1;
  [tmp,idx] = min(abs(time));
  time_zero = time(idx);
end;

% set landmarks for reference
regstruct.landmarks_avg = zeros(1,nlm);
m = 1;
if parms.zero_landmark_flag
  regstruct.landmarks_avg(1,m) = time_zero;
  m = m + 1;  
end;
% find minima and maxima in average
[maxtab,mintab] = mmil_peakdet(data_avg,parms.delta,time);
% choose minima and maxima between range_tmin and range_tmax
tmp_minima = get_minmax([parms.range_tmin,parms.range_tmax],mintab);
tmp_maxima = get_minmax([parms.range_tmin,parms.range_tmax],maxtab);
tmp_minima = get_deflection(tmp_minima,tmp_maxima);
tmp_maxima = get_deflection(tmp_maxima,tmp_minima);
latency = cat(1,tmp_minima.latency',tmp_maxima.latency');
amplitude = cat(1,tmp_minima.amplitude',tmp_maxima.amplitude');
deflection = cat(1,tmp_minima.deflection',tmp_maxima.deflection');
for i=1:parms.npeaks
  tmin = parms.peak_tmin(i);
  tmax = parms.peak_tmax(i);
  pol = parms.peak_pol(i);
  % identify peaks between tmin and tmax
  idm = find(latency>=tmin & latency<=tmax);
  if isempty(idm)
    lat = mean([tmin,tmax]);
    fprintf('%s: WARNING: no peaks between %0.1f and %0.1f\n',...
      mfilename,tmin,tmax);
  else
    % choose biggest peak
    amp = amplitude(idm);
    [tmp,idq] = max(pol*amp);
    if sign(amp) ~= sign(pol)
      lat = mean([tmin,tmax]);
      fprintf('%s: WARNING: no peaks with correct polarity between %0.1f and %0.1f\n',...
        mfilename,tmin,tmax);
    else
      lat = latency(idm(idq));
    end;
  end;
  regstruct.landmarks_avg(1,m) = lat;
  m = m + 1;
end;

% set landmarks for each waveform
regstruct.landmarks = zeros(regstruct.nconds,nlm); 
for s=1:regstruct.nconds
  m = 1;
  if parms.zero_landmark_flag
    regstruct.landmarks(s,m) = time_zero;
    m = m + 1;
  end;
  v = data(:,s);
  [maxtab,mintab] = mmil_peakdet(v,parms.delta,time);
  % choose minima and maxima between range_tmin and range_tmax
  tmp_minima = get_minmax([parms.range_tmin,parms.range_tmax],mintab);
  tmp_maxima = get_minmax([parms.range_tmin,parms.range_tmax],maxtab);
  tmp_minima = get_deflection(tmp_minima,tmp_maxima);
  tmp_maxima = get_deflection(tmp_maxima,tmp_minima);
  latency = cat(1,tmp_minima.latency',tmp_maxima.latency');
  amplitude = cat(1,tmp_minima.amplitude',tmp_maxima.amplitude');
  deflection = cat(1,tmp_minima.deflection',tmp_maxima.deflection');
  for i=1:parms.npeaks
    tmin = parms.peak_tmin(i);
    tmax = parms.peak_tmax(i);
    pol = parms.peak_pol(i);
    % identify peaks between tmin and tmax
    idm = find(latency>=tmin & latency<=tmax);
    if isempty(idm)
      lat = mean([tmin,tmax]);
    else
      % choose biggest peak
      amp = amplitude(idm);
      [tmp,idq] = max(pol*amp);
      if sign(amp) ~= sign(pol)
        lat = mean([tmin,tmax]);
      else
        lat = latency(idm(idq));
      end;
    end;
    regstruct.landmarks(s,m) = lat;
    m = m + 1;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% landmark registration
if parms.monotone_reg_flag
  Wfd0 = fd(zeros(nbasis,1),regstruct.bspline_obj);
  WfdPar = fdPar(Wfd0,1,parms.reg_lambda);
  mono_flag = 1;
else
  WfdPar = [];
  mono_flag = 0;
end;
if parms.quiet_landmark_flag
  [T,data_reg_fd,warpfun_fd] = evalc('landmarkreg(data_fd,regstruct.landmarks,regstruct.landmarks_avg,WfdPar,mono_flag)');
else
  [data_reg_fd,warpfun_fd] = landmarkreg(data_fd,regstruct.landmarks,...
    regstruct.landmarks_avg,WfdPar,mono_flag);
end;

% extract data matrices from fda objects
data_sp = eval_fd(time,data_fd); % spline-fitted original data
data_reg = eval_fd(time,data_reg_fd); % registered values
warpfun = eval_fd(time,warpfun_fd); % warping function

% add results to regstruct
regstruct.data_reg_fd = data_reg_fd;
regstruct.warpfun_fd = warpfun_fd;
regstruct.data_sp = data_sp(npad+1:npad+regstruct.ntpoints,:);
regstruct.data_reg = data_reg(npad+1:npad+regstruct.ntpoints,:);
regstruct.warpfun = warpfun(npad+1:npad+regstruct.ntpoints,:);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax = get_minmax(peak_range,mtab)
  minmax = struct('latency',[],'amplitude',[],'deflection',[]);
  j = 1;
  for i=1:size(mtab,1)
    latency = mtab(i,1);
    amplitude = mtab(i,2);
    if latency<peak_range(1) ||...
       latency>peak_range(2)
      continue;
    end;
    minmax.latency(j) = latency;
    minmax.amplitude(j) = amplitude;
    j = j + 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax1 = get_deflection(minmax1,minmax2)
  % calculate deflection for minima or maxima
  for i=1:length(minmax1.latency)
    latency = minmax1.latency(i);
    % find preceding minimum
    tmp = latency - minmax2.latency;
    tmp(tmp<0)=Inf;
    [tmp,ind1] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp1 = 0;
    else
      amp1 = minmax2.amplitude(ind1);
    end;

    % find subsequent maxima or minima
    tmp = minmax2.latency - latency;
    tmp(tmp<0)=Inf;
    [tmp,ind2] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp2 = 0;
    else
      amp2 = minmax2.amplitude(ind2);
    end;

    minmax1.deflection(i) =...
      abs(minmax1.amplitude(i) - mean([amp1 amp2]));
  end;
return;

