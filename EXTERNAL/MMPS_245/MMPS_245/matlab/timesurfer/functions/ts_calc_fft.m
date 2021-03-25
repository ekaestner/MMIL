function results = ts_calc_fft(data,sfreq,varargin)
%function results = ts_calc_fft(data,sfreq,[options])
%
% Usage:
%  results = ts_calc_fft(data,sfreq,'key1', value1,...);
%
% Required Input:
%   data: matrix with size = [nchans,ntpoints] or [nepochs,nchans,ntpoints]
%     or cell array with length nepochs,
%       containing matrices of [nchans,ntpoints_var]
%   sfreq: sampling frequency
%
% Optional parameters
%  'mask': vector of 0 and 1 corresponding to time points to exclude from epochs
%    {default = []}
%  'epoch_dur': duration of epochs on which to calculate FFT
%      if Inf, use all time points, excluding masked
%    {default = Inf}
%  'norm_flag': [0|1] normalize FFT power spectra by frequency
%    {default = 0}
%  'norm_exp': exponent applied to frequency when normalizing FFT spectra
%    {default = 1.0}
%
% Output:
%   results: struct containing various spindle metrics
%
% Created:  01/27/16 by Don Hagler
% Last Mod: 02/12/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];

% parse input parameters
[parms,data] = check_input(data,sfreq,varargin);

% initialize results struct
results = init_results(parms,data);

% calculate amplitude and phase from FFT
results = calc_fft(data,parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data] = check_input(data,sfreq,options)
  parms = mmil_args2parms(options,{...
    'sfreq',sfreq,[],...
  ...
    'mask',[],[],...
    'epoch_dur',Inf,[],...
    'norm_flag',false,[false true],...
    'norm_exp',1.0,[0,10],...
  ...
    'fft_freq',[0:0.1:200],[],...
    'min_epoch_frac',1,[],... % skip epochs with any masking if 1
    'min_chan_nepochs',10,[],... % minimum number of epochs or replace with NaNs
  });
  if iscell(data)
    parms.nepochs = length(data);
    parms.nchans = size(data{1},1);
    parms.ntpoints = 1;
  else
    if length(size(data))==3
      [parms.nepochs,parms.nchans,parms.epoch_ntpoints] = size(data);
      parms.ntpoints = parms.nepochs * parms.epoch_ntpoints;
    else
      [parms.nchans,parms.ntpoints] = size(data);
      if isempty(parms.mask)
        parms.mask = ones(1,parms.ntpoints);
      end;
      parms.mask = mmil_rowvec(parms.mask);
      if length(parms.mask) ~= parms.ntpoints
        error('mask does not match tpoints');
      end;
      if isinf(parms.epoch_dur);
        parms.nepochs = 1;
        parms.epoch_ntpoints = parms.ntpoints;
        data = reshape(data,[1,parms.nchans,parms.ntpoints]);
      else
        parms.epoch_ntpoints = parms.epoch_dur * parms.sfreq;
        parms.nepochs = ceil(parms.ntpoints / parms.epoch_ntpoints);
        fprintf('%s: %d epochs\n',mfilename,parms.nepochs);
        ntpoints = parms.nepochs * parms.epoch_ntpoints;
        ntpoints_extra = ntpoints - parms.ntpoints;
        data = cat(2,data,zeros(parms.nchans,ntpoints_extra));
        parms.mask = cat(2,parms.mask,zeros(1,ntpoints_extra));
        data = reshape(data,[parms.nchans,parms.epoch_ntpoints,parms.nepochs]);
        data = permute(data,[3,1,2]);
        parms.mask = reshape(parms.mask,[1,parms.epoch_ntpoints,parms.nepochs]);
        parms.mask = permute(parms.mask,[3,1,2]);
        parms.ntpoints = ntpoints;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms,spindles);
  results = [];
  fnames = {'sfreq','nepochs','nchans','ntpoints'};
  for i=1:length(fnames)
    results.(fnames{i}) = parms.(fnames{i});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_fft(data,parms,results)
  % initialize output values
  results.epochs = [];
  % loop over epochs, calculate fft for each channel
  results.nskips = 0;
  for i=1:parms.nepochs
    if ~iscell(data)
      tdata = data(i,:,:);
      ntpoints = size(tdata,3);
      if ~isempty(parms.mask)
        tmask = parms.mask(i,:,:);
      else
        tmask = ones(1,ntpoints);
      end;
    else
      tdata = data{i};
      ntpoints = size(tdata,2);
      tmask = ones(1,ntpoints);
    end;
    tdata = reshape(tdata,[results.nchans,ntpoints]);
    tmask = reshape(tmask,[1,ntpoints]);
    if length(find(tmask)) / ntpoints < parms.min_epoch_frac
      results.epochs(i).chans = [];
      results.nskips = results.nskips + 1;
      fprintf('%s: skipping epoch %d because of masking...\n',mfilename,i);
      continue;
    end;
    % exclude chans if there are nans in this epoch
    chans = find(~isnan(sum(tdata,2)));
    tdata = tdata(chans,:);
    % set the number of zero-padded time points for FFT
    nfft = 2^nextpow2(ntpoints);
    % calculate Fourier transform
    fft_cx = fft(tdata,nfft,2)/ntpoints;
    fft_cx = fft_cx(:,1:nfft/2+1);
    % calculate frequencies
    freq = parms.sfreq/2*linspace(0,1,nfft/2+1);
    % calculate power
    amp = 2*abs(fft_cx);
    % scale power by frequency
    if parms.norm_flag
      amp = bsxfun(@times,amp,freq.^parms.norm_exp);
    end;
    % calculate phase (in cycles, not radians)
    phase = angle(fft_cx)/(2*pi);
    % save results
    results.epochs(i).chans = chans;
    results.epochs(i).freq = freq;
    results.epochs(i).amp = amp;
    results.epochs(i).phase = phase;
  end;
  fprintf('%s: %d of %d epochs skipped\n',...
    mfilename,results.nskips,parms.nepochs);
  % calculate average Fourier spectra for each channel
  nfreq = length(parms.fft_freq);  
  epochs = results.epochs;
  power_sum = zeros(results.nchans,nfreq);
  power_ssq = zeros(results.nchans,nfreq);
  power_n = zeros(results.nchans,1);
  for i=1:results.nepochs
    chans = epochs(i).chans;
    if isempty(chans), continue; end;
    tmp_freq = epochs(i).freq;
    tmp_amp = epochs(i).amp;
    amp = spline(tmp_freq,tmp_amp,parms.fft_freq);
    power_sum(chans,:) = power_sum(chans,:) + amp;
    power_ssq(chans,:) = power_ssq(chans,:) + amp.^2;
    power_n(chans) = power_n(chans) + 1;
  end;
  n = repmat(power_n,[1,nfreq]);
  results.spectra_freq = parms.fft_freq;
  results.spectra_mean = power_sum ./ n;
  results.spectra_std = sqrt(n.*power_ssq - power_sum.^2)./(n.*(n-1));
  results.spectra_n = power_n;
  % set spectra for channels with too few epochs to nan
  ind_excl = find(power_n<parms.min_chan_nepochs);
  spectra_mean(ind_excl,:) = nan;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

