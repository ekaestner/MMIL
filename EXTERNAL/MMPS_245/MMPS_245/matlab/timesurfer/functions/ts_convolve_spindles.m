function wavecoef = ts_convolve_spindles(data,sfreq,varargin)
%function wavecoef = ts_convolve_spindles(data,sfreq,[options])
%
% Usage:
%  results = ts_convolve_spindles(data,'key1', value1,...);
%
% Required Input:
%   data: data matrix with size = [nchans,ntpoints]
%   sfreq: sampling frequency (in Hz)
%
% Optional parameters for spindle detection:
%  'wavelet_freqs': vector of wavelet frequencies
%    {default = [10:16]}
%  'wavelet_dur': duration of the wavelet (sec)
%    {default = 1}
%  'wavelet_width': wavelet width (sec)
%    {default = 0.6}
%  'wavelet_win': moving average window length (sec) of wavelet convolution
%    {default = 0.3}
%  'offset_flag': [0|1] offset output by subtracting median (for each channel)
%    {default = 0}
%  'norm_flag': [0|1] normalize output by dividing
%    by median absolute deviation (for each channel)
%    {default = 0}
%  'excl_perc': for calculating median and MAD, exclude values higher
%    than indicated percentile
%    {default = 100}
%  'combine_flag': [0|1|2] combine convolved waveforms across frequencies
%    0: return convolution time course for each frequency
%    1: average convolution time courses across frequencies
%    2: maximum convolution time courses across frequencies
%    {default = 1}
%
% Output:
%   wavecoef: matrix of data X wavelet convolution time courses
%     if combine_flag = 1, size of wavecoef is [nchans,ntpoints]
%     if combine_flag = 0, size of wavecoef is [nchans,nfreqs,ntpoints]
%
% Created:  07/24/14 by Don Hagler
% Last mod: 01/28/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavecoef = [];
if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters
parms = check_input(data,sfreq,varargin);

% convolve data with wavelets
wavecoef = convolve_data(parms,data);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(data,sfreq,options)
  parms = mmil_args2parms(options,{...
    'sfreq',sfreq,[0 Inf],...
  ... % spindle detection
    'wavelet_freqs',[10:16],[1,1000],...
    'wavelet_dur',1,[0 Inf],...
    'wavelet_width',0.6,[0.01,100],...
    'wavelet_win',0.3,[0,10],...
    'offset_flag',false,[false true],...
    'norm_flag',false,[false true],...
    'excl_perc',100,[1,100],...
    'combine_flag',1,[0:2],...
  });

  % check data matrix
  if numel(size(data))~=2
    error('data matrix must be 2-dimensional');
  end;
  [parms.nchans,parms.ntpoints] = size(data);
  parms.nfreqs = length(parms.wavelet_freqs);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wavecoef = convolve_data(parms,data);
  if parms.combine_flag
    wavecoef = zeros(parms.nchans,parms.ntpoints);
  else
    wavecoef = zeros(parms.nchans,parms.nfreqs,parms.ntpoints);
  end;
%  fprintf('%s: convolving data with wavelets...\n',mfilename);
  % create wavelets
  wavelets = create_wavelets(parms.wavelet_freqs,parms.wavelet_dur,...
    parms.wavelet_width,parms.sfreq);
  for i=1:parms.nchans
    data_chan = data(i,:);
    % convolve data with wavelets
    nwin = round(parms.wavelet_win*parms.sfreq);
    wavecoef_chan = conv_data(data_chan,wavelets,nwin);
    % combine across frequencies
    switch parms.combine_flag
      case {0,1}
        % average across frequencies
        comb_wavecoef_chan = mean(wavecoef_chan,1);
      case 2
        % maximum across frequencies
        comb_wavecoef_chan = max(wavecoef_chan,[],1);
    end;
    % normalize time courses
    if parms.offset_flag || parms.norm_flag
      if parms.excl_perc<100
        % calculate exclusion cut-off
        excl_thresh = prctile(comb_wavecoef_chan,parms.excl_perc);
        % restrict to values below cut-off
        vals = comb_wavecoef_chan(comb_wavecoef_chan <= excl_thresh);
      else
        vals = comb_wavecoef_chan;
      end;
      % calculate median absolute deviation (mad) and median
      [wmad,wmed] = mmil_wtd_mad(vals,[],2);
      % subtract median
      if parms.offset_flag
        if parms.combine_flag
          comb_wavecoef_chan = comb_wavecoef_chan - wmed;
        else
          wavecoef_chan = bsxfun(@minus,wavecoef_chan,wmed);
        end;
      end;
      % normalize by mad
      if parms.norm_flag
        if parms.combine_flag
          comb_wavecoef_chan = comb_wavecoef_chan/wmad;
        else
          wavecoef_chan = bsxfun(@rdivide,wavecoef_chan,wmad);
        end;
      end;
    end;
    if parms.combine_flag
      wavecoef(i,:) = comb_wavecoef_chan;
    else
      wavecoef(i,:,:) = wavecoef_chan;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wavelets = create_wavelets(freqs,dur,width,sfreq)
  nwavelets = length(freqs);
  wavelet_ntpoints = ceil(dur*sfreq);
  wavelets = zeros(nwavelets,wavelet_ntpoints);
  x = [0:wavelet_ntpoints-1]/sfreq;
  x0 = x(round(length(x)/2));
  for i=1:nwavelets
    freq = freqs(i);
    % make sine function
    %% todo: use complex wavelets?
    y = cos(x.*2*pi*freq);
    % make gaussian function
    g = exp(-(pi*(x-x0).^2)/width^2);
    wavelets(i,:) = y.*g;   
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wavecoef = conv_data(data,wavelets,nwin)
  % pad data with nwin zeros on either end
  ntpoints = length(data);
  data_pad = zeros(1,ntpoints+2*nwin);
  data_pad(nwin+1:nwin+ntpoints) = data;
  nwavelets = size(wavelets,1);
  if nwin>0
    avg_window = tukeywin(nwin)'/nwin;
  end;
  wavecoef = zeros(nwavelets,ntpoints);
  for i=1:nwavelets
    % calculate wavelet coefficients by convolving data with wavelets
    tmp = abs(conv(data_pad,wavelets(i,:),'same'));
    % calculate moving average by convolving with averaging window
    if nwin>0
      tmp = conv(tmp,avg_window,'same');
    end;
    % trim zero-padded regions
    wavecoef(i,:) = tmp(nwin+1:nwin+ntpoints);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

