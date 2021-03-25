function [data_env,data_filt] = ts_calc_fband_envelope(data,sfreq,varargin)
%function [data_env,data_filt] = ts_calc_fband_envelope(data,sfreq,[options])
%
% Usage:
%  [data_env,data_filt] = ts_calc_fband_envelope(data,sfreq,'key1', value1,...);
%
% Required Input:
%   data: data matrix with size = [nchans,ntpoints]
%   sfreq: sampling frequency (in Hz)
%
% Optional parameters:
%  'smooth_win': moving average window length (sec)
%     for smoothing band-pass filtered envelopes
%    {default = 0.3}
%  'fband': bounds of frequency range (Hz)
%    {default = [10,16]}
%  'filt_tfrac': width of filter transition bands relative to cut-off frequency
%    {default = 0.3}
%
% Created:  06/18/15 by Don Hagler
% Last mod: 06/18/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_filt = []; data_env = [];
if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters
parms = check_input(data,sfreq,varargin);

% apply band-pass filter to data
if parms.verbose
  fprintf('%s: filtering data...\n',mfilename);
end;
data_filt = filter_data(parms,data,parms.fband);

% calculate envelope
if parms.verbose
  fprintf('%s: calculating amplitude envelope...\n',mfilename);
end;
data_env = calc_envelope(parms,data_filt);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(data,sfreq,options)
  parms = mmil_args2parms(options,{...
    'sfreq',sfreq,[0 Inf],...
  ...
    'smooth_win',0.3,[0,10],...
    'fband',[10,16],[0.5,100],...
    'filt_tfrac',0.3,[0.01,1],...
  ...
    'verbose',false,[false true],...
  });

  % check data matrix
  if numel(size(data))~=2
    error('data matrix must be 2-dimensional');
  end;
  [parms.nchans,parms.ntpoints] = size(data);

  % check filtering band
  if numel(parms.fband)~=2
    error('fband must have 2 elements');
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = filter_data(parms,data,cf_vec)
  tb_vec = parms.filt_tfrac*cf_vec;
  for i=1:parms.nchans
    tmp_data = data(i,:)';
    % remove mean in case there is a problem with detrending
    tmp_data = detrend(tmp_data,'constant');
    % try to detrend (may fail in some cases)
    warning('off');
    tmp_data = detrend(tmp_data);
    warning('on');
    tmp_data = ts_freq_filt(tmp_data,...
      parms.sfreq,cf_vec,tb_vec,'bandpass')';
    data(i,:) = tmp_data';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = calc_envelope(parms,data)
  data = smooth_data(parms,abs(data),parms.smooth_win);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = smooth_data(parms,data,win)
  if win==0, return; end;
  ntpoints = size(data,2);
  nwin = round(win*parms.sfreq);
  avg_window = tukeywin(nwin)'/nwin;
  for i=1:parms.nchans
    data_chan = data(i,:);
    data_tmp = zeros(1,ntpoints+2*nwin);
    data_tmp(nwin+1:nwin+ntpoints) = data_chan;
    data_tmp = conv(data_tmp,avg_window,'same');
    data_chan = data_tmp(nwin+1:nwin+ntpoints);
    data(i,:) = data_chan;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

