function [synth_data,G_norm] = ts_synth_sensors_from_dSPM(stc,varargin)
%function [synth_data,G_norm] = ts_synth_sensors_from_dSPM(stc,[options])
%
% Usage:
%  [synth_data,G_norm] = ts_synth_sensors_from_dSPM(stc, 'key1', value1,...);
%
% Required Input:
%   stc: source time courses (nsources x time points)
%
% Optional Inupt:
%   'prefix': prefix of ts_dSPM output files
%     {default = 'dSPM'}
%   'prefix_data': prefix of processed data file (e.g. 'proc_avg_data')
%     if not supplied, will use avg_data saved by ts_dSPM
%     {default = []}
%   'time': time vector (should be in seconds)
%     if empty or ommitted, will use time vector in avg_data struct
%     saved by ts_dSPM
%     {default = []}
%   'noisefact': stdev of noise added to synth_data,
%      relative to stdev of noise in avg_data
%     {default = 0}
%   'grad_noise': stdev of gradiometer noise, in fT/cm
%     may use instead of or in addition to noisefact
%     {default = 0}
%   'mag_noise': stdev of magnetometer noise, in fT
%     {default = 0}
%   'EEG_noise': stdev of EEG noise, in uV
%     {default = 0}
%   'rootdir': directory containing matfiles dir
%     {default = pwd}
%   'bandpass_flag': [0|1] bandpass filter
%     {default = 0}
%   'bandpass_low_cf': low cutoff frequency (high-pass filter) (Hz)
%     {default = 0}
%   'bandpass_low_tb': low cutoff transition band (Hz)
%     {default = 0}
%   'bandpass_high_cf': high cutoff frequency (low-pass filter) (Hz)
%     {default = 100}
%   'bandpass_high_tb': high cutoff transition band (Hz)
%     {default = 0}
%   'notch_flag': [0|1] notch filter
%     {default = 0}
%   'notch_cf': notch center frequency (notch filter) (Hz)
%     {default = 0}
%   'notch_tb': notch transition band (Hz)
%     {default = 0}
%   'G_norm': normal vector gain matrix
%     {default = []}
%   'verbose': [0|1] display status messages
%     {default = 1}
%
% Output:
%   synth_data: synthesized data in avg_data structure
%   G_norm: normal vector gain matrix
%     NOTE: may use this as input for subsequent calls of this function
%       to avoid having to recalculate this each time
%
% Created:  08/01/07 by Don Hagler
% Last Mod: 03/06/14 by Don Hagler
%

synth_data = []; G_norm = [];

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix','dSPM',[],...
  'prefix_data',[],[],...
  'time',[],[],...
  'noisefact',0,[0,Inf],...
  'grad_noise',0,[0,Inf],...
  'mag_noise',0,[0,Inf],...
  'EEG_noise',0,[0,Inf],...
  'rootdir',pwd,[],...
  'bandpass_flag',false,[false,true],...
  'bandpass_low_cf',0,[0,Inf],...
  'bandpass_low_tb',0,[0,Inf],...
  'bandpass_high_cf',100,[0,Inf],...
  'bandpass_high_tb',0,[0,Inf],...
  'notch_flag',false,[false,true],...
  'notch_cf',0,[0,Inf],...
  'notch_tb',0,[0,Inf],...  
  'G_norm',[],[],...
  'verbose',true,[false true],...
...
  'pad_nsamples',1000,[0,Inf],...
});

if ~isempty(parms.time)
  parms.time = mmil_rowvec(parms.time);
  if size(parms.time,2)~=size(stc,2)
    error('length of time vector (%d) does not match time points in stc',...
      size(parms.time,2),size(stc,2));
  end;
end;

matfile=sprintf('%s/matfiles/%s_parms.mat',parms.rootdir,parms.prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
tmp_parms = parms;
load(matfile);
fields = fieldnames(tmp_parms);
for f=1:length(fields)
  parms = setfield(parms,fields{f},getfield(tmp_parms,fields{f}));
end;

matfile=sprintf('%s/matfiles/%s_forward_prep.mat',parms.rootdir,parms.prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

if isempty(parms.prefix_data)
  matfile=sprintf('%s/matfiles/%s_avg_data.mat',parms.rootdir,parms.prefix);
else
  matfile=sprintf('%s/matfiles/%s.mat',parms.rootdir,parms.prefix_data);
end;
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

synth_data = avg_data;
synth_data.averages = synth_data.averages(1);
if isempty(parms.time)
  parms.time = synth_data.averages.time(1:size(stc,2));
end;

if isempty(parms.G_norm)
  if parms.verbose
    fprintf('%s: calculating G_norm from G_xyz...\n',mfilename);
    tic;
  end;
  G_norm = ts_gain_xyz2norm(G_xyz,...
    parms.lh_dip_info,parms.rh_dip_info,...
    parms.lh_dec_dips,parms.rh_dec_dips,...
    parms.trans);
  if parms.verbose, toc; end;
else
  G_norm = parms.G_norm;
end;

if size(G_norm,2) ~= size(stc,1)
  error('number of sources in stc (%d) does not match G_norm (%d)\n',...
    size(stc,1),size(G_norm,2));
end;

tmp_data = G_norm*stc;
[tmp,ind_grad] = intersect(parms.goodchans,parms.grad_chans);
[tmp,ind_mag] = intersect(parms.goodchans,parms.mag_chans);
[tmp,ind_EEG] = intersect(parms.goodchans,parms.EEG_chans);
tmp_data(ind_grad,:) = tmp_data(ind_grad,:)/parms.grad_scalefact;
tmp_data(ind_mag,:) = tmp_data(ind_mag,:)/parms.mag_scalefact;
tmp_data(ind_EEG,:) = tmp_data(ind_EEG,:)/parms.EEG_scalefact;

% add noise
if parms.noisefact~=0
  noiselevels = sqrt(diag(parms.noisecovar))*parms.noisefact;
  tmp_noise = randn(size(tmp_data)).*(noiselevels*ones(1,length(parms.time)));
  tmp_noise(ind_grad,:) = tmp_noise(ind_grad,:)/parms.grad_scalefact;
  tmp_noise(ind_mag,:) = tmp_noise(ind_mag,:)/parms.mag_scalefact;
  tmp_noise(ind_EEG,:) = tmp_noise(ind_EEG,:)/parms.EEG_scalefact;
  tmp_data = tmp_data + tmp_noise;
end;
if parms.grad_noise
  tmp_noise = zeros(size(tmp_data));
  tmp_noise(ind_grad,:) = randn(length(ind_grad),size(tmp_data,2)) *...
                          parms.grad_noise/parms.grad_scalefact;
  tmp_data = tmp_data + tmp_noise;
end;
if parms.mag_noise
  tmp_noise = zeros(size(tmp_data));
  tmp_noise(ind_mag,:) = randn(length(ind_mag),size(tmp_data,2)) *...
                          parms.mag_noise/parms.mag_scalefact;
  tmp_data = tmp_data + tmp_noise;
end;
if parms.EEG_noise
  tmp_noise = zeros(size(tmp_data));
  tmp_noise(ind_EEG,:) = randn(length(ind_EEG),size(tmp_data,2)) *...
                          parms.EEG_noise/parms.EEG_scalefact;
  tmp_data = tmp_data + tmp_noise;
end;

synth_data.averages.data = zeros(size(synth_data.averages.data,1),...
  size(stc,2));
synth_data.averages.stdev = zeros(size(synth_data.averages.data,1),...
  size(stc,2));
synth_data.averages.time = parms.time;
synth_data.averages.data(parms.goodchans,:) = tmp_data;
synth_data.sfreq = 1000 / (parms.time(2) - parms.time(1));

% filter data
if parms.bandpass_flag
  if parms.verbose
    fprintf('%s: filtering synthesized data...\n',mfilename);
    tic;
  end;
  synth_data = ts_postprocess_avg(synth_data,...
    'baseline_flag',0,...
    'bandpass_flag',parms.bandpass_flag,...
    'bandpass_low_cf',parms.bandpass_low_cf,...
    'bandpass_low_tb',parms.bandpass_low_tb,...
    'bandpass_high_cf',parms.bandpass_high_cf,...
    'bandpass_high_tb',parms.bandpass_high_tb,...
    'notch_flag',parms.notch_flag,...
    'notch_cf',parms.notch_cf,...
    'notch_tb',parms.notch_tb,...
    'pad_nsamples',parms.pad_nsamples);
  if parms.verbose, toc; end;
end;

