function vol = mmil_bandpass_vol(vol,varargin)
%function vol = mmil_bandpass_vol(vol,[options])
%
% Required Input:
%   vol: 4D matrix with brain time series
%
% Optional Input:
%  'TR': sample duration (sec)
%    {default = 1}
%  'bandpass_low': low frequency cut-off (Hz)
%    {default = 0.01}
%  'bandpass_high': high frequency cut-off (Hz)
%    {default = 0.08}
%  'bandpass_tb': bandpass filter transition band (Hz)
%    {default = 0.001}
%  'bandpass_zpad_flag': [0|1] pad time series with zeros on each end to
%     supress edge artifacts
%     {default = 1}
%
% Created:  04/05/12 by Don Hagler
% Last Mod: 11/01/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'TR',1,[],...
  'bandpass_low',0.01,[],...
  'bandpass_high',0.08,[],...
  'bandpass_tb',0.001,[],...
  'bandpass_zpad_flag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% band-pass filter volume time course

volsz = size(vol);
if size(volsz)~=4, error('input vol is not 4D'); end;

% loop over slices to save memory
for z=1:volsz(3)
  % reshape to nvox x ntpoints
  vec = reshape(vol(:,:,z,:),[prod(volsz(1:2)),volsz(4)]);
  % FFT-based filter
  vec = ts_freq_filt(vec',1/parms.TR,...
    [parms.bandpass_low,parms.bandpass_high],...
    [parms.bandpass_tb,parms.bandpass_tb],...
    'bandpass',parms.bandpass_zpad_flag)';
  % reshape back to original size
  vol(:,:,z,:) = reshape(vec,[volsz(1),volsz(2),1,volsz(4)]);
end;

