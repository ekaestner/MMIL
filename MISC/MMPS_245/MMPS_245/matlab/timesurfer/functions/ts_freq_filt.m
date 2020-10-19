function F_filt=ts_freq_filt(F,fs,f0,df,filt_type,zpad_flag)
% function F_filt=freq_filt(F,fs,f0,df,filt_type,[zpad_flag])
%
% Purpose: perform zero phase shift frequency domain filter
%
% Required Input:
%   F: num_time by num_channel (column oriented)
%   fs: sampling frequency
%   f0: center frequency for notch or pointpass
%     For lowpass and highpass, f0 is the -6db point
%     For bandpass option, f0 has two elements [f0_low, f0_high]
%     For bandpassnotch option, f0 has three elements [f0_low, f0_high, f0_notch]
%   df: the width between pass and stop bands. For notch the width from
%     stop-to-stop will be 2*df; for pointpass, the width from pass-to-pass
%     will be 2*df
%     For bandpass option, df has two elements [df_low, df_high]
%     For bandpassnotch option, df has three elements [df_low, df_high, df_notch]
%   filt_type: 'lowpass','highpass','notch','pointpass','bandpass','bandpassnotch'
%
% Optional Input:
%   zpad_flag: [0|1] pad time series with zeros on each end to supress
%     edge artifacts
%     {default = 0}
%
% Notes:
%   The number of freq bins Nfft is chosen to be the smallest power of two
%     which is larger than the the number of time points in F
%   f0+df should be less than fs/2 and f0-df should be greater than 0
%     if not, values will be reset to those limits
%   no artifacts are surpressed at ends of the filtered signal
%     unless zpad_flag = 1
%
% Created:   09/12/06 by Mingxiong Huang, Ph.D.
% Rcnt Mod:  11/23/09 by Don Hagler
% Last Mod:  08/17/12 by Don Hagler
%

F_filt = [];

if ~mmil_check_nargs(nargin,5), return; end;
if ~exist('zpad_flag','var') || isempty(zpad_flag), zpad_flag = 0; end;

m=size(F,1);
if m==1 % one row vector
  F=F';
end
[m,n]=size(F);

% pad with zeros on each end to avoid edge artifacts
if zpad_flag
  ind_orig = [m+1:2*m];
  F = [zeros(size(F));F;zeros(size(F))];
  [m,n]=size(F);
end;

Nfft=2^(ceil(log2(m))); % number of FFT terms

freq=linspace(0,1,Nfft)*fs;

% check input
if length(f0)~=length(df)
  error('length of f0 does not match df');
end;

% set frequency limits to avoid errors
flo = 0;
fhi = fs/2;

for i=1:length(f0)
  if f0(i) - df(i)/2.0 < flo
    fprintf('%s: WARNING: f0(%d)-df(%d)/2 is less than %0.3f\n',...
      mfilename,i,i,flo);
  elseif f0(i) + df(i)/2.0 > fhi
    fprintf('%s: WARNING: f0(%d)+df(%d)/2 is greater than %0.3f\n',...
      mfilename,i,i,fhi);
  end;
end;

switch lower(filt_type)
  case 'lowpass'
    Wp=min(max(f0-df/2.0,flo),fhi); % pass point
    Ws=min(max(f0+df/2.0,flo),fhi); % stop point
    ind=[ceil(Wp*Nfft/fs) ceil(Ws*Nfft/fs)]; % get the index number
    mask=zeros(Nfft/2,1); % the mask function in freq domain
    mask(1:ind(1))=1; % low freq side pass
    w=hanning(2*(ind(2)-ind(1))+1); % hanning window for transition
    mask((ind(1)+1):ind(2))=flipud(w(1:length((ind(1)+1):ind(2)))); % replace a portion of mask with hanning window
  case 'highpass'
    Wp=min(max(f0+df/2.0,flo),fhi);
    Ws=min(max(f0-df/2.0,flo),fhi);
    ind=[ceil(Ws*Nfft/fs) ceil(Wp*Nfft/fs)];
    mask=zeros(Nfft/2,1); % the mask function in freq domain
    mask((ind(2)+1):length(mask))=1;
    w=hanning(2*(ind(2)-ind(1))+1);
    mask((ind(1)+1):ind(2))=w(1:length((ind(1)+1):ind(2)));
  case 'notch'
    Wp1=min(max(f0-df,flo),fhi); % low freq side, pass point, not divided by 2
    Wp2=min(max(f0+df,flo),fhi); % high freq side, pass point
    ind=[ceil(Wp1*Nfft/fs) ceil(Wp2*Nfft/fs)];
    mask=zeros(Nfft/2,1); % the mask function in freq domain
    mask=mask+1;
    w=hanning(ind(2)-ind(1)+1);
    mask(ind(1):ind(2))=mask(ind(1):ind(2))-w;    
  case 'pointpass'
    Ws1=min(max(f0-df,flo),fhi); % low freq side, stop point
    Ws2=min(max(f0+df,flo),fhi); % high freq side, stop point
    ind=[ceil(Ws1*Nfft/fs) ceil(Ws2*Nfft/fs)];
    w=hanning(ind(2)-ind(1)+1);
    mask=zeros(Nfft/2,1); % the mask function in freq domain
    mask(ind(1):ind(2))=w;
  case 'bandpass'
    if length(f0) ~=2 | length(df)~=2,
        error('For bandpass filter, f0 and df must have 2 elements\n');
    end    
    % first the high pass part
    Wp_high=min(max(f0(1)+df(1)/2.0,flo),fhi);
    Ws_high=min(max(f0(1)-df(1)/2.0,flo),fhi);
    ind_high=[ceil(Ws_high*Nfft/fs) ceil(Wp_high*Nfft/fs)];
    mask_high=zeros(Nfft/2,1); % the mask function in freq domain
    mask_high((ind_high(2)+1):length(mask_high))=1;
    w=hanning(2*(ind_high(2)-ind_high(1))+1);
    mask_high((ind_high(1)+1):ind_high(2))=w(1:length((ind_high(1)+1):ind_high(2)));
    % now the low pass part
    Wp_low=min(max(f0(2)-df(2)/2.0,flo),fhi); % pass point
    Ws_low=min(max(f0(2)+df(2)/2.0,flo),fhi); % stop point
    ind_low=[ceil(Wp_low*Nfft/fs) ceil(Ws_low*Nfft/fs)]; % get the index number
    mask_low=zeros(Nfft/2,1); % the mask function in freq domain
    mask_low(1:ind_low(1))=1; % low freq side pass
    w=hanning(2*(ind_low(2)-ind_low(1))+1); % hanning window for transition
    mask_low((ind_low(1)+1):ind_low(2))=flipud(w(1:length((ind_low(1)+1):ind_low(2)))); % replace a portion of mask with hanning window
    mask=mask_high.*mask_low;
  case 'bandpassnotch'
    if length(f0) ~=3 | length(df)~=3,
        error('For bandpassnotch filter, f0 and df must have 3 elements\n');
    end    
    % first the high pass part
    Wp_high=min(max(f0(1)+df(1)/2.0,flo),fhi);
    Ws_high=min(max(f0(1)-df(1)/2.0,flo),fhi);
    ind_high=[ceil(Ws_high*Nfft/fs) ceil(Wp_high*Nfft/fs)];
    mask_high=zeros(Nfft/2,1); % the mask function in freq domain
    mask_high((ind_high(2)+1):length(mask_high))=1;
    w=hanning(2*(ind_high(2)-ind_high(1))+1);
    mask_high((ind_high(1)+1):ind_high(2))=w(1:length((ind_high(1)+1):ind_high(2)));
    % now the low pass part
    Wp_low=min(max(f0(2)-df(2)/2.0,flo),fhi); % pass point
    Ws_low=min(max(f0(2)+df(2)/2.0,flo),fhi); % stop point
    ind_low=[ceil(Wp_low*Nfft/fs) ceil(Ws_low*Nfft/fs)]; % get the index number
    mask_low=zeros(Nfft/2,1); % the mask function in freq domain
    mask_low(1:ind_low(1))=1; % low freq side pass
    w=hanning(2*(ind_low(2)-ind_low(1))+1); % hanning window for transition
    mask_low((ind_low(1)+1):ind_low(2))=flipud(w(1:length((ind_low(1)+1):ind_low(2)))); % replace a portion of mask with hanning window
    % now the notch part
    Wp1=min(max(f0(3)-df(3),flo),fhi); % low freq side, pass point, not divided by 2
    Wp2=min(max(f0(3)+df(3),flo),fhi); % high freq side, pass point
    ind=[ceil(Wp1*Nfft/fs) ceil(Wp2*Nfft/fs)];
    mask_notch=zeros(Nfft/2,1); % the mask function in freq domain
    mask_notch=mask_notch+1;
    w=hanning(ind(2)-ind(1)+1);
    mask_notch(ind(1):ind(2))=mask_notch(ind(1):ind(2))-w;    
    % combine the masks
    mask=mask_high.*mask_low.*mask_notch;
end

mask=[mask;flipud(mask)];

F_fft=fft(F,Nfft).*(mask*ones(1,n));
F_filt_temp=real(ifft(F_fft,Nfft));
F_filt=F_filt_temp(1:m,:);

if zpad_flag
  F_filt=F_filt(ind_orig,:);
end;


