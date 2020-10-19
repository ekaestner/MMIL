function F_filt=freq_filt(F,fs,f0,df,filt_type)
% function F_filt=freq_filt(F,fs,f0,df,filt_type)
% perform zero phase shift frequency domain filter
% no artifacts are surpressed at both ends of the filtered signal
%
% F: num_time by num_channel (column oriented)
% fs: sampling frequency
% f0: center frequency for notch or pointpass
%     For lowpass and highpass, f0 is the -6db point
%     For bandpass option, f0 has two elements [f0_low, f0_high]
% df: the width between pass and stop bands. For notch the width from
%     stop-to-stop will be 2*df; for pointpass, the width from pass-to-pass
%     will be 2*df
%     For bandpass option, df has two elements [df_low, df_high]
% file_type: 'lowpass','highpass','notch','pointpass','bandpass'
%
% The number of freq bins Nfft is chosen to be the smallest power of two
% which is larger than the the number of time points in F
% 
% last modification: 09/12/06
% Mingxiong Huang, Ph.D.

m=size(F,1);
if m==1 % one row vector
    F=F';
end
[m,n]=size(F);

Nfft=2^(ceil(log2(m))); % number of FFT terms

%F_attach=[zeros(size(F));F;zeros(size(F))]; % flip the signal and attach to bath ends
%id_orig=(m+1):2*m;

freq=linspace(0,1,Nfft)*fs;


switch lower(filt_type)
    case 'lowpass'
        Wp=f0-df/2.0; % pass point
        Ws=f0+df/2.0; % stop point
        ind=[ceil(Wp*Nfft/fs) ceil(Ws*Nfft/fs)]; % get the index number
        mask=zeros(Nfft/2,1); % the mask function in freq domain
        mask(1:ind(1))=1; % low freq side pass
        w=hanning(2*(ind(2)-ind(1))+1); % hanning window for transition
        mask((ind(1)+1):ind(2))=flipud(w(1:length((ind(1)+1):ind(2)))); % replace a portion of mask with hanning window
    case 'highpass'
        Wp=f0+df/2.0;
        Ws=f0-df/2.0;
        ind=[ceil(Ws*Nfft/fs) ceil(Wp*Nfft/fs)];
        mask=zeros(Nfft/2,1); % the mask function in freq domain
        mask((ind(2)+1):length(mask))=1;
        w=hanning(2*(ind(2)-ind(1))+1);
        mask((ind(1)+1):ind(2))=w(1:length((ind(1)+1):ind(2)));
    case 'notch'
        Wp1=f0-df; % low freq side, pass point, not divided by 2
        Wp2=f0+df; % high freq side, pass point
        ind=[ceil(Wp1*Nfft/fs) ceil(Wp2*Nfft/fs)];
        mask=zeros(Nfft/2,1); % the mask function in freq domain
        mask=mask+1;
        w=hanning(ind(2)-ind(1)+1);
        mask(ind(1):ind(2))=mask(ind(1):ind(2))-w;    
    case 'pointpass'
        Ws1=f0-df; % low freq side, stop point
        Ws2=f0+df; % high freq side, stop point
        ind=[ceil(Ws1*Nfft/fs) ceil(Ws2*Nfft/fs)];
        w=hanning(ind(2)-ind(1)+1);
        mask=zeros(Nfft/2,1); % the mask function in freq domain
        mask(ind(1):ind(2))=w;
    case 'bandpass'
        if length(f0) ~=2 | length(df)~=2,
            error('For bandpass filter, f0 and df must have 2 elements\n');
        end    
        % first the high pass part
        Wp_high=f0(1)+df(1)/2.0;
        Ws_high=f0(1)-df(1)/2.0;
        ind_high=[ceil(Ws_high*Nfft/fs) ceil(Wp_high*Nfft/fs)];
        mask_high=zeros(Nfft/2,1); % the mask function in freq domain
        mask_high((ind_high(2)+1):length(mask_high))=1;
        w=hanning(2*(ind_high(2)-ind_high(1))+1);
        mask_high((ind_high(1)+1):ind_high(2))=w(1:length((ind_high(1)+1):ind_high(2)));
        % now the low pass part
        Wp_low=f0(2)-df(2)/2.0; % pass point
        Ws_low=f0(2)+df(2)/2.0; % stop point
        ind_low=[ceil(Wp_low*Nfft/fs) ceil(Ws_low*Nfft/fs)]; % get the index number
        mask_low=zeros(Nfft/2,1); % the mask function in freq domain
        mask_low(1:ind_low(1))=1; % low freq side pass
        w=hanning(2*(ind_low(2)-ind_low(1))+1); % hanning window for transition
        mask_low((ind_low(1)+1):ind_low(2))=flipud(w(1:length((ind_low(1)+1):ind_low(2)))); % replace a portion of mask with hanning window
        mask=mask_high.*mask_low;
end

mask=[mask;flipud(mask)];

F_fft=fft(F,Nfft);
F_fft=F_fft.*(mask*ones(1,n));

F_filt_temp=real(ifft(F_fft,Nfft));
F_filt=F_filt_temp(1:m,:);
        
        