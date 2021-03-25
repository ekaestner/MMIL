% ts_filtfft() -  (high|low|band)-pass filter data using inverse fft 
%                 (without using the Matlab signal processing toolbox)
% Usage:
%  >> [filtdata] = filtfft(data,srate,locutoff,hicutoff);
%  >> [filtdata] = filtfft(data,srate,locutoff,hicutoff,loband,hiband);
%
% Inputs:
%   data        = (channels,frames) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {if 0, lowpass only}
%   hicutoff    = high-edge frequency in pass band (Hz) {if 0, highpass only}
%
% Optional Inputs:
%   loband      = low-edge frequency band (Hz)  {0}
%   hiband      = high-edge frequency band (Hz) {0}
%
% Outputs:
%    filtdata = smoothed data
%
% Original author: Tao Song
% Modified by Don Hagler 03/24/06
% Last modifiied: 07/31/06 DH

function filtdata = ts_filtfft_tband(data,srate,lowcut,highcut,lowband,highband);
    
if nargin < 4
  help filtfft;
end;
[chans,frames] = size(data);

try, lowband;  catch, lowband = 0;  end;
try, highband; catch, highband = 0; end;

freq=linspace(0, 1, frames)*srate;

ind_pass_low=min(max(1,ceil((lowcut+lowband/2)*length(freq)/srate)),length(freq));
ind_stop_low=min(max(1,ceil((lowcut-lowband/2)*length(freq)/srate)),length(freq));
ind_pass_high=min(max(1,ceil((highcut-highband/2)*length(freq)/srate)),length(freq));
ind_stop_high=min(max(1,ceil((highcut+highband/2)*length(freq)/srate)),length(freq));

% construct frequency mask
mask=zeros(size(freq));
mask(ind_pass_low:ind_pass_high)=1;
w=hanning(2*(ind_stop_high-ind_pass_high)+1);
mask((ind_stop_high+1):ind_pass_high)=w(1:length((ind_stop_high+1):ind_pass_high));
w=hanning(2*(ind_pass_low-ind_stop_low)+1);
mask((ind_pass_low+1):-1:(ind_stop_low+1))=w(ceil(length(w)/2):length(w));
data_freq=fft(data',length(freq))'.*(mask'*ones(1,chans))';
filtdata=real(ifft((data_freq)', length(freq))');
filtdata=2*filtdata(:,1:frames);

return;
