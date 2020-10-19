function TF_complex = traces2TF_complex(S,freqVec,Fs,width);
% function TF_complex = traces2TF_complex(S,freqVec,Fs,width);
%
% Calculates the time-frequency convolusion coefficient
% multiple trials using a Morlet wavelet method.                            
%
% Input
% -----
% S             : 1 by n_time      
% freqVec       : frequencies over which to calculate TF energy        
% Fs            : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
%
% Output
% ------
% TF_complex    : complex time-freq convolution  
%
%     
%
%------------------------------------------------------------------------
% (c) Mingxiong Huang 03/15/06 
%------------------------------------------------------------------------


timeVec = (1:length(S))/Fs;  

TF_complex = zeros(length(freqVec),length(S)); 

for j=1:length(freqVec)
    TF_complex(j,:) = conv_coef(freqVec(j),detrend(S),Fs,width);
end
     

function y = conv_coef(f,s,Fs,width)
% function y = conv_coef(f,s,Fs,width)
%
% Return a vector containing the complex convolution coefficient
% as a function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
%
% 

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);

y = conv(s,m);

y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));



function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
