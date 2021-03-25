function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width);
% function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width);
%
% Calculates the average of a time-frequency energy representation of
% multiple trials using a Morlet wavelet method.                            
%
% Input
% -----
% S    : signals = time x Trials      
% freqVec    : frequencies over which to calculate TF energy        
% Fs   : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
%
% Output
% ------
% t    : time
% f    : frequency
% B    : phase-locking factor = frequency x time
%
%     
%
%------------------------------------------------------------------------
% Ole Jensen, November 1999
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

S = S';

timeVec = (1:size(S,2))/Fs;  

B = zeros(length(freqVec),size(S,2)); 

for i=1:size(S,1)          
%    fprintf(1,'%d ',i); 
    for j=1:length(freqVec)
        B(j,:) = energyvec(freqVec(j),detrend(S(i,:)),Fs,width) + B(j,:);
    end
end
%fprintf('\n'); 
TFR = B/size(S,1);     


function y = energyvec(f,s,Fs,width)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
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

y = abs(y).^2;
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
