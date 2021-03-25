% ts_filtfft() -  (high|low|band)-pass filter data using inverse fft 
%                 (without using the Matlab signal processing toolbox)
% Usage:
%  >> [filtdata] = filtfft(data,srate,locutoff,hicutoff);
%  >> [filtdata] = filtfft(data,srate,locutoff,hicutoff,notchfilt);
%
% Inputs:
%   data        = (channels,frames) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {if 0, lowpass only}
%   hicutoff    = high-edge frequency in pass band (Hz) {if 0, highpass only}
%
% Optional inputs:
%   notchfilt   = [0|1] notch filter instead of bandpass {0}
%
% Outputs:
%    filtdata = smoothed data
%
% Known problems:
%    The signal drop off is much smaller compared to standard filtering methods
%
% Original author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2003
% Modified by Don Hagler 03/24/06
%
% See also: eegfiltfft()

%% original copyright info from EEGLAB's eegfiltfft.m:
% inspired from a ggogle group message
% http://groups.google.com/groups?q=without+%22the+signal+processing+toolbox%22&hl=en&lr=&ie=UTF-8&oe=UTF-8&selm=f56893ae.0311141025.3069d4f8%40posting.google.com&rnum=8
% Copyright (C) Arnaud Delorme, SCCN/INC/UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% Revision 1.7  2004/03/04 19:30:41  arno
% Revision 1.1  2003/12/03 18:26:17  arno

function filtdata = ts_filtfft(data, fs, lowcut, highcut, notchfilt);
    
if nargin < 4
  help(mfilename);
end;
[chans frames] = size(data);

try, notchfilt; catch, notchfilt = 0; end;

fv=reshape([0:frames-1]*fs/frames,frames,1); % Frequency vector for plotting    

% find closest freq in fft decomposition
% --------------------------------------
if lowcut ~= 0
  [tmp idxl]=min(abs(fv-lowcut));  % Find the entry in fv closest to 5 kHz
else
  idxl = 0;
end;
if highcut ~= 0        
  [tmp idxh]=min(abs(fv-highcut));  % Find the entry in fv closest to 5 kHz    
else 
  idxh = length(fv)/2;
end;

% filter the data
% ---------------
filtdata = zeros(chans,frames);
for c=1:chans
  X=fft(data(c,:));
  if notchfilt
    X(idxl+1:idxh-1)=0;
    X(end/2:end)=0;
  else
    X(1:idxl)=0;
    X(end-idxl:end)=0;
    X(idxh:end)=0;
  end;                
  filtdata(c,:) = 2*real(ifft(X));
end

return;
