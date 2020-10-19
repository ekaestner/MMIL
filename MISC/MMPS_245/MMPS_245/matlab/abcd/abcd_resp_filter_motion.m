function motion_tseries = abcd_resp_filter_motion(motion_tseries,TR,resp_low,resp_high,order)
%function motion_tseries = abcd_resp_filter_motion(motion_tseries,TR,resp_low,resp_high,order)
%
% Purpose: apply notch filter to remove signals from respiration
%
% Required Input:
%   motion_series: matrix of motion estimates with size = [numTRs,6]
%
% Optional Input:
%   TR: repetition time (s)
%     {default = 0.8}
%  resp_low: low frequency cut-off (respirations/minute)
%    {default = 18.6}
%  resp_high: high frequency cut-off (respirations/minute)
%    {default = 25.7}
%  order: filter order
%    {default = 4}
%
% Created:  10/18/17 by Don Hagler
% Last Mod: 10/19/17 by Don Hagler
%
% Based on code from Eric Earl at OHSU (earl@ohsu.edu)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('TR','var') || isempty(TR), TR = 0.8; end;
if ~exist('resp_low','var') || isempty(resp_low), resp_low = 18.6; end;
if ~exist('resp_high','var') || isempty(resp_high), resp_high = 25.7; end;
if ~exist('order','var') || isempty(order), order = 4; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rr=[resp_low,resp_high]/60;
fs = 1/TR;
fNy=fs/2;
fa=abs(rr-floor((rr+fNy)/fs)*fs);
W_resp = fa/fNy;
Wn=mean(W_resp);
bw=diff(W_resp);
[b_filt,a_filt]=iirnotch(Wn,bw);
num_f_apply=floor(order/2); % if order<4 apply filter 1x, if order=4 2x, if order=6 3x

motion_tseries = filtfilt(b_filt,a_filt,motion_tseries);
for i=1:num_f_apply-1
  motion_tseries = filtfilt(b_filt,a_filt,motion_tseries);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

