function R = rc_calc_ref_corr(S,S_ref,sample_range)
%function R = rc_calc_ref_corr(S,S_ref,sample_range)
%
% Purpose: calculate correlation between two waveforms
%
% Required Input:
%   S: source waveform matrix
%   S_ref: references source waveform matrix
%   sample_range: vector of start and end samples from which to calculate diff
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 04/29/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

fprintf('%s: calculating correlation with reference waveforms...\n',mfilename);
R = 0;
n = 0;
for s=1:size(S,2)
  for c=1:size(S,3)
    wform = S(sample_range(1):sample_range(2),s,c);
    wform_ref = S_ref(sample_range(1):sample_range(2),s,c);
    R = R + corr(wform,wform_ref);
    n = n + 1;
  end;
end;
R = R/n;

