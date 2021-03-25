function acc = ts_calc_bs_acceleration(jstat)
%function acc = ts_calc_bs_acceleration(jstat)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acc = [];
if ~mmil_check_nargs(nargin,1), return; end;

% compute acceleration (see DiCiccio and Efron (1996))
N = size(jstat,2);
weights = repmat(1/N,N,1);
mjstat = mean(jstat,2);
score = bsxfun(@minus,mjstat,jstat); % score function at stat
iszer = all(score==0,2);
skew = sum(score.^3,2)./(sum(score.^2,2).^1.5); % skewness of score function
skew(iszer) = 0;
acc = skew/6;  % acceleration

return;
