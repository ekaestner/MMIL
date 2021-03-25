function W = rc_calc_ret_inverse_operator_noprior(G,SNR)
%function W = rc_calc_ret_inverse_operator_noprior(G,SNR)
%
% Purpose: calculate minimum norm inverse from gain matrix
%          assumes source and noise covariance R = I, C = I
%
% Required Input:
%   G: forward gain matrix of sensor amplitudes for each decimated dipole
%      possibly with 3 orientations each
%   SNR: estimated signal to noise ratio
%
% Created:  12/10/05 by Don Hagler
% Last Mod: 02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

W=[];

if isempty(G), error('gain matrix G is empty'); end

[num_sensors,num_sources] = size(G);

if(num_sources>num_sensors)
%  fprintf('%s: num_sources (%d) > num_sensors (%d)\n',...
%          mfilename,num_sources,num_sensors);
  % this is the typical underdetermined inverse
  temp = G*G';
  lambda = mean(diag(temp))/SNR;
  I = eye(num_sensors);
  W = G'*inv(G*G' + lambda*I);
else
%  fprintf('%s: num_sources (%d) < num_sensors (%d)\n',...
%          mfilename,num_sources,num_sensors);
  % special case with more sensors (measurements) than sources
  temp = G'*G;
  lambda = mean(diag(temp))/SNR;
  I = eye(num_sources);
  W = inv(G'*G + lambda*I)*G';
end

return;

