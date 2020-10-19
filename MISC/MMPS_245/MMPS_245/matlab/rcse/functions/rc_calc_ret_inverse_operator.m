function W = rc_calc_ret_inverse_operator(G,SNR,R,C);
%function W = rc_calc_ret_inverse_operator(G,SNR,R,C);
%
% Purpose: calculate minimum norm inverse operator from gain matrix
%   using supplied source covariance matrix
%
% G: forward gain matrix of sensor amplitudes for each decimated dipole
%      possibly with 3 orientations each
% SNR: estimated signal to noise ratio
% R: source covariance matrix
% C: sensor noise covariance matrix
%
% created:  12/10/05 by Don Hagler
% last mod: 03/08/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,4), return; end;

W=[];

[num_sensors,num_sources] = size(G);

[temp1,temp2] = size(R);
if(temp1 ~= temp2)
  error('source covariance matrix is not square');
end
if(temp1 ~= num_sources)
  error('source covariance matrix has wrong number of elements (%d)',temp1);
end

[temp1,temp2] = size(C);
if(temp1 ~= temp2)
  error('sensor noise covariance matrix is not square');
end
if(temp1 ~= num_sensors)
  error('sensor noise covariance matrix has wrong number of elements (%d)',temp1);
end

tic
%tmp = G*R*G';
tmp = R*G';
tmp = G*tmp; % doing this in two steps uses half as much memory
lambdaC = C*(mean(diag(tmp))/mean(full(diag(C))))/SNR^2; % power SNR

if(num_sources>num_sensors)
%  fprintf('%s: num_sources (%d) > num_sensors (%d)\n',...
%          mfilename,num_sources,num_sensors);
  % this is the typical underdetermined inverse
  W = R*G'*inv(tmp + lambdaC);
else
%  fprintf('%s: num_sources (%d) < num_sensors (%d)\n',...
%          mfilename,num_sources,num_sensors);
  % special case with more sensors (measurements) than sources
  invlambdaC = inv(lambdaC);
  invR = inv(R);
  W = inv(G'*invlambdaC*G + invR)*G'*invlambdaC;
end
toc;

return;

