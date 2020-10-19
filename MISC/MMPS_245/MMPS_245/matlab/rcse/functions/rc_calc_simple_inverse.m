function W = rc_calc_simple_inverse(G);
%function W = rc_calc_simple_inverse(G);
%
% Purpose: calculate minimum norm inverse from gain matrix
%
% G: forward gain matrix of sensor amplitudes for each decimated dipole
%      possibly with 3 orientations each
%
% created:  07/28/09 by Don Hagler
% last mod: 03/08/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

W=[];

if isempty(G), error('gain matrix G is empty'); end

[num_sensors,num_sources] = size(G);

if(num_sources>num_sensors)
%  fprintf('%s: num_sources (%d) > num_sensors (%d)\n',...
%          mfilename,num_sources,num_sensors);
  % this is the typical underdetermined inverse
  W = G'*inv(G*G');
else
%  fprintf('%s: num_sources (%d) < num_sensors (%d)\n',...
%          mfilename,num_sources,num_sensors);
  % special case with more sensors (measurements) than sources
  W = inv(G'*G)*G';
end

return;

