function T = dti_components2tensor(vec)
%function T = dti_components2tensor(vec)
%
% Purpose: reshape vector of tensor components to 3x3 tensor matrix
% 
% vec should be either 6 or 9 element vector
%
%   For 6 elements:
%    assumes order used for tensor calculations: xx,yy,zz,xy,xz,yz
%    vec = [1 2 3 4 5 6]
%      T = [1 4 5 ; 4 2 6 ; 5 6 3]
%
%   For 9 elements:
%    vec = [1 2 3 4 5 6 7 8 9]
%      T = [1 2 3 ; 4 5 6 ; 7 8 9]
%
%
% Created:  08/24/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

vec = reshape(vec,[numel(vec),1]);
switch length(vec)
  case 9
    T = reshape(vec,[3,3])';
  case 6
    T = diag(vec(1:3));
    T(1,2) = vec(4); T(2,1) = T(1,2);
    T(1,3) = vec(5); T(3,1) = T(1,3);
    T(2,3) = vec(6); T(3,2) = T(2,3);
  otherwise
    error('vec length must be 6 or 9 (is %d)',length(vec));
end;


return;
