function vec = dti_tensor2components(T,num_tensor_components)
%function vec = dti_tensor2components(T,[num_tensor_components])
%
% Purpose: reshape 3x3 tensor matrix to vector of tensor components
% 
% if num_tensor_components = 9, will return 9 element vector
% otherwise, returns 6 element vector, assuming T is symmetric
%
% Example:
%   T = [1 2 3 ; 4 5 6 ; 7 8 9]
%
%    num_tensor_components = 6:
%      assumes order used for tensor calculations: xx,yy,zz,xy,xz,yz
%      vec = [1 5 9 2 3 6]
%
%    num_tensor_components = 9:
%      vec = [1 2 3 4 5 6 7 8 9]
%
% Created:  08/24/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('num_tensor_components','var') || isempty(num_tensor_components)
  num_tensor_components = 6;
end;

if any(size(T)~=[3,3])
  error('T must be 3x3 matrix');
end;
T = reshape(T',[1,9]);

switch num_tensor_components
  case 9
    vec = T;
  case 6
    vec = T([1 5 9 2 3 6]);
  otherwise
    error('num_tensor_components must be 6 or 9 (is %d)',num_tensor_components);
end;

return;
