function vec = mmil_M_mat2vec(M)
%function vec = mmil_M_mat2vec(M)
%
% Purpose: Converts 4*4 transformation matrix
%   to vector of rotation and translation
%
% Required Input:
%   M: 4x4 transformation matrix
%
% Output:
%   vec: [transx, transy, transz, rotx, roty, rotz]
%      translation in mm and rotation in degrees clockwise
%      rotations would be applied in order [x,y,z]
%
% Created:  07/26/11 by Vijay Venkatraman
% Last Mod: 08/28/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

vec = zeros(6,1);
vec(1:3) = M(1:3,4);
vec(4) = atan2(M(3,2),M(3,3));
vec(5) = -asin(M(3,1));
vec(6) = atan2(M(2,1),M(1,1));
vec(4:6) = -vec(4:6)*180/pi;

return;

