function trans = rc_rotate_trans(trans,rot)
%function trans = rc_rotate_trans(trans,rot)
%
% Purpose: rotate a transformation matrix
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/04/11 by Don Hagler
%
%

if ~mmil_check_nargs(nargin,2), return; end;

rot = rot*pi/180; % expect degrees, convert to radians
Rx=eye(3);
Ry=eye(3);
Rz=eye(3);
if(rot(1))
  c = cos(rot(1));
  s = sin(rot(1));
  Rx(2,2)=c;
  Rx(2,3)=s;
  Rx(3,2)=-s;
  Rx(3,3)=c;
end
if(rot(2))
  c = cos(rot(2));
  s = sin(rot(2));
  Ry(1,1)=c;
  Ry(1,3)=-s;
  Ry(3,1)=s;
  Ry(3,3)=c;
end
if(rot(3))
  c = cos(rot(3));
  s = sin(rot(3));
  Rz(1,1)=c;
  Rz(1,2)=s;
  Rz(2,1)=-s;
  Rz(2,2)=c;
end
R = Rz*Ry*Rx;

% adjust trans matrix so normal vectors get rotated in gain_xyz2norm
tmp_trans = trans(1:3,1:3);
tmp_trans = inv(R*inv(tmp_trans));
trans(1:3,1:3) = tmp_trans;

