function forward = rc_rot_forward(forward,rot)
%function forward = rc_rot_forward(forward,rot)
%
% Early Mod: 06/29/09 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
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

ndips = size(forward.G_xyz,2)/3;
for i=1:ndips
  j = (i-1)*3 + 1;
  k = j+2;
  xyz = forward.G_xyz(:,j:k);
  new_xyz = (R*xyz')';
  forward.G_xyz(:,j:k) = new_xyz;
end;
