function new_coords = rc_rot_coords(coords,rotx,roty,rotz);
%function new_coords = rc_rot_coords(coords,rotx,roty,rotz);
%
% coords: nx3 matrix of 3D coordinates
%
% Early Mod: 01/20/09 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('rotx','var') || isempty(rotx), rotx=0; end;
if ~exist('roty','var') || isempty(roty), roty=0; end;
if ~exist('rotz','var') || isempty(rotz), rotz=0; end;
new_coords = [];

Rx=eye(3);
Ry=eye(3);
Rz=eye(3);
% make same as tksurfer so it is easier to set initial values
if(rotx)
  c = cos(rotx);
  s = sin(rotx);
  Rx(2,2)=c;
  Rx(2,3)=s;
  Rx(3,2)=-s;
  Rx(3,3)=c;
end

if(roty)
  c = cos(roty);
  s = sin(roty);
  Ry(1,1)=c;
  Ry(1,3)=-s;
  Ry(3,1)=s;
  Ry(3,3)=c;
end
  
if(rotz)
  c = cos(rotz);
  s = sin(rotz);
  Rz(1,1)=c;
  Rz(1,2)=s;
  Rz(2,1)=-s;
  Rz(2,2)=c;
end

R = Rz*Ry*Rx;

new_coords = (R*coords')';

