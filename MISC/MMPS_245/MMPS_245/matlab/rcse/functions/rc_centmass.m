function [x,y,z] = rc_centmass(X,Y,Z,V);
%function [x,y,z] = rc_centmass(X,[Y],[Z],[V]);
%
% input:
%  X,Y,Z: vectors of x,y,z coordinates
%  V: vector of values for each point
%
% output:
%  x,y,z: x,y,z coordinates of center of mass
%                  
% Last Mod:  06/12/10 by Don Hagler
% Early Mod: 02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

npoints = length(X);

if ~exist('Y','var')
  Y = zeros(size(X));
end;
if ~exist('Z','var')
  Z = zeros(size(X));
end;
if ~exist('V','var')
  V = ones(size(X));
end;

if length(Y) ~= npoints | length(Z) ~= npoints  | length(V) ~= npoints
  error('X, Y, Z, and V must have same number of points');
end;

sumV = sum(V);
x = sum(X.*V)/sumV;
y = sum(Y.*V)/sumV;
z = sum(Z.*V)/sumV;

return;

