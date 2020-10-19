function [G_x,G_y,G_z]=ts_gain_xyz_split(G_xyz);
%function [G_x,G_y,G_z]=ts_gain_xyz_split(G_xyz);
%
% Split xyz gain matrix into x, y, and z components
%
% Input:
%  G_xyz: gain matrix forward solution with three rows (x,y,z) for each dipole
%
% Output:
%  G_x: forward matrix of sensor amplitudes along x axis
%  G_y: forward matrix of sensor amplitudes along y axis
%  G_z: forward matrix of sensor amplitudes along z axis
%
%
% Created:  10/05/07 by Don Hagler
% Last Mod: 10/05/07 DH
%

G_x = [];
G_y = [];
G_z = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_sensors,num_sources]=size(G_xyz);
n_grid = num_sources/3;
G_x = zeros(num_sensors,n_grid);
G_y = zeros(num_sensors,n_grid);
G_z = zeros(num_sensors,n_grid);
C = eye(3);
for i=1:n_grid;
  % get xyz gain for this dipole
  j = (i-1)*3 + 1;
  k = j + 2;
  tmp_G_xyz = G_xyz(:,j:k);

  % separate components
  G_x(:,i) = tmp_G_xyz(:,1);
  G_y(:,i) = tmp_G_xyz(:,2);
  G_z(:,i) = tmp_G_xyz(:,3);
end

