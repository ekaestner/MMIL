function [G_norm,G_tang1,G_tang2]=ts_gain_xyz2norm(G_xyz,...
  lh_dip_info,rh_dip_info,lh_dec_dips,rh_dec_dips,T_mri2head);
%function [G_norm,G_tang1,G_tang2]=ts_gain_xyz2norm(G_xyz,...
%  lh_dip_info,rh_dip_info,lh_dec_dips,rh_dec_dips,T_mri2head);
%
% calculate normal and tangential components from xyz gain matrix
%
% Input:
%  G_xyz: gain matrix forward solution with three rows (x,y,z) for each dipole
%  lh_dip_info: matrix containing locations and orientations for each dipole
%               in left hemisphere (6 rows, column for each vertex)
%               generated with ts_read_dip_file.m
%  rh_dip_info: matrix containing locations and orientations for each dipole
%               in right hemisphere (6 rows, column for each vertex)
%  lh_dec_dips: vector containing zeros and ones specifying which left
%              hemisphere dipoles to include
%              can be generated with ts_read_dec_file.m
%  rh_dec_dips: vector containing zeros and ones specifying which right
%              hemisphere dipoles to include
%  T_mri2head: transformation matrix specifying registration between
%              mri (freesurfer brain space) and head (MEG/EEG sensor space)
%              can be generated with loadtrans function (fiff access)
%              or with ts_pointreg
%
% Output:
%  G_norm: forward matrix of sensor amplitudes for each decimated dipole
%          (along normal vector)
%  G_tang1: forward matrix of sensor amplitudes for each decimated dipole
%          (along first of two orthogonal tangential componenets)
%  G_tang2: forward matrix of sensor amplitudes for each decimated dipole
%          (along second of two orthogonal tangential componenets)
%
%
% Created:  08/01/06 by Don Hagler
% Last Mod: 12/28/13 by Don Hagler
%

G_norm = [];
G_tang1 = [];
G_tang2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices of select dipoles (non-zero in lh_dec_dips and rh_dec_dips)
id_lh_dip=find(lh_dec_dips==1);
id_rh_dip=find(rh_dec_dips==1);

% combine info for select dipoles
grid_mri=[lh_dip_info(1:3,id_lh_dip)';rh_dip_info(1:3,id_rh_dip)'];
n_grid=size(grid_mri,1);
%fprintf('%s: total number of dipoles in forward = %d\n',mfilename,n_grid);

% get normal vector info
norm_mri=[lh_dip_info(4:6,id_lh_dip)';rh_dip_info(4:6,id_rh_dip)'];

% apply mri2head coordinate transformation to get to "head"-space
grid_head=T_mri2head*[grid_mri'/1000;ones(1,n_grid)]; % mm to meters
grid_head=grid_head(1:3,:)'; % reshape, still in meters
source_locs = grid_head;

% apply mri2head rotation for normal vectors
T_mri2head_rot = T_mri2head(1:3,1:3);
norm_head=(T_mri2head_rot*norm_mri')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate norm and tangential components
%fprintf('%s: calculating gain matrix for normal and tangential components\n',mfilename);
[num_sensors,num_sources]=size(G_xyz);
G_norm = zeros(num_sensors,n_grid);
G_tang1 = zeros(num_sensors,n_grid);
G_tang2 = zeros(num_sensors,n_grid);
C = eye(3);
for i=1:n_grid;
  % find tangential components (2 vectors orthogonal to normal vector)
  N = norm_head(i,:);
  dp = abs(N*C); % dot product with each of three orthogonal vectors
  mindp = find(dp==min(dp)); % min(dp) is farthest from parallel to normal
  if length(mindp)>1
    mindp = mindp(1);
  end;
  T1 = C(mindp,:);  % relatively arbitrary, non-parallel vector to start
  T1 = cross(N,T1); % orthogonal to N (and T1)
  T2 = cross(N,T1); % orthogonal to N and T1

  % get xyz gain for this dipole
  j = (i-1)*3 + 1;
  k = j + 2;
  tmp_G_xyz = G_xyz(:,j:k);

  % caculate normal and tangential components
  G_norm(:,i) = tmp_G_xyz*N';
  G_tang1(:,i) = tmp_G_xyz*T1';
  G_tang2(:,i) = tmp_G_xyz*T2';
end

