function rc_plot_MEG_sensors_sources(prefix,alldips_flag)
%function rc_plot_MEG_sensors_sources(prefix,[alldips_flag])
%
% Purpose: plot locations of MEG sensors relative to sources
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   alldips_flag: [0|1] show all possible dipoles,
%               otherwise, only dipoles in forward
%
% Early Mod: 05/11/07 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

% todo: varargin, allow rootdir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('alldips_flag','var')  || isempty(alldips_flag), alldips_flag = 1; end;
view_angle = [0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_ret_mapping.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_dip_info',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

T_mri_head = parms.trans;
fprintf('%s: T_mri_head:\n',mfilename);
disp(T_mri_head);

% load sensor info
MEG_info = ts_prep_MEG_info(avg_data);
sensor_locs = MEG_info.intpnt_loc;

% get information about source locations
num_dips_lh=length(dip_info_lh);
num_dips_rh=length(dip_info_rh);
if alldips_flag
  % display all dipoles
  dec_dip_lh = ones(length(dip_info_lh),1);
  dec_dip_rh = ones(length(dip_info_rh),1);
else
  % create dec_dip_lh and dec_dip_rh from retmap.uniq_verts_lh and retmap.uniq_verts_rh
  dec_dip_lh = zeros(length(dip_info_lh),1);
  dec_dip_lh(retmap.uniq_verts_lh) = 1;
  dec_dip_rh = zeros(length(dip_info_rh),1);
  dec_dip_rh(retmap.uniq_verts_rh) = 1;
end;

% get indices of select dipoles (non-zero in dec_dip_lh and dec_dip_rh)
id_lh_dip=find(dec_dip_lh==1);
id_rh_dip=find(dec_dip_rh==1);
% combine info for select dipoles
grid_mri=[dip_info_lh(1:3,id_lh_dip)';dip_info_rh(1:3,id_rh_dip)'];
n_grid=size(grid_mri,1);
fprintf('%s: total number of dipoles plotted = %d\n',mfilename,n_grid);
% apply coordinate transformation to get to sensor-space
grid_head=T_mri_head*[grid_mri'/1000;ones(1,n_grid)]; % in meters
grid_head=grid_head(1:3,:)'*1000; % back to mm and reshape
grid_head_cen=grid_head-ones(n_grid,1)*parms.cen_sph;
source_locs = grid_head;

hold on
rc_plot_3d(source_locs,'b*');
rc_plot_3d(sensor_locs*1000,'ro')
xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
view(view_angle);
hold off
drawnow


