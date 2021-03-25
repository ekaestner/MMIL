function ts_plot3d_dSPM_sensors_and_sources(prefix)
%function ts_plot3d_dSPM_sensors_and_sources(prefix)
%
% Created:  12/10/08 by Don Hagler
% Last Mod: 12/11/08 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
view_angle = [0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matfile=sprintf('matfiles/%s_avg_data.mat',prefix);
if ~exist(matfile,'file'), error('file %s not found',matfile); end;
load(matfile);

matfile=sprintf('matfiles/%s_parms.mat',prefix);
if ~exist(matfile,'file'), error('file %s not found',matfile); end;
load(matfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: plotting sensor and source coordinates in 3D...\n',mfilename);

T_mri_head = parms.trans;
fprintf('%s: T_mri_head:\n',mfilename);
disp(T_mri_head);

% load sensor info (in head space)
MEG_info = ts_prep_MEG_info(avg_data);
meg_locs = MEG_info.intpnt_loc;

% get EEG locs (in head space)
EEG_info = ts_prep_EEG_info(avg_data);
eeg_locs = EEG_info.intpnt;

% get information about source locations
num_dips_lh=length(parms.lh_dip_info);
num_dips_rh=length(parms.rh_dip_info);
dec_dip_lh = ones(length(parms.lh_dip_info),1);
dec_dip_rh = ones(length(parms.rh_dip_info),1);

% get indices of select dipoles (non-zero in dec_dip_lh and dec_dip_rh)
id_lh_dip=find(dec_dip_lh==1);
id_rh_dip=find(dec_dip_rh==1);
% combine info for select dipoles
grid_mri=[parms.lh_dip_info(1:3,id_lh_dip)';parms.rh_dip_info(1:3,id_rh_dip)'];
n_grid=size(grid_mri,1);
fprintf('%s: total number of dipoles in forward = %d\n',mfilename,n_grid);
% apply coordinate transformation to get to sensor-space
grid_head=T_mri_head*[grid_mri'/1000;ones(1,n_grid)]; % in meters
grid_head=grid_head(1:3,:)'*1000; % back to mm and reshape
source_locs = grid_head;

hold on
plot3d(source_locs,'b*');
plot3d(meg_locs*1000,'ro')
plot3d(eeg_locs*1000,'gh')
xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
view(view_angle);
hold off
drawnow


%print('-dtiff',sprintf('%s-sensors_and_sources.tif',prefix));

