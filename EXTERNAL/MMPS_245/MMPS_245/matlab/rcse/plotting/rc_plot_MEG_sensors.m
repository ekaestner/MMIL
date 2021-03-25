function rc_plot_MEG_sensors(prefix,view_angle)
%function rc_plot_MEG_sensors(prefix,[view_angle])
%
% Purpose: plot 3D locations of MEG sensors
%
% Required Input:
%   prefix: proc or RCSE prefix
%
% Optional Input:
%   view_angle: azimuth and elevation angles
%     {default = [0 0]}
%
% Early Mod:  03/11/07 by Don Hagler
% Last Mod:   02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('view_angle','var') || isempty(view_angle), view_angle = [0 0]; end;

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

% load sensor info
MEG_info = ts_prep_MEG_info(avg_data);
sensor_locs = MEG_info.intpnt_loc;

hold on
rc_plot_3d(sensor_locs*1000,'ro')
xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
view(view_angle);
hold off
drawnow


