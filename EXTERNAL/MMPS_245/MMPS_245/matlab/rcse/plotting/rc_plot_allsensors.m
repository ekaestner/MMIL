function rc_plot_allsensors(prefix)
%function rc_plot_allsensors(prefix)
%
% Purpose: plot 3D locations of MEG and EEG sensors
%
% Required Input:
%   prefix: proc or RCSE prefix
%
% Early Mod:  05/11/07 by Don Hagler
% Last Mod:   02/20/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
view_angle = [0 0];
if isempty(prefix), error('empty prefix'); end;

matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

% load sensor info (in head space)
MEG_info = ts_prep_MEG_info(avg_data);
meg_locs = MEG_info.intpnt_loc;

% get EEG locs (in head space)
EEG_info = ts_prep_EEG_info(avg_data);
eeg_locs = EEG_info.intpnt;

hold on
rc_plot_3d(meg_locs*1000,'ro')
rc_plot_3d(eeg_locs*1000,'gh')
xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
view(view_angle);
hold off
drawnow


