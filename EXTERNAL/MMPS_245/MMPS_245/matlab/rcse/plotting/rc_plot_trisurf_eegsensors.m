function rc_plot_trisurf_eegsensors(trifile,prefix)
%function rc_plot_trisurf_eegsensors(trifile,prefix)
%
% Purpose: plot EEG sensor locations relative to outer scalp
%
% Required Input:
%   trifile: name of freesurfer trifile (e.g. outer scalp)
%   prefix: RCSE prefix
%
% Early Mod: 01/08/10 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: rootdir

if ~mmil_check_nargs(nargin,1), return; end;

colorval = 10;
face_alpha = 0.6;
edge_alpha = 0.1;
view_angle = [0 0];

surf = fs_read_trisurf(surffile);

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

% get EEG locs (in head space)
T_mri2head = parms.trans;
T_head2mri = inv(T_mri2head);
T = T_head2mri;
chans = find(strcmp('eeg',{avg_data.sensor_info.typestring}));
eeg_locs = zeros(length(chans),3);
for c=1:length(chans)
  k=chans(c);
  loc = avg_data.sensor_info(k).loc;
  % apply coordinate transform
  loc = T*loc;
  eeg_locs(c,1) = loc(1,4);
  eeg_locs(c,2) = loc(2,4);
  eeg_locs(c,3) = loc(3,4);
end;


% plot surf
clf;
hold on;
colors = colorval*ones(surf.nverts,1);
h=trisurf(surf.faces,...
          surf.vertices(:,1),...
          surf.vertices(:,2),...
          surf.vertices(:,3),...
          colors);
set(h,'EdgeAlpha',edge_alpha,'FaceAlpha',face_alpha);
colormap hsv;
caxis([0 100]);

% plot eeg sensors
rc_plot_3d(eeg_locs*1000,'gh')

xlabel('x');
ylabel('y');
zlabel('z');
view(view_angle);
hold off

return;

