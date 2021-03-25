function rc_plot_MEG_sensor_dips(prefix,areas2plot,view_angle);
%function rc_plot_MEG_sensor_dips(prefix,[areas2plot],[view_angle]);
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   areas2plot: vector of area numbers {default: [1:num_areas]}
%   view_angle: [azimuth elevation]
%
% Early Mod: 05/11/07 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('areas2plot','var'), areas2plot = []; end;
if ~exist('view_angle','var'), view_angle = []; end;
if isempty(view_angle), view_angle = [0 0]; end;

area_colors = ['r','b','g','k','y','c','m','r','b','g','k','y','c','m'];
sensor_marker = 'ro';
dip_scale_fact = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_ret_mapping.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_dip_info.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

if isempty(areas2plot)
  areas2plot = [1:retmap.num_areas];
end;
areas2plot = areas2plot(find(areas2plot<=retmap.num_areas & areas2plot>=1));
if isempty(areas2plot)
  error('no valid areas in areas2plot');
end;

T_mri_head = parms.trans;
fprintf('%s: T_mri_head:\n',mfilename);
disp(T_mri_head);

num_locs = retmap.num_locs;
cond_order = retmap.cond_order;
orig_num_locs = retmap.orig_num_locs;

% load sensor info
MEG_info = ts_prep_MEG_info(avg_data);
sensor_locs = MEG_info.intpnt_loc;

hold on
rc_plot_3d(sensor_locs*1000,sensor_marker)

for a=1:length(areas2plot)
  area=areas2plot(a);
  color=area_colors(area);
  coords = zeros(3,num_locs);
  avg_dip = [];
  for theta = 1:num_locs
    cond = cond_order(theta);
    v_lh = retmap.M(area,theta).v_lh;
    w_lh = retmap.M(area,theta).w_lh;
    v_rh = retmap.M(area,theta).v_rh;
    w_rh = retmap.M(area,theta).w_rh;

    sumw_lh = sum(w_lh);
    sumw_rh = sum(w_rh);
    if isempty(sumw_lh), sumw_lh=0; end;
    if isempty(sumw_rh), sumw_rh=0; end;
    sumw = sumw_lh + sumw_rh;

    avg_dip_lh = dip_info_lh(:,v_lh)*w_lh;
    avg_dip_rh = dip_info_rh(:,v_rh)*w_rh;
    if isempty(avg_dip_lh), avg_dip_lh = zeros(6,1); end;
    if isempty(avg_dip_rh), avg_dip_rh = zeros(6,1); end;
    avg_dip(:,theta) = (avg_dip_lh + avg_dip_rh)/sumw;
  end;
  coords = avg_dip(1:3,:);
  normcoords = coords + dip_scale_fact*avg_dip(4:6,:);
  % apply coordinate transformation (rotation,translation,etc.)
  tmp_coords = [coords/1000; ones(1,length(coords))]; % in meters
  tmp_coords = T_mri_head*tmp_coords;
  coords = tmp_coords(1:3,:)*1000; % back to mm
  tmp_normcoords = [normcoords/1000; ones(1,length(normcoords))];
  tmp_normcoords = T_mri_head*tmp_normcoords;
  normcoords = tmp_normcoords(1:3,:)*1000;
  % plot origins
  plot3(coords(1,:),coords(2,:),coords(3,:),...
    [color '.'],'MarkerSize',16);
  % plot normal vectors and labels
  for theta = 1:num_locs
    % left hemi
    line([coords(1,theta) normcoords(1,theta)],...
         [coords(2,theta) normcoords(2,theta)],...
         [coords(3,theta) normcoords(3,theta)],...
         'Color',color,'LineWidth',2);
    text(coords(1,theta),coords(2,theta),coords(3,theta),sprintf('%d',theta),...
      'FontSize',10);
  end;
end;

xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
view(view_angle);
hold off
drawnow


