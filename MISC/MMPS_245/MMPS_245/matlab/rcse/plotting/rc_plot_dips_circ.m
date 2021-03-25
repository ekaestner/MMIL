function rc_plot_dips_circ(prefix,areas2plot,plane);
% usage: rc_plot_dips_circ(prefix,[areas2plot],[plane]);
%
% Purpose: plot orientations of dipoles, arranged in circle,
%   reflecting stimulus locations
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   areas2plot: vector of area numbers
%     {default = [retmap.num_areas:-1:1]}
%   plane: 'cor', 'sag', or 'hor'
%     {default = 'cor'}
%
% Early Mod: 03/04/08 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: use varargin

if ~mmil_check_nargs(nargin,1), return; end;
labelflag = 0;
if ~exist('areas2plot','var'), areas2plot = []; end;
if ~exist('trans','var'), trans = eye(4); end;
if ~exist('plane','var') || isempty(plane), plane = 'cor'; end;
if ~ismember(plane,{'cor','sag','hor'})
  eror('plane must be cor, sag, or hor');
end;
if ~exist('fontname','var') || isempty(fontname), fontname = 'Arial'; end;

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_forward.mat',prefix);
if ~exist(matfile,'file')
  matfile = parms.forward_matfile;
end;
load(matfile);

if isempty(areas2plot)
  areas2plot = fliplr(parms.use_areas);
  if isempty(areas2plot)
    areas2plot = [retmap.num_areas:-1:1];
  end;
end;
areas2plot = areas2plot(find(areas2plot<=retmap.num_areas & areas2plot>=1));
if isempty(areas2plot)
  error('empty areas2plot');
end;


area_colors = ['b','g','r','k','y','c','m','b','g','r','k','y','c','m'];
area_R = [1:length(areas2plot)];
num_locs = retmap.num_locs;
cond_order = retmap.cond_order;
orig_num_locs = retmap.orig_num_locs;
A_offset = 360/(orig_num_locs*2);

ax_pos_init = [0.48,0.48,0.04,0.04];
ax_pos = [0.3,0.3,0.4,0.4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_mri_head = parms.trans;
fprintf('%s: T_mri_head:\n',mfilename);
disp(T_mri_head);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = axes('position',ax_pos_init);
hold on;

for a=1:length(areas2plot)
  area=areas2plot(a);
  color=area_colors(area);
  R=area_R(a);
  coords = zeros(3,num_locs);
  avg_dip = [];
  for theta = 1:num_locs
    cond = cond_order(theta);
    A = A_offset + (cond-1)*360/orig_num_locs;
    v_lh = retmap.M(area,theta).v_lh;
    w_lh = retmap.M(area,theta).w_lh;
    v_rh = retmap.M(area,theta).v_rh;
    w_rh = retmap.M(area,theta).w_rh;

    sumw_lh = sum(w_lh);
    sumw_rh = sum(w_rh);
    if isempty(sumw_lh), sumw_lh=0; end;
    if isempty(sumw_rh), sumw_rh=0; end;
    sumw = sumw_lh + sumw_rh;

    avg_dip_lh = lh_dip_info(:,v_lh)*w_lh;
    avg_dip_rh = rh_dip_info(:,v_rh)*w_rh;
    if isempty(avg_dip_lh), avg_dip_lh = zeros(6,1); end;
    if isempty(avg_dip_rh), avg_dip_rh = zeros(6,1); end;
    avg_dip(:,theta) = (avg_dip_lh + avg_dip_rh)/sumw;

    if(strcmp(plane,'cor'))
      disp_coords(1,theta) = R*cos(A*pi/180);
      disp_coords(2,theta) = 0;
      disp_coords(3,theta) = R*sin(A*pi/180);
    elseif(strcmp(plane,'sag'))
      disp_coords(1,theta) = 0;
      disp_coords(2,theta) = R*cos(A*pi/180);
      disp_coords(3,theta) = R*sin(A*pi/180);
    elseif(strcmp(plane,'hor'))
      disp_coords(1,theta) = R*cos(A*pi/180);
      disp_coords(2,theta) = R*sin(A*pi/180);
      disp_coords(3,theta) = 0;
    end;
  end;
  coords = avg_dip(1:3,:);
  normcoords = coords + avg_dip(4:6,:);

  % apply coordinate transformation (rotation,translation,etc.)
  tmp_coords = [coords/1000; ones(1,length(coords))]; % in meters
  tmp_coords = T_mri_head*tmp_coords;
  coords = tmp_coords(1:3,:)*1000; % back to mm
  tmp_normcoords = [normcoords/1000; ones(1,length(normcoords))];
  tmp_normcoords = T_mri_head*tmp_normcoords;
  normcoords = tmp_normcoords(1:3,:)*1000;
  % get back normal vector
  normcoords = normcoords - coords;
  % get normal relative to display coordinates
  disp_normcoords = disp_coords + normcoords;

  plot3(disp_coords(1,:),disp_coords(2,:),disp_coords(3,:),...
    [color '.'],'MarkerSize',16);
  for theta = 1:num_locs
    line([disp_coords(1,theta) disp_normcoords(1,theta)],...
         [disp_coords(2,theta) disp_normcoords(2,theta)],...
         [disp_coords(3,theta) disp_normcoords(3,theta)],...
         'Color',color,'LineWidth',2);
  end;
end;

t = 1;
R=max(area_R)+0.7;
for theta = 1:num_locs
  cond = cond_order(theta);
  A = A_offset + (cond-1)*360/orig_num_locs;
  if(strcmp(plane,'cor'))
    disp_coords(1,theta) = R*cos(A*pi/180);
    disp_coords(2,theta) = 0;
    disp_coords(3,theta) = R*sin(A*pi/180);
  elseif(strcmp(plane,'sag'))
    disp_coords(1,theta) = 0;
    disp_coords(2,theta) = R*cos(A*pi/180);
    disp_coords(3,theta) = R*sin(A*pi/180);
  elseif(strcmp(plane,'hor'))
    disp_coords(1,theta) = R*cos(A*pi/180);
    disp_coords(2,theta) = R*sin(A*pi/180);
    disp_coords(3,theta) = 0;
  end;
  if labelflag
    set(gca,'FontName',fontname)
    text(disp_coords(1,theta),disp_coords(2,theta),disp_coords(3,theta),...
      sprintf('%d',cond),'FontSize',10,'HorizontalAlignment','Center');
  end;
end;


axis([-max(area_R)-1 max(area_R)+1 -max(area_R)-1 max(area_R)+1 -max(area_R)-1 max(area_R)+1]);
set(ax,'position',ax_pos);
if(strcmp(plane,'cor'))
  view(0,0);
elseif(strcmp(plane,'sag'))
  view(90,0);
elseif(strcmp(plane,'hor'))
  view(0,90);
end;
axis off;
