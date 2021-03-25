function rc_plot_dips(prefix,varargin)
%function rc_plot_dips(prefix,[options])
%
% Purpose: plot location and orientation of RCSE dipoles overlayed on MRI
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Parameters:
%  'fname_conds': file name of csv file containing condition information
%    {default = 'cond_info.csv'}
%  'retfit_stem': file stem for refit results mat file
%    {default = 'retfit/retfit_results'}
%  'r_offset': vector of r_offsets
%    {default = 0}
%  'th_offset': vector of th_offsets
%    {default = 0}
%  'r_sig': eccentricity width of neighborhood around dipoles (deg. vis. ang.)
%    {default = 2}
%  'th_sig': polar angle width of neighborhood around dipoles (degrees)
%    {default = 5}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'plane': 'cor', 'sag', or 'hor'
%    {default = 'cor'}
%  'slicelist': vector of integers between 0 and 256
%              if 0, show all dipoles
%              otherwise, show only dipoles in those slices
%    {default = 0}
%  'plot_outdir': name of directory in which to create output tif files
%    {default = 'plots_dips'}
%
% Created:  04/14/09 by Don Hagler
% Last Mod: 02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix',prefix,[],...
  'fname_conds','cond_info.csv',[],...
  'retfit_stem','retfit/retfit_results',[],...
  'r_offset',0,[-100,100],...
  'th_offset',0,[-180,180],...
  'r_sig',2,[0,100],...
  'th_sig',5,[0,360],...
  'r_max',12.5,[0,Inf],...
  'plot_outdir','plots_dips',[],...
  'slicelist',0,[0,256],...
  'plane','cor',{'sag','cor','hor'},...
... % undocumented
  'area_names',{'v1','v2','v3'},[],...
  'area_colors',{'b','g','r','k','y','c','m','b','g','r','k','y','c','m'},[],...
  'w_thresh',0.5,[0,100],...
  'multi_area_flag',false,[false true],...
  'visible_flag',true,[false true],...
  'label_flag',true,[false true],...
  'dip_scale_fact',8,[],...
  'label_font_size',14,[],...
  'linewidth',3,[],...
  'crop_flag',true,[false true],...
});

switch parms.plane
  case 'cor'
    viewplane = 2;
  case 'sag'
    viewplane = 1;
  case 'hor'
    viewplane = 3;
end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_conds,'file')
  error('file %s not found',parms.fname_conds);
end;
cond_info = read_cond_info(parms.fname_conds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load parameters from retinv
matfile = sprintf('matfiles/%s_parms.mat',parms.prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
tmp_parms = parms;
load(matfile);
ri_parms = parms;
parms = tmp_parms;

% load forward info from retinv
forward = [];
if isfield(ri_parms,'forward_matfile') & ~isempty(ri_parms.forward_matfile)
  matfile = ri_parms.forward_matfile;
else
  matfile = sprintf('matfiles/%s_forward.mat',parms.prefix);
end;
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);
lh_dip_info = forward.lh_dip_info;
rh_dip_info = forward.rh_dip_info;

if parms.crop_flag
  switch parms.plane
    case 'sag'
      x0 = 1;
      x1 = 256;
      y0 = 26;
      y1 = 80;
      z0 = 101;
      z1 = 155;
    case 'cor'
      x0 = 88;
      x1 = 162;
      y0 = 1;
      y1 = 256;
      z0 = 92;
      z1 = 166;
    case 'hor'
      x0 = 88;
      x1 = 152;
      y0 = 26;
      y1 = 90;
      z0 = 1;
      z1 = 256;
  end;
else
  % no cropping
  x0 = 1;
  x1 = 256;
  y0 = 1;
  y1 = 256;
  z0 = 1;
  z1 = 256;
end;
switch parms.plane
  case 'cor'
    view_angle = [0 0];
  case 'sag'
    view_angle = [90 0];
  case 'hor'
    view_angle = [0 90];
end;
ax_pos = [0.05,0.05,0.9,0.9];
% need offset values to switch between voxels and RAS coordinates
% most likely freesurfer coordinate offset
xoff = 129;
yoff = 129;
zoff = 129;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load retfit files
retfit_results = [];
hemilist = {'lh','rh'};
for h=1:length(hemilist)
  hemi = hemilist{h};
  fname = sprintf('%s-%s.mat',parms.retfit_stem,hemi);
  if ~exist(fname,'file'), error('file %s not found',fname); end;
  if strcmp(hemi,'lh')
    clear fit_data area_masks fit_parms;
    load(fname);
    retfit_results.lh_fit_data = fit_data;
    retfit_results.lh_area_masks = area_masks;
  else
    clear fit_data area_masks fit_parms;
    load(fname);
    retfit_results.rh_fit_data = fit_data;
    retfit_results.rh_area_masks = area_masks;
  end;
end;
parms.lh_verts = find(retfit_results.lh_fit_data.th~=0 |...
                      retfit_results.lh_fit_data.r~=0);
parms.rh_verts = find(retfit_results.rh_fit_data.th~=0 |...
                      retfit_results.rh_fit_data.r~=0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get cortical surface vertices for each stimulus location
retmap = rc_define_retmap_from_retfit(...
  retfit_results,...
  'area_names',parms.area_names,...
  'fname_conds',parms.fname_conds,...
  'r_sig',parms.r_sig,...
  'th_sig',parms.th_sig,...
  'r_offset',parms.r_offset,...
  'th_offset',parms.th_offset,...
  'r_max',parms.r_max,...
  'w_thresh',parms.w_thresh);

tmp_contvec = cell2mat({cond_info.contrast});
ind_conds = find(tmp_contvec);
nconds = length(ind_conds);

if parms.multi_area_flag
  figure; clf; ax_dips = axes('position',ax_pos); hold on;
end;

for a=1:length(retmap.areas)
  area_name = upper(retmap.areas(a).name);
  color=parms.area_colors{a};
  avg_dips = zeros(6,nconds);
  disp_coords = zeros(3,nconds);
  tmp_r_max = 0;
  for i=1:nconds
    j = ind_conds(i);

    % stimulus location
    x = cond_info(j).ecc*cos(cond_info(j).theta*pi/180);
    y = cond_info(j).ecc*sin(cond_info(j).theta*pi/180);
    if cond_info(j).ecc>tmp_r_max
      tmp_r_max = cond_info(j).ecc;
    end;

    % vertex numbers -- indices to dip_info
    v_lh = retmap.areas(a).verts(j).v_lh;
    v_rh = retmap.areas(a).verts(j).v_rh;
    w_lh = retmap.areas(a).verts(j).w_lh;
    w_rh = retmap.areas(a).verts(j).w_rh;

    % get weightings for each vertex
    sumw_lh = sum(w_lh);
    sumw_rh = sum(w_rh);
    if isempty(sumw_lh), sumw_lh=0; end;
    if isempty(sumw_rh), sumw_rh=0; end;
    sumw = sumw_lh + sumw_rh;

    % do a weighted average of the dipole location and orientation
    avg_dip_lh = lh_dip_info(:,v_lh)*w_lh;
    avg_dip_rh = rh_dip_info(:,v_rh)*w_rh;
    if isempty(avg_dip_lh), avg_dip_lh = zeros(6,1); end;
    if isempty(avg_dip_rh), avg_dip_rh = zeros(6,1); end;
    avg_dip(:,i) = (avg_dip_lh + avg_dip_rh)/sumw;
  end;
  coords = avg_dip(1:3,:);
  normcoords = coords + parms.dip_scale_fact*avg_dip(4:6,:);

  if parms.slicelist==0
    dips2plot = [1:nconds];
  else
    dips2plot = [];
    for i=1:length(parms.slicelist)
      slicenum = parms.slicelist(i);
      switch parms.plane
        case 'sag'
          % treat x differently because of left-right flip
          a0 = xoff-slicenum-0.5;
          a1 = xoff-slicenum+0.5;
        case 'cor'
          a0 = slicenum-yoff-0.5;
          a1 = slicenum-yoff+0.5;
        case 'hor'
          a0 = slicenum-zoff-0.5;
          a1 = slicenum-zoff+0.5;
      end;    
      tmp_dips2plot = find(coords(viewplane,:)>a0 &...
                           coords(viewplane,:)<a1);
      dips2plot = union(dips2plot,tmp_dips2plot);
    end;
  end;
  tmp_nconds = length(dips2plot);

  if ~parms.multi_area_flag
    figure(a); clf; ax_dips = axes('position',ax_pos); hold on;
  end;
  plot3(coords(1,dips2plot),coords(2,dips2plot),coords(3,dips2plot),...
    [color '.'],'MarkerSize',16);
  axis([xoff-x1 xoff-x0 y0-yoff y1-yoff z0-zoff z1-zoff]);
  view(view_angle);
  t = 1;
  for t=1:tmp_nconds
    d = dips2plot(t);
    line([coords(1,d) normcoords(1,d)],...
         [coords(2,d) normcoords(2,d)],...
         [coords(3,d) normcoords(3,d)],...
         'Color',color,'LineWidth',parms.linewidth);
    if parms.label_flag
      text(coords(1,d),coords(2,d),coords(3,d),...
        sprintf('%d',d),...
        'FontSize',parms.label_font_size,...
        'HorizontalAlignment','Center','Color','w');
    end;
  end;
  xlabel('left-right');
  zlabel('inf-sup');
  ylabel('pos-ant');
  axis off;
  hold off;
  suptitle(...
    sprintf('%s dipoles at MRI locations, %s view',area_name,parms.plane));
  if ~exist(parms.plot_outdir,'dir'), mkdir(parms.plot_outdir); end;
  fname_out = sprintf('%s/dips-area%s-%s',...
    parms.plot_outdir,area_name,parms.plane);
  if ~parms.visible_flag
    set(gcf,'Visible','off');
  end;
  print('-dtiff',fname_out);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
