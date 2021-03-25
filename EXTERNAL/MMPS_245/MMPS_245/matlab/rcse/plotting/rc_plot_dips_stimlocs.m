function rc_plot_dips_stimlocs(prefix,varargin)
%function rc_plot_dips_stimlocs(prefix,[options])
%
% Purpose: plot dipole orientations for each stimulus location
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Parameters:
%  'cond_info' - struct array containing condition information
%    see read_cond_info
%    If supplied, fname_conds will be ignored
%    {default: []}
%  'fname_conds' - full path of csv file containing condition information
%    If cond_info and fname_conds are empty, will attempt to get cond_info
%    from stored mat files
%    {default: []}
%  'retfit_stem' - file stem for refit results mat file
%    {default = 'retfit/retfit_results'}
%  'r_offset' - vector of r_offsets
%    {default = 0}
%  'th_offset' - vector of th_offsets
%    {default = 0}
%  'r_sig' - eccentricity width of neighborhood around dipoles (deg. vis. ang.)
%    {default = 2}
%  'th_sig' - polar angle width of neighborhood around dipoles (degrees)
%    {default = 5}
%  'r_max' - maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'plot_outdir' - name of directory in which to create output tif files
%    {default = 'plots_dips_stimlocs'}
%
% Created:  04/14/09 by Don Hagler
% Last Mod: 02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix',prefix,[],...
  'cond_info',[],[],...
  'fname_conds',[],[],...
  'retfit_stem','retfit/retfit_results',[],...
  'r_offset',0,[-100,100],...
  'th_offset',0,[-180,180],...
  'r_sig',2,[0,100],...
  'th_sig',5,[0,360],...
  'r_max',12.5,[0,Inf],...
  'plot_outdir','plots_dips_stimlocs',[],...
... % undocumented
  'area_names',{'v1','v2','v3'},[],...
  'area_colors',{'b','g','r','k','y','c','m','b','g','r','k','y','c','m'},[],...
  'w_thresh',0.5,[0,100],...
  'multi_area_flag',false,[false true],...
  'visible_flag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.cond_info)
  if ~isempty(parms.fname_conds)
    if ~exist(parms.fname_conds,'file')
      error('file %s not found',parms.fname_conds);
    end;
    parms.cond_info = read_cond_info(parms.fname_conds);
  else
    parms.cond_info = get_cond_info_from_prefix(parms.prefix);
  end;
  if isempty(parms.cond_info)
    error('failed to load cond_info');
  end;
end;

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
  'cond_info',parms.cond_info,...
  'r_sig',parms.r_sig,...
  'th_sig',parms.th_sig,...
  'r_offset',parms.r_offset,...
  'th_offset',parms.th_offset,...
  'r_max',parms.r_max,...
  'w_thresh',parms.w_thresh);

tmp_contvec = cell2mat({parms.cond_info.contrast});
ind_conds = find(tmp_contvec);
nconds = length(ind_conds);

if parms.multi_area_flag
  figure; clf; hold on;
end;

for a=ri_parms.use_areas
  area_name = upper(retmap.areas(a).name);
  color=parms.area_colors{a};
  avg_dips = zeros(6,nconds);
  disp_coords = zeros(3,nconds);
  tmp_r_max = 0;
  for i=1:nconds
    j = ind_conds(i);

    % stimulus location
    x = parms.cond_info(j).ecc*cos(parms.cond_info(j).theta*pi/180);
    y = parms.cond_info(j).ecc*sin(parms.cond_info(j).theta*pi/180);
    if parms.cond_info(j).ecc>tmp_r_max
      tmp_r_max = parms.cond_info(j).ecc;
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

    % disply coordinates according to stimulus location
    disp_coords(1,i) = x;
    disp_coords(2,i) = y;
    disp_coords(3,i) = 0;
%    plot(x,y,'b.');
%    text(x,y,num2str(parms.cond_info(i).event_code));
  end;
  coords = avg_dip(1:3,:);
  normcoords = avg_dip(4:6,:);

  % get normal relative to display coordinates
  disp_normcoords = disp_coords + normcoords;

  if ~parms.multi_area_flag
    figure(a); clf; hold on;
  end;
  plot3(disp_coords(1,:),disp_coords(2,:),disp_coords(3,:),...
    [color '.'],'MarkerSize',16);
  for i = 1:nconds
    line([disp_coords(1,i) disp_normcoords(1,i)],...
         [disp_coords(2,i) disp_normcoords(2,i)],...
         [disp_coords(3,i) disp_normcoords(3,i)],...
         'Color',color,'LineWidth',2);
  end;

  % plot axes
  line([-tmp_r_max tmp_r_max],[0 0],[0 0],...
         'Color','k','LineWidth',1);
  xlabel('x');

  line([0 0],[-tmp_r_max tmp_r_max],[0 0],...
         'Color','k','LineWidth',1);
  ylabel('y');

  line([0 0],[0 0],[-tmp_r_max tmp_r_max],...
         'Color','k','LineWidth',1);
  zlabel('z');

  suptitle(...
    sprintf('%s dipoles at stimulus locations',area_name));

  if ~exist(parms.plot_outdir,'dir'), mkdir(parms.plot_outdir); end;
  fname_out = sprintf('%s/dips-stimlocs-area%s',parms.plot_outdir,area_name);
  if ~parms.visible_flag
    set(gcf,'Visible','off');
  end;
  print('-dtiff',fname_out);
end;


