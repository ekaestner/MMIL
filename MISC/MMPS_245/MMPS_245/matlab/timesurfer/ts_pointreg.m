function [T,min_error] = ts_pointreg(hptsfile,surffile,varargin)
%function [T,min_error] = ts_pointreg(hptsfile,surffile,varargin)
%
% Purpose: Manual and automatic registration between MEG and MRI reference
%   frames.  Plots head points in relation to scalp surface vertices and
%   allows for rotation and translation in 3 axes
%
% Usage:
%   [T,min_error] = ts_pointreg(hptsfile,surffile,'key1', value1,...);
%
% Required input:
%  hptsfile: text file containing head point coordinates
%  surffile: scalp surface file in "tri" format
%
% Optional parameters:
%  'iskullfile': inner skull surface file (tri)
%    {default = []}
%  'oskullfile': outer skull surface file (tri)
%    {default = []}
%  'intransfile': file name for input 4x4 transformation matrix
%    if empty or not found, use identity matrix
%    {default = []}
%  'outtransfile': file name for output 4x4 transformation matrix
%    {default = ts_pointreg.trans}
%
% Colors and transparency:
%  'scalp_alpha': scalp surface transparency (0: 1)
%    {default = 0.5} 
%  'iskull_alpha': inner skull surface transparency (0: 1)
%    {default = 1.0} 
%  'oskull_alpha': outer skull surface transparency (0: 1)
%    {default = 0.6} 
%  'scalp_color': scalp color value (0: 1)
%    {default = 1.0}
%  'iskull_color': inner skull color value (0: 1)
%    {default = 0.0}
%  'oskull_color': outer skull color value (0: 1)
%    {default = 0.5}
%  'lights_flag': [0|1] create spot lights for better contrast
%    {default = 1}
%  'colormap': color map to use for surface colors
%    {default = 'bone'}
%  'card_point_color': color string for cardinal head points
%    {default = 'r'}
%  'point_color': color string for extra head points
%    {default = 'b'}
%  'point_size': size of head point markers
%    {default = 20}
%
% More options:
%  'fitflag': [0|1|2] automated registration to minimize distance between
%     head points and nearest surfce vertices
%     1 = use fmincon, 2 = use random search
%    {default = 0}
%  'T_init': 4x4 transformation matrix
%    {default = use identity matrix}
%  'oldtrans_flag': [0|1] whether intransfile is old format
%    with mm units for translations
%    {default = 0}
%  'trans_type': ['head2mri' or 'mri2head']
%     if 'head2mri': input and output transformation matrices are head2mri
%       i.e. head  points registered to surface (derived from MRI)
%     if 'mri2head': input and output transformation matrices are mri2head
%    {default = 'mri2head'}
%  'cost_include_percentile': head points with minimum distance greater
%    than this percentile will not contribute to error estimate
%    {default = 100}
%  'verbose': [0|1] whether to output messages to terminal
%    {default = 0} 
%
% example of hpts file format:
% cardinal 001 -0.070086 0.000000 -0.000000
% cardinal 002 -0.000000 0.099511 -0.000000
% cardinal 003 0.064181 0.000000 -0.000000
% hpi      001 0.056898 -0.018247 -0.036789
% hpi      002 0.038970 0.098632 0.046132
% hpi      003 -0.018337 0.102662 0.058440
% hpi      004 -0.060982 -0.012347 -0.029710
% extra    001 0.022850 0.103267 0.048321
% extra    002 0.019972 0.101263 0.061493
% extra    003 0.018250 0.096633 0.074459
% extra    004 0.017615 0.089393 0.085852
%
% the hptsfile can be obtained by running fiff2ascii -p on Neuromag fif file
%
% example of trans file format:
% 1.000000 0.000000 0.000000 0.000000 
% 0.000000 1.000000 0.000000 0.000000 
% 0.000000 0.000000 1.000000 0.000000 
% 0.000000 0.000000 0.000000 1.000000 
%
% units of the transformation matrix should be meters
%   if old .trans file created with tkmedit is used, set
%   oldtrans_flag = 1
%
% Created:  05/02/06 by Don Hagler
% Last Mod: 12/07/11 by Don Hagler
%

%% todo: option to plot EEG locations (in different colors) plus ref

%% todo: why so slow when rotating head compared to simple figure?

%% todo: manually click on points to remove
%% todo: allow fif file as input and use hpipoints.m

%  'excl_point_color': color string for excluded head points
%     {default = 'k'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms( varargin, {...
  'fitflag',0,[0 1 2],...
  'niters',100,[0,Inf],... % fitflag = 2 only (random search)
  'iskullfile',[],[],...
  'oskullfile',[],[],...
  'intransfile',[],[],...
  'outtransfile','pointreg.trans',[],...
  'oldtrans_flag',false,[false true],...
  'T_init',[],[],...
  'trans_type','mri2head',{'head2mri','mri2head'},...
  'scalp_alpha',0.5,[0,1],...
  'iskull_alpha',1,[0,1],...
  'oskull_alpha',0.6,[0,1],...
  'scalp_color',0.7,[0,1],...
  'iskull_color',0.2,[0,1],...
  'oskull_color',0.5,[0,1],...
  'colormap','gray',[],...
  'lights_flag',true,[false true],...
  'verbose',false,[false true],...
  'plot_flag',2,[0,1,2],... % 0=no plot, 1=plot, 2=plot + controls
...
  'trans_stepsize',0.1,[],...
  'rot_stepsize',1,[],...
  'trans_stepsize_manual',2,[],...
  'rot_stepsize_manual',5,[],...
  'min_error',Inf,[],...
  'cost_include_percentile',100,[50,100],...
  'init_lpa',[-0.05,0,0],[],...
  'init_rpa',[0.05,0,0],[],...
  'init_nasion',[0,0.1,0],[],...
  'fit_card_flag',false,[false true],...
  'init_card_flag',false,[false true],...
...
  'view',[190,20],[-Inf,Inf],... % azimuth and elevation
  'point_color','b',[],...
  'excl_point_color','k',[],...
  'card_point_color','r',[],...
  'eeg_point_color','g',[],...
  'eeg_ref_color','c',[],...
  'point_size',20,[],...
  'axes_position',[0.2 0.01 0.6  0.99],[],...
  'fig_color',[0.6 0.6 0.6],[],...
  'point_num_flag',false,[false true],...
...
  'hptsfile',hptsfile,[],...
  'surffile',surffile,[],...
  'best_transvec',zeros(3,1),[],...
  'best_rotvec',zeros(3,1),[],...
  'best_T',[],[],...
...
  'fontsize',12,[],...
  'xyz_list',{'x','y','z'},[],...
  'transrot_list',{'trans','rot'},[],...
});

min_error = parms.min_error;
T = parms.best_T;

%  'axes_position',[0.2 0.1 0.5  0.5],[],...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load tri file(s)
parms = load_surf(parms);

% load head points file
parms = load_hpts(parms);

% initialize transformation matrix
parms = init_trans(parms);
if isempty(parms.best_T), parms.best_T = parms.T; end;

% prepare figure and axes
if parms.plot_flag>0, parms = setup_figure(parms); end;

% setup uicontrols
if parms.plot_flag>1, parms = setup_ui(parms); end;

% automated search for better fit
switch parms.fitflag
  case 1
    parms = fmincon_search(parms);
  case 2
    parms = rand_search(parms);
end;      

if ~parms.plot_flag
  parms = save_trans(parms);
  T = parms.best_T;
  min_error = parms.min_error;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = setup_figure(parms)
  parms.figh = figure('Color',parms.fig_color, 'name', mfilename, ...
     'numbertitle', 'off', 'visible', 'on', 'toolbar','figure');
  parms.ax = axes('tag','UI_plot','parent',parms.figh,...
     'Position',parms.axes_position); 
  axis off;
  parms = calc_error(parms);
  set(gcf, 'userdata', parms);
  % plot head model and points
  redraw;
  view(parms.view);
  drawnow;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = setup_ui(parms)
  % lower left corners of uicontrol groups
  ui_trans_text  = [0.01,   0.76];
  ui_transx      = [0.01,   0.70];
  ui_transy      = [0.01,   0.62];
  ui_transz      = [0.01,   0.54];

  ui_rot_text    = [0.01,   0.41];
  ui_rotx        = [0.01,   0.35];
  ui_roty        = [0.01,   0.27];
  ui_rotz        = [0.01,   0.19];

  ui_fit         = [0.83,   0.80];
  ui_snap        = [0.83,   0.70];
  
  ui_small       = [0.83,   0.50];
  ui_large       = [0.83,   0.40];
  
  ui_revert      = [0.83,   0.20];
  ui_save        = [0.83,   0.10];

  % define commands
  trans_cmd = @translate_points;
  rot_cmd = @rotate_points;
  small_cmd = @set_stepsize;
  large_cmd = @set_stepsize;
  fit_cmd = @fit_points;
  snap_cmd = @snap_points;
  revert_cmd = @revert_trans;
  save_cmd = @overwrite_trans;

  % set positions
  n=1;
  uipos(n,:)  = [ ui_trans_text        0.16     0.06   ]; n=n+1; % translate
  uipos(n,:)  = [ ui_rot_text          0.16     0.06   ]; n=n+1; % rotate

  % create uicontrols
  n=1;
  % trans label, rot label
  u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
   'BackgroundColor',parms.fig_color,'Style','text',...
   'HorizontalAlignment','left',...
   'FontSize',parms.fontsize,...
   'Position',uipos(n,:), ...
   'Tag','UI_trans_label','string','Translation');
  n=n+1;
  u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
   'BackgroundColor',parms.fig_color,'Style','text', ...
   'HorizontalAlignment','left',...
   'FontSize',parms.fontsize,...
   'Position',uipos(n,:), ...
   'Tag','UI_rot_label','string','Rotation');
  n=n+1;

  % trans, rot buttons, text entry
  for i=1:length(parms.xyz_list)
    xyz = parms.xyz_list{i};
    for j=1:length(parms.transrot_list)
      transrot = parms.transrot_list{j};
      ui_offset = eval(sprintf('ui_%s%s',transrot,xyz));
      tag_stem = sprintf('UI_%s%s',transrot,xyz);
      switch transrot
        case 'trans'
          stepsize = parms.trans_stepsize_manual;
          tmp_cmd = trans_cmd;
        case 'rot'
          stepsize = parms.rot_stepsize_manual;
          tmp_cmd = rot_cmd;
      end;

      % set positions
      n=1;
      uipos(n,:)  = [ 0.0        0.00      0.055     0.07   ]; n=n+1; % -
      uipos(n,:)  = [ 0.054      0.001     0.055     0.07   ]; n=n+1; % txt
      uipos(n,:)  = [ 0.111      0.00      0.055     0.07   ]; n=n+1; % +
      uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
      uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

      % create uicontrols
      n=1;
      % +- buttons, text box
      u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
       'FontSize',parms.fontsize,...
       'Position',uipos(n,:), ...
       'Tag',[tag_stem '_neg'],'string',['-' upper(xyz)],...
       'callback',tmp_cmd);
      n=n+1;
      u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
       'BackgroundColor',[1 1 1], ...
       'FontSize',parms.fontsize,...
       'Position',uipos(n,:), ...
       'Style','edit', ...
       'Tag',[tag_stem '_txt'],...
       'string',stepsize);
      n=n+1;
      u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
       'FontSize',parms.fontsize,...
       'Position',uipos(n,:), ...
       'Tag',[tag_stem '_pos'],'string',['+' upper(xyz)],...
       'callback',tmp_cmd);
      n=n+1;

      clear uipos u;
    end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%
  % small button
  ui_offset = ui_small;

  % set position
  n=1;
  uipos(n,:)  = [ ui_small              0.16     0.08   ]; n=n+1;

  % create uicontrol
  n=1;
  % fit button
  u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
   'FontSize',parms.fontsize,...
   'Position',uipos(n,:), ...
   'Tag','UI_small','string','Small Steps',...
   'TooltipString','switch to small step sizes',...
   'callback', small_cmd);
  n=n+1;

  %%%%%%%%%%%%%%%%%%%%%%
  % large button
  ui_offset = ui_large;

  % set position
  n=1;
  uipos(n,:)  = [ ui_large              0.16     0.08   ]; n=n+1;

  % create uicontrol
  n=1;
  % fit button
  u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
   'FontSize',parms.fontsize,...
   'Position',uipos(n,:), ...
   'Tag','UI_large','string','Large Steps',...
   'TooltipString','switch to large step sizes',...
   'callback', large_cmd);
  n=n+1;

  %%%%%%%%%%%%%%%%%%%%%%
  % revert button
  ui_offset = ui_revert;

  % set position
  n=1;
  uipos(n,:)  = [ ui_revert              0.16     0.08   ]; n=n+1;

  % create uicontrol
  n=1;
  % fit button
  u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
   'FontSize',parms.fontsize,...
   'Position',uipos(n,:), ...
   'Tag','UI_revert','string','Revert',...
   'TooltipString','revert to original trans',...
   'callback', revert_cmd);
  n=n+1;

  %%%%%%%%%%%%%%%%%%%%%%
  % save button
  if ~strcmp(parms.intransfile,parms.outtransfile)
    ui_offset = ui_save;

    % set position
    n=1;
    uipos(n,:)  = [ ui_save              0.16     0.08   ]; n=n+1;

    % create uicontrol
    n=1;
    % fit button
    u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
     'FontSize',parms.fontsize,...
     'Position',uipos(n,:), ...
     'Tag','UI_save','string','Save',...
     'TooltipString','overwrite input transfile',...
     'callback', save_cmd);
    n=n+1;
  end;

  %%%%%%%%%%%%%%%%%%%%%%
  % fit points button
  if parms.fitflag>0
    ui_offset = ui_fit;

    % set position
    n=1;
    uipos(n,:)  = [ ui_fit              0.16     0.08   ]; n=n+1;

    % create uicontrol
    n=1;
    % fit button
    u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
     'FontSize',parms.fontsize,...
     'Position',uipos(n,:), ...
     'Tag','UI_fit','string','Fit Points',...
     'TooltipString','fit points with random search',...
     'callback', fit_cmd);
    n=n+1;
  end;

  %%%%%%%%%%%%%%%%%%%%%%
  % snap points button
  if parms.fit_card_flag
    ui_offset = ui_snap;

    % set position
    n=1;
    uipos(n,:)  = [ ui_snap              0.16     0.08   ]; n=n+1;

    % create uicontrol
    n=1;
    % snap button
    u(n) = uicontrol('Parent',parms.figh,'Units', 'normalized',...
     'FontSize',parms.fontsize,...
     'Position',uipos(n,:), ...
     'Tag','UI_snap','string','Snap Points',...
     'TooltipString','snap points using cardinal points',...
     'callback', snap_cmd);
    n=n+1;
  end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = fit_card_points(parms)
  if parms.verbose
    fprintf('%s: reference lpa:    [%s]\n',mfilename,sprintf('%0.3f ',parms.init_lpa));
    fprintf('%s: reference nasion: [%s]\n',mfilename,sprintf('%0.3f ',parms.init_nasion));
    fprintf('%s: reference rpa:    [%s]\n',mfilename,sprintf('%0.3f ',parms.init_rpa));
  end;
  ind_card_points = find(strcmp('cardinal',parms.point_types));
  if length(ind_card_points)~=3
    fprintf('%s: WARNING: wrong number of cardinal points (%d instead of 3)\n',...
      mfilename,length(ind_card_points));
    return;
  end;
  mri_points = [parms.init_lpa;parms.init_nasion;parms.init_rpa];
  head_points = parms.coords(ind_card_points,1:3);
  [d,Z,transform] = procrustes(mri_points,head_points);
  c = mean(mri_points) - mean(head_points) * transform.T; % no scaling
  T = eye(4);
  T(1:3,1:3) = inv(transform.T);
  T(1:3,4) = c;
  parms.T = T;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = save_trans(parms)
  if parms.verbose
    fprintf('%s: current transformation matrix:\n',mfilename);
    disp(parms.T);
    parms = calc_error(parms);
    fprintf('%s: error = %f\n',mfilename,parms.error);
    ind_card_points = find(strcmp('cardinal',parms.point_types));
    Y = parms.coords(ind_card_points,:);
    Z = (parms.T*Y')';
    fprintf('%s: lpa:       [%s]\n',mfilename,sprintf('%0.3f ',Z(1,1:3)));
    fprintf('%s: nasion:    [%s]\n',mfilename,sprintf('%0.3f ',Z(2,1:3)));
    fprintf('%s: rpa:       [%s]\n',mfilename,sprintf('%0.3f ',Z(3,1:3)));
  end;

  T = parms.T;
  if strcmp(parms.trans_type,'mri2head')
    T = inv(T);
  end;
  ts_write_transfile(parms.outtransfile,T);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = rand_search(parms)
  if isfield(parms,'error') && ~isempty(parms.error) && parms.verbose
    fprintf('%s: initial error = %f\n',mfilename,parms.error);
  end;
  num_better_fits = 0;
  for i=1:parms.niters
    parms.T = parms.best_T;
    parms.transvec = parms.best_transvec + parms.trans_stepsize*randn(3,1);
    parms.rotvec = parms.best_rotvec + parms.rot_stepsize*randn(3,1);
    parms = translate_matrix(parms);
    parms = rotate_matrix(parms);
    parms = calc_error(parms);
    if parms.error < parms.min_error
      parms.best_transvec = parms.transvec;
      parms.best_rotvec = parms.rotvec;
      parms.best_T = parms.T;
      parms.min_error = parms.error;
      if parms.verbose
        fprintf('%s: # better fit: trans = [%s], rot = [%s], err = %f, niter %d\n',...
          mfilename,sprintf('%0.3f ',parms.transvec),...
          sprintf('%0.3f ',parms.rotvec),parms.error,i);
      end;
      num_better_fits = 1;
    end;
  end;
  parms.T = parms.best_T;
  if parms.verbose
    fprintf('%s: #### best fit: trans = [%s], rot = [%s], err = %f, %d niters, %d better fits\n',...
      mfilename,sprintf('%0.3f ',parms.best_transvec),...
      sprintf('%0.3f ',parms.best_rotvec),parms.min_error,...
      parms.niters,num_better_fits);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = fmincon_search(parms)
  parms.best_T = ts_fit_points_to_surf(parms.hptsfile,parms.surffile,...
    'trans_type','head2mri',...
    'badpoints',parms.excluded,...
    'T_init',parms.best_T,'verbose',parms.verbose);
  parms.T = parms.best_T;
  parms = calc_error(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_error(parms)
  %% todo: use outer skull too?

  tmp_coords = (parms.T*parms.coords')';
  tmp_coords = tmp_coords(:,1:3);
  errvec = zeros(parms.npoints,1);
  for i=1:parms.npoints
    xyz = squeeze(tmp_coords(i,:));
    dist = 0;
    for j=1:length(xyz)
      dist = dist + (parms.surf.vertices(:,j) - xyz(j)).^2;
    end;
    % find closest surface point
    dist = sqrt(dist);
    [mindist,j] = min(dist);
    errvec(i) = mindist;
  end;

  % ignore points with largest distance (e.g. 95% percentile)
  if parms.cost_include_percentile<100
    nbins = 100;
    [N,X] = hist(errvec,nbins);
    P = cumsum(N); Pnorm = P/P(end);
    [tmp,ind_min] = min(abs(Pnorm-parms.cost_include_percentile/100));
    ind_exclude = find(errvec>X(ind_min));
    errvec(ind_exclude) = 0;
    parms.excluded = ind_exclude;
  else
    parms.excluded = [];
  end;

  parms.error = sum(errvec);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_surf(parms)
  parms.surf = fs_read_trisurf(parms.surffile);
  parms.surf.vertices = parms.surf.vertices/1000; % convert to meters
  % get initial card point coords from surf.vertices
  if parms.init_card_flag
    tmp_coords = parms.surf.vertices;
    meanvec = mean(tmp_coords,1);
    for i=1:3
      tmp_coords(:,i) = tmp_coords(:,i) - meanvec(i);
    end;
    maxvec = max(tmp_coords,[],1);
    minvec = min(tmp_coords,[],1);
    parms.init_lpa = [minvec(1),0,0] + meanvec;
    parms.init_nasion = [0,maxvec(2),0] + meanvec;
    parms.init_rpa = [maxvec(1),0,0] + meanvec;
  end;
  % read inner skull file
  if ~isempty(parms.iskullfile)
    parms.iskull = fs_read_trisurf(parms.iskullfile);
    parms.iskull.vertices = parms.iskull.vertices/1000; % convert to meters
  else
    parms.iskull = [];
  end;
  % read outer skull file
  if ~isempty(parms.oskullfile)
    parms.oskull = fs_read_trisurf(parms.oskullfile);
    parms.oskull.vertices = parms.oskull.vertices/1000; % convert to meters
  else
    parms.oskull = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function parms = load_hpts(parms)
  % load hpts file
  points = [];
  fid=fopen(parms.hptsfile,'rt');
  if fid==-1
    error('failed to open file %s',parms.hptsfile);
  end;
  p = 1;
  while (~feof(fid))
    temp=fgetl(fid);
    if strcmp(temp(1),'#'), continue; end;
    spaces=regexp(temp,' +');
    strings = [];
    k = 1;
    for j=1:length(spaces)
      strings{j} = temp(k:spaces(j)-1);
      k = spaces(j) + 1;
    end;
    strings{j+1} = temp(k:length(temp));
    points(p).name =  strings{1};
    points(p).num = str2num(strings{2});
    points(p).x = str2double(strings{3});
    points(p).y = str2double(strings{4});
    points(p).z = str2double(strings{5});
    p=p+1;
  end
  fclose(fid);
  % reformat point coords
  parms.npoints = length(points);
  parms.coords = [];
  parms.coords(:,1) = cell2mat({points.x});
  parms.coords(:,2) = cell2mat({points.y});
  parms.coords(:,3) = cell2mat({points.z});
  parms.coords(:,4) = ones(parms.npoints,1);
  parms.point_types = {points.name};
  parms.point_nums = cell2mat({points.num});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = init_trans(parms)
  % load trans file
  if parms.verbose
    fprintf('%s: initializing transformation matrix...\n',mfilename);
  end;
  intransfile = parms.intransfile;
  if isempty(intransfile)
    parms.intransfile = parms.outtransfile;
  elseif ~exist(parms.intransfile,'file')
    if parms.verbose
      fprintf('%s: WARNING: input transfile not found\n',...
        mfilename);
    end;
    intransfile = [];
  end;
  if ~isempty(intransfile)
    if parms.verbose
      fprintf('%s: reading transfile %s...\n',mfilename,intransfile);
    end;
    if parms.oldtrans_flag
      parms.T = ts_read_transfile(intransfile,0.001);
    else
      parms.T = ts_read_transfile(intransfile);
    end;
  elseif ~isempty(parms.T_init)
    if parms.verbose
      fprintf('%s: using T_init...\n',mfilename);
    end;
    parms.T = parms.T_init;
  elseif parms.fit_card_flag
    if parms.verbose
      fprintf('%s: fitting for init trans using cardinal points...\n',mfilename);
    end;
    % fit cardinal points to init_lpa, init_rpa, init_nasion
    parms = fit_card_points(parms);
    if strcmp(parms.trans_type,'mri2head')
      parms.T = inv(parms.T);
    end;
  else
    if parms.verbose
      fprintf('%s: using identity matrix...\n',mfilename);
    end;
    parms.T = eye(4);
  end;
  % internal trans_type is always head2mri (applied to head points)
  if strcmp(parms.trans_type,'mri2head')
    parms.T = inv(parms.T);
  end;
  if parms.verbose
    fprintf('%s: initial transformation matrix:\n',mfilename);
    disp(parms.T);
  end;
  parms.T_init = parms.T;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = translate_matrix(parms)
  T = eye(4);
  for i=1:3
    T(i,4) = parms.transvec(i)/1000; % convert to meters
  end;
  parms.T = T*parms.T;
return;

function parms = rotate_matrix(parms)
  rot = parms.rotvec*pi/180; % expect degrees, convert to radians
  Rx=eye(3);
  Ry=eye(3);
  Rz=eye(3);
  if(rot(1))
    c = cos(rot(1));
    s = sin(rot(1));
    Rx(2,2)=c;
    Rx(2,3)=s;
    Rx(3,2)=-s;
    Rx(3,3)=c;
  end
  if(rot(2))
    c = cos(rot(2));
    s = sin(rot(2));
    Ry(1,1)=c;
    Ry(1,3)=-s;
    Ry(3,1)=s;
    Ry(3,3)=c;
  end
  if(rot(3))
    c = cos(rot(3));
    s = sin(rot(3));
    Rz(1,1)=c;
    Rz(1,2)=s;
    Rz(2,1)=-s;
    Rz(2,2)=c;
  end
  R = Rz*Ry*Rx;
  T = eye(4);
  T(1:3,1:3) = R;
  parms.T = T*parms.T;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function redraw(parms)
  if ~exist('parms','var') || isempty(parms)
    parms = get(gcf, 'userdata');
  end;
  cla;
  hold on;
  % plot surf(s)

  % plot inner skull
  if ~isempty(parms.iskull)
    val_iskull = parms.iskull_color*ones(parms.iskull.nverts,1);
    h_iskull = triplot(parms.iskull.vertices,parms.iskull.faces,val_iskull);
    set(h_iskull,'FaceAlpha',parms.iskull_alpha);
  end;

  % plot outer skull
  if ~isempty(parms.oskull)
    val_oskull = parms.oskull_color*ones(parms.oskull.nverts,1);
    h_oskull = triplot(parms.oskull.vertices,parms.oskull.faces,val_oskull);
    set(h_oskull,'FaceAlpha',parms.oskull_alpha);
  end;

  % plot scalp
  val_scalp  = parms.scalp_color*ones(parms.surf.nverts,1);
  h_scalp = triplot(parms.surf.vertices,parms.surf.faces,val_scalp);
  set(h_scalp,'FaceAlpha',parms.scalp_alpha);

  colormap(parms.colormap);
  caxis([0 1]);

  % transform coordinates
  tmp_coords = (parms.T*parms.coords')';
  tmp_coords = tmp_coords(:,1:3);

  % plot cardinal points
  ind = find(strcmp(parms.point_types,'cardinal'));
  if ~isempty(ind)
    plot3(tmp_coords(ind,1),tmp_coords(ind,2),tmp_coords(ind,3),...
      [parms.card_point_color '.'],'MarkerSize',parms.point_size);
  end;
  
  % plot extra points
  ind = find(strcmp(parms.point_types,'extra'));
  if ~isempty(ind)
    included = ind;
    if ~isfield(parms,'excluded'), parms.excluded = []; end;
    if ~isempty(parms.excluded)
      included = setdiff(included,parms.excluded);
    end;
    plot3(tmp_coords(included,1),tmp_coords(included,2),tmp_coords(included,3),...
      [parms.point_color '.'],'MarkerSize',parms.point_size);
    if ~isempty(parms.excluded)
      plot3(tmp_coords(parms.excluded,1),...
            tmp_coords(parms.excluded,2),tmp_coords(parms.excluded,3),...
        [parms.excl_point_color '.'],'MarkerSize',parms.point_size);
    end;
  end;

  % plot eeg points
  ind = find(strcmp(parms.point_types,'eeg'));
  if ~isempty(ind)
    plot3(tmp_coords(ind,1),tmp_coords(ind,2),tmp_coords(ind,3),...
      [parms.eeg_point_color '.'],'MarkerSize',parms.point_size);
  end;

  % plot eeg_ref point
  ind = find(strcmp(parms.point_types,'eeg_ref'));
  if ~isempty(ind)
    plot3(tmp_coords(ind,1),tmp_coords(ind,2),tmp_coords(ind,3),...
      [parms.eeg_ref_color '.'],'MarkerSize',parms.point_size);
  end;
  
  if parms.point_num_flag
    numstrs = [];
    for i=1:parms.npoints
      numstrs{i} = num2str(parms.point_nums(i));
    end;
    text(tmp_coords(:,1),tmp_coords(:,2),tmp_coords(:,3),numstrs);
  end;


  axis image;
  axis off;
  xlabel('x');
  ylabel('y');
  zlabel('z');

  % create lights
  if parms.lights_flag
    h1 = light('Color',[0.5,0.5,0.5]);
    lightangle(h1,0,90);
    
    h2 = light('Color',[0.5,0.5,0.5]);
    lightangle(h2,70,20);

    h3 = light('Color',[0.5,0.5,0.5]);
    lightangle(h3,270,20);
  end;  

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_stepsize(hobj,ed)
  parms = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  if strcmp(tmp_UI,'UI_small')
    small_flag = 1;
  else
    small_flag = 0;
  end;
  for j=1:length(parms.transrot_list)
    transrot = parms.transrot_list{j};
    switch transrot
      case 'trans'
        if small_flag
          stepsize = parms.trans_stepsize;
        else
          stepsize = parms.trans_stepsize_manual;
        end;
      case 'rot'
        if small_flag
          stepsize = parms.rot_stepsize;
        else
          stepsize = parms.rot_stepsize_manual;
        end;
    end;
    for i=1:length(parms.xyz_list)
      xyz = parms.xyz_list{i};
      tag = sprintf('UI_%s%s_txt',transrot,xyz);
      set(findobj(gcbf,'tag',tag),'string',num2str(stepsize));
    end;
  end;
return;

function revert_trans(hobj,ed)
  parms = get(gcbf, 'userdata');
  parms.T = parms.T_init;
  parms = save_trans(parms);
  set(gcbf,'userdata',parms);
  redraw;
return;

function overwrite_trans(hobj,ed)
  parms = get(gcbf, 'userdata');
  intransfile = parms.intransfile;
  outtransfile = parms.outtransfile;
  if exist(intransfile,'file')
    question = sprintf('Are you sure you want to overwrite %s?',...
      intransfile);
    btn = questdlg(question,'Overwrite?','OK','Cancel','Cancel');
    if isempty(btn) || strcmp(btn,'Cancel'), return; end;
  end;
  parms.outtransfile = parms.intransfile;
  parms = save_trans(parms);
  if exist(outtransfile,'file')
    delete(outtransfile);
  end;
return;

function translate_points(hobj,ed)
  transx = 0;
  transy = 0;
  transz = 0;
  parms = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  if strcmp(tmp_UI,'UI_transx_neg')
    transx = -str2num(get(findobj(gcbf,'tag','UI_transx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transx_pos')
    transx = str2num(get(findobj(gcbf,'tag','UI_transx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transy_neg')
    transy = -str2num(get(findobj(gcbf,'tag','UI_transy_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transy_pos')
    transy = str2num(get(findobj(gcbf,'tag','UI_transy_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transz_neg')
    transz = -str2num(get(findobj(gcbf,'tag','UI_transz_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transz_pos')
    transz = str2num(get(findobj(gcbf,'tag','UI_transz_txt'),'string'));
  end;
  parms.transvec = [transx,transy,transz];
  parms = translate_matrix(parms);
  parms = save_trans(parms);
  set(gcbf,'userdata',parms);
  redraw;
return;


function rotate_points(hobj,ed)
  rotx = 0;
  roty = 0;
  rotz = 0;
  parms = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  if strcmp(tmp_UI,'UI_rotx_neg')
    rotx = -str2num(get(findobj(gcbf,'tag','UI_rotx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_rotx_pos')
    rotx = str2num(get(findobj(gcbf,'tag','UI_rotx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_roty_neg')
    roty = -str2num(get(findobj(gcbf,'tag','UI_roty_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_roty_pos')
    roty = str2num(get(findobj(gcbf,'tag','UI_roty_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_rotz_neg')
    rotz = -str2num(get(findobj(gcbf,'tag','UI_rotz_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_rotz_pos')
    rotz = str2num(get(findobj(gcbf,'tag','UI_rotz_txt'),'string'));
  end;
  parms.rotvec = [rotx,roty,rotz];
  parms = rotate_matrix(parms);
  parms = save_trans(parms);
  set(gcbf,'userdata',parms);
  redraw;
return;

function fit_points(hobj,ed)
  parms = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  if parms.verbose
    fprintf('%s: initial transformation matrix:\n',mfilename);
    disp(parms.T);
  end;
  parms.best_T = parms.T;
  parms = calc_error(parms);
  parms.min_error = parms.error;
  switch parms.fitflag
    case 1
      parms = fmincon_search(parms);
    case 2
      parms = rand_search(parms);
  end;
  parms = save_trans(parms);
  set(gcbf,'userdata',parms);
  redraw;
return;

function snap_points(hobj,ed)
  parms = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  parms = fit_card_points(parms);
  parms = save_trans(parms);
  set(gcbf,'userdata',parms);
  redraw;
return;

