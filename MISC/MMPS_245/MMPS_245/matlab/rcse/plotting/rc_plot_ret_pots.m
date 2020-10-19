function rc_plot_ret_pots(avg_data,varargin)
%function rc_plot_ret_pots(avg_data,[options])
%
% Purpose: plot EEG potentials on head surface using
%   modified version of EEGLAB's headplot
%   plots multiple display_locs in a ring
%
% Required input:
%  avg_data - average data structure (see avg_fif_data)
%
% Optional parameters:
%  transfile - file name for 4x4 head2mri transformation matrix
%    {default: use identity matrix}
%  trans - 4x4 head2mri transformation matrix -- 
%    overrules transfile if both are supplied
%    {default: use identity matrix}
%  locfile - file name for saving eeg locations
%    {default: 'matfiles/eeglocs.sfp'}
%  splinefile - file name for saving spherical spline
%    {default: 'matfiles/spline.mat'}
%  meshfile - file name for saving mesh
%    {default: 'matfiles/mesh.mat'}
%  surffile - name of surface file (freesurfer tri file format)
%    {default: 'outer_scalp.tri'}
%  badchanfile - name of text file containing bad channel labels
%    {default: []}
%  time0 - start time of averaged data range (msec)
%    {default: 60}
%  time1 - end time of averaged data range (msec)
%    {default: 70}
%  view_angle - 2x1 vector containing azimuth and elevation angles
%    {default: [0 20]}
%  electrodes - ['on'|'off'] toggle whether electrode positions are
%    displayed on head (useful for checking registration)
%    {default: 'off'}
%  scale_max - max value for color scale (uVolts)
%    {default: 4}
%  plotsize - size of single head plot (relative to figure)
%    {default: 0.12}
%  radius - radius of circle
%    {default: 0.3}
%  angle_offset - angle offset of first plot (relative to x=radius,y=0)
%    {default: 11.25}
%  actual_num_locs - number of stimulus locations
%    {default: [16]}
%  display_locs - vector of condition numbers (not event codes) to display
%    {default: [1:16]}
%
% example of trans file format:
%1.000000 0.000000 0.000000 0.000000 
%0.000000 1.000000 0.000000 0.000000 
%0.000000 0.000000 1.000000 0.000000 
%0.000000 0.000000 0.000000 1.000000 
%
% Created:  05/17/06 by Don Hagler
% Last Mod: 06/16/14 by Don Hagler
%

%% todo: use varargin

if ~mmil_check_nargs(nargin,1), return; end;

DEFAULT_LOCFILE = 'matfiles/eeglocs.sfp';
DEFAULT_SPLINEFILE = 'matfiles/spline.mat';
DEFAULT_MESHFILE = 'matfiles/mesh.mat';
DEFAULT_SURFFILE = 'outer_scalp.tri';
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 60;
DEFAULT_TIME1 = 70;
DEFAULT_VIEW_ANGLE = [0 20];
DEFAULT_ELECTRODES = 'off';
DEFAULT_scale_max = 4;
DEFAULT_PLOTSIZE = 0.12;
DEFAULT_PLOTSIZE_INIT = 0.01;
DEFAULT_RADIUS = 0.3;
DEFAULT_RADIUS_INIT = 0.1;
DEFAULT_ANGLE_OFFSET = 11.25;
DEFAULT_ACTUAL_NUM_LOCS = 16;
DEFAULT_DISPLAY_LOCS = [1:DEFAULT_ACTUAL_NUM_LOCS];

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), opt=struct(options{:}); 
  else opt = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

try, opt.transfile;       catch, opt.transfile = []; end;
try, opt.trans;           catch, opt.trans = []; end;
try, opt.locfile;         catch, opt.locfile = DEFAULT_LOCFILE; end;
try, opt.splinefile;      catch, opt.splinefile = DEFAULT_SPLINEFILE; end;
try, opt.meshfile;        catch, opt.meshfile = DEFAULT_MESHFILE; end;
try, opt.surffile;        catch, opt.surffile = DEFAULT_SURFFILE; end;
try, opt.badchanfile;     catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;           catch, opt.time0 = DEFAULT_TIME0; end;
try, opt.time1;           catch, opt.time1 = DEFAULT_TIME1; end;
try, opt.view_angle;      catch, opt.view_angle = DEFAULT_VIEW_ANGLE; end;
try, opt.electrodes;      catch, opt.electrodes = DEFAULT_ELECTRODES; end;
try, opt.scale_max;       catch, opt.scale_max = DEFAULT_scale_max; end;
try, opt.plotsize;        catch, opt.plotsize = DEFAULT_PLOTSIZE; end;
try, opt.radius;          catch, opt.radius = DEFAULT_RADIUS; end;
try, opt.angle_offset;    catch, opt.angle_offset = DEFAULT_ANGLE_OFFSET; end;
try, opt.display_locs;    catch, opt.display_locs = DEFAULT_DISPLAY_LOCS; end;
try, opt.actual_num_locs; catch, opt.actual_num_locs = DEFAULT_ACTUAL_NUM_LOCS; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'transfile' 'trans' 'locfile' 'splinefile' 'meshfile'...
         'surffile' 'badchanfile' 'time0' 'time1' 'view_angle'...
         'electrodes' 'scale_max' 'plotsize' 'radius' 'angle_offset'...
         'display_locs' 'actual_num_locs'...
   },;
   otherwise, error(['rc_plot_ret_pots: unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
transfile = opt.transfile;
trans = opt.trans;
locfile = opt.locfile;
splinefile = opt.splinefile;
meshfile = opt.meshfile;
surffile = opt.surffile;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
view_angle = opt.view_angle;
electrodes = opt.electrodes;
scale_max = opt.scale_max;
plotsize = opt.plotsize;
radius = opt.radius;
angle_offset = opt.angle_offset;
display_locs = opt.display_locs;
actual_num_locs = opt.actual_num_locs;
clear opt optfields options;

% set start and end samples
sfreq = avg_data.sfreq;
t_trigger = avg_data.averages(1).time(1)*1000;
t0 = round((time0 - t_trigger)*sfreq/1000);
t1 = round((time1 - t_trigger)*sfreq/1000);

radius_init = DEFAULT_RADIUS_INIT;
plotsize_init = DEFAULT_PLOTSIZE_INIT;
center=[0 0 0];

color_range = [-scale_max*10^-6,scale_max*10^-6];

% get badchans from avg_data
bad_channels = find(cell2mat({avg_data.sensor_info.badchan})==1);
% read badchan file
labels = {avg_data.sensor_info.label}; 
if ~isempty(badchanfile)
  badchan_i = ts_read_txt_badchans(badchanfile,labels);
else
  badchan_i = [];
end;
bad_channels = unique([bad_channels,badchan_i]);
chans = find(strcmp('eeg',{avg_data.sensor_info.typestring}));
chans = setdiff(chans,badchan_i);
coords = zeros(length(chans),3);

% head2mri coordinate transformation
if ~isempty(trans)
  T_head2mri = trans;
elseif ~isempty(transfile)
  T_head2mri = read_transfile(transfile);
else
  T_head2mri = eye(4);
end;
T = T_head2mri;

if ~file_exists(splinefile)
  %% create electrode location file
  fid = fopen(locfile,'w');
  for c=1:length(chans)
    k=chans(c);
    loc = avg_data.sensor_info(k).loc;
    % apply coordinate transform
    loc = T*loc;
    x(c) = loc(1,4);
    y(c) = loc(2,4);
    z(c) = loc(3,4);
    label = sprintf('EE%02d',c);
    fprintf(fid,'\t%s\t%f\t%f\t%f\n',label,x(c),y(c),z(c));
  end;
  fclose(fid);
  % load tri file
  surf=fs_read_trisurf(surffile);
  POS = surf.vertices;
  TRI1 = surf.faces;
  save(meshfile,'POS','TRI1','center');
  ts_headplot('setup',locfile,splinefile,'meshfile',meshfile);
end;

hold on;
for j=1:length(display_locs)
  k = display_locs(j);
  A = angle_offset + (k-1)*360/actual_num_locs;
  plotloc(1) = 0.5 + radius_init*cos(A*pi/180) - plotsize_init/2;
  plotloc(2) = 0.5 + radius_init*sin(A*pi/180) - plotsize_init/2;
  plotloc(3) = plotsize_init;
  plotloc(4) = plotsize_init;

  B = mean(avg_data.averages(k).data(chans,t0:t1),2);

  subplot('position',plotloc);
  pp_axis(j) = ts_headplot(B,splinefile,'meshfile',meshfile,'view',view_angle,...
    'maplimits',color_range,'electrodes',electrodes);
  axis off;
end;

colormap(mmil_cmap_blueblackred);


for j=1:length(display_locs)
  k = display_locs(j);
  A = angle_offset + (k-1)*360/actual_num_locs;
  plotloc(1) = 0.5 + radius*cos(A*pi/180) - plotsize/2;
  plotloc(2) = 0.5 + radius*sin(A*pi/180) - plotsize/2;
  plotloc(3) = plotsize;
  plotloc(4) = plotsize;
  set(pp_axis(j),'position',plotloc);
end;

return;

