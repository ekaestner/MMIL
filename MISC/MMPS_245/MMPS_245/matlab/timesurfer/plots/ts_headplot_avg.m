function ts_headplot_avg(avg_data,varargin);
% ts_headplot_avg - plots eeg potentials from avg_data structure on head surface
%   using modified version of EEGLAB's headplot
%
% Usage:
%   ts_headplot_avg(avg_data,'key1', value1,...);
%
% equired input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  condition - condition number (not event code) used to index avg_data.averages
%    {default: 1}
%  transfile - file name for 4x4 transformation matrix
%    {default: use identity matrix}
%  trans - 4x4 transformation matrix -- overrules transfile if both are supplied
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
%    {default: 50}
%  time1 - end time of averaged data range (msec)
%    {default: 60}
%  view_angle - 2x1 vector containing azimuth and elevation angles
%    {default: [0 20]}
%  electrodes - ['on'|'off'] toggle whether electrode positions are
%    displayed on head (useful for checking registration)
%    {default: 'off'}
%  scale_max - max value for color scale (uVolts)
%    {default: 4}
%
% example of trans file format:
%1.000000 0.000000 0.000000 0.000000 
%0.000000 1.000000 0.000000 0.000000 
%0.000000 0.000000 1.000000 0.000000 
%0.000000 0.000000 0.000000 1.000000 
%
% Created:  06/05/06 by Don Hagler
% Last Mod: 06/16/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

DEFAULT_LOCFILE = 'matfiles/eeglocs.sfp';
DEFAULT_SPLINEFILE = 'matfiles/spline.mat';
DEFAULT_MESHFILE = 'matfiles/mesh.mat';
DEFAULT_SURFFILE = 'outer_scalp.tri';
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 50;
DEFAULT_TIME1 = 60;
DEFAULT_VIEW_ANGLE = [0 20];
DEFAULT_ELECTRODES = 'off';
DEFAULT_SCALE_MAX = 4;
DEFAULT_CONDITION = 1;

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

try, opt.transfile;     catch, opt.transfile = []; end;
try, opt.trans;         catch, opt.trans = []; end;
try, opt.locfile;       catch, opt.locfile = DEFAULT_LOCFILE; end;
try, opt.splinefile;    catch, opt.splinefile = DEFAULT_SPLINEFILE; end;
try, opt.meshfile;      catch, opt.meshfile = DEFAULT_MESHFILE; end;
try, opt.surffile;      catch, opt.surffile = DEFAULT_SURFFILE; end;
try, opt.badchanfile;   catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;         catch, opt.time0 = DEFAULT_TIME0; end;
try, opt.time1;         catch, opt.time1 = DEFAULT_TIME1; end;
try, opt.view_angle;    catch, opt.view_angle = DEFAULT_VIEW_ANGLE; end;
try, opt.electrodes;    catch, opt.electrodes = DEFAULT_ELECTRODES; end;
try, opt.scale_max;     catch, opt.scale_max = DEFAULT_SCALE_MAX; end;
try, opt.condition;     catch, opt.condition = DEFAULT_CONDITION; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'transfile' 'trans' 'locfile' 'splinefile' 'meshfile'...
         'surffile' 'badchanfile' 'time0' 'time1' 'view_angle'...
         'electrodes' 'scale_max' 'condition'...
   },;
   otherwise, error([mfilename ': unrecognized option: ''' optfields{index} '''' ]);
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
condition = opt.condition;
clear opt optfields options;

% set start and end samples
sfreq = avg_data.sfreq;
t_trigger = avg_data.averages(1).time(1)*1000;
t0 = round((time0 - t_trigger)*sfreq/1000);
t1 = round((time1 - t_trigger)*sfreq/1000);

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
chans = find(strcmp('eeg',lower({avg_data.sensor_info.typestring})));
chans = setdiff(chans,badchan_i);

% head2mri coordinate transformation
if ~isempty(trans)
  T_head2mri = trans;
elseif ~isempty(transfile)
  T_head2mri = ts_read_transfile(transfile);
else
  T_head2mri = eye(4);
end;

T = T_head2mri;

if ~exist(splinefile,'file')
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
  ts_headplot('setup',locfile,splinefile,'meshfile',meshfile,'orilocs','on');
end;

B = mean(avg_data.averages(condition).data(chans,t0:t1),2);
ts_headplot(B,splinefile,'meshfile',meshfile,'view',view_angle,...
  'maplimits',color_range,'electrodes',electrodes);
axis off;

colormap(mmil_cmap_blueblackred);

return;

