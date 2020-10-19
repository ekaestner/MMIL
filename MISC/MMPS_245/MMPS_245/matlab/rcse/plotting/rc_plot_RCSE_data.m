function rc_plot_RCSE_data(surffile,prefix,plane,start_times,end_times)
%function rc_plot_RCSE_data(surffile,prefix,plane,start_times,end_times)
%
% Purpose: plot RCSE potentials, fields, and dipoles
%   for circularly arranged stimuli
%
% Required Input:
%   surffile: file name of freesurfer compatible tri file
%     e.g. outer scalp surface
%        used for EEG contour plots
%   prefix: RCSE prefix
%
% Optional Input:
%   plane: 'cor', 'sag', or 'hor'
%     {default: 'cor'}
%   start_times: vector of start times (msec)
%     {default: [80]}
%   end_times: vector of end times (msec)
%     {default: [90]}
%
% Creeated: 03/04/08  by Don Hagler
% Last Mod: 08/01/13  by Don Hagler
%

%% todo: use varargin
%% todo: is this function necessary (see rc_plot_RCSE_datafit)

if ~mmil_check_nargs(nargin,2), return; end;
titleflag = 0;
if ~exist('plane','var') || isempty(plane), plane = 'cor'; end;
if ~ismember(plane,{'cor','sag','hor'})
  error('plane must be cor, sag, or hor');
end;
if ~exist('start_times','var') || isempty(start_time), start_times = [80]; end;
if ~exist('end_times','var') || isempty(end_time), end_times = [90]; end;
if (strcmp(plane,'cor'))
  view_angle = [0 0];
elseif (strcmp(plane,'sag'))
  view_angle = [90 0];
elseif (strcmp(plane,'hor'))
  view_angle = [0 90];
end;

if ~exist('fontname','var') | isempty(fontname), fontname = 'Arial'; end;

outfix = sprintf('contplots-%s',plane);

for j=1:length(start_times)
  time0 = start_times(j);
  time1 = end_times(j);
  fprintf('%s: plotting data for times %0.1f - %0.1f\n',...
    mfilename,time0,time1);

  clf;
  fprintf('%s: plotting EEG potentials...\n',mfilename);
  plot_retinv_pots(surffile,prefix,time0,time1,view_angle);
  fprintf('%s: plotting MEG fields...\n',mfilename);
  plot_retinv_fields(prefix,time0,time1,view_angle);
  fprintf('%s: plotting dipole orientations...\n',mfilename);
  rc_plot_dips_circ(prefix,[],plane);

  if titleflag
    ax=axes('Units','Normal','Position',[0.075 0.075 0.89 0.89],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    set(gca,'FontName',fontname)
    title([sprintf('MEG/EEG Retinotopy Contour Plots t=%0.1f-%0.1f',time0,time1)]);
    set(gcf,'position',[50 50 1000 1000]);
  end;

  colormap(mmil_cmap_blueblackred);

  outstem = sprintf('%s_t%03d-%03d',outfix,time0,time1)
  print('-dtiff',[outstem '.tif']);
  mmil_printeps(gcf,[outstem '.eps']);
end;

