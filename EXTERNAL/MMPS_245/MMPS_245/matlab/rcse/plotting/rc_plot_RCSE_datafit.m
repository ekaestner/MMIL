function rc_plot_RCSE_datafit(surffile,prefix,plane,start_times,end_times,plot_flags);
%function rc_plot_RCSE_datafit(surffile,[prefix],[plane],[start_times],[end_times],[plot_flags]);
%
% Purpose: plot potentials, fields, and dipoles for RCSE data and fit
%
% Required Input:
%   surffile: freesurfer compatible tri file  e.g. outer scalp surface
%   prefix: RCSE prefix
%
% Optional Input:
%   plane: 'cor', 'sag', or 'hor'
%     {default = 'cor'}
%   start_times: vector of start times (msec)
%     {default = [80]}
%   end_times: vector of end times (msec)
%     {default = [90]}
%   plot_flags: vector of 4 numbers, 1 or 0, indicating whether to plot the following:
%     data, fit, error, fit areas
%     {default = [1 1 1 1]}
%
% Early Mod: 03/04/08 by Don Hagler
% Last Mod:  08/01/13 by Don Hagler
%

%% todo: use varargin
%% todo: use loop instead of three sections of redundant code

if ~mmil_check_nargs(nargin,2), return; end;
titleflag = 0;
if ~exist('plane','var') || isempty(plane), plane = 'cor'; end;
if ~ismember(plane,{'cor','sag','hor'})
  error('plane must be cor, sag, or hor');
end;

if ~exist('start_times','var') || isempty(start_times), start_times = [80]; end;
if ~exist('end_times','var') || isempty(end_times), end_times = [90]; end;

if (strcmp(plane,'cor'))
  view_angle = [0 0];
elseif (strcmp(plane,'sag'))
  view_angle = [90 0];
elseif (strcmp(plane,'hor'))
  view_angle = [0 90];
end;

if ~exist('plot_flags','var') || isempty(plot_flags), plot_flags = ones(4,1); end;
if length(plot_flags)~=4
  fprintf('%s: plot_flags should have 4 values\n',mfilename);
  return;
end;

if ~exist('fontname','var') | isempty(fontname), fontname = 'Arial'; end;

plot_data_flag = plot_flags(1);
plot_fit_flag = plot_flags(2);
plot_err_flag = plot_flags(3);
plot_areas_flag = plot_flags(4);

imgsize = 500;
vvfpfile = 'vectorview_field_plot.mat';
field_scale_max = 50;
field_radius = 0.43;
field_plotsize = 0.13;
pot_scale_max = 2.5;
pot_radius = 0.3;
pot_plotsize = 0.12;
electrodes = 'off';

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_forward.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

if parms.usemag_flag
  usemags = 'on';
else
  usemags = 'off';
end;

badchanfile = parms.badchanfile;
T_mri2head = parms.trans;
T_head2mri = inv(T_mri2head);

num_locs = retmap.num_locs;
display_locs = retmap.cond_order;
orig_num_locs = retmap.orig_num_locs;
angle_offset = 360/(orig_num_locs*2);

for j=1:length(start_times)
  time0 = start_times(j);
  time1 = end_times(j);

  fprintf('%s: plotting fields for times %0.1f - %0.1f\n',...
    mfilename,time0,time1);

  if plot_data_flag
    clf;
    fprintf('%s: plotting measured MEG fields...\n',mfilename);
    rc_plot_ret_fields(avg_data,'vvfpfile',vvfpfile,'badchanfile',badchanfile,...
      'scale_max',field_scale_max,...
      'time0',time0','time1',time1,'view_angle',view_angle,...
      'radius',field_radius,'angle_offset',angle_offset,...
      'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
      'plotsize',field_plotsize,'usemags',usemags);
    fprintf('%s: plotting measured EEG potentials...\n',mfilename);
    rc_plot_ret_pots(avg_data,'surffile',surffile,'badchanfile',badchanfile,...
      'trans',T_head2mri,'scale_max',pot_scale_max,...
      'time0',time0','time1',time1,'view_angle',view_angle,...
      'electrodes',electrodes,'radius',pot_radius,'angle_offset',angle_offset,...
      'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
      'plotsize',pot_plotsize);
    fprintf('%s: plotting dipole orientations...\n',mfilename);
    rc_plot_dips_circ(prefix,[],plane);
    set(gca,'FontName',fontname)
    if titleflag
      ax=axes('Units','Normal','Position',[0.075 0.075 0.89 0.89],'Visible','off');
      set(get(ax,'Title'),'Visible','on');
      title([sprintf('Measured MEG Retinotopy Contour Plots t=%0.1f-%0.1f',time0,time1)]);
      set(gcf,'position',[50 50 imgsize imgsize]);
    end;
    outstem = sprintf('contplots-%s-data-%s_t%03d-%03d',...
      prefix,plane,time0,time1)
    print('-dtiff',[outstem '.tif']);
    mmil_printeps(gcf,[outstem '.eps']);
  end;

  if plot_fit_flag
    clf;
    fprintf('%s: plotting fitted MEG fields...\n',mfilename);
    datamatfile = ...
      sprintf('matfiles/%s_fiterr.mat',prefix);
    load(datamatfile);
    rc_plot_ret_fields(fit_data,'vvfpfile',vvfpfile,'badchanfile',badchanfile,...
      'scale_max',field_scale_max,...
      'time0',time0','time1',time1,'view_angle',view_angle,...
      'radius',field_radius,'angle_offset',angle_offset,...
      'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
      'plotsize',field_plotsize,'usemags',usemags);
    fprintf('%s: plotting fitted EEG potentials...\n',mfilename);
    rc_plot_ret_pots(fit_data,'surffile',surffile,'badchanfile',badchanfile,...
      'trans',T_head2mri,'scale_max',pot_scale_max,...
      'time0',time0','time1',time1,'view_angle',view_angle,...
      'electrodes',electrodes,'radius',pot_radius,'angle_offset',angle_offset,...
      'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
      'plotsize',pot_plotsize);
    fprintf('%s: plotting dipole orientations...\n',mfilename);
    rc_plot_dips_circ(prefix,[],plane);
    set(gca,'FontName',fontname)
    if titleflag
      ax=axes('Units','Normal','Position',[0.075 0.075 0.89 0.89],'Visible','off');
      set(get(ax,'Title'),'Visible','on');
      title([sprintf('Expected MEG Retinotopy Contour Plots t=%0.1f-%0.1f',time0,time1)]);
      set(gcf,'position',[50 50 imgsize imgsize]);
    end;
    outstem = sprintf('contplots-%s-fit-%s_t%03d-%03d',...
      prefix,plane,time0,time1)
    print('-dtiff',[outstem '.tif']);
    mmil_printeps(gcf,[outstem '.eps']);
  end;

  if plot_err_flag
    clf;
    fprintf('%s: plotting residual error MEG fields...\n',mfilename);
    datamatfile = ...
      sprintf('matfiles/%s_fiterr.mat',prefix);
    load(datamatfile);
    rc_plot_ret_fields(err,'vvfpfile',vvfpfile,'badchanfile',badchanfile,...
      'scale_max',field_scale_max,...
      'time0',time0','time1',time1,'view_angle',view_angle,...
      'radius',field_radius,'angle_offset',angle_offset,...
      'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
      'plotsize',field_plotsize,'usemags',usemags);
    fprintf('%s: plotting residual error EEG potentials...\n',mfilename);
    rc_plot_ret_pots(err,'surffile',surffile,'badchanfile',badchanfile,...
      'trans',T_head2mri,'scale_max',pot_scale_max,...
      'time0',time0','time1',time1,'view_angle',view_angle,...
      'electrodes',electrodes,'radius',pot_radius,'angle_offset',angle_offset,...
      'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
      'plotsize',pot_plotsize);
    fprintf('%s: plotting dipole orientations...\n',mfilename);
    rc_plot_dips_circ(prefix,[],plane);
    set(gca,'FontName',fontname)
    if titleflag
      ax=axes('Units','Normal','Position',[0.075 0.075 0.89 0.89],'Visible','off');
      set(get(ax,'Title'),'Visible','on');
      title([sprintf('Residual Error MEG Retinotopy Contour Plots t=%0.1f-%0.1f',time0,time1)]);
      set(gcf,'position',[50 50 imgsize imgsize]);
    end;
    outstem = sprintf('contplots-%s-err-%s_t%03d-%03d',...
      prefix,plane,time0,time1)
    print('-dtiff',[outstem '.tif']);
    mmil_printeps(gcf,[outstem '.eps']);
  end;  

  if plot_areas_flag
    for a=1:retmap.num_areas
      datamatfile = ...
        sprintf('matfiles/%s_fit_%s.mat',prefix,retmap.areas(a).name);
      load(datamatfile);
      clf;
      fprintf('%s: plotting fitted MEG fields for area %s...\n',...
        mfilename,retmap.areas(a).name);
      rc_plot_ret_fields(fit_data,'vvfpfile',vvfpfile,'badchanfile',badchanfile,...
        'scale_max',field_scale_max,...
        'time0',time0','time1',time1,'view_angle',view_angle,...
        'radius',field_radius,'angle_offset',angle_offset,...
        'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
        'plotsize',field_plotsize,'usemags',usemags);
      fprintf('%s: plotting fitted EEG potentials for area %s...\n',...
        mfilename,retmap.areas(a).name);
      rc_rc_plot_ret_pots(fit_data,'surffile',surffile,'badchanfile',badchanfile,...
        'trans',T_head2mri,'scale_max',pot_scale_max,...
        'time0',time0','time1',time1,'view_angle',view_angle,...
        'electrodes',electrodes,'radius',pot_radius,'angle_offset',angle_offset,...
        'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
        'plotsize',pot_plotsize);
      fprintf('%s: plotting dipole orientations...\n',mfilename);
      rc_plot_dips_circ(prefix,a,plane);
      set(gca,'FontName',fontname)
      if titleflag
        ax=axes('Units','Normal','Position',[0.075 0.075 0.89 0.89],'Visible','off');
        set(get(ax,'Title'),'Visible','on');
        title([sprintf('Expected %s MEG Retinotopy Contour Plots t=%0.1f-%0.1f',...
          retmap.areas(a).name,time0,time1)]);
        set(gcf,'position',[50 50 imgsize imgsize]);
      end;
      outstem = sprintf('contplots-%s-fit-%s-%s_t%03d-%03d.tif',...
        prefix,retmap.areas(a).name,plane,time0,time1);
      print('-dtiff',[outstem '.tif']);
      mmil_printeps(gcf,[outstem '.eps']);
    end;
  end;
end;

close(gcf);
