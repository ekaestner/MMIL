function rc_plot_RCSE_pots(surffile,prefix,time0,time1,view_angle)
%function rc_plot_RCSE_pots(surffile,prefix,[time0],[time1],[view_az, view_el])%
%
% Purpose: wrapper for rc_plot_ret_pots
%
% Required Input:
%   surffile: freesurfer compatible tri file
%     on which to display EEG potentials
%     e.g. outer scalp surface
%   prefix: RCSE prefix
%
% Optional Input:
%   time0: start of time window to average (msec)
%     {default = 80}
%   time1: end of time window to average (msec)
%     {default = 90}
%   view_angle: vector of [aximuth,elevation] for plot view
%     {default = [0 0]}
%
% Early Mod: 10/06/07 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: use varargin?
%% todo: combine with rc_plot_RCSE_fields?

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('time0','var') || isempty(time0), time0 = 80; end;
if ~exist('time1','var') || isempty(time1), time1 = 90; end;
if ~exist('view_angle','var') || isempty(view_angle), view_angle = [0 0]; end

scale_max = 2.5;
radius = 0.3;
plotsize = 0.12;
electrodes = 'off';

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_forward.mat',prefix);
if ~exist(matfile,'file')
  matfile = parms.forward_matfile;
end;
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

badchanfile = parms.badchanfile;
T_mri2head = parms.trans;
T_head2mri = inv(T_mri2head);

num_locs = retmap.num_locs;
display_locs = retmap.cond_order;
orig_num_locs = retmap.orig_num_locs;
angle_offset = 360/(orig_num_locs*2);

rc_plot_ret_pots(avg_data,'surffile',surffile,'badchanfile',badchanfile,...
  'trans',T_head2mri,'scale_max',scale_max,...
  'time0',time0','time1',time1,'view_angle',view_angle,...
  'electrodes',electrodes,'radius',radius,'angle_offset',angle_offset,...
  'display_locs',display_locs,'actual_num_locs',orig_num_locs,...
  'plotsize',plotsize);

