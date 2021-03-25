function [motion_stats,motion_fd,deriv] = mmil_motion_stats(motion_tseries,motion_radius,absflag,thresholds,nodyflag)
%function [motion_stats,motion_fd,deriv] = mmil_motion_stats(motion_tseries,[motion_radius][absflag],[thresholds],[nodyflag])
%
% Purpose: calculate statistics from motion time course
%
% Required Parameters:
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%         also called: dL, dP, dS, pitch, yaw, roll
%                      to match output of mmil_load_motion_1D
%     with translations in mm and rotations in degrees clockwise
%
% Optional parameters:
%  motion_radius: approximate radius of typical head (mm)
%    for calculating distance from an angle of rotation
%    {default = 50}
%  absflag: calculate motion as sum of absolute values of each component
%    otherwise, calculate square root of sum of squares
%    {default = 1}
%  thresholds: vector of motion thresholds (mm)
%    used to calculate percent and time of scan below threshold
%    {default = [0.2,0.3,0.4]}
%  nodyflag: calculate motion ignoring dy (translation in y)
%    {default = 0}
%
% Output:
%   motion_stats: struct containing calculated motion statistics with fields:
%     max_rz:      maximum absolute rotz (relative to baseline)
%     max_rx:      maximum absolute rotx (relative to baseline)
%     max_ry:      maximum absolute roty (relative to baseline)
%     max_dz:      maximum absolute rotz (relative to baseline)
%     max_dx:      maximum absolute rotx (relative to baseline)
%     max_dy:      maximum absolute roty (relative to baseline)
%     mean_trans:  mean frame-to-frame translation
%     max_trans:   maximum frame-to-frame translation
%     mean_rot:    mean frame-to-frame rotation
%     max_rot:     maximum frame-to-frame rotation
%     mean_motion: mean frame-to-frame motion (including translation and rotation)
%     std_motion:  standard deviation of frame-to-frame motion
%     mode_motion: mode of frame-to-frame motion
%     med_motion:  median of frame-to-frame motion
%     min_motion:  minimum of  frame-to-frame motion
%     max_motion:  maximum of  frame-to-frame motion
%     mean_accel:  mean of derivative of frame-to-frame motion
%     subthresh_nvols: number of volumes with subthreshold frame-wise displacement
%     subthresh_perc: percent of volumes with subthreshold frame-wise displacement
%     thresholds: vector of motion thresholds
%     nvols: number of volumes
%
%   motion_fd: vector of framewise displacements in mm [nframes,1]
%
%   deriv: first derivative of motion time courses [nframes,6]
%
% Created:  02/12/13 Don Hagler
% Last Mod: 08/23/17 Dani Cornejo 
% Last Mod: 10/16/17 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

motion_stats = [];
motion_fd = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('motion_radius','var') || isempty(motion_radius)
  motion_radius = 50;
end;
if ~exist('absflag','var') || isempty(absflag), absflag = 1; end;
if ~exist('thresholds','var') || isempty(thresholds)
  thresholds = [0.2 0.3 0.4];
end;
if ~exist('nodyflag','var') || isempty(nodyflag), nodyflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate maximum translations relative to baseline
motion_stats.max_dx = max(abs(motion_tseries(:,1)));
motion_stats.max_dy = max(abs(motion_tseries(:,2)));
motion_stats.max_dz = max(abs(motion_tseries(:,3)));
motion_stats.max_rx = max(abs(motion_tseries(:,4)));
motion_stats.max_ry = max(abs(motion_tseries(:,5)));
motion_stats.max_rz = max(abs(motion_tseries(:,6)));

% calculate difference in head position or rotation
if size(motion_tseries,1) > 1
  d = diff(motion_tseries,1,1);
  d2 = diff(motion_tseries,2,1);
else
  d = motion_tseries;
  d2 = 0;
end;

% extract rotation estimates
drot = d(:,4:6);

% convert rotations in radians to approximate displacements in mm
%   (on sphere about the size of a head)
d(:,4:6) = d(:,4:6)*pi*motion_radius/180;

deriv = d; 

% extract translation estimates
if ~nodyflag
  dtrans = d(:,1:3);
else
  % exclude y translation
  dtrans = d(:,[1,3]);
  d = d(:,[1,3:6]);
end;

% calculate absolute value
trans = calc_stat(dtrans,absflag);

% calculate overall rotation
rot = calc_stat(drot,absflag);

% calculate mean translation
motion_stats.mean_trans = mean(trans);

% calculate maximum translation
motion_stats.max_trans = max(trans);

% calculate mean rotation
motion_stats.mean_rot = mean(rot);

% calculate maximum rotation
motion_stats.max_rot = max(rot);

% calculate overall motion
motion = calc_stat(d,absflag);
motion_fd = cat(1,0,motion);

% calculate overall acceleration
accel = calc_stat(d2,absflag);

% calculate mean motion
motion_stats.mean_motion = mean(motion);

% calculate std of motion
motion_stats.std_motion = std(motion);

% calculate mode motion
motion_stats.mode_motion = mode(motion);

% calculate median motion
motion_stats.med_motion = median(motion);

% calculate minimum motion
motion_stats.min_motion = min(motion);

% calculate maximum motion
motion_stats.max_motion = max(motion);

% calculate mean acceleration
motion_stats.mean_accel = mean(accel);

% calculate time and percent of scan below each threshold
nthresh = length(thresholds);
nvols = length(motion_fd);
motion_stats.subthresh_nvols = zeros(nthresh,1);
for i=1:nthresh
  motion_stats.subthresh_nvols(i) = sum(motion<thresholds(i));
end;
motion_stats.subthresh_perc = 100*motion_stats.subthresh_nvols/nvols;
motion_stats.thresholds = thresholds;
motion_stats.nvols = nvols;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stat = calc_stat(vec,absflag)
  if absflag
    stat = sum(abs(vec),2);
  else
    stat = sqrt(sum(vec.^2,2));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
