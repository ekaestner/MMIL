function [mean_motion,mean_trans,mean_rot] = ...
  mmil_mean_motion(motion_tseries,motion_radius)
%function [mean_motion,mean_trans,mean_rot] = ...
%   mmil_mean_motion(motion_tseries,[motion_radius])
%
% Purpose: calculate mean frame-to-frame motion
%
% Required Parameters:
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%                      to match output of mmil_load_motion_1D
%         also called: dL, dP, dS, pitch, yaw, roll
%
% Optional parameters:
%  motion_radius: approximate radius of typical head (mm)
%    for calculating distance from an angle of rotation
%    {default: 50}
%
% Output:
%   mean_motion: mean frame-to-frame motion (including translation and rotation)
%   mean_trans: mean frame-to-frame translation
%   mean_rot: mean frame-to-frame rotation
%
% Created:  02/12/13 Don Hagler
% Last Mod: 03/23/16 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_motion = [];
mean_trans = [];
mean_rot = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('motion_radius','var') || isempty(motion_radius)
  motion_radius = 50;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(motion_tseries,1) > 1
  % calculate absolute difference in head position or rotation
  d = abs(diff(motion_tseries,1,1));
else
  d = abs(motion_tseries);
end;

% calculate approximate displacements on sphere about the size of a head
d(:,4:6) = d(:,4:6)*pi*motion_radius/180;

% calculate mean motion
mean_motion = mean(sum(d,2));

% calculate mean translation
mean_trans = mean(sum(d(:,1:3),2));

% calculate mean rotation
mean_rot = mean(sum(d(:,4:6),2));

