function censor_tseries = mmil_censor_motion(motion_tseries,censor_thresh,contig_thresh,motion_radius,absflag,nodyflag)
%function censor_tseries = mmil_censor_motion(motion_tseries,[censor_thresh],[contig_thresh],[motion_radius][absflag],[nodyflag])
%
% Purpose: determine time points to censor based on motion
%
% Required Parameters:
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%         also called: dL, dP, dS, pitch, yaw, roll
%                      to match output of mmil_load_motion_1D
%     with translations in mm and rotations in degrees clockwise
%
% Optional parameters:
%   censor_thresh: framewise displacement censoring threshold (mm)
%     {default = 0.2}
%   contig_thresh: mininum number of contiguous frames required
%     {default = 5}
%   motion_radius: approximate radius of typical head (mm)
%     for calculating distance from an angle of rotation
%     {default = 50}
%   absflag: calculate motion as sum of absolute values of each component
%     otherwise, calculate square root of sum of squares
%     {default = 1}
%   nodyflag: calculate motion ignoring dy (translation in y)
%     {default = 0}
%
% Output:
%   censor_tseries: vector of 0s and 1s [nframes,1]
%     with 1 indicating censored frame
%
% Created:  10/16/17 Don Hagler
% Last Mod: 10/16/17 Don Hagler
%
% Based on code from Eric Earl at OHSU (earl@ohsu.edu)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

censor_tseries = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('censor_thresh','var') || isempty(censor_thresh)
  censor_thresh = 0.2;
end;
if ~exist('contig_thresh','var') || isempty(contig_thresh)
  contig_thresh = 5;
end;
if ~exist('motion_radius','var') || isempty(motion_radius)
  motion_radius = 50;
end;
if ~exist('absflag','var') || isempty(absflag), absflag = 1; end;
if ~exist('nodyflag','var') || isempty(nodyflag), nodyflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[motion_stats,motion_fd] = mmil_motion_stats(motion_tseries,motion_radius,absflag,[],nodyflag);

censor_tseries = logical(motion_fd > censor_thresh);

% code from eric earl, earl@ohsu.edu
if contig_thresh>0
  contig_tseries = zeros(size(censor_tseries)); 
  contig_groups = bwlabel(logical((censor_tseries-1)*-1)); 
  for group = 1:max(contig_groups)
    if sum(contig_groups == group) < contig_thresh
      contig_tseries(contig_groups == group) = 1;
    else
      contig_tseries(contig_groups == group) = 0;
    end
  end
  censor_tseries = logical(censor_tseries + contig_tseries);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

