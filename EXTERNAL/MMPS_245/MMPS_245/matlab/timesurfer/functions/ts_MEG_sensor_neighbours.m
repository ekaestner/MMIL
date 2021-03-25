function cfg = ts_MEG_sensor_neighbours (cfg,sens,thresh,zshift)
%function cfg = ts_MEG_sensor_neighbours (cfg,sens,thresh,zshift)
%
% Required input:
%  cfg    -- is a structure specifying the configuration for MC clustering
%  sens   -- is sensor_info from any valid TimeSurfer structure
%  thresh -- is an angular threshold for defining neighbours
%
% Output:
%   cfg is a structure containing:
%     neighbours.label       -- a string specifying the label for each vertex
%     neighbours.neighblabel -- a cell array of strings specifying the
%     labels for each neighbour associated with the vertex defined in
%     neighbours.label
%
% created:        05/05/10 Jason Sherfey
% last modified:  05/05/10 Jason Sherfey

if nargin < 4, zshift = 0; end
if nargin < 3, help(mfilename); return; end
  
theta = ts_BetweenSensorAngles(sens,zshift);
cfg.neighbours = [];
for k = 1:size(theta,1)
  ch  = find(theta(k,:)<thresh);
  cfg.neighbours{k}.label       = sens(k).label;
  cfg.neighbours{k}.neighblabel = {sens(ch).label};
end
return