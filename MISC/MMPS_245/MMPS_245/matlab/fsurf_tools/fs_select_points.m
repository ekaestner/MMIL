function [points_sel,ind_sel] = fs_select_points(points,varargin)
%function [points_sel,ind_sel] = fs_select_points(points,[options])
%
% Purpose: select points based on distance to seed point
%
% Required Input:
%   points: npoints x 3 matrix of Talairach coordinates
%
% Optional Input:
%  'seed': seed point coordinates
%    {default = [0 0 0]}
%  'radius': maximum distance (mm) from seed point for selection
%    {default = 50}
%  'MNIflag': [0|1] whether input points are MNI Talaiarach
%      versus true Talaiarach coordinates
%    may be a vector, with one value for each point
%    {default = 1}
%  'seed_MNIflag': [0|1] whether seed point
%      is MNI or true Talairach coordinates
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 1}
%
% Output:
%   points_sel: matrix of selected points (MNI Talairach)
%   ind_sel: index of selected points: points_sel = points(ind_sel,:)
%
% Created:  01/19/16 by Don Hagler
% Last Mod: 01/22/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'seed',[0 0 0],[],...
  'radius',50,[0,100],...
  'MNIflag',true,[false true],...
  'seed_MNIflag',true,[false true],...
  'verbose',true,[false true],...
});
points_sel = []; ind_sel = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(points,2)~=3, error('size of points must be npoints x 3'); end;
npoints = size(points,1);

parms.seed = mmil_rowvec(parms.seed);
if length(parms.seed)~=3, error('seed must have 3 elements'); end;

if length(parms.MNIflag)==1 && npoints>1
  parms.MNIflag = boolean(ones(npoints,1)*parms.MNIflag);
end;
if length(parms.MNIflag)~=npoints
  error('MNIflag must have 1 or %d values',npoints);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% optionally transform from true Talairach to MNI Talairach using fs_tal2mni
if any(~parms.MNIflag)
  % transform points if MNIflag = 1
  ind_tal = find(~parms.MNIflag);
  points(ind_tal,:) = fs_tal2mni(points(ind_tal,:));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find nearest vertex to each point
pdist = sqrt(sum(((repmat(parms.seed,[npoints,1]) - points).^2),2));
ind_sel = find(pdist<parms.radius);
points_sel = points(ind_sel,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

