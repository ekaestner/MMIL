function [v,vdist] = fs_nearest_verts(points,surf)
%function [v,vdist] = fs_nearest_verts(points,surf)
%
% Required Input:
%   points: Nx3 matrix of RAS coordinates (N = number of points)
%   surf: structure containing:
%     nverts: number of vertices
%     nfaces: number of faces (triangles)
%     vertices: x,y,z coordinates for each vertex (nverts x 3)
%     faces:  vertex numbers for each face (3 corners)
%
% Output:
%   v: vector of 1-based vertex numbers for each point
%   vdist: distance from each point to nearest vertex
%
% see also: fs_read_surf
%
% Created:  11/03/11 by Don Hagler
% Last Mod: 11/03/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
v = [];
vdist = [];

if size(points,2)~=3, error('size of points must be Nx3'); end;
N = size(points,1);

% find nearest vertex to each point
v = zeros(1,N);
vdist = Inf(1,N);
for i=1:N
  point = points(i,:);
  point_dist = sqrt(sum(((ones(surf.nverts,1)*point - surf.vertices ).^2),2));
  [min_dist,nearest_v] = min(point_dist);
  v(i) = nearest_v;
  vdist(i) = min_dist;
end;

