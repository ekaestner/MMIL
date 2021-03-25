function coords = ts_roi_centmass(surf,v,w);
%function coords = ts_roi_centmass(surf,[v],[w]);
%
% Required Input:
%  surf: surface structure containing:
%     nverts: number of vertices
%     nfaces: number of faces (triangles)
%     faces:  vertex numbers for each face (3 corners)
%     vertices: x,y,z coordinates for each vertex
%
% Optional Input:
%   v: vector of vertex numbers for ROI (optional)
%   w: vector of vertex values corresponding to v (optional)
%
% Output:
%  coords: x,y,z coordinates of center of mass of verts
%               optionally weighted by vals
%                  
%
% see also: fs_read_surf, fs_read_wfile
% 
% Created:   01/08/10 by Don Hagler
% Last Mod:  01/08/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
if ~exist('v','var') || isempty(v), v = 1:surf.nverts; end;
if ~exist('w','var') || isempty(w), w = ones(length(v),1); end;
coords = [];

X = surf.vertices(v,1);
Y = surf.vertices(v,2);
Z = surf.vertices(v,3);

sumw = sum(w);
x = sum(X.*w)/sumw;
y = sum(Y.*w)/sumw;
z = sum(Z.*w)/sumw;

coords = [x,y,z];

