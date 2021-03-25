function [x,y,z] = rc_roi_centmass(surf,v,w);
%function [x,y,z] = rc_roi_centmass(surf,v,w);
%
% input:
%   surf is a structure containg:
%     nverts: number of vertices
%     nfaces: number of faces (triangles)
%     faces:  vertex numbers for each face (3 corners)
%     vertices: x,y,z coordinates for each vertex
%  v: vector of vertex numbers for ROI
%  w: vector of vertex values corresponding to v (optional)
%
% output:
%  x,y,z: x,y,z coordinates of ROI center of mass
%
% see also: fs_read_surf, fs_read_wfile
%
% Early Mod: 01/08/10 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('w','var') || isempty(w), w = ones(length(v),1); end;

X = surf.vertices(v,1);
Y = surf.vertices(v,2);
Z = surf.vertices(v,3);
V = w;

[x,y,z] = rc_centmass(X,Y,Z,V);

return;

