function [side,centroid,ncentroid,area,lside,inout,vertex] = ...
    tri_analysis(x,y,z,R,geo,scale,verbose);
% TRI_ANALYSIS analyze triangulear data.
% function [side,centroid,ncentroid,area,lside,inout,vertex] = ...
%     tri_analysis(x,y,z,R,geo,scale,verbose);
% Input: x,y,z are the triangular vertices matrices, 3 x nT
% Optional:
%  R is the vertices matrix, mR x 3
%  geo is the index of vertices in R for each triangle, 3 x nT
% Output: Assume there are mR vertices and nT triangles
%  side: nT x 9, each row corresponds to one triangle, first three columns
%   are the displacement vector of side one (geo2 - geo1), then geo3 - geo2,
%   then geo1 - geo3, so we proceed around the triangle
%  centroid: nT x 3, each row is the center of the triangle
%  ncentroid: nT x 3, each row is the unit vector pointing normal to the
%    triangle 
%  area: nT x 1, the area of each triangle
%  lside: 3 x nT, 1st row is length of side 1 (1 to 2),
%    2nd row of side 2 (2 to 3), 3rd row of side 3 (3 to 1).
%  inout: 1 x nT, + indicates pointed "out", -1 indicates "in"
%  vertex: mR rows, each row gives the indices to the triangles which share
%    that vertex (must give R and geo).


if(exist('verbose')~=1),
  verbose = 0;			% silent running
end
v1 = [x(1,:)' y(1,:)' z(1,:)'];	% vertice # 1
v2 = [x(2,:)' y(2,:)' z(2,:)'];	% vertice # 2
v3 = [x(3,:)' y(3,:)' z(3,:)'];	% vertice # 3
r12 = v2 - v1;			% from 1 to 2
r23 = v3 - v2;
r31 = v1 - v3;
side = [r12 r23 r31];

lside = zeros(size(x));	% length of sides

lside(1,:) = rownorm(r12)'; 	% length of side 1
lside(2,:) = rownorm(r23)'; 	% length of side 2
lside(3,:) = rownorm(r31)'; 	% length of side 3

centroid = [mean(x)' mean(y)' mean(z)'];
% form normal to the surface of the triangle
ncentroid = cross_mm(r12,-r31);
area = rownorm(ncentroid)/2; % area of each triangle
ncentroid = ncentroid./(rownorm(ncentroid)*[1 1 1]); % normalize

% is the orientation of the triangle in or out?
% dot product the centroid with the average of the vertices, if 
%  positive, then vertices and centroid have the same polarity.
inout = sign(sum(centroid.' .* ncentroid.'));

% each vertex probably joins six triangles, initialize as 1, will grow
%  automatically
if(exist('R')==1),		% user gave R and geo as well
  vertex = zeros(size(R,1),1); 	
  for i = 1:size(R,1), 		% for each vertex
    [tmp,tri_num] = find(geo==i); % what triangles have this vertice
    vertex(i,[1:length(tri_num)]) = tri_num';
  end
end

if(verbose),
  narea = area/max(area)*scale;	% normalized
  xv = [centroid(:,1)';centroid(:,1)'+ [ncentroid(:,1).*(narea)]'];
  yv = [centroid(:,2)';centroid(:,2)'+ [ncentroid(:,2).*(narea)]'];
  zv = [centroid(:,3)';centroid(:,3)'+ [ncentroid(:,3).*(narea)]'];
  figure(windfind('Triangle analysis'));
  clf
  fill3(x,y,z,z)
  colormap(hsv)
  hold on
  plot3(xv,yv,zv,'r-')
  hold off
  axis('equal'),axis('square')
  drawnow
end

return
