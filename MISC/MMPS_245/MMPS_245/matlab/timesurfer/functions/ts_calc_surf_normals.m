function normals=ts_calc_surf_normals(surf,angle_flag)
%function normals=ts_calc_surf_normals(surf,angle_flag)
%
% Purpose: calculate normal vector for each vertex of a surface mesh
%
% Required Parameters:
%   surf: struct containing vertices and faces
%
% Optional Parameters:
%   angle_flag: [0|1] weight face normals by face angles
%     {default = 1}
%
% Created:    06/25/15 by Don Hagler
% Last Mod:   06/25/15 by Don Hagler
%

% adapted from:
%   patchnormals
%     by Dirk-Jan Kroon, June 2, 2009, MATLAB Central File Exchange
%   normals
%     by Robert Oostenveld, 2002-2007, FieldTrip 20080624

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normals = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('angle_flag','var') | isempty(angle_flag), angle_flag = 1; end;

normals = zeros(surf.nverts,3);

% vertex indices for each triangle
Fa = surf.faces(:,1);
Fb = surf.faces(:,2);
Fc = surf.faces(:,3);

% calculate edge vectors
e1 = surf.vertices(Fb,:) - surf.vertices(Fa,:);
e2 = surf.vertices(Fc,:) - surf.vertices(Fb,:);
e3 = surf.vertices(Fc,:) - surf.vertices(Fa,:);

% calculate face normals
face_normals = cross(e1,e3);

if angle_flag
  % normalize edge vectors
  e1_norm = bsxfun(@rdivide,e1,sqrt(sum(e1.^2,2)));
  e2_norm = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,2)));
  e3_norm = bsxfun(@rdivide,e3,sqrt(sum(e3.^2,2)));

  % calculate angle of face seen from vertices
  face_angles = [acos(dot(e1_norm',e3_norm'));...
                 acos(dot(-e2_norm',e1_norm'));...
                 acos(dot(e3_norm',e2_norm'))]';

  % replace NaNs with 0s (caused by neighbors having identical coordinates)
  face_angles(isnan(face_angles)) = 0;

  % calculate vertex normals
  for i=1:surf.nfaces
    normals(Fa(i),:) = normals(Fa(i),:)+face_normals(i,:)*face_angles(i,1);
    normals(Fb(i),:) = normals(Fb(i),:)+face_normals(i,:)*face_angles(i,2);
    normals(Fc(i),:) = normals(Fc(i),:)+face_normals(i,:)*face_angles(i,3);
  end
else
  % calculate vertex normals
  for i=1:surf.nfaces
    normals(Fa(i),:) = normals(Fa(i),:)+face_normals(i,:);
    normals(Fb(i),:) = normals(Fb(i),:)+face_normals(i,:);
    normals(Fc(i),:) = normals(Fc(i),:)+face_normals(i,:);
  end
end;

% calculate unit vectors
normals = bsxfun(@rdivide,normals,sqrt(sum(normals.^2, 2)));

