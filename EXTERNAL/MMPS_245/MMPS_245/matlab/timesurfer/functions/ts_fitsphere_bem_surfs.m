function [cen_sph, radii]=ts_fitsphere_bem_surfs(bem_surf_files,T_mri2head);
%function [cen_sph, radii]=ts_fitsphere_bem_surfs(bem_surf_files,[T_mri2head]);
%
% ts_fitsphere_bem_surfs: fit a sphere to FreeSurfer compatible tri files
%
% Required input:
%   bem_surf_files: cell array of 3 names of FreeSurfer compatible tri files
%     files should be in this order:
%       (1) inner skull
%       (2) outer skull
%       (3) outer scalp
%     optionally, can supply a single file name for the inner skull if
%       1-shell BEM is used
%
% Optional input:
%   T_mri2head: 4x4 matrix specifying mri2head transformation
%    {default=identity matrix}
%
% Output:
%   cen_sph: vector of x,y,z coordinates of best fit center of first bem surface
%   radii: vector of radii (in mm) of fitted sphere for each bem surface
%
% Created 08/09/06 by Don Hagler
% Last Modified 08/09/06 DH
%

if nargin<1
  help(mfilename)
  return;
end;

cen_sph = [];
radii = [];

if isempty(bem_surf_files)
  fprintf('%s: error: bem_surf_files is empty\n',mfilename);
  return;
end;
if iscell(bem_surf_files)
  num_bem_surfs=length(bem_surf_files);
else
  bem_surf_files={bem_surf_files};
  num_bem_surfs=1;
end;

if ~exist('T_mri2head','var'), T_mri2head=[]; end;
if isempty(T_mri2head), T_mri2head=eye(4); end;

[verts,faces]=ts_load_bem_surfs(bem_surf_files);

cen_sph_mri = fminsearch(@(cen_sph_mri) sphererr(cen_sph_mri,verts{1}),mean(verts{1})')';
fprintf('%s: center of sphere (MRI space) = [%0.4f,%0.4f,%0.4f]\n',...
  mfilename,cen_sph_mri);

% apply coordinate transformation to get to sensor-space
cen_sph=T_mri2head*[cen_sph_mri'/1000;1]; % in meters
cen_sph=cen_sph(1:3)'*1000; % back to mm and reshape

fprintf('%s: center of sphere (HEAD space) = [%0.4f,%0.4f,%0.4f]\n',...
  mfilename,cen_sph);

for s=1:num_bem_surfs
  [err,r] = sphererr(cen_sph_mri',verts{s});
  fprintf('%s: surf %d: radius=%0.4f mm  error=%0.4f\n',...
    mfilename,s,r,err);
  radii(s) = r;
end;

return;

