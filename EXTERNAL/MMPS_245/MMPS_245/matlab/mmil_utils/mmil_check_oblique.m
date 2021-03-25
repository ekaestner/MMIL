function [oblique_flag,oblique_angle]=mmil_check_oblique(fname)
%function [oblique_flag,oblique_angle]=mmil_check_oblique(fname)
%
% Purpose: determine whether scan is oblique
%          and returns oblique angle (degrees clockwise)
%
% Required Input:
%   fname: name of mgh/mgz file
%
% Created:  08/03/11 by Vijay Venkatraman
% Last Mod: 08/28/12 by Don Hagler
%

if ~exist(fname,'file')
  error('%s: %s File doesnt exist',mfilename,fname);
end;

[vol,M] = fs_load_mgh(fname,[],[],1);
[orient,oblique_flag] = fs_read_orient(fname,M);

voxel_sizes = sqrt(sum(M(1:3,1:3).^2,1));
M0 = mmil_construct_M('orient',orient,'scale',voxel_sizes);

xhat=mmil_M_mat2vec(M*inv(M0));
oblique_angle = xhat(4:6); 

return;

