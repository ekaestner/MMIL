function mmil_warp_freesurfer_surface(fname_surf_in,fname_surf_out,fname_dx,fname_dy,fname_dz,forceflag);
%function mmil_warp_freesurfer_surface(fname_surf_in,fname_surf_out,fname_dx,fname_dy,fname_dz,[forceflag]);
%
% Required Input:
%   fname_surf_in: filename of input surface
%   fname_surf_out: filename of output surface
%   fname_dx: filename of dx (displacement in x) volume
%   fname_dy: filename of dy (displacement in y) volume
%   fname_dz: filename of dz (displacement in z) volume
%
% Optional Input:
%   forceflag: [0|1] toggle overwrite if fname_surf_out exists
%     {default = 0}
%
% Created:  06/20/07 Don Hagler
% Rcnt Mod: 03/23/12 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,5), return; end;

if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;

if exist(fname_surf_out,'file') & ~forceflag, return; end;

[fpath,fstem,fext] = fileparts(fname_surf_in);
fprintf('%s: loading input surface file %s...\n',mfilename,fname_surf_in);
if strcmp(fext,'.tri')
  surf = fs_read_trisurf(fname_surf_in);
else
  surf = fs_read_surf(fname_surf_in);
end;

try
  fprintf('%s: loading displacement file %s...\n',mfilename,fname_dx);
  vol_dx = ctx_load_mgh(fname_dx);
  fprintf('%s: loading displacement file %s...\n',mfilename,fname_dy);
  vol_dy = ctx_load_mgh(fname_dy);
  fprintf('%s: loading displacement file %s...\n',mfilename,fname_dz);
  vol_dz = ctx_load_mgh(fname_dz);
catch
  fprintf('%s: ERROR: failed to load displacement files\n',mfilename);
  return;
end;

fprintf('%s: warping surface coordinates...\n',mfilename);
surf.vertices = warp_coords(surf.vertices,vol_dx,vol_dy,vol_dz);

[fpath,fstem,fext] = fileparts(fname_surf_out);
fprintf('%s: writing output surface file %s...\n',mfilename,fname_surf_out);
if strcmp(fext,'.tri')
  fs_write_trisurf(surf,fname_surf_out);
else
  fs_write_surf(surf,fname_surf_out);
end;

fprintf('%s: finished.\n',mfilename);

return;


function coords = warp_coords(coords,vol_dx,vol_dy,vol_dz)
  tmp = zeros(size(coords,1),4);
  tmp(:,1:3) = coords;
  tmp(:,4)=1;
  % transform RAS to LPH
  tmp = (M_RAS_TO_LPH*tmp')';
  dl = vol_getvxlsval(tmp, vol_dx, eye(4), 2);
  dp = vol_getvxlsval(tmp, vol_dy, eye(4), 2);
  dh = vol_getvxlsval(tmp, vol_dz, eye(4), 2);
  tmp(:,1) = tmp(:,1)+dl;
  tmp(:,2) = tmp(:,2)+dp;
  tmp(:,3) = tmp(:,3)+dh;
  % transform LPH to RAS
  tmp = (M_LPH_TO_RAS*tmp')';
  coords = tmp(:,1:3);
return;

