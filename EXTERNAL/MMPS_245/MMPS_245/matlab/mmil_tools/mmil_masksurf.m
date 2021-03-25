function surf = mmil_masksurf(vol_in,varargin)
%function surf = mmil_masksurf(vol_in,[options])
%
% Purpose: create a surface from a mask volume
%
% Usage: surf = mmil_masksurf(fname_mask,'key1', value1,...);
%
% Required Input:
%   vol_in: input mask volume
%
% Optional Input:
%  'fname_in': full path of input mask volume file (e.g. brain.mgz, aseg.mgz)
%    If supplied, vol_in is ignored and fname_in is loaded instead
%    {default = []}
%  'fname_out': full path of output mask surface file (.tri format)
%    If supplied, masksurf is saved as fname_out
%    {default = []}
%  'nverts': approximate number of vertices in output surface
%    {default = 5000}
%  'convert_flag': [0|1] convert surf struct from ctx/LPH to FS/RAS
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output file
%    ignored if fname_out is empty
%    {default = 0}
%
% Created:  06/08/09 by Don Hagler
% Rcnt Mod: 11/15/11 by Vijay Venkatraman
% Last Mod: 12/29/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

surf = [];
if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_in',[],[],...
  'fname_out',[],[],...
  'nverts',5000,[],... % approximate
  'convert_flag',true,[false true],...
  'forceflag',false,[false true],...
});

if ~isempty(parms.fname_out) &&...
   exist(parms.fname_out,'file') && ~parms.forceflag
  if nargout>0
    surf = fs_read_trisurf(parms.fname_out);
  end;
  return;
end;

if ~isempty(parms.fname_in)
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [vol_in,mr_parms] = ctx_load_mgh(parms.fname_in);
else
  mr_parms = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create surface from mask

fprintf('%s: finding surface...\n',mfilename);
surf = isosurface(vol_in.imgs,0.5);
% check if multiple connected components
fprintf('%s: checking surface...\n',mfilename);
surf = preprocessQ(surf);
if regexp(surf.CC,'MULTIPLE CONNECTED COMPONENTS!')
  fprintf('%s: extracting single connected component...\n',mfilename);
  surf = extractpatchCC(surf);    
end;
% resize tesselation to have nfaces faces
fprintf('%s: resampling surface...\n',mfilename);
R = parms.nverts/surf.V;
surf = reducepatch(surf,R);
% get coordinates in scanner coordinates (from voxel indices)
fprintf('%s: transforming coordinates...\n',mfilename);
surf.vertices = ApplyMvxl2lphCoordlist(surf.vertices,vol_in.Mvxl2lph);
% find neighbors, etc.
fprintf('%s: checking surface...\n',mfilename);
surf = preprocessQ(surf);
if regexp(surf.PREPROCESSED,'NON-MANIFOLD!!!')
  fprintf('%s: fixing topology...\n',mfilename);
  surf = topoFixer(surf);
  surf = preprocessQ(surf);
end;

if parms.convert_flag
  surf = convert_surf(surf,vol_in.lphcent);
end;

% save as tri file for viewing with tkmedit and tksurface
if ~isempty(parms.fname_out)
  fprintf('%s: saving file...\n',mfilename);
  fs_write_trisurf(surf,parms.fname_out);
end;

return;


function fssurfstruct = convert_surf(surfstruct,lphcent)
  fssurfstruct = [];
  fssurfstruct.nverts = surfstruct.V;
  fssurfstruct.nfaces = surfstruct.F;
  fssurfstruct.faces = surfstruct.faces;
  fssurfstruct.vertices = surfstruct.vertices;
  fssurfstruct.vertices(:,1) = -fssurfstruct.vertices(:,1) + lphcent(1);
  fssurfstruct.vertices(:,2) = -fssurfstruct.vertices(:,2) + lphcent(2);
  fssurfstruct.vertices(:,3) = fssurfstruct.vertices(:,3) + lphcent(3);
return;

