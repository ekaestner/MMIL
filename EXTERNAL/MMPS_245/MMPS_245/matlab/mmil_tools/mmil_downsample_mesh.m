function mmil_downsample_mesh(fname_in,fname_out,varargin)
%function mmil_downsample_mesh(fname_in,fname_out,[options])
%
% Purpose: smooth and downsample a surface mesh
%
% Required Input:
%   fname_in: name of input mesh file in tri format
%   fname_out: output file name
%
% Optional Input:
%   'dmesh': number of vertices for downsampled mesh
%     {default = 2000}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  10/03/13 by Don Hagler
% Last Mod: 10/27/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{,...
  'dmesh',2000,[],...
  'forceflag',false,[false true],...
...
  'lambda',0.5,[],...
  'mode',1,[],... % smooth only in normal direction
  'itt',1,[],... % number of iterations
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_out,'file') || parms.forceflag
  tmp = fs_read_trisurf(fname_in);
  nverts = tmp.nverts;
  surf = [];
  surf.vertices = tmp.vertices;
  surf.faces = tmp.faces;
  surf = smoothpatch(surf,parms.mode,parms.itt,parms.lambda);
  r = parms.dmesh / nverts;
  if r < 1
    surf = reducepatch(surf,r);
    surf = smoothpatch(surf,parms.mode,parms.itt,parms.lambda);
  end
  surf.nverts = size(surf.vertices,1);
  surf.nfaces = size(surf.faces,1);
  fs_write_trisurf(surf,fname_out);
end;


