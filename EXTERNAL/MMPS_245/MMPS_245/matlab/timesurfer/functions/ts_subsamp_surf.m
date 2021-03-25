function [subsurf,surf,decdips]=ts_subsamp_surf(surffile,subfact,subsurffile)
%function [subsurf,surf,decdips]=ts_subsamp_surf(surffile,[subfact],[subsurffile])
%
% load surf from surffile and subsample with reducepatch
%
% Required Parameters:
%   surffile: full path name of freesurfer binary surface file
%
% Optional Parameters:
%   subfact: subsampling factor
%     0.02 takes full surface with ~160,000 verts to ~2500 verts
%     {default: 0.02}
%    subsurffile: full path name of output freesurfer binary surface file
%     If empty or ommitted, no output file is created
%     {default: []}
%
% Output:
%   subsurf: subsampled surface struct
%   surf: original surface struct
%   decdips: vector of 0's and 1's, indicating which vertices from surf
%     are found in subsurf
%
% Created:   01/20/08 by Don Hagler
% Last Mod:  01/08/10 by Don Hagler
%

surf = [];
subsurf = [];
decdips = [];

if ~exist('subfact','var') | isempty(subfact), subfact = 0.02; end;
if ~exist('subsurffile','var'), subsurffile = []; end;
smf = 10^-5;

% load surface
if ~exist(surffile,'file')
  error('%s: file %s not found',mfilename,surffile);
end;
surf = fs_read_surf(surffile);
surf = fs_find_neighbors(surf);
surf = fs_calc_triarea(surf);

% subsample surface
[subsurf.faces,subsurf.vertices]=reducepatch(surf.faces,surf.vertices,subfact);
subsurf.nverts = size(subsurf.vertices,1);
subsurf.nfaces = size(subsurf.faces,1);
subsurf = fs_find_neighbors(subsurf);
subsurf = fs_calc_triarea(subsurf);

% get decdips vector
decdips = zeros(surf.nverts,1);
for k=1:surf.nverts
  v = surf.vertices(k,:);
  tmp = sum(abs(subsurf.vertices - ones(subsurf.nverts,1)*v),2);
  if length(find(tmp<smf))>0
    decdips(k) = 1;
  end;
end;    

% save as freesurfer surf files
if ~isempty(subsurffile)
  fs_write_surf(subsurf,subsurffile);
end;
