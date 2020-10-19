function weights = rc_create_smooth_area_weights(surf,hemi,centers,areamask,smooth_fwhm)
%function weights = rc_create_smooth_area_weights(surf,hemi,centers,areamask,smooth_fwhm)
%
% Required Input:
%   surf: freesurfer surface structure containing:
%     nverts: number of vertices
%     nfaces: number of faces (triangles)
%     faces:  vertex numbers for each face (3 corners)
%     coords: x,y,z coords for each vertex
%     (see fs_load_subj and fs_read_surf)
%  hemi: cortical hemisphere ('lh' or 'rh')
%  centers: vector of vertex numbers defining center vertex for each
%    stimulus location for a single visual area
%    Note: vertex numbers should be 1-based
%
% Optional Input:
%  areamask: vector of vertex numbers to mask a single visual area
%    if empty or ommitted, weights will not be confined to a single area
%    Note: vertex numbers should be 1-based
%    {default = []}
%  smooth_fwhm: full width half max blurring kernel on surface (mm)
%    N = round(fwhm/1.25).^2) = integer number of smoothing steps
%    {default = 10 mm}
%
% Output:
%   weights: sparse matrix containing weights at each vertex for each center vertex
%
% Created: 05/11/07 Don Hagler
% Lst Mod: 02/19/11 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: use varagin?

if ~mmil_check_nargs(nargin,3), return; end;

if ~exist('areamask','var'), areamask = []; end;
if ~exist('smooth_fwhm','var') | isempty(smooth_fwhm), smooth_fwhm = 10; end;

smooth_N = round((smooth_fwhm/1.25).^2); % FWHM ~ 1.25*sqrt(N)
initval = 10;

weights = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for errors
if ~ismember(hemi,{'lh','rh'})
  error('hemi must be ''lh'' or ''rh''');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nstims = length(centers);
nverts = surf.nverts;
weights = sparse(nstims,nverts);

for i=1:nstims
  v = centers(i);
  w = zeros(nverts,1);
  w(v) = initval;
  w = fs_smooth(surf,w,smooth_N);
  if ~isempty(areamask)
    w_mask = zeros(nverts,1);
    w_mask(areamask) = 1;
    w = w.*w_mask;
  end;
  v = find(w~=0);
  w = w(v);
  w = rc_norm_weights(w,1);  % norm to max
  weights(i,v) = w;
end;

return;

