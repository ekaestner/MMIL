function [gwtype,gwcorr,gwtypelist] = ctx_get_gwtype(vol,varargin)
%function [gwtype,gwcorr,gwtypelist] = ctx_get_gwtype(vol,[options])
%
% Purpose: determine gradwarp signature from an image
%
% Required Input:
%   vol: input volume (ctx format)
%
% Optional Input:
%   'isoctrflag': [0|1] isocenter scanning
%     {default = 1}
%   'gwtypelist': vector of grad warp types to test
%     {default = [2 3 8]}
%   'min_corr': minimum correlation to declare match
%     {default = 0.9}
%   'thresh': threshold applied to image to detect edge
%     {default = 1}
%   'smooth': smoothing sigma applied to binary image
%     {default = 8}
%   'cropx': number of rows cropped on either edge of image in
%     {default = 3}
%   'cropy': number of columns cropped on either edge of image in
%     {default = 30}
%   'cropz': number of slices cropped on either edge of image in
%     {default = 10}
%
% Output:
%   gwtype: best fitting type (empty if none found)
%   gwcorr: correlation for each gwtype in gwtypelist
%   gwtypelist: vector grad warp types tested
% 
% Created:  12/11/10 by Don Hagler
% Last Mod: 05/21/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'isoctrflag',true,[false true],...
  'gwtypelist',[2,3,8],[],...
  'min_corr',0.9,[0.1,0.99],...
  'thresh',1,[0,1e6],...
  'smooth',8,[],...
  'cropx',3,[],...
  'cropy',30,[],...
  'cropz',10,[],...
});

gwtype = [];
gwcorr = [];
gwtypelist = parms.gwtypelist;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% threshold image to detect edge from gradwarp
vol.imgs = 1.0*(vol.imgs>1);

% smooth mask
vol.imgs = mmil_smooth3d(vol.imgs,...
  parms.smooth,parms.smooth,parms.smooth);

% crop image and reshape to vector 
vec = mat2vect(vol.imgs(...
  parms.cropx:(end-parms.cropx+1),...
  parms.cropy:(end-parms.cropy+1),...
  parms.cropz:(end-parms.cropz+1)));

% binarize
vec(vec<0.5) = 0; vec(vec>0.5) = 1;

% test each gwtype
gwcorr = zeros(1,length(gwtypelist));
for j = 1:length(gwtypelist)
  gwtype = gwtypelist(j);  

  % create volume of 1's and apply grad unwarp for this gwtype
  vol_uw = vol; vol_uw.imgs(:) = 1;
  vol_uw = ctx_unwarp_grad(vol_uw,gwtype,2,parms.isoctrflag,1);

  % smooth warped image
  vol_uw.imgs = mmil_smooth3d(vol_uw.imgs,...
    parms.smooth,parms.smooth,parms.smooth);

  % crop image and reshape to vector 
  vec_uw = mat2vect(vol_uw.imgs(...
    parms.cropx:(end-parms.cropx+1),...
    parms.cropy:(end-parms.cropy+1),...
    parms.cropz:(end-parms.cropz+1)));

  % binarize
  vec_uw(vec_uw<0.5) = 0; vec_uw(vec_uw>0.5) = 1;

  % calculate correlation
  gwcorr(j) = corr(vec,vec_uw);

  fprintf('%s: gwtype %d: corr=%f\n',mfilename,gwtype,gwcorr(j));
end

% choose best fitting gwtype
[maxv,maxi] = max(gwcorr);
if (maxv > parms.min_corr)
  gwtype = gwtypelist(maxi);
  fprintf('%s: best fit gwtype %d (corr = %f)\n',...
    mfilename,gwtype,maxv);
else
  gwtype = [];
  fprintf('%s: WARNING: no fitting gwtype found (max corr = %f for gwtype %d)\n',...
    mfilename,maxv,gwtypelist(maxi));
end;


