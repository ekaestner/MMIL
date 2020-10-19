function dti_plot_quiv(vol,M,volV0,varargin)
%function dti_plot_quiv(vol,M,volV0,[options])
%
% Required:
%   vol: 3D volume displayed as underlay
%   M  : 4x4 vox2ras transform
%   volV0: 4D volume containing vector for each voxel (nx,ny,nz,3)
%   vol,volV0,volmask should be same orientation and resolution 
%
% Optional:
%   volmask: 3D volume of mask image (if empty, plot for entire volume)
%     {default: []}
%   plane: 1, 2, or 3 to indicate which slice plane (x, y, or z)
%     {default: 3}
%   xrange: range of pixels in x
%   yrange: range of pixels in y
%   zrange: range of pixels in z
%     {default: []  (use all slices)}
%   linestyle: linestyle for quiver plot
%     {default: 'r.'}
%   linewidth: width in points of lines
%     {default: 1.5}
%   crange: range of values for scaling vol
%     {default: []}
%   
% Notes: Middle slice of range will be shown for plane specified
%   All volumes are reoriented to LPI
%     and diffusion directions are reoriented assuming they are RAS
%
% Created:  07/30/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'volmask',[],[],...
  'plane',3,[1 2 3],...
  'xrange',[],[],...
  'yrange',[],[],...
  'zrange',[],[],...
  'linestyle','r.',[],...
  'linewidth',1.5,[],...
  'crange',[],[],...
});  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(vol)
  [vol,M_vol]= fs_reorient(vol,M,'LPI');
end;

if ~isempty(volV0)
  [volV0,M_V0]= fs_reorient(volV0,M,'LPI');
end;

if ~isempty(parms.volmask)
  [volmask,M_mask]= fs_reorient(parms.volmask,M,'LPI');
end;

if ~exist('volmask','var'), volmask=[]; end;

volsz = size(vol);
nx = volsz(1);
ny = volsz(2);
nz = volsz(3);
if isempty(parms.xrange), parms.xrange = [1:nx]; end;
if isempty(parms.yrange), parms.yrange = [1:ny]; end;
if isempty(parms.zrange), parms.zrange = [1:nz]; end;
if isempty(volmask), volmask = ones(nx,ny,nz); end;

vol = vol(parms.xrange,parms.yrange,parms.zrange);
volmask = volmask(parms.xrange,parms.yrange,parms.zrange);
volV0 = volV0(parms.xrange,parms.yrange,parms.zrange,:);

volsz = size(vol);
if length(volsz)==2
  snum=1;
else
  snum = round(volsz(parms.plane)/2);
end;

% make image
switch parms.plane
  case 1
    im = squeeze(vol(snum,:,:));
    mask = squeeze(volmask(snum,:,:));
    V0 = -squeeze(volV0(snum,:,:,:));
  case 2
    im = squeeze(vol(:,snum,:));
    mask = squeeze(volmask(:,snum,:));
    V0 = -squeeze(volV0(:,snum,:,:));
  case 3
    im = squeeze(vol(:,:,snum));
    mask = squeeze(volmask(:,:,snum));
    V0 = -squeeze(volV0(:,:,snum,:));
  otherwise
end;

n1=size(im,1);
n2=size(im,2);
nvox=n1*n2;

imask = find(mask>0);
vecV0 = reshape(V0,[nvox,3]);
[quiv_x,quiv_y] = ind2sub([n1,n2],imask);

switch parms.plane
  case 1
    quiv_dx = vecV0(imask,2);
    quiv_dy = vecV0(imask,3);
  case 2
    quiv_dx = vecV0(imask,1);
    quiv_dy = vecV0(imask,3);
  case 3
    quiv_dx = vecV0(imask,1);
    quiv_dy = vecV0(imask,2);
end;

hold on;
if ~isempty(parms.crange)
  imagesc(im',parms.crange);
else
  imagesc(im');
end;
axis image; axis ij;
if ~isempty(imask)
  h=quiver(quiv_x-quiv_dx/2,quiv_y-quiv_dy/2,quiv_dx,quiv_dy,0,parms.linestyle);
  set(h,'LineWidth',parms.linewidth);
end;
hold off;
