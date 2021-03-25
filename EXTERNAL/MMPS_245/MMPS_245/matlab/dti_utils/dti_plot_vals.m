function dti_plot_vals(vol,volmask,plane,xrange,yrange,zrange)
%function dti_plot_vals(vol,[volmask],[plane],[xrange],[yrange],[zrange])
%
% Required:
%   vol: volume values to be plotted (as image)
%
% Optional:
%   volmask: 3D mask volume (if empty, plot for entire volume)
%     {default: []}
%   plane: 1, 2, or 3 to indicate which slice plane (x, y, or z)
%     {default: 3}
%   xrange: range of pixels in x
%   yrange: range of pixels in y
%   zrange: range of pixels in z
%     {default: []  (use all slices)}
%
% Note: Middle slice of range will be shown for plane specified
%
% Created:  07/30/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('volmask','var'), volmask=[]; end;
if ~exist('plane','var') | isempty(plane), plane=3; end;
if ~ismember(plane,[1,2,3])
  error('plane must be 1, 2, or 3');
end;

if ~exist('xrange','var'), xrange=[]; end;
if ~exist('yrange','var'), yrange=[]; end;
if ~exist('zrange','var'), zrange=[]; end;

volsz = size(vol);
if isempty(xrange), xrange = [1:volsz(1)]; end;
if isempty(yrange), yrange = [1:volsz(2)]; end;
if isempty(zrange), zrange = [1:volsz(3)]; end;
if isempty(volmask), volmask = ones(volsz(1:3)); end;

vol = vol(xrange,yrange,zrange);
volmask = volmask(xrange,yrange,zrange);

volsz = size(vol);
if length(volsz)==2
  snum=1;
else
  snum = round(volsz(plane)/2);
end;

switch plane
  case 1
    im = squeeze(vol(snum,:,:));
    mask = squeeze(volmask(snum,:,:));
  case 2
    im = squeeze(vol(:,snum,:));
    mask = squeeze(volmask(:,snum,:));
  case 3
    im = squeeze(vol(:,:,snum));
    mask = squeeze(volmask(:,:,snum));
  otherwise
end;

n1=size(im,1);
n2=size(im,2);
nvox=n1*n2;

im(mask==0)=0;

imagesc(im'); axis image; axis ij;

