function dti_plot_cmap(volb0,volFA,volV0,plane,xrange,yrange,zrange)
%function dti_plot_cmap(volb0,volFA,volV0,[plane],[xrange],[yrange],[zrange])
%
% Required:
%   volb0: volume containing b=0 image
%   volFA: volume containing FA
%   volV0: 4D volume containing principal eigen vector
%
% Optional:
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
% Last Mod: 04/08/15 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

if ~exist('plane','var'), plane=[]; end;
if isempty(plane), plane=3; end;
if ~ismember(plane,[1,2,3])
  error('plane must be 1, 2, or 3');
end;

if ~exist('xrange','var'), xrange=[]; end;
if ~exist('yrange','var'), yrange=[]; end;
if ~exist('zrange','var'), zrange=[]; end;

%thresh_b0 = 120;
thresh_b0 = eps;

volsz = size(volb0);
if isempty(xrange), xrange = [1:volsz(1)]; end;
if isempty(yrange), yrange = [1:volsz(2)]; end;
if isempty(zrange), zrange = [1:volsz(3)]; end;

volb0 = volb0(xrange,yrange,zrange);
volFA = volFA(xrange,yrange,zrange);
volV0 = volV0(xrange,yrange,zrange,:);

volsz = size(volb0);
if length(volsz)==2
  snum=1;
else
  snum = round(volsz(plane)/2);
end;

% make colormap
switch plane
  case 1
    b0 = squeeze(volb0(snum,:,:));
    FA = squeeze(volFA(snum,:,:));
    V0 = squeeze(volV0(snum,:,:,:));
  case 2
    b0 = squeeze(volb0(:,snum,:));
    FA = squeeze(volFA(:,snum,:));
    V0 = squeeze(volV0(:,snum,:,:));
  case 3
    b0 = squeeze(volb0(:,:,snum));
    FA = squeeze(volFA(:,:,snum));
    V0 = squeeze(volV0(:,:,snum,:));
  otherwise
end;

FA(b0<thresh_b0)=0;

% make colormap
cmap(:,:,1)=(abs(squeeze(V0(:,:,1))).*FA)';
cmap(:,:,2)=(abs(squeeze(V0(:,:,2))).*FA)';
cmap(:,:,3)=(abs(squeeze(V0(:,:,3))).*FA)';

imshow(cmap); axis image; axis ij;

