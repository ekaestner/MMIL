function mmil_imagesc(X,varargin)
%function output = mmil_imagesc(X,[options])
%
% Purpose: run imagesc to make plot of image, then save to file
%
% Required Input:
%   X: 2-dimensional matrix
%
% Optional Parameters:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem (without extension)
%     If not supplied, image not saved
%     {default = []}
%   'tif_flag': [0|1] save plot as tif
%     {default = 1}
%   'eps_flag': [0|1] save plot as eps
%     {default = 0}
%   'visible_flag': [0|1] display image on screen
%     {default = 1}
%   'cbar_flag': [0|1] display color bar
%     {default = 0}
%   'clim': image scaling range [low high]
%     if not supplied, use full colormap
%     {default = []}
%   'cmap': color map string
%     {default = 'jet'}
%   'title': string displayed above image
%     {default = []}
%   'fontname'
%     {default = 'Arial'}
%   'fontsize'
%     {default = 12}
%
% Created:  03/14/12 by Don Hagler
% Last Mod: 08/01/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',pwd,[],...
  'outstem',[],[],...
  'tif_flag',false,[false true],...
  'eps_flag',false,[false true],...
  'visible_flag',true,[false true],...
  'cbar_flag',false,[false true],...
  'clim',[],[],...
  'cmap','jet',[],...
  'title',[],[],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
...
  'offset',0,[],...
});

if isempty(parms.outstem)
  parms.tif_flag = 0;
  parms.eps_flag = 0;
else
  if mmil_isrelative(parms.outstem)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
end;

clf;
if ~parms.visible_flag
  set(gcf,'Visible','off');
end;

if parms.offset
  parms.clim(2) = parms.clim(2) + parms.offset;
end;

if isempty(parms.clim)
  imagesc(X);
else
  imagesc(X,parms.clim);
end;

axis off
colormap(parms.cmap);

if parms.cbar_flag
  h_cb = colorbar;
  if parms.offset
    YTick = get(h_cb,'YTick');
    YTickLabel = YTick - parms.offset;
    set(h_cb,'YTick',YTick(YTickLabel>=0));
    set(h_cb,'YTickLabel',YTickLabel(YTickLabel>=0));
  end;
end;

set(gca,'FontSize',parms.fontsize,'FontName',parms.fontname)
if ~isempty(parms.title)
  title(parms.title,'FontSize',parms.fontsize,'FontName',parms.fontname);
end;

if parms.tif_flag
  print('-dtiff',sprintf('%s.tif',parms.outstem));
end;
if parms.eps_flag
  mmil_printeps(gcf,parms.outstem);
end;
