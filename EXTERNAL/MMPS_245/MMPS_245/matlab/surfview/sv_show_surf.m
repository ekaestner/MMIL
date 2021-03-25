function h = sv_show_surf(surf,varargin)
%function h = sv_show_surf(surf,[options])
%
% Purpose: display surface mesh with color data
%
% Required Input:
%   surf: surface mesh struct (with fields vertices and faces)
%
% Optional Input:
%   'cvals': matrix of rgb color values (0 - 1) with size [nverts,3]
%     if empty, display surface with uniform gray
%     see sv_linear_cvals, sv_category_cvals, etc.
%     {default = []}
%   'view': string defining view angle
%     allowed values: 'left','right','pos','ant','sup','ven','infl','infr'
%     {default = 'left'}
%   'zoom': zoom factor
%     {default = 1}
%
% Output:
%   h: figure handle
%
% Created:  09/21/12 by Don Hagler (based on code by Anders Dale)
% Last Mod: 08/14/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
h = [];

parms = mmil_args2parms(varargin,{,...
  'cvals',[],[],...
  'view','left',{'left','right','pos','ant','sup','ven','infl','infr'},...
  'zoom',1,[0.01,100],...
});
%% todo: additional options, e.g. AmbientStrength, lighting, etc.

nverts = length(surf.vertices);
if isempty(parms.cvals)
  parms.cvals = 0.5 * ones(nverts,3);
else
  if size(parms.cvals,1) ~= nverts || size(parms.cvals,2) ~= 3
    error('incorrect size of cvals (must be [%d,3])',nverts);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = patch(surf);
set(h,'FaceColor','interp','EdgeColor','none','LineStyle','none',...
      'SpecularStrength',0,'SpecularStrength',0,'AmbientStrength',0.4,...
      'FaceVertexCData',parms.cvals);

%      'DiffuseStrength',0,'SpecularStrength',0,'AmbientStrength',0.4,...

set(gcf,'Color',[0 0 0]);
set(gca,'Color',[0 0 0]);

switch lower(parms.view)
  case 'left',     camorbit(-90,0); camorbit(0,-90);
  case 'right',    camorbit(90,0);  camorbit(0,-90);
  case 'ant',      camorbit(180,0); camorbit(0,-90);
  case 'pos',      camorbit(0,-90);
  case 'sup',      camorbit(0,0); 
  case 'infl'
                   camorbit(180,0); camorbit(0,180);
                   camorbit(90,0);
  case 'infr'
                   camorbit(180,0); camorbit(0,180);
                   camorbit(-90,0);
end

camlight headlight;
daspect([1 1 1]);
material dull;
lighting phong;

if parms.zoom~=1, camzoom(parms.zoom); end;

