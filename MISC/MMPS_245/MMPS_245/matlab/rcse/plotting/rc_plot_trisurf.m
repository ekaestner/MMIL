function rc_plot_trisurf(trifile)
%function rc_plot_trisurf(trifile)
%
% Required Input:
%   trifile: full path of freesurfer tri format surface file
%
% Early Mod: 01/08/10 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

colorval = 10;
face_alpha = 0.6;
edge_alpha = 0.1;

surf = fs_read_trisurf(surffile);
% plot surf
cla;
hold on;
colors = colorval*ones(surf.nverts,1);
h=trisurf(surf.faces,...
          surf.vertices(:,1),...
          surf.vertices(:,2),...
          surf.vertices(:,3),...
          colors);
colormap hsv;
caxis([0 100]);

set(h,'EdgeAlpha',edge_alpha,'FaceAlpha',face_alpha);

%axis image;
%axis off;
xlabel('x');
ylabel('y');
zlabel('z');

return;

