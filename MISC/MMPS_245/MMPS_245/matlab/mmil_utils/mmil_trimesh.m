function mmil_trimesh(faces,X,Y,Z,C,crange,cmap);
%function mmil_trimesh(faces,X,Y,Z,C,[crange],[cmap]);
%
%  Last Mod: 09/27/12   by Don Hagler
%

if ~exist('crange','var'), crange=[]; end;

if isempty(crange)
  crange(1) = min(C);
  crange(2) = max(C);
end;

if ~exist('cmap','var'), cmap=[]; end;
if isempty(cmap), cmap = mmil_cmap_redblackblue; end;

C2 = 1+(C-crange(1))*(size(cmap,1)-1)/(crange(end)-crange(1));
C2 = min(size(cmap,1),max(1,C2));
colormap(cmap);
trimesh(faces,X,Y,Z,C2,'CDataMapping','direct');


