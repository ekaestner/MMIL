function bem_mesh_plot(vertices,faces)

% function bem_mesh_plot(vertices,faces)
% plot bem meshes
% vertices: cell array containing 3 or 1 meshes
% faces: cell array containing 3 or 1 geo matrix
%
% (c) M. Huang, Ph.D.

color_edge=['r','g','b'];

clf
hold on
for i=1:length(vertices)
    mesh=vertices{i};
    geo=faces{i};
    x=reshape(mesh(geo,1),size(geo,1),3)';
    y=reshape(mesh(geo,2),size(geo,1),3)';
    z=reshape(mesh(geo,3),size(geo,1),3)';
    fill3(x,y,z,ones(size(z)),'edgecolor',color_edge(i),'facecolor','none');
end
rotate3d
