function sum_solid_angle=check_tri_order(mesh,geo,dist);
% function sum_solid_angle=check_tri_order(mesh,geo,dist);
% check to see if the trigular mesh is in the right order
% mesh: n_node x 6 [x y z nx ny nz]
% geo: geometry matrix n_tri x 3
% dist: 1 x 1 the distance from inside observation
% sum_solid_angle: n_node by 1, the total sum of solid angle from each
% observation point, should be 4*pi or -4*pi

N_node=size(mesh,1);
r_node=mesh(:,1:3);
r_normal=mesh(:,4:6);
r_node_inside=r_node-r_normal*dist; % observation is a few mm inside each the vertex

r1=mesh(geo(:,1),1:3);
r2=mesh(geo(:,2),1:3);
r3=mesh(geo(:,3),1:3);

sum_solid_angle=zeros(N_node,1);
for i=1:N_node
    sum_solid_angle(i)=sum(solid_angle2(r_node_inside(i,:),r1,r2,r3));
end

