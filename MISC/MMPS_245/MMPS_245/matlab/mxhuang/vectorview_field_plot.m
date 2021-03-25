function vectorview_field_plot(B,chann_index,color_range,view_angle,vv_field_plot);
% function vectorview_field_plot(B,chann_index,color_range,view_angle,vv_field_plot);
% This program plot the MEG field on the Elekta/Neuromag Helmet Surface
% Input variables:
%   B: 306 x 1 MEG field vector, chann 1, 2, 4, 5, ... are gradiometers in (T/m),
%       chann 3, 6, 9,... are magnetometers in (T)
%   chann_index: 306 x 1 with "1" (good channel) and "0" (bad channel)
%   color_range: 1 x 2: [Cmin Cmax] for the range of the color, e.g. 
%       [-300 300]
%   view_angle: the angle to view the field map, same as the one used in "view"
%       e.g., [0 0] is back-view, [90 0] is the R-view, [-90 0] is the L-view, 
%             [180 0] is the front-view, [0 90] is the top-view, etc
%   vv_field_plot: a structure loaded from vectorview_field_plot.mat

id_grad=(1:306)';
id_grad(3:3:306)=[]; % index for gradiometers
id_mag=(3:3:306)'; % index for magnetometers
B(id_grad)=B(id_grad)*1e13; % change from T/m to fT/cm
B(id_mag)=B(id_mag)*1e15; % change from T to fT

id_keep=find(chann_index==1);
B_interp=vv_field_plot.T_306_helmet(:,id_keep)*B(id_keep);
[n,m]=size(vv_field_plot.x_helmet);

fill3(vv_field_plot.x_helmet,vv_field_plot.y_helmet,vv_field_plot.z_helmet,...
    reshape(B_interp(vv_field_plot.geo_helmet(:)),n,m),'FaceColor','interp','EdgeColor','flat');
%axis([-0.15 0.15 -0.15 0.15 -0.15 0.15]), 
axis equal
set(gca,'Clim',color_range), colormap(bluehot), shading interp
xlabel('x')
ylabel('y')
zlabel('z')
view(view_angle)
colorbar('vert')

