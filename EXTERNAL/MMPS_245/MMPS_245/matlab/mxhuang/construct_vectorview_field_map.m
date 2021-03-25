% This script construct the transformation for MEG field map
% of the Elekta/Neuromag VectorView system
% 
% (c) Mingxiong Huang, Ph.D. April 2006

% load the sensor configuration and accurate integration points
load sensor306_accurate_default.mat

id_grad=(1:306)';
id_grad(3:3:306)=[];
id_mag=(3:3:306)';

% Fit the sensor helmet with 6th order spherical harmonics and 
% tessellate the surface with 8mm size triangles 
[x,y,z,RgUp,geoUp,rUp] = helmet2tess(intpnt_loc,6,0.012,90*pi/180,1);

% find those bad vertices that are too far away (>2.5cm) from the sensors
id_bad=[]; % index of bad vertices 
for i=1:size(RgUp,1)
    temp=rownorm(ones(3264,1)*RgUp(i,:)-intpnt_loc);
    if min(temp)>=0.025, id_bad=[id_bad;i]; end
end
R_helmet=RgUp;
R_helmet(id_bad,:)=[]; % x y z of the good vertices on the helmet


id_bad_geo=[]; % index of bad triangles
for i=1:length(id_bad),
    [temp1,temp2]=find(geoUp==id_bad(i));
    id_bad_geo=[id_bad_geo;temp2];
end
x_helmet=x;
x_helmet(:,id_bad_geo)=[];
y_helmet=y;
y_helmet(:,id_bad_geo)=[];
z_helmet=z;
z_helmet(:,id_bad_geo)=[];
r_helmet=rUp;
r_helmet(:,id_bad_geo)=[];

% reorganize the geo matrix of the triangles on the helmet
geo_helmet=geoUp;
geo_helmet(:,id_bad_geo)=[];
id_shift=zeros(3*size(geo_helmet,2),1);
for i=1:length(id_bad)
    id_temp=find(geo_helmet(:)>id_bad(i));
    if isempty(id_temp)==0,
        id_shift(id_temp)=id_shift(id_temp)+1;
    end
end
geo_helmet(:)=geo_helmet(:)-id_shift;

% now estimate the orientations of the helmet 
[side_hmt,centroid_hmt,ncentroid_hmt] =tri_analysis(x_helmet,y_helmet,z_helmet,R_helmet,geo_helmet);

for i=1:size(R_helmet,1)
    temp=rownorm(ones(size(centroid_hmt,1),1)*R_helmet(i,:)-centroid_hmt);
    [temp2,id_temp]=sort(temp);
    O_helmet(i,:)=mean(ncentroid_hmt(id_temp(1:3),:)); % average across the most close 3
end
O_helmet=O_helmet./(rownorm(O_helmet)*ones(1,3)); % normalize the orientation

% now create 7 cm dupole grid
source_grid_struct=sphere_tri('ico',4,0.07);

% prepare for the gain matrix with vector view sensors
P.sensor=intpnt_loc';
P.orient=intpnt_ori';
P.weight=0.1;
G=sarvas5(source_grid_struct.vertices',P,-1);
G=gain_intpnt2chan(Coil(1:306),G); % convert the integration points into 306 channels

% remove the bias between the gradiometer and magnetometer
G_norm_grad=mean(rownorm(G(id_grad,:)));
G_norm_mag=mean(rownorm(G(id_mag,:)));
G_norm_fact=ones(306,1);
G_norm_fact(id_grad)=G_norm_grad;
G_norm_fact(id_mag)=G_norm_mag;

% now the pinv of the gain matrix
[V,S,U]=svd((G./(G_norm_fact*ones(1,size(G,2))))',0); % SVD after removing the bias
S_inv=zeros(306,306);
for i=1:80, % pick up the largest 80 singular values
    S_inv(i,i)=1.0/S(i,i);
end
Ginv=V*S_inv*U';
Ginv=Ginv./(ones(size(Ginv,1),1)*G_norm_fact'); % rescale back, note ./ not .* for inverse

% now the second gain matrix for the tess helmet
P_helmet.sensor=R_helmet';
P_helmet.orient=O_helmet';
P_helmet.weight=0.1;
G_helmet=sarvas5(source_grid_struct.vertices',P_helmet,-1);

% finally the transformation matrix from 306 to tess-helmet 
T_306_helmet=G_helmet*Ginv;

% form a structure
vv_field_plot.x_helmet=x_helmet;
vv_field_plot.y_helmet=y_helmet;
vv_field_plot.z_helmet=z_helmet;
vv_field_plot.r_helmet=r_helmet;
vv_field_plot.R_helmet=R_helmet;
vv_field_plot.geo_helmet=geo_helmet;
vv_field_plot.T_306_helmet=T_306_helmet;

%save vectorview_field_plot.mat vv_field_plot

