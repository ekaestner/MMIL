function ts_plot_vv_fields(B,varargin)
%function ts_plot_vv_fields(B,[options])
%
% Purpose: plot MEG field on the Elekta/Neuromag Helmet Surface
% 
% Required Input:
%   B: 306 x 1 MEG field vector
%     channels 1, 2, 4, 5, ... are gradiometers in (T/m)
%     channels 3, 6, 9,... are magnetometers in (T)
%
% Optional Parameters:
%  'vv_field_plot: structure loaded from vectorview_field_plot.mat
%     if empty, will load vvfpfile
%     {default = []}
%  'vvfpfile': file name for vectorview_field_plot mat file
%     {default = 'vectorview_field_plot.mat'}
%  'chans': 306 x 1 vector of 0s (exclude) and 1s (include)
%     if empty, include all
%     {default = []}
%  'crange': 1 x 2 vector [Cmin Cmax] for the range of color values
%     {default = [-100 100]}
%  'view_angle': angle to view the field map
%       e.g., [0 0] is back-view, [90 0] is the R-view, [-90 0] is the L-view, 
%             [180 0] is the front-view, [0 90] is the top-view, etc
%     {default = [0 20]}
%  'axis_flag': [0|1] draw and label exes
%     {default = 0}
%  'colorbar_flag': [0|1] draw colorbar
%     {default = 0}
%  'colormap': name of colormap
%     {default = 'mmil_cmap_blueblackred'}
%
% Created:  08/20/12 by Don Hagler (based on code from M.X. Huang)
% Last Mod: 09/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'vv_field_plot',[],[],...
  'vvfpfile','vectorview_field_plot.mat',[],...
  'chans',[],[],...
  'crange',[],[],...
  'view_angle',[0 20],[],...
  'axis_flag',false,[false true],...
  'colorbar_flag',false,[false true],...
  'colormap','mmil_cmap_blueblackred',[],...
});

if isempty(parms.vv_field_plot)
  load(parms.vvfpfile);
else
  vv_field_plot = parms.vv_field_plot;
end;

id_grad=(1:306)';
id_grad(3:3:306)=[]; % index for gradiometers
id_mag=(3:3:306)'; % index for magnetometers
B(id_grad)=B(id_grad)*1e13; % change from T/m to fT/cm
B(id_mag)=B(id_mag)*1e15; % change from T to fT

id_keep=find(parms.chans==1);
B_interp=vv_field_plot.T_306_helmet(:,id_keep)*B(id_keep);
[n,m]=size(vv_field_plot.x_helmet);

fill3(vv_field_plot.x_helmet,vv_field_plot.y_helmet,vv_field_plot.z_helmet,...
    reshape(B_interp(vv_field_plot.geo_helmet(:)),n,m),...
    'FaceColor','interp','EdgeColor','flat');

axis equal;
set(gca,'Clim',parms.crange);
colormap(parms.colormap);
shading interp;
if parms.axis_flag
  xlabel('x');
  ylabel('y');
  zlabel('z');
else
  axis off;
end;
view(parms.view_angle)
if parms.colorbar_flag
  colorbar('vert');
end;
