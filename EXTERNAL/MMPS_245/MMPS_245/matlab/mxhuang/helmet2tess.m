function [x,y,z,RgUp,geoUp,rUp] = ...
    helmet2tess(Rsensor,order,spacing,ANG_WIDTH,verbose);
%HELMET2TESS Tesselate based on helmet shaped sensor locations
% function [x,y,z,RgUp,geoUp,rUp,rest,shell,a] = ...
%    hsdig2tess(Rsensor,order,spacing,ANG_WIDTH,verbose);
% Given Rsensor, sensor location data, one three-d Cartesian point per
%              row of Rsensor, in nominal head-coordinate system
%              (x-axis through nasion, y-axis near left prearic, z-axis up);
%       order, the harmonic order to fit (0 = sphere);
%       spacing, the nominal length of one side of the triangles;
% (optional) verbose, to tell you what's happening and plot a picture.
% Output:
%  x,y,z,rUp are suitable for fill3(x,y,z,rUp);  All are 3 x # of triangles,
%   each column of x represents the x-coordinates of the ith triangle,
%   similarly for y and z.  rUp is the distance from the origin to the
%   triangle vertice, makes for nice plots.
%  rtrue is the distance from the translated origin NEW_ORG to the given
%   data in hsdig.
%  rest is the estimated distance to the same data.  plot([rtrue rest])
%   gives an idea of the quality of the fit.
%  a is the coefficients used to fit the spherical harmonics
%  new_org is the new origin, found as the best fitting sphere.
%  nhsdig is the hsdig data, translated by new_org.
%  Try fill3(x,y,z,rUp); 
%      hold on, plot3(nhsdig(:,1),nhsdig(:,2),nhsdig(:,3),'go'), hold off
%   to see how well the triangles fit the 3-d data.



if(exist('verbose') ~= 1),
  verbose = 0;			% silent running
end

X = fmins('sphererr',[0;0;0],[],[],Rsensor); % returns best center
[err,shell] = sphererr(X,Rsensor); % get final fit error and radius
%new_org = X.'; 			% row vector for new origin
%if(verbose),
%  disp(['Best fitting sphere origin:' sprintf(' %.2f',new_org)])
%  disp(sprintf('Best fitting sphere radius: %.2f',shell));
%  disp(sprintf('Arbitrary error in fitting: %f',err));
%end

%Rsensor = Rsensor - new_org(ones(size(Rsensor,1),1),:); % subtract
new_org=[0 0 0]; % using the center of helmet for the center of sphere

rtp = cart2rtp(Rsensor); 		% convert to spherical

rtp(:,2:3) = rtp(:,2:3)*pi/180;	% degrees to radian

t_range = [min(rtp(:,2)) max(rtp(:,2))]; % range of thetas ('elev')
p_range = [min(rtp(:,3)) max(rtp(:,3))]; % range of phi ('azimuth')

% triangles on a sphere, upto max theta
[xs,ys,zs,Rs,geos] = tesselate(shell,spacing,'full'); % get full coverage

% as a function of phi (azimuth), find the largest theta (smallest
%  elevation).
[data_ph,tmp] = sort(rtp(:,3));	% sort the azimuth
data_th = rtp(tmp,2);		% same order for elevation

ndx = find(~diff(data_ph)); 	% where I have repetitions
while(any(ndx)),
  data_ph(ndx) = [];
  % take max of the two and drop the other
  data_th(ndx+1) = max([data_th(ndx) data_th(ndx+1)]')';
  data_th(ndx) = [];
  ndx = find(~diff(data_ph)); 	% where I have repetitions still
end

mean_ph = mean(diff(data_ph));	% average increment in azimuth
win_pts = round(ANG_WIDTH/mean_ph);	% width of window
if(~rem(win_pts,2)),		% it's even
  win_pts = win_pts + 1;	% want odd points
end
win_pts2 = (win_pts-1)/2;	% one side of window
tot_pts = length(data_ph);	% number of data points
data_slide = hankel(data_th([[-win_pts2:0]+tot_pts [1:win_pts2]]),...
    data_th([[win_pts2:tot_pts] [1:(win_pts2-1)]]));
%mx_th = max(data_slide).';	% gives maximum theta as func of phi
mx_th = max(data_slide).'+8*pi/180;	% gives maximum theta as func of phi
% but now have to find where I have multiple identical phi's

xyz_sphere = [xs(:) ys(:) zs(:)]; % all points on the tesselated sphere
rtp_sphere = cart2rtp(xyz_sphere); % in spherical
triangle_s = size(xs,2);	% number of triangles
rs = reshape(rtp_sphere(:,1),3,triangle_s);
ts = reshape(rtp_sphere(:,2),3,triangle_s)/180*pi; % and convert to rads
ps = reshape(rtp_sphere(:,3),3,triangle_s)/180*pi;

% interpolate the corresponding theta val, pad front and back for wrap
n_ph = length(data_ph);
data_ph = [data_ph(n_ph)-2*pi;data_ph;2*pi+data_ph(1)]; % wrap both ends
mx_th = [mx_th(n_ph);mx_th;mx_th(1)];
tsi = interp1(data_ph,mx_th,min(ps)).'; 
trim = find(max(ts) <= tsi'); 	% values in bounds, JCM change min to max 4/4
xUp = xs(:,trim);		% keep only the upper regions
yUp = ys(:,trim);
zUp = zs(:,trim);
geoUp = geos(:,trim);		% so we have geoUp, geos

while(0)			% disable
% find those vertices on the interface between the two regions
gfull = sort(geos(:));		% all of the vertices
gfull = diff([0;find(diff([gfull;0]))]); % number of connections per vertice
% we don't have all vertices here, so have to remap
tmp = sort(geoBot(:));
gBoti = find(diff([tmp;0]));	% index of where changes occured
gBot = zeros(size(gfull));	% default is no connection
gBot(tmp(gBoti)) = diff([0;gBoti]);	% number of connections
% now find where we have non-zero connections in gBot
ndx = find(gBot);
ndx_edge = ndx(find(gBot(ndx) ~= gfull(ndx))); % number of vertices has changed
% so now we have the vertices on the edge between the two regions
% Later, we will adjust these values to the harmonic fit
end				% disable

xyzUp = [xUp(:) yUp(:) zUp(:)]; % all points on the tesselated sphere
rtpUp = cart2rtp(xyzUp); % in spherical
thetaUp = rtpUp(:,2)*pi/180;	% the thetas
phiUp = rtpUp(:,3)*pi/180;	% the phis

cols = sum(1:2:(2*order+1));	% number of terms given the order
PUp = zeros(length(thetaUp),cols); 	% polynomials for interpolation
Prtp = zeros(size(rtp,1),cols);	% polynomials of data

ndx = 0;
for l = 0:order,
  for m = -l:l,
    ndx = ndx + 1;
    PUp(:,ndx) = spherharm(thetaUp,phiUp,l,m); % the upper region
    Prtp(:,ndx) = spherharm(rtp(:,2),rtp(:,3),l,m); % the original data
  end
end

rtrue = rtp(:,1);		% true distances

a = Prtp\rtrue; 		% fit of distances

rest = abs(Prtp*a); 		% estimated distances

rUp = abs(PUp*a); 		% refit to interpolations

% all of the interpolated points
xyzUp = rtp2cart([rUp [thetaUp phiUp]*180/pi]); 

% reshape suitable for fill3
triUp = size(xyzUp,1)/3;	% the number of triangles in upper

x = reshape(xyzUp(:,1),3,triUp);
y = reshape(xyzUp(:,2),3,triUp);
z = reshape(xyzUp(:,3),3,triUp);
rUp = reshape(rUp,3,triUp);	% suitable for fill3 color

if(verbose)
  disp(sprintf('Generated %.0f triangles interpolated through your data',...
      triUp));
  figure(windfind('Tesselation Results'));
  windclf
  % give almost black edges (so they don't invert on print)
  colormap(hot(64));
  cmap = colormap;
  % set edge color to halfway.
  % ht=fill3(x,y,z,rUp,'face','interp','edge',cmap(size(cmap,1)/2,:));
  ht=fill3(x,y,z,rUp,'FaceColor','interp','EdgeColor','flat');
  cax = caxis;
  caxis = [cax(1)*.9 cax(2)*1.1];	% not quite from black to hot
  hold on
  hp=plot3(Rsensor(:,1),Rsensor(:,2),Rsensor(:,3),'go');
  hold off
  title(sprintf(...
      'Tesselation Results, order %.0f, spacing of %.2f on a %.2f sphere',...
      order,spacing,shell));
  xlabel('X-axis')
  ylabel('Y-axis')
  zlabel('Z-axis')
  axis('square')
  axis('equal')
  view([135 5])			% somewhat off to the subject's right side
  axis_cube;
  hh = gca;
  hhc = colorbar('v');
  axes(hhc);
  ylabel('Distance from centroid to surface')
  axes(hh);			% return to original axis
end

RgUp = zeros(max(geoUp(:)),3) + NaN; % initialize
RgUp(geoUp(:),:) = [x(:) y(:) z(:)]; % load valid indices

% now reconvert everything back to the original coordinates
x = x+new_org(1);
y = y+new_org(2);
z = z+new_org(3);
RgUp = RgUp + ones(size(RgUp,1),1)*new_org;

return
