function rc_plot_thetas_normerr(prefix)
%function rc_plot_thetas_normerr(prefix)
%
% Purpose: plotting normalized residual error for
%   each stimulus location
%
% Required Input:
%   prefix: RCSE prefix
%
% Early Mod: 06/21/10 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: generalize for any stimulus locations, not just ring
%% todo: merge with rc_plot_thetas_err

if ~mmil_check_nargs(nargin,1), return; end;
radius = 0.8;
width = 0.1;
height = 0.1;
scale_max = 10;

fprintf('%s: plotting error for %s\n',mfilename,prefix);

matname=sprintf('matfiles/%s_retsources.mat',prefix);
load(matname);
matname=sprintf('matfiles/%s_results.mat',prefix);
load(matname);

[num_tpoints num_sources] = size(S);
[num_tpoints num_measurements] = size(err);
num_locs = retmap.num_locs;
cond_order = retmap.cond_order;
num_sensors = num_measurements/num_locs;
orig_num_locs = retmap.orig_num_locs;
tv = retmap.areas(1).time_vector*1000;
A_offset = 360/(orig_num_locs*2);

hold on;
clf;
for theta=1:num_locs
  cond = cond_order(theta);
  A = A_offset + (cond-1)*360/orig_num_locs;
  px = 0.5+radius*cos(A*pi/180)/2-width/2;
  py = 0.5+radius*sin(A*pi/180)/2-height/2;  
  j = 1+(theta-1)*num_sensors;
  k = j+num_sensors-1;
  h=subplot('position',[px,py,width,height]);
  stdvY = std(Y,0,2);
  plot(h,tv,err(:,j:k)./(stdvY*ones(1,num_sensors)));
  axis(h,[-100 300 -scale_max scale_max]);
  axis(h,'manual');
  axis off;
  title(h,sprintf('loc %d',cond));
end

hold off;

print('-dtiff',sprintf('%s-thetas_normerr.tif',prefix));
