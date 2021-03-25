function rc_plot_RCSE_summary(prefix)
%function rc_plot_RCSE_summary(prefix)
%
% Purpose: plot RCSE source waveforms, data and err variance
%
% Required Input:
%   prefix: RCSE prefix
%
% Early Mod: 10/06/07 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

fprintf('%s: plotting inverse summary for %s\n',mfilename,prefix);

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matname=sprintf('matfiles/%s_retsources.mat',parms.prefix);
load(matname);
matname=sprintf('matfiles/%s_results.mat',parms.prefix);
load(matname);

num_locs = retmap.num_locs;
[num_tpoints num_measurements] = size(Y);
num_channels = num_measurements/num_locs;
[num_tpoints num_sources] = size(S);
num_areas = retmap.num_areas;
num_locs = retmap.num_locs;

indy_locs_flag = 0;
loose_flag = 0;
if exist('parms','var')
  if isfield(parms,'indy_locs_flag')
    indy_locs_flag = parms.indy_locs_flag;
  end;
  if isfield(parms,'loose_flag')
    loose_flag = parms.loose_flag;
  end;
end;

if ~exist('channel','var')
  channel = 0;
end;
if(channel > 0)
  indices = channel:num_channels:num_measurements;
else
  indices = 1:num_measurements;
end

%figure

titles = {'Y' 'Yfit' 'err' 'var Y' 'var err' 'norm var err' ...
          'S(V1-V2)' 'S(other visual areas)' 'S(extra dips)'};

for i=1:9
  h(i) = subplot(3,3,i);
end

tv = retmap.areas(1).time*1000;

plot(h(1),tv,Y(:,indices));
plot(h(2),tv,Yfit(:,indices));
plot(h(3),tv,E(:,indices));
%% todo: average var_Y across channel types
plot(h(4),tv,var_Ygrad);
plot(h(5),tv,var_Egrad);
plot(h(6),tv,norm_var_E);

if loose_flag
  if(num_areas>=3)
    plot(h(7),tv,S(:,1:2*num_locs*3));
    plot(h(8),tv,S(:,2*num_locs*3+1:num_areas*num_locs*3));
  else
    plot(h(7),tv,S(:,1:num_areas*num_locs*3));
  end
  if(num_sources>num_areas*num_locs*3)
    plot(h(9),tv,S(:,num_areas*num_locs+1:num_sources));
  end
elseif indy_locs_flag
  if(num_areas>=3)
    plot(h(7),tv,S(:,1:2*num_locs));
    plot(h(8),tv,S(:,2*num_locs+1:num_areas*num_locs));
  else
    plot(h(7),tv,S(:,1:num_areas*num_locs));
  end
  if(num_sources>num_areas*num_locs)
    plot(h(9),tv,S(:,num_areas*num_locs+1:num_sources));
  end
else
  if(num_areas>=3)
    plot(h(7),tv,S(:,1:2));
    plot(h(8),tv,S(:,3:num_areas));
  else
    plot(h(7),tv,S(:,1:num_areas));
  end
  if(num_sources>num_areas)
    plot(h(9),tv,S(:,num_areas+1:num_sources));
  end
end;

for i=1:9
  title(h(i),titles{i});
  axis(h(i),'tight');
end
set(h(6),'YLim',[0 1]); % normalized residual variance

print('-dtiff',sprintf('%s-inverse_summary.tif',prefix));
