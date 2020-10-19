function rc_plot_areas_legend(prefix)
%function rc_plot_areas_legend(prefix)
%
% Purpose: plot waveforms for each visual area
%   and include legend
%   particularly for independent locations (indy_locs_flag)
%
% Required Inout:
%   prefix: RCSE prefix
%
% Early Mod: 03/08/08 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: use varargin to add more options?
if ~mmil_check_nargs(nargin,1), return; end;

fprintf('%s: plotting areas for %s\n',mfilename,prefix);

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matname=sprintf('matfiles/%s_retsources.mat',parms.prefix);
load(matname);
matname=sprintf('matfiles/%s_results.mat',parms.prefix);
load(matname);

[num_tpoints num_sources] = size(S);
num_areas = retmap.num_areas;
num_locs = retmap.num_locs;

locstrings=[];
for loc=1:floor(num_locs/2)
  locstrings{loc} = num2str(loc);
end

loose_flag = 0;
if exist('parms','var')
  if isfield(parms,'loose_flag')
    loose_flag = parms.loose_flag;
  end;
  if isfield(parms,'indy_locs_flag')
    indy_locs_flag = parms.indy_locs_flag;
  end;
end;

tv = retmap.areas(1).time*1000;

nrows = num_areas;
ncols = 2;


clf;
p = 1;
for i=1:num_areas
  if loose_flag
    h=subplot(nrows,ncols,p); p=p+1;
    j = 1+(i-1)*num_locs*3;
    k = j+(num_locs/2)*3-1;
    plot(h,tv,S(:,j:3:k));
    axis tight;
    tmp_title = sprintf('%s lh',retmap.areas(i).name);
    title(h,tmp_title);
    legend(locstrings,'Location','EastOutside');

    h=subplot(nrows,ncols,p); p=p+1;
    j = k+1;
    k = j+(num_locs/2)*3-1;
    plot(h,tv,S(:,j:3:k));
    axis tight;
    set(h,'XGrid','on','XTick',[-100 0 100 200 300 400])
    tmp_title = sprintf('%s rh',retmap.areas(i).name);
    title(h,tmp_title);
    legend(locstrings,'Location','EastOutside');
  elseif indy_locs_flag
    h=subplot(nrows,ncols,p); p=p+1;
    j = 1+(i-1)*num_locs;
    k = j+(num_locs/2)-1;
    plot(h,tv,S(:,j:k));
    axis tight;
    set(h,'XGrid','on','XTick',[-100 0 100 200 300 400])
    tmp_title = sprintf('%s lh',retmap.areas(i).name);
    title(h,tmp_title);
    legend(locstrings,'Location','EastOutside');

    h=subplot(nrows,ncols,p); p=p+1;
    j = k+1;
    k = j+(num_locs/2)-1;
    plot(h,tv,S(:,j:k));
    set(h,'XGrid','on','XTick',[-100 0 100 200 300 400])
    axis tight;
    tmp_title = sprintf('%s rh',retmap.areas(i).name);
    title(h,tmp_title);
    legend(locstrings,'Location','EastOutside');
  else
    h=subplot(nrows,ncols,p); p=p+1;
    plot(h,tv,S(:,i));
    axis tight;
    set(h,'XGrid','on','XTick',[-100 0 100 200 300 400])
    tmp_title = sprintf('%s',retmap.areas(i).name);
    title(h,tmp_title);

    p=p+1;
  end;
end

print('-dtiff',sprintf('%s-areas-legend.tif',prefix));
