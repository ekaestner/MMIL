function rc_plot_loose_areas(prefix)
%function rc_plot_loose_areas(prefix)
%
% Purpose: plot RCSE sources waveforms for each area
%   assumes loose orientation constraint
%   plots normal and tangential components
%
% Required Input:
%   prefix: RCSE prefix
%
% Early Mod: 04/29/09 by Don Hagler
% Last Mod:  04/29/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

fprintf('%s: plotting areas for %s\n',mfilename,prefix);

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matname=sprintf('matfiles/%s_results.mat',parms.prefix);
load(matname);

[num_tpoints num_sources num_contrasts] = size(results.S);
num_areas = results.retmap.num_areas;
num_locs = results.retmap.num_locs;

loose_flag = 0;
if exist('parms','var')
  if isfield(parms,'loose_flag')
    loose_flag = parms.loose_flag;
  end;
end;

tv = results.retmap.areas(1).time*1000;

nrows = num_areas;
ncols = 2;

hold on;
p = 1;
for i=1:num_areas
  for c=1:num_contrasts
    if loose_flag
      j = 1+(i-1)*num_locs*3;
      k = j+num_locs*3-1;

      h=subplot(nrows,ncols,p); p=p+1;
      plot(h,tv,results.S(:,j:3:k,c));
      axis(h,'tight');
      tmp_title = sprintf('%s normal',results.retmap.areas(i).name);
      title(h,tmp_title);

      h=subplot(nrows,ncols,p); p=p+1;
      plot(h,tv,results.S(:,[j+1:3:k+1,j+2:3:k+2],c));
      tmp_title = sprintf('%s tangential',results.retmap.areas(i).name);
      title(h,tmp_title);
      axis(h,'tight');
      set(h,'XGrid','on','XTick',[-100 0 100 200 300 400])
    else
      j = 1+(i-1)*num_locs;
      k = j+num_locs-1;
      h=subplot(nrows,ncols,p); p=p+2;
      plot(h,tv,results.S(:,j:k,c));
      axis(h,'tight');
      set(h,'XGrid','on','XTick',[-100 0 100 200 300 400])
      tmp_title = sprintf('%s normal',results.retmap.areas(i).name);
      title(h,tmp_title);
    end;
  end;
end

print('-dtiff',sprintf('%s-loose-areas.tif',prefix));
