function rc_plot_retdips(prefix)
%function rc_plot_retdips(prefix)
%
% Purpose: plot waveforms for extra retinotopy dipoles
%
% Required Input:
%   prefix: RCSE prefix
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

colr = ['b' 'g' 'r' 'c' 'm' 'y' 'k'];
colr = [colr colr colr];

fprintf('%s: plotting areas for %s\n',mfilename,prefix);

matname=sprintf('matfiles/%s_retsources.mat',prefix);
load(matname);

num_ret_dips = retmap.num_ret_dips;

if ~num_ret_dips
  error('no ret dips to plot');
end;

tv = retmap.areas(1).time*1000;

nrows = floor(sqrt(num_ret_dips));
ncols = ceil(num_ret_dips/nrows);

clf;
for d=1:num_ret_dips
  h=subplot(nrows,ncols,d);
  hold on
  for q=1:retmap.num_ret_dip_locs
    plot(h,tv,retmap.ret_dips(d).a(:,q),colr(q));
  end;
  axis(h,'tight');
  name = sprintf('%s-%s',retmap.ret_dips(d).name,retmap.ret_dips(d).hemisphere);
  title(h,name);
end

print('-dtiff',sprintf('%s-retdips.tif',prefix));
