function rc_plot_nonretdips(prefix,ampflag,linewidth,fontsize,ylim)
%function rc_plot_nonretdips(prefix,[ampflag],[linewidth],[fontsize],[ylim])
%
% Purpose: plot waveforms for non-retinotopic dipoles included in RCSE fit
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   ampflag: [0|1] whether to plot overall amplitude
%     otherwise, plot normal and tangential components
%     {default = 0}
%   linedith
%     {default = 1.5}
%   fontsize
%     {default = 'Arial'}
%   ylim: if empty, auto-scale y-axis
%     {default = []}
%
% Early Mod: 03/05/08 by Don Hagler
% Last Mod:  08/01/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('ampflag','var') | isempty(ampflag), ampflag = 0; end;
if ~exist('linewidth','var') | isempty(linewidth), linewidth = 1.5; end;
if ~exist('fontsize','var') | isempty(fontsize), fontsize = 12; end;
if ~exist('ylim','var'), ylim = []; end;
if ~exist('fontname','var') | isempty(fontname), fontname = 'Arial'; end;

fprintf('%s: plotting areas for %s\n',mfilename,prefix);

matname=sprintf('matfiles/%s_retsources.mat',prefix);
load(matname);

num_nonret_dips = retmap.num_nonret_dips;

if ~num_nonret_dips
  error('no non-ret dips to plot');
end;

tv = retmap.areas(1).time*1000;

for d=1:num_nonret_dips
  clf;
  if ampflag
    plot(tv,retmap.nonret_dips(d).a,'LineWidth',linewidth);
  else
    hold on;
    plot(tv,retmap.nonret_dips(d).t1,'g','LineWidth',linewidth);
    plot(tv,retmap.nonret_dips(d).t2,'m','LineWidth',linewidth);
    plot(tv,retmap.nonret_dips(d).n,'k','LineWidth',linewidth);
  end;    
  axis(gca,'tight');
  if ~isempty(ylim), set(gca,'YLim',ylim); end;
  set(gca,'FontSize',fontsize)
  set(gca,'FontName',fontname)
  xlabel('Time (msec)','FontSize',fontsize)
  ylabel('Source Amplitude (nA M)','FontSize',fontsize)
  legendstrings = {'tangential 1','tangential 2','normal'};
  legend(legendstrings,'Location','SouthEast','FontSize',fontsize,'FontName',fontname);

  outstem = sprintf('%s-nonretdip-%s-%s',...
    prefix,retmap.nonret_dips(d).name,retmap.nonret_dips(d).hemisphere);;
  print('-dtiff',[outstem '.tif']);
  mmil_printeps(gcf,[outstem '.eps']);
end

