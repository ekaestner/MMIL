function fname = rc_plot_areas(prefix,areas,linewidth,fontsize,loc,ylim,linestyle)
%function fname = rc_plot_areas(prefix,[areas],[linewidth],[fontsize],[loc],[ylim],[linestyle])
%
% Purpose: plot waveforms for each visual area
%
% Early Mod: 05/27/09 by Don Hagler
% Last Mod:  08/01/13 by Don Hagler
%

%% todo: use varargin

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('areas','var'), areas = []; end;
if ~exist('linewidth','var') || isempty(linewidth), linewidth = 1.5; end;
if ~exist('fontsize','var') || isempty(fontsize), fontsize = 12; end;
if ~exist('loc','var') || isempty(loc), loc = []; end;
if ~exist('ylim','var'), ylim = []; end;
if ~exist('linestyle','var') || isempty(linestyle), linestyle = []; end;
if ~exist('fontname','var') || isempty(fontname), fontname = 'Arial'; end;

fprintf('%s: plotting areas for %s\n',mfilename,prefix);

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matname=sprintf('matfiles/%s_results.mat',prefix);
load(matname);

if ~isfield(parms,'ncontrasts')
  parms.ncontrasts = 1;
end;

[num_tpoints num_sources] = size(results.S);
num_areas = results.retmap.num_areas;
num_locs = results.retmap.num_locs;

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

tv = results.retmap.areas(1).time*1000;

areas2plot = 1:num_areas;
if ~isempty(areas)
  areas2plot = intersect(areas,areas2plot);
end;
if isempty(areas2plot)
  fprintf('%s: no areas to plot!\n',mfilename);
  return;
end;

num_areas2plot = length(areas2plot);

ncols = floor(sqrt(num_areas2plot));
nrows = ceil(num_areas2plot/ncols);

clf;
namelist = [];
for i=1:num_areas2plot
  a=areas2plot(i);
  namelist = sprintf('%s%s',namelist,results.retmap.areas(a).name);
  h=subplot(nrows,ncols,i);
  tmpdata = [];
  for c=1:parms.ncontrasts
    if loose_flag
      j = 1+(a-1)*num_locs*3 + (c-1)*num_locs*3*num_areas;
      k = j+num_locs*3-1;
      srange = [j:3:k];
    elseif indy_locs_flag
      if isempty(loc)
        j = 1+(a-1)*num_locs + (c-1)*num_locs*num_areas;
        k = j+num_locs-1;
        srange = [j:k];
      else
%% todo: this will break if there are extra ret_dips or nonret_dips
        j = (a-1)*num_locs;
        srange = j+loc;
      end;
    else
      srange = a + (c-1)*num_areas;
    end;
    tmpdata = [tmpdata results.S(:,srange)];
  end;
  
  if isempty(linestyle)
    plot(h,tv,tmpdata,'LineWidth',linewidth);
  else
    plot(h,tv,tmpdata,linestyle,'LineWidth',linewidth);
  end;
  axis(h,'tight');
  if ~isempty(ylim), set(gca,'YLim',ylim); end;
%  title(h,retmap.areas(a).name,'FontSize',fontsize);
%  set(h,'XGrid','on','XTick',[-100 0 100 200 300 400],'GridLineStyle','-');
  set(h,'FontSize',fontsize)
  set(h,'FontName',fontname)
  xlabel('Time (msec)','FontSize',fontsize)
  ylabel('Source Amplitude (nA M)','FontSize',fontsize)
end

outstem = [prefix '-areas'];
if ~isempty(areas)
  outstem = [outstem '-' namelist];
end;
if ~isempty(loc)
  outstem = [outstem '-loc'];
  for i=1:length(loc)
    if (i>1) outstem = [outstem '_']; end;
    outstem = sprintf('%s%d',outstem,loc(i));
  end;
end;

fname = sprintf('%s.tif',outstem);
%print('-dtiff',fname);
%mmil_printeps(gcf,[outstem '.eps']);
