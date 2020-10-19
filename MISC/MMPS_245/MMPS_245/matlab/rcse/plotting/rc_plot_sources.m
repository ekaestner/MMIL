function rc_plot_sources(varargin)
%function rc_plot_sources([options])
%
% Required Parameters:
%   'prefix': prefix of output files from a run of retinv
%     if supplied, will load parms and results from matfiles
%     {default = []}
%   'parms': parms struct from retinv
%     Required if prefix is empty
%     {default = []}
%   'results': results struct from retinv
%     Required if prefix is empty
%     {default = []}
%
% Optional Parameters:
%   'rootdir': root directory for retinv analysis (containing matfiles dir)
%     {default = pwd}
%   'area_colors': cell array of colors for each area waveform (V1,V2,V3)
%     {default = {'b','g','r'}}
%   'area_names': cell array of area names
%     {default = {'V1','V2','V3'}}
%   'linewidth': plotted line width
%     {default = 1.5}
%   'areas': vector of areas to include in plot
%     If empty, plot all
%     {default = []}
%   'contrasts': vector contrast levels to include in plot
%     If empty, plot all
%     {default = []}
%   'ylim': vector of y-axis limits
%     If empty, auto-scale y-axis
%     {default = []}
%   'xlim': vector of x-axis limits
%     If empty, auto-scale x-axis
%     {default = []}
%   'slopeflag': [0|1] whether to plot derivative of waveforms
%     {default = 0}
%   'normflag': [0|1] whether to normalize waveforms to peak
%     {default = 0}
%   'avgindyflag': [0|1] if retinv done with indy locs,
%     whether to plot individual source estimates (0) or average (1)
%     {default = 0}
%   'labelflag': [0|1] whether to draw text labels
%     {default = 1}
%   'units': amplitude units
%     {default = 'nA m'}
%   'fontname': font name for labels
%     {default = 'Arial'}
%   'fontsize': font size for labels
%     {default = 12}
%   'err_samp_dur': number of samples between error bars
%     {default = 15}
%
% Created:   06/04/10 by Don Hagler
% Last Mod:  12/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get parms from varargin

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix',[],[],...
  'parms',[],[],...
  'results',[],[],...
  'rootdir',pwd,[],...
  'area_colors',{'b','g','r'},[],...
  'area_names',{'v1','v2','v3'},[],...
  'linewidth',1.5,[],...
  'areas',[],[],...
  'contrasts',[],[],...
  'ylim',[],[],...
  'xlim',[],[],...
  'slopeflag',false,[false true],...
  'normflag',false,[false true],...
  'avgindyflag',false,[false true],...
  'labelflag',true,[false true],...
  'units','nA m',[],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'err_samp_dur',15,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load retinv parms and results if prefix supplied

if ~isempty(parms.prefix)
  matfile = sprintf('%s/matfiles/%s_parms.mat',parms.rootdir,parms.prefix);
  if ~exist(matfile,'file')
    error('file %s not found',matfile);
  end;
  tmp = load(matfile);
  parms.parms = tmp.parms;
  matfile = sprintf('%s/matfiles/%s_results.mat',parms.rootdir,parms.prefix);
  if ~exist(matfile,'file')
    error('file %s not found',matfile);
  end;
  tmp = load(matfile);
  parms.results = tmp.results;
elseif isempty(parms.parms) | isempty(parms.results)
  error('prefix not supplied and parms or results are empty');
end;

if parms.normflag
  parms.units = 'arbitrary units';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check some retinv parms

if isfield(parms.parms,'indy_locs_flag')
  parms.indy_locs_flag = parms.parms.indy_locs_flag;
else
  parms.indy_locs_flag = 0;
end;
if isfield(parms.parms,'loose_flag')
  parms.loose_flag = parms.parms.loose_flag;
else
  parms.loose_flag = 0;
end;

if isfield(parms.parms,'ncontrasts')
  parms.ncontrasts = parms.parms.ncontrasts;
else
  parms.ncontrasts = 1;
end;

if isempty(parms.contrasts)
  parms.contrasts = [1:parms.ncontrasts];
elseif any(~ismember(parms.contrasts,[1:parms.ncontrasts]))
  error('specified contrast numbers must be less than %d',...
    parms.ncontrasts);
end;
[num_tpoints num_sources num_contrasts] = size(parms.results.S);
if num_contrasts ~= parms.ncontrasts
  error('number of source contrasts (%d) does not match ncontrasts (%d)',...
    num_contrasts,parms.ncontrasts);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

retmap = parms.results.retmap;
num_areas = retmap.num_areas;
num_locs = retmap.num_locs;
tv = parms.results.time*1000;
timestep = tv(2)-tv(1);

areas2plot = 1:num_areas;
if ~isempty(parms.areas)
  areas2plot = intersect(parms.areas,areas2plot);
end;
if isempty(areas2plot)
  error('no areas to plot');
end;

num_areas2plot = length(areas2plot);

legendstrings=[];
parms.new_area_colors = cell(num_areas2plot,1);
for i=1:num_areas2plot
  a=areas2plot(i);
  area_name = upper(retmap.areas(a).name);
  legendstrings{i} = area_name;
  ind_area = find(strcmp(area_name,upper(parms.area_names)));
  if ~isempty(ind_area)
    parms.new_area_colors{i} = parms.area_colors{ind_area};
  else
    parms.new_area_colors{i} = parms.area_colors{i};
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
hold on;
for c=parms.contrasts
  if length(parms.linewidth)>1
    tmp_linewidth = parms.linewidth(c);
  else
    tmp_linewidth = parms.linewidth;
  end;
  data = [];
  stdv = [];
  for i=1:num_areas2plot
    a=areas2plot(i);
    if parms.loose_flag
      j = 1+(a-1)*num_locs*3;
      k = j+num_locs*3-1;
      srange = [j:3:k];
    elseif parms.indy_locs_flag
      j = 1+(a-1)*num_locs;
      k = j+num_locs-1;
      srange = [j:k];
    else
      srange = a;
    end;
    tmpdata = squeeze(parms.results.S(:,srange,c));
    if isfield(parms.results,'S_stdv')
      tmpstdv = squeeze(parms.results.S_stdv(:,srange,c));
    else
      tmpstdv = [];
    end;
    if parms.indy_locs_flag & parms.avgindyflag
      tmpdata = mean(tmpdata,2);
    end;
    if parms.slopeflag
      tmpdata = [0;diff(tmpdata)/timestep];
    end;
    if parms.normflag
      tmpdata = tmpdata / max(abs(tmpdata(:)));
    end;
    if parms.indy_locs_flag & ~parms.avgindyflag
      if isempty(data)
        data = zeros([size(tmpdata),num_areas2plot]);
      end;
      data(:,:,i) = tmpdata;
      if ~isempty(tmpstdv)
        if isempty(stdv)
          stdv = zeros([size(tmpstdv),num_areas2plot]);
        end;
        stdv(:,:,i) = tmpstdv;
      end;
    else
      data = [data tmpdata];
      stdv = [stdv tmpstdv];
    end;
  end

  if parms.indy_locs_flag & ~parms.avgindyflag
    for s=1:num_locs
      for i=1:num_areas2plot
        a=areas2plot(i);
        plot(tv,data(:,s,i),'Color',parms.new_area_colors{a},...
          'LineWidth',tmp_linewidth);
      end;
      if ~isempty(stdv)
        for i=1:num_areas2plot
          a=areas2plot(i);
          errorbar(tv(1:parms.err_samp_dur:end),...
            data(1:parms.err_samp_dur:end,s,i),...
            stdv(1:parms.err_samp_dur:end,s,i),'.',...
            'Color',parms.new_area_colors{a},'LineWidth',tmp_linewidth);
        end;
      end;
    end;
  else
    for i=1:num_areas2plot
      a=areas2plot(i);
      plot(tv,data(:,i),'Color',parms.new_area_colors{a},'LineWidth',tmp_linewidth);
      if ~isempty(stdv)
        errorbar(tv(1:parms.err_samp_dur:end),...
          data(1:parms.err_samp_dur:end,i),...
          stdv(1:parms.err_samp_dur:end,i),'.',...
          'Color',parms.new_area_colors{a},'LineWidth',tmp_linewidth);
      end;  
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis('tight');
if ~isempty(parms.ylim), set(gca,'YLim',parms.ylim); end;
if ~isempty(parms.xlim), set(gca,'XLim',parms.xlim); end;

if parms.labelflag
  set(gca,'FontSize',parms.fontsize,'FontName',parms.fontname)
  legend(legendstrings,'Location','SouthEast',...
    'FontSize',parms.fontsize,'FontName',parms.fontname);
  xlabel('Time (msec)','FontSize',parms.fontsize,'FontName',parms.fontname)
  title('Source Waveforms','FontSize',parms.fontsize,'FontName',parms.fontname)
  if parms.slopeflag
    ylabel(sprintf('Slope of Source Amplitude (%s / ms)',parms.units),...
      'FontSize',parms.fontsize,'FontName',parms.fontname)
  else
    ylabel(sprintf('Source Amplitude (%s)',parms.units),...
      'FontSize',parms.fontsize,'FontName',parms.fontname)
  end;
end;

