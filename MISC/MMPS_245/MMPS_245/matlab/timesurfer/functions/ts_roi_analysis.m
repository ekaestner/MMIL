function ts_roi_analysis(varargin)
%function ts_roi_analysis([options])
%
% Purpose: extract and plot source waveforms for surface ROIs
%
% Optional Parameters:
%  'rootdir': root directory containing stcfiles directory
%    {default: pwd} (current working directory)
%  'indir': input directory, relative to rootdir unless full path is given
%    {default: 'stcfiles'}
%  'prefix': prefix of stc files
%    If full path is given, rootdir is used only for output
%    If cell array, will loop over multiple prefixes like different conditions
%    {default: 'dSPM'}
%  'conditions' - vector or cell array of conditions to analyze
%    If vector of condition numbers, will look for stcfiles like
%      <prefix>_cond<#>-<hemi>.stc, with # formatted like %02d
%    If cell array of strings, will look for stcfiles like
%      <prefix>_<cond>-<hemi>.stc
%    If empty, use all stc files found in indir
%    {default: []}
%  'outdir': output directory, relative to rootdir unless full path is given
%    {default: 'roi_analysis'}
%  'outstem': output file stem
%    {default: 'roi_results'}
%  'roidir': directory containing label files
%    full path or relative to subjdir/subjname
%    {default = 'label'}
%  'subjname': name of FreeSurfer subject with label files (ROIs)
%    {default = fsaverage}
%  'subjdir': FreeSurfer root directory containing subjname
%    {default = $SUBJECTS_DIR}
%  'ico': icosahedron order number:
%      Order     Number of Vertices
%        1              42
%        2             162
%        3             642
%        4            2562
%        5           10242
%        6           40962
%        7          163842
%     If 0, assumes stcfiles are in native subject vertices
%    {default: 4}
%  'roinames': cell array of ROI names
%      Label files have this format: <hemi>.<roiname>.label
%      If empty, will search for label files in subjdir/subjname/label
%      If none found, will quit with error
%    {default: []}
%  'ico_infix_flag':[0|1] if ico>4, expect label files to have
%    this format: <hemi>.<roiname>.ico<ico>.label
%    {default: 0}
%  'overlay_roi_flag': [0|1] plot with ROIs overlayed on same plot
%    {default: 0}
%  'hemilist': cell array of cortical hemispheres
%    {default: {'lh','rh'}}
%  'plotflag': [0|1] plot waveforms
%    {default: 1}
%  'condnames' - cell array of condition names (must correspond to conditions)
%    If empty, will use conditions
%    {default: []}
%  'prenames' - cell array of shortened names for each prefix
%    Will be prepended to condnames for legend
%    {default: []}
%  'plot_type': 'jpeg' or 'epsc'  or {'jpeg','epsc'}
%    {default: 'jpeg'}
%  'plot_xlim': x axis (time) plot limits (vector of min/max)
%    {default: []}
%  'plot_ylim' y axis (amplitude) plot limits (vector of min/max)
%    {default: []}
%  'plot_offset': subtract this value from plotted waveforms
%    {default: 1}
%  'linewidth': plot line width
%    {default: 1.5}
%  'fontname': font name for labels
%    {default = 'Arial'}
%  'fontsize': font size for labels
%    {default = 12}
%  'label_flag': [0|1] whether to draw text labels
%     {default = 1}
%  'units': amplitude units
%     {default = 'sqrt(F)-1'}
%  'roi_colors' cell array of colors for each roi
%    (only applicable if overlay_roi_flag = 1)
%    {default: {'b','c','g','y','r','m','k'}}
%  'roi_linewidths' vector of linewidhts for each roi
%    (only applicable if overlay_roi_flag = 1)
%    If empty, use linewidth for all
%    {default: []}
%  'legend_flag': [0|1] whether to include legend on plot
%    {default: 1}
%  'legend_loc': location of legend relative to plot
%    {default: 'NorthEastOutside'}
%  'visible_flag': [0|1] make plots visible (otherwise plot in background)
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
%
% Created:   04/14/10 by Don Hagler
% Last Mod:  04/26/13 by Don Hagler
%

%% todo: different linewidths for different ROIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'rootdir',pwd,[],...
  'indir','stcfiles',[],...
  'prefix','dSPM',[],...
  'conditions',[],[],...
  'outdir','roi_analysis',[],...
  'outstem','roi_results',[],...
  'roidir','label',[],...
  'subjname','fsaverage',[],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'ico',4,[0:7],...
  'roinames',[],[],...
  'ico_infix_flag',false,[false true],...
  'overlay_roi_flag',false,[false true],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'plotflag',true,[false true],...
  'condnames',[],[],...
  'prenames',[],[],...
  'plot_type','jpeg',{'jpeg','epsc'},...
  'plot_xlim',[],[],...
  'plot_ylim',[],[],...
  'plot_offset',1,[-Inf,Inf],...
  'visible_flag',false,[false true],...
  'linewidth',1.5,[],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'label_flag',true,[false true],...
  'units','sqrt(F)-1',[],...
  'roi_colors',{'b','c','g','y','r','m','k'},[],...
  'roi_linewidths',[],[],...
  'legend_flag',true,[false true],...
  'legend_loc','NorthEastOutside',[],...
  'avg_hemis_flag',false,[false true],...
  'forceflag',false,[false true],...
});

if parms.plotflag
  if ~iscell(parms.plot_type)
    switch parms.plot_type
      case 'epsc'
        parms.plotext = {'eps'};
      case 'jpeg'
        parms.plotext = {'jpg'};
    end;
    parms.plot_type = {parms.plot_type};
  else
    parms.plotext = cell(size(parms.plot_type));
    for i=1:length(parms.plot_type)
      plot_type = parms.plot_type{i};
      switch plot_type
        case 'epsc'
          parms.plotext{i} = 'eps';
        case 'jpeg'
          parms.plotext{i} = 'jpg';
      end;
    end;
  end;
end;

% check prefix
if ~iscell(parms.prefix)
  parms.prefix_list = {parms.prefix};
else
  parms.prefix_list = parms.prefix;
end;

% check that prenames match prefix_list
if ~isempty(parms.prenames)
  if ~iscell(parms.prenames)
    parms.prenames = {parms.prenames};
  end;
  if length(parms.prenames) ~= length(parms.prefix_list)
    error('length of prenames does not match length of prefix');
  end;
end;

% if not a vector of numbers, make sure is a cell array
if ~isempty(parms.conditions) & ~isnumeric(parms.conditions) &...
   ~iscell(parms.conditions)
  parms.conditions = {parms.conditions};
end;

% check that condnames match conditions
if ~isempty(parms.condnames)
  if ~iscell(parms.condnames)
    parms.condnames = {parms.condnames};
  end;
  if length(parms.condnames) ~= length(parms.conditions)
    error('length of condnames does not match length of conditions');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check indir
if isempty(regexp(parms.indir,'^/')) % relative
  parms.indir = [parms.rootdir '/' parms.indir];
end;
if ~exist(parms.indir,'dir')
  error('input dir %s not found',parms.indir);
end;

% find stcfiles
stcfiles = {};
legend_strs = {};
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  s=1;
  for p=1:length(parms.prefix_list)
    prefix = parms.prefix_list{p};
    if isempty(parms.conditions)
      flist = dir(sprintf('%s/%s*.stc',parms.indir,prefix));
      if isempty(flist)
        fprintf('%s: WARNING: no stc files found in indir %s with prefix %s and hemi %s\n',...
          mfilename,parms.indir,prefix,hemi);
        continue;
      end;
      tmp_stcfiles = {flist.name};
      pat = sprintf('^%s_(?<condname>.+)-%s.stc$',prefix,hemi);
      n = regexp(tmp_stcfiles,pat,'names');
      ind_keep = find(~cellfun(@isempty,n));
      for k=1:length(ind_keep)
        stcfiles{h}{s} = [parms.indir '/' tmp_stcfiles{ind_keep(k)}];
        tmpstr = n{ind_keep(k)}.condname;
        if ~isempty(parms.prenames)
          tmpstr = [parms.prenames{p} ' ' tmpstr];
        end;
        legend_strs{h}{s} = regexprep(tmpstr,'_',' ');
        s=s+1;
      end;
    elseif iscell(parms.conditions)
      for c=1:length(parms.conditions)
        cond = parms.conditions{c};
        if isnumeric(cond), cond = num2str(cond); end;
        tmp_stcfile = sprintf('%s/%s_%s-%s.stc',...
          parms.indir,prefix,cond,hemi);
        if ~exist(tmp_stcfile,'file')
          fprintf('%s: WARNING: stc file %s not found\n',...
            mfilename,tmp_stcfile);
          continue;
        end;        
        stcfiles{h}{s} = tmp_stcfile;
        if isempty(parms.condnames)
          tmpstr = cond;
        else
          tmpstr = parms.condnames{c};
        end;
        if ~isempty(parms.prenames)
          tmpstr = [parms.prenames{p} ' ' tmpstr];
        end;
        legend_strs{h}{s} = regexprep(tmpstr,'_',' ');
        s=s+1;
      end;
    else
      for c=1:length(parms.conditions)
        cond = parms.conditions(c);
        tmp_stcfile = sprintf('%s/%s_cond%02d-%s.stc',...
          parms.indir,prefix,cond,hemi);
        if ~exist(tmp_stcfile,'file')
          fprintf('%s: WARNING: stc file %s not found\n',...
            mfilename,tmp_stcfile);
          continue;
        end;        
        stcfiles{h}{s} = tmp_stcfile;
        if isempty(parms.condnames)
          tmpstr = sprintf('cond%02d',cond);
        else
          tmpstr = parms.condnames{c};
        end;
        if ~isempty(parms.prenames)
          tmpstr = [parms.prenames{p} ' ' tmpstr];
        end;
        legend_strs{h}{s} = regexprep(tmpstr,'_',' ');
        s=s+1;
      end;
    end;
  end;
end;
if isempty(stcfiles)
  error('no stc files found in indir %s',parms.indir);
end;

% check roi dir
if mmil_isrelative(parms.roidir)
  parms.roidir = [parms.subjdir '/' parms.subjname '/' parms.roidir];
end;
if ~exist(parms.roidir,'dir')
  error('dir %s not found',parms.roidir);
end;

% find roi names
if isempty(parms.roinames)
  flist = dir(sprintf('%s/*.label',parms.roidir));
  if isempty(flist)
    error('no label files found in roidir %s',parms.roidir);
  end;
  tmp_roinames = {flist.name};
  if parms.ico>0 && parms.ico_infix_flag
    pat = sprintf('^[lr]h.(?<roiname>.+).ico%d.label$',parms.ico);
  else
    pat = '^[lr]h.(?<roiname>.+).label$';
  end;
  n = regexp(tmp_roinames,pat,'names');
  ind_keep = find(~cellfun(@isempty,n));
  for k=1:length(ind_keep)
    parms.roinames{end+1} = n{ind_keep(k)}.roiname;
  end;
  parms.roinames = unique(parms.roinames);
end;

% create roi file list
roifiles = {};
roinames = {};
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  t=1;
  for r=1:length(parms.roinames)
    if parms.ico>0 && parms.ico_infix_flag
      tmp_roifile = sprintf('%s/%s.%s.ico%d.label',...
        parms.roidir,hemi,parms.roinames{r},parms.ico);
    else
      tmp_roifile = sprintf('%s/%s.%s.label',...
        parms.roidir,hemi,parms.roinames{r});
    end;
    if ~exist(tmp_roifile,'file')
      fprintf('%s: WARNING: label file %s not found\n',...
        mfilename,tmp_roifile);
      continue;
    end;
    roifiles{h}{t} = tmp_roifile;
    roinames{h}{t} = parms.roinames{r};
    t=t+1;
  end;
end;

% create outdir
if mmil_isrelative(parms.outdir)
  parms.outdir = [parms.rootdir '/' parms.outdir];
end;
parms.outstem = [parms.outdir '/' parms.outstem];
mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract waveforms for ROIs
data = [];
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  matfile = sprintf('%s-%s.mat',parms.outstem,hemi);
  matfile_avg = sprintf('%s-avg.mat',parms.outstem);
  if ~exist(matfile,'file') || parms.forceflag
    tmp_stcfiles = stcfiles{h};
    tmp_roifiles = roifiles{h};
    fprintf('%s: extracting waveforms for hemi %s...\n',mfilename,hemi);
    [time,wforms] = ts_stc2roi(tmp_stcfiles,tmp_roifiles);
    if ~isempty(wforms)
      save_results(matfile,time,wforms,roinames{h},parms);
    end;
  end;
end;

if parms.avg_hemis_flag & (~exist(matfile_avg,'file') || parms.forceflag)
  if length(union(roinames{1},roinames{2}))~=length(roinames{1})
    fprintf('%s: WARNING: ROIs do not match for left and right hemispheres... cannot average across hemispheres\n',...
      mfilename);
  else
    data = [];
    for h=1:length(parms.hemilist)
      [time,wforms] = load_results(matfile); % don't get parms
      data(h).time = time;
      data(h).wforms = wforms;
    end;
    time = data(1).time;
    wforms = (data(1).wforms + data(2).wforms)/2;
    if ~isempty(wforms)
      save_results(matfile_avg,time,wforms,roinames{1},parms);
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot waveforms
if parms.plotflag
  plotext = parms.plotext{end};
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    if parms.avg_hemis_flag
      matfile = sprintf('%s-avg.mat',parms.outstem);
    else
      matfile = sprintf('%s-%s.mat',parms.outstem,hemi);
    end;
    [time,wforms] = load_results(matfile); % do not get parms
    if parms.overlay_roi_flag
      num_rois = length(roinames{h});
      num_strs = length(legend_strs{h});
      tmp_legend_strs = cell(num_rois*num_strs,1);
      for r=1:length(tmp_legend_strs)
        for s=1:num_strs
          j = r + (s-1)*num_strs;
          tmp_legend_strs{j} = regexprep([roinames{h}{r} '_' legend_strs{h}{s}],...
            '_',' ');
        end;
      end;
      if parms.avg_hemis_flag
        fname = sprintf('%s-avg.%s',...
          parms.outstem,plotext);
      else
        fname = sprintf('%s-%s.%s',...
          parms.outstem,hemi,plotext);
      end;
      if ~exist(fname,'file') || parms.forceflag
        if parms.avg_hemis_flag
          fprintf('%s: plotting waveforms for avg of hemis...\n',mfilename);
        else
          fprintf('%s: plotting waveforms for hemi %s...\n',mfilename,hemi);
        end;
        clf; hold on;
        if ~parms.visible_flag
          set(gcf,'Visible','Off');
        end;
        num_linewidths = length(parms.roi_linewidths);
        if num_linewidths >= num_rois
          tmp_linewidths = parms.roi_linewidths(1:num_rois);
        else
          tmp_linewidths = parms.linewidth*ones(num_rois,1);
          if num_linewidths>0
            tmp_linewidths(1:num_linewidths) = parms.roi_linewidths(1:num_linewidths);
          end;
        end;
        k = 1;
        for r=1:num_rois
          if k>length(parms.roi_colors)
            k=1;
          end;
          roi_color = parms.roi_colors{k}; k=k+1;
          if size(wforms,3)>1
            wform = squeeze(wforms(:,r,:))-parms.plot_offset;
          else
            wform = squeeze(wforms(r,:))-parms.plot_offset;
          end;
          plot(time,wform,roi_color,'LineWidth',tmp_linewidths(r));
        end;
        save_plot(parms,fname,tmp_legend_strs);
      end;
    else
      for r=1:length(roinames{h})
        if parms.avg_hemis_flag
          fname = sprintf('%s-%s-avg.%s',...
            parms.outstem,roinames{h}{r},plotext);
        else
          fname = sprintf('%s-%s-%s.%s',...
            parms.outstem,roinames{h}{r},hemi,plotext);
        end;
        tmp_legend_strs = legend_strs{h};
        if ~exist(fname,'file') || parms.forceflag
          if parms.avg_hemis_flag
            fprintf('%s: plotting waveforms for avg of hemis...\n',mfilename);
          else
            fprintf('%s: plotting waveforms for hemi %s...\n',mfilename,hemi);
          end;
          if size(wforms,3)>1
            wform = squeeze(wforms(:,r,:))-parms.plot_offset;
          else
            wform = squeeze(wforms(r,:))-parms.plot_offset;
          end;
          clf;
          if ~parms.visible_flag
            set(gcf,'Visible','Off');
          end;
          %% todo: loop over conditions to set color for each according
          %%  to some input parameter
          plot(time,wform,'LineWidth',parms.linewidth);
          save_plot(parms,fname,tmp_legend_strs);
          if ~parms.visible_flag
            close(gcf);
          end;
        end;
      end;
    end;
    if parms.avg_hemis_flag, break; end;
  end;
end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function save_plot(parms,fname,legend_strs)
  if ~isempty(parms.plot_xlim)
    set(gca,'XLim',parms.plot_xlim);
  end;
  if ~isempty(parms.plot_ylim)
    set(gca,'YLim',parms.plot_ylim);
  end;
  if parms.legend_flag
    legend(legend_strs,'Location',parms.legend_loc,...
      'FontSize',parms.fontsize,'FontName',parms.fontname);
  end;
  if parms.label_flag
    set(gca,'FontSize',parms.fontsize,'FontName',parms.fontname)
    xlabel('Time (msec)','FontSize',parms.fontsize,'FontName',parms.fontname)
    title('Source Waveforms','FontSize',parms.fontsize,'FontName',parms.fontname)
    ylabel(sprintf('Source Amplitude (%s)',parms.units),...
      'FontSize',parms.fontsize,'FontName',parms.fontname)
  end;
  for i=1:length(parms.plot_type)
    if ~parms.visible_flag
      set(gcf,'Visible','Off');
    end;
    plot_type = parms.plot_type{i};
    [fpath,fstem,fext]=fileparts(fname);
    tmp_fname = [fpath '/' fstem '.' parms.plotext{i}];
    switch plot_type
      case 'epsc'
        mmil_printeps(gcf,tmp_fname);
      otherwise
        print(['-d' plot_type],tmp_fname);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,wforms] = load_results(matfile)
  time = [];
  wforms = [];
  load(matfile);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_results(matfile,time,wforms,roinames,parms)
  save(matfile,'time','wforms','roinames','parms');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

