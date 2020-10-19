function [results,wforms_fit] = rc_fit_wforms(wforms,varargin)
%function [results,wforms_fit] = rc_fit_wforms(wforms,[options])
%
% Purpose: fit source estimate waveforms with multiple components
%
% Required Input:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nareas,nconditions]
%
% Optional Input:
%   'condition_values': vector of x-axis values for each condition
%     number must equal nconditions; if empty, will use 1:nconditions
%     {default = []}
%   'area_names': cell array of area names (must match nareas)
%     {default = {'V1','V2','V3'}}
%   'area_colors': cell array of colors for each area (must match nareas)
%     {default = {'b','g','r'}}
%   'search_type': algorithm for search
%       'rand' = random search   'fmincon': = fmincon MATLAB function
%     {default = fmincon}
%   'niters': number of iterations for random or fmincon search
%     {default = 500}
%   'stepsize': size of random search step relative to range of bounds
%     {default = 0.05}
%   'rand_init_flag': [0|1] use random starting values for parameter estimates
%     otherwise use middle of bounds
%   'ncomponents': number of components to model
%     {default = 1}
%   'polarity': waveform component polarity (may be vector with ncomponents)
%     {default = -1}
%   'latency': initial estimate(s) for component latency (msec)
%     if ncomponents>1, should be matrix with nareas x ncomponents
%     {default = []}
%   'amplitude': initial estimate(s) for component amplitude
%     if ncomponents>1, should be matrix with nareas x ncomponents
%     {default = []}
%   'rise_tc':initial estimate(s) for exponential rise time constant (msec)
%     if ncomponents>1, should be matrix with nareas x ncomponents
%     {default = []}
%   'fall_tc': initial estimate(s) for exponential fall time constant (msec)
%     if ncomponents>1, should be matrix with nareas x ncomponents
%     {default = []}
%   'latency_bounds': vector of lower and upper bounds for latency
%     May be matrix with size = [nareas,ncomponents,2]
%     {default = [40,120]}
%   'amplitude_bounds': lower and upper bounds for amplitude
%     May be matrix with size = [nareas,ncomponents,2]
%     {default = [0.5,25]}
%   'rise_tc_bounds': exponential rise time constant (msec)
%     May be matrix with size = [nareas,ncomponents,2]
%     {default = [1,20]}
%   'fall_tc_bounds': exponential fall time constant (msec)
%     May be matrix with size = [nareas,ncomponents,2]
%     {default = [1,20]}
%   'sfreq': sampling frequency; used to calculate time vector
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 400}
%   'time': time vector (msec)
%     if supplied, sfreq, t0, and t1 are ignored
%     length of time vector must match length of input wform
%     {default = []}
%   'contrast_latency_flag': [0|1] require that higher condition
%     levels have equal or shorter latency
%     {default = 0}
%   'outdir': output directory; if none supplied, no output
%     otherwise, will save plots as tifs and results in mat file
%     {default = []}
%   'outstem': output file stem (relative to outdir or full path)
%     {default = 'wfit'}
%
% Created:  05/10/11 by Don Hagler
% Last Mod: 11/19/15 by Don Hagler
%

%% todo: allow bounds to be [ncomponents,2]

%% todo: allow restricted time range for fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'condition_values',[],[],...
  'area_names',{'V1','V2','V3'},[],...
  'area_colors',{'b','g','r'},[],...
  'niters',500,[0,1e10],...
  'search_type','rand',{'rand','fmincon'},...
  'stepsize',0.05,[],...
  'rand_init_flag',true,[false true],...
  'ncomponents',1,[],...
  'polarity',-1,[],...
...
  'latency',[],[],...
  'amplitude',[],[],...
  'rise_tc',[],[],...
  'fall_tc',[],[],...
...
  'latency_bounds',[40,120],[],...
  'amplitude_bounds',[0.5,25],[],...
  'rise_tc_bounds',[1,30],[],...
  'fall_tc_bounds',[1,30],[],...
...
  'wfit_func','sigmoid',{'sigmoid','gamma'},...
  'single_tc_flag',false,[false true],...
...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
  'contrast_latency_flag',false,[false true],...
  'outdir',[],[],...
  'outstem','wfit',[],...
...
  'delay_sf',0,[],...
...
  'DiffMinChange',0.01,[],...
  'TolFun',1e-5,[],...
  'TolX',1e-6,[],...
...
  'plotflag',true,[false true],...
  'xlim',[-100,300],[],...
  'ylim',[-25,10],[],...
  'use_areas',[],[],...
  'linewidth_data',1,[],...
  'linewidth_fit',2,[],...
...
  'visible_flag',true,[false true],...
  'forceflag',false,[false true],...
};

excl_tags = {'plotflag','xlim','ylim','use_areas','condition_values',...
  'area_names','area_colors','linewidth_data','linewidth_fit',...
  'outdir','outstem','visible_flag','contrast_latency_flag',...
  'beta_names','nbeta','ntpoints','nareas','nconditions','forceflag'};

results = [];
wforms_fit = [];

parms = check_input(wforms,varargin,parms_filter);

tags = setdiff(fieldnames(parms),excl_tags);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mat = [parms.outstem '_results.mat'];
if isempty(parms.outdir) || (~exist(fname_mat,'file') || parms.forceflag)
  fprintf('%s: fitting waveforms...\n',mfilename);
  wforms_fit = zeros(size(wforms));
  clear results;
  for c=parms.nconditions:-1:1
    if parms.plotflag
      figure(c); clf; hold on;
      legstrs = [];
    end;
    for s=parms.use_areas
      wform = wforms(:,s,c);
      tmp_parms = parms;
      for b=1:parms.nbeta
        bname = parms.beta_names{b};
        tmp_parms.(bname) = parms.(bname)(s,:);
        bname = [parms.beta_names{b} '_bounds'];
        tmp_parms.(bname) = ...
          reshape(parms.(bname)(s,1:parms.ncomponents,:),[parms.ncomponents,2]);
      end;

      % reset latency_bounds so lower contrast conditions have longer latency
      if parms.contrast_latency_flag && c<parms.nconditions
        tmp_bounds = tmp_parms.latency_bounds;
        tmp_latency = results(s,c+1).latency;
        tmp_bounds(:,1) = tmp_latency;
        tmp = diff(tmp_bounds,[],2);
        ind = find(tmp<1);
        if ~isempty(ind)
          tmp_bounds(ind,2) = tmp_bounds(ind,2) + 1;
        end;
        tmp_parms.latency_bounds = tmp_bounds;
      end;

      args = mmil_parms2args(tmp_parms,tags);
      tmp_results = rc_fit_wform_components(wform,args{:});
      tmp_results.area_name = parms.area_names{s};
      tmp_results.condition = parms.condition_values(c);
      results(s,c) = tmp_results;
      wforms_fit(:,s,c) = tmp_results.wform_fit;
      fprintf('%s: area %s, condition %0.2f, error = %0.3f\n',...
        mfilename,parms.area_names{s},parms.condition_values(c),tmp_results.min_err);
      if parms.plotflag
        plot(tmp_results.time,wform,parms.area_colors{s},'LineWidth',parms.linewidth_data);
        legstrs{end+1} = sprintf('%s data',parms.area_names{s});
        plot(tmp_results.time,tmp_results.wform_fit,parms.area_colors{s},'LineWidth',parms.linewidth_fit);
        legstrs{end+1} = sprintf('%s fit',parms.area_names{s});
      end;
    end;
    if parms.plotflag
      title(sprintf('condition = %0.2f',parms.condition_values(c)));
      legend(legstrs,'Location','SouthEast');
      set(gca,'XLim',parms.xlim,'YLim',parms.ylim);
      if ~parms.visible_flag
        set(gcf,'Visible','Off');
      end;
      if ~isempty(parms.outdir)
        fname_plot = sprintf('%s_cond%0.2f.tif',...
          parms.outstem,parms.condition_values(c));
        print(gcf,'-dtiff',fname_plot);
      end;
      if ~parms.visible_flag
        close(gcf);
      end;
    end;
  end;
  if ~isempty(parms.outdir)
    save(fname_mat,'results');
  end;
else
  load(fname_mat);
  for c=1:parms.nconditions
    for s=parms.use_areas
      wforms_fit(:,s,c) = results(s,c).wform_fit;
    end;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.ntpoints = size(wforms,1);
  parms.nareas = size(wforms,2);
  parms.nconditions = size(wforms,3);
  parms.beta_names = {'latency','amplitude','rise_tc','fall_tc'};
  parms.nbeta = length(parms.beta_names);
  if isempty(parms.use_areas), parms.use_areas = [1:parms.nareas]; end;
  if isempty(parms.condition_values)
    parms.condition_values = 1:parms.nconditions;
  end;
  if ~isempty(parms.outdir)
    mmil_mkdir(parms.outdir);
  end;
  if mmil_isrelative(parms.outstem)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
  for b=1:parms.nbeta
    bname = [parms.beta_names{b} '_bounds'];
    bsz = size(parms.(bname));
    if numel(parms.(bname))==2
      bounds = zeros(parms.nareas,parms.ncomponents,2);
      for j=1:2
        bounds(:,:,j) = parms.(bname)(j);
      end;
      parms.(bname) = bounds;
    elseif any(bsz~=[parms.nareas,parms.ncomponents,2])
      error('%s has wrong size [ %s] instead of [ %d %d 2 ]',...
        bname,sprintf('%d ',bsz),parms.nareas,parms.ncomponents);
    end;
  end;
return;

