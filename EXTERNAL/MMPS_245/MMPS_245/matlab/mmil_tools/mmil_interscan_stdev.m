function results = mmil_interscan_stdev(fnamelist,varargin)
%function results = mmil_interscan_stdev(fnamelist,[options])
%
% Required input:
%   fnamelist: cell array of file names (mgh/mgz format)
%
% Optional input:
%   'outdir': output directory
%     {default = [pwd '/output]}
%   'outstem': output file stem
%     {default: 'data'}
%   'plotflag': [0|1] plot images of inter- and intra-scan variability
%     {default = 0}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  08/14/2017 by Don Hagler
% Last Mod: 08/15/2017 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'outdir',[pwd '/output'],[],...
  'outstem','data',[],...
  'plotflag',false,[false true],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
  ...
  'fig_size',[],[],...
  'tif_dpi',300,[10,10000],...
  'planestrings',{'HOR','SAG','COR'},[],...
  'slice_fract',[0.5,0.6,0.5],[0,1],...
  ... % hidden for epi_brainmask
  'log_flag',true,[false true],...
  'thresh',0.75,[],... %% todo: may need way to choose best value for
          ...          %%       different levels of intensity inhomogeneity
          ...          %%   or, need to correct for intensity inhomogeneity
  'fill1_smooth1',10,[],...
  'fill1_thresh1',0.95,[],...
  'fill1_smooth2',30,[],... % different from default
  'fill1_thresh2',0.1,[],...
  'fill1_smooth3',0,[],...
  'fill2_smooth1',20,[],...
  'fill2_thresh1',0.95,[],...
  'fill2_smooth2',30,[],... % different from default
  'fill2_thresh2',0.1,[],...
  'fill2_smooth3',20,[],...
  'fill2_thresh3',0.1,[],... % different from default
  'binary_flag',true,[false true],...
  ...
  'plot_tags',{'outdir','outstem','border_flag','fig_size','tif_dpi',...
    'frame','planestrings','slice_fract','clim','visible_flag','forceflag',...
    'plim_flag','plim'},[],...
  'epi_tags',{'log_flag','thresh',...
              'fill1_smooth1','fill1_thresh1',...
              'fill1_smooth2','fill1_thresh2',...
              'fill1_smooth3','fill1_erode_flag',...
              'clip_edges_flag','clip_edges_width',...
              'fill2_smooth1','fill2_thresh1',...
              'fill2_smooth2','fill2_thresh2',...
              'fill2_smooth3','fill2_thresh3','fill2_erode_flag',...
              'binary_flag','forceflag'},[],...
});

if ~iscell(fnamelist), fnamelist = {fnamelist}; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for missing files
ind_missing = find(~cellfun(@(x) exist(x,'file'),fnamelist));
if parms.verbose
  nummissing = length(ind_missing);
  if nummissing>0
    fprintf('%s: WARNING: excluding %d/%d files that do not exist:\n',...
      mfilename,nummissing,length(fnamelist));
    for i=1:nummissing
      fprintf('%s\n',fnamelist{ind_missing(i)});
    end;
  end;
end;
ind_exist = find(cellfun(@(x) exist(x,'file'),fnamelist));
if isempty(ind_exist)
  fprintf('%s: ERROR: all files in fnamelist are missing\n',mfilename);
  return;
end;
fnamelist = fnamelist(ind_exist);
nfiles = length(fnamelist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

% save first frame of each volume timeseries
fnamelist_first = {};
for i=1:nfiles
  fname_in = fnamelist{i};
  [fpath,fstem,fext] = fileparts(fname_in);
  fname_out = sprintf('%s/%s_first.mgz',parms.outdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: saving first frame of %s...\n',mfilename,fname_in);
    end;
    [vol,M] = fs_load_mgh(fname_in,[],1);
    fs_save_mgh(vol,fname_out,M);
  end;
  clear vol;
  fnamelist_first{i} = fname_out;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_first_mean = sprintf('%s/%s_interscan_first_mean.mgz',parms.outdir,parms.outstem);
fname_first_stdev = sprintf('%s/%s_interscan_first_stdev.mgz',parms.outdir,parms.outstem);
if ~exist(fname_first_mean,'file') || ~exist(fname_first_stdev,'file') || parms.forceflag
  vol_first = [];
  for i=1:nfiles
    fname_in = fnamelist_first{i};
    if parms.verbose
      fprintf('%s: loading %s...\n',mfilename,fname_in);
    end;
    [vol,M] = fs_load_mgh(fname_in);
    vol_first = cat(4,vol_first,vol);
  end;
  if parms.verbose
    fprintf('%s: calculating interscan mean of first frame...\n',mfilename);
  end;
  vol_first_mean = mean(vol_first,4);
  fs_save_mgh(vol_first_mean,fname_first_mean,M);
  if parms.verbose
    fprintf('%s: calculating interscan stdev of first frame...\n',mfilename);
  end;
  vol_first_stdev = std(vol_first,[],4);
  fs_save_mgh(vol_first_stdev,fname_first_stdev,M);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create brain mask from mean image
[fpath,fstem,fext] = fileparts(fname_first_mean);
fname_mask = sprintf('%s/%s_mask%s',parms.outdir,fstem,fext);
args = mmil_parms2args(parms,parms.epi_tags);
epi_brainmask(fname_first_mean,fname_mask,args{:});

% summarize results in mat and json files
fname_mat = sprintf('%s/%s_results.mat',parms.outdir,parms.outstem);
fname_json = sprintf('%s/%s_results.json',parms.outdir,parms.outstem);
if ~exist(fname_mat,'file') || ~exist(fname_json,'file') || parms.forceflag
  % load mask
  if parms.verbose
    fprintf('%s: loading %s...\n',mfilename,fname_mask);
  end;
  vol_mask = fs_load_mgh(fname_mask);
  ind_mask = find(vol_mask>0);

  % calculate average stdev value within brain mask
  if parms.verbose
    fprintf('%s: calculating average stdev within brain mask...\n',mfilename);
  end;
  [vol,M] = fs_load_mgh(fname_first_stdev,[],1);
  vals = vol(ind_mask);
  results.stdev_mean = mean(vals);
  results.stdev_median = median(vals);

  % calculate average stdev value within brain mask
  if parms.verbose
    fprintf('%s: calculating average mean within brain mask...\n',mfilename);
  end;
  [vol,M] = fs_load_mgh(fname_first_mean,[],1);
  vals = vol(ind_mask);
  results.mean_mean = mean(vals);
  results.mean_median = median(vals);

  % calculate coefficient of variation
  results.coef_var_mean = results.stdev_mean / results.mean_mean;
  results.coef_var_median = results.stdev_median / results.mean_median;

  save(fname_mat,'results');
  savejson('',results,'filename',fname_json);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot images for each output file
if parms.plotflag
  if parms.verbose
    fprintf('%s: plotting interscan mean of first frame...\n',mfilename);
  end;
  tp = parms;
  tp.outstem = sprintf('%s_interscan_first_mean',parms.outstem);
  args = mmil_parms2args(tp,parms.plot_tags);
  mmil_subplots(fname_first_mean,args{:});

  if parms.verbose
    fprintf('%s: plotting mask of interscan mean of first frame...\n',mfilename);
  end;
  tp = parms;
  tp.outstem = sprintf('%s_interscan_first_mean_mask',parms.outstem);
  args = mmil_parms2args(tp,parms.plot_tags);
  mmil_subplots(fname_mask,args{:});
  
  if parms.verbose
    fprintf('%s: plotting interscan stdev of first frame...\n',mfilename);
  end;
  tp = parms;
  tp.outstem = sprintf('%s_interscan_first_stdev',parms.outstem);
  args = mmil_parms2args(tp,parms.plot_tags);
  mmil_subplots(fname_first_stdev,args{:});
  
  for i=1:nfiles
    fname_in = fnamelist_first{i};
    [fpath,fstem,fext] = fileparts(fname_in);
    if parms.verbose
      fprintf('%s: plotting first frame of %s...\n',mfilename,fname_in);
    end;
    tp = parms;
    tp.outstem = fstem;
    args = mmil_parms2args(tp,parms.plot_tags);
    mmil_subplots(fname_in,args{:});
  end;  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

