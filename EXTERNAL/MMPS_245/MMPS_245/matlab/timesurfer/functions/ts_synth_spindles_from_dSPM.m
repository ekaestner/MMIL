function G_norm = ts_synth_spindles_from_dSPM(fspath,varargin)
%function G_norm = ts_synth_spindles_from_dSPM(fspath,[options])
%
% Usage:
%  G_norm = ts_synth_spindles_from_dSPM(fspath,'key1', value1,...);
%
% Required Input:
%   fspath: full path of FreeSurfer recon
%     NOTE: this function assumes that dSPM forward solution was
%       calculated from ALL surface vertices in fspath surfaces
%
% Optional parameters for input dSPM:
%   'rootdir': directory containing matfiles dir
%     {default = pwd}
%   'prefix': prefix of ts_dSPM output files
%     {default = 'dSPM'}
%   'prefix_data': prefix of processed data file (e.g. 'proc_avg_data')
%     if not supplied, will use avg_data saved by ts_dSPM
%     {default = []}
%   'G_norm': normal vector gain matrix
%     NOTE: may supply this to avoid having to calculate from G_xyz
%     {default = []}
%
% Optional parameters for output:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'spindles'}
%   'source_mat_flag': [0|1] save source time courses in matlab file format
%     otherwise use mgz
%     {default = 1}
%   'save_data_flag': [0|1] save source and sensor data
%     otherwise just create plots and calculate correlations
%     {default = 1}
%   'niters': number of iterations
%     {default = 1}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Optional parameters for simulation time window:
%   'time_dur': duration of simulation time window (msec)
%     {default = 1000}
%   'sf': sampling frequency (Hz)
%     {default = 500}
%
% Optional parameters for spindle characteristics:
%   'uniform_flag': [0|1] create identical spindle for all dipoles
%     {default = 0}
%   'indy_flag': [0|1] create independent spindles for each dipole
%     {default = 1}
%   'npatches': number of discrete patches (split evenly across hemispheres)
%     if 0, and indy_flag = 1, create independent spindles at every dipole
%     {default = 0}
%   'uniform_patches_flag': [0|1] create identical spindle for all patches
%     ignored if npatches = 0
%     {default = 0}
%   'npatches_pool': number of patches in pool from which patches are drawn
%     if npatches_pool < npatches, will be set to npatches
%     {default = 0}
%   'patch_radius': number of sparse smoothing steps applied to dilate patches
%     {default = 10}
%   'freq_mean': mean spindle frequency
%     {default = 13}
%   'freq_std': standard deviation of spindle frequency
%     {default = 1}
%   'phase_mean': mean spindle phase (cycles)
%     {default = 0}
%   'phase_std': st.dev. of spindle phase (cycles)
%     {default = 1}
%   'phase_gauss_flag': [0|1] use Gaussian distribution for spindle phases
%     otherwise use uniform distribution with range of phase_std
%     {default = 0}
%   'amp_mean': mean spindle amplitude (nA m)
%     {default = 1}
%   'amp_std': st.dev. of spindle amplitude (nA m)
%     {default = 0}
%   'amp_uniform': spindle amplitude (nA m) if uniform_flag = 1
%     {default = 0.1}
%   't0_min': minimum spindle start time (msec)
%     {default = 100}
%   't0_mean': mean spindle start time (msec)
%     {default = 250}
%   't0_std': st.dev. of spindle start time (msec)
%     {default = 50}
%   'td_min': minimum spindle duration (msec)
%     {default = 200}
%   'td_mean': mean spindle duration (msec)
%     {default = 500}
%   'td_std': st.dev. of spindle duration (msec)
%     {default = 50}
%   'source_noise': st.dev. of Gaussian noise added to each dipole (nA m)
%     {default = 0.1}
%
% Optional parameters for plotting source waveforms:
%   'source_plot_flag': [0|1] create plots of source time courses
%     {default = 1}
%   'source_plot_ylim': scale limits for source plots
%     {default = [-2,2]}
%   'source_image_flag': [0|1] create images of source time courses
%     {default = 1}
%   'source_image_ylim': scale limits for source images
%     {deafult = [-1,1]}
%
% Optional parameters for simulating sensor data:
%   'sensor_noise': st.dev.of Gaussian noise added to each sensor
%     relative to st.dev. of noise in avg_data (raw noise covariance)
%     {default = 0}
%   'grad_noise': stdev of gradiometer noise, in fT/cm
%     may use instead of or in addition to sensor_noise
%     {default = 0}
%   'mag_noise': stdev of magnetometer noise, in fT
%     {default = 0}
%   'EEG_noise': stdev of EEG noise, in uV
%     {default = 0}
%   'bandpass_flag': [0|1] bandpass filter
%     {default = 1}
%   'bandpass_low_cf': low cutoff frequency (high-pass filter) (Hz)
%     {default = 8}
%   'bandpass_low_tb': low cutoff transition band (Hz)
%     {default = 0.1}
%   'bandpass_high_cf': high cutoff frequency (low-pass filter) (Hz)
%     {default = 18}
%   'bandpass_high_tb': high cutoff transition band (Hz)
%     {default = 0.1}
%   'notch_flag': [0|1] notch filter
%     {default = 0}
%   'notch_cf': notch center frequency (notch filter) (Hz)
%     {default = 0}
%   'notch_tb': notch transition band (Hz)
%     {default = 0}
%
% Optional parameters for plotting sensor waveforms:
%   'sensor_plot_flag': [0|1] create plots of sensor time courses
%     {default = 1}
%   'sensor_plot_grad_ylim': scale limits for gradiometer sensor plots
%     {default = [-200,200]}
%   'sensor_plot_mag_ylim': scale limits for magnetometer sensor plots
%     {default = [-400,400]}
%   'sensor_plot_EEG_ylim': scale limits for EEG sensor plots
%     {default = [-30,30]}
%   'sensor_image_flag': [0|1] create images of sensor time courses
%     {default = 1}
%   'sensor_image_grad_ylim': scale limits for gradiometer sensor images
%     {default = [-400,400]}
%   'sensor_image_mag_ylim': scale limits for magnetometer sensor images
%     {default = [-800,800]}
%   'sensor_image_EEG_ylim': scale limits for EEG sensor images
%     {default = [-60,60]}
%   'sensor_topo_flag': [0|1] create topography plot of sensors at the
%     time of greatest magnitude
%     {default = 1}
%   'grad_scalefact': scaling factor applied to gradiometer data
%     purpose of scaling factors is to get data from different channel types
%     into roughly the same scale
%     {default = 10^13}
%   'mag_scalefact': scaling factor applied to magnetometer data
%     {default = 10^15}
%   'EEG_scalefact': scaling factor applied to EEG data
%     {default = 10^6}
%   'sensor_types': sensor types to plot
%     {default = {'grad','mag','EEG'}}
%
% Optional parameters for plotting sensor correlations:
%   'sensor_corr_image_flag': [0|1] create images of sensor correlations matrices
%     {default = 1}
%
% Optional parameters for correlation calculations:
%   'corr_t0': start time (msec) of correlation time window
%     {default = 100}
%   'corr_t1': stop time (msec) of correlation time window
%     {default = 1000}
%   'corr_source_subsamp': number of subsampled sources
%     for source correlation calculation
%     {default = 500}
%
% Output:
%   G_norm: normal vector gain matrix
%     NOTE: may use this as input for subsequent calls of this function
%       to avoid having to recalculate this each time
%
% Created:  12/28/13 by Don Hagler
% Last Mod: 08/04/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: smoothing option
%% todo: multiple spindles per dipole with refractory period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_norm = [];

% parse input parameters
parms = check_input(fspath,varargin);

% check if output files exist
all_exist = check_output(parms);
if all_exist, return; end;

% create patches
if parms.npatches
  parms = create_patches(parms);
end;

source_corr_arr = cell(parms.niters,1);
sensor_corr_arr = cell(parms.niters,1);
for iter=1:parms.niters
  if parms.niters>1 && parms.verbose
    fprintf('%s: iteration %d...\n',mfilename,iter);
  end;
  fname_out = set_corr_mat_fname(parms,iter);
  if ~exist(fname_out,'file') || parms.forceflag
    % simulate source waveforms
    [stc,parms] = sim_sources(parms,iter);

    % create source plots and images
    plot_sources(stc,parms,iter);

    % synthesize sensor waveforms
    [synth_data,parms] = synth_sensors(stc,parms,iter);

    % plot sensor waveforms
    plot_sensors(synth_data,parms,iter);

    % calculate measures of correlation across sources
    source_corr = calc_source_metrics(stc,parms);

    % calculate measures of correlation across sensors
    sensor_corr = calc_sensor_metrics(synth_data,parms);

    save(fname_out,'source_corr','sensor_corr');
  else
    load(fname_out);
  end;
  source_corr_arr{iter} = source_corr;
  sensor_corr_arr{iter} = sensor_corr;
end;

% write results to csv file
write_results_csv(source_corr_arr,sensor_corr_arr,parms);

G_norm = parms.G_norm;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fspath,options)
  parms = mmil_args2parms(options,{...
    'fspath',fspath,[],...
    ... % input dSPM
    'rootdir',pwd,[],...
    'prefix','dSPM',[],...
    'prefix_data',[],[],...
    'G_norm',[],[],...
    ... % output
    'outdir',pwd,[],...
    'outstem','spindles',[],...
    'source_mat_flag',true,[false true],...
    'save_data_flag',true,[false true],...
    'niters',1,[1,1e4],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
    ... % time window
    'time_dur',1000,[1e2,1e10],...
    'sf',500,[1e2,1e4],...
    ... % spindle characteristics
    'uniform_flag',false,[false true],...
    'indy_flag',true,[false true],...
    'npatches',0,[0,1e6],...
    'uniform_patches_flag',false,[false true],...
    'npatches_pool',0,[],...
    'patch_radius',10,[0,1e4],...
    'freq_mean',13,[0,1e2],...
    'freq_std',1,[0,1e2],...
    'phase_mean',0,[0,1],...
    'phase_std',1,[0,1],... 
    'phase_gauss_flag',true,[false true],...
    'amp_mean',1,[-100,100],...
    'amp_std',0,[0,100],...
    'amp_uniform',0.1,[-100,100],...
    't0_min',100,[-1e4,1e4],...
    't0_mean',250,[-1e4,1e4],...
    't0_std',50,[0,1e4],...
    'td_min',200,[-1e4,1e4],...
    'td_mean',500,[-1e4,1e4],...
    'td_std',50,[0,1e4],...
    'source_noise',0.1,[0,1e100],...
    ... % sensor data    
    'sensor_noise',0,[0,1e100],...
    'grad_noise',0,[0,1e100],...
    'mag_noise',0,[0,1e100],...
    'EEG_noise',0,[0,1e100],...
    'bandpass_flag',true,[false,true],...
    'bandpass_low_cf',8,[0,Inf],...
    'bandpass_low_tb',0.1,[0,Inf],...
    'bandpass_high_cf',18,[0,Inf],...
    'bandpass_high_tb',0.1,[0,Inf],...
    'notch_flag',false,[false,true],...
    'notch_cf',0,[0,Inf],...
    'notch_tb',0,[0,Inf],...  
    ... % plotting
    'source_plot_flag',true,[false true],...
    'source_plot_ylim',[-2,2],[],...
    'source_image_flag',true,[false true],...
    'source_image_ylim',[-1,1],[],...
    'sensor_plot_flag',true,[false true],...
    'sensor_plot_grad_ylim',[-200,200],[],...
    'sensor_plot_mag_ylim',[-400,400],[],...
    'sensor_plot_EEG_ylim',[-30,30],[],...
    'sensor_image_flag',true,[false true],...
    'sensor_image_grad_ylim',[-400,400],[],...
    'sensor_image_mag_ylim',[-800,800],[],...
    'sensor_image_EEG_ylim',[-60,60],[],...
    'sensor_topo_flag',true,[false true],...
    'grad_scalefact',10^13,[-Inf Inf],...
    'mag_scalefact',10^15,[-Inf Inf],...
    'EEG_scalefact',10^6,[-Inf Inf],...
    'sensor_types',{'grad','mag','EEG'},[],...
    'sensor_corr_image_flag',true,[false true],...
    ... % correlation calculation
    'corr_t0',100,[0,Inf],...
    'corr_t1',1000,[0,Inf],...
    'corr_source_subsamp',500,[10,10000],...
    ... % undocumented
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'surfname','white',[],...
  });

  parms.nhemi = length(parms.hemilist);
  parms.ntpoints = round(parms.time_dur * parms.sf / 1000);
  parms.time = 1000*[0:parms.ntpoints-1]/parms.sf;
  parms.corr_t1 = max(parms.corr_t0,parms.corr_t1);
  [tmp,parms.corr_s0] = min(abs(parms.time-parms.corr_t0));
  [tmp,parms.corr_s1] = min(abs(parms.time-parms.corr_t1));
  fname_parms = sprintf('%s/matfiles/%s_parms.mat',...
    parms.rootdir,parms.prefix);
  dSPM = load(fname_parms);
  parms.grad_chans = dSPM.parms.grad_chans;
  parms.mag_chans = dSPM.parms.mag_chans;
  parms.EEG_chans = dSPM.parms.EEG_chans;
  parms.nverts_lh = [];
  parms.nverts_rh = [];
  parms.v_ctx_lh = [];
  parms.v_ctx_rh = [];
  if parms.npatches_pool < parms.npatches
    parms.npatches_pool = parms.npatches;
  end;
  % divide npatches across hemispheres
  parms.npatches_pool = ceil(parms.npatches_pool/parms.nhemi);

  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_exist = check_output(parms)
  all_exist = 1;
  if parms.forceflag
    all_exist = 0;
    return;
  end;
  for iter=1:parms.niters
    if parms.save_data_flag
      % check source files
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fname_out = set_source_fname(parms,hemi,iter);
        if ~exist(fname_out,'file'), all_exist = 0; return; end;
      end;
      % check sensor files
      fname_out = set_sensor_fname(parms,iter);
      if ~exist(fname_out,'file'), all_exist = 0; return; end;
    end;
    % check plots
    if parms.source_plot_flag
      fname_out = set_source_plot_fname(parms,iter);
      if ~exist(fname_out,'file'), all_exist = 0; return; end;
    end;
    if parms.source_image_flag
      fname_out = set_source_image_fname(parms,iter);
      if ~exist(fname_out,'file'), all_exist = 0; return; end;
    end;
    for t=1:length(parms.sensor_types)
      sensor_type = parms.sensor_types{t};
      if parms.sensor_plot_flag
        fname_out = set_sensor_plot_fname(parms,sensor_type,iter);
        if ~exist(fname_out,'file'), all_exist = 0; return; end;
      end;
      if parms.sensor_image_flag  
        fname_out = set_sensor_image_fname(parms,sensor_type,iter);
        if ~exist(fname_out,'file'), all_exist = 0; return; end;
      end;
      if parms.sensor_topo_flag
        fname_out = set_sensor_topo_fname(parms,sensor_type,iter);
        if ~exist(fname_out,'file'), all_exist = 0; return; end;
      end;
      if parms.sensor_corr_image_flag
        fname_out = set_sensor_corr_image_fname(parms,sensor_type,iter);
        if ~exist(fname_out,'file'), all_exist = 0; return; end;
      end;
    end;
    % check corr mat file
    fname_out = set_corr_mat_fname(parms,iter);
    if ~exist(fname_out,'file'), all_exist = 0; return; end;
  end;
  % check corr csv file
  fname_out = set_corr_csv_fname(parms);
  if ~exist(fname_out,'file'), all_exist = 0; return; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_patches(parms)
  parms.fname_surfs = ...
    sprintf('%s/%s_surfs.mat',parms.outdir,parms.outstem);  
  parms.fname_patches = ...
    sprintf('%s/%s_patches.mat',parms.outdir,parms.outstem);
  if exist(parms.fname_surfs,'file') &&...
     exist(parms.fname_patches,'file') && ~parms.forceflag
    return;
  end;
  % load surfaces
  if ~exist(parms.fname_surfs,'file') || parms.forceflag
    clear surfs;
    for h=1:parms.nhemi;
      hemi = parms.hemilist{h};
      % load surface file
      fname_surf = [parms.fspath '/surf/' hemi '.' parms.surfname];
      if parms.verbose
        fprintf('%s: loading surface file %s...\n',mfilename,fname_surf);
        tic;
      end;
      surf = fs_read_surf(fname_surf);
      surf = fs_find_neighbors(surf,parms.verbose);
      if parms.verbose, toc; end;
      % load cortex label file
      fname_label = [parms.fspath '/label/' hemi '.cortex.label'];
      if parms.verbose
        fprintf('%s: loading label file %s...\n',mfilename,fname_label);
        tic;
      end;
      surf.v_ctx = fs_read_label(fname_label);
      if parms.verbose, toc; end;
      surfs(h) = surf;
    end;
    save(parms.fname_surfs,'surfs');  
  else
    load(parms.fname_surfs);
  end;
  for h=1:parms.nhemi;
    hemi = parms.hemilist{h};
    parms.(['v_ctx_' hemi]) = surfs(h).v_ctx;
  end;
  % create patches
  if ~exist(parms.fname_patches,'file') || parms.forceflag
    patches = cell(1,2);
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      surf = surfs(h);
      % select a subset of evenly distributed vertices
      if parms.verbose
        fprintf('%s: selecting %s vertices...\n',mfilename,hemi);
        tic;
      end;
      r = parms.npatches_pool / length(surf.v_ctx);
      surf_red = reducepatch(surf,r);
      % find original vertex numbers
      v_red = zeros(size(surf_red.vertices,1),1);
      for i=1:length(v_red)
        d = sqrt(sum(bsxfun(@minus,...
                            surf.vertices,surf_red.vertices(i,:)).^2,2));
        [mind,ind] = min(d);
        v_red(i) = ind;
      end;
      % exclude non-cortical vertices
      v_ctx_red = intersect(v_red,surf.v_ctx);
      % randomly select extra vertices to exclude
      if length(v_ctx_red) > parms.npatches_pool
        ind = randperm(length(v_ctx_red));
        v_ctx_red = v_ctx_red(ind(1:parms.npatches_pool));    
      end;
      if parms.verbose, toc; end;

      % set values
      if parms.verbose
        fprintf('%s: creating %d %s patches...\n',...
          mfilename,parms.npatches_pool,hemi);
        tic;
      end;
      vals = zeros(surf.nverts,parms.npatches_pool);
      for p=1:parms.npatches_pool
        v = v_ctx_red(p);
        tvals = zeros(surf.nverts,1);
        tvals(v) = 1;
        % dilate patches by sparse smoothing
        if parms.patch_radius>0
          tvals = fs_smooth_sparse(surf,tvals,parms.patch_radius);
        end;
        vals(:,p) = tvals;
      end;
      if parms.verbose, toc; end;
      patches{h} = vals;
    end;
    save(parms.fname_patches,'patches');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outstem = set_outstem(parms,iter)
  if parms.niters==1
    outstem = parms.outstem;
  else
    outstem = sprintf('%s_iter%d',parms.outstem,iter);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_source_fname(parms,hemi,iter)
  if parms.source_mat_flag
    fname_out = sprintf('%s/%s_sources-%s.mat',...
      parms.outdir,set_outstem(parms,iter),hemi);
  else
    fname_out = sprintf('%s/%s_sources-%s.mgz',...
      parms.outdir,set_outstem(parms,iter),hemi);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_sensor_fname(parms,iter)
  fname_out = sprintf('%s/%s_synth_sensors.mat',...
    parms.outdir,set_outstem(parms,iter));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_source_plot_fname(parms,iter)
  fname_out = sprintf('%s/%s_source_plot.tif',...
    parms.outdir,set_outstem(parms,iter));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_source_image_fname(parms,iter)
  fname_out = sprintf('%s/%s_source_image.tif',...
    parms.outdir,set_outstem(parms,iter));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_sensor_plot_fname(parms,sensor_type,iter)
  fname_out = sprintf('%s/%s_%s_sensor_plot.tif',...
    parms.outdir,set_outstem(parms,iter),sensor_type);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_sensor_image_fname(parms,sensor_type,iter)
  fname_out = sprintf('%s/%s_%s_sensor_image.tif',...
    parms.outdir,set_outstem(parms,iter),sensor_type);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_sensor_topo_fname(parms,sensor_type,iter)
  fname_out = sprintf('%s/%s_%s_sensor_topo.tif',...
    parms.outdir,set_outstem(parms,iter),sensor_type);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_sensor_corr_image_fname(parms,sensor_type,iter)
  fname_out = sprintf('%s/%s_%s_sensor_corr_image.tif',...
    parms.outdir,set_outstem(parms,iter),sensor_type);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_corr_mat_fname(parms,iter)
  fname_out = sprintf('%s/%s_corr.mat',...
    parms.outdir,set_outstem(parms,iter));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_corr_csv_fname(parms)
  fname_out = sprintf('%s/%s_corr.csv',parms.outdir,parms.outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stc,parms] = sim_sources(parms,iter)
  stc = [];
  % simulate spindle
  if parms.uniform_flag
    uniform_spindle = sim_spindle(parms,1);
  end;
  % load patches
  if parms.npatches
    if parms.uniform_patches_flag
      uniform_patches_spindle = sim_spindle(parms);
    end;
    load(parms.fname_patches);
    % set number of patches for each hemisphere
    npatches = ones(parms.nhemi,1)*floor(parms.npatches/2);
    % if npatches is odd, randomly chose which hemisphere has the extra one
    if sum(npatches)<parms.npatches 
      extra = round(rand);
      npatches(1) = npatches(1) + extra;
      npatches(2) = npatches(2) + (1-extra);
    end;
  else
    npatches = zeros(parms.nhemi,1);
  end;
  % set values for each vertex
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    npatches_hemi = npatches(h);
    fname_out = set_source_fname(parms,hemi,iter);
    if ~exist(fname_out,'file') || parms.forceflag
      if isempty(parms.(['nverts_' hemi]))
        % read number of vertices
        fname_surf = [parms.fspath '/surf/' hemi '.' parms.surfname];
        nverts = fs_read_surf_nverts(fname_surf);
        parms.(['nverts_' hemi]) = nverts;
      else
        nverts = parms.(['nverts_' hemi]);
      end;
      if isempty(parms.(['v_ctx_' hemi]))
        % load cortex label file
        fname_label = [parms.fspath '/label/' hemi '.cortex.label'];
        if parms.verbose
          fprintf('%s: loading label file %s...\n',mfilename,fname_label);
          tic;
        end;
        v_ctx = fs_read_label(fname_label);
        if parms.verbose, toc; end;
        parms.(['v_ctx_' hemi]) = v_ctx;
      else
        v_ctx = parms.(['v_ctx_' hemi]);
      end;
      if npatches_hemi
        % select patches from pool
        tmp_patches = patches{h};
        if npatches_hemi && parms.npatches_pool > npatches_hemi
          p_order = randperm(parms.npatches_pool);
          tmp_patches = tmp_patches(:,p_order(1:npatches_hemi));
        end;
        % normalize each vertex to sum across patches
        %   so that maximum amplitude is constant (in case of overlap)
        tmp_patches = bsxfun(@rdivide,tmp_patches,max(1,sum(tmp_patches,2)));
      end;
      nverts_ctx = length(v_ctx);
      % set values
      vals = zeros(nverts,parms.ntpoints);
      if parms.uniform_flag
        vals(v_ctx,:) = repmat(uniform_spindle,[nverts_ctx,1]);
      end;
      if parms.indy_flag
        if parms.npatches
          if npatches_hemi
            if parms.uniform_patches_flag
              v_patches = find(sum(tmp_patches,2));
              vals(v_patches,:) = vals(v_patches,:) + ...
                                  repmat(uniform_patches_spindle,[length(v_patches),1]);
            else
              if parms.verbose
                fprintf('%s: simulating %d spindles for %s...\n',...
                  mfilename,npatches_hemi,hemi);
                tic
              end;
              for k=1:npatches_hemi
                v = find(tmp_patches(:,k));
                spindle = sim_spindle(parms);
                vals = vals + tmp_patches(:,k) * spindle;
              end;
            end;
          end;
        else
          if parms.verbose
            fprintf('%s: simulating %d spindles for %s...\n',...
              mfilename,nverts_ctx,hemi);
            tic
          end;
          for k=1:nverts_ctx
            v = v_ctx(k);
            spindle = sim_spindle(parms);
            vals(v,:) = vals(v,:) + spindle;
          end;
        end;
        if parms.verbose, toc; end;
      end;
      % add noise and save output file
      if parms.source_noise>0
        vals = vals + parms.source_noise*randn(size(vals));
      end;
      if parms.save_data_flag
        % save output
        if parms.verbose
          fprintf('%s: saving output to %s...\n',mfilename,fname_out);
          tic;
        end;
        if parms.source_mat_flag
          save(fname_out,'vals');
        else
          fs_save_mgh(reshape(vals,[nverts,1,1,parms.ntpoints]),fname_out);
        end;
        if parms.verbose, toc; end;
      end;
    else
      if parms.verbose
        fprintf('%s: loading source time courses from %s...\n',...
          mfilename,fname_out);
        tic;
      end;
      if parms.source_mat_flag
        load(fname_out);
      else
        vals = squeeze(fs_load_mgh(fname_out));
      end;
      if parms.verbose, toc; end;
    end;
    stc = cat(1,stc,vals);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spindle = sim_spindle(parms,uniform_flag)
  if ~exist('uniform_flag','var') || isempty(uniform_flag)
    uniform_flag = 0;
  end;
  freq = parms.freq_mean + parms.freq_std*randn;
  if parms.phase_gauss_flag
    phas = parms.phase_mean + parms.phase_std*randn;
  else
    phas = parms.phase_mean + parms.phase_std*rand;
  end;
  if uniform_flag
    amp = parms.amp_uniform;
  else
    amp = parms.amp_mean + parms.amp_std*randn;
  end;
  t0 = max(parms.t0_min,parms.t0_mean + parms.t0_std*randn);
  td = max(parms.td_min,parms.td_mean + parms.td_std*randn);
  spindle = ts_sim_spindle(parms.ntpoints,freq,phas,t0,td,amp,parms.sf);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sources(stc,parms,iter)
  if parms.source_plot_flag
    fname_out = set_source_plot_fname(parms,iter);
    if ~exist(fname_out,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: creating source plot...\n',mfilename);
        tic;
      end;
      % create source plot
      figure(1); clf;
      plot(parms.time,stc);
      axis tight;
      ylim(parms.source_plot_ylim);
      xlabel('time');
      set(gcf,'visible','off');
      print(gcf,'-dtiff',fname_out);
      close(gcf);
    end;
  end;
  if parms.source_image_flag  
    fname_out = set_source_image_fname(parms,iter);
    if ~exist(fname_out,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: creating source image...\n',mfilename);
        tic;
      end;
      % create source image
      figure(1); clf;
      imagesc(stc,parms.source_image_ylim);
      axis off;
      xlabel('time');
      set(gcf,'visible','off');
      print(gcf,'-dtiff',fname_out);
      close(gcf);
      if parms.verbose, toc; end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [synth_data,parms] = synth_sensors(stc,parms,iter)
  fname_out = set_sensor_fname(parms,iter);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: synthesizing sensor data...\n',mfilename);
      tic;
    end;
    [synth_data,G_norm] = ts_synth_sensors_from_dSPM(stc,...
      'rootdir',parms.rootdir,...
      'prefix',parms.prefix,...
      'prefix_data',parms.prefix_data,...
      'time',parms.time,...
      'noisefact',parms.sensor_noise,...
      'grad_noise',parms.grad_noise,...
      'mag_noise',parms.mag_noise,...
      'EEG_noise',parms.EEG_noise,...
      'bandpass_flag',parms.bandpass_flag,...
      'bandpass_low_cf',parms.bandpass_low_cf,...
      'bandpass_low_tb',parms.bandpass_low_tb,...
      'bandpass_high_cf',parms.bandpass_high_cf,...
      'bandpass_high_tb',parms.bandpass_high_tb,...
      'notch_flag',parms.notch_flag,...
      'notch_cf',parms.notch_cf,...
      'notch_tb',parms.notch_tb,...
      'G_norm',parms.G_norm,...
      'verbose',parms.verbose);
    if parms.verbose, toc; end;
    if parms.save_data_flag
      save(fname_out,'synth_data','G_norm');
    end;
  else
    if parms.verbose
      fprintf('%s: loading synthesized sensor data from %s...\n',...
        mfilename,fname_out);
      tic;
    end;
    load(fname_out);
    if parms.verbose, toc; end;
  end;
  parms.G_norm = G_norm;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sensors(synth_data,parms,iter)
  for t=1:length(parms.sensor_types)
    sensor_type = parms.sensor_types{t};
    switch sensor_type
      case 'grad'
        plot_ylim = parms.sensor_plot_grad_ylim;
        image_ylim = parms.sensor_image_grad_ylim;
        scalefact = parms.grad_scalefact;
        chans = parms.grad_chans;
      case 'mag'
        plot_ylim = parms.sensor_plot_mag_ylim;
        image_ylim = parms.sensor_image_mag_ylim;
        scalefact = parms.mag_scalefact;
        chans = parms.mag_chans;
      case 'EEG'
        plot_ylim = parms.sensor_plot_EEG_ylim;
        image_ylim = parms.sensor_image_EEG_ylim;
        scalefact = parms.EEG_scalefact;
        chans = parms.EEG_chans;
    end;
    if isempty(chans), continue; end;
    % get data for this sensor type
    data = scalefact*synth_data.averages.data(chans,:);
    % create plot
    if parms.sensor_plot_flag
      fname_out = set_sensor_plot_fname(parms,sensor_type,iter);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: creating %s sensor plot...\n',mfilename,sensor_type);
          tic;
        end;
        % create sensor plot
        figure(1); clf;
        plot(parms.time,data);
        axis tight;
        ylim(plot_ylim);
        xlabel('time');
        set(gcf,'visible','off');
        print(gcf,'-dtiff',fname_out);
        close(gcf);
      end;
    end;
    % create image
    if parms.sensor_image_flag  
      fname_out = set_sensor_image_fname(parms,sensor_type,iter);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: creating %s sensor image...\n',mfilename,sensor_type);
          tic;
        end;
        % create sensor image
        figure(1); clf;
        imagesc(data,image_ylim);
        axis off;
        xlabel('time');
        set(gcf,'visible','off');
        print(gcf,'-dtiff',fname_out);
        close(gcf);
        if parms.verbose, toc; end;
      end;
    end;
    % create topoplot
    if parms.sensor_topo_flag
      fname_out = set_sensor_topo_fname(parms,sensor_type,iter);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: creating %s sensor topography...\n',mfilename,sensor_type);
          tic;
        end;
        figure(1); clf;
        set(gcf,'visible','off');
        % find time of max val
        maxvals = max(abs(data),[],1);
        [maxval,ind_max] = max(maxvals);
        t_max = synth_data.averages(1).time(ind_max);
        if strcmp(sensor_type,'EEG')
          % check for bad EEG sensor locations
          nchans = length(parms.EEG_chans);
          EEG_loc = zeros(nchans,3);
          for i=1:nchans
            i_chan = parms.EEG_chans(i);
            EEG_loc(i,:) = synth_data.sensor_info(i_chan).loc(1:3,4);
          end;
          d = sqrt(sum(EEG_loc.^2,2));
          ind_bad = find(d<eps);
          for i=1:length(ind_bad)
            i_chan = parms.EEG_chans(ind_bad(i));
            synth_data.sensor_info(i_chan).badchan = 1;
          end;
        end;
        % create topoplot
        ts_topoplot_avg(synth_data,'chantype',lower(sensor_type),...
          'time0',t_max,'time1',t_max);
        set(gcf,'visible','off');
        print(gcf,'-dtiff',fname_out);
        close(gcf);
        if parms.verbose, toc; end;
      end;
    end;
    % create sensor corr image
    if parms.sensor_corr_image_flag
      fname_out = set_sensor_corr_image_fname(parms,sensor_type,iter);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: creating %s sensor corr image...\n',mfilename,sensor_type);
          tic;
        end;
        % calculate correlation between sensors
        tmp_data = data(:,parms.corr_s0:parms.corr_s1);
        R = corr(tmp_data');
        % create correlation matrix image
        figure(1); clf;
        imagesc(R,[-1,1]);
        colorbar
        axis off;
        xlabel('channel');
        ylabel('channel');
        set(gcf,'visible','off');
        print(gcf,'-dtiff',fname_out);
        close(gcf);
        if parms.verbose, toc; end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source_corr = calc_source_metrics(stc,parms)
  source_corr = [];
  fprintf('%s: calculating source corr...\n',mfilename);
  nsources = size(stc,1);
  parms.corr_source_subsamp = min(nsources,parms.corr_source_subsamp);
  ind_sources = round(linspace(1,nsources,parms.corr_source_subsamp));
  data = stc(ind_sources,parms.corr_s0:parms.corr_s1);

  % correlation between sources
  R = corr(data');
  R2 = R.^2;
  % average across upper half of matrix
  [X,Y] = meshgrid(1:parms.corr_source_subsamp,1:parms.corr_source_subsamp);
  mask_upper = (X>Y);
  source_corr.r = ts_nan_mean(abs(R(mask_upper)));
  source_corr.r2 = ts_nan_mean(R2(mask_upper));

  % calculate mean and cv across sensors
  %   of standard deviation across time
  std_time = std(data,1,2);
  source_corr.ms = mean(std_time);
  source_corr.cv = std(std_time)/source_corr.ms;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sensor_corr = calc_sensor_metrics(synth_data,parms)
  sensor_corr = [];
  fprintf('%s: calculating sensor corr...\n',mfilename);
  for t=1:length(parms.sensor_types)
    sensor_type = parms.sensor_types{t};
    switch sensor_type
      case 'grad'
        plot_ylim = parms.sensor_plot_grad_ylim;
        image_ylim = parms.sensor_image_grad_ylim;
        scalefact = parms.grad_scalefact;
        chans = parms.grad_chans;
      case 'mag'
        plot_ylim = parms.sensor_plot_mag_ylim;
        image_ylim = parms.sensor_image_mag_ylim;
        scalefact = parms.mag_scalefact;
        chans = parms.mag_chans;
      case 'EEG'
        plot_ylim = parms.sensor_plot_EEG_ylim;
        image_ylim = parms.sensor_image_EEG_ylim;
        scalefact = parms.EEG_scalefact;
        chans = parms.EEG_chans;
    end;
    if isempty(chans), continue; end;
    % get data for this sensor type
    data = scalefact*synth_data.averages.data(chans,parms.corr_s0:parms.corr_s1);
    nchans = size(data,1);
    % calculate max sensor magnitude
    sensor_corr(t).maxval = max(abs(data(:)));
    % correlation between sensors
    R = corr(data');
    R2 = R.^2;
    % average across upper half of matrix
    [X,Y] = meshgrid(1:nchans,1:nchans);
    mask_upper = (X>Y);
    sensor_corr(t).r = mean(abs(R(mask_upper)));
    sensor_corr(t).r2 = mean(R2(mask_upper));
    % calculate mean and cv across sensors
    %   of standard deviation across time
    std_time = std(data,1,2);
    sensor_corr(t).ms = mean(std_time);
    sensor_corr(t).cv = std(std_time)/sensor_corr(t).ms;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_results_csv(source_corr,sensor_corr,parms)
  fname_out = set_corr_csv_fname(parms);
  if ~exist(fname_out,'file') || parms.forceflag
    % set column labels
    col_labels = {};    
    j = 1;
    col_labels{j} = 'iter'; j = j + 1;
    fnames = fieldnames(source_corr{1});
    for f=1:length(fnames)
      col_labels{j} = ['source_' fnames{f}];
      j = j + 1;
    end;
    for t=1:length(parms.sensor_types)
      sensor_type = parms.sensor_types{t};
      fnames = fieldnames(sensor_corr{1});
      for f=1:length(fnames)
        col_labels{j} = [sensor_type '_' fnames{f}];
        j = j + 1;
      end;
    end;
    % set output values
    output = cell(parms.niters,length(col_labels));
    for iter=1:parms.niters
      j = 1;
      output{iter,j} = iter; j = j + 1;
      fnames = fieldnames(source_corr{1});
      for f=1:length(fnames)
        col_labels{j} = ['source_' fnames{f}];
        output{iter,j} = source_corr{iter}.(fnames{f});
        j = j + 1;
      end;
      for t=1:length(parms.sensor_types)
        sensor_type = parms.sensor_types{t};
        fnames = fieldnames(sensor_corr{1});
        for f=1:length(fnames)
          col_labels{j} = [sensor_type '_' fnames{f}];
          output{iter,j} = sensor_corr{iter}(t).(fnames{f});
          j = j + 1;
        end;
      end;
    end;
    mmil_write_csv(fname_out,output,'col_labels',col_labels);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

