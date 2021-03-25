function ts_sim_dSPM(varargin)
%function ts_sim_dSPM([options])
%
% NOTE: this function relies on the output of ts_dSPM
%
% Usage:
%  ts_sim_dSPM('key1', value1,...);
%
% Optional Inupt:
%   'prefix': prefix of ts_dSPM output files
%     {default = 'dSPM'}
%   'prefix_data': prefix of processed data file (e.g. 'proc_avg_data')
%     if not supplied, will use avg_data saved by ts_dSPM
%     {default = []}
%   'prefix_forward': prefix of dSPM output files for forward solution
%     if empty, use prefix
%     {default = []}
%   'resp_time': latency of modeled evoked response
%     {default = 100}
%   'time': vector of time points (msec)
%     if empty, will use time vector in avg_data.averages(1).time
%     {default = []}
%   'stc': source time course matrix
%     if not empty, size must be [num_sources,num_tpoints]
%       and other source parameters will be ignored
%     {default = []}
%   'source_roiname': name of label file roi to be source
%     Can be a cell array of multiple ROIs
%     {default = []}
%   'source_vert': vertex number for source
%     Can be a vector of vertex numbers
%     Ignored if source_roiname is not empty
%     {default = 1}
%   'source_heminum': 1 = 'lh', 2 = 'rh'
%     Can be a vector to match number of sources
%       (in source_roiname or source_vert) otherwise, use same for all sources
%     {default = 1}
%   'source_dur': duration of evoked response for each source
%     Can be a vector to match number of sources
%     {default = 100}
%   'source_dt': source-specific delta time added to resp_time
%     Can be a vector to match number of sources
%     {default = 0}
%   'sourcefact': amplitude of modeled source
%     {default = 1}
%   'noisefact': stdev of noise added to synth_data,
%      relative to stdev of noise in avg_data
%     {default = 0.05}
%   'SNR': assumed signal to noise ratio (for regularization)
%     {default = 3}
%   'rootdir': directory containing matfiles dir
%     {default = pwd}
%   'outdir': output directory
%     {default = pwd/sim_dSPM}
%   'outstem': output file stem
%     {default = sim}
%   'source_infix': output file infix for sources
%     if not supplied, will be a concatenated list of source names
%     {default = []}
%   'fif_flag': [0|1] whether to save synthesized sensor data to fif file
%     {default = 0}
%   'template_fif': template fif file (e.g. raw or online avg)
%     Required for saving synthesized data to fif file
%     {default = []}
%   'roisubj': freesurfer subject with labels subdir containing ROI label files
%     {default = 'fsaverage'}
%   'roisubjdir': freesurfer subjects dir containing roisubj
%     {default = $SUBJECTS_DIR}
%   'roinames': cell array of ROI names
%      Label files have this format: hemi.roiname.label or 
%                                    hemi.roiname.ico4.label
%      If empty, will search for label files in roisubjdir/roisubj/label
%      If none found, will skip ROI analysis
%      {default = []}
%   'roi_hemilist': cell array of hemispheres for ROIs
%      {default = {'lh','rh'}}
%   'ico': icosahedral order number
%     If 0, assume native subject space
%     {default = 0}
%
% Created:  01/06/10 by Don Hagler
% Last Mod: 11/26/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix','dSPM',[],...
  'prefix_data',[],[],...
  'prefix_forward',[],[],...
  'resp_time',100,[0,Inf],...
  'time',[],[],...
  'stc',[],[],...
  'source_roiname',[],[],...
  'source_vert',1,[1,Inf],...
  'source_heminum',1,[1,2],...
  'source_dur',100,[0,Inf],...
  'source_dt',0,[-Inf,Inf],...
  'sourcefact',1,[0,Inf],...
  'noisefact',0.05,[0,Inf],...
  'SNR',3,[eps,Inf],...
  'rootdir',pwd,[],...
  'outdir',[pwd '/sim_dSPM'],[],...
  'outstem','sim',[],...
  'source_infix',[],[],...
  'fif_flag',false,[false true],...
  'template_fif',[],[],...
  'roisubj','fsaverage',[],...
  'roisubjdir',getenv('SUBJECTS_DIR'),[],...
  'roinames',[],[],...
  'roi_hemilist',{'lh','rh'},{'lh','rh'},...
  'ico',0,[0:7],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'ncov_type',2,[1:3],...
  'noisenorm_flag',true,[false true],...
  'noisenorm_identity_flag',true,[false true],...
  'depthweight_flag',false,[false true],...
  'depthweight_p',0.5,[0,1],...
  'bandpass_flag',false,[false,true],...
  'bandpass_low_cf',0,[0,Inf],...
  'bandpass_low_tb',0,[0,Inf],...
  'bandpass_high_cf',100,[0,Inf],...
  'bandpass_high_tb',0,[0,Inf],...
  'notch_flag',false,[false,true],...
  'notch_cf',0,[0,Inf],...
  'notch_tb',0,[0,Inf],...  
  'write_mgh_flag',false,[false true],...
  'overwrite_forward_flag',false,[false true],...
  'overwrite_inverse_flag',false,[false true],...
  'overwrite_output_flag',true,[false true],...
  'save_avg_flag',true,[false true],...
  'plotflag',true,[false true],...
  'plot_type','tiff',{'tiff','jpeg','epsc'},...
  'plot_xlim',[],[],...
  'plot_ylim',[],[],...
  'plot_offset',1,[-Inf,Inf],...
  'visible_flag',false,[false true],...
...
  'dSPM_tags',{'subjdir','rootoutdir','prefix','conditions',...
               'calc_dipinfo_flag','lh_dip_info','rh_dip_info',...
               'lh_dip_file','rh_dip_file','nodec_flag','lh_dec_dips',...
               'rh_dec_dips','lh_dec_file','rh_dec_file','ncov_type',...
               'ncov_conditions','raw_ncov_lambda','calc_scalefacts_flag',...
               'baseline_start','baseline_end','baseline_flag',...
               'ssp_projmat','forceflag','SNR','noisenorm_flag',...
               'noisenorm_identity_flag','depthweight_flag',...
               'depthweight_p','bem_flag','openmeeg_flag',...
               'radii','conductivities',...
               'EEG_gain_scalefact','nlayers','badchans','badchanfile',...
               'usegrad_flag','usemag_flag','useEEG_flag','grad_scalefact',...
               'mag_scalefact','EEG_scalefact','write_stc_flag',...
               'stc_scalefact','write_mgh_flag','sparsesmooth',...
               'postsmooth','mbmask_flag','resamp2ico_flag','icolevel',...
               'icosmooth','write_fif_flag','template_fif','bem_surf_files',...
               'cen_sph','trans','transfile','alignment_fif',...
               'lh_sourcecov_file','rh_sourcecov_file','sourcecov_thresh',...
               'sourcecov_thresh_abs_flag','sourcecov_maxvar',...
               'sourcecov_minvar','forward_matfile','inverse_matfile',...
               'refEEG_coords','forward_only_flag','prewhiten_flag',...
               'orient_constr_flag',...
               'orient_tang','smooth_constr_flag','smooth_constr_nsmooth',...
               'signed_sources_flag','datatype','save_results_flag',...
               'save_avg_flag','save_fiterr_flag','save_fiterr_struct_flag',...
               'hemilist'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cwd = pwd;

parms.roidir = [parms.roisubjdir '/' parms.roisubj '/label'];
if ~exist(parms.roidir,'dir')
  fprintf('%s: WARNING: roidir %s not found, skipping roi analysis\n',...
    mfilename,parms.roidir);
  parms.roidir = [];
end;

if isempty(parms.roinames)
  if ~isempty(parms.roidir)
    flist = dir(sprintf('%s/*.label',parms.roidir));
    if isempty(flist)
      fprintf('%s: WARNING: no label files found in roidir %s, skipping roi analysis\n',...
        mfilename,parms.roidir);
    else
      tmp_roinames = {flist.name};
      if parms.ico>0
        pat = sprintf('^[lr]h.(?<roiname>.+).ico%d.label$',parms.ico);
      else
        pat = '^[lr]h.(?<roiname>.+).label$';
      end;
      n = regexp(tmp_roinames,pat,'names');
      keep_flags = ~cellfun(@isempty,n);
      ind_keep = find(keep_flags);
      for k=1:length(ind_keep)
        parms.roinames{end+1} = n{ind_keep(k)}.roiname;
      end;
      parms.roinames = unique(parms.roinames);
    end;
  end;
end;

if ~isempty(parms.roinames)
  for h=1:length(parms.roi_hemilist)
    hemi = parms.roi_hemilist{h};
    % check that rois exist
    tmp_roifiles = [];
    for r=1:length(parms.roinames)
      if parms.ico>0
        fname = sprintf('%s/%s.%s.ico%d.label',...
          parms.roidir,hemi,parms.roinames{r},parms.ico);
      else
        fname = sprintf('%s/%s.%s.ico%d.label',...
          parms.roidir,hemi,parms.roinames{r},parms.ico);
      end;
      if ~exist(fname,'file')
        fprintf('%s: WARNING: label file %s not found, skipping\n',...
          mfilename,fname);
      else
        tmp_roifiles{end+1} = fname;
      end;
    end;
    switch hemi
      case 'lh'
        parms.lh_roifiles = tmp_roifiles;
      case 'rh'
        parms.rh_roifiles = tmp_roifiles;
    end;
  end;
end;

if parms.plotflag
  switch parms.plot_type
    case 'tiff'
      parms.plotext = 'tif';
    case 'epsc'
      parms.plotext = 'eps';
    case 'jpeg'
      parms.plotext = 'jpg';
  end;
end;

% check source roi files
if ~isempty(parms.source_roiname)
  % check if source_roiname is a cell
  if ~iscell(parms.source_roiname)
    parms.source_roiname = {parms.source_roiname};
  end;
  parms.num_sources = length(parms.source_roiname);
  parms.source_roifiles = cell(parms.num_sources,1);
  for r=1:parms.num_sources
    hemi = parms.hemilist{parms.source_heminum(r)};
    if parms.ico>0
      fname = sprintf('%s/%s.%s.ico%d.label',...
        parms.roidir,hemi,parms.source_roiname{r},parms.ico);
    else
      fname = sprintf('%s/%s.%s.ico%d.label',...
        parms.roidir,hemi,parms.source_roiname{r},parms.ico);
    end;
    if ~exist(fname,'file')
      error('source ROI file %s not found',fname);
    end;
    parms.source_roifiles{r} = fname;
  end;
else
  parms.num_sources = length(parms.source_vert);
end;

% check that source_dur, source_dt match
if length(parms.source_dur)==1
  parms.source_dur = parms.source_dur*ones(parms.num_sources,1);
elseif length(parms.source_dur)~=parms.num_sources
  error('length of source_dur (%d) does not match number of sources (%d)',...
    length(parms.source_dur),parms.num_sources);
end;
if length(parms.source_dt)==1
  parms.source_dt = parms.source_dt*ones(parms.num_sources,1);
elseif length(parms.source_dt)~=parms.num_sources
  error('length of source_dt (%d) does not match number of sources (%d)',...
    length(parms.source_dt),parms.num_sources);
end;
if length(parms.source_heminum)==1
  parms.source_heminum = parms.source_heminum*ones(parms.num_sources,1);
elseif length(parms.source_heminum)~=parms.num_sources
  error('length of source_heminum (%d) does not match number of sources (%d)',...
    length(parms.source_heminum),parms.num_sources);
end;

parms.prefix_orig = parms.prefix;
if isempty(parms.prefix_forward)
  parms.prefix_forward = parms.prefix_orig;
end;

if isempty(parms.source_infix)
  sourcelist = [];
  for r=1:parms.num_sources
    hemi = parms.hemilist{parms.source_heminum(r)};
    if isempty(parms.source_roiname)
      roiname = sprintf('v%d',parms.source_vert(r));
    else
      roiname = parms.source_roiname{r};
    end;
    if ~isempty(sourcelist)
      sourcelist = [sourcelist '+' roiname '-' hemi];
    else
      sourcelist = [roiname '-' hemi];
    end;
  end;
  parms.source_infix = sourcelist;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.outdir,'dir'), mkdir(parms.outdir); end;

% load parms
fname = sprintf('%s/matfiles/%s_parms.mat',parms.rootdir,parms.prefix_orig);
if ~exist(fname,'file'), error('file %s not found',fname); end;
tmp_parms = parms;
load(fname);

if mmil_isrelative(parms.forward_matfile)
  parms.forward_matfile = sprintf('%s/matfiles/%s',...
    tmp_parms.rootdir,parms.forward_matfile);
end;
if ~exist(parms.forward_matfile,'file')
  fprintf('%s: WARNING: forward_matfile %s not found\n',...
    mfilename,parms.forward_matfile);
  parms.forward_matfile = sprintf('%s/matfiles/%s_forward.mat',...
    tmp_parms.outdir,tmp_parms.prefix_orig);
end;

if tmp_parms.SNR~=parms.SNR
  parms.inverse_matfile = sprintf('%s/matfiles/%s_inverse.mat',...
    tmp_parms.rootdir,tmp_parms.prefix_orig);
else
  parms.inverse_matfile = [];
end;
fields = fieldnames(tmp_parms);
for f=1:length(fields)
  parms = setfield(parms,fields{f},getfield(tmp_parms,fields{f}));
end;

% load avg_data
if isempty(parms.prefix_data)
  fname = sprintf('%s/matfiles/%s_avg_data.mat',parms.rootdir,parms.prefix_orig);
else
  fname = sprintf('%s/matfiles/%s.mat',parms.rootdir,parms.prefix_data);
end;
if ~exist(fname,'file'), error('file %s not found',fname); end;
load(fname);

if isempty(parms.time)
  parms.time = avg_data.averages(1).time;
end;

% create simulated source waveforms
if isempty(parms.stc)
  fprintf('%s: creating simulated source waveforms...\n',mfilename);
  lh_stc = zeros(parms.num_dec_dips_lh,length(parms.time));
  rh_stc = zeros(parms.num_dec_dips_rh,length(parms.time));
  for r=1:parms.num_sources
    hemi = parms.hemilist{parms.source_heminum(r)};
    if isempty(parms.source_roiname)
      roi = parms.source_vert(r);
      roiname = sprintf('v%d',roi);
    else
      roiname = parms.source_roiname{r};
      fname_roi = parms.source_roifiles{r};
      roi = fs_read_label(fname_roi);
    end;
    % create source time course matrix
    time_on = parms.resp_time + parms.source_dt(r);
    time_off = time_on + parms.source_dur(r);
    [tmp,ind_on] = min(abs(parms.time-time_on/1000));
    [tmp,ind_off] = min(abs(parms.time-time_off/1000));
    source_wf = zeros(1,length(parms.time));
    source_wf(ind_on:ind_off) = parms.sourcefact;
    switch hemi
      case 'lh'
        lh_stc(roi,:) = ones(length(roi),1)*source_wf;
      case 'rh'
        rh_stc(roi,:) = ones(length(roi),1)*source_wf;
    end;
  end;
  parms.stc = [lh_stc;rh_stc];
else
  [num_dips,ntpoints] = size(parms.stc);
  if num_dips ~= parms.num_sources
    error('stc has wrong number of dipoles (%d instead of %d)\n',...
      num_dips,parms.num_sources);
  end;
  if ntpoints ~= length(parms.time)
    error('stc has wrong number of time points (%d instead of %d)\n',...
      ntpoints,length(parms.time));
  end;
end;

% simulate source waveforms
fprintf('%s: synthesizing sensor data...\n',mfilename);
%% todo: simulate conditions other than 1st?
synth_data = ts_synth_sensors_from_dSPM(parms.stc,...
  'prefix',parms.prefix_forward,...
  'prefix_data',parms.prefix_data,...
  'time',parms.time,...
  'rootdir',parms.rootdir,...
  'noisefact',parms.noisefact,...
  'bandpass_flag',parms.bandpass_flag,...
  'bandpass_low_cf',parms.bandpass_low_cf,...
  'bandpass_low_tb',parms.bandpass_low_tb,...
  'bandpass_high_cf',parms.bandpass_high_cf,...
  'bandpass_high_tb',parms.bandpass_high_tb,...
  'notch_flag',parms.notch_flag,...
  'notch_cf',parms.notch_cf,...
  'notch_tb',parms.notch_tb);
parms.conditions = 1;
parms.ncov_conditions = 1;

if parms.fif_flag
  fif_outdir = [parms.outdir '/fifs'];
  mmil_mkdir(fif_outdir);
  tmp_outstem = sprintf('%s/%s_%s',...
    fif_outdir,parms.outstem,parms.source_infix);
  ts_avg2fif(synth_data,parms.template_fif,tmp_outstem);
end;

% run dSPM on synth_data
fprintf('%s: running dSPM on synthesized sensor data...\n',mfilename);
parms.prefix = sprintf('%s_%s_%s',...
  parms.outstem,parms.source_infix,parms.prefix_orig);
parms.rootoutdir = parms.outdir;
dSPM_args = mmil_parms2args(parms,parms.dSPM_tags);
cd(parms.rootdir)
ts_dSPM(synth_data,parms.subjname,dSPM_args{:});
cd(cwd);
% todo: just apply inverse and get stc directly?
%       calculating inverse and file i/o are slow

if ~isempty(parms.roinames)
  % load stc files and extract roi waveforms
  for h=1:length(parms.roi_hemilist)
    hemi = parms.roi_hemilist{h};

    % create stc file name
    stcfile = sprintf('%s/stcfiles/%s_cond%02d-%s.stc',...
      parms.rootoutdir,parms.prefix,1,hemi);

    % create label list
    switch hemi
      case 'lh'
        labelfiles = parms.lh_roifiles;
      case 'rh'
        labelfiles = parms.rh_roifiles;
    end;

    % extract waveforms for ROIs
    matfile = sprintf('%s/matfiles/%s-roi-stc-%s.mat',...
      parms.rootoutdir,parms.prefix,hemi);
    fprintf('%s: extracting waveforms for hemi %s...\n',mfilename,hemi);
    [time,wforms] = ts_stc2roi(stcfile,labelfiles);
    if ~isempty(wforms)
      save(matfile,'time','wforms','stcfile','labelfiles');
    end;

    if parms.plotflag
      fprintf('%s: plotting waveforms for hemi %s...\n',mfilename,hemi);
      plotdir = [parms.outdir '/plots'];
      mmil_mkdir(plotdir);
      figure;
      for r=1:length(parms.roinames)
        wform = squeeze(wforms(r,:))-parms.plot_offset;
        clf; plot(time,wform);
        if ~isempty(parms.plot_xlim)
          set(gca,'XLim',parms.plot_xlim);
        end;
        if ~isempty(parms.plot_ylim)
          set(gca,'YLim',parms.plot_ylim);
        end;
        if ~parms.visible_flag
          set(gcf,'Visible','Off');
        end;
        fname = sprintf('%s/%s_%s_roi_stc_%s-%s.%s',...
          plotdir,parms.outstem,parms.source_infix,...
          parms.roinames{r},hemi,parms.plotext);
        title([hemi ' ' parms.roinames{r}]);
        print(['-d' parms.plot_type],fname);
        if ~parms.visible_flag
          close(gcf);
        end;
      end;
    end;
  end;
end;
