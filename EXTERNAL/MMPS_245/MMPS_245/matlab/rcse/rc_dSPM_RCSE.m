function rc_dSPM_RCSE(varargin)
%function rc_dSPM_RCSE([options])
%
% Purpose: Create dSPM using RCSE results
%
% Usage:
%  rc_dSPM_RCSE('key1', value1,...);
%
% Optional Inupt:
%  'rootdir': directory containing matfiles dir for dSPM
%     {default = pwd}
%  'rootoutdir': root output directory
%     {default = pwd}
%  'prefix': prefix of dSPM-RCSE output files
%     if empty, will construct from outstem, dSPM_prefix, and RCSE_prefix
%     {default = []}
%  'outstem': output file stem
%     {default = dSPM_RCSE}
%  'dSPM_prefix': prefix of ts_dSPM output files
%     {default = 'dSPM'}
%  'RCSE_prefix': prefix of RCSE output files for source waveforms
%     {default = 'RCSE'}
%  'RCSE_rootdir': directory containing matfiles dir for RCSE sources
%     If empty, rootdir
%     {default = rootdir}
%  'RCSE_forward_prefix': prefix of RCSE output files for forward solution
%     If empty, RCSE_prefix
%     {default = []}
%  'RCSE_forward_rootdir': directory containing matfiles dir for RCSE forward
%     If empty, rootdir
%     {default = []}
%  'fiterr_flag': which type of sensor data
%    0: run dSPM on RCSE fitted data plus noise
%    1: run dSPM on residual error of RCSE fit
%    2: run dSPM on RCSE fit projected into sensor data
%    3: run dSPM on data minus RCSE fit projected into sensor data
%     {default = 0}
%  'areas': when fiterr_flag=0, use only these area numbers
%    if empty, use all areas
%    {default = []}
%  'conditions': vector of condition numbers (index to avg_data.averages)
%     defining subset of stimulus locations in retmap to use
%     If empty, use all conditions in retmap.cond_order
%       or conditions in cond_info with non-zero contrast
%     {default = []}
%  'sourcefact': amplitude of modeled sources
%     {default = 1}
%  'baselineflag': [0|1] use repeated baseline as noise added to synth_data
%     {default = 1}
%  'save_synth_flag': [0|1] save synthesized sensor data to mat file
%     {default = 0}
%  'fif_flag': [0|1] save synthesized sensor data to fif file
%     {default = 0}
%  'template_fif': template fif file (e.g. raw or online avg)
%     Required for saving synthesized data to fif file
%     {default = []}
%  'write_mgh_flag': [0|1] save output maps as mgh files
%     {default = 0}
%  'ncov_type': type of noise covariance matrix
%     0 = identity matrix, 1 = average prestimulus, 2 = single trials
%     {default = 2}
%  'calc_scalefacts_flag': [0|1] calculate scaling factors
%     for different channel types (standard deviation of avg baseline)
%     if 1, calculate these scaling factors and apply to data and forward matrix
%     if 0, use default or user specified scaling factors
%     NOTE: ignored if prewhiten_flag=1
%     {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% NOTE: this function relies on having pre-run rc_RCSE and ts_dSPM
%
% Created:  03/12/10 by Don Hagler
% Last Mod: 11/01/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'rootdir',pwd,[],...
  'rootoutdir',pwd,[],...
  'prefix',[],[],...
  'outstem','dSPM_RCSE',[],...
  'dSPM_prefix','dSPM',[],...
  'RCSE_prefix','RCSE',[],...
  'RCSE_rootdir',[],[],...
  'RCSE_forward_prefix',[],[],...
  'RCSE_forward_rootdir',[],[],...
  'fiterr_flag',0,[0:3],...
  'areas',[],[],...
  'conditions',[],[],...
  'sourcefact',1,[0,Inf],...
  'baselineflag',true,[false true],...
  'save_synth_flag',false,[false true],...
  'fif_flag',false,[false true],...
  'template_fif',[],[],...
  'write_mgh_flag',false,[false true],...
  'ncov_type',2,[0:2],...
  'calc_scalefacts_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'dSPM_tags',{'subjdir','rootoutdir','prefix','conditions',...
               'calc_dipinfo_flag','lh_dip_info','rh_dip_info',...
               'lh_dip_file','rh_dip_file','nodec_flag',...
               'lh_dec_dips','rh_dec_dips','lh_dec_file','rh_dec_file',...
               'ncov_type','ncov_conditions','raw_ncov_lambda',...
               'calc_scalefacts_flag','baseline_start','baseline_end',...
               'baseline_flag','ssp_projmat','forceflag','SNR',...
               'noisenorm_flag','noisenorm_identity_flag',...
               'depthweight_flag','depthweight_p','bem_flag','radii',...
               'conductivities','EEG_gain_scalefact','nlayers',...
               'badchans','badchanfile','usegrad_flag','usemag_flag',...
               'useEEG_flag','grad_scalefact','mag_scalefact',...
               'EEG_scalefact','write_stc_flag','stc_scalefact',...
               'write_mgh_flag','sparsesmooth','postsmooth',...
               'mbmask_flag','resamp2ico_flag','icolevel',...
               'icosmooth','write_fif_flag','template_fif',...
               'bem_surf_files','cen_sph','trans','transfile',...
               'alignment_fif','lh_sourcecov_file','rh_sourcecov_file',...
               'sourcecov_thresh','sourcecov_thresh_abs_flag',...
               'sourcecov_maxvar','sourcecov_minvar','forward_matfile',...
               'inverse_matfile','refEEG_coords','prewhiten_flag',...
               'orient_constr_flag','orient_tang','smooth_constr_flag',...
               'smooth_constr_nsmooth','signed_sources_flag','datatype',...
               'save_results_flag','save_avg_flag','save_fiterr_flag',...
               'save_fiterr_struct_flag','hemilist'},[],...
});

cwd = pwd;

if isempty(parms.RCSE_rootdir)
  parms.RCSE_rootdir = parms.rootdir;
end;
if isempty(parms.RCSE_forward_rootdir)
  parms.RCSE_forward_rootdir = parms.RCSE_rootdir;
end;
if isempty(parms.RCSE_forward_prefix)
  parms.RCSE_forward_prefix = parms.RCSE_prefix;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load parms for dSPM
fname = sprintf('%s/matfiles/%s_parms.mat',...
  parms.rootdir,parms.dSPM_prefix);
if ~exist(fname,'file'), error('file %s not found',fname); end;
tmp = load(fname);
dSPM_parms = tmp.parms;

% load parms for RCSE source
fname = sprintf('%s/matfiles/%s_parms.mat',...
  parms.RCSE_rootdir,parms.RCSE_prefix);
if ~exist(fname,'file'), error('file %s not found',fname); end;
tmp = load(fname);
RCSE_source_parms = tmp.parms;

% load parms for RCSE forward
fname = sprintf('%s/matfiles/%s_parms.mat',...
  parms.RCSE_forward_rootdir,parms.RCSE_forward_prefix);
if ~exist(fname,'file'), error('file %s not found',fname); end;
tmp = load(fname);
RCSE_forward_parms = tmp.parms;

% merge parms with dSPM parms so we could override existing settings if desired
%dSPM_tags = fieldnames(dSPM_parms);
tags = fieldnames(parms);
for t=1:length(tags)
  dSPM_parms.(tags{t}) = parms.(tags{t});
end;

if isempty(parms.areas)
  parms.areas = [1:RCSE_forward_parms.nareas];
end;

% construct output prefix
if ~isempty(parms.prefix)
  dSPM_parms.prefix = parms.prefix;
elseif strcmp(parms.RCSE_prefix,parms.RCSE_forward_prefix)
  dSPM_parms.prefix = sprintf('%s_%s_%s',...
    parms.outstem,parms.RCSE_prefix,parms.dSPM_prefix);
else
  dSPM_parms.prefix = sprintf('%s_%s_%s_%s',...
    parms.outstem,parms.RCSE_prefix,parms.RCSE_forward_prefix,parms.dSPM_prefix);
end;

% synthesize data from RCSE source estimates + noise
fname = sprintf('%s/matfiles/%s_synth_avg_data.mat',...
  parms.rootoutdir,dSPM_parms.prefix);
if ~exist(fname,'file') || parms.forceflag
  tmp_outdir = [parms.rootoutdir '/matfiles'];
  mmil_mkdir(tmp_outdir);
  tags = {'fiterr_flag','areas',...
          'sourcefact','baselineflag'};
  args = mmil_parms2args(parms,tags);
  avg_data = rc_synth_sensors_from_RCSE(...
    'rootdir',parms.RCSE_rootdir,...
    'prefix',parms.RCSE_prefix,...
    'forward_rootdir',parms.RCSE_forward_rootdir,...
    'forward_prefix',parms.RCSE_forward_prefix,...
    args{:});
  save(fname,'avg_data');
else
  load(fname);
end;

if parms.fif_flag
  tmp_outdir = [parms.rootoutdir '/fifs'];
  mmil_mkdir(tmp_outdir);
  tmp_outstem = sprintf('%s/%s',tmp_outdir,dSPM_parms.prefix);
  tmp_data = avg_data;
  if ~isempty(parms.conditions)
    tmp_data.averages = tmp_data.averages(parms.conditions);
  end;
  ts_avg2fif(tmp_data,parms.template_fif,tmp_outstem,parms.forceflag);
  clear tmp_data;
end;

% run dSPM on synthesized avg_data
dSPM_args = mmil_parms2args(dSPM_parms,parms.dSPM_tags);
cd(dSPM_parms.rootdir)
ts_dSPM(avg_data,dSPM_parms.subjname,dSPM_args{:});
cd(cwd);

