function MEG_MMIL_Synth_RCSE(ContainerPath,varargin)
%function MEG_MMIL_Synth_RCSE(ContainerPath,varargin)
%
% Purpose: Synthesize sensor data using RCSE results
%
% Required Input:
%   ContainerPath: full path containing processed MEG data
%
% Optional Parameters:
%  'rootoutdir': root output directory
%     {default = ContainerPath}
%  'outstem': prefix of output files
%     {default = synth_RCSE}
%  'prefix': prefix of RCSE files to be used
%    {default = 'RCSE'}
%  'forward_prefix': prefix of RCSE output files with different forward
%    If empty, use prefix
%    {default = []}
%  'rootdir': root input directory with matfiles subdir
%    {default = ContainerPath}
%  'forward_rootdir': directory containing matfiles dir for different forward
%    If empty, use rootdir
%    {default = []}
%  'fiterr_flag': which type of sensor data
%    0: RCSE fitted data plus noise
%    1: residual error of RCSE fit
%    2: RCSE fit projected into sensor data plus noise
%    3: data minus RCSE fit projected into sensor data
%     {default = 0}
%  'areas': when fiterr_flag=0, use only these area numbers
%    if empty, use all areas
%    {default = []}
%  'sourcefact': amplitude of modeled sources
%    {default = 1}
%  'baselineflag': [0|1] use repeated baseline as noise added to synth_data
%    {default = 0}
%  'sources': matrix of source time courses
%     if empty, will use RCSE estimated source time courses
%     size must be [ntpoints,narea,ncontrasts]
%     only used if fiterr_flag = 0
%    {default = []}
%  'fif_flag': [0|1] save synthesized sensor data to fif file
%     {default = 0}
%  'template_fif': template fif file (e.g. raw or online avg)
%     Required for saving synthesized data to fif file
%     {default = []}
%  'conditions': vector of condition numbers (index to avg_data.averages)
%     to be saved as fif files
%     {default = []}
%  'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  10/25/13 by Don Hagler
% Last Mod: 10/25/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'rootoutdir',ContainerPath,[],...
  'outstem','synth_RCSE',[],...
...
  'prefix','RCSE',[],...
  'forward_prefix',[],[],...
  'rootdir',ContainerPath,[],...
  'forward_rootdir',[],[],...
  'fiterr_flag',0,[0:3],...
  'areas',[],[],...
  'sourcefact',1,[-Inf,Inf],...
  'baselineflag',false,[false true],...
  'sources',[],[],...
...
  'fif_flag',false,[false true],...
  'template_fif',[],[],...
  'conditions',[],[],...
  'forceflag',false,[false true],...
...
  'synth_tags',{'prefix','forward_prefix','rootdir','forward_rootdir',...
                'fiterr_flag','areas','sourcefact','baselineflag',...
                'sources'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% synthesize sensor data
fname = sprintf('%s/matfiles/%s_synth_avg_data.mat',...
  parms.rootoutdir,parms.outstem);
if ~exist(fname,'file') || parms.forceflag
  fprintf('%s: synthesizing sensor data from RCSE results...\n',mfilename);
  outdir = [parms.rootoutdir '/matfiles'];
  mmil_mkdir(outdir);
  args = mmil_parms2args(parms,parms.synth_tags);
  avg_data = rc_synth_sensors_from_RCSE(args{:});
  save(fname,'avg_data');
elseif parms.fif_flag
  load(fname);
end;

% save fif files
if parms.fif_flag
  outdir = [parms.rootoutdir '/fifs'];
  mmil_mkdir(outdir);
  outstem = sprintf('%s/%s',outdir,parms.outstem);
  tmp_data = avg_data;
  if ~isempty(parms.conditions)
    tmp_data.averages = tmp_data.averages(parms.conditions);
  end;
  ts_avg2fif(tmp_data,parms.template_fif,outstem,parms.forceflag);
  clear tmp_data;
end;


fprintf('%s: finished\n',mfilename);

