function MEG_MMIL_dSPM_RCSE(ContainerPath,varargin)
%function MEG_MMIL_dSPM_RCSE(ContainerPath,varargin)
%
% Purpose: Create dSPM using RCSE results
%
% Required Input:
%   ContainerPath: full path containing processed MEG data
%
% Optional Parameters:
%  'rootoutdir': root output directory
%     {default = ContainerPath}
%  'prefix': prefix of dSPM_RCSE output files
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
%  'baselineflag': use repeated baseline as noise added to synth_data
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
% For a complete list of options with detailed descriptions, see rc_dSPM_RCSE
%
% Created:  02/21/11 by Don Hagler
% Last Mod: 11/01/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'rootdir',ContainerPath,[],...
  'rootoutdir',ContainerPath,[],...
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
  'forceflag',0,{0,1},...
...
  'ROI_outdir',[],[],... % not used, just accepted from MMIL_Analyze_MEG_Exam
  'ROI_outstem',[],[],...
  'dSRC_tags',{'rootdir','rootoutdir','prefix','outstem',...
               'dSPM_prefix','RCSE_prefix','RCSE_rootdir',...
               'RCSE_forward_prefix','RCSE_forward_rootdir','fiterr_flag',...
               'areas','conditions','sourcefact','baselineflag',...
               'save_synth_flag''fif_flag',...
               'template_fif','write_mgh_flag','ncov_type',...
               'calc_scalefacts_flag','forceflag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: running dSPM-RCSE...\n',mfilename);
args = mmil_parms2args(parms,parms.dSRC_tags);
rc_dSPM_RCSE(args{:});
fprintf('%s: finished\n',mfilename);

