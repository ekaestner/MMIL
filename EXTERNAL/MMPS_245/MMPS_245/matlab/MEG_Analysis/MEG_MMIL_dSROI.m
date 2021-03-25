function errcode = MEG_MMIL_dSPM_ROI(ContainerPath,varargin)
%function errcode = MEG_MMIL_dSPM_ROI(ContainerPath,[options])
%
% Purpose: run Dynamic Statistical Parametric Mapping
%
% Required Input:
%   ContainerPath: full path containing processed MEG data
%
% Optional Parameters:
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
%    {default: 'dSPM_ROI_analysis'}
%  'outstem': output file stem
%    {default: 'dSPM_ROI_results'}
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
%  'overlay_roi_flag' [0|1] plot with ROIs overlayed on same plot
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
%    {default: 'NorthEast'}
%  'visible_flag': [0|1] make plots visible (otherwise plot in background)
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% For a complete list of options with detailed descriptions, see ts_dSPM
%
% Created:  02/25/11 by Don Hagler
% Last Mod: 04/26/13 by Don Hagler
%

%% todo: rootoutdir in ContainerPath/dSPM ?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'rootdir',ContainerPath,[],...
  'indir','stcfiles',[],...
  'prefix','dSPM',[],...
  'conditions',[],[],...
  'outdir','dSPM_ROI_analysis',[],...
  'outstem','dSPM_ROI_results',[],...
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
  'legend_loc','NorthEast',[],...
  'avg_hemis_flag',false,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: running dSPM ROI analysis...\n',mfilename);
args = mmil_parms2args(parms);
ts_roi_analysis(args{:});
fprintf('%s: finished\n',mfilename);

