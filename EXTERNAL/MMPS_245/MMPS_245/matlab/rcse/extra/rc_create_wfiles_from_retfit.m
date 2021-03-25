function rc_create_wfiles_from_retfit(subjname,varargin)
%function rc_create_wfiles_from_retfit(subjname,[options])
%
% Purpose: get patch of vertices from retfit, save as w files
%
% Usage:
%  rc_create_wfiles_from_retfit(subjname,'key1', value1,...);
%
% Required Input:
%  subjname - freesurfer recon subject name
%
% Optional parameters:
%  'retfit_stem' - file stem for refit results mat file
%    {default = 'retfit/retfit_results'}
%  'fname_conds' - full path of csv file containing condition information
%    {default: 'cond_info.csv'}
%  'ecc_width' - eccentricity width of stimuli (deg. vis. ang.)
%    will be ignored if fname_conds has ecc_width column
%    {default = 1}
%  'theta_width' - polar angle width of stimuli (degrees)
%    will be ignored if fname_conds has theta_width column or if stim_type=1
%    {default = 10}
%  'r_offset' - eccentricity offset
%    {default = 0}
%  'th_offset' - polar angle offset
%    {default = 0}
%  'w_thresh' - threshold applied to weights
%    {default = 0.05}
%  'single_vertex_flag' - [0|1] select one vertex for each stimulus location
%    0: use all vertices
%    1: vertex with maximum weight
%    2: vertex closest to center of mass
%    {default = 0}
%  'forward' - struct containing RCSE forward matrix and other things
%    including lh_dip_info and rh_dip_info
%    required if single_vertex_flag = 2
%    {default = []}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'subjdir' - root directory containing freesurfer recons
%    {default: $SUBJECTS_DIR}
%  'surfname' - surface file to load (for smoothing)
%    {default: 'white'}
%  'hemilist' - cell array containing 'lh' and/or 'rh'
%    {default: {'lh' 'rh'}}
%  'outdir' - output directory for w files
%    {default: './nbrhoods'}
%  'matdir' - output directory for mat file
%    {default: './matfiles'}
%  'rf_sizes': vector of receptive field sizes for each visual area
%    Will be used as initial estimate if rf_niters>0
%    {default = [1,1,1]}
%  'rf_slopes': vector of slopes of linear trend of receptive field sizes
%    w.r.t. ecc for each visual area
%    Intercept is assumed to be half of r_max
%    {default = [0.1,0.1,0.1]}
%  'restrict_hemi_flag': [0|1] whether to restrict dipole clusters
%    so no ipsilateral
%    {default = 0}
%  'restrict_uplow_flag': [0|1] whether to restrict dipole clusters
%    so no upper-lower cross-over
%    {default = 0}
%  'vf2ctx': struct containing lh and rh vf2ctx sparce matrices
%    {default = []}
%  'stim_type' - type of stimulus to model
%    0: point - generation of retmap is very fast, but less accurate
%    1: circle - slower, but more accurate
%    2: wedge - slower, but more accurate
%    {default = 2}
%  'forceflag' - [0|1] whether to overwrite existing output
%    {default: 0}
%
% Created:  07/13/09  by Don Hagler
% Last Mod: 12/08/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'retfit_stem','retfit/retfit_results',[],...
  'fname_conds','cond_info.csv',[],...
  'ecc_width',1,[0,100],...
  'theta_width',10,[0,360],...
  'r_offset',0,[-100,100],...
  'th_offset',0,[-180,180],...
  'w_thresh',0.05,[0,100],...
  'single_vertex_flag',0,[0 1 2],...
  'forward',[],[],...
  'r_max',12.5,[0,Inf],...
  'subjdir',[],[],...
  'surfname','white',[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'outdir','./nbrhoods',[],...
  'matdir','./matfiles',[],...
  'rf_sizes',[1,1,1],[],...
  'rf_slopes',[0.1,0.1,0.1],[0,10],...
  'restrict_hemi_flag',false,[false true],...
  'restrict_uplow_flag',false,[false true],...
  'restrict_flag',false,[false true],...
  'vf2ctx',[],[],...
  'stim_type',2,[0,1,2],... 
  'forceflag',false,[false true],...
...
  'contrasts',1,[],...
  'nbrhood_flag',false,[false true],...
  'smoothed_areas_flag',true,[false true],...
  'polar_flag',true,[false true],...
});

if parms.restrict_flag
  parms.restrict_hemi_flag = 1;
  parms.restrict_uplow_flag = 1;
end;

if isempty(parms.subjdir)
  parms.subjdir = getenv('SUBJECTS_DIR');
end;
if isempty(parms.subjdir)
  error('no subjdir specified',parms.subjdir);
end;
if ~exist(parms.subjdir,'dir')
  error('subjdir %s not found',parms.subjdir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load retfit files
retfit_results = [];
hemilist = {'lh','rh'};
for h=1:length(hemilist)
  hemi = hemilist{h};
  fname = sprintf('%s-%s.mat',parms.retfit_stem,hemi);
  if ~exist(fname,'file'), error('file %s not found',fname); end;
  if strcmp(hemi,'lh')
    clear fit_data area_masks fit_parms;
    load(fname);
    retfit_results.lh_fit_data = fit_data;
    retfit_results.lh_area_masks = area_masks;
  else
    clear fit_data area_masks fit_parms;
    load(fname);
    retfit_results.rh_fit_data = fit_data;
    retfit_results.rh_area_masks = area_masks;
  end;
end;

fprintf('%s: creating wfiles for %s...\n',...
  mfilename,parms.fname_conds);

% define retmap
retmap = rc_define_retmap_from_retfit(...
  retfit_results,...
  'fname_conds',parms.fname_conds,...
  'contrasts',parms.contrasts,...
  'r_offset',parms.r_offset,...
  'th_offset',parms.th_offset,...
  'ecc_width',parms.ecc_width,...
  'theta_width',parms.theta_width,...
  'r_max',parms.r_max,...
  'rf_sizes',parms.rf_sizes,...
  'rf_slopes',parms.rf_slopes,...
  'smoothed_areas_flag',parms.smoothed_areas_flag,...
  'restrict_hemi_flag',parms.restrict_hemi_flag,...
  'restrict_uplow_flag',parms.restrict_uplow_flag,...
  'stim_type',parms.stim_type,...
  'vf2ctx',parms.vf2ctx,...
  'polar_flag',parms.polar_flag,...
  'single_vertex_flag',parms.single_vertex_flag,...
  'forward',parms.forward,...
  'w_thresh',parms.w_thresh);

rc_create_wfiles_from_retmap_verts(subjname,retmap,...
  'nbrhood_flag',parms.nbrhood_flag,'outdir',parms.outdir,...
  'forceflag',parms.forceflag);


