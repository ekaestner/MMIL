function rc_retfit(subj,varargin)
%function rc_retfit(subj,[options])
%
% Required Parameters:
%   subj: freesurfer subject (directory name)
%
% Optional Parameters:
%  'subjdir': freesurfer root subject dir
%    {default = SUBJECTS_DIR} (environment variable)
%  'hemilist': cell array of cortical hemispheres to fit
%    {default = {'lh','rh'}}
%  'prereg_niter': pre-registration iterations
%    (random search with global map parameters)
%    {default = 100}
%  'prereg_nruns': number of runs of pre-registration
%    if greater than 1,  will randomly choose seed point from parameter ranges
%    {default = 200}
%  'nruns': number of runs of nonlinear fitting
%    (updated fit parameters saved after each run)
%    {default = 1}
%  'run_reset_flag': [0|1] whether to set du and dv to 0 at start of each 
%    fine fitting run (and update u0 and v0)
%    {default = 0}
%  'niter': maximum number of iterations within nonlinear fit
%    {default = 1000}
%  'cost_include_percentile': map nodes with error cost greater than this
%    percentile will not contribute to error cost for fine scale fit
%    {default = 100}
%  'area_smooth_steps': number of sparse smoothing steps to apply to
%    area masks
%    {default = 2}
%  'save_results_flag': [0|1] whether to save results to mat file
%    This also involves resampling grid to vertices, which is slow
%    Set to 0 if you want to plot results and interactively adjust parameters
%    NOTE: fit parameters are saved regardless, allowing results to be
%     saved on subsequent runs
%    {default = 1}
%  'save_mgh_flag': [0|1] whether to save fitted data to mgh files
%    Ignored if 'save_results_flag' = 0
%    {default = 1}
%  'save_label_flag': [0|1] whether to save area masks as label files
%    Ignored if 'save_results_flag' = 0
%    {default = 1}
%  'resample_ico_flag': [0|1] whether to save area masks as ico4 label files
%    Ignored if 'save_results_flag' = 0 or 'save_label_flag' = 0
%    {default = 1}
%  'plot_flag': [0|1|2|3] whether to plot images
%    0: no plots
%    1: plots of fit results
%    2: plots of initial registration (pauses before fitting) and fit results
%    3: plots of initial and fitted registration and more
%    {default = 1}
%  'outdir': output directory
%    {default = pwd}
%  'outstem':  file stem for output files
%    {default = 'retfit'}
%  'forceflag': [0|1] whether to overwrite files (e.g. data sampled to grid)
%    {default = 0}
%
% Optional Data Parameters:
%  'polstem': file stem of polar angle retinotopy data (mgh format)
%    May include full path
%    {default = 'polret'}
%  'eccstem': file stem of eccentricity retinotopy data (mgh format)
%    {default = 'eccret'}
%  'data_r_min': minimum radius (degrees visual angle)
%     used for eccentricity mapping - determines phase for a given eccentricity
%    {default = 0.25}
%  'data_r_max': maximum radius (degrees visual angle)
%     used for eccentricity mapping - determines phase for a given eccentricity
%    {default = 12.5}
%  'data_logtrans_flag': [0|1] whether eccentricity mapping data was collected
%    with log transform
%    {default = 0}
%   Any changes to these parameters require running with forceflag=1
%
% Optional ROI Parameters:
%   For initial resampling of roi to unit grid coordinates (u and v)
%   Each of these should be vectors with lengths corresponding to number
%     of elements in hemilist (one or each cortical hemisphere)
%   Entire V1-V2-V3 complex must fit unit grid (0 - 1)
%   Ideally, horizontal meridian of V1 should be lined up with u=0.5
%     and mid eccentricity should line up with v=0.5
%  'roi_name': file stem for ROI label file (e.g. lh.roi_name.label)
%    {default = 'v123'}
%  'roi_excl_name': file stem for exclusion ROI label file
%    i.e. region to ignore in calculating error
%    note: this label file is NOT required to exist
%    {default = 'exclude'}
%  'roi_indir': input directory containing ROI label files
%    {default = pwd}
%  'roi_dilate_niters': number of smoothing iterations to dilate ROI
%    this will create a border around map fit area
%    {default = 0}
%  'roi_rotation': degrees of clockwise rotation
%    {default = [0 0]}
%  'roi_shift_u': horizontal shift on unit grid (after rotation and scale)
%    {default = [0 0]}
%  'roi_shift_v': vertical shift on unit grid (after rotation and scale)
%    {default = [0 0]}
%  'roi_scale_u': horizontal scaling factor
%    {default = [0.7 0.7]}
%  'roi_scale_v': vertical scaling factor
%    {default = [0.7 0.7]}
%   Any changes to these parameters require running with forceflag=1
%   NOTE: If number of hemis in hemilist=1, but number of parameters=2
%     second parameter will be ignored
%         If number of hemis in hemilist=2, but number of parameters=1
%     that parameter will be used for both hemis
%
% Optional Map Template Parameters:
%   Each of these except 'map_v123_flag', 'map_poly_flag', 'map_poly_order',
%     map_poly_coef', 'map_model_type', 'map_area_name', and 'map_rev_polar_flag'
%     should be vectors with lengths corresponding to number of elements
%     in hemilist (one for each cortical hemisphere)
%  'map_v123_flag': [0|1] whether to model V1-V2-V3 complex or
%    a single area mapping entire hemifield
%    {default = 1}
%  'map_area_name': area label if map_v123_flag=0
%    {default = 'v'}
%  'map_poly_flag': [0|1] use polynomial function to deform template
%    {default = 1}
%  'map_poly_order': order of polynomial function (n+1 additional parameters)
%    used to deform template
%    {default = 4}
%  'map_model_type': [0|1|2] model used for initial estimates of u and v
%    0: rectangle
%    1: wedge
%    2: radial wedge
%    {default = 2}
%  'map_rev_polar_flag': [0|1] whether to reverse direction of polar angle
%    {default = 0}
%  'map_v1p_width': width (along theta axis) of v1+ map
%    {default = [1 1]}
%  'map_v1m_width': width (along theta axis) of v1- map
%    {default = [1 1]}
%  'map_v2p_width': width (along theta axis) of v2+ map
%    {default = [1 1]}
%  'map_v2m_width': width (along theta axis) of v2- map
%    {default = [1 1]}
%  'map_v3p_width': width (along theta axis) of v3+ map
%    {default = [1 1]}
%  'map_v3m_width': width (along theta axis) of v3- map
%    {default = [1 1]}
%  'map_v1p_length': length (along ecc axis) of v1+ map
%    {default = [1 1]}
%  'map_v1m_length': length (along ecc axis) of v1- map
%    {default = [1 1]}
%  'map_v2p_length': length (along ecc axis) of v2+ map
%    {default = [1 1]}
%  'map_v2m_length': length (along ecc axis) of v2- map
%    {default = [1 1]}
%  'map_v3p_length': length (along ecc axis) of v3+ map
%    {default = [1 1]}
%  'map_v3m_length': length (along ecc axis) of v3- map
%    {default = [1 1]}
%  'map_poly_coef': vector of polynomial function coefficients
%     must have number of elements equal to map_poly_order + 1
%    {default = zeros(1,map_poly_order+1)}
%  'map_wedge_fact': how much narrower map's foveal edge is compared to peripheral
%    {default = [1 1]}
%  'map_radial_wedge_fact': how horizontal extent of map translates to polar
%    angle (smaller value means less bending)
%    {default = [0.2 0.2]}
%  'map_radial_offset': value added to vertical map coordinates before
%    radial wedge transform (smaller value means sharper point)
%    {default = [1 1]}
%  'map_scale_u': horizontal scaling factor for initial u coordinates
%     for template map
%    {default = [0.7 0.7]}
%  'map_scale_v': vertical scaling factor for initial v coordinates
%     for template map
%    {default = [0.4 0.4]}
%  'map_rotation': degrees of clockwise rotation for initial template map
%    {default = [0 0]}
%  'map_shift_u': horizontal shift on unit grid for initial template map
%    {default = [0 0]}
%  'map_shift_v': vertical shift on unit grid for initial template map
%    {default = [0 0]}
%  'map_r_min': minimum radius (degrees visual angle) in template map
%    {default = [2 2]}
%  'map_r_max': maximum radius (degrees visual angle) in template map
%    {default = [12 12]}
%   Any changes to these parameters require running with forceflag=1
%     or deleting *parm* mat files
%
% Optional Map Template Preregistration Parameters:
%  'map_v1_width_range': vector of lower and upper bounds for map parameter
%    {default = [1 1]}
%  'map_v2_width_range'
%    {default = [0.6 1]}
%  'map_v3_width_range'
%    {default = [0.5 1]}
%  'map_v1_length_range'
%    {default = [0.8 1.2]}
%  'map_v2_length_range'
%    {default = [0.8 1.2]}
%  'map_v3_length_range'
%    {default = [0.8 1.2]}
%  'map_poly_coef_range'
%    {default = [-5,5]}
%  'map_wedge_fact_range'
%    {default = [0.7 1.3]}
%  'map_radial_wedge_fact_range'
%    {default = [0.05 0.4]}
%  'map_radial_offset_range'
%    {default = [1 4]}
%  'map_scale_u_range'
%    {default = [0.6 0.8]}
%  'map_scale_v_range'
%    {default = [0.3 0.6]}
%  'map_rotation_range'
%    {default = [-10 10]}
%  'map_shift_u_range'
%    {default = [-0.1,0.1]}
%  'map_shift_v_range'
%    {default = [-0.1,0.1]}
%  'map_r_min_range'
%    {default = [1 1]}
%  'map_r_max_range'
%    {default = [12 12]}
%
% Output:
%   Creates results mat file containing fit_data, area_masks, and fit_parms
%
% Created:  01/15/09 by Don Hagler
% Last Mod: 10/26/13 by Don Hagler
%

%% todo: expect cell array with two elements containing vectors for
%%       map_poly_coef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'subjdir',[],[],...
  'nruns',1,[0,Inf],...
  'run_reset_flag',false,[false true],...
  'niter',1000,[1,Inf],...
  'prereg_nruns',200,[0,Inf],...
  'prereg_niter',100,[0,Inf],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'save_results_flag',true,[false true],...
  'save_mgh_flag',true,[false true],...
  'save_label_flag',true,[false true],...
  'resample_ico_flag',true,[false true],...
  'plot_flag',1,[0,1,2,3],...
  'outdir',pwd,[],...
  'outstem','retfit',[],...
  'forceflag',false,[false true],...
  'polstem','polret',[],...
  'eccstem','eccret',[],...
  'data_r_min',0.25,[0,Inf],...
  'data_r_max',12.5,[0,Inf],...
  'data_logtrans_flag',false,[false true],...
  'roi_name','v123',[],...
  'roi_excl_name','exclude',[],...
  'roi_indir',pwd,[],...
  'roi_dilate_niters',0,[0,1000],...
  'roi_rotation',[0 0],[-180,180],...
  'roi_shift_u',[0 0],[-1,1],...
  'roi_shift_v',[0 0],[-1,1],...
  'roi_scale_u',[0.7 0.7],[0.1,2],...
  'roi_scale_v',[0.7 0.7],[0.1,2],...
  'map_r_min',[2 2],[0,Inf],...
  'map_r_max',[12 12],[0,Inf],...
  'map_v1p_width',[1 1],[0.01,2],...
  'map_v1m_width',[1 1],[0.01,2],...
  'map_v2p_width',[1 1],[0.01,2],...
  'map_v2m_width',[1 1],[0.01,2],...
  'map_v3p_width',[1 1],[0.01,2],...
  'map_v3m_width',[1 1],[0.01,2],...
  'map_v1p_length',[1 1],[0.1,10],...
  'map_v1m_length',[1 1],[0.1,10],...
  'map_v2p_length',[1 1],[0.1,10],...
  'map_v2m_length',[1 1],[0.1,10],...
  'map_v3p_length',[1 1],[0.1,10],...
  'map_v3m_length',[1 1],[0.1,10],...
  'map_scale_u',[0.7 0.7],[0.01,1.0],...
  'map_scale_v',[0.4 0.4],[0.01,1.0],...
  'map_rotation',[0 0],[-180,180],...
  'map_shift_u',[0 0],[-1,1],...
  'map_shift_v',[0 0],[-1,1],...
  'map_poly_coef',[],[],...
  'map_wedge_fact',[1 1],[0.01,2],...
  'map_radial_wedge_fact',[0.2 0.2],[0,1],...
  'map_radial_offset',[1 1],[0,10],...
  'map_v1_width_range',[1 1],[0.01,2],...
  'map_v2_width_range',[0.6 1],[0.01,2],...
  'map_v3_width_range',[0.5 1],[0.01,2],...
  'map_v1_length_range',[0.8 1.2],[0.1 10],...
  'map_v2_length_range',[0.8 1.2],[0.1 10],...
  'map_v3_length_range',[0.8 1.2],[0.1 10],...
  'map_poly_coef_range',[-5,5],[-100,100],...
  'map_wedge_fact_range',[0.7 1.3],[0.01,2],...
  'map_radial_wedge_fact_range',[0.05 0.4],[0,1],...
  'map_radial_offset_range',[1 4],[0,10],...
  'map_scale_u_range',[0.6 0.8],[0.01,1.0],...
  'map_scale_v_range',[0.3 0.6],[0.01,1.0],...
  'map_rotation_range',[-10 10],[-180,180],...
  'map_shift_u_range',[-0.1,0.1],[-1,1],...
  'map_shift_v_range',[-0.1,0.1],[-1,1],...
  'map_r_min_range',[1 1],[0,Inf],...
  'map_r_max_range',[12 12],[0,Inf],...
  'map_v123_flag',true,[false true],...
  'map_poly_flag',true,[false true],...
  'map_poly_order',4,[1,10],...
  'map_model_type',2,[0,1,2],...
  'map_area_name','v',[],...
  'map_rev_polar_flag',false,[false true],...
  'cost_include_percentile',100,[50,100],...
  'area_smooth_steps',2,[0,100],...
... % undocumented
  'prereg_search_type','rand',{'rand','fmincon'},...
  'map_logtrans_flag',true,[false true],...
  'data_gridsize_u',500,[100,1e10],...
  'data_gridsize_v',500,[100,1e10],...
  'map_gridsize_r',40,[2,1000],...
  'map_gridsize_th',40,[2,1000],...
  'step_size',1,[0,100],...
  'max_step_size',1,[0,1000],...
  'prereg_step_size',0.1,[0.001,1],...
  'ms_trans',12,[0,100],...
  'ms_scale',0.6,[0,10],...
  'err_fact',1,[0,Inf],...
  'data_smooth_sigma',0.15,[0,1],...
  'grad_smooth_sigma',0,[0,1],...
  'ecc_fact',1,[0,1],...
  'smooth_fact',1,[0,Inf],...
  'fold_fact',15,[0 Inf],...
  'vacancy_fact',0,[0,Inf],...
  'max_outbound_penalty',20,[0 10000],...
  'step_type','msgrad',{'rand','grad','msgrad'},...
  'interp_method','linear',{'linear','cubic','nearest'},...
  'area_code_smf',1e-5,[0,1],...
  'check_grad_flag',false,[false true],...
  'ico',4,[1:7],...
  'roi_excl',[],[],...
});

excl_tags = {'subjdir' 'ico' ...
  'polstem' 'eccstem' 'hemilist' ...
  'save_results_flag' 'save_mgh_flag' 'save_label_flag' ...
  'resample_ico_flag' 'roi_name' 'roi_excl_name' 'roi_indir' ...
  'fit_outdir' 'label_outdir'};
hemi_tags = {'roi_rotation' 'roi_scale_u' 'roi_scale_v' ...
  'roi_shift_u' 'roi_shift_v' 'map_wedge_fact' ...
  'map_radial_wedge_fact' 'map_radial_offset' 'map_rotation'...
  'map_scale_u' 'map_scale_v' 'map_shift_u' 'map_shift_v' ...
  'map_v1p_width' 'map_v1m_width'...
  'map_v2p_width' 'map_v2m_width'...
  'map_v3p_width' 'map_v3m_width'...
  'map_v1p_length' 'map_v1m_length' ...
  'map_v2p_length' 'map_v2m_length'...
  'map_v3p_length' 'map_v3m_length'...
  'map_r_min' 'map_r_max'};
tags = setdiff(fieldnames(parms),excl_tags);

if isempty(parms.subjdir)
  parms.subjdir = getenv('SUBJECTS_DIR');
end;
if isempty(parms.subjdir)
  error('SUBJECTS_DIR not specified');
end;
if ~exist(parms.subjdir,'dir')
  error('subjdir %s not found',parms.subjdir);
end;
subjpath = [parms.subjdir '/' subj];
if ~exist(subjpath,'dir')
  error('freesurfer recon dir %s not found',subjpath);
end;

parms.fit_outdir = [parms.outdir '/fit'];
parms.label_outdir = [parms.outdir '/label'];
parms.outdir = [parms.outdir '/matfiles'];

switch parms.step_type
  case 'grad'
    parms.max_nsteps_inc_cost = 0;
  case 'msgrad'
    parms.max_nsteps_inc_cost = parms.niter;
end;

if ~iscell(parms.hemilist)
  parms.hemilist = {parms.hemilist};
end;
nhemis = length(parms.hemilist);
if nhemis > 2
  error('number of hemispheres in hemilist (%d) greater than 2',nhemis);
end;
% for roi_ and map_ parameters, allow to be vector of two or one value
%   if one value but nhemis==2, use for both hemis (but negative for lh)
if length(parms.roi_rotation)<nhemis
  parms.roi_rotation = [-parms.roi_rotation,parms.roi_rotation];
end;
if length(parms.roi_shift_u)<nhemis
  parms.roi_shift_u = [parms.roi_shift_u,parms.roi_shift_u];
end;
if length(parms.roi_shift_v)<nhemis
  parms.roi_shift_v = [parms.roi_shift_v,parms.roi_shift_v];
end;
if length(parms.roi_scale_u)<nhemis
  parms.roi_scale_u = [parms.roi_scale_u,parms.roi_scale_u];
end;
if length(parms.roi_scale_v)<nhemis
  parms.roi_scale_v = [parms.roi_scale_v,parms.roi_scale_v];
end;
if length(parms.map_r_min)<nhemis
  parms.map_r_min = [parms.map_r_min,parms.map_r_min];
end;
if length(parms.map_r_max)<nhemis
  parms.map_r_max = [parms.map_r_max,parms.map_r_max];
end;
if length(parms.map_v1p_width)<nhemis
  parms.map_v1p_width = [parms.map_v1p_width,parms.map_v1p_width];
end;
if length(parms.map_v1m_width)<nhemis
  parms.map_v1m_width = [parms.map_v1m_width,parms.map_v1m_width];
end;
if length(parms.map_v2p_width)<nhemis
  parms.map_v2p_width = [parms.map_v2p_width,parms.map_v2p_width];
end;
if length(parms.map_v2m_width)<nhemis
  parms.map_v2m_width = [parms.map_v2m_width,parms.map_v2m_width];
end;
if length(parms.map_v3p_width)<nhemis
  parms.map_v3p_width = [parms.map_v3p_width,parms.map_v3p_width];
end;
if length(parms.map_v3m_width)<nhemis
  parms.map_v3m_width = [parms.map_v3m_width,parms.map_v3m_width];
end;
if length(parms.map_v1p_length)<nhemis
  parms.map_v1p_length = [parms.map_v1p_length,parms.map_v1p_length];
end;
if length(parms.map_v1m_length)<nhemis
  parms.map_v1m_length = [parms.map_v1m_length,parms.map_v1m_length];
end;
if length(parms.map_v2p_length)<nhemis
  parms.map_v2p_length = [parms.map_v2p_length,parms.map_v2p_length];
end;
if length(parms.map_v2m_length)<nhemis
  parms.map_v2m_length = [parms.map_v2m_length,parms.map_v2m_length];
end;
if length(parms.map_v3p_length)<nhemis
  parms.map_v3p_length = [parms.map_v3p_length,parms.map_v3p_length];
end;
if length(parms.map_v3m_length)<nhemis
  parms.map_v3m_length = [parms.map_v3m_length,parms.map_v3m_length];
end;
if length(parms.map_scale_u)<nhemis
  parms.map_scale_u = [parms.map_scale_u,parms.map_scale_u];
end;
if length(parms.map_scale_v)<nhemis
  parms.map_scale_v = [parms.map_scale_v,parms.map_scale_v];
end;
if length(parms.map_rotation)<nhemis
  parms.map_rotation = [parms.map_rotation,parms.map_rotation];
end;
if length(parms.map_shift_u)<nhemis
  parms.map_shift_u = [parms.map_shift_u,parms.map_shift_u];
end;
if length(parms.map_shift_v)<nhemis
  parms.map_shift_v = [parms.map_shift_v,parms.map_shift_v];
end;
if length(parms.map_wedge_fact)<nhemis
  parms.map_wedge_fact = [parms.map_wedge_fact,parms.map_wedge_fact];
end;
if length(parms.map_radial_wedge_fact)<nhemis
  parms.map_radial_wedge_fact = [parms.map_radial_wedge_fact,parms.map_radial_wedge_fact];
end;
if length(parms.map_radial_offset)<nhemis
  parms.map_radial_offset = [parms.map_radial_offset,parms.map_radial_offset];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  % load surfaces
  matfile = sprintf('%s/%s_surfs-%s.mat',parms.outdir,parms.outstem,hemi);
  if ~exist(matfile,'file')
    fprintf('%s: loading %s cortical surfaces...\n',mfilename,hemi);
    % load surface
    fname_sph = sprintf('%s/%s/surf/%s.sphere',...
      parms.subjdir,subj,hemi);
    fname_wm = sprintf('%s/%s/surf/%s.white',...
      parms.subjdir,subj,hemi);
    if ~exist(fname_sph,'file')
      error('file %s not found',fname_sph);
    end;
    if ~exist(fname_wm,'file')
      error('file %s not found',fname_wm);
    end;
    surf = fs_read_surf(fname_sph);
    surf = fs_find_neighbors(surf);
    surf_wm = fs_read_surf(fname_wm);
    if surf.nverts ~= surf_wm.nverts
      error('mismatch between %s and %s',fname_sph,fname_wm);
    end;
    surf.vertices_wm = surf_wm.vertices;
    save(matfile,'surf');
  else
    load(matfile);
  end;
  % load data
  matfile = sprintf('%s/%s_data-%s.mat',parms.outdir,parms.outstem,hemi);
  if ~exist(matfile,'file')
    fprintf('%s: loading %s data...\n',mfilename,hemi);
    % load map data
    fname = sprintf('%s_r-%s.mgh',parms.polstem,hemi);
    pol_r = mmil_colvec(fs_load_mgh(fname));
    fname = sprintf('%s_i-%s.mgh',parms.polstem,hemi);
    pol_i = mmil_colvec(fs_load_mgh(fname));

    fname = sprintf('%s_r-%s.mgh',parms.eccstem,hemi);
    ecc_r = mmil_colvec(fs_load_mgh(fname));
    fname = sprintf('%s_i-%s.mgh',parms.eccstem,hemi);
    ecc_i = mmil_colvec(fs_load_mgh(fname));
    clear data;
    data.pol_r = pol_r;
    data.pol_i = pol_i;
    data.ecc_r = ecc_r;
    data.ecc_i = ecc_i;
    save(matfile,'data');
  else
    load(matfile);
  end;
  % load ROI
  fname = sprintf('%s/%s.%s.label',parms.roi_indir,hemi,parms.roi_name);
  fprintf('%s: loading roi %s...\n',mfilename,fname);
  mask = fs_read_label(fname);
  % load excl ROI
  fname = sprintf('%s/%s.%s.label',parms.roi_indir,hemi,parms.roi_excl_name);
  if exist(fname,'file')
    fprintf('%s: loading roi_excl %s...\n',mfilename,fname);
    mask_excl = fs_read_label(fname);
  else
    mask_excl = [];
  end;
  % run rc_retfit_func
  tmp_parms = parms;
  for i=1:length(hemi_tags)
    tmp_parms.(hemi_tags{i}) = parms.(hemi_tags{i})(h);
  end;
  tmp_parms.roi_excl = mask_excl;
  args = mmil_parms2args(tmp_parms,tags);
  if parms.save_results_flag
    matfile = sprintf('%s/%s_results-%s.mat',parms.outdir,parms.outstem,hemi);
    if ~exist(matfile,'file') || parms.forceflag
      [fit_data,area_masks,fit_parms,data] = ...
         rc_retfit_func(surf,data,mask,hemi,args{:});
      save(matfile,'fit_parms','fit_data','area_masks','data');
    else
      load(matfile);
    end;
    if parms.save_mgh_flag
      mmil_mkdir(parms.fit_outdir);
      fname = sprintf('%s/%s_pol_r-%s.mgh',parms.fit_outdir,parms.outstem,hemi);
      if ~exist(fname,'file') || parms.forceflag
        fs_save_mgh(full(fit_data.pol_r),fname);
      end;
      fname = sprintf('%s/%s_pol_i-%s.mgh',parms.fit_outdir,parms.outstem,hemi);
      if ~exist(fname,'file') || parms.forceflag
        fs_save_mgh(full(fit_data.pol_i),fname);
      end;
      fname = sprintf('%s/%s_ecc_r-%s.mgh',parms.fit_outdir,parms.outstem,hemi);
      if ~exist(fname,'file') || parms.forceflag
        fs_save_mgh(full(fit_data.ecc_r),fname);
      end;
      fname = sprintf('%s/%s_ecc_i-%s.mgh',parms.fit_outdir,parms.outstem,hemi);
      if ~exist(fname,'file') || parms.forceflag
        fs_save_mgh(full(fit_data.ecc_i),fname);
      end;
    end;
    if parms.save_label_flag
      mmil_mkdir(parms.label_outdir);
      for a=1:length(area_masks)
        fname = sprintf('%s/%s.%s_%s.label',...
          parms.label_outdir,hemi,parms.outstem,area_masks(a).name);
        if ~exist(fname,'file') || parms.forceflag
          if parms.area_smooth_steps
            verts = area_masks(a).vertices_smoothed;
          else
            verts = area_masks(a).vertices;
          end;
          fs_write_label(verts,fname,subj);
        end;
      end;
      if parms.resample_ico_flag
        setenv('SUBJECTS_DIR',parms.subjdir);
        for a=1:length(area_masks)
          fname_in = sprintf('%s/%s.%s_%s.label',...
            parms.label_outdir,hemi,parms.outstem,area_masks(a).name);
          fname_out = sprintf('%s/%s.%s_%s.ico%d.label',...
            parms.label_outdir,hemi,parms.outstem,area_masks(a).name,parms.ico);
          if ~exist(fname_in,'file'), warning('%s not found',fname_in); end;
          if ~exist(fname_out,'file') || parms.forceflag
            cmd = ['mri_label2label --srclabel ' fname_in ...
                   ' --trglabel ' fname_out ...
                   ' --regmethod surface'...
                   ' --hemi ' hemi ...
                   ' --srcsubject ' subj ...
                   ' --trgsubject ico'...
                   ' --trgicoorder ' num2str(parms.ico)];
            [status,result] = unix(cmd);
            if status
              fprintf('%s: WARNING: cmd %s failed:\n%s\n',cmd,result);
            end;
          end
        end;
      end;
    end;
  else
    rc_retfit_func(surf,data,mask,hemi,args{:});
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
