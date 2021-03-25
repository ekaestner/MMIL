function [fit_data,area_masks,fit_parms,data] = rc_retfit_func(surf,data,roi,hemi,varargin)
%function [fit_data,area_masks,fit_parms,data] = rc_retfit_func(surf,data,roi,hemi,[options])
%
% Required Parameters:
%  surf: freesurfer sphere surface structure containing:
%    nverts: number of vertices
%    nfaces: number of faces (triangles)
%    faces:  vertex numbers for each face (3 corners)
%    vertices: x,y,z coords for each vertex (from sphere surface)
%    vertices_wm:  x,y,z coords for each vertex 
%                               (from folded white matter surface)
%     (see fs_load_subj and fs_read_surf)
%  data: structure containing:
%    pol_r: real component of polar angle map
%    pol_i: imaginary component of polar angle map
%    ecc_r: real component of eccentricity map
%    ecc_i: imaginary component of eccentricity map
%      these should be vectors of values for each vertex (can be sparse)
%      (see fs_read_wfile or fs_load_mgh)
%  roi: vector of vertex numbers defining region of interest
%     including retinotopic map(s) to be fit
%  hemi: cortical hemisphere ('lh' or 'rh')
%
% Optional Control Parameters:
%  'prereg_nruns': number of runs of pre-registration
%    if greater than 1,  will randomly choose seed point from parameter ranges
%    {default = 1}
%  'prereg_niter': pre-registration iterations
%    (random search with global map parameters)
%    {default = 1000}
%  'prereg_search_type':  string specifying which type of search to do
%    allowed values: 'rand', 'fmincon'
%    'rand' = random steps
%    'fmincon' = matlab's constrained nonlinear optimization
%    {default = 'rand'}
%  'nruns': number of runs of fine fitting
%    (updated fit parameters saved after each run)
%    {default = 1}
%  'run_reset_flag': [0|1] whether to set du and dv to 0 at start of each 
%    fine fitting run (and update u0 and v0)
%    {default = 0}
%  'niter': maximum number of iterations within nonlinear fit
%    {default = 100}
%  'step_type': string specifying which type of stepping to do
%    allowed values: 'rand' , 'grad', 'msgrad'
%    'rand' = random steps
%    'grad' = gradient descent (scaled to max_step_size)
%    'msgrad' = multi-start gradient descent
%    {default = 'rand'}
%  'fit_parms': fit parameters
%     parameters for initial fit or fit parameters from earlier runs
%    {default = []}
%  'outdir': output directory
%    {default = pwd}
%  'outstem': file stem for output mat files
%    fit_parms matfile is ignored if fit_parms supplied as parameter
%    {default = 'retfit'}
%  'step_size': standard deviation of random step (in grid units)
%    (also used in calculating gradients)
%    {default = 1}
%  'max_step_size': maximum size of random or gradient step (in grid units)
%    {default = 10}
%  'prereg_step_size': random step size for pre-registration
%    specify a vector of sizes to do multi-scale search
%    (relative to ranges for each map parameter below)
%    {default = 0.1}
%  'ms_trans': standard deviation of random translation of entire map grid
%    (units of data grid nodes)
%    (used for step_type='msgrad' when cost increases)
%    {default = 1}
%  'ms_scale': standard deviation of random scale factor applied to map
%    (used for step_type='msgrad' when cost increases)
%    {default = 0.1}
%  'max_nsteps_inc_cost': total number of steps allowed with increasing cost
%    (relevant for step_type='grad' or 'msgrad')
%    {default = 1000}
%  'area_smooth_steps': number of sparse smoothing steps to apply to
%    area masks
%    {default = 0}
%  'plot_flag': [0|1|2|3] whether to plot images
%     0: no plots
%     1: plots of fit results
%     2: plots of initial registration (pauses before fitting) and fit results
%     3: plots of initial and fitted registration and more
%     {default = 1}
%
% Optional Data Parameters:
%  'data_r_min': minimum radius (degrees visual angle)
%     used for eccentricity mapping - determines phase for a given eccentricity
%    {default = 0.25}
%  'data_r_max': maximum radius (degrees visual angle)
%     used for eccentricity mapping - determines phase for a given eccentricity
%    {default = 12.5}
%  'data_logtrans_flag': [0|1] whether eccentricity mapping data was collected
%    with log transform
%    {default = 0}
%  'interp_method': interpolation method for resampling data to grid
%    allowed values: 'linear', 'cubic', 'nearest'
%    {default = 'linear'}
%  'data_smooth_sigma': smoothing kernel applied resampling data to grid
%    Relative to unit coordinates, so values > 0.1 are very large
%    {default = 0.05}
%  'grad_smooth_sigma': smoothing kernel applied to gradients of data
%    Relative to unit coordinates, so values > 0.1 are very large
%    {default = 0.05}
%  'data_gridsize_r': number of different eccentricity data grid nodes
%    {default = 500}
%  'data_gridsize_th': number of different polar angle data grid nodes
%    {default = 500}
%
% Optional ROI Parameters
%   For initial resampling of roi to unit grid coordinates (u and v)
%   Entire V1-V2-V3 complex must fit unit grid (0 - 1)
%   Ideally, horizontal meridian of V1 should be lined up with u=0.5
%     and mid eccentricity should line up with v=0.5
%  'roi_dilate_niters': number of smoothing iterations to dilate ROI
%    this will create a border around map fit area
%    {default = 0}
%  'roi_rotation': degrees of clockwise rotation
%     {default = 0}
%  'roi_scale_u': horizontal scaling factor
%     {default = 1}
%  'roi_scale_v': vertical scaling factor
%     {default = 1}
%  'roi_shift_u': horizontal shift on unit grid (after rotation and scale)
%     {default = 0}
%  'roi_shift_v': vertical shift on unit grid (after rotation and scale)
%     {default = 0}
%  'roi_excl': vector of vertex numbers to be excluded from error calculations
%     {default = []}
%
% Optional Map Template Parameters:
%  'map_v123_flag': [0|1] whether to model V1-V2-V3 complex or
%     a single area mapping entire hemifield
%     {default = 1}
%  'map_area_name': area label if map_v123_flag=0
%     {default = 'v'}
%  'map_poly_flag': [0|1] use polynomial function to deform template
%     {default = 1}
%  'map_poly_order': order of polynomial function (n+1 additional parameters)
%      used to deform template
%     {default = 4}
%  'map_model_type': [0|1|2] model used for initial estimates of u and v
%    0: rectangle
%    1: wedge
%    2: radial wedge
%    {default = 2}
%  'map_rev_polar_flag': [0|1] whether to reverse direction of polar angle
%     {default = 0}
%  'map_v1p_width': width (along theta axis) of v1+ map
%     {default = 1}
%  'map_v1m_width': width (along theta axis) of v1- map
%     {default = 1}
%  'map_v2p_width': width (along theta axis) of v2+ map
%     {default = 1}
%  'map_v2m_width': width (along theta axis) of v2- map
%     {default = 1}
%  'map_v3p_width': width (along theta axis) of v3+ map
%     {default = 1}
%  'map_v3m_width': width (along theta axis) of v3- map
%     {default = 1}
%  'map_v1p_length': length (along ecc axis) of v1+ map
%     {default = 1}
%  'map_v1m_length': length (along ecc axis) of v1- map
%     {default = 1}
%  'map_v2p_length': length (along ecc axis) of v2+ map
%     {default = 1}
%  'map_v2m_length': length (along ecc axis) of v2- map
%     {default = 1}
%  'map_v3p_length': length (along ecc axis) of v3+ map
%     {default = 1}
%  'map_v3m_length': length (along ecc axis) of v3- map
%     {default = 1}
%  'map_poly_coef': vector of polynomial function coefficients
%     must have number of elements equal to map_poly_order + 1
%     {default = zeros(1,map_poly_order+1)}
%  'map_wedge_fact': narrowness of map's foveal edge is compared to peripheral
%    Only applies for model types 1 and 2
%    {default = 1}
%  'map_radial_wedge_fact': how horizontal extent of map translates to polar
%    angle (smaller value means less bending)
%    Only apples for model type 2
%    {default = 0.2}
%  'map_radial_offset': value added to vertical map coordinates before
%    radial wedge transform (smaller value means sharper point)
%    Only apples for model type 2
%    {default = 1}
%  'map_scale_u': horizontal scaling factor for initial u coordinates
%     for template map
%     {default = 0.7}
%  'map_scale_v': vertical scaling factor for initial v coordinates
%     for template map
%     {default = 0.4}
%  'map_rotation': degrees of clockwise rotation for initial template map
%     {default = 0}
%  'map_shift_u': horizontal shift on unit grid for initial template map
%     {default = 0}
%  'map_shift_v': vertical shift on unit grid for initial template map
%     {default = 0}
%  'map_r_min': minimum radius (degrees visual angle) in template map
%    {default = 2}
%  'map_r_max': maximum radius (degrees visual angle) in template map
%    {default = 12}
%  'map_logtrans_flag': [0|1] whether to use log transform for
%    eccentricity map grid spacing
%    {default = 1}
%  'map_gridsize_r': number of different eccentricity map grid nodes
%    {default = 10}
%  'map_gridsize_th': number of different polar map grid nodes
%    (for each visual area)
%    {default = 10}
%
% Optional Map Template Preregistration Parameters:
%  'map_v1_width_range': vector of lower and upper bounds for map parameter
%     {default = [1 1]}
%  'map_v2_width_range'
%     {default = [0.6 1]}
%  'map_v3_width_range'
%     {default = [0.5 1]}
%  'map_v1_length_range'
%     {default = [0.8 1.2]}
%  'map_v2_length_range'
%     {default = [0.8 1.2]}
%  'map_v3_length_range'
%     {default = [0.8 1.2]}
%  'map_poly_coef_range'
%     {default = [-5,5]}
%  'map_wedge_fact_range'
%    {default = [0.7 1.3]}
%  'map_radial_wedge_fact_range'
%    {default = [0.05 0.4]}
%  'map_radial_offset_range'
%    {default = [1 4]}
%  'map_scale_u_range'
%     {default = [0.6 0.8]}
%  'map_scale_v_range'
%     {default = [0.3 0.6]}
%  'map_rotation_range'
%     {default = [-10 10]}
%  'map_shift_u_range'
%     {default = [-0.1,0.1]}
%  'map_shift_v_range'
%     {default = [-0.1,0.1]}
%  'map_r_min_range'
%    {default = [1 1]}
%  'map_r_max_range'
%    {default = [12 12]}
%
% Optional Cost Function Parameters:
%  'err_fact': scalar multiplier for fit error constraint in cost function
%    {default = 1}
%  'ecc_fact': scalar multiplier for eccentricity error
%    {default = 1}
%  'smooth_fact': scalar multiplier for smoothness constraint in cost function
%    weighting relative to fit error
%    {default = 1}
%  'fold_fact': scalar cost multiplier for folding constraint in cost function
%    negative smoothness determinant
%    {default = 15}
%  'vacancy_fact': scalar cost multiplier for vacancy penaly in cost function
%     weighting factor for fraction of unoccupied grid nodes
%     applies only for prereg
%    {default = 0}
%  'max_outbound_penalty': maximum penalty for a map node being outside data ROI
%    penalty will range from 0 at ROI boundary to the max penalty at edge
%    of data grid
%    {default = 20}
%  'cost_include_percentile': map nodes with error cost greater than this
%    percentile will not contribute to error cost for fine scale fit
%    {default = 100}
%  'deformation_flag': calculate deformation between flat and folded surfaces
%    {default = 0}
%
% Output:
%  fit_data: structure containing:
%    u: u coordinate for each vertex in roi (sparse)
%    v: v coordinate for each vertex in roi
%    th: polar angle for each vertex in roi (sparse vector)
%    r: eccentricity for each vertex in roi
%    pol_r: real component of polar angle
%    pol_i: imaginary component of polar angle
%    ecc_r: real component of eccentricity
%    ecc_i: imaginary component of eccentricity
%    cost: minimum cost for each vertex in roi
%  area_masks: struct array containing:
%    name: name of each visual area
%    vertices: list of vertices for each visual area
%  fit_parms: structure containing:
%    u0: initial u coordinate for each map node
%    v0: initial v coordinate for each map node
%    du: change in u coordinate for each map node
%    dv: change in u coordinate for each map node
%    u: estimated u coordinate for each map node
%    v: estimated v coordinate for each map node
%    cost: minimum cost for each map node
%    cost_err: fit error part of cost
%    cost_smooth_u: u smoothness part of cost
%    cost_smooth_v: v smoothness part of cost
%    cost_fold: cost for folding (negative determinant of du and dv)
%    cost_outbound: out of bounds penalty
%  data: structure containing:
%    th: polar angle for each vertex in roi (sparse vector)
%    r: eccentricity for each vertex in roi
%    pol_r: real component of polar angle
%    pol_i: imaginary component of polar angle
%    ecc_r: real component of eccentricity
%    ecc_i: imaginary component of eccentricity
%
% Created:  09/18/08 by Don Hagler
% Last Mod: 12/04/13 by Don Hagler
%

%% TODO: allow lsqnonlin as prereg_search_type
%  'prereg_search_type':  string specifying which type of search to do
%    allowed values: 'rand', 'fmincon', 'lsqnonlin'
%    'rand' = random steps
%    'fmincon' = matlab's constrained nonlinear optimization
%    'lsqnonlin' = matlab's least squares nonlinear optimization
%    {default = 'rand'}

fit_data = [];
area_masks = [];
fit_parms = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

% initialize parameters, parse input
parms = init_parms(hemi,varargin);

% check that surf and data match
check_surf_data(surf,data);

% initialize fit parameters
fit_parms = init_fit_parms(parms);

% prepare data grid
[datagrid,roi,parms] = prep_data(parms,data,surf,roi);

% search for best initial starting point
[mapgrid,parms] = run_coarse_fit(parms,datagrid,fit_parms);

% calculate initial cost
[mapvec,parms] = calc_init_cost(parms,datagrid,mapgrid);

% plot data, template, and error for initial parameters
if parms.plot_flag>=2
  mapgrid_full = plot_coarse_fit(parms,mapvec,datagrid);
end;

% identify map nodes with highest cost_err
if ~isfield(mapvec,'ind_exclude')
  mapvec = exclude_nodes(parms,mapvec);
end;

% small scale fit
[mapvec,fit_parms,parms] = run_fine_fit(parms,datagrid,mapvec);

if (parms.plot_flag>=1 && ~parms.plotted_results) ||...
   (nargout>0 && ~exist('mapgrid_full','var'))
  mapgrid_full = resamp_mapgrid_full(parms,mapvec,datagrid);
end;

% plot results if not already done
if parms.plot_flag && ~parms.plotted_results
  parms = plot_final_fit(parms,datagrid,mapgrid_full);
end;

if nargout==0, return; end;

% sample fitted data and area labels to vertices
[fit_data,area_masks,data] = resamp_results_surf(parms,...
                                                 mapgrid_full,roi,surf,data);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = init_parms(hemi,options)
  parms = mmil_args2parms(options, { ...
    'prereg_nruns',1,[0,Inf],...
    'prereg_niter',1000,[0,Inf],...
    'nruns',1,[0,Inf],...
    'run_reset_flag',false,[false true],...
    'niter',100,[],...
    'outdir',pwd,[],...
    'outstem','retfit',[],...
    'step_type','rand',{'rand','grad','msgrad'},...
    'prereg_search_type','rand',{'rand','fmincon'},... %% TODO: lsqnonlin
    'interp_method','linear',{'linear','cubic','nearest'},...
    'roi_dilate_niters',0,[0,1000],...
    'roi_rotation',0,[-180,180],...
    'roi_scale_u',1,[0.1,2],...
    'roi_scale_v',1,[0.1,2],...
    'roi_shift_u',0,[-1,1],...
    'roi_shift_v',0,[-1,1],...
    'roi_excl',[],[],...
    'map_v123_flag',true,[false true],...
    'map_area_name','v',[],...
    'map_rev_polar_flag',false,[false true],...
    'map_v1p_width',1,[0.01,2],...
    'map_v1m_width',1,[0.01,2],...
    'map_v2p_width',1,[0.01,2],...
    'map_v2m_width',1,[0.01,2],...
    'map_v3p_width',1,[0.01,2],...
    'map_v3m_width',1,[0.01,2],...
    'map_v1p_length',1,[0.1,10],...
    'map_v1m_length',1,[0.1,10],...
    'map_v2p_length',1,[0.1,10],...
    'map_v2m_length',1,[0.1,10],...
    'map_v3p_length',1,[0.1,10],...
    'map_v3m_length',1,[0.1,10],...
    'map_poly_flag',true,[false true],...
    'map_poly_order',4,[1,10],...
    'map_poly_coef',[],[],...
    'map_model_type',2,[0,1,2],...
    'map_wedge_fact',1,[0.01,2],...
    'map_radial_wedge_fact',0.2,[0,1],...
    'map_radial_offset',1,[0,10],...
    'map_scale_u',0.7,[0.01,1.0],...
    'map_scale_v',0.4,[0.01,1.0],...
    'map_rotation',0,[-180,180],...
    'map_shift_u',0,[-1,1],...
    'map_shift_v',0,[-1,1],...
    'map_gridsize_r',10,[2,1000],...
    'map_gridsize_th',10,[2,1000],...
    'map_r_min',2,[0,Inf],...
    'map_r_max',12,[0,Inf],...
    'map_v1_width_range',[1.0 1.0],[0.01,2],...
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
    'data_gridsize_u',500,[100,1e10],...
    'data_gridsize_v',500,[100,1e10],...
    'fit_parms',[],[],...
    'plot_flag',0,[0,1,2,3],...
    'data_r_min',0.25,[0,Inf],...
    'data_r_max',12.5,[0,Inf],...
    'data_logtrans_flag',false,[false true],...
    'map_logtrans_flag',true,[false true],...
    'step_size',1,[0,100],...
    'max_step_size',10,[0,1000],...
    'prereg_step_size',0.1,[0.001,1],...
    'ms_trans',1,[0,100],...
    'ms_scale',0.1,[0,100],...
    'err_fact',1,[0,Inf],...
    'ecc_fact',1,[0,1],...
    'smooth_fact',1,[0,Inf],...
    'fold_fact',15,[0 Inf],...
    'vacancy_fact',0,[0,Inf],...
    'data_smooth_sigma',0.05,[0,1],...
    'grad_smooth_sigma',0.05,[0,1],...
    'cost_tol',1e-6,[1e-100,1e-2],...
    'check_grad_flag',false,[false true],...
    'max_nsteps_inc_cost',1000,[0 Inf],...
    'area_smooth_steps',0,[0,100],...
    'max_outbound_penalty',20,[0 10000],...
    'min_outbound_penalty',1e-5,[1e-100,1e-2],...
    'penalty_smooth_sigma',0.2,[0,1],...
    'penalty_smooth_niters',10,[0,100],...
    'area_code_smf',1e-5,[0,1],...
    'cost_include_percentile',100,[50,100],...
    'deformation_flag',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'smooth_sparse_niters',2,[],...
    'fit_data_smf',0.01,[0,1],...
    'area_smooth_sigma',0.08,[],...
    'area_smooth_thresh',0.1,[],...
    'area_smooth_sigma2',0.05,[],...
    'area_smooth_thresh2',0.98,[],...
    'vacancy_subsamp',10,[],...
  ...
    'hemi',hemi,[],...
  ...
    'prereg_tags',{...
      'map_v1p_width' 'map_v1m_width'...
      'map_v2p_width' 'map_v2m_width'...
      'map_v3p_width' 'map_v3m_width'...
      'map_v1p_length' 'map_v1m_length'...
      'map_v2p_length' 'map_v2m_length'...
      'map_v3p_length' 'map_v3m_length'...
      'map_wedge_fact' 'map_poly_coef' ...
      'map_radial_wedge_fact' 'map_radial_offset'...
      'map_scale_u' 'map_scale_v' 'map_rotation' 'map_shift_u'...
      'map_shift_v' 'map_r_min' 'map_r_max'},[],...
  }); 

  parms.grid_step_size = parms.step_size/...
          (mean([parms.data_gridsize_v,parms.data_gridsize_u]));
  parms.max_grid_step_size = parms.max_step_size/...
          (mean([parms.data_gridsize_v,parms.data_gridsize_u]));
  parms.grid_ms_trans = parms.ms_trans/...
          (mean([parms.data_gridsize_v,parms.data_gridsize_u]));

  parms.map_poly_ncoef = parms.map_poly_order + 1;
  if isempty(parms.map_poly_coef)
    parms.map_poly_coef = zeros(1,parms.map_poly_ncoef);
  else
    parms.map_poly_coef = mmil_rowvec(parms.map_poly_coef);
    if length(parms.map_poly_coef)~=parms.map_poly_ncoef
      error('map_poly_coef has wrong number of elements (%d not %d)',...
        length(parms.map_poly_coef),parms.map_poly_ncoef);
    end;
  end;

  parms = check_prereg_ranges(parms);

  parms = check_prereg_parms(parms);

  parms = define_areas(parms);

  parms.plotted_results = 0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = define_areas(parms)
  % for retinotopic map template grid
  if parms.map_v123_flag
    % V1-V2-V3 complex
    switch parms.hemi % cortical hemisphere
      case 'lh'
        parms.area_labels = {'v3-','v2-','v1-','v1+','v2+','v3+'};
        parms.area_bound_th = [-pi/2,0,-pi/2,0,pi/2,0,pi/2];
      case 'rh'
        parms.area_labels = {'v3+','v2+','v1+','v1-','v2-','v3-'};
        parms.area_bound_th = [pi/2,pi,pi/2,pi,3*pi/2,pi,3*pi/2];
    end;
  else
    % single area, entire hemifield (treat as two areas)
    parms.area_labels = {[parms.map_area_name '-'],[parms.map_area_name '+']};
    switch parms.hemi % cortical hemisphere
      case 'lh'
        parms.area_bound_th = [-pi/2,0,pi/2];
      case 'rh'
        parms.area_bound_th = [pi/2,pi,3*pi/2];
    end;
  end;
  if parms.map_rev_polar_flag
    parms.area_bound_th = fliplr(parms.area_bound_th);
    parms.area_labels = fliplr(parms.area_labels);
  end;

  % define retinotopic map areas
  parms.area_vec = [];
  for i=1:length(parms.area_labels)
    for j=1:parms.map_gridsize_th
      parms.area_vec{end+1} = parms.area_labels{i};
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_prereg_ranges(parms)
  % set v2 and v3 width and length ranges to [1 1] if map_v123_flag = 0
  if ~parms.map_v123_flag
    parms.map_v2_width_range = [1 1];
    parms.map_v3_width_range = [1 1];
    parms.map_v2_length_range = [1 1];
    parms.map_v3_length_range = [1 1];
  end;

  % set width ranges for individual p amd m halves
  parms.map_v1p_width_range = parms.map_v1_width_range;
  parms.map_v1m_width_range = parms.map_v1_width_range;
  parms.map_v2p_width_range = parms.map_v2_width_range;
  parms.map_v2m_width_range = parms.map_v2_width_range;
  parms.map_v3p_width_range = parms.map_v3_width_range;
  parms.map_v3m_width_range = parms.map_v3_width_range;

  % set length ranges for individual p amd m halves
  parms.map_v1p_length_range = parms.map_v1_length_range;
  parms.map_v1m_length_range = parms.map_v1_length_range;
  parms.map_v2p_length_range = parms.map_v2_length_range;
  parms.map_v2m_length_range = parms.map_v2_length_range;
  parms.map_v3p_length_range = parms.map_v3_length_range;
  parms.map_v3m_length_range = parms.map_v3_length_range;

  % set bounds equal for map_poly_coef if unused
  if parms.map_poly_flag==0
    parms.map_poly_coef_range = [0 0];
  end;
  % set bounds equal for map_wedge_fact if unused
  if parms.map_model_type==0
    parms.map_wedge_fact_range = [1 1];
  end;
  % set bounds equal for map_radial_wedge_fact if unused
  if parms.map_model_type<2
    parms.map_radial_wedge_fact_range = [0 0];
    parms.map_radial_offset_range = [0 0];
  end;

  % parameter to indicate which prereg_parms are variable
  parms.prereg_vars = zeros(1,length(parms.prereg_tags));
  for t=1:length(parms.prereg_tags)
    parm_range = parms.(sprintf('%s_range',parms.prereg_tags{t}));
    if range(parm_range)~=0
      parms.prereg_vars(t) = 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_surf_data(surf,data)
  % check that dimensions of data match
  if any(size(data.pol_r)~=size(data.pol_i)) ||...
     any(size(data.ecc_r)~=size(data.ecc_i)) ||...
     any(size(data.pol_r)~=size(data.ecc_r))
    error('dimensions of data (pol_r, pol_i, ecc_r, ecc_i) must match');
  end;
  nverts = length(data.pol_r);
  if nverts~=surf.nverts
    error('number of elements in data.pol_r (%d) does not match surf.nverts (%d)',...
      nverts,surf.nverts);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fit_parms = init_fit_parms(parms)
  % check for existing fit_parms
  fit_parms = [];
  if ~isempty(parms.fit_parms)
    fit_parms = parms.fit_parms;
  else
    matfile = sprintf('%s/%s_fit_parms-%s.mat',...
      parms.outdir,parms.outstem,parms.hemi);
    if exist(matfile,'file') && ~parms.forceflag, load(matfile); end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [datagrid,roi,parms] = prep_data(parms,data,surf,roi)
  % resample data to grid with interpolation
  matfile = sprintf('%s/%s_datagrid-%s.mat',...
    parms.outdir,parms.outstem,parms.hemi);
  if ~exist(matfile,'file') || parms.forceflag
    [datagrid,roi,roi_orig] = resamp_data(data,parms,surf,roi);
    save(matfile,'datagrid','roi','roi_orig');
  else
    fprintf('%s: loading existing data grid\n',mfilename);
    fprintf('\t(ignoring any changes to data and roi parameters)...\n');
    load(matfile);
  end;
  if parms.vacancy_fact > 0
    parms.sub_datagrid = subsamp_datagrid(datagrid,parms.vacancy_subsamp);
  end;

  % renormalize outbound penalty
  if parms.max_outbound_penalty>0
    if isempty(datagrid.outbound_penalty)
      error('outbound_penalty was already initialized as emptpy');
    else
      datagrid.outbound_penalty = datagrid.outbound_penalty *...
        parms.max_outbound_penalty / max(datagrid.outbound_penalty(:));
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sub_datagrid = subsamp_datagrid(datagrid,subsamp_fact)
  ind_X = subsamp_fact:subsamp_fact:length(datagrid.X);
  ind_Y = subsamp_fact:subsamp_fact:length(datagrid.Y);
  sub_datagrid = [];
  sub_datagrid.X = datagrid.X(ind_X);
  sub_datagrid.Y = datagrid.Y(ind_Y);
  sub_datagrid.origX = sub_datagrid.X;
  sub_datagrid.origY = sub_datagrid.Y;
  sub_datagrid.roi_mask = datagrid.roi_mask(ind_Y,ind_X);
  sub_datagrid.size = size(sub_datagrid.roi_mask);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mapgrid,parms] = run_coarse_fit(parms,datagrid,fit_parms)
  matfile = sprintf('%s/%s_prereg_parms-%s.mat',...
    parms.outdir,parms.outstem,parms.hemi);
  if isempty(fit_parms)
    if parms.plot_flag>=2 && ~exist(matfile,'file')
      if parms.prereg_nruns<=1
        basis_flag = 0;
      else
        basis_flag = 1;
      end;
      [beta,parms] = select_prereg_betas(parms,basis_flag);
      prereg_parms = parms2prereg_parms(parms);
      fprintf('%s: initial pre-registration parameters:\n',mfilename);
      disp(prereg_parms);
      plot_init_fit(parms,datagrid);
    else
      prereg_parms = parms2prereg_parms(parms);
    end;
    if parms.prereg_niter>0 && parms.prereg_nruns>0
      if ~exist(matfile,'file') || parms.forceflag
        fprintf('%s: large scale search with global parameters...\n',mfilename);
        tic;
        best_parms = parms;
        best_min_sumcost = Inf;
        for r=1:parms.prereg_nruns
          if parms.prereg_nruns>1 && r<parms.prereg_nruns
            if r==1
              basis_flag = 1;
            else
              basis_flag = 2;
            end;
            [beta_cur,parms] = select_prereg_betas(parms,basis_flag);
            fprintf('%s: prereg run %d initial min sum cost:\n',...
              mfilename,r);
            min_sumcost = prereg_costfunc(beta_cur,datagrid,parms,1);
          elseif parms.prereg_nruns==1
            parms = best_parms;
            beta_cur = prereg_parms2beta(parms);
            fprintf('%s: prereg initial min sum cost:\n',mfilename);
            min_sumcost = prereg_costfunc(beta_cur,datagrid,parms,1);
          elseif r==parms.prereg_nruns
            parms = best_parms;
            min_sumcost = best_min_sumcost;
            fprintf('%s: prereg run %d initial min sum cost = %0.4f\n',...
              mfilename,r,min_sumcost);
          end;
          switch parms.prereg_search_type
            case 'rand'
              min_sumcost_last = min_sumcost;
              for p=1:length(parms.prereg_step_size)
                tmp_parms = parms;
                tmp_parms.prereg_step_size = parms.prereg_step_size(p);
                for i=1:parms.prereg_niter
                  tmp2_parms = randstep_prereg_parms(tmp_parms);
                  beta_cur = prereg_parms2beta(tmp2_parms);
                  sumcost = prereg_costfunc(beta_cur,datagrid,parms);
                  if sumcost < min_sumcost
                    tmp_parms = tmp2_parms;
                    min_sumcost = sumcost;
                    fprintf('%s: iter %d: min sum cost:\n',mfilename,i);
                    sumcost = prereg_costfunc(beta_cur,datagrid,parms,1);
                  end;
                end;
                if min_sumcost < min_sumcost_last
                  tmp_parms.prereg_step_size = parms.prereg_step_size;
                  parms = tmp_parms;
                  min_sumcost_last = min_sumcost;
                  fprintf('%s: prereg scale %d: min sum cost = %0.4f\n',...
                    mfilename,p,min_sumcost);
                end;
              end;
            case 'fmincon'
              beta_init = prereg_parms2beta(parms);
              [lbounds,ubounds] = prereg_parms2bounds(parms);
              %% TODO: opts = {...}
              try % for R2009b
                fit_options = optimset('Display','iter',...
                  'Algorithm','active-set' ,...
                  'MaxFunEvals',Inf,...
                  'MaxIter',parms.prereg_niter,...
                  'DiffMinChange',0.01,'TolFun',1e-3,'TolX',1e-4);
              catch % for R2007a
                fit_options = optimset('Display','iter',...
                  'LargeScale','off',...
                  'MaxFunEvals',Inf,...
                  'MaxIter',parms.prereg_niter,...
                  'DiffMinChange',0.01,'TolFun',1e-3,'TolX',1e-4);
              end;
              [parms.beta_fit,fval,exitflag,output] = ...
                fmincon(@(beta) prereg_costfunc(beta,datagrid,parms),beta_init,...
                  [],[],[],[],...
                  lbounds,ubounds,[],fit_options);
              if exitflag~=1
                fprintf('%s: WARNING: fmincon returned with exitflag %d\n',mfilename,exitflag);
              end;
              min_sumcost = fval;
              parms = prereg_beta2parms(parms.beta_fit,parms);
            %% TODO: lsqnonlin
          end;
          if min_sumcost < best_min_sumcost
            if parms.prereg_nruns>1
              fprintf('%s: prereg run %d: min sum cost = %0.4f\n',...
                mfilename,r,min_sumcost);
            end;
            best_min_sumcost = min_sumcost;
            best_parms = parms;
          end;
        end;
        fprintf('%s: prereg min sum cost = %0.4f\n',...
          mfilename,min_sumcost);
        parms = best_parms;
        toc;
        % save prereg parms
        prereg_parms = parms2prereg_parms(parms);
        save(matfile,'prereg_parms');
      else
        % load prereg parms
        fprintf('%s: loading global parameters from previous pre-reg...\n',mfilename);
        load(matfile);
        parms = prereg_parms2parms(prereg_parms,parms);
      end;
    end;
    mapgrid = [];
  else
    fprintf('%s: using fit parameters from previous fit...\n',mfilename);
    mapgrid.u0 = fit_parms.u0;
    mapgrid.v0 = fit_parms.v0;
    mapgrid.du = fit_parms.du;
    mapgrid.dv = fit_parms.dv;
    mapgrid.u = fit_parms.u;
    mapgrid.v = fit_parms.v;
    if isfield(fit_parms,'ind_exclude')
      mapgrid.ind_exclude = fit_parms.ind_exclude;
    end;
    if exist(matfile,'file')
      % load prereg parms
      fprintf('%s: loading global parameters from previous pre-reg...\n',mfilename);
      load(matfile);
      parms = prereg_parms2parms(prereg_parms,parms);
    end;
  end;
  fprintf('%s: pre-registration parameters:\n',mfilename);
  disp(prereg_parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_init_fit(parms,datagrid)
  mapvec = init_ret_grid(parms,datagrid);
  % calculate each type of cost separately, attach to mapvec
  min_cost_err = parms.err_fact*costfunc_err(mapvec,datagrid,parms.ecc_fact);
  [cost_sm,cost_fold] = costfunc_smooth(mapvec);
  min_cost_smooth = parms.smooth_fact*cost_sm;
  min_cost_fold = parms.fold_fact*cost_fold;
  min_cost_vacancy = parms.vacancy_fact*costfunc_vacancy(mapvec,parms);
  min_cost_outbound = costfunc_outbound(mapvec,datagrid);
  mapvec.cost_err = min_cost_err(1:mapvec.num_nodes);
  mapvec.cost_smooth_u = min_cost_smooth(1:mapvec.num_nodes);
  mapvec.cost_smooth_v = min_cost_smooth(mapvec.num_nodes+1:end);
  mapvec.cost_fold = min_cost_fold(1:mapvec.num_nodes);
  mapvec.cost_outbound = min_cost_outbound(1:mapvec.num_nodes);
  mapgrid_full = resamp_mapgrid_full(parms,mapvec,datagrid);
  plot_datafit(datagrid,mapgrid_full,parms);
  plot_areas(mapgrid_full);
  plot_cost(datagrid,mapgrid_full,parms);
  fprintf('%s: press a key to start rough fitting\n',mfilename);
  drawnow;
  pause;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta,parms] = select_prereg_betas(parms,basis_flag)
  if basis_flag == 0 % use input parameters
    beta = prereg_parms2beta(parms);
  else
    [lbounds,ubounds] = prereg_parms2bounds(parms);
    beta = zeros(size(lbounds));
    if basis_flag == 1 % middle of specified ranges
      for b=1:length(beta)
        beta(b) = (lbounds(b) + ubounds(b))/2;
      end;
    elseif basis_flag == 2 % random values inside ranges
      for b=1:length(beta)
        beta(b) = lbounds(b) + rand*(ubounds(b)-lbounds(b));
      end;
    end;
    parms = prereg_beta2parms(beta,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mapvec,parms] = calc_init_cost(parms,datagrid,mapgrid)
  [mapvec,parms]=init_ret_grid(parms,datagrid,mapgrid);
  parms.beta_fit = [mapvec.du;mapvec.dv];
  min_sumcost = sum(costfunc(parms.beta_fit,datagrid,mapvec,parms));
  min_cost_err = parms.err_fact*costfunc_err(mapvec,datagrid,parms.ecc_fact);
  [cost_sm,cost_fold] = costfunc_smooth(mapvec);
  min_cost_smooth = parms.smooth_fact*cost_sm;
  min_cost_fold = parms.fold_fact*cost_fold;
  min_cost_vacancy = parms.vacancy_fact*costfunc_vacancy(mapvec,parms);
  min_cost_outbound = costfunc_outbound(mapvec,datagrid);
  mapvec.cost_err = min_cost_err(1:mapvec.num_nodes);
  mapvec.cost_smooth_u = min_cost_smooth(1:mapvec.num_nodes);
  mapvec.cost_smooth_v = min_cost_smooth(mapvec.num_nodes+1:end);
  mapvec.cost_fold = min_cost_fold(1:mapvec.num_nodes);
  mapvec.cost_vacancy = min_cost_vacancy(1:mapvec.num_nodes);
  mapvec.cost_outbound = min_cost_outbound(1:mapvec.num_nodes);
  fprintf('%s: initial min sum cost = %0.4f\n',mfilename,min_sumcost);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mapgrid_full = plot_coarse_fit(parms,mapvec,datagrid)
  mapgrid_full = resamp_mapgrid_full(parms,mapvec,datagrid);
  plot_datafit(datagrid,mapgrid_full,parms);
  plot_areas(mapgrid_full);
  plot_cost(datagrid,mapgrid_full,parms);
  fprintf('%s: press a key to start fine fitting\n',mfilename);
  drawnow;
  pause;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mapvec = exclude_nodes(parms,mapvec)
  % identify map nodes with highest cost_err
  nbins = 100;
  [N,X] = hist(mapvec.cost_err,nbins);
  P = cumsum(N); Pnorm = P/P(end);
  [tmp,ind_min] = min(abs(Pnorm-parms.cost_include_percentile/100));
  if ind_min == nbins
    mapvec.ind_exclude = [];
  else
    mapvec.ind_exclude = find(mapvec.cost_err>X(ind_min));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mapvec,fit_parms,parms] = run_fine_fit(parms,datagrid,mapvec)
  fprintf('%s: multi-dimensional small scale fit...\n',mfilename);
  if parms.check_grad_flag && ismember(parms.step_type,{'grad','msgrad'})
    check_grad(parms,datagrid,mapvec);
  end;
  parms.vacancy_fact = 0;
  fit_parms = [];
  for r=1:parms.nruns
    tic;
    fprintf('%s: ####### run %d ########\n',mfilename,r);
    if parms.run_reset_flag
      mapvec.u0 = mapvec.u0 + mapvec.du;
      mapvec.v0 = mapvec.v0 + mapvec.dv;
      mapvec.du = zeros(size(mapvec.u0));
      mapvec.dv = zeros(size(mapvec.v0));
    end;
    % nonlin fit
    beta_cur = [mapvec.du;mapvec.dv];
    [min_cost,g] = costfunc(beta_cur,datagrid,mapvec,parms);
    min_sumcost = sum(min_cost);
    beta_last = parms.beta_fit;
    fprintf('%s: run %d initial min sum cost = %0.4f\n',...
      mfilename,r,min_sumcost);  
    if r==1, min_sumcost = Inf; end;
    min_sumcost_last = min_sumcost;
    nsteps_inc_cost = 0;
    for i=1:parms.niter
      if strcmp(parms.step_type,'rand')
        beta_cur = rand_step(parms.beta_fit,...
          parms.grid_step_size,parms.max_grid_step_size);
        cost = costfunc(beta_cur,datagrid,mapvec,parms);
      else
        beta_cur = grad_step(beta_last,g,parms.max_grid_step_size);
        [cost,g] = costfunc(beta_cur,datagrid,mapvec,parms);
      end;
      sumcost = sum(cost);
  %    fprintf('%s: iter %d: norm cost = %0.4f\n',mfilename,i,sumcost);  
      if (sumcost < min_sumcost)
        min_cost = cost;
        min_sumcost = sumcost;
        parms.beta_fit = beta_cur;
        beta_last = parms.beta_fit;
        mapvec.du = parms.beta_fit(1:mapvec.num_nodes);
        mapvec.dv = parms.beta_fit(mapvec.num_nodes+1:end);
        mapvec.u = mapvec.u0 + mapvec.du;
        mapvec.v = mapvec.v0 + mapvec.dv;
        min_cost_err = ...
          parms.err_fact*costfunc_err(mapvec,datagrid,parms.ecc_fact);
        [cost_sm,cost_fold] = costfunc_smooth(mapvec);
        min_cost_smooth = parms.smooth_fact*cost_sm;
        min_cost_fold = parms.fold_fact*cost_fold;
        min_cost_vacancy = costfunc_vacancy(mapvec,parms);
        min_cost_outbound = costfunc_outbound(mapvec,datagrid);
        sumcost_err = sum(min_cost_err);
        sumcost_smooth = sum(min_cost_smooth);
        sumcost_fold = sum(min_cost_fold);
        sumcost_vacancy = sum(min_cost_vacancy);
        sumcost_outbound = sum(min_cost_outbound);
        fprintf('%s: iter %d: min sum cost = %0.4f\n',...
          mfilename,i,min_sumcost);
        fprintf('    err = %0.4f, smooth = %0.4f, fold = %0.4f, outbound = %0.4f\n',...
          sumcost_err,sumcost_smooth,sumcost_fold,sumcost_outbound);  
      elseif strcmp(parms.step_type,'grad')
  %      fprintf('%s: iter %d: cost increased to %0.4f\n',...
  %        mfilename,i,sumcost);
        parms.max_grid_step_size = parms.max_grid_step_size/2;
        nsteps_inc_cost = nsteps_inc_cost + 1;
        if nsteps_inc_cost > parms.max_nsteps_inc_cost
          break;
        end;
      elseif strcmp(parms.step_type,'msgrad')
  %      fprintf('%s: iter %d: cost increased to %0.4f\n',...
  %        mfilename,i,sumcost);
        beta_last = (1+parms.ms_scale*randn)*...
                   (parms.beta_fit + parms.grid_ms_trans*randn);
        %% todo: rotation?
        [cost,g] = costfunc(beta_last,datagrid,mapvec,parms);
        nsteps_inc_cost = nsteps_inc_cost + 1;
        if nsteps_inc_cost > parms.max_nsteps_inc_cost
          break;
        end;
      end;
      if (min_sumcost_last-min_sumcost<parms.cost_tol) &&...
         (i>parms.max_nsteps_inc_cost || strcmp(parms.step_type,'grad'))
        fprintf('%s: difference in cost (%0.1e) less than tolerance (%0.1e)\n',...
          mfilename,min_sumcost_last-min_sumcost,parms.cost_tol);
        break;
      end;
    end;

    fprintf('%s: ## run %d: final min sum cost = %0.4f ##\n',...
      mfilename,r,min_sumcost);
    toc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if min_sumcost < min_sumcost_last || isempty(fit_parms)
                                                %|| ~strcmp(parms.step_type,'msgrad')
      % reshape deformed map vec to grid
      mapvec.du = parms.beta_fit(1:mapvec.num_nodes);
      mapvec.dv = parms.beta_fit(mapvec.num_nodes+1:end);
      mapvec.u = mapvec.u0 + mapvec.du;
      mapvec.v = mapvec.v0 + mapvec.dv;
      mapvec.cost_err = min_cost_err(1:mapvec.num_nodes);
      mapvec.cost_smooth_u = min_cost_smooth(1:mapvec.num_nodes);
      mapvec.cost_smooth_v = min_cost_smooth(mapvec.num_nodes+1:end);
      mapvec.cost_fold = min_cost_fold(1:mapvec.num_nodes);
      mapvec.cost_outbound = min_cost_outbound(1:mapvec.num_nodes);
      fit_parms.u0 = reshape(mapvec.u0,mapvec.grid_size);
      fit_parms.v0 = reshape(mapvec.v0,mapvec.grid_size);
      fit_parms.du = reshape(mapvec.du,mapvec.grid_size);
      fit_parms.dv = reshape(mapvec.dv,mapvec.grid_size);
      fit_parms.u = reshape(mapvec.u,mapvec.grid_size);
      fit_parms.v = reshape(mapvec.v,mapvec.grid_size);
      fit_parms.cost_err = reshape(mapvec.cost_err,mapvec.grid_size);
      fit_parms.cost_smooth_u = reshape(mapvec.cost_smooth_u,mapvec.grid_size);
      fit_parms.cost_smooth_v = reshape(mapvec.cost_smooth_v,mapvec.grid_size);
      fit_parms.cost_fold = reshape(mapvec.cost_fold,mapvec.grid_size);
      fit_parms.cost_outbound = reshape(mapvec.cost_outbound,mapvec.grid_size);
      fit_parms.ind_exclude = mapvec.ind_exclude;

      % save fit parameters
      matfile = sprintf('%s/%s_fit_parms-%s.mat',...
        parms.outdir,parms.outstem,parms.hemi);
      save(matfile,'fit_parms');

      if parms.plot_flag>=1
        mapgrid_full = resamp_mapgrid_full(parms,mapvec,datagrid);
        fprintf('%s: plotting fit results...\n',mfilename);
        plot_datafit(datagrid,mapgrid_full,parms);
        if parms.plot_flag>=3
          plot_datafit_sub(datagrid,mapgrid_full,parms);
        end;
        plot_cost(datagrid,mapgrid_full,parms);
        plot_areas(mapgrid_full);
        parms.plotted_results = 1;
      end;
    end;
    min_sumcost_last = min_sumcost;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_final_fit(parms,datagrid,mapgrid_full)
  fprintf('%s: plotting fit results...\n',mfilename);
  if parms.plot_flag>=2
    plot_datafit(datagrid,mapgrid_full,parms);
  end;
  if parms.plot_flag>=3
    plot_datafit_sub(datagrid,mapgrid_full,parms);
  end;
  plot_cost(datagrid,mapgrid_full,parms);
  plot_areas(mapgrid_full);
  parms.plotted_results = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [datagrid,roi,roi_orig] = resamp_data(data,parms,surf,roi)
  roi_orig = roi;
  if parms.roi_dilate_niters
    fprintf('%s: dilating ROI with %d iterations...\n',...
      mfilename,parms.roi_dilate_niters);
    roi = dilate_roi(roi,surf,parms.roi_dilate_niters);
  end;

  fprintf('%s: transforming spherical coordinates to unit cartesian...\n',...
    mfilename);
  coords_uv = transform_coords(surf,roi,roi_orig,...
    parms.roi_rotation,...
    parms.roi_shift_u,parms.roi_shift_v,...
    parms.roi_scale_u,parms.roi_scale_v);
  X = coords_uv(:,1);
  Y = coords_uv(:,2);
  [tmp,ind_roi_orig] = intersect(roi,roi_orig);
  X_roi = X(ind_roi_orig);
  Y_roi = Y(ind_roi_orig);
  datagrid = [];
  datagrid.X = ((1:parms.data_gridsize_u)-1)/parms.data_gridsize_u;
  datagrid.Y = (((1:parms.data_gridsize_v)-1)/parms.data_gridsize_v)';
  datagrid.origX = X;
  datagrid.origY = Y;
  datagrid.origX_roi = X_roi;
  datagrid.origY_roi = Y_roi;

  % sample original region of interest to grid
  fprintf('%s: sampling ROI to grid (%s interpolation)...\n',...
    mfilename,parms.interp_method);
  mask = zeros(size(roi));
  mask(ind_roi_orig) = 1;
  datagrid.roi_mask = griddata(X,Y,mask,...
    datagrid.X,datagrid.Y,parms.interp_method);
  datagrid.roi_mask = 1.0*(datagrid.roi_mask>0.5);
  datagrid.size = size(datagrid.roi_mask);

  % create mask of region, excluding defined region
  if ~isempty(parms.roi_excl)
    fprintf('%s: excluding defined region...\n',mfilename);
    [tmp,ind_excl] = intersect(roi,parms.roi_excl);
    mask = ones(size(roi));
    if ~isempty(ind_excl)
      mask(ind_excl) = 0;
    end;
    datagrid.err_mask = griddata(X,Y,mask,...
      datagrid.X,datagrid.Y,parms.interp_method);
    datagrid.err_mask = 1.0*(datagrid.err_mask>0.5);
  else
    datagrid.err_mask = ones(datagrid.size);
  end;

  if parms.max_outbound_penalty>0
    fprintf('%s: initializing boundary penalty...\n',mfilename);
    sig_u = parms.penalty_smooth_sigma*parms.data_gridsize_u;
    sig_v = parms.penalty_smooth_sigma*parms.data_gridsize_v;
    datagrid.outbound_penalty = zeros(datagrid.size);
    ind_grid_roi = find(datagrid.roi_mask);
    for i=1:parms.penalty_smooth_niters
      datagrid.outbound_penalty(ind_grid_roi) = parms.max_outbound_penalty;
      datagrid.outbound_penalty = mmil_smooth2d(datagrid.outbound_penalty,...
        sig_u,sig_v);
    end;
    datagrid.outbound_penalty = ...
      parms.max_outbound_penalty - datagrid.outbound_penalty;
    datagrid.outbound_penalty(datagrid.outbound_penalty<...
      parms.min_outbound_penalty) = 0;
    fprintf('%s: calculating gradients of boundary penalty...\n',mfilename);
    [datagrid.dpenalty_du,datagrid.dpenalty_dv] =...
      gradient(datagrid.outbound_penalty,...
      1/parms.data_gridsize_u,1/parms.data_gridsize_v);
  else
    datagrid.outbound_penalty = [];
  end;

  % polar angle
  fprintf('%s: sampling pol data to grid (%s interpolation)...\n',...
    mfilename,parms.interp_method);
  pol_r = griddata(X,Y,double(full(data.pol_r(roi))),...
    datagrid.X,datagrid.Y,parms.interp_method);
  pol_i = griddata(X,Y,double(full(data.pol_i(roi))),...
    datagrid.X,datagrid.Y,parms.interp_method);
  if parms.data_smooth_sigma
    fprintf('%s: smoothing data with sigma = %0.2f...\n',...
      mfilename,parms.data_smooth_sigma);
    ind_bad_nodes = find(isnan(pol_r) | isnan(pol_i));
    pol_r(ind_bad_nodes)=0;
    pol_i(ind_bad_nodes)=0;
    sig_u = parms.data_smooth_sigma*parms.data_gridsize_u;
    sig_v = parms.data_smooth_sigma*parms.data_gridsize_v;
    pol_r = mmil_smooth2d(pol_r,sig_u,sig_v);
    pol_i = mmil_smooth2d(pol_i,sig_u,sig_v);
    pol_r(ind_bad_nodes)=0;
    pol_i(ind_bad_nodes)=0;
  end;

  % eccentricity
  fprintf('%s: sampling ecc data to grid (%s interpolation)...\n',...
    mfilename,parms.interp_method);
  ecc_r = griddata(X,Y,double(full(data.ecc_r(roi))),...
    datagrid.X,datagrid.Y,parms.interp_method);
  ecc_i = griddata(X,Y,double(full(data.ecc_i(roi))),...
    datagrid.X,datagrid.Y,parms.interp_method);
  if parms.data_smooth_sigma
    fprintf('%s: smoothing data with sigma = %0.2f...\n',...
      mfilename,parms.data_smooth_sigma);
    ind_bad_nodes = find(isnan(ecc_r) | isnan(ecc_i));
    ecc_r(ind_bad_nodes)=0;
    ecc_i(ind_bad_nodes)=0;
    sig_u = parms.data_smooth_sigma*parms.data_gridsize_u;
    sig_v = parms.data_smooth_sigma*parms.data_gridsize_v;
    ecc_r = mmil_smooth2d(ecc_r,sig_u,sig_v);
    ecc_i = mmil_smooth2d(ecc_i,sig_u,sig_v);
    ecc_r(ind_bad_nodes)=0;
    ecc_i(ind_bad_nodes)=0;
  end;

  % create mask of grid points with data
  datagrid.data_mask = ones(datagrid.size);
  datagrid.data_mask((pol_r==0 & pol_i==0) |...
                     ~isfinite(atan2(pol_i,pol_r))) = 0;

  % calculate theta from r and i
  datagrid.th = atan2(pol_i,pol_r).*datagrid.data_mask;
  datagrid.pol_r = cos(datagrid.th).*datagrid.data_mask;
  datagrid.pol_i = sin(datagrid.th).*datagrid.data_mask;

  % calculate standard deviation of polar angle data
  datagrid.std_pol_r = std(datagrid.pol_r(datagrid.data_mask==1));
  datagrid.std_pol_i = std(datagrid.pol_i(datagrid.data_mask==1));
  % make sure that they are not too small
  datagrid.std_pol_r = max(datagrid.std_pol_r,datagrid.std_pol_i/2);
  datagrid.std_pol_i = max(datagrid.std_pol_i,datagrid.std_pol_r/2);
  
  % calculate eccentricity phase from r and i
  datagrid.ecc_phase = atan2(ecc_i,ecc_r)/(2*pi);
  datagrid.ecc_phase(datagrid.ecc_phase<0) =...
    datagrid.ecc_phase(datagrid.ecc_phase<0)+1;
  % convert eccentricity phase to deg visual angle
  if parms.data_logtrans_flag
    datagrid.r = parms.data_r_min*...
      exp(datagrid.ecc_phase*log(parms.data_r_max/parms.data_r_min));
  else
    datagrid.r = parms.data_r_min + ...
      (parms.data_r_max - parms.data_r_min)*datagrid.ecc_phase;
  end;
  datagrid.r(~isfinite(datagrid.r)) = 0;
  tmp_r = (datagrid.r - parms.data_r_min)*...
    2*pi/(parms.data_r_max - parms.data_r_min);
  datagrid.ecc_r = cos(tmp_r);
  datagrid.ecc_i = sin(tmp_r);
  datagrid.ecc_r(tmp_r<=0) = 0;
  datagrid.ecc_i(tmp_r<=0) = 0;

  % calculate standard deviation of eccentricity data
  datagrid.std_ecc_r = std(datagrid.ecc_r(datagrid.data_mask==1));
  datagrid.std_ecc_i = std(datagrid.ecc_i(datagrid.data_mask==1));
  % make sure that they are not too small
  datagrid.std_ecc_r = max(datagrid.std_ecc_r,datagrid.std_ecc_i/2);
  datagrid.std_ecc_i = max(datagrid.std_ecc_i,datagrid.std_ecc_r/2);

  % calculate gradients of pol_r and pol_i wrt u and v
  if ismember(parms.step_type,{'rand','grad','msgrad'})
    fprintf('%s: calculating pol gradients...\n',mfilename);
    [datagrid.dpr_du,datagrid.dpr_dv] =...
      gradient(datagrid.pol_r,...
      1/parms.data_gridsize_u,1/parms.data_gridsize_v);
    [datagrid.dpi_du,datagrid.dpi_dv] =...
      gradient(datagrid.pol_i,...
      1/parms.data_gridsize_u,1/parms.data_gridsize_v);
    if parms.grad_smooth_sigma
      fprintf('%s: smoothing gradients with sigma = %0.2f...\n',...
        mfilename,parms.grad_smooth_sigma);
      sig_u = parms.grad_smooth_sigma*parms.data_gridsize_u;
      sig_v = parms.grad_smooth_sigma*parms.data_gridsize_v;
      datagrid.dpr_du = mmil_smooth2d(datagrid.dpr_du,sig_u,sig_v);
      datagrid.dpr_dv = mmil_smooth2d(datagrid.dpr_dv,sig_u,sig_v);
      datagrid.dpi_du = mmil_smooth2d(datagrid.dpi_du,sig_u,sig_v);
      datagrid.dpi_dv = mmil_smooth2d(datagrid.dpi_dv,sig_u,sig_v);
    end;
  end;

  % calculate gradients of ecc_r and ecc_i wrt u and v
  if ismember(parms.step_type,{'rand','grad','msgrad'})
    fprintf('%s: calculating ecc gradients...\n',mfilename);
    [datagrid.der_du,datagrid.der_dv] =...
      gradient(datagrid.ecc_r,...
      1/parms.data_gridsize_u,1/parms.data_gridsize_v);
    [datagrid.dei_du,datagrid.dei_dv] =...
      gradient(datagrid.ecc_i,...
      1/parms.data_gridsize_u,1/parms.data_gridsize_v);
    if parms.grad_smooth_sigma
      fprintf('%s: smoothing gradients with sigma = %0.2f...\n',...
        mfilename,parms.grad_smooth_sigma);
      sig_u = parms.grad_smooth_sigma*parms.data_gridsize_u;
      sig_v = parms.grad_smooth_sigma*parms.data_gridsize_v;
      datagrid.der_du = mmil_smooth2d(datagrid.der_du,sig_u,sig_v);
      datagrid.der_dv = mmil_smooth2d(datagrid.der_dv,sig_u,sig_v);
      datagrid.dei_du = mmil_smooth2d(datagrid.dei_du,sig_u,sig_v);
      datagrid.dei_dv = mmil_smooth2d(datagrid.dei_dv,sig_u,sig_v);
    end;
  end;

  if parms.deformation_flag
    % calculate deformation for each vertex in ROI
    datagrid.deformation = calc_deformation(surf,roi,coords_uv);

    fprintf('%s: sampling deformation det to grid (%s interpolation)...\n',...
      mfilename,parms.interp_method);
    datagrid.defdet = griddata(X,Y,datagrid.deformation.det,...
      datagrid.X,datagrid.Y,parms.interp_method);
    datagrid.defdet(isnan(datagrid.defdet))=0;

    datagrid.dudxy = griddata(X,Y,datagrid.deformation.dudxy,...
      datagrid.X,datagrid.Y,parms.interp_method);
    datagrid.dudxy(isnan(datagrid.dudxy))=0;

    datagrid.dvdxy = griddata(X,Y,datagrid.deformation.dvdxy,...
      datagrid.X,datagrid.Y,parms.interp_method);
    datagrid.dvdxy(isnan(datagrid.dvdxy))=0;

    if 1
      figure(1); clf;
      imagesc(datagrid.defdet,[0,0.1]); colorbar;
      title('deformation determinant');

      figure(2); clf;
      imagesc(datagrid.dudxy); colorbar;
      title('deformation dudxy');

      figure(3); clf;
      imagesc(datagrid.dvdxy); colorbar;
      title('deformation dvdxy');

      fprintf('%s: smoothing gradients with sigma = %0.2f...\n',...
        mfilename,parms.data_smooth_sigma);
      sig_u = parms.data_smooth_sigma*parms.data_gridsize_u;
      sig_v = parms.data_smooth_sigma*parms.data_gridsize_v;
      datagrid.dudxy = mmil_smooth2d(datagrid.dudxy,sig_u,sig_v);
      datagrid.dvdxy = mmil_smooth2d(datagrid.dvdxy,sig_u,sig_v);

      figure(4); clf;
      imagesc(datagrid.dudxy); colorbar;
      title('smoothed deformation dudxy');

      figure(5); clf;
      imagesc(datagrid.dvdxy); colorbar;
      title('smoothed deformation dvdxy');

      keyboard
    end;
  end;

  if parms.plot_flag>=3
    plot_data(datagrid,parms);
    if ismember(parms.step_type,{'grad','msgrad'})
      plot_gradients(datagrid,parms);
    end;
  end;

  clear pol_r pol_i ecc_r ecc_i;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deformation=calc_deformation(surf,roi,coords_uv);
  fprintf('%s: calculating deformations...\n',mfilename);
  tic
  
  % initialize output
  deformation = [];
  nroi = length(roi);
  deformation.dudx = zeros(nroi,1);
  deformation.dudy = zeros(nroi,1);
  deformation.dvdx = zeros(nroi,1);
  deformation.dvdy = zeros(nroi,1);
  deformation.det  = zeros(nroi,1);

  % consistently reorders neighbor vertices
  msurf = [];
  msurf.vertices = surf.vertices_wm;
  msurf.faces = surf.faces;
  msurf = preprocessQ(msurf);

  % calculate deformation for each vertex in ROI
  for k=1:nroi
    v=roi(k);
    if (msurf.NNoV(v)~=msurf.NTaV(v)) || (msurf.NNoV(v)<=1)
      continue;
    end;

    % get vertex numbers for neighbors
    num_nbrs = msurf.NNoV(v);
    v_nbrs = msurf.Neighbors(msurf.Index2(v) + [0:num_nbrs-1]);

    % calculate distances between neighbors and center
    dists = vertexDist(msurf,v*ones(num_nbrs,1),v_nbrs);

    % calculate angles between each neighbor
    angles = vertexAngle(msurf,v_nbrs,v*ones(num_nbrs,1),[v_nbrs(2:end);v_nbrs(1)]);

    % scale angles so total is 2*pi
    angles = angles*2*pi/sum(angles);

    % calculate new locally flat coordinates
    cum_angles = cumsum(angles)-angles(1);
    coords_flat = [dists.*cos(cum_angles),dists.*sin(cum_angles)];

    % calculate partial derivatives (adapted from fieldsign code)
    m11=0; m12=0; m13=0; m22=0; m23=0; z1=0; z2=0; z3=0; z1b=0; z2b=0; z3b=0;
    for m=1:num_nbrs
      vn = v_nbrs(m);
      j = find(vn==roi);
      if isempty(j), continue; end;
      du = coords_uv(j,1) - coords_uv(k,1);
      dv = coords_uv(j,2) - coords_uv(k,2);
      dx = coords_flat(m,1);
      dy = coords_flat(m,2);
      m11 = m11 + dx*dx;
      m12 = m12 + dx*dy;
      m13 = m13 + dx;
      m22 = m22 + dy*dy;
      m23 = m23 + dy;
      z1  = z1 + dx*du;
      z2  = z2 + dy*du;
      z3  = z3 + du;
      z1b = z1b + dx*dv;
      z2b = z2b + dy*dv;
      z3b = z3b + dv;
    end;
    dudx = m22*z1-m23*m23*z1-m12*z2+m13*m23*z2-m13*m22*z3+m12*m23*z3;
    dvdx = m22*z1b-m23*m23*z1b-m12*z2b+m13*m23*z2b-m13*m22*z3b+m12*m23*z3b;
    dudy = -m12*z1+m13*m23*z1+m11*z2-m13*m13*z2+m12*m13*z3-m11*m23*z3;
    dvdy = -m12*z1b+m13*m23*z1b+m11*z2b-m13*m13*z2b+m12*m13*z3b-m11*m23*z3b;
    
    if 0
      denom = -m12*m12+m11*m22-m13*m13*m22+2*m12*m13*m23-m11*m23*m23;
      if (denom~=0)
        defdet = (dudx*dvdy-dvdx*dudy)/(denom*denom);
      else
        defdet = 0;
      end;
    else
      tmp = [dudx dudy; dvdx dvdy];
      defdet = det(tmp);
    end;

    if defdet>0.5 && 0
      
      figure(1); clf; hold on;
      for m=1:num_nbrs
        vn = v_nbrs(m);
        j = find(vn==roi);
        if isempty(j), continue; end;
        x = surf.vertices_wm(v_nbrs(m),1) - surf.vertices_wm(v,1);
        y = surf.vertices_wm(v_nbrs(m),2) - surf.vertices_wm(v,2);
        z = surf.vertices_wm(v_nbrs(m),3) - surf.vertices_wm(v,3);
        line([0,x],[0,y],[0,z]);
        text(x,y,z,num2str(m));
      end;
      title('coords wm');
      
      figure(2); clf; hold on;
      for m=1:num_nbrs
        vn = v_nbrs(m);
        j = find(vn==roi);
        if isempty(j), continue; end;
        x = surf.vertices(v_nbrs(m),1) - surf.vertices(v,1);
        y = surf.vertices(v_nbrs(m),2) - surf.vertices(v,2);
        z = surf.vertices(v_nbrs(m),3) - surf.vertices(v,3);
        line([0,x],[0,y],[0,z]);
        text(x,y,z,num2str(m));
      end;
      title('coords sph');

      figure(3); clf; hold on;
      for m=1:num_nbrs
        vn = v_nbrs(m);
        j = find(vn==roi);
        if isempty(j), continue; end;
        x = coords_flat(m,1);
        y = coords_flat(m,2);
        line([0,x],[0,y]);
        text(x,y,num2str(m));
      end;
      title('coords flat');

      figure(4); clf; hold on;
      for m=1:num_nbrs
        vn = v_nbrs(m);
        j = find(vn==roi);
        if isempty(j), continue; end;
        x = coords_uv(j,1) - coords_uv(k,1);
        y = coords_uv(j,2) - coords_uv(k,2);
        line([0,x],[0,y]);
        text(x,y,num2str(m));
      end;
      title('coords uv');

      keyboard
    end;

    deformation.dudx(k) = dudx;
    deformation.dvdx(k) = dvdx;
    deformation.dudy(k) = dudy;
    deformation.dvdy(k) = dvdy;
    deformation.det(k) = defdet;    
  end;

  deformation.dudxy = sqrt(deformation.dudx.^2 + deformation.dudy.^2);
  deformation.dvdxy = sqrt(deformation.dvdx.^2 + deformation.dvdy.^2);

  toc
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mapvec,parms] = init_ret_grid(parms,datagrid,mapgrid)
  if ~exist('mapgrid','var') | isempty(mapgrid)
    [tmp_mapgrid,parms]=init_map_grid(parms,0);
    [mapgrid,parms]=init_map_grid(parms,1);
    mapgrid.du = mapgrid.u0 - tmp_mapgrid.u0;
    mapgrid.dv = mapgrid.v0 - tmp_mapgrid.v0;
    mapgrid.u0 = tmp_mapgrid.u0;
    mapgrid.v0 = tmp_mapgrid.v0;
  else
    mapgrid.size = size(mapgrid.u);
    mapgrid.num_nodes = prod(mapgrid.size);
  end;

  % define eccentricity and polar angle
  r_vec = linspace(parms.map_r_min,parms.map_r_max,parms.map_gridsize_r);
  if parms.map_logtrans_flag
    phase_vec = (r_vec-parms.map_r_min)/(parms.map_r_max-parms.map_r_min);
    r_vec = parms.map_r_min*(parms.map_r_max/parms.map_r_min).^(phase_vec);
  end;
  th_vec = [];
  for i=1:length(parms.area_labels)
    tmp_th_vec = linspace(parms.area_bound_th(i),parms.area_bound_th(i+1),...
      parms.map_gridsize_th);
    th_vec = [th_vec tmp_th_vec];
  end;
  ind = find(th_vec>pi); % deal with wrap (left hemifield only)
  if ~isempty(ind)
    % left hemifield is from first point to the end
    %   even though vertical meridian of V2/V3 is pi, so won't be in ind
    th_vec(ind(1):end) = th_vec(ind(1):end) - 2*pi;
  end;
  [mapgrid.th,mapgrid.r] = meshgrid(th_vec,r_vec);

  % assign area codes
  mapgrid.area_codes = zeros(mapgrid.size);
  mapgrid.area_labels = parms.area_labels;
  for i=1:mapgrid.size(2)
    area_code = find(strcmp(parms.area_vec{i},parms.area_labels));
    mapgrid.area_codes(:,i) = area_code;
  end;
  mapgrid.pol_r = cos(mapgrid.th);
  mapgrid.pol_i = sin(mapgrid.th);
  tmp_r = (mapgrid.r - parms.data_r_min)*...
    2*pi/(parms.data_r_max - parms.data_r_min);
  mapgrid.ecc_r = cos(tmp_r);
  mapgrid.ecc_i = sin(tmp_r);
  mapgrid.ecc_r(tmp_r<=0) = 0;
  mapgrid.ecc_i(tmp_r<=0) = 0;

  % reshape grid to vector
  mapvec = [];
  mapvec.th = reshape(mapgrid.th,[mapgrid.num_nodes,1]);
  mapvec.r  = reshape(mapgrid.r, [mapgrid.num_nodes,1]);
  mapvec.pol_r = reshape(mapgrid.pol_r,[mapgrid.num_nodes,1]);
  mapvec.pol_i  = reshape(mapgrid.pol_i, [mapgrid.num_nodes,1]);
  mapvec.ecc_r = reshape(mapgrid.ecc_r,[mapgrid.num_nodes,1]);
  mapvec.ecc_i  = reshape(mapgrid.ecc_i, [mapgrid.num_nodes,1]);
  mapvec.u0  = reshape(mapgrid.u0, [mapgrid.num_nodes,1]);
  mapvec.v0  = reshape(mapgrid.v0, [mapgrid.num_nodes,1]);
  mapvec.du  = reshape(mapgrid.du, [mapgrid.num_nodes,1]);
  mapvec.dv  = reshape(mapgrid.dv, [mapgrid.num_nodes,1]);
  mapvec.u  = reshape(mapgrid.u, [mapgrid.num_nodes,1]);
  mapvec.v  = reshape(mapgrid.v, [mapgrid.num_nodes,1]);
  mapvec.area_codes = reshape(mapgrid.area_codes, [mapgrid.num_nodes,1]);
  mapvec.area_labels = parms.area_labels;
  mapvec.grid_size = mapgrid.size;
  mapvec.num_nodes = mapgrid.num_nodes;
  mapvec.num_betas = 2*mapgrid.num_nodes;
  if isfield(mapgrid,'ind_exclude')
    mapvec.ind_exclude = mapgrid.ind_exclude;
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mapgrid,parms] = init_map_grid(parms,length_flag)
  mapgrid = [];
  if parms.map_v123_flag
    switch parms.hemi
      case 'lh'
        parms.area_spacing = [...
          parms.map_v3m_width,...
          parms.map_v2m_width,...
          parms.map_v1m_width,...
          parms.map_v1p_width,...
          parms.map_v2p_width,...
          parms.map_v3p_width,...
        ];
      case 'rh'
        parms.area_spacing = [...
          parms.map_v3p_width,...
          parms.map_v2p_width,...
          parms.map_v1p_width,...
          parms.map_v1m_width,...
          parms.map_v2m_width,...
          parms.map_v3m_width,...
        ];
    end;      
    if length_flag
      switch parms.hemi
        case 'lh'
          parms.area_lengths = [...
            parms.map_v3m_length,...
            parms.map_v2m_length,...
            parms.map_v1m_length,...
            parms.map_v1p_length,...
            parms.map_v2p_length,...
            parms.map_v3p_length,...
          ];
        case 'rh'
          parms.area_lengths = [...
            parms.map_v3p_length,...
            parms.map_v2p_length,...
            parms.map_v1p_length,...
            parms.map_v1m_length,...
            parms.map_v2m_length,...
            parms.map_v3m_length,...
          ];
      end;
    else
      parms.area_lengths = [1 1 1 1 1 1];
    end;
  else
    parms.area_spacing = [1 1];
    parms.area_lengths = [1 1];
  end;

  % create vector of widths for each column of nodes
  parms.area_spacing_vec = [];
  for i=1:length(parms.area_labels)
    for j=1:parms.map_gridsize_th
      parms.area_spacing_vec(end+1) = parms.area_spacing(i);
    end;
  end;
  % normalize spacing so total width = 1
  parms.area_spacing_vec = parms.area_spacing_vec/sum(parms.area_spacing_vec);

  % create vector of lengths for each column of nodes
  parms.area_length_vec = [];
  for i=1:length(parms.area_labels)
    for j=1:parms.map_gridsize_th
      parms.area_length_vec(end+1) = parms.area_lengths(i);
    end;
  end;
  % normalize lengths to max
  parms.area_length_vec = parms.area_length_vec/max(parms.area_length_vec);

  mapgrid.size = [parms.map_gridsize_r,length(parms.area_spacing_vec)];
  mapgrid.num_nodes = prod(mapgrid.size);
  u_vec = [0,cumsum(parms.area_spacing_vec(1:end-1))];
  v_vec = linspace(0,1,mapgrid.size(1));
  [mapgrid.u,mapgrid.v] = meshgrid(u_vec,v_vec);

  % adjust grid spacing for different area lengths
  if any(parms.area_lengths~=1)
    tmp = ones(size(mapgrid.u,1),1)*parms.area_length_vec;
    mapgrid.v = mapgrid.v .* tmp;
  end;

  % apply wedge transform
  if parms.map_model_type>=1 && parms.map_wedge_fact~=1 % wedge
    %   scale u coordinates so that lower edge of map (foveal) is narrower
    %                     while upper edge of map (peripheral) is wider
    min_v = min(mapgrid.v(:));
    max_v = max(mapgrid.v(:));
    b = parms.map_wedge_fact;
    slope = 1-parms.map_wedge_fact;
    x = (mapgrid.v-min_v)/(max_v-min_v);
    y = slope*x + b;
    mapgrid.u = 0.5 + (mapgrid.u - 0.5).*y;
  end;

  % apply polynomial transform
  if parms.map_poly_flag
    u_range = max(mapgrid.u(:)) - min(mapgrid.u(:));
    v_range = max(mapgrid.v(:)) - min(mapgrid.v(:));
    % transform coordinates
    x = mmil_rowvec(mapgrid.u);
    y = mmil_rowvec(mapgrid.v);    
    y = y + polyval(parms.map_poly_coef,x);
    mapgrid.u = reshape(x,size(mapgrid.u));
    mapgrid.v = reshape(y,size(mapgrid.v));
    % shift and rescale to fit on unit grid
    min_u = min(mapgrid.u(:));
    min_v = min(mapgrid.v(:));
    u_range_new = max(mapgrid.u(:)) - min_u;
    v_range_new = max(mapgrid.v(:)) - min_v;
    mapgrid.u = (mapgrid.u - min_u)*u_range/u_range_new;
    mapgrid.v = (mapgrid.v - min_v)*v_range/v_range_new;
  end;

  % apply radial wedge transform
  if parms.map_model_type>=2 && parms.map_radial_wedge_fact~=0
    u_range = max(mapgrid.u(:)) - min(mapgrid.u(:));
    v_range = max(mapgrid.v(:)) - min(mapgrid.v(:));
    % transform coordinates
    r = mapgrid.v + parms.map_radial_offset;
    th = 2*pi*(mapgrid.u - 0.5)*parms.map_radial_wedge_fact;
    mapgrid.u = r.*sin(th);
    mapgrid.v = r.*cos(th);
    % shift and rescale to fit on unit grid
    min_u = min(mapgrid.u(:));
    min_v = min(mapgrid.v(:));
    u_range_new = max(mapgrid.u(:)) - min_u;
    v_range_new = max(mapgrid.v(:)) - min_v;
    mapgrid.u = (mapgrid.u - min_u)*u_range/u_range_new;
    mapgrid.v = (mapgrid.v - min_v)*v_range/v_range_new;
  end;

  % scale map coordinates
  mapgrid.u = (mapgrid.u-0.5)*parms.map_scale_u;
  mapgrid.v = (mapgrid.v-0.5)*parms.map_scale_v;

  % rotate map coordinates
  if parms.map_rotation~=0
    rot_rad = parms.map_rotation*pi/180;
    coords = zeros(mapgrid.num_nodes,3);
    coords(:,1) = reshape(mapgrid.u, [mapgrid.num_nodes,1]);
    coords(:,2) = reshape(mapgrid.v, [mapgrid.num_nodes,1]);
    coords = rc_rot_coords(coords,0,0,rot_rad);
    mapgrid.u = reshape(coords(:,1),mapgrid.size);
    mapgrid.v = reshape(coords(:,2),mapgrid.size);
  end;

  % shift map coordinates
  mapgrid.u = mapgrid.u + 0.5 + parms.map_shift_u;
  mapgrid.v = mapgrid.v + 0.5 + parms.map_shift_v;

  mapgrid.u0 = mapgrid.u;
  mapgrid.v0 = mapgrid.v;
  mapgrid.du = zeros(size(mapgrid.u));
  mapgrid.dv = zeros(size(mapgrid.v));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coords_uv = transform_coords(surf,roi,roi_orig,...
  rotation,shift_u,shift_v,scale_u,scale_v);

  % identify vertices in dilated ROI that are also in the original ROI
  [tmp,ind_roi_orig] = intersect(roi,roi_orig);

  % get coordinates for vertices in roi
  coords = surf.vertices(roi,:);

  % calculate center of mass of roi
  [xc,yc,zc]=rc_roi_centmass(surf,roi);
  [thc,phc,rc]=cart2sph(xc,yc,zc);

  % rotate to get to equator and prime meridian
  coords = rc_rot_coords(coords,0,0,thc);
  coords = rc_rot_coords(coords,0,-phc,0);

  % additional rotation with rotation
  if rotation~=0
    rot_rad = rotation*pi/180;
    coords = rc_rot_coords(coords,rot_rad,0,0);
  end;

  % convert to spherical coordinates, change theta and phi to u and v
  coords_sph = roi_cart2sph(coords);
  coords_uv = coords_sph(:,1:2);

  % scale width of original ROI to 1
  avg_range = mean([range(coords_uv(ind_roi_orig,1)),...
                    range(coords_uv(ind_roi_orig,2))]);
  coords_uv(:,1:2) = coords_uv(:,1:2)/avg_range;

  % additional scaling with scale_u and scale_v
  if scale_u~=1
    coords_uv(:,1) = coords_uv(:,1)*scale_u;
  end;
  if scale_v~=1
    coords_uv(:,2) = coords_uv(:,2)*scale_v;
  end;

  % shift so that u and v range from 0 to 1
  coords_uv = coords_uv + 0.5;

  % additional shifts with shift_u and shift_v
  if shift_u~=0
    coords_uv(:,1) = coords_uv(:,1) + shift_u;
  end;
  if shift_v~=0
    coords_uv(:,2) = coords_uv(:,2) + shift_v;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_coords_sph = roi_cart2sph(roi_coords_cart);
  roi_coords_sph = zeros(size(roi_coords_cart));
  [roi_th,roi_ph,roi_r] = cart2sph(roi_coords_cart(:,1),roi_coords_cart(:,2),roi_coords_cart(:,3));
  roi_coords_sph(:,1)=roi_th;
  roi_coords_sph(:,2)=roi_ph;
  roi_coords_sph(:,3)=roi_r;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_coords_cart = roi_sph2cart(roi_coords_sph);
  roi_coords_cart = zeros(size(roi_coords_sph));
  [roi_x,roi_y,roi_z] = sph2cart(roi_coords_sph(:,1),roi_coords_sph(:,2),roi_coords_sph(:,3));
  roi_coords_cart(:,1)=roi_x;
  roi_coords_cart(:,2)=roi_y;
  roi_coords_cart(:,3)=roi_z;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi = dilate_roi(roi,surf,niters)
  vals = zeros(surf.nverts,1);
  vals(roi) = 1;
  vals = fs_smooth_sparse(surf,vals,niters);
  roi = find(vals);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cost,g] = costfunc(betavec,datagrid,mapvec,parms)
  % cost = fit error + weighted smoothness constraint
  mapvec.du = betavec(1:mapvec.num_nodes);
  mapvec.dv = betavec(mapvec.num_nodes+1:end);
  mapvec.u = mapvec.u0 + mapvec.du;
  mapvec.v = mapvec.v0 + mapvec.dv;
  if nargout>1
    [cost_err,g_err] = costfunc_err(mapvec,datagrid,parms.ecc_fact);
    [cost_sm,cost_fold,g_sm,g_fold] = costfunc_smooth(mapvec);
    [cost_vacancy,g_vacancy] = costfunc_vacancy(mapvec,parms);
    [cost_outbound,g_outbound] = costfunc_outbound(mapvec,datagrid);
    g = parms.err_fact*g_err + parms.smooth_fact*g_sm +...
        parms.fold_fact*g_fold + parms.vacancy_fact*g_vacancy + g_outbound;
  else
    cost_err = costfunc_err(mapvec,datagrid,parms.ecc_fact);
    [cost_sm,cost_fold] = costfunc_smooth(mapvec);
    cost_vacancy = costfunc_vacancy(mapvec,parms);
    cost_outbound = costfunc_outbound(mapvec,datagrid);
  end;
  cost = parms.err_fact*cost_err + parms.smooth_fact*cost_sm +...
         parms.fold_fact*cost_fold + parms.vacancy_fact*cost_vacancy +...
         cost_outbound;
return;

function [cost,g] = costfunc_err(mapvec,datagrid,ecc_fact);
  % fit error: difference of r and th between data and fit
  if ~isfield(mapvec,'ind_exclude'), mapvec.ind_exclude = []; end;

  % look up data th and r based on u and v for each map node
  i_u_float = mapvec.u*datagrid.size(2);
  i_v_float = mapvec.v*datagrid.size(1);
  i_u = 1+floor(mapvec.u*datagrid.size(2));
  i_v = 1+floor(mapvec.v*datagrid.size(1));
  delta_u = (i_u_float - i_u)/datagrid.size(2);
  delta_v = (i_v_float - i_v)/datagrid.size(1);
  % prevent out of bounds access of datagrid
  ind_bad_nodes = find(i_u<=0 | i_u>datagrid.size(2) |...
                       i_v<=0 | i_v>datagrid.size(1));
  i_u(ind_bad_nodes) = 1;
  i_v(ind_bad_nodes) = 1;
  ind = sub2ind(datagrid.size,i_v,i_u);
  pol_r = datagrid.pol_r(ind);
  pol_i = datagrid.pol_i(ind);
  ecc_r = datagrid.ecc_r(ind);
  ecc_i = datagrid.ecc_i(ind);

  if isfield(datagrid,'dpr_du')
    pol_r = pol_r + delta_u.*datagrid.dpr_du(ind) +...
                    delta_v.*datagrid.dpr_dv(ind);
    pol_i = pol_i + delta_u.*datagrid.dpi_du(ind) +...
                    delta_v.*datagrid.dpi_dv(ind);
    ecc_r = ecc_r + delta_u.*datagrid.der_du(ind) +...
                    delta_v.*datagrid.der_dv(ind);
    ecc_i = ecc_i + delta_u.*datagrid.dei_du(ind) +...
                    delta_v.*datagrid.dei_dv(ind);
  end;
  pol_r(ind_bad_nodes) = 0;
  pol_i(ind_bad_nodes) = 0;
  ecc_r(ind_bad_nodes) = 0;
  ecc_i(ind_bad_nodes) = 0;

  % calculate fit error
  pol_r_err = (pol_r - mapvec.pol_r)/datagrid.std_pol_r;
  pol_i_err = (pol_i - mapvec.pol_i)/datagrid.std_pol_i;
  ecc_r_err = (ecc_r - mapvec.ecc_r)/datagrid.std_ecc_r;
  ecc_i_err = (ecc_i - mapvec.ecc_i)/datagrid.std_ecc_i;
  % set error to 0 for regions in err_mask (based on roi_excl)
  ind_exclude = find(datagrid.err_mask(ind)==0);
  ind_exclude = union(ind_exclude,mapvec.ind_exclude);
  pol_r_err(ind_exclude) = 0;
  pol_i_err(ind_exclude) = 0;
  ecc_r_err(ind_exclude) = 0;
  ecc_i_err(ind_exclude) = 0;
  % set error to max for regions outside of data_mask
  ind_nondata = find(datagrid.data_mask(ind)==0);
  pol_r_err(ind_nondata) = 1/datagrid.std_pol_r;
  pol_i_err(ind_nondata) = 1/datagrid.std_pol_i;
  ecc_r_err(ind_nondata) = 1/datagrid.std_ecc_r;
  ecc_i_err(ind_nondata) = 1/datagrid.std_ecc_i;
  % combine pol and ecc error
  cost = ((2-ecc_fact)*(pol_r_err.^2 + pol_i_err.^2) +...
              ecc_fact*(ecc_r_err.^2 + ecc_i_err.^2))/2;
  cost = [cost;cost]; % cost for u and v same at a given node
  nbeta = length(cost);
  cost = cost/nbeta;

  if nargout>1
    % calculate gradient of cost
    dpr_du = datagrid.dpr_du(ind)/datagrid.std_pol_r;
    dpr_dv = datagrid.dpr_dv(ind)/datagrid.std_pol_r;
    dpi_du = datagrid.dpi_du(ind)/datagrid.std_pol_i;
    dpi_dv = datagrid.dpi_dv(ind)/datagrid.std_pol_i;
    der_du = datagrid.der_du(ind)/datagrid.std_ecc_r;
    der_dv = datagrid.der_dv(ind)/datagrid.std_ecc_r;
    dei_du = datagrid.dei_du(ind)/datagrid.std_ecc_i;
    dei_dv = datagrid.dei_dv(ind)/datagrid.std_ecc_i;
    ind_const = union(union(ind_exclude,ind_nondata),ind_bad_nodes);
    dpr_du(ind_const) = 0;
    dpr_dv(ind_const) = 0;
    dpi_du(ind_const) = 0;
    dpi_dv(ind_const) = 0;
    der_du(ind_const) = 0;
    der_dv(ind_const) = 0;
    dei_du(ind_const) = 0;
    dei_dv(ind_const) = 0;
    % multiply by 2 because cost for each node is added twice into sumcost
    g = 2*[((2-ecc_fact)*(pol_r_err.*dpr_du + pol_i_err.*dpi_du) +...
               ecc_fact*(ecc_r_err.*der_du + ecc_i_err.*dei_du));...
            ((2-ecc_fact)*(pol_r_err.*dpr_dv + pol_i_err.*dpi_dv) +...
                ecc_fact*(ecc_r_err.*der_dv + ecc_i_err.*dei_dv))];
    % normalize by number of parameters
    g = g/nbeta;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cost_sm,cost_fold,g_sm,g_fold] = costfunc_smooth(mapvec);
  %% todo: use datagrid.deformation
  
  % smoothness: gradient squared
  [ddu_dr,ddu_dth] = gradient(reshape(mapvec.du,mapvec.grid_size));
  ddu_dr = reshape(ddu_dr,[mapvec.num_nodes,1]);
  ddu_dth = reshape(ddu_dth,[mapvec.num_nodes,1]);
  [ddv_dr,ddv_dth] = gradient(reshape(mapvec.dv,mapvec.grid_size));
  ddv_dr = reshape(ddv_dr,[mapvec.num_nodes,1]);
  ddv_dth = reshape(ddv_dth,[mapvec.num_nodes,1]);
  cost_sm = [ddu_dr.^2 + ddu_dth.^2 ; ddv_dr.^2 + ddv_dth.^2]/2;

  % negative determinant indicates folding
  D = ddu_dr.*ddv_dth - ddu_dth.*ddv_dr;
  cost_fold = -[D;D];

  if nargout>2
    % calculate gradient of cost_sm
    % gradient at each node depends on ddu and ddv at neighboring nodes
    % with chain rule, grad(cost_sm) =
    %      (-neighbor weights from gradient) * (ddudr + ...)
    [g_ur_r,g_ur_th] = gradient(reshape(ddu_dr,mapvec.grid_size));
    [g_uth_r,g_uth_th] = gradient(reshape(ddu_dth,mapvec.grid_size));
    [g_vr_r,g_vr_th] = gradient(reshape(ddv_dr,mapvec.grid_size));
    [g_vth_r,g_vth_th] = gradient(reshape(ddv_dth,mapvec.grid_size));
    g_u = g_ur_r + g_uth_th;
    g_v = g_vr_r + g_vth_th;
    g_u = reshape(g_u,[mapvec.num_nodes,1]);
    g_v = reshape(g_v,[mapvec.num_nodes,1]);
    g_sm = -[g_u;g_v];

    % calculate gradient of cost_fold
    g_u = g_uth_r - g_vth_r + g_vr_th - g_ur_th;    
    g_v = g_vth_r - g_uth_r + g_ur_th - g_vr_th;    
    g_u = reshape(g_u,[mapvec.num_nodes,1]);
    g_v = reshape(g_v,[mapvec.num_nodes,1]);
    g_fold = [g_u;g_v];
    g_fold(cost_fold<0) = 0; 
  end;

  % only penalize nodes with negative determinant
  cost_fold(cost_fold<0) = 0; 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cost,g] = costfunc_vacancy(mapvec,parms);
  cost = zeros(mapvec.num_betas,1); g = zeros(mapvec.num_betas,1);
  if parms.vacancy_fact==0, return; end;
  % ratio of grid points with no th value to all grid points
  mapgrid_sub = resamp_mapgrid_sub(parms,mapvec);
  ind_map = find(mapgrid_sub.area_mask);
  ind_mask = find(parms.sub_datagrid.roi_mask);
  ind_inter = intersect(ind_map,ind_mask);
  n_inter = length(ind_inter);
  n_mask = length(ind_mask);
  cost = (n_mask - n_inter) / n_mask;
  cost = cost * ones(mapvec.num_betas,1) / mapvec.num_betas;
  %% todo: is it possible to calculate g?
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cost,g] = costfunc_outbound(mapvec,datagrid);
  cost = zeros(mapvec.num_betas,1); g = zeros(mapvec.num_betas,1);
  % penalty for map node outside of ROI
  if isempty(datagrid.outbound_penalty), return; end;
  
  % look up datagrid.outbound_penalty based on u and v for each map node
  i_u_float = mapvec.u*datagrid.size(2);
  i_v_float = mapvec.v*datagrid.size(1);
  i_u = 1+floor(mapvec.u*datagrid.size(2));
  i_v = 1+floor(mapvec.v*datagrid.size(1));
  delta_u = (i_u_float - i_u)/datagrid.size(2);
  delta_v = (i_v_float - i_v)/datagrid.size(1);
  % prevent out of bounds access of datagrid
  ind_bad_nodes = find(i_u<=0 | i_u>datagrid.size(2) |...
                       i_v<=0 | i_v>datagrid.size(1));
  i_u(ind_bad_nodes) = 1;
  i_v(ind_bad_nodes) = 1;
  ind = sub2ind(datagrid.size,i_v,i_u);
  penalty = datagrid.outbound_penalty(ind);
  penalty(ind_bad_nodes) = max(datagrid.outbound_penalty(:));

  cost = (penalty.^2)/2;
  cost = [cost;cost]; % cost for u and v same at a given node
  nbeta = length(cost);
  cost = cost/nbeta;

  if nargout>1
    % calculate gradient of cost
    dpenalty_du = datagrid.dpenalty_du(ind);
    dpenalty_dv = datagrid.dpenalty_dv(ind);
    dpenalty_du(ind_bad_nodes) = 0;
    dpenalty_dv(ind_bad_nodes) = 0;
    % multiply by 2 because cost for each node is added twice into sumcost
    g = 2*[(penalty.*dpenalty_du);(penalty.*dpenalty_dv)];
    % normalize by number of parameters
    g = g/nbeta;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out_betas=rand_step(in_betas,step_size,max_step_size)
  dx = step_size*randn(size(in_betas));
  dx(dx<-max_step_size) = -max_step_size;
  dx(dx>max_step_size) = max_step_size;
  out_betas = in_betas + dx;
return;

function out_betas=grad_step(in_betas,g,max_step_size)
  dx = -max_step_size*g/max(abs(g));
  dx(~isfinite(dx)) = 0;
  out_betas = in_betas + dx;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g=cost_grad_b(beta,b,datagrid,mapvec,parms);
  delta = zeros(size(beta));
  delta(b) = parms.step_size;
  cost_plus = sum(costfunc(beta+delta,datagrid,mapvec,parms));
  cost_minus = sum(costfunc(beta-delta,datagrid,mapvec,parms));
  g = (cost_plus-cost_minus)/(2*parms.step_size);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mapgrid_sub = resamp_mapgrid_sub(parms,mapvec)
%  fprintf('%s: resampling map to sub grid...\n',mfilename);
  % resample deformed map to sub grid
  mapgrid_sub = [];
  X = mapvec.u;
  Y = mapvec.v;
  gridsize_u = size(parms.sub_datagrid.roi_mask,2);
  gridsize_v = size(parms.sub_datagrid.roi_mask,1);

  % area labels
  mapgrid_sub.area_labels = parms.area_labels;
  mapgrid_sub.area_codes = griddata(X,Y,mapvec.area_codes,...
    parms.sub_datagrid.X,parms.sub_datagrid.Y,parms.interp_method);
  % get rid of fractional area codes caused by interpolation
  tmp = round(mapgrid_sub.area_codes);
  d = abs(tmp - mapgrid_sub.area_codes);
  mapgrid_sub.area_codes(d>parms.area_code_smf) = 0;
  mapgrid_sub.area_codes(d<=parms.area_code_smf) =...
                      tmp(d<=parms.area_code_smf);

  % create mask around grid points with defined area codes
  area_mask = zeros(parms.sub_datagrid.size);
  area_mask(mapgrid_sub.area_codes>0) = 1;

  % dilate mask (fill gaps)
  sig_u = parms.area_smooth_sigma*gridsize_u;
  sig_v = parms.area_smooth_sigma*gridsize_v;
  area_mask = mmil_smooth2d(area_mask,sig_u,sig_v);
  area_mask(area_mask>=parms.area_smooth_thresh) = 1;
  area_mask(area_mask<parms.area_smooth_thresh) = 0;
  area_mask(isnan(area_mask)) = 0;

  % erode mask
  sig_u = parms.area_smooth_sigma2*gridsize_u;
  sig_v = parms.area_smooth_sigma2*gridsize_v;
  area_mask = mmil_smooth2d(area_mask,sig_u,sig_v);
  area_mask(area_mask>=parms.area_smooth_thresh2) = 1;
  area_mask(area_mask<parms.area_smooth_thresh2) = 0;
  area_mask(isnan(area_mask)) = 0;

  mapgrid_sub.area_mask = area_mask;  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mapgrid_full = resamp_mapgrid_full(parms,mapvec,datagrid)
  fprintf('%s: resampling map to full grid...\n',mfilename);
  % resample deformed map to full grid
  X = mapvec.u;
  Y = mapvec.v;
  mapgrid_full = [];

  % polar angle
  mapgrid_full.pol_r = griddata(X,Y,mapvec.pol_r,...
    datagrid.X,datagrid.Y,parms.interp_method);
  mapgrid_full.pol_i = griddata(X,Y,mapvec.pol_i,...
    datagrid.X,datagrid.Y,parms.interp_method);
  ind = find(mapgrid_full.pol_r==0 & mapgrid_full.pol_i==0);
  mapgrid_full.pol_r(ind) = NaN;
  mapgrid_full.pol_i(ind) = NaN;
  mapgrid_full.th = atan2(mapgrid_full.pol_i,mapgrid_full.pol_r);

  % eccentricity
  mapgrid_full.r  = griddata(X,Y,mapvec.r,...
    datagrid.X,datagrid.Y,parms.interp_method);

  % error
  mapgrid_full.cost_err = griddata(X,Y,mapvec.cost_err,...
    datagrid.X,datagrid.Y,parms.interp_method);
  mapgrid_full.cost_smooth_u = griddata(X,Y,mapvec.cost_smooth_u,...
    datagrid.X,datagrid.Y,parms.interp_method);
  mapgrid_full.cost_smooth_v = griddata(X,Y,mapvec.cost_smooth_v,...
    datagrid.X,datagrid.Y,parms.interp_method);
  mapgrid_full.cost_fold = griddata(X,Y,mapvec.cost_fold,...
    datagrid.X,datagrid.Y,parms.interp_method);
  mapgrid_full.cost_outbound = griddata(X,Y,mapvec.cost_outbound,...
    datagrid.X,datagrid.Y,parms.interp_method);

  % area labels
  mapgrid_full.area_labels = parms.area_labels;
  mapgrid_full.area_codes = griddata(X,Y,mapvec.area_codes,...
    datagrid.X,datagrid.Y,parms.interp_method);
  % get rid of fractional area codes caused by interpolation
  tmp = round(mapgrid_full.area_codes);
  d = abs(tmp - mapgrid_full.area_codes);
  mapgrid_full.area_codes(d>parms.area_code_smf) = 0;
  mapgrid_full.area_codes(d<=parms.area_code_smf) =...
                      tmp(d<=parms.area_code_smf);

  % create mask around grid points with defined area codes
  area_mask = zeros(datagrid.size);
  area_mask(mapgrid_full.area_codes>0) = 1;

  % dilate mask (fill gaps)
  sig_u = parms.area_smooth_sigma*parms.data_gridsize_u;
  sig_v = parms.area_smooth_sigma*parms.data_gridsize_v;
  area_mask = mmil_smooth2d(area_mask,sig_u,sig_v);
  area_mask(area_mask>=parms.area_smooth_thresh) = 1;
  area_mask(area_mask<parms.area_smooth_thresh) = 0;
  area_mask(isnan(area_mask)) = 0;

  % erode mask
  sig_u = parms.area_smooth_sigma2*parms.data_gridsize_u;
  sig_v = parms.area_smooth_sigma2*parms.data_gridsize_v;
  area_mask = mmil_smooth2d(area_mask,sig_u,sig_v);
  area_mask(area_mask>=parms.area_smooth_thresh2) = 1;
  area_mask(area_mask<parms.area_smooth_thresh2) = 0;
  area_mask(isnan(area_mask)) = 0;

  mapgrid_full.area_codes(find(area_mask)) = tmp(find(area_mask));
  
  mapgrid_full.pol_r(~area_mask) = NaN;
  mapgrid_full.pol_i(~area_mask) = NaN;
  mapgrid_full.th(~area_mask) = NaN;
  mapgrid_full.r(~area_mask) = NaN;
  mapgrid_full.cost_err(~area_mask) = 0;
  mapgrid_full.cost_smooth_u(~area_mask) = 0;
  mapgrid_full.cost_smooth_v(~area_mask) = 0;
  mapgrid_full.cost_fold(~area_mask) = 0;
  mapgrid_full.cost_outbound(~area_mask) = 0;

  % save info needed for mapping back to surface
  mapgrid_full.origX = datagrid.origX;
  mapgrid_full.origY = datagrid.origY;
  mapgrid_full.size = datagrid.size;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fit_data,area_masks,data] = resamp_results_surf(parms,...
                                                    mapgrid_full,roi,surf,data)
  fprintf('%s: resampling deformed map to cortical surface vertices...\n',mfilename);
  u = mapgrid_full.origX;
  v = mapgrid_full.origY;
  i_u = 1+floor(u*mapgrid_full.size(2));
  i_v = 1+floor(v*mapgrid_full.size(1));
  % prevent out of bounds access of datagrid
  ind_bad_nodes = find(i_u<=0 | i_u>mapgrid_full.size(2) |...
                       i_v<=0 | i_v>mapgrid_full.size(1));
  i_u(ind_bad_nodes) = 1;
  i_v(ind_bad_nodes) = 1;
  ind = sub2ind(mapgrid_full.size,i_v,i_u);
  pol_r = mapgrid_full.pol_r(ind);
  pol_i = mapgrid_full.pol_i(ind);
  th = mapgrid_full.th(ind);
  r = mapgrid_full.r(ind);
  cost_err = mapgrid_full.cost_err(ind);
  cost_smooth_u = mapgrid_full.cost_smooth_u(ind);
  cost_smooth_v = mapgrid_full.cost_smooth_v(ind);
  cost_fold = mapgrid_full.cost_fold(ind);
  area_codes = mapgrid_full.area_codes(ind);

  th(ind_bad_nodes) = NaN;
  ind_bad_nodes = find(isnan(th) | isnan(r)); % outside of grid
  pol_r(ind_bad_nodes) = 0;
  pol_i(ind_bad_nodes) = 0;
  th(ind_bad_nodes) = 0;
  r(ind_bad_nodes) = 0;
  cost_u(ind_bad_nodes) = 0;
  cost_v(ind_bad_nodes) = 0;
  cost_err(ind_bad_nodes) = 0;
  cost_smooth_u(ind_bad_nodes) = 0;
  cost_smooth_v(ind_bad_nodes) = 0;
  cost_fold(ind_bad_nodes) = 0;
  area_codes(ind_bad_nodes) = 0;

  area_masks = struct('name',parms.area_labels,'vertices',[],'vertices_smoothed',[]);
  area_codes = round(area_codes);
  area_verts = [];
  for a=1:length(area_masks)
    vertlist = [];
    for i=1:length(mapgrid_full.area_labels) % includes borders
      if ~isempty(findstr(area_masks(a).name,mapgrid_full.area_labels{i}))
        verts = find(area_codes==i);
        vertlist = [vertlist;verts];
      end;  
    end;
    area_masks(a).vertices = sort(roi(vertlist));
    area_verts = [area_verts;area_masks(a).vertices];
    % smooth area_masks
    if parms.area_smooth_steps
      tmp_vals = zeros(surf.nverts,1);
      tmp_vals(area_masks(a).vertices) = 1;
      tmp_vals = fs_smooth_sparse(surf,tmp_vals,parms.area_smooth_steps);
      area_masks(a).vertices_smoothed = find(tmp_vals);
    end;
  end;
  % find vertices with no area label
  [bad_verts,ind] = setdiff(roi,area_verts);

  % create fit_data output
  fit_data.u = sparse(surf.nverts,1);
  fit_data.v = sparse(surf.nverts,1);
  fit_data.pol_r = sparse(surf.nverts,1);
  fit_data.pol_i = sparse(surf.nverts,1);
  fit_data.th = sparse(surf.nverts,1);
  fit_data.r = sparse(surf.nverts,1);
  fit_data.cost_err = sparse(surf.nverts,1);
  fit_data.cost_smooth_u = sparse(surf.nverts,1);
  fit_data.cost_smooth_v = sparse(surf.nverts,1);
  fit_data.cost_fold = sparse(surf.nverts,1);
  fit_data.u(roi) = u;
  fit_data.v(roi) = v;
  fit_data.pol_r(roi) = pol_r;
  fit_data.pol_i(roi) = pol_i;
  fit_data.th(roi) = th;
  fit_data.r(roi) = r;
  fit_data.pol_r(bad_verts) = 0; % replace values at bad verts with 0s
  fit_data.pol_i(bad_verts) = 0;
  fit_data.th(bad_verts) = 0;
  fit_data.r(bad_verts) = 0;
  tmp_r = (fit_data.r - parms.data_r_min)*...
    2*pi/(parms.data_r_max - parms.data_r_min);
  tmp_r(union(bad_verts,find(fit_data.r==0))) = 0;    
  tmp_r = sparse(tmp_r);
  fit_data.ecc_r = cos(tmp_r);
  fit_data.ecc_i = sin(tmp_r);
  fit_data.ecc_r(tmp_r<=0) = 0;
  fit_data.ecc_i(tmp_r<=0) = 0;
  fit_data.cost_err(roi) = cost_err;
  fit_data.cost_smooth_u(roi) = cost_smooth_u;
  fit_data.cost_smooth_v(roi) = cost_smooth_v;
  fit_data.cost_fold(roi) = cost_fold;

  if parms.smooth_sparse_niters
    % sparsely smooth fit_data
    fit_data.pol_r = fs_smooth_sparse(surf,fit_data.pol_r,parms.smooth_sparse_niters);
    fit_data.pol_i = fs_smooth_sparse(surf,fit_data.pol_i,parms.smooth_sparse_niters);
    fit_data.ecc_r = fs_smooth_sparse(surf,fit_data.ecc_r,parms.smooth_sparse_niters);
    fit_data.ecc_i = fs_smooth_sparse(surf,fit_data.ecc_i,parms.smooth_sparse_niters);

    % recalculate th
    fit_data.th = atan2(fit_data.pol_i,fit_data.pol_r);
    fit_data.th(~isfinite(fit_data.th)) = 0;

    % recalculate r
    ecc_phase = atan2(fit_data.ecc_i,fit_data.ecc_r)/(2*pi);
    ecc_phase(ecc_phase<0) = ecc_phase(ecc_phase<0)+1;
    ind_ecc = find(ecc_phase);
    tmp_r = parms.data_r_min +...
      ecc_phase*(parms.data_r_max - parms.data_r_min);
    fit_data.r = sparse(surf.nverts,1);
    fit_data.r(ind_ecc) = tmp_r(ind_ecc);
    fit_data.r(~isfinite(fit_data.r)) = 0;
  end;

  % adjust ecc_r and ecc_i so phase reflects logtrans_flag
  ind_ecc = find(fit_data.r);
  if parms.data_logtrans_flag
    ecc_phase = ...
      log(fit_data.r/(parms.data_r_min*log(parms.data_r_max/parms.data_r_min)));
  else
    ecc_phase = ...
      (fit_data.r - parms.data_r_min)/(parms.data_r_max-parms.data_r_min);
  end;
  tmp = sparse(surf.nverts,1);
  tmp(ind_ecc) = 2*pi*ecc_phase(ind_ecc);
  ecc_phase = tmp;
  fit_data.ecc_r = cos(ecc_phase);
  fit_data.ecc_i = sin(ecc_phase);
  fit_data.ecc_r(fit_data.r==0) = 0;
  fit_data.ecc_i(fit_data.r==0) = 0;

  % apply threshold
  ind = find(fit_data.ecc_r~=0 &...
             fit_data.ecc_i~=0 &...
             abs(fit_data.ecc_r)<parms.fit_data_smf &...
             abs(fit_data.ecc_i)<parms.fit_data_smf);
  fit_data.ecc_r(ind) = 0;
  fit_data.ecc_i(ind) = 0;

  ind = find(fit_data.pol_r~=0 &...
             fit_data.pol_i~=0 &...
             abs(fit_data.pol_r)<parms.fit_data_smf &...
             abs(fit_data.pol_i)<parms.fit_data_smf);
  fit_data.pol_r(ind) = 0;
  fit_data.pol_i(ind) = 0;

  ind = find(fit_data.r | fit_data.th);

  % calculate th for data
  tmp_th = atan2(data.pol_i,data.pol_r);
  data.th = sparse(surf.nverts,1);
  data.th(ind) = tmp_th(ind);

  % calculate r for data
  ecc_phase = atan2(data.ecc_i,data.ecc_r)/(2*pi);
  ecc_phase(ecc_phase<0) = ecc_phase(ecc_phase<0)+1;
  % take into account log transformation
  if parms.data_logtrans_flag
    tmp_r = parms.data_r_min *...
      exp(ecc_phase*log(parms.data_r_max/parms.data_r_min));
  else
    tmp_r = parms.data_r_min +...
      (parms.data_r_max - parms.data_r_min)*ecc_phase;
  end;
  data.r = sparse(surf.nverts,1);
  data.r(ind) = tmp_r(ind);
  data.r(~isfinite(data.r)) = 0;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_grad(parms,datagrid,mapvec)
  beta = [mapvec.du;mapvec.dv];
  if parms.err_fact==0 && parms.smooth_fact~=0
%    beta = rand_step(beta,parms.max_grid_step_size,2*parms.max_grid_step_size);
    x = 1:mapvec.grid_size(2);
    y = 1:mapvec.grid_size(1);
    [X,Y] = meshgrid(x,y);
    Zu = cos(X/5) + Y.^2/100;
    Zv = sin(X/4) - Y/20;
    if 0
      clf;
      subplot(2,1,1);
      mesh(X,Y,Zu);
      subplot(2,1,2);
      mesh(X,Y,Zv);
    end;
    step_u = parms.max_grid_step_size*reshape(Zu,[mapvec.num_nodes,1]);
    step_v = parms.max_grid_step_size*reshape(Zv,[mapvec.num_nodes,1]);
    beta = beta + [step_u;step_v];
  end;
  nbeta = length(beta);

  [cost_init,g_init] = costfunc(beta,datagrid,mapvec,parms);

  fprintf('%s: checking gradients...\n',mfilename);
  g = zeros(nbeta,1);
  for b=1:nbeta
    g(b) = cost_grad_b(beta,b,datagrid,mapvec,parms);
  end;

  med_ratio = median(g./(g_init+eps));
  fprintf('%s: median ratio between g and g_init = %0.3f\n',...
    mfilename,med_ratio);
  if isnan(med_ratio) || med_ratio>1.1  || med_ratio<0.9
    fprintf('%s: WARNING: calculated gradients do not match prediction\n',...
      mfilename);
    keyboard
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_prereg_parms(parms)
  for t=1:length(parms.prereg_tags)
    tag = parms.prereg_tags{t};
    val = parms.(tag);
    nval = length(val);
    range_tag = sprintf('%s_range',tag);
    range_val = parms.(range_tag);
    val = check_bounds(val,range_val);
    parms.(tag) = val;
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = randstep_prereg_parms(parms)
  beta = prereg_parms2beta(parms);
  [lbounds,ubounds] = prereg_parms2bounds(parms);
  for i=1:length(beta)
    tmp_range = [lbounds(i),ubounds(i)];
    beta(i) = check_bounds(beta(i) +...
      parms.prereg_step_size*range(tmp_range)*randn,...
      tmp_range);
  end;
  parms = prereg_beta2parms(beta,parms);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val=check_bounds(val,bounds);
  val = max(val,bounds(1));
  val = min(val,bounds(2));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set beta vector from parameters struct
function beta=prereg_parms2beta(parms)
  ind_vars = find(parms.prereg_vars);
  nvars = length(ind_vars);
  beta = [];
  for n=1:nvars
    t = ind_vars(n);
    tag = parms.prereg_tags{t};
    val = parms.(tag);
    nval = length(val);
    beta(end+1:end+nval) = val;
  end;
return;

% set fields in parameters struct from beta vector
function parms=prereg_beta2parms(beta,parms)
  ind_vars = find(parms.prereg_vars);
  nvars = length(ind_vars);
  j = 1;
  for n=1:nvars
    t = ind_vars(n);
    tag = parms.prereg_tags{t};
    nval = length(parms.(tag));
    val = beta(j:j+nval-1);
    parms.(tag) = val;
    j = j + nval;
  end;
return;

% extract prereg parameters from parms struct
function prereg_parms = parms2prereg_parms(parms)
  prereg_parms = [];
  for t=1:length(parms.prereg_tags)
    prereg_parms.(parms.prereg_tags{t}) = parms.(parms.prereg_tags{t});
  end;
return;

% add prereg parameters to parms struct
function parms = prereg_parms2parms(prereg_parms,parms)
  for t=1:length(parms.prereg_tags)
    parms.(parms.prereg_tags{t}) = prereg_parms.(parms.prereg_tags{t});
  end;
return;

% set bounds vectors from parameter ranges
function [lbounds,ubounds]=prereg_parms2bounds(parms)
  ind_vars = find(parms.prereg_vars);
  nvars = length(ind_vars);
  lbounds = [];
  ubounds = [];
  for n=1:nvars
    t = ind_vars(n);
    tag = parms.prereg_tags{t};
    val = parms.(tag);
    nval = length(val);
    tag = sprintf('%s_range',tag);
    val = parms.(tag);
    lbounds(end+1:end+nval) = val(1)*ones(1,nval);
    ubounds(end+1:end+nval) = val(2)*ones(1,nval);
  end;
return;

%% TODO: change prereg_costfunc to return vector of differences
%   if prereg_search_type is lsqnonlin
function cost=prereg_costfunc(beta,datagrid,parms,msg_flag);
  if ~exist('msg_flag','var') || isempty(msg_flag), msg_flag = 0; end;
  parms = prereg_beta2parms(beta,parms);
  mapvec = init_ret_grid(parms,datagrid);
  cost_err = costfunc_err(mapvec,datagrid,parms.ecc_fact);
  cost_vacancy = costfunc_vacancy(mapvec,parms);
  cost_outbound = costfunc_outbound(mapvec,datagrid);
  cost = sum(parms.err_fact*cost_err + ...
             parms.vacancy_fact*cost_vacancy + cost_outbound);
  if msg_flag
    fprintf(' total cost = %0.4f; err = %0.4f, vacancy = %0.4f, outbound = %0.4f\n',...
      cost,parms.err_fact*sum(cost_err),...
      parms.vacancy_fact*sum(cost_vacancy),sum(cost_outbound));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_data(datagrid,parms)
  fprintf('%s: plotting data on grid...\n',mfilename);

  figure(1); clf;
  subplot(3,1,1);
  crange = [-pi,pi];
  imagesc(datagrid.th,crange);
  colorbar; axis xy;
  title('pol phase');
  subplot(3,1,2);
  crange = [0,2*pi];
  imagesc(2*pi*datagrid.ecc_phase,crange);
  colorbar; axis xy;
  title('ecc phase');
  subplot(3,1,3);
  crange = [0,parms.data_r_max];
  imagesc(datagrid.r,crange);
  colorbar; axis xy;
  title('ecc vang');

  if parms.max_outbound_penalty>0
    figure(2); clf;
    subplot(3,1,1);
    imagesc(datagrid.outbound_penalty);
    colorbar; axis xy;
    title('outbound penalty');
  end;
  drawnow;
return

function plot_gradients(datagrid,parms)
  fprintf('%s: plotting data gradients...\n',mfilename);
  figure(3); clf;
  subplot(2,2,1);
  imagesc(datagrid.dpr_du);
  colorbar; axis xy;
  title('dpr/du');
  subplot(2,2,2);
  imagesc(datagrid.dpi_du);
  colorbar; axis xy;
  title('dpi/du');
  subplot(2,2,3);
  imagesc(datagrid.dpr_dv);
  colorbar; axis xy;
  title('dpr/dv');
  subplot(2,2,4);
  imagesc(datagrid.dpi_dv);
  colorbar; axis xy;
  title('dpi/dv');

  figure(4); clf;
  subplot(2,2,1);
  imagesc(datagrid.der_du);
  colorbar; axis xy;
  title('der/du');
  subplot(2,2,2);
  imagesc(datagrid.dei_du);
  colorbar; axis xy;
  title('dei/du');
  subplot(2,2,3);
  imagesc(datagrid.der_dv);
  colorbar; axis xy;
  title('der/dv');
  subplot(2,2,4);
  imagesc(datagrid.dei_dv);
  colorbar; axis xy;
  title('dei/dv');

  if parms.max_outbound_penalty>0
    figure(5); clf;
    subplot(1,2,1);
    imagesc(datagrid.dpenalty_du);
    colorbar; axis xy;
    title('dpenalty/du');
    subplot(1,2,2);
    imagesc(datagrid.dpenalty_dv);
    colorbar; axis xy;
    title('dpenalty/dv');
  end;
  drawnow;
return

function plot_datafit(datagrid,mapgrid_full,parms)
  figure(11); clf;
  crange = [-pi,pi];
  tmp = datagrid.th;
  tmp(datagrid.roi_mask==0) = NaN;
  imagesc(tmp,crange);
  colorbar; axis xy;
  title('data pol phase');

  figure(12); clf;
  imagesc(mapgrid_full.th,crange);
  colorbar; axis xy;
  title('map th');

  figure(21); clf;
  crange = [0,parms.data_r_max];
  tmp = datagrid.r;
  tmp(datagrid.roi_mask==0) = NaN;
  imagesc(tmp,crange);
  colorbar; axis xy;
  title('data ecc vang');

  figure(22); clf;
  crange = [0,parms.data_r_max];
  imagesc(mapgrid_full.r,crange);
  colorbar; axis xy;
  title('map r');
  drawnow;
return

function plot_datafit_sub(datagrid,mapgrid_full,parms)
  figure(31); clf;
  subplot(3,1,1)
  crange = [-pi,pi];
  imagesc(datagrid.th,crange);
  colorbar; axis xy;
  title('data pol phase');
  subplot(3,1,2)
  imagesc(mapgrid_full.th,crange);
  colorbar; axis xy;
  title('map th');
  subplot(3,1,3)
  tmp = datagrid.th;
  tmp(~isnan(mapgrid_full.th)) = -10;
  imagesc(tmp,crange);
  colorbar; axis xy;
  title('overlay');

  figure(32); clf;
  subplot(1,3,1)
  crange = [0,parms.data_r_max];
  imagesc(datagrid.r,crange);
  colorbar; axis xy;
  title('data ecc vang');
  subplot(1,3,2)
  crange = [0,parms.data_r_max];
  imagesc(mapgrid_full.r,crange);
  colorbar; axis xy;
  title('map r');
  subplot(1,3,3)
  tmp = datagrid.r;
  tmp(~isnan(mapgrid_full.r)) = -10;
  imagesc(tmp,crange);
  colorbar; axis xy;
  title('overlay');
  drawnow;
return

function plot_cost(datagrid,mapgrid_full,parms)
  figure(51); clf;
  subplot(2,3,1)
  imagesc(datagrid.roi_mask .* datagrid.err_mask);
  colorbar; axis xy;
  title('roi mask');
  if parms.err_fact>0
    subplot(2,3,2)
    imagesc(mapgrid_full.cost_err);
    colorbar; axis xy;
    title('cost err');
  end;
  if parms.max_outbound_penalty>0
    subplot(2,3,3)
    imagesc(mapgrid_full.cost_outbound);
    colorbar; axis xy;
    title('cost outbound');
  end;
  if parms.smooth_fact>0
    subplot(2,3,4)
    imagesc(mapgrid_full.cost_smooth_u);
    colorbar; axis xy;
    title('cost smooth u');
    subplot(2,3,5)
    imagesc(mapgrid_full.cost_smooth_v);
    colorbar; axis xy;
    title('cost smooth v');
  end;
  if parms.fold_fact>0
    subplot(2,3,6)
    imagesc(mapgrid_full.cost_fold);
    colorbar; axis xy;
    title('cost fold');
  end;
  drawnow;
return

function plot_areas(mapgrid_full)
  figure(41); clf;
  area_codes = mapgrid_full.area_codes;
  area_codes(isnan(area_codes)) = 0;
  imagesc(area_codes);
  colorbar; axis xy;
  title('area codes');
  drawnow;
return

