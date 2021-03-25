function [vf2ctx,retfit_results] = rc_calc_vf2ctx(retfit_results,varargin)
%function [vf2ctx,retfit_results] = rc_calc_vf2ctx(retfit_results,[options])
%
% Purpose: calculate visual field to cortex mapping
%
% Required Input:
%   retfit_results: struct containing results from retinotopy fitting
%     (see rc_retfit, rc_load_retfit_results)
%
% Optional Input:
%  'w_thresh': threshold applied to weights (relative to max)
%    {default = 0.01}
%  'vfnorm_flag': [0|1] for applying w_thresh for vf2ctx
%     0: normalize to global max
%     1: normalize to max for each cortical location
%     {default = 1}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'rf_sizes': vector of receptive field sizes for each visual area
%    Will be used as initial estimate if rf_niters>0
%    {default = [1,1,1]}
%  'rf_slopes': vector of slopes of linear trend of receptive field sizes
%    w.r.t. ecc for each visual area
%    Will be used as initial estimate if rf_niters>0
%    {default = [0.1,0.1,0.1]}
%  'rf_r_inter': radius value used as intercept at which rf size is rf_sizes
%    {default = 6}
%  'surround_flag': [0|1] model center and surround with difference of Gaussians
%    {default = 0}
%  'surround_rf_fact': size of surround rf relative to center rf
%    {default = 2}
%  'surround_amp_fact': amplitude of surround relative to center
%    {default = 0.2}
%  'retfit_data_flag': [0|1] whether to use refit area ROIs but original
%    retinotopy data instead of template values for selection of dipole clusters
%    {default = 0}
%
% Output:
%   vf2ctx: struct array with size = [nareas,nhemis] with fields:
%     verts: vector of vertices for each area
%     map: matrix of size [nverts,nstimres*nstimres]
%     
% Created:  04/18/11 by Don Hagler
% Last Mod: 12/08/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vf2ctx = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(retfit_results,varargin);

% create visual field grid
[parms.stim_x,parms.stim_y] = define_vf_grid(parms);

% define visual field to cortex mappings
vf2ctx = define_vf2ctx_maps(retfit_results,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(retfit_results,options)
  parms = mmil_args2parms(options, { ...
    'w_thresh',0.01,[0,100],...
    'vfnorm_flag',true,[false true],...
    'r_max',12.5,[0,Inf],...
    'rf_sizes',[1,1,1],[0.01,10],...
    'rf_slopes',[0.1,0.1,0.1],[0,10],...
    'rf_r_inter',6,[0,Inf],...
    'surround_flag',false,[false true],...
    'surround_rf_fact',2,[0,100],...
    'surround_amp_fact',0.2,[0,100],...
    'retfit_data_flag',false,[false true],...
    'area_names',{'v1','v2','v3'},[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'stimres',100,[50,1000],...
  });
  parms.nareas = length(parms.area_names);
  parms.nhemis = length(parms.hemilist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stim_x,stim_y] = define_vf_grid(parms)
  % initialize grid
  [X,Y] = meshgrid(1:parms.stimres);
  stim_x = reshape(X,[parms.stimres.^2,1]);
  stim_y = reshape(Y,[parms.stimres.^2,1]);
  % subtract center, scale
  c = mean(1:parms.stimres);
  maxR = parms.stimres/2;
  stim_x = (stim_x - c)*parms.r_max/maxR;
  stim_y = (c - stim_y)*parms.r_max/maxR;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vf2ctx = define_vf2ctx_maps(retfit_results,parms)
  vf2ctx = [];
  for a=1:parms.nareas
    for h=1:parms.nhemis
      [vf2ctx(a,h).verts,vf2ctx(a,h).map] = ...
        define_vf2ctx_map(retfit_results,parms,a,h);
      vf2ctx(a,h).area_name = parms.area_names{a};
      vf2ctx(a,h).hemi = parms.hemilist{h};
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [verts,map] = define_vf2ctx_map(retfit_results,parms,a,h)
  verts = []; map = [];
  area_name = parms.area_names{a};
  hemi = parms.hemilist{h};

  % get area masks
  area_masks = retfit_results.([hemi '_area_masks']);

  % get data
  if parms.retfit_data_flag
    fit_data = retfit_results.([hemi '_data']);
  else
    fit_data = retfit_results.([hemi '_area_data']);
    if isempty(fit_data)
      fit_data = retfit_results.([hemi '_fit_data']);
    end;
  end;

  % identify vertices for this area
  verts = []; th = []; r = [];
  for i=1:length(area_masks)
    mask_name = area_masks(i).name;
    if isempty(regexp(mask_name,[area_name,'[+-]'])), continue; end;
    tmp_verts = area_masks(i).vertices; 
    if length(fit_data)>1
      tmp_th = full(fit_data(i).th(tmp_verts));
      tmp_r = full(fit_data(i).r(tmp_verts));
    else
      tmp_th = full(fit_data.th(tmp_verts));
      tmp_r = full(fit_data.r(tmp_verts));
    end;
    verts = [verts;tmp_verts];
    th = [th;tmp_th];
    r = [r;tmp_r];
  end;
  
  % calculate preferred location cartesian coordinates for each mapped vertex
  x = r.*cos(th);
  y = r.*sin(th);

  % calculate receptive field sizes
  rf_sizes = parms.rf_sizes(a)*ones(length(verts),1);
  if parms.rf_slopes(a)~=0
    % rf_size depends on r
    rf_sizes = rf_sizes + parms.rf_slopes(a)*(r-parms.rf_r_inter);
  else
    rf_sizes = parms.rf_sizes(a)*ones(length(verts),1);
  end;

  % calculate weights for each vertex for each location in visual field
  map = zeros(length(verts),parms.stimres.^2);
  for i=1:length(parms.stim_x)
    dx = parms.stim_x(i) - x;
    dy = parms.stim_y(i) - y;
    w = exp(-(dx.^2 + dy.^2)./(2*(rf_sizes.^2)));
    if parms.surround_flag
      rf_sizes2 = parms.surround_rf_fact*rf_sizes;
      w2 = parms.surround_amp_fact*exp(-(dx.^2 + dy.^2)./...
                                       (2*((rf_sizes2).^2)));
      w = w./(rf_sizes*sqrt(2*pi)) - w2./(rf_sizes2*sqrt(2*pi));
    end;    
    map(:,i) = map(:,i) + w;
  end;

  % normalize so positive values in each row (across visual field) sum to one
  tmp = map;
  tmp(map<0) = 0;
  rowsums = sum(tmp,2);
  rowsums(rowsums<eps) = 1;
  clear tmp;
  map = bsxfun(@rdivide,map,rowsums);

  % apply threshold
  if parms.w_thresh>0
    if parms.vfnorm_flag
      maxvals = max(eps,max(map,[],2)); % max value across visual field
    else
      maxvals = max(eps,max(map(:))); % max value for all
    end;
    tmp = bsxfun(@rdivide,map,max(eps,maxvals));
    tmp = 1.0*(abs(tmp)>=parms.w_thresh);
    map = tmp .* map;
    clear tmp;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

