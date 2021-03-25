function retmap=rc_define_retmap_from_retfit(retfit_results,varargin)
%function retmap=rc_define_retmap_from_retfit(retfit_results,[options])
%
% Usage:
%  retmap=rc_define_retmap_from_retfit(retfit_results,'key1', value1,...);
%
% Required parameters:
%   retfit_results: struct containing these fields (output from retfit)
%     lh_data: data struct for left hemisphere
%     rh_data: data struct for right hemisphere
%     lh_fit_data: fitted data struct for left hemisphere
%     rh_fit_data: fitted data struct for right hemisphere
%     lh_area_masks: area masks struct for left hemisphere
%     rh_area_masks: area masks struct for right hemisphere
%
% Optional parameters:
%  'cond_info': struct array containing condition information
%    see rc_read_cond_info
%    If supplied, fname_conds will be ignored
%    {default: []}
%  'fname_conds': full path of csv file containing condition information
%    {default: []}
%  'contrasts': vector of acceptable contrast levels
%    if empty, use all contrasts specified in fname_conds
%    ignored if fname_conds is empty
%    {default: []}
%  'r_vec': vector of eccentricities (degrees visual angle)
%    ignored if fname_conds is not empty
%    {default = [5,7]}
%  'th_vec' vector of polar angles (degrees)
%    ignored if fname_conds is not empty
%    {default = [45,135,225,315]}
%  'ecc_width': eccentricity width of stimuli (deg. vis. ang.)
%    will be ignored if fname_conds has ecc_width column
%    {default = 1}
%  'theta_width': polar angle width of stimuli (degrees)
%    will be ignored if fname_conds has theta_width column or if stim_type=1
%    {default = 10}
%  'area_names': cell array containing area names
%    {default = {'v1','v2','v3'}}
%  'w_thresh': threshold applied to weights relative to max
%    {default = 0.01}
%  'single_vertex_flag': [0|1] select one vertex for each stimulus location
%    0: use all vertices
%    1: vertex with maximum weight
%    2: vertex closest to center of mass
%    {default = 0}
%  'forward': struct containing RCSE forward matrix and other things
%    including lh_dip_info and rh_dip_info
%    required if single_vertex_flag = 2
%    {default = []}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'r_offset': value (degrees visual angle) added to
%    eccentricities in fname_conds or r_vec
%    {default = 0}
%  'th_offset': degrees polar angle added to all
%    polar angles in fname_conds or th_vec
%    {default = 0}
%  'grid_offset_flag': instead of r and th, make offets to grid
%    coordinates u and v
%    {default = 0}
%  'rf_sizes': vector of receptive field sizes (degrees visual angle)
%    -- one for each visual area
%    {default = [1,1,1]}
%  'rf_slopes': vector of slopes of linear trend of receptive field sizes
%    w.r.t. ecc for each visual area
%    Intercept is assumed to be half of r_max
%    {default = [0.1,0.1,0.1]}
%  'smoothed_areas_flag': [0|1] whether to use smoothed area masks
%    {default = 0}
%  'restrict_hemi_flag': [0|1] whether to restrict dipole clusters
%    so no ipsilateral
%    {default = 0}
%  'restrict_uplow_flag': [0|1] whether to restrict dipole clusters
%    so no upper-lower cross-over
%    {default = 0}
%  'data_flag': [0|1] whether to use refit area ROIs but original
%    retinotopy data instead of template values for selection of dipole clusters
%    {default = 0}
%  'vf2ctx': struct array containing lh and rh mapping matrices
%    {default = []}
%  'stim_type': type of stimulus to model
%    0: point: generation of retmap is very fast, but less accurate
%    1: circle: slower, but more accurate
%    2: wedge: slower, but more accurate
%    {default = 2}
%  'point_stim_scalefact': fractional scaling of point stimulus size
%    {default = 0.67}
%
%
% Created:  12/18/08 by Don Hagler
% Last Mod: 12/08/13 by Don Hagler
%

%% todo: use_areas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms( varargin, {...
  'cond_info',[],[],...
  'fname_conds',[],[],...
  'contrasts',[],[],...
  'r_vec',[5,7],[],...
  'th_vec',[45,135,225,315],[],...
  'ecc_width',1,[0,100],...
  'theta_width',10,[0,360],...
  'w_thresh',0.01,[0,100],...
  'single_vertex_flag',0,[0 1 2],...
  'forward',[],[],...
  'r_max',12.5,[0,Inf],...
  'r_offset',0,[-100,100],...
  'th_offset',0,[-180,180],...
  'grid_offset_flag',false,[false true],...
  'rf_sizes',[1,1,1],[],...
  'rf_slopes',[0.1,0.1,0.1],[0,10],...
  'smoothed_areas_flag',false,[false true],...
  'restrict_hemi_flag',false,[false true],...
  'restrict_uplow_flag',false,[false true],...
  'data_flag',false,[false true],...
  'vf2ctx',[],[],...
  'stim_type',2,[0,1,2],... 
  'point_stim_scalefact',0.05,[0.001,1],... % only affects stim_type = 0
...
  'area_names',{'v1','v2','v3'},[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'stimres',100,[50,1000],...
  'polar_flag',false,[false true],... % only affects stim_type = 0
...
  'nbrhds',[],[],...
});

if ~isempty(parms.vf2ctx) && parms.stim_type==0
  parms.stim_type = 1;
end;

if parms.stim_type>0
  parms.polar_flag = false;
end;

if parms.single_vertex_flag==2
  if isempty(parms.forward)
    error('forward struct must be supplied if single_vertex_flag=2');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

retmap = [];
cond_info = []; % should contain struct array with ecc, theta, for each condition
areas = []; % should have names for each visual area, verts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms
if ~iscell(parms.area_names), parms.area_names = {parms.area_names}; end;
num_areas = length(parms.area_names);

if ~isempty(parms.fname_conds)
  if ~exist(parms.fname_conds)
    error('file %s not found',parms.fname_conds);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get condition information

if ~isempty(parms.cond_info) || ~isempty(parms.fname_conds)
  if ~isempty(parms.cond_info)
    cond_info = parms.cond_info;
  else
    cond_info = rc_read_cond_info(parms.fname_conds);
  end;
  if ~isempty(parms.contrasts) && isfield(cond_info,'contrast')
    contrasts = cell2mat({cond_info.contrast});
    cond_order = find(ismember(contrasts,parms.contrasts));
    cond_info = cond_info(cond_order);
  else
    cond_order = 1:length(cond_info);
  end;
else
  c=1;
  for i=1:length(parms.r_vec)
    r = parms.r_vec(i);
    for j=1:length(parms.th_vec)
      th = parms.th_vec(j);
      if th<0, th=th+360; end;
      cond_info(c).ecc = r;
      cond_info(c).theta = th;
      c=c+1;
    end;
  end;
  cond_order = 1:length(cond_info);
end;
nconds = length(cond_info);

offset_infix_list = {''};
if ~parms.restrict_hemi_flag
  offset_infix_list{end+1} = '_ipsi';
end;
if ~parms.restrict_uplow_flag
  offset_infix_list{end+1} = '_cross';
end;
if ~parms.restrict_hemi_flag & ~parms.restrict_uplow_flag
  offset_infix_list{end+1} = '_ipsi_cross';
end;
npatches = length(offset_infix_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offset_types = {'r','th'};
for t=1:length(offset_types)
  fieldname = [offset_types{t} '_offset'];
  tmp_offsets = parms.(fieldname);
  num_offsets = numel(tmp_offsets);
  if num_offsets==1
    parms.(fieldname) = tmp_offsets*ones(num_areas,nconds,npatches);
  elseif num_offsets==nconds
    tmp_offsets = reshape(tmp_offsets,[1,nconds]);
    parms.(fieldname) = repmat(tmp_offsets,[num_areas,1,npatches]);
  elseif num_offsets==(num_areas*nconds)
    tmp_offsets = reshape(tmp_offsets,[num_areas,nconds]);
    parms.(fieldname) = repmat(tmp_offsets,[1,1,npatches]);
  elseif num_offsets==(num_areas*nconds*npatches)
    if numel(size(tmp_offsets))~=3
      parms.(fieldname) = reshape(tmp_offsets,[num_areas,nconds,npatches]);
    end;
  else
    error('size of %s should be 1x%d or %dx%d or %dx%dx%d',...
      fieldname,nconds,num_areas,nconds,num_areas_nconds,npatches);
  end;
end;

% area, condition, and patch specific offsets from cond_info
for a=1:num_areas
  for t=1:length(offset_types)      
    dest_fieldname = [offset_types{t} '_offset'];
    for p=1:npatches
      fieldname = sprintf('%s_offset_%s%s',...
        offset_types{t},parms.area_names{a},offset_infix_list{p});
      if isfield(parms.cond_info,fieldname)
        for i=1:nconds
          tmp_offset = mmil_getfield(cond_info(i),fieldname);
          if ~isempty(tmp_offset)
            parms.(dest_fieldname)(a,i,p) = ...
              parms.(dest_fieldname)(a,i,p) + tmp_offset;
          end;
        end;
      end;
    end;
  end;
end;

% condition specific ecc_width and theta_width
parms.ecc_widths = zeros(nconds,1);
fieldname = 'ecc_width';
for i=1:nconds
  tmp_width = mmil_getfield(cond_info(i),fieldname);
  if isempty(tmp_width), tmp_width = parms.ecc_width; end;
  parms.ecc_widths(i) = tmp_width;
end;

parms.theta_widths = zeros(nconds,1);
fieldname = 'theta_width';
for i=1:nconds
  tmp_width = mmil_getfield(cond_info(i),fieldname);
  if isempty(tmp_width), tmp_width = parms.theta_width; end;
  parms.theta_widths(i) = tmp_width;
end;

% find unique locations from cond_info
parms.same_location_conds = [1:nconds]';
if ~isempty(cond_info)
  if isfield(cond_info,'contrast')
    parms.contrasts = [cond_info.contrast];
  else
    parms.contrasts = ones(1,length(cond_info));
  end;
  parms.unique_contrasts = unique(parms.contrasts);
  parms.ncontrasts = length(parms.unique_contrasts);
  if parms.ncontrasts>1
    for i=nconds:-1:2
      r = cond_info(i).ecc;
      th = cond_info(i).theta;
      for j=1:i-1
        tmp_r = cond_info(j).ecc;
        tmp_th = cond_info(j).theta;
        if tmp_r==r & tmp_th==th
          parms.same_location_conds(i) = j;
          break;
        end;
      end;
    end;
  end;
  parms.unique_location_conds = unique(parms.same_location_conds);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define retmap

fprintf('%s: defining dipoles for all areas...\n',mfilename);
maxR = parms.stimres/2;
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  % fitted data on cortical surface
  if parms.data_flag
    fit_data = retfit_results.([hemi '_data']);
  else
    fit_data = retfit_results.([hemi '_fit_data']);
  end;
  if ~isfield(fit_data,'x') && ~parms.polar_flag
    fit_data.x = fit_data.r.*cos(fit_data.th);
    fit_data.y = fit_data.r.*sin(fit_data.th);
  end;
  area_masks = retfit_results.([hemi '_area_masks']);
  area_names = {area_masks.name};
  area_data = retfit_results.([hemi '_area_data']);
  if parms.grid_offset_flag
    roi = find(fit_data.u | fit_data.v);
    if isempty(area_data)
      all_u = full(fit_data.u);
      all_v = full(fit_data.v);
      roi_u = all_u(roi);
      roi_v = all_v(roi);
    end;
  end;

  % substitute smoothed area masks
  if parms.smoothed_areas_flag
    for i=1:length(area_masks)
      area_masks(i).vertices = area_masks(i).vertices_smoothed;
    end;
  end;

  % visual field to cortex mapping matrix
  if isempty(parms.vf2ctx)
    % calculate rf sizes for vertices in each area mask
    for i=1:length(area_masks)
      tmp_verts = area_masks(i).vertices;
      area_name = area_masks(i).name;
      tmp_name = regexprep(area_name,'[+-]','');
      a = find(strcmp(tmp_name,parms.area_names));
      if isempty(a), continue; end;
      area_masks(i).rf_sizes = parms.rf_sizes(a)*...
                               ones(size(tmp_verts));
      if parms.rf_slopes(a)~=0
        % rf_size depends on fit_data.r
        r0 = parms.r_max/2;
        r = fit_data.r(tmp_verts);
        area_masks(i).rf_sizes = area_masks(i).rf_sizes +...
                                 parms.rf_slopes(a)*(r-r0);
      end;
    end;
  end;

  for a=1:num_areas
    if ~isempty(parms.vf2ctx)
      vf2ctx = parms.vf2ctx(a,h).map;
      verts = parms.vf2ctx(a,h).verts;
    else
      vf2ctx = [];
    end;

    % get u and v from area_data if it exists
    if parms.grid_offset_flag && ~isempty(area_data)
      all_u = full(area_data(a).u);
      all_v = full(area_data(a).v);
      roi_u = all_u(roi);
      roi_v = all_v(roi);
    end;

    % match area names
    areas(a).name = parms.area_names{a};
    ind = find(strcmp(areas(a).name,area_names));

    % get vertex numbers and rf_sizes from area masks
    if ~isempty(ind)
      mask_lower = area_masks(ind).vertices;
      mask_upper = mask_lower;
      if isempty(vf2ctx)
        rf_sizes_lower = area_masks(ind).rf_sizes;
        rf_sizes_upper = rf_sizes_lower;
      end;
    else % check for names with - or + in them
      mask_name = [areas(a).name '-'];
      ind = find(strcmp(mask_name,area_names));
      if isempty(ind)
        error('mask name %s not found in retfit_results.%s_area_masks',...
          mask_name,hemi);
      end;
      mask_lower = area_masks(ind).vertices;
      if isempty(vf2ctx)
        rf_sizes_lower = area_masks(ind).rf_sizes;
      end;
      mask_name = [areas(a).name '+'];
      ind = find(strcmp(mask_name,area_names));
      if isempty(ind)
        error('mask name %s not found in retfit_results.%s_area_masks',...
          mask_name,hemi);
      end;
      mask_upper = area_masks(ind).vertices;
      if isempty(vf2ctx)
        rf_sizes_upper = area_masks(ind).rf_sizes;
      end;
    end;

    for i=1:nconds
      if ~ismember(i,parms.unique_location_conds)
        j = parms.same_location_conds(i);
        areas(a).verts(i).(['v_' hemi]) = areas(a).verts(j).(['v_' hemi]);
        areas(a).verts(i).(['w_' hemi]) = areas(a).verts(j).(['w_' hemi]);
        continue;
      else
        areas(a).verts(i).(['v_' hemi]) = [];
        areas(a).verts(i).(['w_' hemi]) = [];
        r = cond_info(i).ecc;
        th = cond_info(i).theta;
        dr = parms.ecc_widths(i);
        dth = parms.theta_widths(i);
        if dr==0, continue; end;
        if th>90 && th<270 % left visual field
          hemifield = 'lh';
        else
          hemifield = 'rh';
        end;
        % disallow cross-over between hemispheres
        if parms.restrict_hemi_flag & strcmp(hemi,hemifield), continue; end;

        if th>0 && th<=180
          stim_uplow = 1;
        else
          stim_uplow = 2;;
        end;

        for uplow=1:2
          % disallow cross-over between upper and lower sub areas
          if parms.restrict_uplow_flag & uplow~=stim_uplow, continue; end;

          switch uplow
            case 1
              mask = mask_upper;
              if isempty(vf2ctx)
                rf_sizes = rf_sizes_upper;
              end;
            case 2
              mask = mask_lower;
              if isempty(vf2ctx)
                rf_sizes = rf_sizes_lower;
              end;
          end;
          
          % determine which patch we are creating
          if ~strcmp(hemi,hemifield) & uplow==stim_uplow
            p=1;
          elseif strcmp(hemi,hemifield) & uplow==stim_uplow
            p=2;
          elseif parms.restrict_hemi_flag & uplow~=stim_uplow
            p=2;
          elseif ~strcmp(hemi,hemifield) & uplow~=stim_uplow
            p=3;
          else
            p=4;
          end;

          % apply r and th offsets to nominal visual field location
          if ~parms.grid_offset_flag
            r = r + parms.r_offset(a,i,p);
            th = th + parms.th_offset(a,i,p);
          end;
          % define stimulus mask
          if parms.stim_type
            stim_mask = define_stimulus(parms,r,th,dr,dth);
          end;

          % calculate Gaussian weights
          if parms.stim_type==0 % point
            % convert from polar degrees to degrees visual angle
            tmp_dth = 2*r*sin(dth*(pi/180)/2);
                             % /2 for right triangles *2 for both

            % account for rf sizes for each vertex
            tmp_dr = parms.point_stim_scalefact*dr + rf_sizes;
            tmp_dth = parms.point_stim_scalefact*tmp_dth + rf_sizes;

            if parms.polar_flag
              % convert back to polar angle
              tmp_dth = 2*asin(tmp_dth/(2*r))*180/pi;
              % calculate difference in polar angle in degrees
              th_diff = th - fit_data.th(mask)*180/pi;
              % correct for wrap-around
              th_diff(th_diff>180) = th_diff(th_diff>180)-360;
              th_diff(th_diff<-180) = th_diff(th_diff<-180)+360;
              % calculate difference in eccentricity in visual angle
              r_diff = r-fit_data.r(mask);
              w = exp((-th_diff.^2 ./ (2*tmp_dth.^2)) + ...
                      (-r_diff.^2  ./ (2*tmp_dr.^2 )));
            else
              % calculate difference in cartesian coordinates
              x = r*cos(th*pi/180);
              y = r*sin(th*pi/180);
              x_diff = x - fit_data.x(mask);
              y_diff = y - fit_data.y(mask);
              w = exp(-(x_diff.^2 + y_diff.^2)./(2*tmp_dr.^2));
            end;
            % apply threshold
            tmp_thresh = parms.w_thresh*max(w);
            w(abs(w)<tmp_thresh) = 0;
          elseif ~isempty(vf2ctx)
            stim_vec = reshape(stim_mask,[1,prod(size(stim_mask))]);
            w = (stim_vec * vf2ctx')';
            [tmp,ind_mask] = intersect(verts,mask);
            w = w(ind_mask);
          else
            ind_stim_mask = find(stim_mask);
            c = mean(1:parms.stimres);
            w = zeros(size(mask));
            for s=1:length(ind_stim_mask)
              [y,x] = ind2sub(size(stim_mask),ind_stim_mask(s));
              x = (x - c)*parms.r_max/maxR;
              y = (c - y)*parms.r_max/maxR;
              % calculate difference in cartesian coordinates
              dx = x - fit_data.x(mask);
              dy = y - fit_data.y(mask);
              % calculate Gaussian weights
              tmp_w = exp(-(dx.^2 + dy.^2)./(2*(rf_sizes.^2)));
              tmp_w(abs(tmp_w)<parms.w_thresh)=0;
              w = w + tmp_w;
            end;
            w = w / length(ind_stim_mask);
            % NOTE: this is different from how it is done with vf2ctx
            %  where normalization is across entire visual field
            %  convolved with receptive field
          end;

          % exclude ipsilateral phases
          if strcmp(hemi,'rh')
            ind = find(fit_data.th(mask)<=pi/2 & fit_data.th(mask)>=-pi/2);
          else
            ind = find(fit_data.th(mask)>=pi/2 | fit_data.th(mask)<=-pi/2);
          end;
          w(ind) = 0;

          % get vertices with non-zero weights
          v = mask(w~=0);
          w = w(w~=0);

          if ~isempty(w)
            % select single vertex
            switch parms.single_vertex_flag
              case 1 % vertex with max weight
                [tmp,ind] = max(w);
                v = v(ind(1));
                w = 1;
              case 2 % vertex closest to center of mass
                switch hemi
                  case 'lh'
                    dip_locs = squeeze(parms.forward.lh_dip_info(1:3,:));
                  case 'rh'
                    dip_locs = squeeze(parms.forward.rh_dip_info(1:3,:));
                end;
                [x,y,z] = rc_centmass(dip_locs(1,v),dip_locs(2,v),dip_locs(3,v),w');
                cmass = [x,y,z]';
                dist = sqrt(sum(((cmass*ones(1,length(v)) - dip_locs(:,v)).^2),1));
                [tmp,ind] = min(dist);
                v = v(ind(1));
                w = 1;
            end;

            % apply grid offset
            if parms.grid_offset_flag
              du = parms.r_offset(a,i,p);
              dv = parms.th_offset(a,i,p);
              if abs(du)>0 | abs(dv)>0
                patch_u = all_u(v) + du;
                patch_v = all_v(v) + dv;
                % find nearest roi vertices
                new_v = v;
                if ~isempty(parms.nbrhds)
                  % use pre-calculated neighbors if supplied
                  for k=1:length(v)
                    ind_roi = find(v(k)==roi);
                    if ~isempty(ind_roi)
                      nbrs = parms.nbrhds(h).n{ind_roi};
                      nbr_u = all_u(nbrs);
                      nbr_v = all_v(nbrs);
                      dist = hypot(nbr_u - patch_u(k),nbr_v - patch_v(k));                
                      [tmp,ind] = min(dist);
                      ind = find(roi==nbrs(ind));
                      new_v(k) = roi(ind);
                    else
                      fprintf('%s: WARNING: grid offset to nowhere\n',mfilename);
                    end;
                  end;
                else
                  for k=1:length(v)
                    dist = hypot(roi_u - patch_u(k),roi_v - patch_v(k));                
                    [tmp,ind] = min(dist);
                    new_v(k) = roi(ind);
                  end;
                end;
                v = new_v;
              end;
            end;

            % combine across upper and lower
            old_v = areas(a).verts(i).(['v_' hemi]);
            old_w = areas(a).verts(i).(['w_' hemi]);
            if isempty(old_v)
              new_v = v;
              new_w = w;
            else
              new_v = [v;old_v];
              new_w = [w;old_w];
              [new_v,ind] = unique(new_v);
              new_w = new_w(ind);
            end;
            areas(a).verts(i).(['v_' hemi]) = new_v;
            areas(a).verts(i).(['w_' hemi]) = new_w;
          end;
        end;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

retmap.cond_info = cond_info;
retmap.cond_order = cond_order;
retmap.num_locs = nconds;
retmap.areas = areas;
retmap.num_areas = num_areas;
retmap.ret_dips = [];
retmap.nonret_dips = [];
retmap.areas_vbase = 1;
 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stim_mask=define_stimulus(parms,r,th,dr,dth)
  maxR = parms.stimres/2;

  % define stimulus
  switch parms.stim_type
    case 1 % circle
      % calculate stim locations in cartesian coordinates
      % stim_width = dr;
      x = r*cos(th*pi/180);
      y = r*sin(th*pi/180);
      % setup grid
      [X,Y] = meshgrid(1:parms.stimres);
      dx = ((X-1) - maxR)*parms.r_max/maxR - x;
      dy = (parms.stimres - (Y-1) - maxR)*parms.r_max/maxR - y;
      stim_mask = exp(-(dx.^2 + dy.^2)/(2*dr.^2))/(2*exp(-0.5));
      stim_mask(stim_mask<0.5) = 0;
      stim_mask(stim_mask>=0.5) = 1;
    case 2 % wedge
      tmp_dr = dr;
      tmp_dth = dth;
      % calculate stim locations in polar coordinates
      rA = r - dr/2;
      rB = r + dr/2;
      thA = th - dth/2;
      thB = th + dth/2;
      % convert radius from degrees visual angle to fractional distance
      rA = rA/parms.r_max;
      rB = rB/parms.r_max;
      % convert theta from degrees to radians
      thA = (pi/180) * thA;
      thB = (pi/180) * thB;
      % correct for wrap around
      thA = mod(thA,2*pi);
      thB = mod(thB,2*pi);
      % setup grid
      [X,Y] = meshgrid(1:parms.stimres);
      dX = X-maxR; dY = Y-maxR;
      Theta = mod(atan2(dX,dY) - pi/2,2*pi);
      R = hypot(dX,dY);
      normR = (R)./maxR;
      % make wedge
      if thB <= thA   % wrapped around 0 degrees
          wedge_th = Theta >= thA | Theta <= thB;
      else
          wedge_th = Theta > thA & Theta < thB;
      end
      if rB < rA   % wrapped around 0 degrees
          wedge_r = (normR > rA & normR<1);
      else
          wedge_r = normR > rA & normR < rB;
      end
      stim_mask = wedge_th & wedge_r;
  end;
  % imagesc(stim_mask); colormap gray; drawnow; pause(0.4);

return;
