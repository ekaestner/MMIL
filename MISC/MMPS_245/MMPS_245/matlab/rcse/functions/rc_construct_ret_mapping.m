function new_retmap = rc_construct_ret_mapping(retmap,varargin);
%function new_retmap = rc_construct_ret_mapping(retmap,[options]);
%
% Purpose: Construct retinotopy mapping structure for RCSE
%   defines stimulus location maps to dipoles
%
% Usage:
%  new_retmap = rc_construct_ret_mapping(retmap,'key1',value1,...);
%
% Required Input:
%  retmap - specifies vertices for visual areas and additional dipoles
%           visual area vertices can be defined explicitly with vertex numbers
%           or with wfiles (specifying a cluster of vertices)
%           (see examp_define_retmap.m and examp_define_retmap_wfiles
%            for how to generate)
%
% Optional Input:
%  'norm_weights_flag': [0|1|2] whether and how to normalize vertex weights
%     0: do not normalize
%     1: normalize so sum of weights for each stimulus location equals 1
%     2: normalize so average sum of weights (within a visual area) equals 1
%     {default: 2}
%  'use_areas': vector of area numbers defining subset of visual areas
%     defined in retmap to use
%     {default: []} (if empty, use all areas in retmap)
%  'conditions': vector of condition numbers (index to avg_data.averages)
%     defining subset of stimulus locations in retmap to use
%     If empty, use all stim locs in retmap.cond_order
%     {default: []}
%  'event_codes': vector of event codes
%     defining subset of stimulus conditions in cond_info to use
%     If specified, 'conditions' will be ignored
%     {default: []}
%  'wfiledir': directory (relative to running dir) containing FreeSurfer wfiles
%     used if retmap is defined using wfiles
%     {default: './weights'}
%  'w_thresh_patch': threshold applied to cortical patch weights
%       relative to max val for each patch
%     this reduces the number of dipoles that included in the forward matrix
%       cutting off tails of smooth, gaussian-like cluster
%     {default: 0}
%  'ret_dips_lh_dec_file': name of left hemi FreeSurfer dec file (from tksurfer)
%     for extra retinotopic dipoles -- these can be instead of or in addition
%     to extra dipoles defined in retmap
%     dec file contains 0's and 1's indicating subset of dipoles to use
%     {default: []} (relative to subject dir)
%  'ret_dips_rh_dec_file': name of right hemi dec file for retinotopic dipoles
%     {default: []} (relative to subject dir)
%  'nonret_dips_lh_dec_file': name of left hemi FreeSurfer dec file (from tksurfer)
%     for extra non-retinotopic dipoles -- these can be instead of or in addition
%     to extra dipoles defined in retmap
%     {default: []} (relative to subject dir)
%  'nonret_dips_rh_dec_file': name of right hemi dec file for non-ret dipoles
%     {default: []} (relative to subject dir)
%  'ret_dip_qfield_flag': [1|0] Toggle quarterfield constraint
%     if 1, extra retinotopic dipoles are constrained to have the
%     same source strength and orientation within a visual quarter field
%     (reduces number of unknowns and prevents these dipoles from absorbing too
%      much variance)
%     {default: 1}
%
% Created:  03/07/06 by Don Hagler
% Last Mod: 03/09/15 by Don Hagler
%

new_retmap = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'norm_weights_flag',2,[0 1 2 3],...
  'use_areas',[],[],...
  'conditions',[],[],...
  'event_codes',[],[],...
  'wfiledir','./weights',[],...
  'w_thresh_patch',0,[0,1],...
  'ret_dips_lh_dec_file',[],[],...
  'ret_dips_rh_dec_file',[],[],...
  'ret_dip_qfield_flag',true,[false true],...
  'nonret_dips_lh_dec_file',[],[],...
  'nonret_dips_rh_dec_file',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sure retmap fields exist, check for problems

areas = retmap.areas;
num_areas = length(areas);
if isfield(retmap,'orig_areas')
  orig_areas = retmap.orig_areas;
else  
  orig_areas = areas;
end;
orig_num_areas = length(orig_areas);
if isfield(retmap,'orig_cond_info')
  orig_cond_info = retmap.orig_cond_info;
  orig_cond_order = cell2mat({retmap.orig_cond_info.cond_number});
elseif isfield(retmap,'cond_info');
  orig_cond_info = retmap.cond_info;
  orig_cond_order = cell2mat({retmap.cond_info.cond_number});
else
  orig_cond_info = [];
  orig_cond_order = retmap.cond_order;
end;
num_locs = retmap.num_locs;
if isfield(retmap,'orig_num_locs')
  orig_num_locs = retmap.orig_num_locs;
else  
  orig_num_locs = num_locs;
end;

if ~isfield(retmap,'areas_vbase')
  areas_vbase = 0;
else
  areas_vbase = retmap.areas_vbase;
end;
if ~isfield(retmap,'ret_dip_vbase')
  ret_dip_vbase = 0;
else
  ret_dip_vbase = retmap.ret_dip_vbase;
end;
if ~isfield(retmap,'nonret_dip_vbase')
  nonret_dip_vbase = 0;
else
  nonret_dip_vbase = retmap.nonret_dip_vbase;
end;

if (isfield(retmap,'ret_dips'))
  ret_dips = retmap.ret_dips;
  num_ret_dips = length(ret_dips);
  if ret_dip_vbase==0
    % matlab indices start at 1 but freesurfer starts at 0
    for d=1:num_ret_dips
      ret_dips(d).v = ret_dips(d).v + 1;
    end
    ret_dip_vbase = 1;
  end;
else
  fprintf('%s: no ret_dips defined in retmap...ok\n',mfilename);
  ret_dips = [];
  num_ret_dips = 0;
end

if (isfield(retmap,'nonret_dips'))
  nonret_dips = retmap.nonret_dips;
  num_nonret_dips = length(nonret_dips);
  if nonret_dip_vbase==0
    % matlab indices start at 1 but freesurfer starts at 0
    for d=1:num_nonret_dips
      nonret_dips(d).v = nonret_dips(d).v + 1;
    end
    nonret_dip_vbase = 1;
  end;
else
  fprintf('%s: no nonret_dips defined in retmap...ok\n',mfilename);
  nonret_dips = [];
  num_nonret_dips = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modify retmap based on areas in use_areas

if isempty(parms.use_areas)
  parms.use_areas = 1:orig_num_areas;
end;
use_areas = find(ismember(1:orig_num_areas,parms.use_areas));
if isempty(use_areas)
  error('use_areas is empty');
end;
areas = orig_areas(use_areas);
num_areas = length(areas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modify retmap based on conditions or event codes

if ~isempty(parms.event_codes)
  if isempty(orig_cond_info)
    error('must supply cond_info to use event_codes option');
  end;
  event_codes = cell2mat({orig_cond_info.event_code});
  [parms.event_codes,use_conds] = ...
    intersect(event_codes,parms.event_codes);
else
  if isempty(parms.conditions)
    parms.conditions = retmap.cond_order;
  end;
  use_conds = find(ismember(retmap.cond_order,parms.conditions));
  if isempty(use_conds)
    error('condition numbers do not match those in cond_order');
  end;
end;

cond_order = orig_cond_order(use_conds);
if ~isempty(orig_cond_info)
  cond_info = orig_cond_info(use_conds);
end;

if isfield(areas, 'verts')
  for a=1:num_areas
    verts = areas(a).verts;
    areas(a).verts = areas(a).verts(use_conds);
  end
elseif isfield(areas, 'wfiles')
  for a=1:num_areas
    areas(a).wfiles = areas(a).wfiles(use_conds);
  end
  areas_vbase = 1;
else
  error('retmap lacks both verts and wfiles!');
end;
num_locs = length(cond_order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uniq_verts_lh = [];
uniq_verts_rh = [];

% create M matrix/structure
for a=1:num_areas
  for theta=1:num_locs

    if isfield(areas,'verts') % use pre-determined vertex numbers
      % left hemisphere dipoles for this location
      v_lh = areas(a).verts(theta).v_lh;
      w_lh = areas(a).verts(theta).w_lh;
      % right hemisphere dipoles for this location
      v_rh = areas(a).verts(theta).v_rh;
      w_rh = areas(a).verts(theta).w_rh;
    elseif isfield(areas,'wfiles') % read wfiles
      % left hemisphere dipoles for this location
      filename = areas(a).wfiles(theta).lh_fname;
      v_lh=[]; w_lh=[];
      if ~isempty(filename)
        fname = sprintf('%s/%s',parms.wfiledir,filename);
        [w_lh,v_lh]=fs_read_wfile(fname);
      end
      % right hemisphere dipoles for this location
      filename = areas(a).wfiles(theta).rh_fname;
      v_rh=[]; w_rh=[];
      if ~isempty(filename)
        fname = sprintf('%s/%s',parms.wfiledir,filename);
        [w_rh,v_rh]=fs_read_wfile(fname);
      end
    end;

    % threshold relative to max value
    if parms.w_thresh_patch>0
      % normalize to max across both hemispheres
      w = rc_norm_weights([w_lh;w_rh],1);
      w_lh = w(1:length(w_lh));
      w_rh = w(length(w_lh)+1:end);
      % apply threshold
      [w_lh,v_lh]=rc_thresh_weights(w_lh,v_lh,parms.w_thresh_patch);
      [w_rh,v_rh]=rc_thresh_weights(w_rh,v_rh,parms.w_thresh_patch);
    end;

    % normalize w values so sum = 1
    if parms.norm_weights_flag==1
      w=rc_norm_weights([w_lh;w_rh]);
      w_lh = w(1:length(w_lh));
      w_rh = w(length(w_lh)+1:end);
    end;

    % add weights and vertex numbers of M for lh
    M(a,theta).w_lh = w_lh;
    if areas_vbase==0
      M(a,theta).v_lh = v_lh+1;
      areas(a).verts(theta).v_lh = areas(a).verts(theta).v_lh+1;
    else
      M(a,theta).v_lh = v_lh;
    end;

    % add weights and vertex numbers of M for rh
    M(a,theta).w_rh = w_rh;
    if areas_vbase==0
      M(a,theta).v_rh = v_rh+1;
      areas(a).verts(theta).v_rh = areas(a).verts(theta).v_rh+1;
    else
      M(a,theta).v_rh = v_rh;
    end;

    % add new verts to uniq_verts lists
    uniq_verts_lh=union(uniq_verts_lh,M(a,theta).v_lh);
    uniq_verts_rh=union(uniq_verts_rh,M(a,theta).v_rh);
  end;
end;
areas_vbase = 1;

if parms.norm_weights_flag==2
  fprintf('%s: normalizing vertex weights to average of sums\n',mfilename);
  for a=1:num_areas
    w_sums = zeros(1,num_locs);
    for theta=1:num_locs
      w_sums(theta) = sum([M(a,theta).w_lh;M(a,theta).w_rh]);
    end;
    w_sum_avg = mean(w_sums); %% todo: abs to allow negative weights
    if w_sum_avg==0
      fprintf('%s: WARNING: average of vertex weight sums = 0\n',mfilename);
      w_sum_avg = 1;
    end;
    for theta=1:num_locs
      M(a,theta).w_lh = M(a,theta).w_lh/w_sum_avg;
      M(a,theta).w_rh = M(a,theta).w_rh/w_sum_avg;
    end;
  end;
end;

% add extra retinotopic dipoles
if ~isempty(ret_dips)
  lh_ret_dips = find(strcmp({ret_dips.hemisphere},'lh'));
  rh_ret_dips = find(strcmp({ret_dips.hemisphere},'rh'));
  lh_ret_dips_v = cell2mat({ret_dips(lh_ret_dips).v});
  rh_ret_dips_v = cell2mat({ret_dips(rh_ret_dips).v});
else
  lh_ret_dips_v = [];
  rh_ret_dips_v = [];
end;

if ~isempty(parms.ret_dips_lh_dec_file)
  dec_dip_lh=ts_read_dec_file(parms.ret_dips_lh_dec_file);
  tmp_v = find(dec_dip_lh);
  lh_ret_dips_v = union(lh_ret_dips_v,tmp_v);
end;
if ~isempty(parms.ret_dips_rh_dec_file)
  dec_dip_rh=ts_read_dec_file(parms.ret_dips_rh_dec_file);
  tmp_v = find(dec_dip_rh);
  rh_ret_dips_v = union(rh_ret_dips_v,tmp_v);
end;

% exclude dips that are already in uniq_verts (i.e. in area models)
keep_dips = ~ismember(lh_ret_dips_v,uniq_verts_lh);
lh_ret_dips_v = lh_ret_dips_v(keep_dips);
keep_dips = ~ismember(rh_ret_dips_v,uniq_verts_rh);
rh_ret_dips_v = rh_ret_dips_v(keep_dips);

% add to list of vertices
num_lh_ret_dips = length(lh_ret_dips_v);
num_rh_ret_dips = length(rh_ret_dips_v);
uniq_verts_lh=union(uniq_verts_lh,lh_ret_dips_v);
uniq_verts_rh=union(uniq_verts_rh,rh_ret_dips_v);

% add extra nonretinotopic dipoles
if ~isempty(nonret_dips)
  lh_nonret_dips = find(strcmp({nonret_dips.hemisphere},'lh'));
  rh_nonret_dips = find(strcmp({nonret_dips.hemisphere},'rh'));
  lh_nonret_dips_v = cell2mat({nonret_dips(lh_nonret_dips).v});
  rh_nonret_dips_v = cell2mat({nonret_dips(rh_nonret_dips).v});
else
  lh_nonret_dips_v = [];
  rh_nonret_dips_v = [];
end;
if ~isempty(parms.nonret_dips_lh_dec_file)
  dec_dip_lh=read_dec_file(parms.nonret_dips_lh_dec_file);
  tmp_v = find(dec_dip_lh);
  lh_nonret_dips_v = union(lh_nonret_dips_v,tmp_v);
end;
if ~isempty(parms.nonret_dips_rh_dec_file)
  dec_dip_rh=read_dec_file(parms.nonret_dips_rh_dec_file);
  tmp_v = find(dec_dip_rh);
  rh_nonret_dips_v = union(rh_nonret_dips_v,tmp_v);
end;

% exclude dips that are already in uniq_verts (i.e. in area models)
keep_dips = ~ismember(lh_nonret_dips_v,uniq_verts_lh);
lh_nonret_dips_v = lh_nonret_dips_v(keep_dips);
keep_dips = ~ismember(rh_nonret_dips_v,uniq_verts_rh);
rh_nonret_dips_v = rh_nonret_dips_v(keep_dips);

% add to list of vertices
num_lh_nonret_dips = length(lh_nonret_dips_v);
num_rh_nonret_dips = length(rh_nonret_dips_v);
uniq_verts_lh=union(uniq_verts_lh,lh_nonret_dips_v);
uniq_verts_rh=union(uniq_verts_rh,rh_nonret_dips_v);

% check that qfield constraint is valid
if (num_ret_dips)
  if parms.ret_dip_qfield_flag
    num_ret_dip_locs = num_locs/4;
  else
    num_ret_dip_locs = num_locs;
  end
  if (~mmil_isint(num_ret_dip_locs))
    fprintf('%s: since number of locations is not a multiple of 4,\n',mfilename);
    fprintf('      qfield constraint disabled\n',mfilename);
    num_ret_dip_locs = num_locs;
  end
else
  num_ret_dip_locs = num_locs;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect into a single structure
new_retmap.cond_order = cond_order;
new_retmap.cond_info = cond_info;
new_retmap.orig_cond_order = orig_cond_order;
new_retmap.orig_cond_info = orig_cond_info;
new_retmap.use_areas = parms.use_areas;
new_retmap.areas = areas;
new_retmap.orig_areas = orig_areas;
new_retmap.num_areas = num_areas;
new_retmap.orig_num_areas = orig_num_areas;
new_retmap.num_locs = num_locs;
new_retmap.orig_num_locs = orig_num_locs;
new_retmap.ret_dips = ret_dips;
new_retmap.lh_ret_dips_v = lh_ret_dips_v;
new_retmap.rh_ret_dips_v = rh_ret_dips_v;
new_retmap.num_lh_ret_dips = num_lh_ret_dips;
new_retmap.num_rh_ret_dips = num_rh_ret_dips;
new_retmap.num_ret_dips = num_lh_ret_dips + num_rh_ret_dips;
new_retmap.num_ret_dip_locs = num_ret_dip_locs;
new_retmap.nonret_dips = nonret_dips;
new_retmap.lh_nonret_dips_v = lh_nonret_dips_v;
new_retmap.rh_nonret_dips_v = rh_nonret_dips_v;
new_retmap.num_lh_nonret_dips = num_lh_nonret_dips;
new_retmap.num_rh_nonret_dips = num_rh_nonret_dips;
new_retmap.num_nonret_dips = num_lh_nonret_dips + num_rh_nonret_dips;
new_retmap.uniq_verts_lh = reshape(uniq_verts_lh,[prod(size(uniq_verts_lh)),1]);
new_retmap.uniq_verts_rh = reshape(uniq_verts_rh,[prod(size(uniq_verts_rh)),1]);
new_retmap.M = M;
new_retmap.areas_vbase = areas_vbase;
new_retmap.ret_dip_vbase = ret_dip_vbase;
new_retmap.nonret_dip_vbase = nonret_dip_vbase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
