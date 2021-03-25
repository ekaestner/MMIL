function [tmp_retmap] = rc_choose_retdips(retmap);
%function [tmp_retmap] = rc_choose_retdips(retmap);
%
% Note: must have run rc_construct_ret_mapping on retmap first
%
% created:  06/26/06 by Don Hagler
% modified: 02/19/11 by Don Hagler
%
%% todo: optionally input best_retmap and surf
%%   for each dipole, pick a nearest neighbor within cluster
%%   if a vertex has no neighbors in cluster, then pick any vertex in cluster   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select a single vertex from 

areas = retmap.orig_areas;
num_areas = length(areas);
num_locs = retmap.num_locs;
M = retmap.M;

for a=1:num_areas
  for theta=1:num_locs
    num_neighbors = length(M(a,theta).v_lh);
    i_v = 1 + floor(rand(1)*num_neighbors);
    if isempty(M(a,theta).v_lh)
      v = [];
      w = [];
    else
      v = M(a,theta).v_lh(i_v);
      w = 1;
    end;
    areas(a).verts(theta).v_lh = v;
    areas(a).verts(theta).w_lh = w;
        
    num_neighbors = length(M(a,theta).v_rh);
    i_v = 1 + floor(rand(1)*num_neighbors);
    if isempty(M(a,theta).v_rh)
      v = [];
      w = [];
    else
      v = M(a,theta).v_rh(i_v);
      w = 1;
    end;
    areas(a).verts(theta).v_rh = v;
    areas(a).verts(theta).w_rh = w;
  end;
end;

tmp_retmap.areas = areas;
tmp_retmap.orig_areas = areas;
tmp_retmap.cond_order = retmap.cond_order;
tmp_retmap.cond_info = retmap.cond_info;
tmp_retmap.num_locs = retmap.num_locs;
tmp_retmap.orig_num_locs = retmap.orig_num_locs;
tmp_retmap.orig_num_areas = num_areas;
tmp_retmap.ret_dips = retmap.ret_dips;
tmp_retmap.nonret_dips = retmap.nonret_dips;
tmp_retmap.areas_vbase = retmap.areas_vbase;
tmp_retmap.ret_dip_vbase = retmap.ret_dip_vbase;
tmp_retmap.nonret_dip_vbase = retmap.nonret_dip_vbase;
tmp_retmap.num_areas = length(areas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

status = 1;
return;
