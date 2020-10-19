function F = rc_construct_ret_forward(G_norm,G_tang1,G_tang2,...
                                    retmap, num_dips_lh, num_dips_rh,...
                                    indy_locs_flag,loose_flag,near_nbr_weight);
%function F = rc_construct_ret_forward(G_norm,G_tang1,G_tang2,...
%                                    retmap, num_dips_lh, num_dips_rh,
%                                    indy_locs_flag,loose_flag,near_nbr_weight);
%
% Purpose:
%  Construct a retinotopy-constrained forward matrix that makes assumptions
%   about the relationships between responses within a given visual area
%   to different stimulus locations
%
% Required Input:
%  G_norm: forward matrix of sensor amplitudes for each decimated dipole
%          (along normal vector)
%  G_tang1: forward matrix of sensor amplitudes for each decimated dipole
%          (along first of two orthogonal tangential componenets)
%  G_tang2: forward matrix of sensor amplitudes for each decimated dipole
%          (along second of two orthogonal tangential componenets)
%  retmap: output of rc_construct_ret_mapping
%  num_dips_lh: number of left hemisphere dipoles in surface model
%  num_dips_rh: number of right hemisphere dipoles in surface model
%
% Optional Input:
%  indy_locs_flag: [1|0] Toggle calculation of independent source estimates for
%    each stimulus location
%    {default: 0}
%  loose_flag: [1|0] Toggle loose orientation constraint
%    indy_locs_flag must = 1
%    {default: 0}
%  near_nbr_weight: if non-zero, sum visual area dipoles across nearest
%    neighbors and apply this weighting to the nearest neighbors (central
%    dipole will get weight=1)
%    {default: 0}
%
% Created:  03/07/06 by Don Hagler
% Last Mod: 03/09/15 by Don Hagler
%

if ~mmil_check_nargs(nargin,6), return; end;

F=[];

if ~exist('indy_locs_flag','var') || isempty(indy_locs_flag), indy_locs_flag=0; end;
if ~exist('loose_flag','var') || isempty(loose_flag), loose_flag=0; end;
if ~exist('near_nbr_weight','var'), near_nbr_weight=[]; end;
if isempty(near_nbr_weight), near_nbr_weight=0; end;

if loose_flag & ~indy_locs_flag,
  error('loose_flag=1 so indy_locs_flag must=1 too');
end;

if loose_flag ||...
   retmap.num_lh_ret_dips>0 || retmap.num_rh_ret_dips>0 ||...
   retmap.num_lh_nonret_dips>0 || retmap.num_rh_nonret_dips>0
  % only spend time on tangential components if needed
  tang_flag = 1;
else
  tang_flag = 0;
end;

fprintf('%s: creating retinotopy forward matrix...\n',mfilename);
[num_sensors,num_sources]=size(G_norm);
fprintf('%s: size(G_norm) = [%d sensors,%d sources]\n',...
  mfilename,num_sensors,num_sources);

if loose_flag
  comp_per_source = 3;
else
  comp_per_source = 1;
end;

if indy_locs_flag
  num_sources = comp_per_source*retmap.num_areas*retmap.num_locs;
else
  num_sources = retmap.num_areas;
end;

if isfield(retmap,'orig_uniq_verts_lh')
  [c,ia_lh,ib_lh] = intersect(retmap.uniq_verts_lh,retmap.orig_uniq_verts_lh);
  [c,ia_rh,ib_rh] = intersect(retmap.uniq_verts_rh,retmap.orig_uniq_verts_rh);
  ind_verts = [ib_lh;length(retmap.orig_uniq_verts_lh)+ib_rh];
  % ia_lh and ia_rh are indices to used forward matrix
  % ib_lh and ib_rh are indices to original forward matrix
  % ind_verts is index to all vertices in original forward matrix
  tmp_G_norm = G_norm(:,ind_verts);
  if tang_flag
    tmp_G_tang1 = G_tang1(:,ind_verts);
    tmp_G_tang2 = G_tang2(:,ind_verts);
  end;
else
  tmp_G_norm = G_norm;
  if tang_flag
    tmp_G_tang1 = G_tang1;
    tmp_G_tang2 = G_tang2;
  end;
end;

% for ret and extra dips: multiply by 3 to allow dipole orientation to freely vary
% for ret dips: multiply by number of constrained locations
num_sources = num_sources + ...
  3*retmap.num_ret_dip_locs*retmap.num_ret_dips + ...
  3*retmap.num_nonret_dips;

num_measurements = num_sensors*retmap.num_locs;
F = zeros(num_measurements,num_sources);
% populate F matrix matrix
i=1;
k=0;
for a=1:retmap.num_areas
  % calculate forward solution for each location
  for theta=1:retmap.num_locs
    w_lh = zeros(num_dips_lh,1);
    w_rh = zeros(num_dips_rh,1);
    w_lh(retmap.M(a,theta).v_lh) = retmap.M(a,theta).w_lh;
    w_rh(retmap.M(a,theta).v_rh) = retmap.M(a,theta).w_rh;

    % nearest neighbor weighting
    if near_nbr_weight
      nbrw = [near_nbr_weight 1 near_nbr_weight];
      nbrw = nbrw/sum(nbrw);
      % add nearest neighbors
      if(theta==1)
        t = [retmap.num_locs,theta+1];
      elseif(theta==retmap.num_locs)
        t = [theta-1,1];
      else
        t = [theta-1,theta+1];
      end;

      % left hemisphere vertices
      w_lh = w_lh*nbrw(2);
      w_lh(retmap.M(a,t(1)).v_lh) = ...
        w_lh(retmap.M(a,t(1)).v_lh) + nbrw(1)*retmap.M(a,t(1)).w_lh;
      w_lh(retmap.M(a,t(2)).v_lh) = ...
        w_lh(retmap.M(a,t(2)).v_lh) + nbrw(3)*retmap.M(a,t(2)).w_lh;
      % right hemisphere vertices
      w_rh = w_rh*nbrw(2);
      w_rh(retmap.M(a,t(1)).v_rh) = ...
        w_rh(retmap.M(a,t(1)).v_rh) + nbrw(1)*retmap.M(a,t(1)).w_rh;
      w_rh(retmap.M(a,t(1)).v_rh) = ...
        w_rh(retmap.M(a,t(1)).v_rh) + nbrw(3)*retmap.M(a,t(1)).w_rh;
    end;

    if isfield(retmap,'orig_uniq_verts_lh')
      w_dec_dip = [w_lh(retmap.uniq_verts_lh(ia_lh));...
                   w_rh(retmap.uniq_verts_rh(ia_rh))];
    else
      w_dec_dip = [w_lh(retmap.uniq_verts_lh);w_rh(retmap.uniq_verts_rh)];
    end;
    
    % multiply by forward matrix to get sensor amplitudes
    % for this area and location
    if loose_flag
      f = zeros(num_sensors,3);
      %% todo: what if tmp_G_norm has more dips because of orig_areas in retmap?
      f(:,1) = tmp_G_norm*w_dec_dip;
      f(:,2) = tmp_G_tang1*w_dec_dip;
      f(:,3) = tmp_G_tang2*w_dec_dip;
    else
      f = tmp_G_norm*w_dec_dip;
    end;

    % insert f into F matrix
    row1 = 1 + (theta-1)*num_sensors;
    row2 = row1 + num_sensors - 1;

    if loose_flag
      F(row1:row2,i:i+2) = f;
    else
      F(row1:row2,i) = f;
    end;

    if indy_locs_flag
      i=i+comp_per_source; % different source for each stim loc
    end;
  end;
  if ~indy_locs_flag
    i=i+1; % different source for each area
  end;   
end

% setup forward solution for ret dipoles
% left hemisphere
for d=1:retmap.num_lh_ret_dips
  v = retmap.lh_ret_dips_v(d);
  j = find(retmap.uniq_verts_lh==v);
  
  if(isempty(j))
    error('bad vertex number for lh ret dip %d (v=%d)',...
      d,v-1);
  end
  f = zeros(num_sensors,3);
  f(:,1) = G_norm(:,j);
  f(:,2) = G_tang1(:,j);
  f(:,3) = G_tang2(:,j);

  % insert f into F matrix for every theta
  for theta=1:retmap.num_locs
    row1 = 1 + (theta-1)*num_sensors;
    row2 = row1 + num_sensors - 1;
    F(row1:row2,i:i+2) = f;
    % new column for each quarterfield
    if(mod(theta,round(retmap.num_locs/retmap.num_ret_dip_locs))==0)
      % increment column index for F
      i = i + 3;
    end
  end
end;
% right hemisphere
for d=1:retmap.num_rh_ret_dips
  v = retmap.rh_ret_dips_v(d);
  j = length(retmap.uniq_verts_lh) + find(retmap.uniq_verts_rh==v);
  if(isempty(j))
    error('bad vertex number for rh ret dip %d (v=%d)',d,v-1);
  end
  f = zeros(num_sensors,3);
  f(:,1) = G_norm(:,j);
  f(:,2) = G_tang1(:,j);
  f(:,3) = G_tang2(:,j);

  % insert f into F matrix for every theta
  for theta=1:retmap.num_locs
    row1 = 1 + (theta-1)*num_sensors;
    row2 = row1 + num_sensors - 1;
    F(row1:row2,i:i+2) = f;
    % new column for each quarterfield
    if(mod(theta,round(retmap.num_locs/retmap.num_ret_dip_locs))==0)
      % increment column index for F
      i = i + 3;
    end
  end
end;


% setup forward solution for nonret dipoles
% left hemisphere
for d=1:retmap.num_lh_nonret_dips
  v = retmap.lh_nonret_dips_v(d);
  j = find(retmap.uniq_verts_lh==v);
  if(isempty(j))
    error('bad vertex number for lh nonret dip %d (v=%d)',d,v-1);
  end
  f = zeros(num_sensors,3);
  f(:,1) = G_norm(:,j);
  f(:,2) = G_tang1(:,j);
  f(:,3) = G_tang2(:,j);

  % insert f into F matrix for every theta
  for theta=1:retmap.num_locs
    row1 = 1 + (theta-1)*num_sensors;
    row2 = row1 + num_sensors - 1;
    F(row1:row2,i:i+2) = f;
    % new column for each quarterfield
    if(mod(theta,round(retmap.num_locs))==0)
      % increment column index for F
      i = i + 3;
    end
  end
end;
% right hemisphere
for d=1:retmap.num_rh_nonret_dips
  v = retmap.rh_nonret_dips_v(d);
  j = length(retmap.uniq_verts_lh) + find(retmap.uniq_verts_rh==v);
  if(isempty(j))
    error('bad vertex number for rh nonret dip %d (v=%d)',d,v-1);
  end
  f = zeros(num_sensors,3);
  f(:,1) = G_norm(:,j);
  f(:,2) = G_tang1(:,j);
  f(:,3) = G_tang2(:,j);

  % insert f into F matrix for every theta
  for theta=1:retmap.num_locs
    row1 = 1 + (theta-1)*num_sensors;
    row2 = row1 + num_sensors - 1;
    F(row1:row2,i:i+2) = f;
    % new column for each quarterfield
    if(mod(theta,round(retmap.num_locs))==0)
      % increment column index for F
      i = i + 3;
    end
  end
end;

fprintf('%s: size(F) = [%d measurements,%d sources]\n',...
  mfilename,num_measurements,num_sources);

return;
