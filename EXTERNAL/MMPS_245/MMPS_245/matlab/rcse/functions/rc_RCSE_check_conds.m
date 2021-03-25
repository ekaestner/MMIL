function parms = rc_RCSE_check_conds(parms,avg_data)
%function parms = rc_RCSE_check_conds(parms,avg_data)
%
% Purpose: check conditions, etc.
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 05/06/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

% conditions to use
if ~isempty(parms.event_codes)
  if isempty(parms.cond_info)
    error('must supply cond_info to use event_codes option');
  end;
  event_codes = cell2mat({parms.cond_info.event_code});
  [parms.event_codes,parms.conditions] = ...
    intersect(event_codes,parms.event_codes);
elseif ~isempty(parms.conditions)
  conditions = 1:length(avg_data.averages);
  parms.conditions = intersect(parms.conditions,conditions);
else
  parms.conditions = 1:length(avg_data.averages);
end;

% reduce cond_info to selected conditions / event_codes
if ~isempty(parms.retfit_results)
  if isfield(parms.cond_info,'contrast')
    contrasts = cell2mat({parms.cond_info.contrast});
    parms.conditions = intersect(find(contrasts>0),parms.conditions);
    parms.cond_info = parms.cond_info(parms.conditions);
  end;
  parms.nconds = length(parms.cond_info);
  parms.event_codes = cell2mat({parms.cond_info.event_code});

  % get initial offsets from cond_info if present
  parms.num_r_offsets = length(parms.r_offset);
  parms.num_th_offsets = length(parms.th_offset);
  % independent offsets for hemi or uplow patches

  parms.offset_infix_list = {''};
  if ~parms.restrict_hemi_flag
    parms.offset_infix_list{end+1} = '_ipsi';
  end;
  if ~parms.restrict_uplow_flag
    parms.offset_infix_list{end+1} = '_cross';
  end;
  if ~parms.restrict_hemi_flag & ~parms.restrict_uplow_flag
    parms.offset_infix_list{end+1} = '_ipsi_cross';
  end;
  parms.npatches = length(parms.offset_infix_list);
  parms.init_r_offsets = zeros(parms.nareas,parms.nconds,parms.npatches);
  parms.init_th_offsets = zeros(parms.nareas,parms.nconds,parms.npatches);
  offset_types = {'r','th'};
  for a=1:parms.nareas
    for t=1:length(offset_types)
      for p=1:parms.npatches
        fieldname = sprintf('%s_offset_%s%s',...
          offset_types{t},parms.area_names{a},parms.offset_infix_list{p});
        if isfield(parms.cond_info,fieldname)
          for i=1:parms.nconds
            tmp_offset = getfield(parms.cond_info,{i},fieldname);
            if ~isempty(tmp_offset)
              parms.(['init_' offset_types{t} '_offsets'])(a,i,p) = tmp_offset;
            end;
          end;
        end;
      end;
    end;
  end;
end;

% find unique locations from cond_info
parms.same_location_conds = [1:parms.nconds]';
if ~isempty(parms.cond_info)
  parms.contrasts = cell2mat({parms.cond_info.contrast});
  parms.unique_contrasts = unique(parms.contrasts);
  parms.ncontrasts = length(parms.unique_contrasts);
  if parms.ncontrasts>1
    for i=parms.nconds:-1:2
      r = parms.cond_info(i).ecc;
      th = parms.cond_info(i).theta;
      for j=1:i-1
        tmp_r = parms.cond_info(j).ecc;
        tmp_th = parms.cond_info(j).theta;
        if tmp_r==r & tmp_th==th
          parms.same_location_conds(i) = j;
          break;
        end;
      end;
    end;
  end;
  parms.unique_location_conds = unique(parms.same_location_conds);

  if ~isfield(parms.cond_info,'loc_num')
    % set location number for each cond, save in cond_info
    for i=1:parms.nconds
      parms.cond_info(i).loc_num = ...
        find(parms.unique_location_conds==...
             parms.same_location_conds(i));
    end;
  end;

  % require each contrast level to have all locations
  for i=1:parms.ncontrasts
    ind = find(parms.contrasts == parms.unique_contrasts(i));
    if length(parms.unique_location_conds)~=...
       length(intersect(parms.same_location_conds(ind),...
                        parms.unique_location_conds))
      error('stimuli with different contrast levels must be represented at each stimulus location');
    end;
  end;
else
  parms.ncontrasts = 1;
  parms.contrasts = 1;
  parms.unique_contrasts = 1;
  parms.unique_location_conds = parms.same_location_conds;
end;

if parms.offset_niters>0
  parms.offset_group_conds = zeros(parms.nconds,1);

  if parms.offset_group_flag
    % find which quadrant each condition is in
    quad_nums = zeros(parms.nconds,1);
    for i=1:parms.nconds
      th = parms.cond_info(i).theta;
      quad_nums(i) = 1+floor(th/90);
    end;

    % find which hemisphere each condition is in
    hemi_nums = ones(parms.nconds,1);
    ind = find(ismember(quad_nums,[2,3]));
    hemi_nums(ind) = 2;

    % find which ecc ring each condition is in
    eccs = cell2mat({parms.cond_info.ecc});
%      uniq_eccs = unique(eccs(eccs>0));
    uniq_eccs = unique(eccs);
    ecc_nums = zeros(parms.nconds,1);
    for i=1:parms.nconds
      ind = find(uniq_eccs==parms.cond_info(i).ecc);
      if isempty(ind), ind = 0; end;
      ecc_nums(i) = ind;
    end;
  end;

  switch parms.offset_group_flag
    case 0
      parms.offset_group_conds = parms.same_location_conds;
    case 1 % all conditions within an ecc ring and quadrant have same offset
      for i=1:parms.nconds
        for j=i:parms.nconds
          if ecc_nums(i)==ecc_nums(j) &&...
             quad_nums(i)==quad_nums(j) &&...
             ~parms.offset_group_conds(j)
            parms.offset_group_conds(j) = i;
          end;
        end;
      end;
    case 2 % all conditions within a quadrant have same offset
      [uniq_quads,first_conds] = unique(quad_nums,'first');
      for i=1:parms.nconds
        ind = find(uniq_quads==quad_nums(i));
        parms.offset_group_conds(i) = first_conds(ind);
      end;
    case 3 % all conditions within an hemisphere have same offset
      [uniq_hemis,first_conds] = unique(hemi_nums,'first');
      for i=1:parms.nconds
        ind = find(uniq_hemis==hemi_nums(i));
        parms.offset_group_conds(i) = first_conds(ind);
      end;
    case 4 % all conditions have same offset
      parms.offset_group_conds = ones(parms.nconds,1);
    otherwise
      error('unrecognized offset_group_flag %d',parms.offset_group_flag);
  end;

  if parms.offset_mstarts==0
    parms.offset_mstarts = 1;
    parms.mstart_flag = 0;
  else
    parms.mstart_flag = 1;
    if parms.offset_niters_last > 0
      parms.offset_mstarts = parms.offset_mstarts + 1;
      % an extra iteration to search longer for a better fit in the
      %   best neighborhood
    end;
  end;
else
  parms.mstart_flag = 0;
end;

