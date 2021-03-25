function cond_info = rc_RCSE_set_cond_info_offsets(parms,r_offsets,th_offsets)
%function cond_info = rc_RCSE_set_cond_info_offsets(parms,r_offsets,th_offsets)
%
% Purpose: set cond r and th offsets in parms.cond_info struct
%
% Required Input:
%   parms: RCSE parms struct
%   r_offsets: eccentricity offset vector
%   th_offsets: polar angle offset vector
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/04/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

cond_info = parms.cond_info;
r_offsets = rc_check_bounds(r_offsets,parms.r_offset_range);
th_offsets = rc_check_bounds(th_offsets,parms.th_offset_range);
offset_types = {'r','th'};
for a=1:parms.nareas
  if ~ismember(a,parms.use_areas), continue; end;
  for i=1:parms.nconds
    for t=1:length(offset_types)
      for p=1:parms.npatches
        switch offset_types{t}
          case 'r'
            tmp_offset = r_offsets(a,i,p);
          case 'th'
            tmp_offset = th_offsets(a,i,p);
        end;
        fieldname = sprintf('%s_offset_%s%s',...
          offset_types{t},parms.area_names{a},parms.offset_infix_list{p});
        cond_info = setfield(cond_info,{i},fieldname,{1},tmp_offset);
      end;
    end;
  end;
  fprintf('%s: %s r_offsets:\n  ',mfilename,parms.area_names{a});
  fprintf('%0.4f ',squeeze(r_offsets(a,:,:)));
  fprintf('\n');
  fprintf('%s: %s th_offsets:\n  ',mfilename,parms.area_names{a});
  fprintf('%0.4f ',squeeze(th_offsets(a,:,:)));
  fprintf('\n');
end;
