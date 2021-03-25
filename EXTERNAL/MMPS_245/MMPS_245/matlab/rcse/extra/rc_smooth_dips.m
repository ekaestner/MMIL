function rc_smooth_dips(subj,lh_dip_name,rh_dip_name,smoothsteps)
%function rc_smooth_dips(subj,lh_dip_name,rh_dip_name,[smoothsteps])
%
% subj: freesurfer subject name (must exist in $SUBJECTS_DIR)
% lh_dip_name:
% rh_dip_name:
% smoothsteps: 
%
% Early Mod: 04/03/09 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

prefix = 'smdips';
surfname = 'white';
if ~exist('smoothsteps','var') || isempty(smoothsteps)
  smoothsteps = 0;
end;

% load surface and dipole info
matname=sprintf('matfiles/%s_surf.mat',prefix);
if exist(matname)
  load(matname);
else
  fprintf('%s: loading surfaces...\n',mfilename);
  surf_lh = fs_load_subj(subj,'lh',surfname);
  surf_rh = fs_load_subj(subj,'rh',surfname);

  fprintf('%s: finding neighbor relations...\n',mfilename);
  surf_lh = fs_find_neighbors(surf_lh);
  surf_rh = fs_find_neighbors(surf_rh);

  fprintf('%s: saving surface...\n',mfilename);
  save(matname,'surf_lh','surf_rh');
end;

% load original dipole info
matname=sprintf('matfiles/%s_orig.mat',prefix);
if exist(matname)
  load(matname);
else
  % get information about source locations
  fprintf('%s: loading dipole information...\n',mfilename);
  dip_info_lh=ts_read_dip_file(lh_dip_name);
  dip_info_rh=ts_read_dip_file(rh_dip_name);

  fprintf('%s: saving dipole info...\n',mfilename);
  save(matname,'dip_info_lh', 'dip_info_rh');
end;

orig_dip_info_lh = dip_info_lh;
orig_dip_info_rh = dip_info_rh;

% smooth dipoles
for s=1:length(smoothsteps)
  n = smoothsteps(s);
  fprintf('%s: smoothing dipoles %d steps...\n',mfilename,n);
  % smooth coordinates and normal vector
  dip_info_lh = orig_dip_info_lh;
  dip_info_rh = orig_dip_info_rh;
  for i=1:6
    dip_info_lh(i,:) = fs_smooth(surf_lh,dip_info_lh(i,:)',n)';
    dip_info_rh(i,:) = fs_smooth(surf_rh,dip_info_rh(i,:)',n)';
  end;
  fprintf('%s: saving smoothed dipoles...\n',mfilename);
  matname=sprintf('matfiles/%s_sm%02d.mat',prefix,n);
  save(matname,'dip_info_lh', 'dip_info_rh');
end

fprintf('%s: finished.\n',mfilename);
