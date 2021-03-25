function rc_set_cond_offsets(fname_conds,fname_offsets,fname_out,forceflag)
%function rc_set_cond_offsets(fname_conds,fname_offsets,fname_out,forceflag)
%
% Created:  02/05/09 by Don Hagler
% Last Mod: 02/19/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag=0; end;

area_names = {'v1','v2','v3'};
cond_info = read_cond_info(fname_conds);
eccs = cell2mat({cond_info.ecc});
eccs = eccs(eccs>0);
uniq_eccs = unique(eccs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_conds,'file'), error('file %s not found',fname_conds); end;
if ~exist(fname_offsets,'file'), error('file %s not found',fname_offsets); end;
if exist(fname_out,'file') & ~forceflag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load offsets, convert to struct array

raw_offsets = mmil_readtext(fname_offsets);
[nrows,ncols]=size(raw_offsets);
collabels = {raw_offsets{1,:}};
offsets = [];
for i=2:nrows
  tmp = [];
  for j=1:ncols
    fieldname = regexprep(collabels{j},' ','_');
    tmp = setfield(tmp,fieldname,raw_offsets{i,j});
  end;
  if isempty(offsets)
    offsets = tmp;
  else
    offsets(i-1) = tmp;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add offsets to cond_info

hemi_inds = cell2mat({offsets.hemi});
uplow_inds = cell2mat({offsets.uplow});
ecc_inds = cell2mat({offsets.ecc});
areas = cell2mat({offsets.area});

if all(hemi_inds==0), hemi_inds = []; end;
if all(uplow_inds==0), uplow_inds = []; end;
if all(ecc_inds==0), ecc_inds = []; end;
if all(areas==0), areas = []; end;

for i=1:length(cond_info)
  if cond_info(i).contrast==0
    cond_info(i).hemifield = 'null';
    cond_info(i).uplowfield = 'null';
    for a=1:length(area_names)
      fieldname = sprintf('r_offset_%s',area_names{a});
      cond_info = setfield(cond_info,{i},fieldname,{1},0);
      fieldname = sprintf('th_offset_%s',area_names{a});
      cond_info = setfield(cond_info,{i},fieldname,{1},0);
    end;
    continue;
  end;

  th = cond_info(i).theta;
  if th<90 || th>270
    hemifield = 'right';
    hemi_ind = 1;
  else
    hemifield = 'left';
    hemi_ind = 2;
  end;
  cond_info(i).hemifield = hemifield;

  if th<180
    uplowfield = 'upper';
    uplow_ind = 1;
  else
    uplowfield = 'lower';
    uplow_ind = 2;
  end;
  cond_info(i).uplowfield = uplowfield;

  ecc = cond_info(i).ecc;
  ecc_ind = find(ecc==uniq_eccs);

  if isempty(hemi_inds)
    ind_hemi = 1:length(offsets);
  else
    ind_hemi = find(hemi_ind==hemi_inds);
  end;
  if isempty(uplow_inds)
    ind_uplow = 1:length(offsets);
  else
    ind_uplow = find(uplow_ind==uplow_inds);
  end; 
  if isempty(ecc_inds)
    ind_ecc = 1:length(offsets);
  else
    ind_ecc = find(ecc_ind==ecc_inds);
  end; 
  ind = intersect(ind_hemi,ind_uplow);
  ind = intersect(ind,ind_ecc);

  if isempty(areas) && length(ind)==1
    k = ind;
    for a=1:length(area_names)
      fieldname = sprintf('r_offset_%s',area_names{a});
      cond_info = setfield(cond_info,{i},fieldname,{1},offsets(k).r_offset);
      fieldname = sprintf('th_offset_%s',area_names{a});
      cond_info = setfield(cond_info,{i},fieldname,{1},offsets(k).th_offset);
    end;
  else
    for j=1:length(ind)
      k = ind(j);
      a = areas(k);
      fieldname = sprintf('r_offset_%s',area_names{a});
      cond_info = setfield(cond_info,{i},fieldname,{1},offsets(k).r_offset);
      fieldname = sprintf('th_offset_%s',area_names{a});
      cond_info = setfield(cond_info,{i},fieldname,{1},offsets(k).th_offset);
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rc_write_cond_info(cond_info,fname_out,forceflag);


