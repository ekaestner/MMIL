function retmap=rc_define_retmap_from_csv(subjname,fname_dips,fname_conds,varargin)
%function retmap=rc_define_retmap_from_csv(subjname,fname_dips,fname_conds,[options])
%
% Purpose: create retmap struct using dipole and condition info from csv files
%
% Usage:
%  retmap=rc_define_retmap_from_csv(subjname,fname_dips,fname_conds,'key1', value1,...);
%
% Required Input:
%  subjname - freesurfer recon subject name
%  fname_dips - full path of csv file containing dipole vertex numbers
%  fname_conds - full path of csv file containing condition information
%
% Optional parameters:
%  'subjdir' - root directory containing freesurfer recons
%    {default: $SUBJECTS_DIR}
%  'vert_colname' - name of column containing vertex numbers
%    {default: 'v'}
%  'nbrhood_flag' - [0|1] whether to use wfiles in nbrhood_dir
%    {default: 0}
%  'nbrhood_dir' - input directory containing area masks (w files)
%    {default: './nbrhoods'}
%  'nbrhood_smooth' - fwhm (mm) smoothing kernel used to create nbrhood w files
%    {default: 2}
%  'contrasts' - vector of acceptable contrast levels
%    if empty, use all contrasts specified in fname_conds
%    {default: []}
%
% Created:  12/10/08  by Don Hagler
% Last Mod: 02/19/11  by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

retmap = [];

if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin, { ...
  'subjdir',[],[],...
  'vert_colname','v',[],...
  'nbrhood_flag',false,[false true],...
  'nbrhood_dir','./nbrhoods',[],...
  'nbrhood_smooth',2,[0,1000],...
  'contrasts',[],[],...
});

if isempty(parms.subjdir)
  parms.subjdir = getenv('SUBJECTS_DIR');
end;
if isempty(parms.subjdir)
  error('no subjdir specified',parms.subjdir);
end;
if ~exist(parms.subjdir,'dir')
  error('subjdir %s not found',parms.subjdir);
end;

if parms.nbrhood_flag && ~exist(parms.nbrhood_dir,'dir')
  error('nbrhood_dir %s not found',parms.nbrhood_dir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load dips file
raw_dipinfo = mmil_readtext(fname_dips);
% convert dip info to struct array
[nrows,ncols]=size(raw_dipinfo);
collabels = {raw_dipinfo{1,:}};
dipinfo = [];
for i=2:nrows
  tmp = [];
  for j=1:ncols
    fieldname = regexprep(collabels{j},' ','_');
    tmp = setfield(tmp,fieldname,raw_dipinfo{i,j});
  end;
  if isempty(dipinfo)
    dipinfo = tmp;
  else
    dipinfo(i-1) = tmp;
  end;
end;

% replace dipinfo.v with dipinfo.v1, dipinfo.v2, or dipinfo.v3
if ~strcmp(parms.vert_colname,'v')
  if ~isfield(dipinfo,parms.vert_colname)
    error('vert_colname is invalid -- no column with label "%s"',...
      parms.vert_colname);
  end;
  for i=1:length(dipinfo)
    dipinfo(i).v = getfield(dipinfo,{i},parms.vert_colname);
  end;
elseif ~isfield(dipinfo,'v')
  error('missing column with label "v"');
end;

% load cond file
raw_cond_info = mmil_readtext(fname_conds);
% convert cond info to struct array
[nrows,ncols]=size(raw_cond_info);
collabels = {raw_cond_info{1,:}};
cond_info = [];
for i=2:nrows
  tmp = [];
  for j=1:ncols
    fieldname = regexprep(collabels{j},' ','_');
    tmp = setfield(tmp,fieldname,raw_cond_info{i,j});
  end;
  if isempty(cond_info)
    cond_info = tmp;
  else
    cond_info(i-1) = tmp;
  end;
end;

if ~isempty(parms.contrasts)
  contrasts = cell2mat({cond_info.contrast});
  cond_order = ismember(contrasts,parms.contrasts);
  cond_info = cond_info(cond_order);
else
  cond_order = 1:length(cond_info);
end;
nconds = length(cond_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define visual areas

dip_areas = {dipinfo.area};
dip_thetas = cell2mat({dipinfo.theta});
dip_eccs = cell2mat({dipinfo.ecc});
area_names = unique(dip_areas);
areas = [];
for i=1:length(cond_info)
  i_theta = find(cond_info(i).theta == dip_thetas);
  i_ecc = find(cond_info(i).ecc == dip_eccs);
  for a=1:length(area_names)
    areas(a).name = area_names{a};
    i_areas = find(strcmp(areas(a).name,dip_areas));
    i_match = intersect(i_areas,i_theta);
    i_match = intersect(i_match,i_ecc);
    if isempty(i_match)
      error('dipole not found with area %s, theta = %0.1f, and ecc = %0.1f',...
        areas(a).name,cond_info(i).theta,cond_info(i).ecc);
    elseif length(i_match)>1
      i_match = i_match(1);
      fprintf('%s: WARNING: multiple dipoles for area %s, theta = %0.1f, and ecc = %0.1f\n',...
        mfilename,areas(a).name,cond_info(i).theta,cond_info(i).ecc);
    end;
    hemisphere = dipinfo(i_match).hemisphere;
    v = dipinfo(i_match).v;
    if parms.nbrhood_flag
      fname = sprintf('%s-theta%0.1f-ecc%0.1f-v%d-smooth%0.1fmm-%s.w',...
        areas(a).name,dipinfo(i_match).theta,dipinfo(i_match).ecc,...
        dipinfo(i_match).v,parms.nbrhood_smooth,hemisphere);
      full_fname = [parms.nbrhood_dir '/' fname];
      if ~exist(full_fname)
        error('file %s not found',full_fname);
      end;
      if strcmp(hemisphere,'rh')
        areas(a).wfiles(i).lh_fname = [];
        areas(a).wfiles(i).rh_fname = fname;
      else
        areas(a).wfiles(i).lh_fname = fname;
        areas(a).wfiles(i).rh_fname = [];
      end;
    else
      if strcmp(hemisphere,'rh')
        areas(a).verts(i).v_lh = []; % ipsilateral
        areas(a).verts(i).w_lh = [];
        areas(a).verts(i).v_rh = v; % contralateral
        areas(a).verts(i).w_rh = 1;
      else
        areas(a).verts(i).v_lh = v; % contralateral
        areas(a).verts(i).w_lh = 1;
        areas(a).verts(i).v_rh = []; % ipsilateral
        areas(a).verts(i).w_rh = [];
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

retmap.cond_info = cond_info;
retmap.cond_order = cond_order;
retmap.num_locs = nconds;
retmap.areas = areas;
retmap.ret_dips = [];
retmap.nonret_dips = [];

return;

