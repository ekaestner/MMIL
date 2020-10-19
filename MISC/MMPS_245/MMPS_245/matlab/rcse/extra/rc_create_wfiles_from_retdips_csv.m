function rc_create_wfiles_from_retdips_csv(subjname,fname_retdips,varargin)
%function rc_create_wfiles_from_retdips_csv(subjname,fname_retdips,[options])
%
% Purpose: get vertices from csv file, smooth, save as w files
%
% Usage:
%  rc_create_wfiles_from_retdips_csv(subjname,fname_retdips,'key1', value1,...);
%
% Required Input:
%  subjname - freesurfer recon subject name
%  fname_retdips - full path of csv file containing dipole vertex numbers
%
% Optional parameters:
%  'subjdir' - root directory containing freesurfer recons
%    {default: $SUBJECTS_DIR}
%  'surfname' - surface file to load (for smoothing)
%    {default: 'white'}
%  'hemilist' - cell array containing 'lh' and/or 'rh'
%    {default: {'lh' 'rh'}}
%  'outdir' - output directory for w files
%    {default: './nbrhoods'}
%  'matdir' - output directory for mat file
%    {default: './matfiles'}
%  'smooth_fwhm' - smoothing kernel full width half max (mm)
%    {default: 2}
%  'maskflag' - [0|1] whether to use area masks
%    {default: 0}
%  'maskdir' - input directory containing area masks (w files)
%    {default: './areamasks'}
%  'vert_colname' - name of column containing vertex numbers
%    {default: 'v'}
%  'forceflag' - [0|1] whether to overwrite existing output
%    {default: 0}
%
% Created:  12/10/08  by Don Hagler
% Last Mod: 02/19/11  by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'subjdir',[],[],...
  'surfname','white',[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'outdir','./nbrhoods',[],...
  'matdir','./matfiles',[],...
  'maskdir','./areamasks',[],...
  'smooth_fwhm',2,[0,100],...
  'maskflag',false,[false true],...
  'vert_colname','v',[],...
  'forceflag',false,[false true],...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load dips file
raw_dipinfo = mmil_readtext(fname_retdips);
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
area_names = unique({dipinfo.area});
thetas = unique(cell2mat({dipinfo.theta}));
eccs = unique(cell2mat({dipinfo.ecc}));
areas = [];
for a=1:length(area_names)
  areas(a).name = area_names{a};
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

[success,msg] = mkdir(parms.outdir);

% load surface, data, masks
areamasks=[];
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};

  matfile = sprintf('%s/surf-%s.mat',parms.matdir,hemi);
  if ~exist(matfile,'file') || parms.forceflag
    % load surface
    surf = fs_load_subj(subjname,hemi,parms.surfname,[],parms.subjdir);
    save(matfile,'surf');
  else
    load(matfile);
  end;

  if parms.maskflag
    matfile = sprintf('%s/areamasks-%s.mat',parms.matdir,hemi);
    if ~exist(matfile,'file') || parms.forceflag
      % load masks
      for a=1:length(areas)
        areaname = areas(a).name;
        fname = sprintf('%s/%s_mask-%s.w',parms.maskdir,areaname,hemi);
        [w,mask] = fs_read_wfile(fname);
        if isempty(w)
          return;
        end;
        mask = mask(find(w));
        areamasks(a).hemi(h).mask = mask;
      end;
      save(matfile,'areamasks');
    else
      load(matfile);
    end;
  end;
end;

dip_areas = {dipinfo.area};

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  matfile = sprintf('%s/surf-%s.mat',parms.matdir,hemi);
  load(matfile);
  if parms.maskflag
    matfile = sprintf('%s/areamasks-%s.mat',parms.matdir,hemi);
    load(matfile);
  end;
  for a=1:length(areas)
    areaname = areas(a).name;
    fprintf('%s: creating weights files for area %s...\n',mfilename,areaname);
    if parms.maskflag
      areamask = areamasks(a).hemi(h).mask;
    else
      areamask = [];
    end;
    verts = cell2mat({dipinfo.v})+1; % adjust for zero-indexing
    hemis = {dipinfo.hemisphere};
    [uniq_verts,i_uniq] = unique(verts);
    i_areas = find(strcmp(areaname,dip_areas));
    i_hemi = find(strcmp(hemi,hemis));
    i_match = intersect(i_areas,i_uniq);
    i_match = intersect(i_match,i_hemi);
    verts = verts(i_match);
    weights = [];
    for i=1:length(verts)
      ind = i_match(i);
      theta = dipinfo(ind).theta;
      ecc = dipinfo(ind).ecc;
      vert = dipinfo(ind).v;
      fname = sprintf('%s/%s-theta%0.1f-ecc%0.1f-v%d-smooth%0.1fmm-%s.w',...
        parms.outdir,areaname,theta,ecc,vert,parms.smooth_fwhm,hemi);
      if ~exist(fname,'file') || parms.forceflag
        if isempty(weights)
          weights = rc_create_smooth_area_weights(surf,hemi,verts,...
            areamask,parms.smooth_fwhm);
        end;
        if isempty(weights)
          error('failed to create smooth area weights');
        end;
        if parms.maskflag
          w = squeeze(weights(i,areamask));
          v = areamask;
        else
          w = squeeze(weights(i,:));
          v = find(w);
          w = w(v);
        end;
        fs_write_wfile(fname,w,v);
      end;
    end;
  end;
end;

