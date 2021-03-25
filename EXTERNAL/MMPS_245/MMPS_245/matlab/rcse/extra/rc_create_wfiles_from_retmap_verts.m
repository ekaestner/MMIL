function rc_create_wfiles_from_retmap_verts(subjname,retmap,varargin)
%function rc_create_wfiles_from_retmap_verts(subjname,retmap,[options])
%
% Purpose: get vertices from retmap struct, smooth, save as w files
%
% Usage:
%  rc_create_wfiles_from_retmap_verts(subjname,retmap,'key1', value1,...);
%
% Required Input:
%  subjname - freesurfer recon subject name
%  retmap - retmap struct
%           (see rc_define_retmap_from_csv.m or examp_define_retmap.m)
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
%  'nbrhood_flag' - [0|1] whether to generate neighborhood around a single
%     vertex (1) or simply include the specified cluster of vertices (0)
%     {default: 1}
%  'maskflag' - [0|1] whether to use area masks to constrain the neighborhoods
%    {default: 0}
%  'maskdir' - input directory containing area masks (w files)
%    {default: './areamasks'}
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
  'smooth_fwhm',2,[0,100],...
  'nbrhood_flag',true,[false true],...
  'maskflag',false,[false true],...
  'maskdir','./areamasks',[],...
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

if isfield(retmap,'cond_info')
  contrasts = cell2mat({retmap.cond_info.contrast});
else
  contrasts = [];
end;

if ~isfield(retmap,'areas_vbase')
  areas_vbase = 0;
else
  areas_vbase = retmap.areas_vbase;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[success,msg] = mkdir(parms.outdir);

if parms.nbrhood_flag
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
end;

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  if parms.nbrhood_flag
    matfile = sprintf('%s/smw_surf-%s.mat',parms.matdir,hemi);
    load(matfile);
    if parms.maskflag
      matfile = sprintf('%s/smw_areamasks-%s.mat',parms.matdir,hemi);
      load(matfile);
    end;
  end;
  for a=1:length(retmap.areas)
    areaname = retmap.areas(a).name;
    fprintf('%s: creating weights files for area %s...\n',mfilename,areaname);
    tmp_verts = retmap.areas(a).verts;
    if parms.nbrhood_flag
      if parms.maskflag
        areamask = areamasks(a).hemi(h).mask;
      else
        areamask = [];
      end;
      verts = [];
      dips = [];
      j = 1;
      for i=1:length(tmp_verts)
        if ~isempty(contrasts) && contrasts(i)==0
          continue;
        end;
        switch hemi
          case 'lh'
            tmp_vert = tmp_verts(i).v_lh;
          case 'rh'
            tmp_vert = tmp_verts(i).v_rh;
        end;
        if areas_vbase==0
          tmp_vert = tmp_vert + 1;
        end;
        if ~isempty(tmp_vert)
          verts(j) = tmp_vert(1);
          dips(j) = i;
          j = j + 1;
        end;
      end;
      weights = rc_create_smooth_area_weights(surf,hemi,verts,areamask,smooth_fwhm);
      if isempty(weights), return; end;
      for j=1:length(verts)
        fname = sprintf('%s/%s-retmap-dip%d-v%d-smooth%0.1fmm-%s.w',...
          parms.outdir,areaname,dips(j),verts(j),smooth_fwhm,hemi);
        if ~exist(fname,'file') || parms.forceflag
          if parms.maskflag
            w = squeeze(weights(j,areamask));
            v = areamask;
          else
            w = squeeze(weights(j,:));
            v = find(w);
            w = w(v);
          end;
          fs_write_wfile(fname,w,v);
        end;
      end;
    else
      for i=1:length(tmp_verts)
        if ~isempty(contrasts) && contrasts(i)==0
          continue;
        end;
        switch hemi
          case 'lh'
            v = tmp_verts(i).v_lh;
            w = tmp_verts(i).w_lh;
          case 'rh'
            v = tmp_verts(i).v_rh;
            w = tmp_verts(i).w_rh;
        end;
        if areas_vbase==0
          v = v + 1;
        end;
        if isempty(v), continue; end;
        fname = sprintf('%s/%s-retmap-loc%02d-%s.w',...
          parms.outdir,areaname,i,hemi);
        if ~exist(fname,'file') || parms.forceflag
          fs_write_wfile(fname,w,v);
        end;
      end;
    end; % if nbrhood_flag
  end;
end;

