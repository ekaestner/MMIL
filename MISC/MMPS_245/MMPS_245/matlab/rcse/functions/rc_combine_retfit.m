function rc_combine_retfit(subj,varargin)
%function rc_combine_retfit(subj,varargin)
%
% Purpose: combine results from multiple retfit directories
%   e.g. v123 + v4
%
% Required Parameters:
%   subj: FreeSurfer recon subject name
%
% Optional Parameters:
%   'outdir': output directory (full path or relative to rootdir)
%     {default = 'retfit_combo'}
%   'indirs': cell array of input directories
%     {default = {'retfit','retfit_V4'}}
%   'rootdir': root directory containing 'outdir' and 'indirs'
%     ignored if 'outdir' and 'indirs' are given as full path
%     {default = pwd}
%   'retfit_stem': retfit file stem
%     {default = 'retfit'}
%   'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%     subjdir/subj should contain the FreeSurfer subject directory
%     {default = $SUBJECTS_DIR}
%   'allow_overlap_flag': [0|1] allow visual areas in 'indirs' to overlap
%     (i.e. share vertices)
%     {default = 1}
%   'save_label_flag': [0|1] whether to save area masks as label files
%     {default = 1}
%   'resample_ico_flag': [0|1] whether to save area masks as ico4 label files
%     ignored if save_label_flag = 0
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/14/11 by Don Hagler
% Last Mod: 10/25/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin, { ...
  'outdir','retfit',[],...
  'indirs',{'retfit','retfit_V4'},[],...
  'rootdir',pwd,[],...
  'retfit_stem','retfit',[],...
  'subjdir',[],[],...
  'allow_overlap_flag',true,[false true],...
  'save_label_flag',true,[false true],...
  'resample_ico_flag',true,[false true],...
  'forceflag',false,[false true],...
... % retinotopy data / tksurfer
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'suffixlist',{'_r','_i'},{'_r','_i'},...
  'stemlist',{'pol','ecc'},[],...
  'surf','sphere',{'white','pial','inflated','sphere'},...
  'smooth',10,[0,100],...
  'fthresh',0,[],...
  'fit_fthresh',0.01,[],...
  'fmid',1.5,[],...
  'fit_fmid',0.5,[],...
  'fslope',3,[],...
  'fit_fslope',3,[],...
  'revflag',false,[false true],...
  'sph_rot',{[45 0 90],[45 -20 -90]},[],...
...
  'default_hemi','lh',[],...
  'coord_fields',{'u','v'},[],...
  'data_fields',{'pol_r','pol_i','ecc_r','ecc_i','th','r'},[],...
  'ico',4,[1:7],...
});

parms.fit_fields = cat(2,parms.coord_fields,parms.data_fields);

if isempty(parms.subjdir)
  parms.subjdir = deblank(getenv('SUBJECTS_DIR'));
  if isempty(parms.subjdir)
    error('Cannot find SUBJECTS_DIR environment variable');
  end
end;

for r=1:length(parms.indirs)
  if mmil_isrelative(parms.indirs{r})
    parms.indirs{r} = [parms.rootdir '/' parms.indirs{r}];
  end;
end;
if mmil_isrelative(parms.outdir)
  parms.outdir = [parms.rootdir '/' parms.outdir];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% copy data from first indir
indir = [parms.indirs{1} '/data'];
if ~exist(indir,'dir'), error('directory %s not found',indir); end;
outdir = [parms.outdir '/data'];
if ~exist(outdir,'dir') || parms.forceflag
  mmil_mkdir(outdir);
  cmd = sprintf('cp -p %s/* %s',indir,outdir);
  [s,r] = unix(cmd);
  if s, error('cmd %s failed:\n%s',cmd,r); end;
end;

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  [data,fit_data,area_masks,area_data] = init_results(parms);

  % create combined data, fit_data, area_masks, and area_data
  outdir = [parms.outdir '/matfiles'];
  fname_out = sprintf('%s/%s_results-%s.mat',outdir,parms.retfit_stem,hemi);
  if ~exist(fname_out,'file') || parms.forceflag
    for r=1:length(parms.indirs)
      indir = [parms.indirs{r} '/matfiles'];
      if ~exist(indir,'dir'), error('directory %s not found',indir); end;
      fname = sprintf('%s/%s_results-%s.mat',indir,parms.retfit_stem,hemi);
      if ~exist(fname,'file'), error('file %s not found',fname); end;
      results = load(fname);
      nverts = length(results.fit_data.(parms.coord_fields{1}));

      % identify new vertices with visual field coordinates
      v_curr = find(fit_data.pol_r~=0 | fit_data.pol_i~=0);                      
      v_tmp = find(results.fit_data.pol_r~=0 | results.fit_data.pol_r~=0);
      v_new = setdiff(v_tmp,v_curr); % exclude overlap

      % add data and fit_data for new vertices
      for f=1:length(parms.data_fields)
        tag = parms.data_fields{f};
        if isempty(data.(tag))
          data.(tag) = sparse(nverts,1);
        end;
        data.(tag)(v_new) = results.data.(tag)(v_new);
        if isempty(fit_data.(tag))
          fit_data.(tag) = sparse(nverts,1);
        end;
        fit_data.(tag)(v_new) = results.fit_data.(tag)(v_new);
      end;

      % identify new vertices with grid coordinates
      v_curr = find(fit_data.u | fit_data.v);
      v_tmp = find(results.fit_data.u | results.fit_data.v);
      v_new = setdiff(v_tmp,v_curr); % exclude overlap

      % add coordinates for new vertices
      for f=1:length(parms.coord_fields)
        tag = parms.coord_fields{f};
        if isempty(fit_data.(tag))
          fit_data.(tag) = sparse(nverts,1);
        end;
        fit_data.(tag)(v_new) = results.fit_data.(tag)(v_new);
      end;

      % combine area_masks and area_data
      nareas = length(area_masks);
      v_curr = [];
      if ~parms.allow_overlap_flag
        for i=1:nareas
          v_curr = [v_curr; area_masks(i).vertices];
        end;
      end;
      for i=1:length(results.area_masks);
        nareas = nareas + 1;
        % combine area_masks
        v_area = setdiff(results.area_masks(i).vertices,v_curr);
        area_masks(nareas).name = results.area_masks(i).name;
        area_masks(nareas).vertices = v_area;
        area_masks(nareas).vertices_smoothed =...
          results.area_masks(i).vertices_smoothed;
        % combine area_data
        v_all = find(results.fit_data.u | results.fit_data.v);
        for f=1:length(parms.coord_fields)
          area_data(nareas).(parms.coord_fields{f}) = sparse(nverts,1);
          area_data(nareas).(parms.coord_fields{f})(v_all) =...
            results.fit_data.(parms.coord_fields{f})(v_all);
        end;
        for f=1:length(parms.data_fields)
          area_data(nareas).(parms.data_fields{f}) = sparse(nverts,1);
          area_data(nareas).(parms.data_fields{f})(v_area) =...
            results.fit_data.(parms.data_fields{f})(v_area);
        end;
      end;
    end;

    % save results as matfile
    mmil_mkdir(outdir);
    save(fname_out,'data','fit_data','area_masks','area_data');
  else
    load(fname_out);
  end;

  % save fit_data as mgh files
  outdir = [parms.outdir '/fit'];
  mmil_mkdir(outdir);
  fname = sprintf('%s/%s_pol_r-%s.mgh',outdir,parms.retfit_stem,hemi);
  if ~exist(fname,'file') || parms.forceflag
    fs_save_mgh(full(fit_data.pol_r),fname);
  end;
  fname = sprintf('%s/%s_pol_i-%s.mgh',outdir,parms.retfit_stem,hemi);
  if ~exist(fname,'file') || parms.forceflag
    fs_save_mgh(full(fit_data.pol_i),fname);
  end;
  fname = sprintf('%s/%s_ecc_r-%s.mgh',outdir,parms.retfit_stem,hemi);
  if ~exist(fname,'file') || parms.forceflag
    fs_save_mgh(full(fit_data.ecc_r),fname);
  end;
  fname = sprintf('%s/%s_ecc_i-%s.mgh',outdir,parms.retfit_stem,hemi);
  if ~exist(fname,'file') || parms.forceflag
    fs_save_mgh(full(fit_data.ecc_i),fname);
  end;

  % save labels for each area
  outdir = [parms.outdir '/label'];
  mmil_mkdir(outdir);
  for a=1:length(area_masks)
    fname = sprintf('%s/%s.%s_%s.label',...
      outdir,hemi,parms.retfit_stem,area_masks(a).name);
    if ~exist(fname,'file') || parms.forceflag
      if ~isempty(area_masks(a).vertices_smoothed)
        verts = area_masks(a).vertices_smoothed;
      else
        verts = area_masks(a).vertices;
      end;
      fs_write_label(verts,fname,subj);
    end;
  end;
  if parms.resample_ico_flag
    setenv('SUBJECTS_DIR',parms.subjdir);
    for a=1:length(area_masks)
      fname_in = sprintf('%s/%s.%s_%s.label',...
        outdir,hemi,parms.retfit_stem,area_masks(a).name);
      fname_out = sprintf('%s/%s.%s_%s.ico%d.label',...
        outdir,hemi,parms.retfit_stem,area_masks(a).name,parms.ico);
      if ~exist(fname_in,'file'), warning('%s not found',fname_in); end;
      if ~exist(fname_out,'file') || parms.forceflag
        cmd = ['mri_label2label --srclabel ' fname_in ...
               ' --trglabel ' fname_out ...
               ' --regmethod surface'...
               ' --hemi ' hemi ...
               ' --srcsubject ' subj ...
               ' --trgsubject ico'...
               ' --trgicoorder ' num2str(parms.ico)];
        [status,result] = unix(cmd);
        if status
          fprintf('%s: WARNING: cmd %s failed:\n%s\n',cmd,result);
        end;
      end
    end;
  end;

  %% todo: create combined label with all vertices?
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create csh scripts for viewing retinoptopy data and fit

for i=1:length(parms.stemlist)
  fstem = parms.stemlist{i};
  indir = [parms.outdir '/data'];
  fname_view = sprintf('%s/view_%s.csh',parms.outdir,fstem);
  write_view_script(parms,fname_view,subj,fstem,indir)

  fstem = ['retfit_' parms.stemlist{i}];
  indir = [parms.outdir '/fit'];
  fname_view = sprintf('%s/view_%s.csh',parms.outdir,fstem);
  tmp_parms = parms;
  tmp_parms.fthresh = parms.fit_fthresh;
  tmp_parms.fmid = parms.fit_fmid;
  tmp_parms.fslope = parms.fit_fslope;
  write_view_script(tmp_parms,fname_view,subj,fstem,indir)
end;

return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,fit_data,area_masks,area_data] = init_results(parms)
  data = [];
  for f=1:length(parms.data_fields)
    data.(parms.data_fields{f}) = [];
  end;
  fit_data = [];
  for f=1:length(parms.fit_fields)
    fit_data.(parms.fit_fields{f}) = [];
  end;
  area_masks = [];
  area_data = [];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_view_script(parms,fname_out,subj,fstem,indir)
  tags = {'fname_out','fstem','indir','roi_name','label_dir','subjdir',...
    'forceflag','surf','smooth','fthresh','fmid','fslope','revflag',...
    'sph_rot','default_hemi'};

  parms.fstem = fstem;
  parms.fname_out = fname_out;
  parms.indir = indir;
  parms.label_dir = parms.outdir;
  
  args = mmil_parms2args(parms,tags);
  rc_write_viewscript(subj,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

