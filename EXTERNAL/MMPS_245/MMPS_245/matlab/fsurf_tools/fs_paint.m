function fnames_out = fs_paint(subj,fname_in,varargin);
%function fnames_out = fs_paint(subj,fname_in,[options]);
%
% Purpose:
%   a wrapper around freesurfer binaries for sampling
%   volume data to surface, sampling to sphere, and smoothing
%
% Usage:
%  fs_paint(subj,fname_in,'key1', value1,...);
%
% Required parameters:
%  subj is a string specifying the subject name
%  fname_in is a string specifying the full or relative path
%    of the functional volume (must be mgh format)
%    (can be empty if 'meas' is supplied -- see below)
%
% Optional parameters:
%  'hemi': should be either 'lh' for left hemisphere or 'rh' for right hemi
%    {default = both}
%  'meas': surface measure such as 'thickness', 'sulc', 'area', or 'icoarea'
%     Note: for 'icoarea', must be an "ico" subject
%    {default = []}
%  'outstem' : output file stem (omit extension, hemi)
%    {default = fname_in or meas}
%  'outdir': output directory
%    if empty, will use pwd or subjdir/analysis if meas supplied
%    {default = []}
%  'infix' : add extra suffix to outstem before hemi and extension
%    {default = []}
%  'outtype': output file type ('w', 'mgh', or 'mgz')
%    {default = 'mgh'}
%  'intype': input file type (e.g. 'mgh', 'analyze', 'bfloat')
%    {default = 'mgh'}
%  'regfile': register.dat file containing 4x4 registration matrix
%    {default = []}
%  'regmat': 4x4 registration matrix (ignored if 'regfile' is specified)
%    {default = identity matrix}
%  'inplane': inplane resolution (mm) (ignored if 'regfile' is specified)
%    {default = 1}
%  'slicethick': slice thickness (mm) (ignored if 'regfile' is specified)
%    {default = 1}
%  'sphere_flag': [0|1] whether to sample to icosohedral sphere
%    {default = 0}
%  'sphere_ico': spherical icosahedral order (0-7)
%       Order  Number of Vertices
%         0              12
%         1              42
%         2             162
%         3             642
%         4            2562
%         5           10242
%         6           40962
%         7          163842
%    {default = 7}
%  'projfrac_flag': [0|1] whether to use projdist (0) or projfract (1)
%    {default = 0}
%  'projdist': distance (mm) to project along surface vertex normal
%    {default = 1}
%  'projfrac': fractional distance to project along surface vertex normal
%    relative to cortical thickness
%    {default = 0.5}
%  'projfrac_avg': vector of [min max del] for averaging multiple samples
%     with mri_surf2surf projfract-avg option
%    If empty, use projfrac instead if projfrac_flag=1
%    {default = []}
%  'projdist_avg': vector of [min max del] for averaging multiple samples
%     with mri_surf2surf projdist-avg option
%    If empty, use projdist instead
%    {default = []}
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 0}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default = 0}
%  'cortex_flag': apply smoothing only inside cortex mask
%    only applies if $FREESURFER_VER >= 400
%    {default = 1}
%  'mask_midbrain_flag': [0|1] mask out thalamus, corpus collosum
%    {default = 0}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'surfname': surface to paint onto
%    {default = white}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Output:
%   fnames_out: cell array of output file names (e.g. left and right)
%     with multiple steps (e.g. sphere, smoothing, mask), only the final
%     output files are included in fnames_out
%
% Created:  10/19/06 by Don Hagler
% Last Mod: 08/21/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

fnames_out = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'hemi',[],{'lh','rh'},...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'outstem',[],[],...
  'outdir',[],[],...
  'outtype','mgh',{'w','mgh','mgz'},...
  'regmat',eye(4),[],...
  'regfile',[],[],...
  'inplane',[],[],...
  'slicethick',[],[],...
  'sphere_flag',false,[false true],...
  'sphere_ico',7,[0 7],...
  'subjdir',[],[],...
  'surfname','white',[],...
  'projdist',1,[-10,10],...
  'projfrac',0.5,[-2,2],...
  'projfrac_flag',false,[false true],...
  'projdist_avg',[],[],...
  'projfrac_avg',[],[],...
  'intype','mgh',{'mgh','bfloat' 'analyze'},...
  'forceflag',false,[false true],...
  'overwrite_flag',false,[false true],...
  'infix',[],[],...
  'smoothsteps',[],[],...
  'sphsmoothsteps',[],[],...
  'cortex_flag',true,[false true],...
  'mask_midbrain_flag',false,[false true],...
  'meas',[],{'thickness','sulc','area','icoarea'},...
...
  'verbose',false,[false true],...
  'mask_roinums',[1,5],[],... % 'unknown' and 'corpuscallosum'
});

if parms.overwrite_flag, parms.forceflag = true; end;
if ~isempty(parms.hemi), parms.hemilist = {parms.hemi}; end;

if isempty(parms.subjdir)
  parms.subjdir = getenv('SUBJECTS_DIR');
  if isempty(parms.subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  surffile = sprintf('%s/%s/surf/%s.%s',parms.subjdir,subj,hemi,parms.surfname);
  if ~exist(surffile,'file')
    error('surface file %s not found',surffile);
  end
end;
setenv('SUBJECTS_DIR',parms.subjdir);

if ~isempty(parms.regfile) & ~exist(parms.regfile,'file')
  error('regfile %s not found',parms.regfile);
elseif any(size(parms.regmat)~=[4 4])
  error('regmat is wrong size');
end;

if parms.cortex_flag
  fs_ver = getenv('FREESURFER_VER');
  if isempty(fs_ver)
    fprintf('%s: WARNING: FREESURFER_VER not set, setting cortex_flag=0\n',...
      mfilename);
    parms.cortex_flag = 0;
  else
    fs_ver = str2num(fs_ver);
    if isempty(fs_ver)
      fprintf('%s: WARNING: non-numeric FREESURFER_VER, setting cortex_flag=0\n',...
        mfilename);
      parms.cortex_flag = 0;
    elseif fs_ver < 400
      parms.cortex_flag = 0;
    end;
  end;
end;

if isempty(parms.outstem)
  if isempty(parms.meas)
    [tmp_path,tmp_fstem,tmp_ext] = fileparts(fname_in);
    if isempty(parms.outdir)
      parms.outstem = [tmp_path '/' tmp_fstem];
    else
      parms.outstem = [parms.outdir '/' tmp_fstem];
    end;
  else
    if isempty(parms.outdir)
      parms.outstem = sprintf('%s/%s/analysis/%s',...
        parms.subjdir,subj,parms.meas);
    else
      parms.outstem = [parms.outdir '/' parms.meas];
    end;
  end;
elseif mmil_isrelative(parms.outstem)
  if ~isempty(parms.outdir)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
end;
[outdir,tmp_fstem] = fileparts(parms.outstem);
if isempty(outdir)
  outdir = pwd;
  parms.outstem = [outdir '/' parms.outstem];
end;
mmil_mkdir(outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.meas) % i.e. painting from volume to surface
  if ~isempty(parms.regfile)
   [parms.regmat,tmp_subj,parms.inplane,parms.slicethick] = fs_read_regdat(parms.regfile);
  end;

  % create tmp register.dat file (in case parms.regfile not specified or subj is different)
  [tmp,tmp_stem] = fileparts(tempname);
  tmp_regfile = sprintf('%s/%s-register.dat',outdir,tmp_stem);

  fs_write_regdat(tmp_regfile,...
    'M',parms.regmat,...
    'subj',subj,...
    'inplane',parms.inplane,...
    'slicethick',parms.slicethick,...
    'forceflag',1);
end;

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};

  % sample from volume to surface
  if isempty(parms.infix)
    prefix_out = parms.outstem;
  else
    prefix_out = sprintf('%s%s',parms.outstem,parms.infix);
  end;
  outfile = sprintf('%s-%s.%s',prefix_out,hemi,parms.outtype);
  fnames_out{h} = outfile;
  if ~exist(outfile,'file') | parms.forceflag
    % create mri_surf2surf or mri_vol2surf cmd string
    if ~isempty(parms.meas)
      if parms.verbose
        fprintf('%s: painting cortical %s\n',mfilename,parms.meas);
      end;
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_surf2surf --srcsubject %s --hemi %s',...
        parms.subjdir,subj,hemi);
      switch parms.meas
        case {'thickness','sulc'}
          cmd = sprintf('%s --srcsurfval %s --src_type curv',cmd,parms.meas);
        case 'area'
          cmd = sprintf('%s --srcsurfval area --jac --src_type curv',cmd);
        case 'icoarea'
          cmd = sprintf('%s --sval %s/%s/surf/%s.white.avg.area.mgh',...
            cmd,parms.subjdir,subj,hemi);
      end;
      cmd = sprintf('%s --tval %s --trgsubject %s',cmd,outfile,subj);
    else
      if parms.verbose
        fprintf('%s: painting func volume from %s for hemi %s\n',...
          mfilename,fname_in,hemi);
      end;
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_vol2surf --src %s --src_type %s --srcreg %s --hemi %s',...
        parms.subjdir,fname_in,parms.intype,tmp_regfile,hemi);
      cmd = sprintf('%s --surf %s --out %s --out_type %s --trgsubject %s',...
        cmd,parms.surfname,outfile,parms.outtype,subj);
      if ~parms.projfrac_flag
        if ~isempty(parms.projdist_avg)
          cmd = sprintf('%s --projdist-avg %s',...
            cmd,sprintf('%0.3f ',parms.projdist_avg));
        else
          cmd = sprintf('%s --projdist %0.3f',...
            cmd,parms.projdist);
        end;
      else
        if ~isempty(parms.projfrac_avg)
          cmd = sprintf('%s --projfrac-avg %s',...
            cmd,sprintf('%0.3f ',parms.projfrac_avg));
        else
          cmd = sprintf('%s --projfrac %0.3f',...
            cmd,parms.projfrac);
        end;
      end;
    end;
    cmd = sprintf('%s --noreshape',cmd);
    % run cmd
    [status,result]=unix(cmd);
    if status
      error('cmd %s failed:\n%s',cmd,result);
    end;
  end;

  % mask out midbrain
  if parms.mask_midbrain_flag
    prefix_in = prefix_out;
    prefix_out = sprintf('%s-mbmask',prefix_in);
    infile = sprintf('%s-%s.%s',prefix_in,hemi,parms.outtype);
    outfile = sprintf('%s-%s.%s',prefix_out,hemi,parms.outtype);
    fnames_out{h} = outfile;
    %% todo: if parms.outtype is w, then need to load file and mask surfstats
    fs_mask_surfmgh_with_aparc(subj,hemi,infile,outfile,...
      parms.subjdir,parms.mask_roinums);
  end;

  % smooth on surface
  if parms.smoothsteps
    prefix_in = prefix_out;
    prefix_out = sprintf('%s-sm%d',prefix_in,parms.smoothsteps);
    infile = sprintf('%s-%s.%s',prefix_in,hemi,parms.outtype);
    outfile = sprintf('%s-%s.%s',prefix_out,hemi,parms.outtype);
    fnames_out{h} = fs_surf2surf(infile,subj,...
      'fname_out',outfile,...
      'hemi',hemi,...
      'outtype',parms.outtype,...
      'smooth_out',parms.smoothsteps,...
      'cortex_flag',parms.cortex_flag',...
      'subjdir',parms.subjdir,...
      'verbose',parms.verbose,...
      'forceflag',parms.forceflag);
  end;

  % sample to sphere
  if parms.sphere_flag
    prefix_in = prefix_out;
    infile = sprintf('%s-%s.%s',prefix_in,hemi,parms.outtype);
    if parms.sphere_ico == 7
      sph_infix = 'sphere';
    else
      sph_infix = sprintf('ico%d',parms.sphere_ico);
    end;
    prefix_out = sprintf('%s-%s',prefix_in,sph_infix);

    if ~parms.cortex_flag || parms.sphere_ico~=7
      %% note: cortex_flag does not work with sphmoothsteps if trgsubj is ico
      %%    because trgsubj must have cortex label
      %%    similarly will not work if ico level is not 7
      if parms.sphsmoothsteps
        prefix_out = sprintf('%s-sm%d',prefix_out,parms.sphsmoothsteps);
      end;    
      outfile = sprintf('%s-%s.%s',prefix_out,hemi,parms.outtype);
      fnames_out{h} = fs_surf2surf(infile,subj,...
        'fname_out',outfile,...
        'hemi',hemi,...
        'trgsubj','ico',...
        'icolevel',parms.sphere_ico,...
        'outtype',parms.outtype,...
        'smooth_out',parms.sphsmoothsteps,...
        'subjdir',parms.subjdir,...
...        'verbose',parms.verbose,...
        'forceflag',parms.forceflag);
    else
      infile = sprintf('%s-%s.%s',prefix_in,hemi,parms.outtype);
      outfile = sprintf('%s-%s.%s',prefix_out,hemi,parms.outtype);
      fnames_out{h} = fs_surf2surf(infile,subj,...
        'fname_out',outfile,...
        'hemi',hemi,...
        'trgsubj','ico',...
        'icolevel',parms.sphere_ico,...
        'outtype',parms.outtype,...
        'subjdir',parms.subjdir,...
...        'verbose',parms.verbose,...
        'forceflag',parms.forceflag);
    end;

    if parms.sphsmoothsteps && parms.cortex_flag
      subjdir = [getenv('FREESURFER_HOME') '/subjects'];
      prefix_out = sprintf('%s-sm%d',prefix_out,parms.sphsmoothsteps);
      infile = outfile;
      outfile = sprintf('%s-%s.%s',prefix_out,hemi,parms.outtype);
      fnames_out{h} = fs_surf2surf(infile,'fsaverage',...
        'fname_out',outfile,...
        'hemi',hemi,...
        'outtype',parms.outtype,...
        'smooth_out',parms.sphsmoothsteps,...
        'cortex_flag',parms.cortex_flag',...
        'subjdir',subjdir,...
...        'verbose',parms.verbose,...
        'forceflag',parms.forceflag);
    end;
  end;
end;

if exist('tmp_regfile','var') & exist(tmp_regfile,'file')
  delete(tmp_regfile);
end;
