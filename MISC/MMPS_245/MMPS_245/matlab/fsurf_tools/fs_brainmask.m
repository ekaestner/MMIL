function fname_out = fs_brainmask(fname_in,varargin)
%function fname_out = fs_brainmask(fname_in,varargin)
%
% Purpose: create brain mask volume and surfaces from T1-weighted image
%   using mri_watershed
%
% Usage: fs_brainmask(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_out': output file name
%    {default = 'brainmask.mgz'}
%  'conform_flag': [0|1] convert input image to be 256^3 and uchar type
%    {default = 1}
%  'nu_flag': [0|1] run N3 non-uniformity intensity correction
%    {default = 1}
%  'T1_flag': [0|1] whether to perform T1 intensity normalization
%    {default = 1}
%  'surfs_flag': [0|1] create mesh files
%    for inner skull, outer skull, and outer scalp
%    {default = 1}
%  'nu3T_flag': [0|1] whether to use nu intensity correction
%    with settings optimized for 3T
%    irrelevant if nu_flag = 0
%    {default = 0}
%  'tal_flag': [0|1] register to Talairach space
%    for improved non-uniformity correction
%    irrelevant if nu_flag = 0
%    {default = 1}
%  'res2raw_flag': [0|1] resample final brain mask to slicing of input image
%    {default = 1}
%  'conform_options': options used when calling mri_convert to conform input
%    {default = '--conform -rt cubic -odt uchar -ns 1'}
%  'fname_talxfm': [0|1] existing Talairach transform file
%    if supplied, will skip Talairach registration
%    {default = []}
%  'surfs_outdir': output directory for mesh files
%    irrelevant if surfs_flag = 0
%    {default = 'brainmask_surfs'}
%  'surfs_nverts': approximate number of vertices for reduced resolution meshes
%    {default = 3000}
%  'watershed_options': additional options for mri_watershed
%    {default = []}
%  'tmpdir': temporary directory containing intermediate output
%    {default = 'tmp_brainmask'}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  12/29/14 by Don Hagler
% Last Mod: 02/19/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output file
if exist(parms.fname_out) && ~parms.forceflag
  fname_out = parms.fname_out;
  return;
end;

% create temporary output directory
mmil_mkdir(parms.tmpdir);

% conform input to FreeSurfer standard
if parms.conform_flag
  parms.fname_in = conform_input(parms,'orig');
end;

% N3 nonuniformity correction
if parms.nu_flag
  parms.fname_in = nu_correction(parms,'nu');
end;

% T1 intensity normalization
if parms.T1_flag
  parms.fname_in = T1_correction(parms,'T1');
end;

% create brain mask
parms.fname_in = create_brainmask(parms,'brainmask_cor');

% resample to original slicing and resolution
if parms.res2raw_flag
  parms.fname_in = resamp_to_raw(parms,'brainmask_raw');
end;

% copy to final output
if ~exist(parms.fname_out,'file') || parms.forceflag
  fs_copy_mgh(parms.fname_in,parms.fname_out);
end;
fname_out = parms.fname_out;

% remove temporary files
if parms.cleanup_flag
  parms = cleanup_tmpdir(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_out','brainmask.mgz',[],...
    'conform_flag',true,[false true],...
    'nu_flag',true,[false true],...
    'T1_flag',true,[false true],...
    'surfs_flag',true,[false true],...
    'nu3T_flag',false,[false true],...
    'tal_flag',true,[false true],...
    'res2raw_flag',true,[false true],...
    'conform_options','--conform -rt cubic -odt uchar -ns 1',[],...
    'fname_talxfm',[],[],...
    'surfs_outdir','brainmask_surfs',[],...
    'surfs_nverts',3000,[100,10000],...
    'watershed_options',[],[],...
    'tmpdir','tmp_brainmask',[],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'nu_tags',{'fname_out','nu3T_flag','tal_flag','fname_talxfm',...
               'tmpdir','cleanup_flag','verbose','forceflag'},[],...
    'T1_tags',{'fname_out','fname_mask','verbose','forceflag'},[],...
    'surf_names',{'lh.bem_inner_skull_surface',...
                  'lh.bem_outer_skull_surface',...
                  'lh.bem_outer_skin_surface'},[],...
    'tri_names',{'innner_skull','outer_skull','outer_scalp'},[],...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  parms.fname_raw = parms.fname_in;
  if parms.tal_flag &&...
     ~isempty(parms.fname_talxfm) && ~exist(parms.fname_talxfm,'file')
    error('file %s not found',parms.fname_talxfm);
  end;
  if mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [pwd '/' parms.tmpdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = conform_input(parms,outstem)
  fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  fs_mri_convert(parms.fname_in,fname_out,...
    'options',parms.conform_options,...
    'verbose',parms.verbose,...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = nu_correction(parms,outstem)
  parms.fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  parms.cleanup_flag = 0;
  args = mmil_parms2args(parms,parms.nu_tags);
  fname_out = fs_nu_corr(parms.fname_in,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = T1_correction(parms,outstem)
  parms.fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  parms.cleanup_flag = 0;
  args = mmil_parms2args(parms,parms.T1_tags);
  fname_out = fs_T1_corr(parms.fname_in,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = create_brainmask(parms,outstem)
  fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  if parms.surfs_flag
    fname_tri1 = sprintf('%s/%s.tri',parms.surfs_outdir,parms.tri_names{1});
    fname_tri2 = sprintf('%s/%s.tri',parms.surfs_outdir,parms.tri_names{2});
    fname_tri3 = sprintf('%s/%s.tri',parms.surfs_outdir,parms.tri_names{3});
  end;
  if ~exist(fname_out,'file') ||...
      (parms.surfs_flag &&...
       (~exist(fname_tri1,'file') ||...
        ~exist(fname_tri2,'file') || ~exist(fname_tri3,'file'))) ||...
      parms.forceflag
    if parms.verbose
      fprintf('%s: creating brain mask...\n',mfilename);
      tic;
    end;
    cmd = 'mri_watershed';
    if parms.surfs_flag
      mmil_mkdir(parms.surfs_outdir);
      cmd = sprintf('cd %s; %s -surf bem -useSRAS',parms.tmpdir,cmd);
    end;
    if ~isempty(parms.watershed_options)
      cmd = sprintf('%s %s',cmd,parms.watershed_options);
    end;
    cmd = sprintf('%s %s %s',cmd,parms.fname_in,fname_out);
    [status,msg] = unix(cmd);
    if status || ~exist(fname_out,'file')
      error('brain mask creation failed:\n%s\n%s\n',cmd,msg);
    end;
    if parms.verbose
      toc;
    end;
    % convert surfaces to tri format
    if parms.surfs_flag
      if parms.verbose
        fprintf('%s: converting surface files to tri format...\n',mfilename);
        tic;
      end;
      for i=1:length(parms.surf_names)
        fname_surf = sprintf('%s/%s',parms.tmpdir,parms.surf_names{i});
        fname_tri = sprintf('%s/%s.tri',parms.surfs_outdir,parms.tri_names{i});
        fs_convert_fsurf2tri(fname_surf,fname_tri,parms.surfs_nverts);
      end;
      if parms.verbose
        toc;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = resamp_to_raw(parms,outstem)
  fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  fs_mri_convert(parms.fname_in,fname_out,...
    'options',sprintf('-rl %s',parms.fname_raw),...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = cleanup_tmpdir(parms)
  if exist(parms.tmpdir,'dir')
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,msg] = unix(cmd);
    if status
      error('tmpdir cleanup failed:\n%s\n%s\n',cmd,msg);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

