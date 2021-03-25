function fname_out = fs_wmseg(fname_in,varargin)
%function fname_out = fs_wmseg(fname_in,varargin)
%
% Purpose: create white matter segmentation volume from T1-weighted image
%
% Usage: fs_wmseg(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_out': output file name
%    {default = 'wmseg.mgz'}
%  'fname_mask': input brain mask file name
%    if not supplied, will run mri_watershed to create brainmask
%    {default = []}
%  'conform_flag': [0|1] convert input image to be 256^3 and uchar type
%    {default = 1}
%  'nu_flag': [0|1] run N3 non-uniformity intensity correction
%    {default = 1}
%  'nu3T_flag': [0|1] whether to use nu intensity correction
%    with settings optimized for 3T
%    irrelevant if nu_flag = 0
%    {default = 0}
%  'tal_flag': [0|1] register to Talairach space
%    for improved non-uniformity correction
%    irrelevant if nu_flag = 0
%    {default = 1}
%  'conform_options': options used when calling mri_convert to conform input
%    {default = '--conform -rt cubic -odt uchar -ns 1'}
%  'fname_talxfm': [0|1] existing Talairach transform file
%    if supplied, will skip Talairach registration
%    {default = []}
%  'thicken_flag': [0|1] make thin strands of white matter thicker
%    {default = 0}
%  'fillbg_flag': [0|1] fill basal ganglia
%    {default = 0}
%  'fillv_flag': [0|1] fill ventricles
%    {default = 0}
%  'tmpdir': temporary directory containing intermediate output
%    {default = 'tmp_wmseg'}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 1}
%  'verbose': [0|1] display progress messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  12/21/14 by Don Hagler
% Last Mod: 05/12/15 by Don Hagler
% Last Mod: 07/13/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output file
if exist(parms.fname_out) && ~parms.forceflag, return; end;

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

% create brain mask
if isempty(parms.fname_mask)
  parms.fname_mask = create_brainmask(parms,'brainmask');
end;

% T1 intensity normalization
parms.fname_in = T1_correction(parms,'brain');

% white matter segmentation
wm_segmentation(parms);

% remove temporary files
if parms.cleanup_flag
  parms = cleanup_tmpdir(parms);
end;

fname_out = parms.fname_out;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_out','wmseg.mgz',[],...
    'fname_mask',[],[],...
    'conform_flag',true,[false true],...
    'nu_flag',true,[false true],...
    'nu3T_flag',false,[false true],...
    'tal_flag',true,[false true],...
    'conform_options','--conform -rt cubic -odt uchar -ns 1',[],...
    'fname_talxfm',[],[],...
    'thicken_flag',false,[false true],...
    'fillbg_flag',false,[false true],...
    'fillv_flag',false,[false true],...
    'tmpdir','tmp_wmseg',[],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'nu_tags',{'fname_out','nu3T_flag','tal_flag','fname_talxfm',...
               'tmpdir','cleanup_flag','verbose','forceflag'},[],...
    'T1_tags',{'fname_out','fname_mask','verbose','forceflag'},[],...
    'brainmask_tags',{'fname_out','conform_flag','nu_flag','T1_flag',...
                      'surfs_flag','nu3T_flag','tal_flag',...
                      'conform_options','fname_talxfm','surfs_outdir',...
                      'surfs_nverts','watershed_options','tmpdir',...
                      'cleanup_flag','verbose','forceflag'},[]...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if ~isempty(parms.fname_mask) && ~exist(parms.fname_mask,'file')
    error('file %s not found',parms.fname_mask);
  end;
  if parms.tal_flag &&...
     ~isempty(parms.fname_talxfm) && ~exist(parms.fname_talxfm,'file')
    error('file %s not found',parms.fname_talxfm);
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

function fname_out = create_brainmask(parms,outstem)
  parms.fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  parms.cleanup_flag = 0;
  parms.conform_flag = 0;
  parms.nu_flag = 0;
  parms.T1_flag = 1;
  parms.surfs_flag = 0;
  args = mmil_parms2args(parms,parms.brainmask_tags);
  fname_out = fs_brainmask(parms.fname_in,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = T1_correction(parms,outstem)
  parms.fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  args = mmil_parms2args(parms,parms.T1_tags);
  fname_out = fs_T1_corr(parms.fname_in,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wm_segmentation(parms)
  fname_out = parms.fname_out;
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: segmenting white matter...\n',mfilename);
      tic;
    end;
    outdir = fileparts(parms.fname_out);
    cmd = sprintf('cd %s; mri_segment',outdir);
    if ~parms.thicken_flag
      cmd = sprintf('%s -thicken 0',cmd);
    end;
    if parms.fillbg_flag
      cmd = sprintf('%s -fillbg',cmd);
    end;
    if parms.fillv_flag
      cmd = sprintf('%s -fillv',cmd);
    end;
    cmd = sprintf('%s %s %s',cmd,parms.fname_in,fname_out);
    [status,msg] = unix(cmd);
    if status || ~exist(fname_out,'file')
      error('segmentation failed:\n%s\n%s\n',cmd,msg);
    end;
    if parms.verbose
      toc;
    end;
    % move segment.dat file to tmpdir (for later delection or not)
    fname_segdat = sprintf('%s/segment.dat',outdir);
    if exist(fname_segdat,'file') & ~strcmp(outdir,parms.tmpdir)
      cmd = sprintf('mv %s %s',fname_segdat,parms.tmpdir);
      [status,msg] = unix(cmd);
      if status
        error('failed to move segement dat file:\n%s\n%s\n',cmd,msg);
      end;
    end;
  end;
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

