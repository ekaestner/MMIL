function fname_out = fs_surf2surf(fname_in,subject,varargin)
%function fname_out = fs_surf2surf(fname_in,subject,[options])
%
% Purpose: resample a surface stats file (mgh or w formats)
%   (wrapper for freesurfer's mri_surf2surf)
%
% Usage:
%  fname_out = fs_surf2surf(fname_in,subject,'key1', value1,...); 
%
% Required Parameters:
%   fname_in:  full pathname of input file
%   subject: freesurfer recon subject name
%
% Optional parameters:
%  'trgsubj': target subject (e.g. fsaverage, ico)
%     if empty, use subject
%     {default = []}
%  'fname_out' : full path of output file
%     if ommitted, fname_out will be constructed from fname_in
%     {default = []}
%  'outdir' : directory for output file
%     if ommitted, will use same directory as input file
%     if fname_out is supplied, outdir is ignored
%     {default = []}
%  'hemi': cortical hemisphere (lh or rh)
%     necessary only if input file name does not have hemi tag at end
%       e.g. 'stem-lh.mgh' or 'stem-rh.mgh'
%     if file name does have hemi tag, this parameter will be ignored
%     {default = []}
%  'outtype': output file type ('mgh' or 'mgz')
%     if empty, will be set to intype
%     {default = []}
%  'icolevel': icosahedron order number used if trgsubj = 'ico':
%               Order  Number of Vertices
%                 0              12
%                 1              42
%                 2             162
%                 3             642
%                 4            2562
%                 5           10242
%                 6           40962
%                 7          163842
%    {default = 7}
%  'smooth_in': smoothing steps on surface before sampling to trgsubject
%    {default = 0}
%  'smooth_out': smoothing steps on surface after sampling to trgsubject
%    {default = 0}
%  'cortex_flag': apply smoothing only inside cortex mask
%    only applies if $FREESURFER_VER >= 400
%    {default = 0}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'options': string containing additional command-line options for mri_surf2surf
%    {default = []}
%  'verbose': [0|1] display mri_surf2surf output
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Created:  05/30/07 by Ben Cipollini
% Last Mod: 10/14/14 by Don Hagler
%

% created as fs_surf2ico 05/30/07 by Ben Cipoolini

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin, 2), return; end;
parms = mmil_args2parms(varargin,{...
  'trgsubj',[],[],...
  'fname_out',[],[],...
  'outdir',[],[],...            
  'hemi',[],[{'lh','rh'}],...
  'intype','mgh',{'mgh','mgz','w'},...
  'outtype',[],{'mgh','mgz','w'},...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'smooth_in',0,[0 Inf],...
  'smooth_out',0,[0 Inf],...
  'cortex_flag',false,[false true],...
  'icolevel',7,[0 7],...
  'options',[],[],...
  'forceflag',false,[false true],...
  'verbose',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check fname_in
if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

% extract the file stem
[indir,inname,inext] = fileparts(fname_in);
if strcmp(parms.intype,'mgh') && strcmp(inext,'.mgz')
  parms.intype = 'mgz';
end;
if ~strcmp(['.' parms.intype],inext)
  error('fname_in %s missing correct file extension for intype %s',...
    fname_in,parms.intype);
end;
if isempty(parms.outtype)
  parms.outtype = parms.intype;
  parms.outext = inext;
else
  parms.outext = ['.' parms.outtype];
end;

% check subject
fspath = [parms.subjdir '/' subject];
if ~exist(fspath,'dir'), error('fs recon directory %s not found',fspath); end;

% check trgsubj
if isempty(parms.trgsubj), parms.trgsubj = subject; end;
if ~ismember(parms.trgsubj,{subject,'ico'})
  fspath = [parms.subjdir '/' parms.trgsubj];
  if ~exist(fspath,'dir'), error('fs recon directory %s not found',fspath); end;
end;
if strcmp(parms.trgsubj,'ico')
  parms.ico_flag = 1;
else
  parms.ico_flag = 0;
end;

pat = sprintf('(?<instem>.+)-(?<hemi>[lr]h)');
n = regexp(inname,pat,'names');
if isempty(n)
  if isempty(parms.hemi)
    error('input file name %s missing hemi tag (lh or rh) at end of name and hemi parameter not specified',inname);
  else
    fprintf('%s: WARNING: input file name %s missing hemi tag (lh or rh) at end of name',inname);
  end;
  instem = inname;
else
  instem = n.instem;
  parms.hemi = n.hemi;
end;
if isempty(parms.outdir)
  parms.outdir = indir;
end;

% construct the output file name
if isempty(parms.fname_out)
  outstem = instem;
  if (parms.smooth_in > 0)
    outstem = sprintf('%s-sm%d',outstem,parms.smooth_in);
  end;
  if parms.ico_flag
    outstem = sprintf('%s-ico%d',outstem,parms.icolevel);
  end;
  if (parms.smooth_out > 0)
    outstem = sprintf('%s-sm%d',outstem,parms.smooth_out);
  end;
  parms.fname_out = sprintf('%s-%s%s',outstem,parms.hemi,parms.outext);
  parms.fname_out = fullfile(parms.outdir,parms.fname_out);
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

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_out,'file') || parms.forceflag
  % construct the unix command
  cmd = sprintf('setenv SUBJECTS_DIR %s\n',parms.subjdir);
  cmd = sprintf('%smri_surf2surf --srcsubject %s --trgsubject %s --hemi %s',...
    cmd,subject,parms.trgsubj,parms.hemi);
  if parms.ico_flag
    cmd = sprintf('%s --trgicoorder %d',cmd,parms.icolevel);
  end;
  cmd = sprintf('%s --sval %s --tval %s',cmd,...
                fname_in,parms.fname_out);
  cmd = sprintf('%s --sfmt %s --tfmt %s',cmd,...
                parms.intype,parms.outtype);
  if parms.smooth_in>0
    cmd = sprintf('%s --nsmooth-in %d',cmd,parms.smooth_in);
  end;
  if parms.smooth_out>0
    cmd = sprintf('%s --nsmooth-out %d',cmd,parms.smooth_out);
  end;
  if parms.cortex_flag
    cmd = sprintf('%s --cortex',cmd);
  end;
  cmd = sprintf('%s --noreshape',cmd);
  if ~isempty(parms.options)
    cmd = sprintf('%s %s',cmd,parms.options);
  end;
  if parms.verbose, display(cmd); end;
  % run the unix command
  [status,result]=unix(cmd);
  % check for success or failure
  if status || ~isempty(findstr(result,'could not open')) ||...
     ~exist(parms.fname_out,'file')
    error('%s\n\ncmd:\n%s',result,cmd);
  end;
  if parms.verbose, display(result); end;
end;

fname_out = parms.fname_out;


