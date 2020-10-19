function fname_out = fs_ico2ico(fname_in,varargin)
%function fname_out = fs_ico2ico(fname_in,[options])
%
% Purpose: resample a surface stats file (mgh format)
%   from one icoshedral order to another
%
% Usage:
%  fs_ico2ico(fname_in,'key1', value1,...); 
%
% Required Parameters:
%   fname_in:  full pathname of input file
%
% Optional parameters:
%  'ico_subj': name of icosahedral freesurfer subject
%     {default: 'fsaverage'}
%  'fname_out': full path of output file
%     if ommitted, fname_out will be constructed from fname_in
%     {default: []}
%  'presmooth': number of smooth steps applied before resampling
%    {default = 0}
%  'sparsesmooth': number of sparse smoothing steps
%    {default = 0}
%    note: sparse smoothing is a fast way to fill
%          in gaps between sparse vertices
%  'postsmooth': number of normal smoothing steps
%    {default = 0}
%    note: postsmoothing is additional nearest-neighbor average
%          smoothing applied after sparse smoothing
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the avgerage subject directory
%    {default = $FREESURFER_HOME/subjects}
%  'forceflag': [0|1] toggle overwrite existing output files
%    {default: 0}
%
%
% NOTES:
%   Input file must have number of values corresponding to one of the
%     following icosahedral orders:
%
%   Order  Number of Vertices
%     1              42
%     2             162
%     3             642
%     4            2562
%     5           10242
%     6           40962
%     7          163842
%
%  Output file will have ico order corresponding to ico_subj  
%
% Created:       09/20/08 by Don Hagler
% Last Modified: 11/05/12 by Don Hagler
%

% 02/23/09: Andrei Irimia: changed vol to invals, sparsesmooth to parms.sparsesmooth

% 10/27/10: Don Hagler: use mri_surf2surf to go from one ico to another

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns = [42 162 642 2562 10242 40962 163842];
fname_out = [];

if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin, {...
  'ico_subj','fsaverage',[],...
  'fname_out',[],[],...
  'presmooth',0,[0 1000],...
  'sparsesmooth',0,[0 1000],...
  'postsmooth',0,[0 1000],...
  'subjdir',[],[],...
  'forceflag',false,[false true],...
... % hidden
  'hemi','lh',{'lh','rh'},...
});

if isempty(parms.subjdir)
  tmp = getenv('FREESURFER_HOME');
  if isempty(tmp)
    error('FREESURFER_HOME environment variable not set');
  end;
  parms.subjdir = [tmp '/subjects'];
end;

fullsubj = [parms.subjdir '/' parms.ico_subj];
if ~exist(fullsubj,'dir')
  error('ico subj %s not found',fullsubj);
end;

surf = fs_load_subj(parms.ico_subj,parms.hemi,[],1,parms.subjdir); % read nverts only
ico_out = find(surf.nverts==Ns);
if isempty(ico_out)
  error('ico subj %s has non-ico number of vertices (%d)',fullsubj,surf.nverts);
end;

% set fname_out
if isempty(parms.fname_out)
  [indir,inname,inext] = fileparts(fname_in);
  pat = sprintf('(?<instem>.+)-(?<hemi>[lr]h)');
  n = regexp(inname,pat,'names');
  if isempty(n)
    fprintf('%s: WARNING: input file name %s missing hemi tag (lh or rh) at end of name',inname);
    instem = inname;
    hemi = [];
  else
    instem = n.instem;
    hemi = n.hemi;
    parms.hemi = hemi;
  end;
  parms.fname_out = sprintf('%s/%s-ico%d',indir,instem,ico_out);
  if parms.presmooth>0
    parms.fname_out = sprintf('%s-presm%d',parms.fname_out,parms.presmooth);
  end;
  if parms.sparsesmooth>0
    parms.fname_out = sprintf('%s-spsm%d',parms.fname_out,parms.sparsesmooth);
  end;
  if parms.postsmooth>0
    parms.fname_out = sprintf('%s-sm%d',parms.fname_out,parms.postsmooth);
  end;
  if ~isempty(hemi)
    parms.fname_out = sprintf('%s-%s',parms.fname_out,hemi);
  end;
  parms.fname_out = [parms.fname_out '.mgh'];
end;

if exist(parms.fname_out) && ~parms.forceflag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find number of vertices in input file
[tmp,volsz] = fs_read_header(fname_in);
if volsz(2)>1 || volsz(3)>1
  error('input file %s appears to contain volume data',fname_in);
end;
nverts = volsz(1);
nframes = volsz(4);
ico_in = find(nverts==Ns);
if isempty(ico_in)
  error('input file has non-ico number of vertices (%d)',nverts);
end;

% set tmpfile
if parms.sparsesmooth || parms.postsmooth
  [fpath,fstem,fext] = fileparts(parms.fname_out);
  tmpfile = [fpath '/' fstem '-tmp' fext];
else
  tmpfile = parms.fname_out;
end;

% use mri_surf2surf to resample
fprintf('%s: resampling from ico%d to ico%d...\n',mfilename,ico_in,ico_out);
cmd = sprintf('setenv SUBJECTS_DIR %s\n',parms.subjdir); 
cmd = sprintf('%smri_surf2surf',cmd);
cmd = sprintf('%s --srcsubject ico --srcicoorder %d',cmd,ico_in);
%cmd = sprintf('%s --trgsubject ico --trgicoorder %d',...
%  cmd,ico_out); %% setting trgicoorder does not work for some reason
cmd = sprintf('%s --trgsubject %s',cmd,parms.ico_subj);
cmd = sprintf('%s --sval %s --tval %s',cmd,fname_in,tmpfile);
if parms.presmooth
  cmd = sprintf('%s --nsmooth-in %d',cmd,parms.presmooth);
end;
cmd = sprintf('%s --hemi %s',cmd,parms.hemi);

[s,r] = unix(cmd);
if s, error('cmd %s failed:\n%s',cmd,r); end;

fname_out = tmpfile;
if strcmp(tmpfile,parms.fname_out), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load input file
[invals,tmp,tmp,volsz] = fs_load_mgh(tmpfile);
nverts = volsz(1);
nframes = volsz(4);
invals = reshape(invals,[nverts,nframes]);

% load ico_subj
surf = fs_load_subj(parms.ico_subj,parms.hemi,[],[],parms.subjdir);

% create output
outvals = zeros(nverts,nframes);

% apply smoothing
fprintf('%s: smoothing %d frames...\n',mfilename,nframes);
for f=1:nframes
  tmpvals = invals(:,f);
  if parms.sparsesmooth>0
    tmpvals = fs_smooth_sparse(surf,tmpvals,parms.sparsesmooth);
  end;
  if parms.postsmooth>0
    tmpvals = fs_smooth(surf,tmpvals,parms.postsmooth);
  end;
  outvals(:,f) = tmpvals;
end;

% save output
outvals = reshape(outvals,[nverts,1,1,nframes]);
fs_save_mgh(outvals,parms.fname_out);
fname_out = parms.fname_out;

delete(tmpfile);


