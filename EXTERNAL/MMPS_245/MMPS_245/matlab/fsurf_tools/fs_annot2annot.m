function fname_out = fs_annot2annot(fname,varargin)
%function fname_out = fs_annot2annot(fname,[options])
%
% Purpose: resample an annotation file from one subject to another
%   by first converting to labels, resampling those labels,
%   then creating a new annotation file
%
% Required Input:
%   fname: full or relative path of annotation file
%
% Optional Input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     if not supplied, will use name from input file
%     {default = []}
%   'source_subj': FreeSurfer subject from whom labels were made
%     {default = 'fsaverage'}
%   'subj': FreeSurfer subject for whom annotation will be made
%     {default = 'fsaverage'}
%   'subjdir': FreeSurfer subject root directory
%     {default = $FREESURFER_HOME/subjects}
%   'hemi': hemisphere string ('lh' or 'rh')
%     not required if fname has form {hemi}.{name}.annot
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_out: output file name
%
% Created:  11/29/12 by Don Hagler
% Prev Mod: 11/29/12 by Don Hagler
% Last Mod: 10/19/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
fname_out = [];

parms = mmil_args2parms(varargin,{...
  'outdir',pwd,[],...
  'outstem',[],[],...
  'source_subj','fsaverage',[],...
  'subj','fsaverage',[],...
  'subjdir',[],[],...
  'hemi',[],[],...
  'forceflag',false,[false true],...
...
  'verbose',true,[false true],...
  'cleanup_flag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.subjdir)
  parms.subjdir = sprintf('%s/subjects',getenv('FREESURFER_HOME'));
end;

% determine hemi using regexp
[tpath,tstem,text] = fileparts(fname);
n = regexp(tstem,'(?<hemi>[lr]h)\.(?<name>.+)','names');
if isempty(n)
  fprintf('%s: WARNING: fname %s has unexpected pattern...\n',...
    mfilename,fname);
else
  parms.hemi = n.hemi;
end;
if isempty(parms.hemi)
  error('hemi not specified');
end;
if isempty(parms.outstem)
  parms.outstem = n.name;
end;

if mmil_isrelative(parms.outdir)
  parms.outdir = [pwd '/' parms.outdir];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = sprintf('%s/%s.%s.annot',parms.outdir,parms.hemi,parms.outstem);

if ~exist(fname_out,'file') || parms.forceflag
  tmpdir = [parms.outdir '/tmp_labels_source_subj'];
  mmil_mkdir(tmpdir);
  % convert annotation file to labels
  [label_fnames,ctab] = fs_annot2labels(fname,...
    'outdir',tmpdir,'subj',parms.source_subj,'hemi',parms.hemi,...
    'verbose',parms.verbose,'forceflag',parms.forceflag);
  % create ctab file
  fname_ctab = sprintf('%s/%s.%s.ctab',parms.outdir,parms.hemi,parms.outstem);
  if ~exist(fname_ctab,'file') || parms.forceflag
    roinames = ctab.struct_names;
    roicolors = ctab.table(:,1:3);
    fs_write_ctab(roinames,fname_ctab,roicolors);
  end;
  % resample labels and create annotation file
  [annot_fnames,label_fnames] = fs_labels2annot(tmpdir,...
    'fname_ctab',fname_ctab,...
    'outdir',parms.outdir,...
    'source_subj',parms.source_subj,...
    'subj',parms.subj,...
    'subjdir',parms.subjdir,...
    'annotname',parms.outstem,...
    'hemilist',{parms.hemi},...
    'verbose',parms.verbose,...
    'forceflag',parms.forceflag);
  % remove temporary files
  if parms.cleanup_flag
    cmd = sprintf('rm -r %s',tmpdir);
    [s,r] = unix(cmd);
    if s, error('cmd %s failed:\n%s',cmd,r); end;
  end;
end;
