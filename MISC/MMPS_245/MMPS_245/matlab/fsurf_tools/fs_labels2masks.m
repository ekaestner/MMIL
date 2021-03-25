function fs_labels2masks(subj,varargin)
%function fs_labels2masks(subj,varargin)
%
% Usage:
%  fs_labels2masks(subj,'key1', value1,...);
%
% Required parameters:
%   subj:  a string specifying the subject name (freesurfer recon dir)
%
% Optional parameters:
%  'indir' - input directory containing label files
%     {default = pwd}
%  'outdir' - output directory to contain mask files
%     {default = pwd}
%  'outtype' - output file type ('mgh' or 'w')
%    {default: 'mgh'}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'forceflag' - [0|1] whether to overwrite existing output files
%    {default: 0}
%
% Created:  10/09/08 by Don Hagler
% Last Mod: 10/09/08 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms( varargin, {...
  'indir',pwd,[],...
  'outdir',pwd,[],...
  'outtype','mgh',{'mgh','w'},...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false true],...
... % hidden parms
  'surf','white',[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if ~isempty(parms.subjdir)
  setenv('SUBJECTS_DIR',parms.subjdir)
else
  error('FreeSurfer SUBJECTS_DIR environment variable not set');
end;

if ~exist(parms.indir,'dir'), error('directory %s not found',parms.indir); end;
if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to create output dir %s:\n%s',parms.outdir,msg);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  flist = dir(sprintf('%s/%s.*.label',parms.indir,hemi));
  if isempty(flist), continue; end;
  if strcmp(parms.outtype,'mgh')
    surf = fs_load_subj(subj,hemi,parms.surf,1,parms.subjdir);
  end;
  for f=1:length(flist)
    name = flist(f).name;
    pat = sprintf('^%s\\.(?<stem>.+)\\.label$',hemi);
    n = regexp(flist(f).name,pat,'names');
    if isempty(n)
      fprintf('%s: WARNING: failed to match pattern %s to name %s',...
        pat,name);
      continue;
    end;
    fname_in = sprintf('%s/%s',parms.indir,name);
    fname_out = sprintf('%s/%s-%s.%s',parms.outdir,n.stem,hemi,parms.outtype);
    fprintf('%s: converting %s to %s...\n',mfilename,fname_in,fname_out);
    if ~exist(fname_out,'file') || parms.forceflag
      v = fs_read_label(fname_in);
      if strcmp(parms.outtype,'mgh')
        vals = zeros(surf.nverts,1);
        vals(v) = 1;
        fs_save_mgh(vals,fname_out);
      else
        w = ones(length(v),1);
        fs_write_wfile(fname_out,w,v);
      end;
    end;
  end;
end;
