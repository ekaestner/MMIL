function fs_label2mgh(fname_in,fname_out,subject,varargin)
%function fs_label2mgh(fname_in,fname_out,subject,varargin)
%
% Purpose: convert label file to mgh file
%
% Required Input:
%   fname_in: input filename (label format)
%   fname_out: output filename (mgh format)
%   subject: freesurfer recon dir
%
% Optional parameters:
%  'hemi': cortical hemisphere for label files
%    if empty, will attempt to determine hemisphere from first input label file name
%      for that to work, label files should be like "lh.name.label" or "rh.name.label"
%    {default = []}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'forceflag': [0|1] whether to overwrite existing output files
%    {default: 0}
%
% Created:  03/16/11 by Don Hagler
% Last Mod: 08/14/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms( varargin, {...
  'hemi',[],{'lh','rh'},...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false true],...
... % hidden
  'surf','white',[],...
});

if ~exist(fname_in,'file'), error('file %s not found',fname_in); end;

if isempty(parms.subjdir)
  error('must specify subjdir or set SUBJECTS_DIR environment variable');
end;

if isempty(parms.hemi)
  n = regexp(fname_in,'(?<hemi>[lr]h)\.(?<name>.+)\.label','names');
  if isempty(n)
    error('fname_in (%s) does not match pattern hemi.name.label and hemi was not specified',...
      fname_in);
  else
    parms.hemi = n.hemi;
  end;
end;

if ~exist(fname_out,'file') || parms.forceflag
  surf = fs_load_subj(subject,parms.hemi,parms.surf,1,parms.subjdir);
  v = fs_read_label(fname_in);
  vals = zeros(surf.nverts,1);
  vals(v) = 1;
  fs_save_mgh(vals,fname_out);
end;

