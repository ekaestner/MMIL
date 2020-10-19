function fname_out = fs_label_groupavg(labels,subjnames,varargin)
%function fname_out = fs_label_groupavg(labelist,subjnames,[options])
%
% Required Input:
%   labels: cell array of label file names
%   subjnames: cell array of subject names
%     NOTE: length of these arrays must match
%
% Optional Parameters:
%  'outdir': output directory
%    will contain labels for each subject, resampled to ico
%    {default = pwd}
%  'outstem': output file stem (do not include hemi or extension)
%    {default = 'groupavg'}
%  'outtype': output file type ('mgh' or 'label')
%    {default = 'mgh'}
%  'hemi': cortical hemisphere for label files
%    if empty, will attempt to determine hemisphere from first input label file name
%      for that to work, label files should be like "lh.name.label" or "rh.name.label"
%    {default = []}
%  'icolevel': icosahedron order number:
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
%  'smooth_subj': smoothing steps on native subject surface before sampling to ico
%    {default = 0}
%  'smooth_ico': smoothing steps on surface after sampling to ico
%    {default = 0}
%  'cortex_flag': apply smoothing only inside cortex mask
%    only applies if $FREESURFER_VER >= 400
%    {default = 0}
%  'thresh': threshold applied to average of binary masks
%    {default = 0}
%  'subjdir': root directory for subject directories
%    {default = getenv(SUBJECTS_DIR)}
%  'verbose': [0|1] display mri_surf2surf output
%   {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Output:
%   fname_out: output file name for average label
%
% Created:  03/15/11 by Don Hagler
% Last Mod: 07/19/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',pwd,[],...
  'outstem','groupavg',[],...
  'outtype','mgh',{'mgh','label'},...
  'hemi',[],{'lh','rh'},...
  'smooth_subj',0,[0 Inf],...
  'smooth_ico',0,[0 Inf],...
  'cortex_flag',false,[false true],...
  'icolevel',7,[0 7],...
  'thresh',0,[0,1],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false true],...
  'verbose',false,[false true],...
...
  'fpat','(?<hemi>[l,r]h)\.(?<name>.+)\.label',[],...
});

if ~iscell(labels)
  error('labels must be cell array of label file names');
end;
if ~iscell(subjnames)
  error('subjnames must be cell array of subject names');
end;

nlabels = length(labels);
if length(subjnames) ~= nlabels
  error('number of labels (%d) must match number of subjects (%d)',...
    nlabels,length(subjnames));
end;

if isempty(parms.subjdir)
  error('must specify subjdir or set SUBJECTS_DIR environment variable');
end;

if mmil_isrelative(parms.outdir)
  parms.outdir = [pwd '/' parms.outdir];
end;
parms.indivdir = [parms.outdir '/individuals'];

if isempty(parms.hemi)
  [name,hemi] = check_label_name(labels{1},parms.fpat);
  if isempty(hemi)
    error('label file %s does not match pattern hemi.name.label and hemi was not specified',...
      labels{1});
  else
    parms.hemi = hemi;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create group average label

[succ,msg] = mkdir(parms.outdir);
if ~succ, error('failed to create output directory %s\n:%s',parms.outdir,msg); end;
[succ,msg] = mkdir(parms.indivdir);
if ~succ, error('failed to create output directory %s\n:%s',parms.indivdir,msg); end;

switch parms.outtype
  case 'mgh'
    fname_out = [parms.outdir '/' parms.outstem '-' parms.hemi '.' parms.outtype];
  case 'label'
    fname_out = [parms.outdir '/' parms.hemi '.' parms.outstem '.' parms.outtype];
end;
if ~exist(fname_out,'file') || parms.forceflag
  fnamelist = cell(1,nlabels);
  for i=1:nlabels
    [name,hemi] = check_label_name(labels{i},parms.fpat);
    if ~isempty(hemi) && ~strcmp(hemi,parms.hemi)
      error('hemi (%s) for label file %s does not match specified hemi or that of first label (%s)',...
        hemi,labels{i},parms.hemi);
    end;
    fname_mgh = sprintf('%s/%s-%s-%s.mgh',...
      parms.indivdir,subjnames{i},name,parms.hemi);

    % convert from label file to mgh
    fs_label2mgh(labels{i},fname_mgh,subjnames{i},...
      'hemi',parms.hemi,'subjdir',parms.subjdir,'forceflag',parms.forceflag);

    % resample to ico
    fnamelist{i} = fs_surf2surf(fname_mgh,subjnames{i},'trgsubj','ico',...
      'smooth_in',parms.smooth_subj,...
      'smooth_out',parms.smooth_ico,...
      'icolevel',parms.icolevel,...
      'subjdir',parms.subjdir,...
      'verbose',parms.verbose,...
      'forceflag',parms.forceflag);
  end;
%      'cortex_flag',parms.cortex_flag',...

  % call fs_groupavg
  results = fs_groupavg(fnamelist,'stats_flag',0);

  % apply threhsold
  vals = results.groups.mean;
  if parms.thresh, vals(vals<parms.thresh) = 0; end;
  
  % save output file as mgh or label
  if strcmp(parms.outtype,'label')
    fs_write_label(find(vals),fname_out);
  else
    fs_save_mgh(vals,fname_out);
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [name,hemi] = check_label_name(fname_label,fpat)
  name = [];
  hemi = [];
  if isempty(fname_label), error('fname_label is empty'); end;
  [tpath,tstem,text] = fileparts(fname_label);
  tname = [tstem text];
  n = regexp(tname,fpat,'names');
  if ~isempty(n)
    name = n.name;
    hemi = n.hemi;
  else
    name = tstem;
  end;
return;

