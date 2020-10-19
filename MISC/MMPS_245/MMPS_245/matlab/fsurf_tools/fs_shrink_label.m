function fs_shrink_label(fname,varargin)
%function fs_shrink_label(fname,[options])
%
% purpose: make a smaller version of a label
%
% Required Input:
%   fname: file name of label file
%
% Optional Input:
%  'fname_out': output file name
%    If empty, will append -sm to stem of fname
%    {default: []}
%  'surf' surface struct (see fs_load_subj)
%    If empty, will load using 'subj' and 'subjdir'
%    {default: []}
%  'subj': FreeSurfer subject name
%    {default: 'fsaverage'}
%  'subjdir': FreeSurfer subject root dir
%    {default: $FREESURFER_HOME/subjects}
%  'smooth': number of smoothing steps
%    {default: 10}
%  'thresh': threshold applied after smoothing
%    {default: 0.5}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% created:  10/20/10 Don Hagler
% last mod: 10/20/10 Don Hagler
%

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_out',[],[],...
  'surf',[],[],...
  'subj','fsaverage',[],...
  'subjdir',[],[],...
  'hemi','lh',[],...
  'smooth',10,[0,1000],...
  'thresh',0.5,[0,1],...
  'forceflag',false,[false true],...
});

if isempty(parms.surf)
  if isempty(parms.subjdir)
    parms.subjdir = [getenv('FREESURFER_HOME') '/subjects'];
  end;

  if ~exist(parms.subjdir,'dir')
    error('directory %s not found',parms.subjdir);
  end;

  subjpath = [parms.subjdir '/' parms.subj];
  if ~exist(subjpath,'dir')
    error('directory %s not found',subjpath);
  end;

  parms.surf = fs_load_subj(parms.subj,parms.hemi,[],[],parms.subjdir);
end;

if isempty(parms.fname_out)
  [tpath,tstem,text]=fileparts(fname);
  parms.fname_out = sprintf('%s/%s-sm%d-thresh%0.2f%s',...
    tpath,tstem,parms.smooth,parms.thresh,text);
end;

vertices = fs_read_label(fname);

vals = zeros(parms.surf.nverts,1);
vals(vertices) = 1;

vals = fs_smooth(parms.surf,vals,parms.smooth);

vals = vals/max(eps,max(vals));

vertices = find(vals>parms.thresh);

fs_write_label(vertices,parms.fname_out,parms.subj);

