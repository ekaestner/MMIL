function fs_thresh_surfstats(fname_in,fname_out,varargin);
%function fs_thresh_surfstats(fname_in,fname_out,[options]);
%
% Usage:
%  fs_thresh_surfstats(fname_in,fname_out,'key1', value1,...);
%
% Required input:
%  fname_in: full/relative path name for input surface data file (mgh format)
%  fname_out: full/relative path name for output surface data file (mgh format)
%
% Optional parameters:
%  'fname_tstat' - full/relative path name for surface file containing
%     t-stats (mgh format)
%    {default = []}
%  'p_thresh' - probability threshold applied to t-stats in fname_tstat
%    (only applicable if fname_tstat is supplied)
%    {default = 0.01}
%  'dof' - degrees of freedom for stats
%    (only applicable if fname_tstat is supplied)
%    {default = 100}
%  'tails' - [1|2] tails for statistical test
%    (only applicable if fname_tstat is supplied)
%    {default = 2}
%  'val_thresh' - threshold applied to values in fname_in
%    {default = 0}
%  'thresh_abs_flag' [0|1] - toggle thresholding of absolute values
%    if 0, values < val_thresh are set to 0
%    if 1, abs(values) < abs(val_thresh) are set to 0
%    {default = 1}
%  'clust_thresh' - cluster size threshold (mm^2)
%    {default = 0}
%  'subj' - freesurfer subject name (required if clust_thresh>0)
%    {default = []}
%  'hemi' - cortical hemisphere (required if clust_thresh>0)
%    {default = []}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    (required if clust_thresh>0)
%    {default = $SUBJECTS_DIR}
%
% Created:  02/20/07 by Don Hagler
% Rcnt Mod: 08/10/09 by Don Hagler
% Last Mod: 09/13/12 by Don Hagler
%

%% todo: use mmil_args2parms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  error('calling convention {''key'', value, ... } error');
end;

if ~isfield(g,'fname_tstat'), g.fname_tstat = []; end;
if ~isfield(g,'p_thresh'), g.p_thresh = 0.01; end;
if ~isfield(g,'dof'), g.dof = 100; end;
if ~isfield(g,'tails'), g.tails = 2; end;
if ~isfield(g,'val_thresh'), g.val_thresh = 0; end;
if ~isfield(g,'thresh_abs_flag'), g.thresh_abs_flag = 1; end;
if ~isfield(g,'clust_thresh'), g.clust_thresh = 0; end;
if ~isfield(g,'subj'), g.subj = []; end;
if ~isfield(g,'hemi'), g.hemi = []; end;
if ~isfield(g,'subjdir'), g.subjdir = []; end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'fname_tstat' 'p_thresh' 'dof' 'tails' 'val_thresh'...
         'thresh_abs_flag' 'clust_thresh' 'subj' 'hemi' 'subjdir'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
fname_tstat = g.fname_tstat;
p_thresh = g.p_thresh;
dof = g.dof;
tails = g.tails;
val_thresh = g.val_thresh;
thresh_abs_flag = g.thresh_abs_flag;
clust_thresh = g.clust_thresh;
subj = g.subj;
hemi = g.hemi;
subjdir = g.subjdir;
clear g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parameters
if nargin<2, help(mfilename); return; end;  

if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;
if ~isempty(fname_tstat) & ~exist(fname_tstat,'file')
  error('file %s not found',fname_tstat);
end;

if dof<2
  error('dof (%d) must be >= 2',dof);
end;

if ~ismember(tails,[1,2])
  error('tails (tails) must be 1 or 2',tails);
end;

if clust_thresh>0
  if isempty(subj)
    error('subj must be specified if clust_thresh>0');
  end;
  if isempty(hemi)
    error('hemi must be specified if clust_thresh>0');
  end;
  if ~ismember(hemi,{'lh','rh'})
    error('hemi must be lh or rh (is %s)',hemi);
  end;
  if isempty(subjdir)
    subjdir = getenv('SUBJECTS_DIR');
    if isempty(subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  else
    setenv('SUBJECTS_DIR',subjdir);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
vals = mmil_rowvec(fs_load_mgh(fname_in,[],1)); % first frame only
if ~isempty(fname_tstat)
  if tails==1, p_thresh = 2*p_thresh; end;
  t_thresh = abs(tinv(p_thresh,dof));
  tstats = mmil_rowvec(fs_load_mgh(fname_tstat,[],1)); % first frame only
  if size(tstats)~=size(vals)
    error('size of tstats does not match input vals');
  end;
end;

if clust_thresh
  if ~isempty(fname_tstat)
    tmp_vals = tstats;
    tmp_thresh = t_thresh;
  else
    tmp_vals = vals;
    tmp_thresh = val_thresh;
  end;

  [fname_out_path,fname_out_fstem] = fileparts(fname_out);
  [tmp,tmp_fstem] = fileparts(tempname);
  tmp_fstem = fullfile(fname_out_path,[fname_out_fstem '_' tmp_fstem]);

  % convert to w format
  tmp_fname_in = sprintf('%s-%s.w',tmp_fstem,hemi);
  v=[1:length(tmp_vals)];
  fs_write_wfile(tmp_fname_in,tmp_vals,v);
  if ~exist(tmp_fname_in,'file')
    error('failed to create temp file %s for clust_thresh',...
      tmp_fname_in);
  end;
  
  % run mri_surfcluster
  tmp_fname_out = sprintf('%s-clust-%s.w',tmp_fstem,hemi);
  cmd = 'mri_surfcluster';
  cmd = sprintf('%s --in %s',cmd,tmp_fname_in);
  cmd = sprintf('%s --thmin %0.6f',cmd,tmp_thresh);
  cmd = sprintf('%s --subject %s --hemi %s',cmd,subj,hemi);
  cmd = sprintf('%s --sd %s',cmd,subjdir);
  cmd = sprintf('%s --o %s',cmd,tmp_fname_out);
  cmd = sprintf('%s --minarea %0.6f',cmd,clust_thresh);
  fprintf('%s: cmd = %s\n',mfilename,cmd);
  [status,result] = unix(cmd);
  if status
    error('cmd %s failed:\n%s',cmd,result);
  end;
  fprintf('%s\n',result);
  
  % load mask
  if ~exist(tmp_fname_out,'file')
    error('mri_surfcluster output %s not found',...
      tmp_fname_out);
  end;
  [w,v] = fs_read_wfile(tmp_fname_out);
  v = v(find(w)); % exclude 0 vals (shouldn't be there anyway)
  mask = ones(size(vals));
  mask(v) = 1;

  % apply mask to vals
  tmp_vals = zeros(size(vals));
  tmp_vals(v) = vals(v);
  vals = tmp_vals;

  % clean up tmp files
  delete(tmp_fname_in);
  delete(tmp_fname_out);
elseif ~isempty(fname_tstat)
  vals(abs(tstats)<t_thresh)=0;
end;

% apply value threshold
if thresh_abs_flag
  vals(abs(vals)<abs(val_thresh))=0;
else
  vals(vals<val_thresh)=0;
end;

fs_save_mgh(vals,fname_out,eye(4));

