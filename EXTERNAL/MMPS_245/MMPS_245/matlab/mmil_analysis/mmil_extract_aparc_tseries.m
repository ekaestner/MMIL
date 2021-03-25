function [data,roinames] = mmil_extract_aparc_tseries(subj,fname,varargin);
%function [data,roinames] = mmil_extract_aparc_tseries(subj,fname,[options]);
%
% Purpose:
%   paint 4D volume to surface and extract time series for
%     each ROI in annotation (e.g. aparc)
%
% Required parameters:
%  subj:  string specifying the subject name
%  fname: full or relative path of 4D functional volume
%    (must be mgh/mgz format)
%
% Optional parameters:
%  'fnames_surf': cell array of surface data files
%     must have two elements, first for left hemi, second for right hemi
%     If supplied, fname will be ignored (may be empty)
%     {default = []}
%  'outdir' : output directory for painted surface files
%    If empty, will attempt to write to path containing fname
%    {default = []}
%  'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     if empty, will use ?h.aparc.annot files in subjdir/subj/label
%    {default = []}
%  'skipTRs': number of frames to remove from beginning of fname
%    {default = 0}
%  'regfile': register.dat file containing 4x4 registration matrix
%    If not supplied, fname should be resampled to structural space
%    {default = []}
%  'projfrac_flag': [0|1] whether to use projdist (0) or projfract (1)
%    {default = 0}
%  'projdist': distance (mm) to project along surface vertex normal
%    {default = 1}
%  'projfrac': fractional distance to project along surface vertex normal
%    relative to cortical thickness
%    {default = 0.5}
%  'projfrac_avg': vector of [min max del] for averaging multiple samples
%     with mri_surf2surf projfract-avg option
%    If empty, use projfrac instead if projfrac_flag=1
%    {default = []}
%  'projdist_avg': vector of [min max del] for averaging multiple samples
%     with mri_surf2surf projdist-avg option
%    If empty, use projdist instead
%    {default = []}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'surfname': name of surface onto which to sample volume data
%    {default = white}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Output:
%   data: matrix of time x ROI
%   roinames: cell array of aparc ROI names
%
% Created:  03/06/12 by Don Hagler
% Prev Mod: 05/14/13 by Don Hagler
% Last Mod: 06/30/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

data = []; roinames = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'subj',subj,[],...
  'fname_in',fname,[],...
  'fnames_surf',[],[],...
  'outdir',[],[],...
  'fnames_aparc',[],[],...
  'skipTRs',0,[0,1000],...
  'regfile',[],[],...
  'projfrac_flag',false,[false true],...
  'projdist',1,[-10,10],...
  'projfrac',0.5,[-2,2],...
  'projdist_avg',[],[],...
  'projfrac_avg',[],[],...
  'subjdir',[],[],...
  'surfname','white',[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
});

paint_tags = {
  'regfile' 'projfrac_flag' 'projdist' 'projfrac' 'projdist_avg'...
  'projfrac_avg' 'subjdir' 'surfname' 'forceflag' 'hemilist'...
  'outstem'...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

parms.nhemi = length(parms.hemilist);

% set subjdir if not supplied
if isempty(parms.subjdir)
  parms.subjdir = getenv('SUBJECTS_DIR');
  if isempty(parms.subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;

% check FS recon exists
parms.fspath = sprintf('%s/%s',parms.subjdir,parms.subj);
if ~exist(parms.fspath,'dir')
  error('FreeSurfer recon dir %s not found',parms.fspath);
end;

% check fnames_aparc
if isempty(parms.fnames_aparc)
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    parms.fnames_aparc{h} = sprintf('%s/label/%s.aparc.annot',...
      parms.fspath,hemi);
  end;
else
  if ~iscell(parms.fnames_aparc), parms.fnames_aparc = {parms.fnames_aparc}; end;
  if length(parms.fnames_aparc) ~= parms.nhemi
    error('must have %d elements in fnames_aparc (have %d)',...
      parms.nhemi,length(parms.fnames_aparc));
  end;
end;
for h=1:parms.nhemi
  if ~exist(parms.fnames_aparc{h},'file')
    error('annot file %s not found',parms.fnames_aparc{h});
  end;
end;

% check fnames_surf
if ~isempty(parms.fnames_surf)
  if ~iscell(parms.fnames_surf) || length(parms.fnames_surf)~=parms.nhemi
    error('fnames_surf must be a cell array with %d elements',parms.nhemi);
  end;
  for h=1:parms.nhemi
    if ~exist(parms.fnames_surf{h},'file')
      error('file %s not found',parms.fnames_surf{h});
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fnames_surf)
  % set outdir and outstem
  [tmp_outdir,tmp_outstem,tmp] = fileparts(parms.fname_in);
  if isempty(parms.outdir), parms.outdir = tmp_outdir; end;
  parms.outstem = [parms.outdir '/' tmp_outstem];

  mmil_mkdir(parms.outdir);
  % sample 4D data to surface
  args = mmil_parms2args(parms,paint_tags);
  parms.fnames_surf = fs_paint(subj,fname,args{:});
end;

% determine number of frames
%[~,volsz] = fs_read_header(parms.fnames_surf{1});
[M,volsz] = mmil_load_mgh_info(parms.fnames_surf{1},parms.forceflag,parms.outdir);

nframes = volsz(4);
if parms.skipTRs>0
  frames = setdiff([1:nframes],[1:parms.skipTRs]);
else
  frames = [1:nframes];
end;
ntpoints = length(frames);

% extract time series for each cortical ROI
nroi = 0;
aparc_results = [];
for h=1:parms.nhemi
  hemi = parms.hemilist{h};
  fname_in = parms.fnames_surf{h};
  fname_annot = parms.fnames_aparc{h};
  aparc_results{h} = mmil_surf_roi(fname_in,...
    'fname_aparc',fname_annot,'frames',frames,'verbose',parms.verbose);
  nroi = nroi + length(aparc_results{h});
end;

% combine lh and rh aparc
k = 1;
roinames = cell(1,nroi);
data = zeros(ntpoints,nroi);
for h=1:parms.nhemi
  hemi = parms.hemilist{h};
  for r=1:length(aparc_results{h})
    tmp_roiname = aparc_results{h}(r).roiname;
    if isempty(regexp(tmp_roiname,sprintf('^ctx-%s',hemi)))
      tmp_roiname = sprintf('ctx-%s-%s',hemi,tmp_roiname);
    end;
    roinames{k} = tmp_roiname;
    data(:,k) = aparc_results{h}(r).avg;
    k = k + 1;
  end;
end;

return;

