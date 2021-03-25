function mmil_extract_tseries(subj,fname,varargin);
%function mmil_extract_tseries(subj,fname,[options]);
%
% Purpose: extract time series for aparc and aseg ROIs
%   from 4D functional volume
%
% Required parameters:
%  subj:  string specifying the subject name
%  fname: full or relative path of 4D functional volume
%    (must be mgh/mgz format)
%
% Optional parameters:
%  'fname_out': name of output csv file containing time series for each ROI
%    If empty, will construct from file stem of funcname
%      with '_avgtseries.csv' appended
%    If relative, will create in outdir
%    {default = []}
%  'outdir' : output directory for painted surface files
%    If empty, will attempt to write to path containing funcname
%    {default = []}
%  'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     if empty, will use ?h.aparc.annot files in subjdir/subj/label
%    {default = []}
%  'fname_aseg': name of aseg file
%     if empty, will use aseg.mgz in subjdir/subj/mri
%    {default = []}
%  'TR': duration (sec) of each time point in funcname
%    {default = 1}
%  'skipTRs': number of frames to remove from beginning of funcname
%    {default = 0}
%  'regfile': register.dat file containing 4x4 registration matrix
%    If not supplied, funcname should be resampled to structural space
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
%  'aparc_flag': extract time series for aparc cortical surface ROIs
%    {default = 1}
%  'aseg_flag': extract time series for aseg volume ROIs
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Created:  11/30/11 by Don Hagler
% Last Mod: 10/29/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'subj',subj,[],...
  'fname_in',fname,[],...
  'fname_out',[],[],...
  'outdir',[],[],...
  'fnames_aparc',[],[],...
  'fname_aseg',[],[],...
  'TR',1,[0.01,100],...
  'skipTRs',0,[0,1000],...
  'regfile',[],[],...
  'projfrac_flag',false,[false true],...
  'projdist',1,[-10,10],...
  'projfrac',0.5,[-2,2],...
  'projdist_avg',[],[],...
  'projfrac_avg',[],[],...
  'subjdir',[],[],...
  'surfname','white',[],...
  'aparc_flag',true,[false true],...
  'aseg_flag',true,[false true],...
  'forceflag',false,[false true],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79],[1,Inf],...
});

%% todo: use aseg_roilist from function

aparc_tags = {'outdir','fnames_aparc',...
  'skipTRs','regfile','projfrac_flag',...
  'projdist','projfrac','projdist_avg','projfrac_avg',...
  'subjdir','surfname','forceflag','hemilist'};

aseg_tags = {'outdir','fname_aseg',...
  'skipTRs','regfile','subjdir','forceflag','aseg_roilist'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if ~parms.aparc_flag && ~parms.aseg_flag
  error('aparc_flag and aseg_flag cannot both be false');
end;

% set outdir and outstem
[tmp_outdir,tmp_outstem,tmp] = fileparts(parms.fname_in);
if isempty(parms.outdir), parms.outdir = tmp_outdir; end;
parms.outstem = [parms.outdir '/' tmp_outstem];

% set fname_out
if isempty(parms.fname_out)
  parms.fname_out = sprintf('%s_avgtseries.csv',parms.outstem);
end;
if mmil_isrelative(parms.fname_out)
  parms.fname_out = [parms.outdir '/' parms.fname_out];
end;

if exist(parms.fname_out,'file') && ~parms.forceflag
  return;
end;

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.aparc_flag
  args = mmil_parms2args(parms,aparc_tags);
  [aparc_data,aparc_roinames] = ...
    mmil_extract_aparc_tseries(parms.subj,parms.fname_in,args{:});
  ntpoints = size(aparc_data,1);
end;

if parms.aseg_flag
  args = mmil_parms2args(parms,aseg_tags);
  [aseg_data,aseg_roinames] = ...
    mmil_extract_aseg_tseries(parms.subj,parms.fname_in,args{:});
  ntpoints = size(aseg_data,1);
end;

% create row labels (time)
row_labels = cell(1,ntpoints);
for t=1:ntpoints
  row_labels{t} = sprintf('%0.2f',parms.TR*(t-1));
end;

if parms.aparc_flag && parms.aseg_flag
  % combine lh and rh aparc and aseg data and create column labels
  data = cat(2,aparc_data,aseg_data);
  col_labels = cat(2,aparc_roinames,aseg_roinames);
elseif parms.aparc_flag
  data = aparc_data;
  col_labels = aparc_roinames;
else
  data = aseg_data;
  col_labels = aseg_roinames;
end;

% write output file
mmil_write_csv(parms.fname_out,data,'col_labels',col_labels,...
  'row_labels',row_labels','firstcol_label','time (sec)');
