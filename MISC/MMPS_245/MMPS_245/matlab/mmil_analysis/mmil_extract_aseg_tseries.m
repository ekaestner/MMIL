function [data,roinames] = mmil_extract_aseg_tseries(subj,fname,varargin)
%function [data,roinames] = mmil_extract_aseg_tseries(subj,fname,[options])
%
% Purpose:
%   extract time series for each ROI in aseg
%
% Required parameters:
%  subj:  string specifying the subject name
%  fname: full or relative path of 4D functional volume
%    (must be mgh/mgz format)
%
% Optional parameters:
%  'outdir' : output directory for aseg eroded and resampled
%      to resolution of fname
%    If empty, will attempt to write to path containing fname
%    {default = []}
%  'fname_aseg': name of aseg file
%     if empty, will use aseg.mgz in subjdir/subj/mri
%    {default = []}
%  'skipTRs': number of frames to remove from beginning of fname
%    {default = 0}
%  'regfile': register.dat file containing 4x4 registration matrix
%    If not supplied, fname should already be resampled to structural space
%    {default = []}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'erode_flag': [0|1] "erode" ROIs by smoothing and thresholding
%    {default = 0}
%  'erode_nvoxels': number of voxels to erode
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Created:  03/06/12 by Don Hagler
% Prev Mod: 10/29/12 by Don Hagler
% Last Mod: 06/30/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'subj',subj,[],...
  'fname_in',fname,[],...
  'outdir',[],[],...
  'fname_aseg',[],[],...
  'skipTRs',0,[0,1000],...
  'regfile',[],[],...
  'subjdir',[],[],...
  'erode_flag',false,[false true],...
  'erode_nvoxels',1,[1:100],...
  'forceflag',false,[false true],...
...
  'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79],[1,Inf],...
});

if parms.erode_flag
  parms.erode_outfix = 'erode';
  if parms.erode_nvoxels>1
    parms.erode_outfix = sprintf('%s%d',parms.erode_outfix,parms.erode_nvoxels);
  end;
end;

%% todo: use aseg_roilist from function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

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

% check fname_aseg
if isempty(parms.fname_aseg)
  parms.fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
end;
if ~exist(parms.fname_aseg,'file')
  error('aseg file %s not found',parms.fname_aseg);
end;

% set outdir and outstem
[tmp_outdir,tmp_outstem,tmp] = fileparts(parms.fname_in);
if isempty(parms.outdir), parms.outdir = tmp_outdir; end;
parms.outstem = [parms.outdir '/' tmp_outstem];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

%[tmp,volsz] = fs_read_header(parms.fname_in);
[M,volsz] = mmil_load_mgh_info(parms.fname_in,parms.forceflag,parms.outdir);
nframes = volsz(4);
if parms.skipTRs>0
  frames = setdiff([1:nframes],[1:parms.skipTRs]);
else
  frames = [1:nframes];
end;
ntpoints = length(frames);

% resample aseg to BOLD resolution
if ~isempty(parms.regfile)
  fname_aseg_res = sprintf('%s_aseg.mgz',parms.outstem);
  fs_resample_aseg(parms.fname_aseg,parms.regfile,parms.fname_in,...
    fname_aseg_res,parms.forceflag);
  parms.fname_aseg = fname_aseg_res;
end;

% erode aseg
if parms.erode_flag
  fname_aseg_erode = sprintf('%s_aseg_%s.mgz',...
    parms.outstem,parms.erode_outfix);
  fs_erode_aseg(parms.fname_aseg,fname_aseg_erode,...
    'nvoxels',parms.erode_nvoxels,'forceflag',parms.forceflag);
  parms.fname_aseg = fname_aseg_erode;
end;

% extract time series for each aseg ROI
aseg_results = mmil_aseg_roi(parms.fname_in,parms.fname_aseg,...
  'roilist',parms.aseg_roilist,'frames',frames);
nroi = length(aseg_results);

% combine aseg data
roinames = cell(1,nroi);
data = zeros(ntpoints,nroi);
for r=1:nroi
  roinames{r} = aseg_results(r).roiname;
  data(:,r) = aseg_results(r).avg;
end;

return;

