function data = mmil_extract_aseg_group_tseries(subj,fname,varargin);
%function data = mmil_extract_aseg_group_tseries(subj,fname,[options]);
%
% Purpose:
%   extract time series for group of aseg ROI codes
%
% Required parameters:
%  subj:  string specifying the subject name
%  fname: full or relative path of 4D functional volume
%    (must be mgh/mgz format)
%
% Optional parameters:
%  'outdir': output directory for aseg eroded and resampled
%      to resolution of fname
%     if empty, will attempt to write to path containing fname
%    {default = []}
%  'outstem': output file stem
%     if empty, will use file stem of fname
%    {default = []}
%  'fname_aseg': name of aseg file
%     if empty, will use aseg.mgz in subjdir/subj/mri
%    {default = []}
%  'aseg_codes': vector of aseg codes to be grouped together
%    note: [2,41] are left and right cerebral white matter
%    {default = [2,41]}
%  'mask_name': name of ROI group to be included in output file names
%    {default = 'wm_mask'}
%  'skipTRs': number of frames to remove from beginning of fname
%    {default = 0}
%  'regfile': register.dat file containing 4x4 registration matrix
%    If not supplied, fname should already be resampled to structural space
%    {default = []}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'erode_flag': [0|1] "erode" ROIs by smoothing and thresholding
%    {default = 1}
%  'erode_nvoxels': number of voxels to erode
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Output:
%   data: time series vector with size = [ntpoints,1]
%
% Created:  03/06/12 by Don Hagler
% Prev Mod: 05/02/13 by Don Hagler
% Last Mod: 06/30/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'subj',subj,[],...
  'fname_in',fname,[],...
...
  'outdir',[],[],...
  'outstem',[],[],...
  'fname_aseg',[],[],...
  'aseg_codes',[2,41],[1,Inf],...
  'mask_name','wm_mask',[],...
  'skipTRs',0,[0,1000],...
  'regfile',[],[],...
  'subjdir',[],[],...
  'erode_flag',true,[false true],...
  'erode_nvoxels',1,[1:100],...
  'forceflag',false,[false true],...
});

if parms.erode_flag
  parms.erode_outfix = 'erode';
  if parms.erode_nvoxels>1
    parms.erode_outfix = sprintf('%s%d',parms.erode_outfix,parms.erode_nvoxels);
  end;
end;

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
if isempty(parms.outdir) || isempty(parms.outstem)
  [tmp_outdir,tmp_outstem,tmp] = fileparts(parms.fname_in);
  if isempty(parms.outdir), parms.outdir = tmp_outdir; end;
  if isempty(parms.outstem), parms.outstem = tmp_outstem; end;
end;
if mmil_isrelative(parms.outstem)
  parms.outstem = [parms.outdir '/' parms.outstem];
end;

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

% create mask
fname_mask = sprintf('%s_%s.mgz',parms.outstem,parms.mask_name);
fs_aseg_mask(parms.fname_aseg,'fname_mask',fname_mask,...
  'aseg_codes',parms.aseg_codes,'forceflag',parms.forceflag);

% erode mask
if parms.erode_flag
  fname_mask_erode = sprintf('%s_%s_%s.mgz',...
    parms.outstem,parms.mask_name,parms.erode_outfix);
  fs_erode_mask(fname_mask,fname_mask_erode,...
    'nvoxels',parms.erode_nvoxels,'forceflag',parms.forceflag);
  fname_mask = fname_mask_erode;
end;

% extract time series for mask
results = mmil_multi_roi(parms.fname_in,fname_mask);

% combine aseg data
data = reshape(results(1).avg,[ntpoints,1]);

return;

