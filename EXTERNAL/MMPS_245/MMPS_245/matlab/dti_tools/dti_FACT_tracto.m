function dti_FACT_tracto(varargin)
%function dti_FACT_tracto([options])
%
% Purpose:
%  Create fiber streamlines using FACT tractography
%     from Susumo Mori's DTI Studio
%
% Usage:
%  dti_FACT_tracto('key1', value1,...);
%
%  Must supply either fname_FA and fname_V0, or vol, M, qmat, and bvals
%
% Optional Parameters for specifying FA and V0:
%   'fname_FA': full or relative path of file containing FA volume
%     {default = []}
%   'fname_V0': full or relative path of file containing V0 volume
%     (1st eigen vector from diffusion tensor fit, with x, y, and z components)
%     {default = []}
%
% Optional Parameters for specifying diffusion MRI data:
%   'vol': 4D volume containing multiple diffusion weighted volumes
%     {default = []}
%   'M': 4x4 vox2ras matrix for vol
%     {default = []}
%   'qmat': matrix of diffusion direction vectors
%     {default = []}
%   'bvals': vector of b values
%     one for all, or one for each diffusion direction
%     {default = []}
%
% Optional Parameters:
%   'outdir': full or relative path of output directory
%     {default = 'FACT_Tracto'}
%   'outstem': output file stem
%     {default = []}
%   'orient': reorient data to specified slice orientation
%     (e.g. 'LAS', 'LPI', etc.)
%     {default = []}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Optional Parameters for FACT tractography:
%  'thresh_FA': FA threshold for starting and ending tracking
%    {default = 0.15}
%  'thresh_angle': turning angle threshold
%    {default = 50}
%  'min_fiberlen': minimm fiber length
%    {default = 20}
%
% Created:  05/16/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,4), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_FA',[],[],...
  'fname_V0',[],[],...
...
  'vol',[],[],...
  'M',[],[],...
  'qmat',[],[],...
  'bvals',[],[],...
...
  'outdir','FACT_Tracto',[],...
  'outstem',[],[],...
  'orient',[],[],...
  'forceflag',false,[false true],...
...
  'thresh_FA',0.15,[],...
  'thresh_angle',50,[],...
  'min_fiberlen',20,[],...
...
  'measlist',{'FA','V0'},[],...
  'orient_ref','LPS',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if ~isempty(parms.fname_FA)
  if ~exist(parms.fname_FA,'file')
    error('file %s not found',parms.fname_FA);
  end;
  if isempty(parms.fname_V0)
    error('fname_V0 not specified');
  end;
  if ~exist(parms.fname_V0,'file')
    error('file %s not found',parms.fname_V0);
  end;
elseif ~isempty(parms.vol)
  if isempty(parms.M), error('M not specified'); end;
  if isempty(parms.qmat), error('qmat not specified'); end;
  if isempty(parms.bvals), error('bvals not specified'); end;
  [nx,ny,nz,nf] = size(parms.vol);
  if length(parms.bvals)==1
    parms.bvals = parms.bvals*ones(nf,1);
  end;
  if length(parms.bvals)~=nf
    error('number of bvals must match frames in vol');
  end;
  parms.bvals = mmil_colvec(parms.bvals);
  if size(parms.qmat,1)~=nf
    error('number of directions in qmat must match frames in vol');
  end;
  if size(parms.qmat,2)~=3
    error('qmat must have 3 columns');
  end;
else
  error('must specify fname_FA and fname_V0 or vol, M, qmat, and bvals');
end;

if mmil_isrelative(parms.outdir)
  parms.outdir = [pwd '/' parms.outdir];
end;
if ~isempty(parms.outstem)
  parms.outstem = [parms.outstem '_'];
end;
parms.outstem = sprintf('%s/%s',parms.outdir,parms.outstem);

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fname_FA)
  % prepare data, qmat, and bvals
  % save as mgh, reorient if necessary
  fname_data_mgh = sprintf('%sdata.mgz',parms.outstem);
  if ~exist(fname_data_mgh,'file') || parms.forceflag
    % reorient data
    if ~isempty(parms.orient)
      [parms.vol,parms.M] = fs_reorient(parms.vol,parms.M,parms.orient);
    end;
    % save to mgh file
    fs_save_mgh(parms.vol,fname_data_mgh,parms.M);
  end;
  % fit data with tensor
  fname_DTfit = sprintf('%sDTfit.mat',parms.outstem);
  made_flag = 0;
  if ~exist(fname_DTfit,'file') || parms.forceflag
    if ~isempty(parms.orient)
      [parms.vol,parms.M] = fs_load_mgh(fname_data_mgh);
    end;
    DTfit = dti_fit_tensor(parms.vol,parms.qmat,...
      'M',parms.M,'bvals',parms.bvals);
    save(fname_DTfit,'DTfit');
    made_flag = 1;
  end;
  if ~made_flag, load(fname_DTfit); end;
  % calculate diffusion tensor measures
  fname_DTmeas = sprintf('%sDTmeas.mat',parms.outstem);
  made_flag = 0;
  if ~exist(fname_DTmeas,'file') || parms.forceflag
    DTmeas = dti_calc_DTmeas(DTfit);
    save(fname_DTmeas,'DTmeas');
    made_flag = 1;
  end;
  if ~made_flag, load(fname_DTmeas); end;
  % convert diffusion tensor measures to mgz files
  fstem = sprintf('%sDT',parms.outstem);
  dti_convert_DTmeas(DTmeas,fstem,'M',DTfit.M,...
    'mgz_flag',0,'measlist',parms.measlist,'forceflag',parms.forceflag);
  parms.fname_FA = [fstem '_FA.mgh'];
  parms.fname_V0 = [fstem '_V0.mgh'];
  parms.M = DTfit.M;
  parms.volsz = DTfit.volsz;
else
  % make local copy of FA and V0, reorient if desired
  fname_FA = sprintf('%sDT_FA.mgh',parms.outstem);
  fname_V0 = sprintf('%sDT_V0.mgh',parms.outstem);  
  if ~isempty(parms.orient)
    if ~exist(fname_FA,'file') ||...
       ~exist(fname_V0,'file') || parms.forceflag
      [vol,M] = fs_load_mgh(fname_V0);
      [vol,M] = fs_reorient(vol,M,parms.orient);
      fs_save_mgh(vol,fname_V0,M);      
      [vol,M] = fs_load_mgh(fname_FA);
      [vol,M] = fs_reorient(vol,M,parms.orient);
      fs_save_mgh(vol,fname_FA,M);
    end;
    parms.fname_FA = fname_FA;
    parms.fname_V0 = fname_V0;
  else
    fs_copy_mgh(parms.fname_FA,fname_FA);
    fs_copy_mgh(parms.fname_V0,fname_V0);  
  end;
  parms.fname_FA = fname_FA;
  parms.fname_V0 = fname_V0;
  [parms.M,parms.volsz] = fs_read_header(parms.fname_FA);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run tractography

if ~isempty(parms.orient)
  orient = parms.orient;
else
  orient = fs_read_orient([],parms.M);
end;
[permvec,flipvec] = fs_compare_orient(orient,parms.orient_ref);
flipflags = flipvec<0;

fname_path = sprintf('%sfiber_path.grp',parms.outstem);
if ~exist(fname_path,'file') || parms.forceflag
  dti_track_fibers(parms.fname_FA,parms.fname_V0,...
    'flipx_flag',flipflags(1),...
    'flipy_flag',flipflags(2),...
    'flipz_flag',flipflags(3),...
    'thresh_FA',parms.thresh_FA,...
    'thresh_angle',parms.thresh_angle,...
    'min_fiberlen',parms.min_fiberlen,...
    'outstem',[parms.outstem 'fiber_path']);
end;

fname_dat = sprintf('%sfiber_count.dat',parms.outstem);
if ~exist(fname_dat,'file') || parms.forceflag
  dti_fiber_path_to_mask(fname_path,fname_dat,1);
end;

fname_mgh = sprintf('%sfiber_count.mgh',parms.outstem);
if ~exist(fname_mgh,'file') || parms.forceflag
  dti_dat2mgh(fname_dat,parms.volsz(1:3),parms.M,[],flipflags(3));
end;

