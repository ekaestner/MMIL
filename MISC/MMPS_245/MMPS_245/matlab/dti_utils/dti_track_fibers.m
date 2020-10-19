function dti_track_fibers(fname_FA,fname_V0,varargin)
%function dti_track_fibers(fname_FA,fname_V0,[options])
%
% Purpose: runs fiber tracking program for DTI data
%
% Usage:
%  dti_track_fibers(fname_FA,fname_V0,'key1', value1,...);
%
% Required Input:
%   fname_FA: full or relative path of file containing FA volume
%   fname_V0: full or relative path of file containing V0 volume
%     (1st eigen vector from diffusion tensor fit, with x, y, and z components)
%
% Optional Parameters:
%  'rois': structure defining rois to be used in fiber tracking
%     example:
%     rois(1).fname = 'CST_L_RoiMap_00.dat'
%     rois(1).logic = 'OR'
%     rois(2).fname = 'CST_L_RoiMap_01.dat'
%     rois(2).logic = 'AND'
%     rois(3).fname = 'CST_L_RoiMap_02.dat'
%     rois(3).logic = 'AND'
%     rois(4).fname = 'CST_L_RoiMap_03.dat'
%     rois(4).logic = 'NOT'
%    If empty, will use fname_prob or entire volume as OR
%     {default = []}
%  'fname_prob' - full or relative path of file containing
%     probability map -- values must be between and including 0 and 1
%     {default = []}
%  'prob_exponent' - exponent applied to probability map
%     must be non-negative
%     zero exponent will binarize probability map (0 or 1)
%     {default = 1}
%  'outstem'  - output file stem for fiber track files
%    {default = 'fibers'}
%  'roidir' - directory containing roi files
%    if empty or ommitted, assumes that rois contains full path names
%    {default = []}
%  'thresh_prob' - probability threshold for starting and ending tracking
%    {default = 0.25}
%  'thresh_FA' - FA threshold for starting and ending tracking
%    {default = 0.15}
%  'thresh_angle' - turning angle threshold
%    {default = 50}
%  'min_fiberlen' - minimm fiber length
%    {default = 7}
%  'flipx_flag' - [0|1] toggle flip (make negative) x component of eigen vector
%    {default = 0}
%  'flipy_flag' - [0|1] toggle flip (make negative) y component of eigen vector
%    {default = 0}
%  'flipz_flag' - [0|1] toggle flip (make negative) z component of eigen vector
%    {default = 0}
%  'verbose': [0|1] display status messages and fact_suabe output
%    {default = 0}
%
% Optional Parameters Relevant with raw input (e.g. dat format)
%  'imgorient' - image orientation (0=coronal, 1=axial, 2=sagittal)
%    {default = 1}
%  'imgseq' - image sequence (slice order) (0=positive, 1=negative)
%    {default = 0}
%  'volsz' - vector of dimension sizes (e.g. [256,256,256])
%    {default = [256,256,256]}
%  'voxsz' - vector of voxel sizes (e.g. [1,1,1])
%    {default = [1,1,1]}
%
% Created:  11/05/06 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'rois',[],[],...
  'fname_prob',[],[],...
  'prob_exponent',1,[],...
  'outstem','fibers',[],...
  'roidir',[],[],...
  'thresh_prob',0.25,[],...
  'thresh_FA',0.15,[],...
  'thresh_angle',50,[],...
  'min_fiberlen',7,[],...
  'flipx_flag',false,[false true],...
  'flipy_flag',false,[false true],...
  'flipz_flag',false,[false true],...
  'verbose',false,[false true],...
  'imgorient',1,[0:2],...
  'imgseq',0,[0,1,1],...
  'volsz',[256,256,256],[],...
  'voxsz',[1,1,1],[],...
...
  'outmask_flag',false,[false true],... % output masks
  'fact_thresh_FA',0.0001,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

if ~mmil_check_nargs(nargin,3), return; end;

smf = 10^-5;

[tmp,tmp_fstem] = fileparts(tempname);
tmp_fstem = [parms.outstem '_' tmp_fstem];

% check that input files exist, check formats
if ~exist(fname_FA,'file')
  error('FA file %s not found',fname_FA);
end;
if ~exist(fname_V0,'file')
  error('V0 file %s not found',fname_V0);
end;
if ~isempty(parms.fname_prob) & ~exist(parms.fname_prob,'file')
  error('prob file %s not found',parms.fname_prob);
end;

if parms.prob_exponent < 0
  error('parms.prob_exponent must be non-negative (is %f)',...
    parms.prob_exponent);
end;

% check input formats
raw_input_flag = 0;
cleanup_flag = 0;

if ~isempty(parms.fname_prob)
  [fpath,fstem,fext] = fileparts(parms.fname_prob);
  switch fext
    case {'.mgh','.mgz'}
      vol_probmap = fs_load_mgh(parms.fname_prob);
    case '.mat'
      vol_probmap = mmil_load_sparse(parms.fname_prob);
    otherwise
      vol_probmap = dti_load_dat(parms.fname_prob,parms.volsz);
      if isempty(vol_probmap)
        fprintf('%s: WARNING: retrying to read probability map from %s as uint8\n',...
          mfilename,parms.fname_prob);
        vol_probmap = dti_load_dat(parms.fname_prob,parms.volsz,[],[],'uint8');
      end;
      if isempty(vol_probmap)
        fprintf('%s: ERROR: failed to read probability map from %s\n',...
          mfilename,parms.fname_prob);
        return;
      end;
      raw_input_flag = 1;
  end;
else
  vol_probmap = [];
end;

% check max val in vol_probmap
if ~isempty(vol_probmap)
  maxval = max(vol_probmap(:));
  minval = min(vol_probmap(:));
  if maxval > 1
    fprintf('%s: WARNING: probability map contains values > 1\n',mfilename);
    vol_probmap = vol_probmap / (maxval + eps);
  elseif minval < 0
    fprintf('%s: WARNING: probability map contains values < 0\n',mfilename);
    vol_probmap(vol_probmap<0)=0;
  end;

  % apply parms.thresh_prob to vol_probmap
  vol_probmap(vol_probmap<parms.thresh_prob)=0;

  % apply parms.prob_exponent
  if parms.prob_exponent == 0
    vol_probmap(vol_probmap>smf)=1;
  else
    vol_probmap = vol_probmap.^parms.prob_exponent;
  end;
end;

% determine format of FA file
[fpath,fstem,fext] = fileparts(fname_FA);
vol_FA = [];
switch fext
  case {'.mgh','.mgz','.mat'}
    if ~isempty(vol_probmap) || ~strcmp(fext,'.mgh')
      if raw_input_flag
        error('file type mismatch between probmask and FA');
      end;
      if strcmp(fext,'.mat')
        [vol_FA,M] = mmil_load_sparse(fname_FA);
      else
        [vol_FA,M] = fs_load_mgh(fname_FA);
      end;
    end;
  otherwise
    if isempty(vol_probmap)
      raw_input_flag = 1;
    elseif ~raw_input_flag
      error('file type mismatch between probmask and FA');
    else
      vol_FA = dti_load_dat(fname_FA,parms.volsz);
    end;
end;

% apply probability map to FA if indicated
if ~isempty(vol_FA)
  % apply parms.thresh_FA to vol_FA
  vol_FA(vol_FA<parms.thresh_FA)=0;
  if ~isempty(vol_probmap)
    if any(size(vol_FA)~=size(vol_probmap))
      error('volume size mismatch between probmask and FA');
    end;
    % combine prob and FA
    vol_FA = vol_FA.*vol_probmap;
  end;
  if raw_input_flag
    fname_FA = sprintf('%s_FA.dat',tmp_fstem);
    dti_save_dat(vol_FA,fname_FA);
  else
    fname_FA = sprintf('%s_FA.mgh',tmp_fstem);
    fs_save_mgh(vol_FA,fname_FA,M);
  end;
  cleanup_flag = 1;
else % otherwise, let fact apply FA threshold
  parms.fact_thresh_FA = parms.thresh_FA;
end;

[fpath,fstem,fext] = fileparts(fname_V0);
switch fext
  case {'.mgh','.mgz','.mat'}
    if raw_input_flag
      error('file type mismatch between FA and V0');
    end;
    if strcmp(fext,'.mgz')
      [vol,M] = fs_load_mgh(fname_V0);
      fname_V0 = sprintf('%s_V0.mgh',tmp_fstem);
      fs_save_mgh(vol,fname_V0,M);
      cleanup_flag = 1;
    elseif strcmp(fext,'.mat')
      fname_V0_mgh = sprintf('%s_V0.mgh',tmp_fstem);
      fname_V0 = mmil_sparse2mgh(fname_V0,fname_V0_mgh,1);
      cleanup_flag = 1;
    end;
  otherwise
    if ~raw_input_flag
      error('file type mismatch between FA and V0');
    end;
end;

if isempty(parms.rois)
  parms.rois.fname = fname_FA; % entire volume unless fname_prob supplied
  parms.rois.logic = 'OR';
end;

for r=1:length(parms.rois)
  parms.rois(r).raw_flag = 0;
  if ~isempty(parms.roidir)
    parms.rois(r).fname = sprintf('%s/%s',parms.roidir,parms.rois(r).fname);
  end;
  if ~exist(parms.rois(r).fname,'file')
    error('ROI file %s not found',parms.rois(r).fname);
  else
    [fpath,fstem,fext] = fileparts(parms.rois(r).fname);
    switch fext
      case {'.mgh','.mgz'}
        if strcmp(fext,'.mgz')
          [vol,M] = fs_load_mgh(parms.rois(r).fname);
          parms.rois(r).fname = sprintf('%s_ROI%d.mgh',tmp_fstem,r);
          fs_save_mgh(vol,parms.rois(r).fname,M);
        end;
      case '.mat'
        [vol,M] = mmil_load_sparse(parms.rois(r).fname);
        parms.rois(r).fname = sprintf('%s_ROI%d.mgh',tmp_fstem,r);
        fs_save_mgh(vol,parms.rois(r).fname,M);
      otherwise
        parms.rois(r).raw_flag = 1;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run fiber tracking program

% create cmd
cmd = 'fact_suabe';
if raw_input_flag
  cmd = sprintf('%s -iFAraw %s -iV0raw %s \\\n',cmd,fname_FA,fname_V0);
  cmd = sprintf('%s -iImgDims %d %d %d \\\n',cmd,parms.volsz(1),parms.volsz(2),parms.volsz(3));
  cmd = sprintf('%s -iVoxSize %0.5f %0.5f \\\n',cmd,parms.voxsz(1),parms.voxsz(2));
  cmd = sprintf('%s -iSliceThickness %0.5f \\\n',cmd,parms.voxsz(3));
  cmd = sprintf('%s -iImageOrientation %d \\\n',cmd,parms.imgorient);
%  cmd = sprintf('%s -iImageSequence %d' \\\n,cmd,parms.imgseq);
else
  cmd = sprintf('%s -iFAmgh %s -iV0mgh %s \\\n',cmd,fname_FA,fname_V0);
end;
cmd = sprintf('%s -iStartFA %0.6f \\\n',cmd,parms.fact_thresh_FA);
cmd = sprintf('%s -iStopFA %0.6f \\\n',cmd,parms.fact_thresh_FA);
cmd = sprintf('%s -iTurnAngleDeg %0.2f \\\n',cmd,parms.thresh_angle);
cmd = sprintf('%s -iFiberLenMin %d \\\n',cmd,parms.min_fiberlen);
cmd = sprintf('%s -iFlipEigenX %d -iFlipEigenY %d -iFlipEigenZ %d \\\n',...
  cmd,parms.flipx_flag,parms.flipy_flag,parms.flipz_flag);
for r=1:length(parms.rois)
  logic = upper(parms.rois(r).logic);
  switch logic
    case {'OR','AND','NOT'}
    case 'CUT0'
      logic = 'CUT_OR';
    case 'CUT1'
      logic = 'CUT_AND';
    otherwise
      fprintf('%s: WARNING: unsupported ROI logic for ROI %d: %s\n',...
        mfilename,r,logic);
      continue;
  end;
  if parms.rois(r).raw_flag
    cmd = sprintf('%s -iRoiraw%s %s \\\n',...
      cmd,upper(parms.rois(r).logic),parms.rois(r).fname);
  else
    cmd = sprintf('%s -iRoimgh%s %s \\\n',...
      cmd,upper(parms.rois(r).logic),parms.rois(r).fname);
  end;
end;
if parms.outmask_flag
  if raw_input_flag
    cmd = sprintf('%s -oSelFiberVolFile %s_mask.dat \\\n',cmd,parms.outstem);
  else
    cmd = sprintf('%s -oSelFiberMghFile %s_mask.mgh \\\n',cmd,parms.outstem);
  end;
end;
cmd = sprintf('%s -oSelFiberFile %s.grp\n',cmd,parms.outstem);

if parms.verbose
  fprintf('%s: cmd = %s\n',mfilename,cmd);
end;

[status,result] = unix(cmd);
if status
  error('fiber tracking failed with these messages:\n%s',result);
elseif parms.verbose
  fprintf('%s: fiber tracking successful with these messages:\n',mfilename);
  disp(result);
end;

%return;

% clean up tmp files
if cleanup_flag
  cmd = sprintf('rm %s*',tmp_fstem);
  [status,result]=unix(cmd);
  if status
    error('failed to remove temp files:\n%s',result);
  end;
end;

