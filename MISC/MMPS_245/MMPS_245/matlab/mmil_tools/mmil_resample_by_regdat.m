function [M,subj,inplane,slicethick] = mmil_resample_by_regdat(fname,varargin)
%function [M,subj,inplane,slicethick] = mmil_resample_by_regdat(fname,varargin)
%
% Usage:
%   mmil_resample_by_regdat(fname,'key1', value1,...);
%
% Required Input:
%   fname: full path or relative file name of input volume file (mgh format)
%     to be registered to fname_ref
%  'fname_regdat': full path or relative file name of input register.dat file
% 
% Optional Input:
%  'fname_ref': full path of reference volume file (mgh format)
%     If empty will assume identity vox2ras matrix and
%     256 cubed voxels
%     {default = []}
%  'fname_out': full path of output registered volume file (mgh format)
%     If empty will generate a name based on input fname
%     {default = []}
%  'frames': selected frames of input fname to resample
%     If empty will resample all frames
%     {default = []}
%  'interpm': interpolation method
%     0: nearest neighbor, 1: linear, 2: cubic,
%     3: Key's spline, 4: cubic spline, 5: Hamming sinc
%     {default = 2}
%  'bclamp': [0|1] set negative values to zero
%     {default = 0}
%  'forceflag': [0|1] force overwrite of fname_out
%     {default = 0}
%
% Output:
%   M: 4x4 matrix specifying transformation from reference to registered image
%   subj: subject name string (freesurfer recon directory name)
%   inplane: in-plane voxel size
%   slicethick: slice thickness
%
%
% Note:  A FreeSurfer register.dat file looks like this
% (minus leading spaces and # comments):
%   bert           # subject name
%   1.875000       # in-plane voxel size
%   2.500000       # slice thickness
%   1.000000       # ?
%   9.969735e-01 -7.713004e-02 -9.725398e-03 1.009433e+00
%   2.020466e-02 3.778771e-01 -9.256357e-01 2.806142e+01
%   7.506926e-02 9.226375e-01 3.782920e-01 9.858908e+00
%   0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00
%   round           # how to handle float to int conversion
%
%  See fs_write_regdat.m, fs_read_regdat.m
%
%  created:       07/18/08   by Don Hagler
%  last modified: 09/23/08   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = []; subj = []; inplane = []; slicethick = [];
if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms( varargin, {...
  'fname_regdat',[],[],...
  'fname_ref',[],[],...
  'fname_out',[],[],...
  'frames',[],[],...
  'interpm',2,[1,2,3,4],...
  'bclamp',false,[false true],...
  'forceflag',false,[false true],...
});

if isempty(parms.fname_regdat)
  error('must specify fname_regdat');
end;

if ~isempty(parms.frames) && min(parms.frames)<1
  error('frame numbers must be > 1');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fname_out)
  [fpath,fstem,fext] = fileparts(fname);
  parms.fname_out = sprintf('%s/%s_res.mgh',fpath,fstem); % todo: allow infix as parameter?
end;

if ~exist(parms.fname_out,'file') || parms.forceflag
  if ~isempty(parms.fname_ref)
    [vol_ref,M_ref,mrparms,nvox_ref] = fs_load_mgh(parms.fname_ref,[],[],1);
  else
    M_ref = eye(4);
    nvox_ref = [256 256 256 1];
  end;  

  [vol_reg,M_reg,mrparms,nvox_reg] = fs_load_mgh(fname,[],[],1);
  nframes = nvox_reg(4);
  if isempty(parms.frames)
    parms.frames = [1:nframes];
  else
    if max(parms.frames)>nframes
      error('only %d frames in file %s',nframes,fname);
    end;
    nframes = length(parms.frames);
  end;

  [M_ref_to_reg,subj,inplane,slicethick] = fs_read_regdat(parms.fname_regdat,...
    'tk2ras_flag',1,...
    'M_ref',M_ref,...
    'M_reg',M_reg,...
    'nvox_ref',nvox_ref,...
    'nvox_reg',nvox_reg);

  vol_ref = ctx_mgh2ctx(zeros(nvox_ref(1:3)),M_ref);

  % loop over multiple frames
  fprintf('%s: resampling reg to ref...\n',mfilename);
  vol_res = zeros([nvox_ref(1:3),nframes],'single');
  for f=parms.frames
    vol = fs_load_mgh(fname,[],f,0,1);
    vol_reg = ctx_mgh2ctx(vol,M_reg);
    vol_res_tmp = vol_resample_pad(vol_reg,vol_ref,M_ref_to_reg,...
                                   parms.interpm,parms.bclamp);
    vol_res(:,:,:,f) = vol_res_tmp.imgs;
  end;
  fs_save_mgh(vol_res,parms.fname_out,M_ref);

  M = M_ref_to_reg;
end;

