function [M,subj,inplane,slicethick] = fs_read_regdat(fname,varargin)
%function [M,subj,inplane,slicethick] = fs_read_regdat(fname,varargin)
%
% Usage:
%   fs_read_regdat(fname,'key1', value1,...);
%
% Required Input:
%   fname: full path or relative file name of input register.dat file
%
% Optional Input:
%   'tk2ras_flag': [0|1] whether to convert from tkregister transformation
%     matrix to ras2ras matrix
%     If 1, be sure M_ref, M_reg, nvox_ref, and nvox_reg are correctly specified
%      {default = 0}
%   'M_ref': vox2ras matrix for reference image
%      Only used if tk2ras_flag = 1
%      {default = eye(4)}
%   'M_reg': vox2ras matrix for registered image
%      Only used if tk2ras_flag = 1
%      {default = eye(4)}
%   'nvox_ref': vector of number of voxels in x, y, and z dimensions
%      for reference image
%      Only used if tk2ras_flag = 1
%      {default = [256 256 256]}
%   'nvox_reg': vector of number of voxels in x, y, and z dimensions
%      for registered image
%      Only used if tk2ras_flag = 1
%      {default = [256 256 256]}
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
%   9.969735e-01 -7.713004e-02 -9.725398e-03 1.009433e+00  # transformation matrix
%   2.020466e-02 3.778771e-01 -9.256357e-01 2.806142e+01
%   7.506926e-02 9.226375e-01 3.782920e-01 9.858908e+00
%   0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00
%   round           # how to handle float to int conversion
%
%  See fs_write_regdat.m
%
%  Created:  07/18/08   by Don Hagler
%  Last Mod: 08/21/15   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = []; subj = []; inplane = []; slicethick = [];
if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms( varargin, {...
  'tk2ras_flag',false,sort([false true]),...
  'M_ref',eye(4),[],...
  'M_reg',eye(4),[],...
  'nvox_ref',[256 256 256],[],...
  'nvox_reg',[256 256 256],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read register.dat file
if ~exist(fname,'file')
  error('file %s not found',fname);
end;
fid=fopen(fname,'rt');
if fid==-1
  error('failed to open file %s for writing',fname);
end;
subj = fscanf(fid,'%s',1);
inplane = fscanf(fid,'%f',1);
slicethick = fscanf(fid,'%f',1);
tmp = fscanf(fid,'%f',1);
M = fscanf(fid,'%f',[4,4])';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.tk2ras_flag
  % modify transformation matrix from tkregister to ras2ras

  M_ref = parms.M_ref;
  voxres_ref = sqrt(sum(M_ref(1:3,1:3).^2,1));
  nx = parms.nvox_ref(1);
  ny = parms.nvox_ref(2);
  nz = parms.nvox_ref(3);
  ps = voxres_ref(1);
  st = voxres_ref(3);
  M_ref_tk = [-ps 0 0 (nx/2+1)*ps; 0 0 st -(nz/2+1)*st; 0 -ps 0 (ny/2+1)*ps; 0 0 0 1];

  M_reg = parms.M_reg;
  voxres_reg = sqrt(sum(M_reg(1:3,1:3).^2,1));
  nx = parms.nvox_reg(1);
  ny = parms.nvox_reg(2);
  nz = parms.nvox_reg(3);
  ps = voxres_reg(1);
  st = voxres_reg(3);
  M_reg_tk = [-ps 0 0 (nx/2+1)*ps; 0 0 st -(nz/2+1)*st; 0 -ps 0 (ny/2+1)*ps; 0 0 0 1];

  % convert from tkras2tkras to ras2ras
  M_ref_to_reg = M_reg * inv(M_reg_tk) * M * M_ref_tk * inv(M_ref);

  % reverse translations
  Y = diag([-1 -1 -1]);
  M_ref_to_reg(1:3,4) = Y*M_ref_to_reg(1:3,4);

  % reverse x rotation (by flipping z)
  X = diag([1 1 -1 1]);
  M_ref_to_reg = X*M_ref_to_reg*X;

  M = M_ref_to_reg;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

