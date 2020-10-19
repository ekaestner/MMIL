function fs_write_regdat(fname,varargin)
%function fs_write_regdat(fname,varargin)
%
% Usage:
%   fs_write_regdat(fname,'key1', value1,...);
%
% Required Input:
%   fname: full path or relative file name of output register.dat file
%
% Optional Input:
%   'M': 4x4 matrix specifying transformation from reference to registered image
%      {default = eye(4)}
%   'subj': subject name string (freesurfer recon directory name)
%      {default = 'unknown'}
%   'inplane': in plane voxel size
%      {default = 1}
%   'slicethick': slice thickness
%      {default = 1}
%   'ras2tk_flag': [0|1] whether to convert to tkregister transformation
%     matrix from ras2ras matrix
%     If 1, inplane and slicethick will be calculated from M_ref and M_reg
%      {default = 0}
%   'M_ref': vox2ras matrix for reference image
%      Only used if ras2tk_flag = 1
%      {default = eye(4)}
%   'M_reg': vox2ras matrix for registered image
%      Only used if ras2tk_flag = 1
%      {default = eye(4)}
%   'nvox_ref': vector of number of voxels in x, y, and z dimensions
%      for reference image
%      Only used if ras2tk_flag = 1
%      {default = [256 256 256]}
%   'nvox_reg': vector of number of voxels in x, y, and z dimensions
%      for registered image
%      Only used if ras2tk_flag = 1
%      {default = [256 256 256]}
%   'forceflag': [0|1] whether to overwrite existing output file
%      {default = 0}
%   
%  Created:  07/17/08   by Don Hagler
%  Last Mod: 08/21/15   by Don Hagler
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'M',eye(4),[],...
  'subj','unknown',[],...
  'inplane',1,[],...
  'slicethick',1,[],...
  'ras2tk_flag',false,sort([false true]),...
  'M_ref',eye(4),[],...
  'M_reg',eye(4),[],...
  'nvox_ref',[256 256 256],[],...
  'nvox_reg',[256 256 256],[],...
  'forceflag',false,sort([false true]),...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(fname,'file') && ~parms.forceflag
  fprintf('%s: WARNING: NOT overwriting file %s\n',mfilename,fname);
  return;
end;

if parms.ras2tk_flag
  % modify transformation matrix for tkregister
  M_ref_to_reg = parms.M;

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

  % reverse x rotation (by flipping z)
  X = diag([1 1 -1 1]);
  M_ref_to_reg = X*M_ref_to_reg*X;

  % reverse translations
  Y = diag([-1 -1 -1]);
  M_ref_to_reg(1:3,4) = Y*M_ref_to_reg(1:3,4);

  % convert from ras2ras to tkras2tkras
  parms.M = M_reg_tk * inv(M_reg) * M_ref_to_reg * M_ref * inv(M_ref_tk);

  parms.inplane = voxres_reg(1);
  parms.slicethick = voxres_reg(3);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(fname,'wt');
if fid==-1
  error('failed to open file %s for writing',fname);
end;
fprintf(fid,'%s\n',parms.subj);
fprintf(fid,'%f\n',parms.inplane);
fprintf(fid,'%f\n',parms.slicethick);
fprintf(fid,'1\n');
fprintf(fid,'%f ',parms.M(1,:));
fprintf(fid,'\n');
fprintf(fid,'%f ',parms.M(2,:));
fprintf(fid,'\n');
fprintf(fid,'%f ',parms.M(3,:));
fprintf(fid,'\n');
fprintf(fid,'%f ',parms.M(4,:));
fprintf(fid,'\n');
fprintf(fid,'round\n');
fclose(fid);

