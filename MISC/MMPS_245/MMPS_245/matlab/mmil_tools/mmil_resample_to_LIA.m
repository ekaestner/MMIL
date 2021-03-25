function M_LIA_to_orig = mmil_resample_to_LIA(fname_in,fname_out,varargin)
%function M_LIA_to_orig = mmil_resample_to_LIA(fname_in,fname_out,[options])
%
% Purpose: resamples a volume to standard space (LIA, 256^3)
%
% Required:
%   fname_in : input file name (mgh format)
%   fname_out : output file name (mgh format)
%
% Optional:
%   'M_reg' : registration matrix
%     {default = eye(4)}
%   'interpm' : interpolation method
%       0: Nearest      1:Linear         2:cubic
%       3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%     {default  = 2}
%   'bclamp' : [0|1] set negative value to zero
%     {default = 1}
%   'forceflag' : [0|1] overwrite output
%     {default = 0}
%
% Output:
%   M_LIA_to_orig : voxel to voxel registration matrix
%     accounting for change in orientation and translation of center
%
% NOTE: if multi-frame volume, only resamples first frame
%
% Created:   04/23/10 by Don Hagler
% Last Mod:  12/15/11 by Vijay Venkatraman
%

if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'M_reg',eye(4),[],...
  'interpm',2,[0:5],...
  'bclamp',true,[false true],...
  'forceflag',false,[false true],...
});

if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

nvox_res = [256 256 256];
M_res =  [-1     0     0   129;...
           0     0     1  -129;...
           0    -1     0   129;...
           0     0     0     1];

if ~exist(fname_out,'file') || parms.forceflag
  [vol,M,tmp,nvox] = fs_load_mgh(fname_in,[],1); % first frame only in case multi-frame
  nvox = nvox(1:3);
  vol = ctx_mgh2ctx(vol,M);
  vol_res = ctx_mgh2ctx(zeros(nvox_res),M_res);
  vol_res = vol_resample_pad(vol,vol_res,parms.M_reg,parms.interpm,parms.bclamp);
  ctx_save_mgh(vol_res,fname_out);
  M_lph = vol.Mvxl2lph;
  M_res_lph = vol_res.Mvxl2lph;
else
  M = fs_read_header(fname_in);
  M_lph = M_RAS_TO_LPH*M;
  M_res_lph = M_RAS_TO_LPH*M_res;
end;

M_LIA_to_orig = inv(M_lph)*parms.M_reg*M_res_lph;

