function fname_out = epi_gradunwarp(fname_in,varargin)
%function fname_out = epi_gradunwarp(fname_in,[options])
%
% Purpose: correction for distortion caused by gradient non-linearities
%
% Required Parameters:
%  fname_in: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Parameters:
%  'fname_out': output file name of motion corrected 4D volume
%    If empty, will append '_gruw' to file stem of fname_in
%    {default = [])
%  'gwtype': gradient type number
%     0:  Siemens Sonata
%     1:  Siemens Allegra
%     2:  GE BRM
%     3:  GE CRM
%     4:  Siemens Avanto
%     5:  Siemens AXXess/Espree
%     6:  Siemens Quantum/Symphony
%     7:  GE Twin Speed Whole Body
%     8:  GE Twin Speed Zoom
%     9:  GE mr450 or mr750
%     10: GE MR750W
%     11: Siemens Skyra
%     12: Siemens Connectome Skyra
%     13: Siemens Prisma
%     {default = 8}
%  'unwarpflag': [0|1|2]
%     0: unwarp 3D
%     1: unwarp through plan only
%     2: unwarp inplane only
%     {default = 0}
%  'isoctrflag': [0|1] whether to adjust for isocenter coordinates
%     {default = 1}
%  'jacobian_flag': [0|1] whether to apply jacobian
%     {default = 0}
%  'interpm': interpolation method
%      0 = nearest neighbor, 1 = linear, 2 = cubic
%      3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%     {default = 2}
%  'forceflag': overwrite existing output file
%     {default = 0}
%
% Output:
%   fname_out: output file name of grad warp corrected 4D volume
%
% Created:  02/23/10 by Don Hagler
% Prev Mod: 10/27/12 by Don Hagler
% Last Mod: 01/08/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms( varargin, { ...
  'fname_out',[],[],...
  'gwtype',8,[0:13],...
  'unwarpflag',0,[],...
  'isoctrflag',true,[false true],...
  'jacobian_flag',false,[false true],...
  'interpm',2,[0:5],...
  'forceflag',false,[false true],...
  'suffix','gruw',[],...
});

[tpath,tstem,text] = fileparts(fname_in);
if isempty(parms.fname_out)
  fname_out = [tpath '/' tstem '_' parms.sufix text];
else
  fname_out = parms.fname_out;
end;
parms.outdir = fileparts(fname_out);

if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_out,'file') || parms.forceflag
  mmil_mkdir(parms.outdir);
  [vol,M,mrparms,volsz] = fs_load_mgh(fname_in);
  nf = volsz(4);
  for f=1:nf
    vol_tmp = ctx_mgh2ctx(vol(:,:,:,f),M);
    vol_tmp = mmil_unwarp_grad(vol_tmp,...
      'gwtype',parms.gwtype,...
      'unwarpflag',parms.unwarpflag,...
      'isoctrflag',parms.isoctrflag,...
      'jacobian_flag',parms.jacobian_flag,...
      'interpm',parms.interpm);
    vol(:,:,:,f) = vol_tmp.imgs;
  end;
  fs_save_mgh(vol,parms.fname_out,M);
end;

