function M_volA_to_volB = epi_register_scans(fnameA,fnameB,varargin)
%function M_volA_to_volB = epi_register_scans(fnameA,fnameB,[options])
%
% Purpose: register one EPI scan to another
%
% Required Parameters:
%   fnameA: full path name of mgh/mgz file containing
%     reference 4D EPI volume
%   fnameB: full path name of mgh/mgz file containing 4D EPI volume
%     to be registered to reference
%
% Optional Parameters:
%   'fname_mask': input/output file name for brain mask of fnameA
%     If not empty and does not exist, will create one and save as fname_mask
%     If empty, will create one temporarily
%     {default=[]}
%   'fname_reg': output file name for registration matrix
%     {default=[]}
%   'rigid_flag': whether to perform rigid reg
%     {default = 1}
%   'affine_flag': whether to perform affine reg
%     for best affine reg, use rigid_flag=1 and affine_flag=1
%     {default=0}
%   'forceflag': overwrite existing output files
%     {default: 0}
%
% Created:  02/27/10 by Don Hagler
% Last Mod: 02/11/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_volA_to_volB = [];
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'fname_reg',[],[],...
  'rigid_flag',true,[false true],...
  'affine_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'tmpdir',[],[],...
  'cleanupflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input files
if ~exist(fnameA,'file')
  error('file %s not found',fnameA);
end;
if ~exist(fnameB,'file')
  error('file %s not found',fnameB);
end;

% check for output file
if ~isempty(parms.fname_reg) && exist(parms.fname_reg,'file') &&...
   ~parms.forceflag
  load(parms.fname_reg);
  return;
end;

if isempty(parms.tmpdir)
  if ~isempty(parms.fname_reg)
    [tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_reg);
  else
    [tmp_path,tmp_fstem,tmp_ext] = fileparts(fnameA);
  end;    
  if isempty(tmp_path), tmp_path = pwd; end;
  parms.tmpdir = [tmp_path '/tmp_EPI_reg'];
elseif exist(parms.tmpdir,'dir')
  parms.cleanupflag = false;
end;

if isempty(parms.fname_mask)
  parms.fname_mask = [parms.tmpdir '/mask.mgz'];
end;

mmil_mkdir(parms.tmpdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create brain mask from fnameA (if does not already exist)
epi_brainmask(fnameA,parms.fname_mask,'forceflag',parms.forceflag);

% run reg program
M_volA_to_volB = mmil_reg(fnameA,fnameB,...
  'fname_maskA',parms.fname_mask,...
  'outdir',parms.tmpdir,...
  'rigid_flag',parms.rigid_flag,...
  'affine_flag',parms.affine_flag);

% save output file
if ~isempty(parms.fname_reg)
  save(parms.fname_reg,'M_volA_to_volB');
end;

% remove tmp dir
if parms.cleanupflag
  cmd = sprintf('rm -r %s',parms.tmpdir);
  [status,result]=unix(cmd);
  if status
    fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n',mfilename,result);
  end;
end;

