function M_volA_to_volB = mmil_reg(fname_volA,fname_volB,varargin)
%function M_volA_to_volB = mmil_reg(fname_volA,fname_volB,varargin)
%
% Purpose: registration between two MRI volumes using Dominic Holland's reg
%
% Usage: M_volA_to_volB = mmil_reg(fname_volA,fname_volB,'key1',val1,...)
%
% Required Input:
%   fname_volA: file name of volume A (mgh/mgz format)
%   fname_volB: file name of volume B (mgh/mgz format)
%
% Optional Input:
%   'outdir': output path
%     {default = pwd}
%   'fname_maskA': full path of brain mask for volume A
%     {default = []}
%   'fname_maskB': full path of brain mask for volume B
%     {default = []}
%   'rigid_flag': [0|1] whether to perform rigid reg
%     {default = 1}
%   'affine_flag': [0|1] whether to perform affine reg
%     NOTE: for best affine reg, use rigid_flag=1 and affine_flag=1
%     {default = 0}
%   'logfile': file name of output log file (relative to outdir)
%     {default = 'reg.log'}
%   'paramfile': full or relative path of parameter file
%     If relative, relative to $MMPS_PARMS/REG
%     {default = []}
%   'binfile': full or relative path of reg binary file
%     {default = 'reg'}
%   'cleanup_flag': [0|1] remove output files/folder
%     of residual step registration
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
% 
% Created:  04/24/10 by Don Hagler
% Last Mod: 10/30/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_volA',fname_volA,[],...
  'fname_volB',fname_volB,[],...
  'outdir',pwd,[],...
  'fname_maskA',[],[],...
  'fname_maskB',[],[],...
  'rigid_flag',true,[false true],...
  'affine_flag',false,[false true],...
  'logfile','reg.log',[],...
  'paramfile',[],[],...
  'binfile','reg',[],...
  'forceflag',false,[false true],...
...
  'sincflag',false,[false true],...
  'interpm',2,[1:5],...
  'fname_regmat','M_volA_to_volB.mat',[],...
  'fname_regtxt','M_volA_to_volB.txt',[],...
  'fname_regout',[],[],...
  'bindir',[],[],...
  'parmsdir',[],[],...
  'cleanup_flag',true,[false true],...
...
  'options',[],[],...
  'fstem_regtxt',[],[],... 
  'fstem_regout',[],[],... 
});

if parms.rigid_flag && parms.affine_flag
  tmp_options = '-rr -ar';
  parms.fstem_regtxt = 'affineRegMatrix';
  parms.fstem_regout = 'AffReg';
elseif parms.rigid_flag
  tmp_options= '-rr';
  parms.fstem_regtxt = 'rigidRegMatrix';% hard-coded in reg, depends on options
  parms.fstem_regout = 'RigidReg';% hard-coded in reg, depends on options
elseif parms.affine_flag
  tmp_options = '-ar';
  parms.fstem_regtxt = 'affineRegMatrix';
  parms.fstem_regout = 'AffReg';
end;

if ~isempty(parms.options)
  parms.options = [tmp_options ' ' parms.options];
else
  parms.options = tmp_options;
end; 

M_volA_to_volB = [];

nvox_LIA = [256 256 256,1];
M_LIA =  [-1     0     0   129;...
           0     0     1  -129;...
           0    -1     0   129;...
           0     0     0     1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

% check parameter file
if ~isempty(parms.paramfile)
  if isempty(parms.parmsdir)
    parms.parmsdir = [getenv('MMPS_PARMS') '/REG'];
  end;
  if isempty(regexp(parms.paramfile,'^/')) % relative path
    parms.paramfile = [parms.parmsdir '/' parms.paramfile];
  end;
  if ~exist(parms.paramfile,'file')
    error('parameter file %s not found',parms.paramfile);
  end;
end;

% check binary
if isempty(parms.bindir)
  parms.bindir = [getenv('MMPS_DIR') '/bin'];
end;
if mmil_isrelative(parms.binfile)
  parms.binfile = [parms.bindir '/' parms.binfile];
end;
if ~exist(parms.binfile,'file')
  error('binary file %s not found',parms.binfile);
end;

% set full log file name
if mmil_isrelative(parms.logfile)
  parms.logfile = [parms.outdir '/' parms.logfile];
end;

% set full regtxt file name
if mmil_isrelative(parms.fname_regtxt)
  parms.fname_regtxt = [parms.outdir '/' parms.fname_regtxt];
end;

% set full regmat file name
if mmil_isrelative(parms.fname_regmat)
  parms.fname_regmat = [parms.outdir '/' parms.fname_regmat];
end;

% set full regout file name
if ~isempty(parms.fname_regout) & mmil_isrelative(parms.fname_regout)
  parms.fname_regout = [parms.outdir '/' parms.fname_regout];
end;

% check input files exist
if ~exist(parms.fname_volA,'file')
  error('file %s not found',parms.fname_volA);
end;
[tmp,tstem_volA,text_volA] = fileparts(parms.fname_volA);
if ~exist(parms.fname_volB,'file')
  error('file %s not found',parms.fname_volB);
end;
[tmp,tstem_volB,text_volB] = fileparts(parms.fname_volB);
if ~isempty(parms.fname_maskA)
  if ~exist(parms.fname_maskA,'file')
    error('file %s not found',parms.fname_maskA);
  end;
  [tmp,tstem_maskA,text_maskA] = fileparts(parms.fname_maskA);
end;
if ~isempty(parms.fname_maskB)
  if ~exist(parms.fname_maskB,'file')
    error('file %s not found',parms.fname_maskB);
  end;
  [tmp,tstem_maskB,text_maskB] = fileparts(parms.fname_maskB);
end;

% create output dir if needed
mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resample volumes to standard space (LIA, 256^3) if necessary

[M_volA,nvox_volA,orient_volA] = fs_read_header(parms.fname_volA);
if ~all([M_volA(:);nvox_volA(:)]==[M_LIA(:);nvox_LIA(:)])
  fname_volA_LIA = [parms.outdir '/' tstem_volA '_LIA' text_volA];
  mmil_resample_to_LIA(parms.fname_volA,fname_volA_LIA,...
    'interpm',parms.interpm,'forceflag',parms.forceflag);
else
  fname_volA_LIA = parms.fname_volA;
end;

[M_volB,nvox_volB,orient_volB] = fs_read_header(parms.fname_volB);
if ~all([M_volB(:);nvox_volB(:)]==[M_LIA(:);nvox_LIA(:)])
  fname_volB_LIA = [parms.outdir '/' tstem_volB '_LIA' text_volB];
  mmil_resample_to_LIA(parms.fname_volB,fname_volB_LIA,...
    'interpm',parms.interpm,'forceflag',parms.forceflag);
else
  fname_volB_LIA = parms.fname_volB;
end;

if ~isempty(parms.fname_maskA)
  [M_maskA,nvox_maskA,orient_maskA] = fs_read_header(parms.fname_maskA);
  if ~all([M_maskA(:);nvox_maskA(:)]==[M_LIA(:);nvox_LIA(:)])
    fname_maskA_LIA = [parms.outdir '/' tstem_maskA '_LIA' text_maskA];
    mmil_resample_to_LIA(parms.fname_maskA,fname_maskA_LIA,...
      'interpm',parms.interpm,'forceflag',parms.forceflag);
  else
    fname_maskA_LIA = parms.fname_maskA;
  end;
else
  fname_maskA_LIA = [];
end;

if ~isempty(parms.fname_maskB)
  [M_maskB,nvox_maskB,orient_maskB] = fs_read_header(parms.fname_maskB);
  if ~all([M_maskB(:);nvox_maskB(:)]==[M_LIA(:);nvox_LIA(:)])
    fname_maskB_LIA = [parms.outdir '/' tstem_maskB '_LIA' text_maskB];
    mmil_resample_to_LIA(parms.fname_maskB,fname_maskB_LIA,...
      'interpm',parms.interpm,'forceflag',parms.forceflag);
  else
    fname_maskB_LIA = parms.fname_maskB;
  end;
else
  fname_maskB_LIA = [];
end;

[tmp,tstem_volB,text_volB] = fileparts(fname_volB_LIA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_regmat,'file') || parms.forceflag
  if ~exist(parms.fname_regtxt,'file') || parms.forceflag
    cmd = parms.binfile;
    cmd = sprintf('%s %s -t %s -s %s -od %s',...
      cmd,parms.options,fname_volA_LIA,fname_volB_LIA,parms.outdir);
    if ~isempty(fname_maskA_LIA)
      cmd = sprintf('%s -tm %s',cmd,fname_maskA_LIA);
    end;
    if ~isempty(fname_maskB_LIA)
      cmd = sprintf('%s -sm %s',cmd,fname_maskB_LIA);
    end;
    if parms.sincflag
      cmd = [cmd ' -sinc']; 
    else
      cmd = [cmd ' -cubic'];
    end;
    if ~isempty(parms.paramfile)
      cmd = [cmd ' -ip ' parms.paramfile];
    end;
    [status,result]=unix(cmd);
    mmil_logstr(parms,result);
    if status, error('registration failed:\n%s\n',result); end;

    % rename reg matrix
    if ~isempty(parms.fstem_regtxt) && ~isempty(parms.fname_regtxt)
      fname_tmp = [parms.outdir '/' parms.fstem_regtxt '.txt'];
      if ~exist(fname_tmp,'file')
        error('reg output %s not found',fname_tmp);
      end;
      cmd = ['mv ' fname_tmp ' ' parms.fname_regtxt];
      [status,result]=unix(cmd);
      if status, error('cmd %s failed:\n%s',cmd,result); end; 
    end;

    if ~isempty(parms.fstem_regout) && ~isempty(parms.fname_regout)
      % rename output file
      fname_tmp = sprintf('%s/%s_%s%s',...
        parms.outdir,tstem_volB,parms.fstem_regout,'.mgz');
      if ~exist(fname_tmp,'file')
        error('reg output %s not found',fname_tmp);
      end;
      cmd = ['mv ' fname_tmp ' ' parms.fname_regout];
      [status,result]=unix(cmd);
      if status, error('cmd %s failed:\n%s',cmd,result); end;
    end; 
  end;
  
  % there is residual registration from reg, so this additional step is added
  % use tmp_mmilreg folder
  tmp_outdir = sprintf('%s/tmp_mmilreg',parms.outdir);
  if ~exist(tmp_outdir,'dir')
    mmil_mkdir(tmp_outdir);
  end;
  
  % Apply M_volA_to_volB on fname_volB
  fname_volB2volA_resid= sprintf('%s/%s_resid%s',tmp_outdir,...
      tstem_volB,text_volB);
  [vol1,mr_parms1,volsz1]= ctx_load_mgh(fname_volA_LIA);
  [vol2,mr_parms2,volsz2]= ctx_load_mgh(fname_volB_LIA);
  M_vol1_to_vol2 = mmil_Mtxt2M(parms.fname_regtxt);% convert rigidRegMatrix to rbreg
  vol2r= vol_resample_pad(vol2,vol1,M_vol1_to_vol2,2,1);% cubic hard-coded      
  ctx_save_mgh(vol2r,fname_volB2volA_resid,mr_parms1);
  clear vol2 vol2r;
  % resample the mask
  if ~isempty(fname_maskB_LIA)
    [vol2_mask,mr_parms2_mask,volsz2_mask]=ctx_load_mgh(fname_maskB_LIA);
    vol2r_mask= vol_resample_pad(vol2_mask,vol1,M_vol1_to_vol2,0,1);% nearest neighbor hard-coded 
    fname_volB2volA_residmask= sprintf('%s/%s_residmask%s',tmp_outdir,...
      tstem_volB,text_volB);
    ctx_save_mgh(vol2r_mask,fname_volB2volA_residmask,mr_parms1);
    clear vol2_mask vol2r_mask;
  end;
  clear vol1;
  
  if isempty(parms.fname_regout)
    fname_res = sprintf('%s/%s_%s.mgz',parms.outdir,tstem_volB,parms.fstem_regout);
  else
    fname_res = parms.fname_regout;
  end;
   
  %  Registration between fname_volB2volA_residual to fname_regout
  cmd = parms.binfile;
  cmd = sprintf('%s -rr -t %s -s %s -od %s',cmd,fname_res,...
        fname_volB2volA_resid,tmp_outdir);
  if ~isempty(fname_maskA_LIA)
    cmd = sprintf('%s -tm %s',cmd,fname_maskA_LIA);
  end;  
  if ~isempty(fname_maskB_LIA)
    cmd = sprintf('%s -sm %s',cmd,fname_volB2volA_residmask);
  end;
  if parms.sincflag
    cmd = [cmd ' -sinc'];
  else
    cmd = [cmd ' -cubic'];
  end;
  if ~isempty(parms.paramfile)
    cmd = [cmd ' -ip ' parms.paramfile];
  end;
  [status,result]= unix(cmd);
  mmil_logstr(parms,result); 
  if status, error('residual registration failed:\n%s\n',result); end;
      
  % Combining the residual reg and initial registration
  fname_residreg = sprintf('%s/rigidRegMatrix.txt',tmp_outdir);
  M_volA_to_volB = mmil_Mtxt2M(parms.fname_regtxt) * mmil_Mtxt2M(fname_residreg);
  save(parms.fname_regmat,'M_volA_to_volB');
  
  % clean up the tmp_outdir directory
  if parms.cleanup_flag
    cmd = sprintf('rm -r %s',tmp_outdir);
    [status,result] = unix(cmd);
    if status, error('temporary directory deletion failed:\n%s\n',result); end;
  end;
else
  load(parms.fname_regmat);
end;
