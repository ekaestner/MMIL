function [M_T1_to_T2,RegInfo,errcode] = MMIL_Register_T2_to_T1(fname_T2,varargin)
% function [M_T1_to_T2,RegInfo,errcode] = MMIL_Register_T2_to_T1(fname_T2,[options])
%
% Usage:
%  MMIL_Register_T2_to_T1(fname_T2,'key1', value1,...);
%
% Required Parameters:
%   'fname_T2': full path name of T2-wegihted MRI volume
%
% Optional Parameters:
%   'fname_T1': full path name of T1-weighted image volume
%     if empty, will look for MPR_res.mgz or hiFA_res.mgz in T1ContainerPath
%     (preference depends on T1type)
%     {default = []}
%   'T1ContainerPath': ContainerPath for T1-weighted image volume
%     if empty, will use directory containing fname_T2
%     {default = []}
%   'T1type': select this type of T1 volume for target
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     Ignored if fname_T1 not empty
%     {default=2}
%   'FSContainerPath': full path of directory containing FreeSurfer recon
%     Will use FSContainerPath/mri/nu.mgz as fname_T1
%     Ignored if fname_T1 is not empty
%     Required if bbregflag=1 (uses mri/orig.mgz)
%     {default = []}
%   'fname_reg': output file name for mat file with RegInfo struct
%     If empty, will generate based on fname_T2
%     {default = []}
%   'bbregflag': [0|1] use FreeSurfer's bbregister (requires FSContainerPath)
%     {default = 0}
%   'T2_type': type of T2-weighted scan (e.g. 'T2','BOLD','DTI')
%      used for output T1_resT2.mgz file name
%    {default = 'T2'}
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%   'ext': file extension (i.e. '.mgh' or '.mgz')
%     {default = '.mgz'}
%   'forceflag': run calculations even if output files exist
%     {default = 0}
%
% Created:  04/21/10 by Don Hagler
% Last Mod: 02/15/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
errcode = 0;
M_T1_to_T2 = [];
RegInfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_T2',fname_T2,[],...
  'fname_T1',[],[],...
  'T1ContainerPath',[],[],...
  'T1type',2,[0,1,2,3],... 
  'FSContainerPath',[],[],...
  'fname_reg',[],[],...
  'bbregflag',false,[false true],...
  'T2_type','T2',[],...
  'interpm',2,[0:5],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
... % for jpdf
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...  
  'cleanup_flag',true,[false true],...
... % for creating brain mask if does not exist
  'mask_smooth',15,[0,Inf],...
});

if isempty(parms.fname_T2), error('fname_T2 is empty'); end;
if ~exist(parms.fname_T2,'file')
  fprintf('%s: ERROR: T2 file %s not found\n',mfilename,parms.fname_T2);
  errcode = 1;
  return;
end;
[ContainerPath,tstem_T2,text_T2] = fileparts(parms.fname_T2);

if isempty(parms.T1ContainerPath), parms.T1ContainerPath = ContainerPath; end;

if parms.bbregflag && isempty(parms.FSContainerPath)
  fprintf('%s: ERROR: bbregflag=1 and FSContainerPath is empty\n',mfilename);
  errcode = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose T1 file
tags = {'fname_T1','T1type','FSContainerPath','bbregflag','ext'};
args = mmil_parms2args(parms,tags);
[parms.fname_T1,errcode] = MMIL_Choose_T1(parms.T1ContainerPath,args{:});
if errcode, return; end;

% check whether T1 mask exist
[T1pathstr,T1stem,T1ext] = fileparts(parms.fname_T1);
fname_T1_mask = sprintf('%s/%s_mask%s',T1pathstr,T1stem,T1ext);
if exist(fname_T1_mask, 'file')
  parms.fname_T1_mask = fname_T1_mask; % this is a binary mask
else
  parms.fname_T1_mask = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get info about input files
[M_T1,volsz_T1] = ...
  mmil_load_mgh_info(parms.fname_T1,parms.forceflag,ContainerPath);
[M_T2,volsz_T2] = ...
  mmil_load_mgh_info(parms.fname_T2,parms.forceflag,ContainerPath);

% set output file names
if isempty(parms.fname_reg)
  parms.fname_reg = [ContainerPath '/' tstem_T2 '_regT1.mat'];
end;
if volsz_T2(4)==1
  fname_f0 = parms.fname_T2;
  fname_f0_res = [ContainerPath '/' tstem_T2 '_resT1' parms.ext];
else
  fname_f0 = [ContainerPath '/' tstem_T2 '_f0' parms.ext];
  fname_f0_res = [ContainerPath '/' tstem_T2 '_f0_resT1' parms.ext];
end;
[tmp,tstem_f0,text_f0] = fileparts(fname_f0);
fname_regdat = [ContainerPath '/' tstem_T2 '_register.dat'];
fname_regcsh = [ContainerPath '/' tstem_T2 '_register.csh'];
fname_tkmcsh = [ContainerPath '/' tstem_T2 '_tkmedit.csh'];
fname_T1_res = [ContainerPath '/T1_res' parms.T2_type parms.ext];

% check if output files exist
output_list = {parms.fname_reg,fname_f0,fname_f0_res,...
  fname_regdat,fname_regcsh,fname_T1_res};
if parms.forceflag
  runflag = 1;
else
  runflag = 0;
  for i=1:length(output_list)
    if ~exist(output_list{i},'file'), runflag = 1; break; end;
  end;
end;

if ~runflag
  load(parms.fname_reg);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save first frame of fname_T2 as fname_f0
vol_T2 = [];
if ~strcmp(parms.fname_T2,fname_f0) &...
   (~exist(fname_f0,'file') || parms.forceflag)
  fprintf('%s: saving first frame of T2 file %s...\n',...
    mfilename,parms.fname_T2);
  [vol,M_T2,tmp,volsz_T2]=fs_load_mgh(parms.fname_T2,[],1);
  vol_T2 = ctx_mgh2ctx(vol,M_T2);
  ctx_save_mgh(vol_T2,fname_f0);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% register T2 to T1
if ~exist(parms.fname_reg,'file') || parms.forceflag
  fprintf('%s: registering T2 to T1 (%s to %s)...\n',...
    mfilename,fname_f0,parms.fname_T1);
  if parms.bbregflag % FreeSurfer's bbregister
    mmil_bbregister(fname_f0,parms.FSContainerPath,...
      'fname_regdat',fname_regdat,'forceflag',1);
    % load register.dat file
    fprintf('%s: loading registration from %s...\n',...
      mfilename,fname_regdat);
    [M_T1_to_T2,subj,inplane,slicethick] = fs_read_regdat(fname_regdat,...
      'tk2ras_flag',1,...
      'M_ref',M_T1,...
      'M_reg',M_T2,...
      'nvox_ref',volsz_T1,...
      'nvox_reg',volsz_T2);
  else % Dominic Holland's jpdf reg
    tmpdir = [ContainerPath '/tmp_jpdf'];
    mmil_mkdir(tmpdir);
    try
      M_T1_to_T2 = mmil_jpdfreg_T1T2(parms.fname_T1,parms.fname_T2,...
        'outdir',tmpdir,'fname_T1_mask',parms.fname_T1_mask,...
        'atlasname',parms.atlasname,'atlasdir',parms.atlasdir,...
        'forceflag',parms.forceflag);
    catch
      fprintf('%s: ERROR: failed to register images:\n%s\n',mfilename,lasterr);
      errcode = 1;
      return;
    end;
    % remove tmp files
    if parms.cleanup_flag 
      cmd = sprintf('rm -r %s',tmpdir);
      [status,result] = unix(cmd);
      if status
        fprintf('%s: WARNING: failed to remove tmp dir %s\n',...
          mfilename,tmpdir);
      end;
    end;
  end;

  % save registration info
  RegInfo=[];
  RegInfo.M_T1_to_T2 = M_T1_to_T2;
  RegInfo.fname_T1 = parms.fname_T1;
  RegInfo.fname_T2 = parms.fname_T2;
  RegInfo.M_T2 = M_T2;
  RegInfo.M_T1 = M_T1;
  RegInfo.volsz_T2 = volsz_T2;
  RegInfo.volsz_T1 = volsz_T1;
  save(parms.fname_reg,'M_T1_to_T2','RegInfo');
else
  load(parms.fname_reg);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample T1 to T2 space
if ~exist(fname_T1_res,'file') || parms.forceflag
  fprintf('%s: resampling %s to T2 space...\n',...
    mfilename,parms.fname_T1);
  [vol_T1,M_T1] = fs_load_mgh(parms.fname_T1);
  [vol_T1_res,M_T1_res] = mmil_resample_vol(vol_T1,M_T1,...
      'M_ref',RegInfo.M_T2,'nvox_ref',RegInfo.volsz_T2(1:3),...
      'interpm',parms.interpm,'bclamp',1,...
      'M_reg',inv(RegInfo.M_T1_to_T2));
  fs_save_mgh(vol_T1_res,fname_T1_res,M_T1_res);
end;

%try
  mmil_write_regscripts(parms.fname_reg);
%catch
%  fprintf('%s: ERROR: unable to write reg scripts:\n%s\n',...
%    mfilename,lasterr);
%  errcode = 1;
%end;

