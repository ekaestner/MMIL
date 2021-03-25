function mmil_write_regscripts(fname_reg,varargin)
%function mmil_write_regscripts(fname_reg,[options])
%
% Usage:
%  mmmil_write_regscripts(fname_reg,'key1', value1,...);
%
% Required Parameters:
%   fname_reg: input file name for mat file with RegInfo struct
%
% Optional Parameters:
%   'fname_T1': full path name of T1-weighted image volume
%     If supplied, override RegInfo
%   'fname_T2': full path name of T2-wegihted MRI volume
%     If supplied, override RegInfo
%   'FSContainerPath': full path of directory containing FreeSurfer recon
%     If supplied, display surfaces
%   'outdir': output directory
%     {default = path of fname_T2}
%   'ext': output file extension (i.e. '.mgh' or '.mgz')
%     {default = '.mgz'}
%   'forceflag': overwrite existing output files
%     {default = 0}
%
% Created:  05/04/10 by Don Hagler
% Prev Mod: 02/16/13 by Don Hagler
% Last Mod: 11/06/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_reg',fname_reg,[],...
  'fname_T2',[],[],...
  'fname_T1',[],[],...
  'FSContainerPath',[],[],...
  'outdir',[],[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
...
  'interpm',2,[1:5],...
});

require_exist_tags = {'fname_T1','fname_T2','fname_reg','FSContainerPath'};
replace_reg_tags = {'fname_T1','fname_T2'};

load(parms.fname_reg);

if ~exist(RegInfo.fname_T1,'file')
  [procpath,~] = fileparts(parms.fname_reg);
  [procpath,~] = fileparts(procpath);
  [projpath,~] = fileparts(procpath);
  [ContainerPath_T1,fstem_T1,ext_T1] = fileparts(RegInfo.fname_T1);
  [procpath_T1,ContainerName_T1,ContainerName_ext_T1] = fileparts(ContainerPath_T1);
  [~,procname_T1] = fileparts(procpath_T1);
  RegInfo.fname_T1 = strcat([projpath '/' procname_T1 '/' ContainerName_T1 ContainerName_ext_T1 '/' fstem_T1 ext_T1]);
  
  [ContainerPath_T2,fstem_T2,ext_T2] = fileparts(RegInfo.fname_T2);
  [procpath_T2,ContainerName_T2,ContainerName_ext_T2] = fileparts(ContainerPath_T2);
  [~,procname_T2] = fileparts(procpath_T2);
  RegInfo.fname_T2 = strcat([projpath '/' procname_T2 '/' ContainerName_T2 ContainerName_ext_T2 '/' fstem_T2 ext_T2]);
end

for i=1:length(replace_reg_tags)
  tag = replace_reg_tags{i};
  parm_val = getfield(parms,tag);
  if isempty(parm_val)
    reg_val = getfield(RegInfo,tag);
    parms = setfield(parms,tag,reg_val);
  end;
end;

for i=1:length(require_exist_tags)
  tag = require_exist_tags{i};
  val = getfield(parms,tag);
  if ~isempty(val) & ~exist(val,'file')
    error('%s not found',val);
  end;    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set output file names
[M_T2,volsz_T2] = ...
  mmil_load_mgh_info(parms.fname_T2,parms.forceflag,parms.outdir);
[tpath,tstem,text] = fileparts(parms.fname_T2);
if isempty(parms.outdir), parms.outdir = tpath; end;
fname_regdat = [parms.outdir '/' tstem '_register.dat'];
fname_regcsh = [parms.outdir '/' tstem '_register.csh'];
fname_tkmcsh = [parms.outdir '/' tstem '_tkmedit.csh'];
if volsz_T2(4)==1
  fname_f0 = parms.fname_T2;
  fname_f0_res = [parms.outdir '/' tstem '_resT1' parms.ext];
else
  fname_f0 = [parms.outdir '/' tstem '_f0' parms.ext];
  fname_f0_res = [parms.outdir '/' tstem '_f0_resT1' parms.ext];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% write register.dat
fprintf('%s: writing %s...\n',mfilename,fname_regdat);
fs_write_regdat(fname_regdat,...
  'M',RegInfo.M_T1_to_T2,...
  'ras2tk_flag',1,...
  'M_ref',RegInfo.M_T1,...
  'M_reg',RegInfo.M_T2,...
  'nvox_ref',RegInfo.volsz_T1,...
  'nvox_reg',RegInfo.volsz_T2,...
  'forceflag',parms.forceflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create script to run tkregister2
if ~exist(fname_regcsh) || parms.forceflag
  fid = fopen(fname_regcsh,'wt');
  if fid==-1
    error('unable to write %s',fname_regcsh);
  end;
  fprintf(fid,'#!/bin/csh -f\n');
  fprintf(fid,'tkregister2 --targ %s \\\n  --mov %s \\\n',...
    parms.fname_T1,fname_f0);
  if ~isempty(parms.FSContainerPath)
    [subjdir,subj,ext] = fileparts(parms.FSContainerPath);
    subj = [subj ext];
    fprintf(fid,'  --surf --s %s --sd %s\\\n',...
      subj,subjdir);
  end;
  fprintf(fid,'  --reg %s\n',fname_regdat);
  fclose(fid);
  % change permissions so group can write and execute
  cmd = sprintf('chmod ug+rwx %s',fname_regcsh);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: failed to set permissions for %s:\n%s\n',...
      mfilename,fname_regcsh,result);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create script to run tkmedit
if ~exist(fname_tkmcsh) || parms.forceflag
  fid = fopen(fname_tkmcsh,'wt');
  if fid==-1
    error('unable to write %s',fname_tkmcsh);
  end;
  fprintf(fid,'#!/bin/csh -f\n');
  fprintf(fid,'tkmedit -f %s \\\n  -aux %s\n',...
    parms.fname_T1,fname_f0_res);
  fclose(fid);
  % change permissions so group can write and execute
  cmd = sprintf('chmod ug+rwx %s',fname_tkmcsh);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: failed to set permissions for %s:\n%s\n',...
      mfilename,fname_tkmcsh,result);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save first frame of fname_T2 as fname_f0
if ~strcmp(parms.fname_T2,fname_f0) &...
   (~exist(fname_f0,'file') || parms.forceflag)
  fprintf('%s: saving first frame of T2 file %s...\n',...
    mfilename,parms.fname_T2);
  [vol,M_T2,tmp,volsz_T2]=fs_load_mgh(parms.fname_T2,[],1);
  vol_T2 = ctx_mgh2ctx(vol,M_T2);
  ctx_save_mgh(vol_T2,fname_f0);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply registration to fname_f0
if ~exist(fname_f0_res,'file') || parms.forceflag
  fprintf('%s: resampling %s to T1 space...\n',...
    mfilename,fname_f0);
  if isempty(RegInfo), load(parms.fname_reg); end;
  vol_T1 = ctx_mgh2ctx(zeros(RegInfo.volsz_T1),RegInfo.M_T1);
  if ~exist('vol_T2','var')
    vol_T2 = ctx_load_mgh(fname_f0);
  end;
  vol_T2_res = vol_resample_pad(vol_T2,vol_T1,...
    RegInfo.M_T1_to_T2,parms.interpm);
  ctx_save_mgh(vol_T2_res,fname_f0_res);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
