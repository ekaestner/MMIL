function errcode = DTI_MMIL_Manual_Register_Override(ContainerPath,varargin)
%function errcode = DTI_MMIL_Manual_Register_Override(ContainerPath,[options])
%
% Usage:
%  DTI_MMIL_Manual_Register_Override(ContainerPath,'key1', value1,...);
%
% Required Input Parameters:
%   ContainerPath: full path of directory containing processed diffusion data (mgh format)
%
% Optional Input Parameters:
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr', 'corr_resDTI'
%     {default = 'corr'}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%
% Created:  05/18/09 Don Hagler
% Last Mod: 02/22/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options, set defaults

% parse input parameters
errcode = 0;
if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'infix',corr,[],...
  'revflag',0,[0,1,2],...
...
  'min_nb0',1,[],...
  'fnamestem','DTI',[],...
});

% load registration info
fprintf('%s: loading registration info...\n',mfilename);
tags = {'infix','min_nb0','revflag'};
args = mmil_parms2args(parms,tags);
[RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,args{:});
if errcode, return; end;
origRegInfo = RegInfo;

% get register.dat file name from regT1.mat file name
fname_regdat = regexprep(fname_reg,'_regT1.mat','_register.dat');
if ~exist(fname_regdat,'file')
  fprintf('%s: ERROR: %s not found\n',mfilename,fname_regdat);
  errcode = 1;
  return;
end;

fname_bak = sprintf('%s/%s_regT1_orig.mat',...
    ContainerPath,ftem);
if ~exist(fname_bak,'file')
  fprintf('%s: copying %s to %s...\n',...
    mfilename,fname_reg,fname_bak);
  cmd = sprintf('cp %s %s',fname_reg,fname_bak);
  [s,r] = unix(cmd);
  if s
    error('cmd %s failed:\n%s',cmd,r);
  end;
end;

% load register.dat file
RegInfo = origRegInfo;
[M_T1_to_T2,subj,inplane,slicethick] = fs_read_regdat(fname_regdat,...
  'tk2ras_flag',1,...
  'M_ref',RegInfo.M_T1,...
  'M_reg',RegInfo.M_T2,...
  'nvox_ref',RegInfo.volsz_T1,...
  'nvox_reg',RegInfo.volsz_T2);

% modify registration info
fprintf('%s: modifying registration info in %s...\n',...
  mfilename,fname_reg);
RegInfo.M_T1_to_T2 = M_T1_to_T2;
save(fname_reg,'M_T1_to_T2','RegInfo');


fprintf('%s: removing stale dependencies...\n',mfilename);
cleanup_list = {...
  'T1_resT2_DTI.mgh'...
  'T1_resDTI_regT1.mgh'...
  'fiber_paths_from_atlas'...
  'fiber_masks_from_atlas'...
  '*resT1*'...
  'DTcalc/*resT1*'...
  'RSIcalc/*resT1*'...
  'DTanalysis'...
  'RSIanalysis'...
};
for i=1:length(cleanup_list)
  if exist(target)
    cmd = sprintf('rm -r %s/%s',ContainerPath,cleanup_list{i});
    [s,r] = unix(cmd);
    if s
      fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
    end;
  end;
end;  

