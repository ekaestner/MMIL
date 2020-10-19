function errcode = BOLD_MMIL_Manual_Register_Override(ContainerPath,varargin)
%function errcode = BOLD_MMIL_Manual_Register_Override(ContainerPath,[options])
%
% Usage:
%  BOLD_MMIL_Manual_Register_Override(ContainerPath,'key1', value1,...);
%
% Required Input Parameters:
%   ContainerPath: full path of directory containing processed BOLD data
%
% Optional Input Parameters:
%   'snums' - vector of scan numbers that were processed (if empty, assume all)
%     {default = []}
%   'refsnum' - reference scan number registered to T1 volume
%     if empty, will use reference number determined by BOLD_MMIL_Get_ScanInfo
%      {default = []}
%   'infix' - BOLD file name infix
%      e.g. '', 'corr', 'corr_resBOLD', 'corr_regT1'
%      {default = 'corr'}
%
% Created:  08/16/10 by Don Hagler
% Last Mod: 03/23/12 by Don Hagler
%

%% todo: should call MMIL_Process_Exam or MMIL_Process_BOLD
%%  using ProjInfo?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options, set defaults

% parse input parameters
errcode = 0;
if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'refsnum',[],[],...
  'infix','corr',[],...
...
  'fnamestem','BOLD',[],...
});

% load registration info
fprintf('%s: loading registration info...\n',mfilename);
tags = {'snums','refsnum','infix','fnamestem'};
args = mmil_parms2args(parms,tags);
[RegInfo,fname_reg,errcode] = BOLD_MMIL_Load_RegInfo(ContainerPath,args{:});
if errcode, return; end;
origRegInfo = RegInfo;

% get register.dat file name from regT1.mat file name
fname_regdat = regexprep(fname_reg,'_regT1.mat','_register.dat');
if ~exist(fname_regdat,'file')
  fprintf('%s: ERROR: %s not found\n',mfilename,fname_regdat);
  errcode = 1;
  return;
end;

fname_bak = regexprep(fname_reg,'_regT1.mat','_regT1_orig.mat');
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
  'BOLD*tkmedit.csh'...
  'view_BOLD*.csh'...
  'T1_resBOLD.mgh'...
  'T1_resT2_BOLD.mgh'...
  'BOLD*regT1*.mg?'...
  'BOLD*resT1*.mg?'...
  'BOLD*sis/*resT1*.mg?'...
  'BOLD*sis/*-?h.mgh'...
};
if strcmp(parms.infix,'corr')
  cleanup_list{end+1} = 'BOLD*corr_resBOLD_regT1.mat';
  cleanup_list{end+1} = 'BOLD*corr_resBOLD_register.*';
end;

for i=1:length(cleanup_list)
  cmd = sprintf('rm -r %s/%s',ContainerPath,cleanup_list{i});
  [s,r] = unix(cmd);
  if s
    fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
  end;
end;  

