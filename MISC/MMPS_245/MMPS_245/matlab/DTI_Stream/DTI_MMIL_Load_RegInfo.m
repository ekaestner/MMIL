function [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,varargin)
%function [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,[options])
%
% Usage:
%  [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,'key1', value1,...);
%
% Required Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%
% Optional Parameters specifying DTI image:
%   'fname_reg': name of mat file with RegInfo struct
%     If empty, will generate based on fname_DTI
%     {default = []}
%   'fname_DTI': full path name of 4d diffusion data or b=0 volume registered to T1
%     If empty, will choose reference DTI file based on the following
%     parameters using DTI_MMIL_Get_ScanInfo's reg_ref_for or reg_ref_rev
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr'
%     {default = []}
%   'revflag': [0|1|2] specify whether to process non-rev or rev data
%     0: register forward phase-encode polarity data
%     1: register reverse phase-encode polarity data
%     2: register either forward or reverse data,
%        depending on which is more common
%     {default = 2}
%
% Output:
%   RegInfo: struct containing these fields
%     M_T1_to_T2
%     M_T1
%     M_T2
%     fname_T1
%     fname_T2
%   fname_reg: name of registration mat file
%   errcode: 0 if success, 1 if error
%
% Created:  01/01/07 by Don Hagler
% Prev Mod: 02/18/13 Don Hagler
% Last Mod: 11/06/17 Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
RegInfo = [];
fname_reg = [];
errcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options, set defaults

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_reg',[],[],...
  'fname_DTI',[],[],...
  'infix',[],[],...
  'revflag',2,[0,1,2],...
...
  'min_nb0',[],[],...
  'fnamestem','DTI',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fname_reg)
  % choose DTI file
  ref_snum = [];
  if isempty(parms.fname_DTI)
    % load ContainerInfo, get DTI scan info, determine valid scans, reference scans
    tags = {'revflag','min_nb0'};
    args = mmil_parms2args(parms,tags);
    [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
    if errcode || isempty(ScanInfo)
      errcode = 1;
      return;
    end;
    if ~SessInfo.revflag
      ref_snum = SessInfo.reg_ref_for;
      fstem_ref = sprintf('%s/%s%d',...
        ContainerPath,parms.fnamestem,ref_snum);
    else
      ref_snum = SessInfo.reg_ref_rev;
      fstem_ref = sprintf('%s/%s%d_rev',...
        ContainerPath,parms.fnamestem,ref_snum);
    end;
    % check for first frame copy
    if ~isempty(parms.infix)
      parms.fname_DTI = [fstem_ref '_f0_' parms.infix parms.ext];
    else
      parms.fname_DTI = [fstem_ref '_f0' parms.ext];
    end;
    if ~exist(parms.fname_DTI,'file')
      % look for full file
      if ~isempty(parms.infix)
        parms.fname_DTI = [fstem_ref '_' parms.infix parms.ext];
      else
        parms.fname_DTI = [fstem_ref parms.ext];
     end;
    end;
  end;
  if ~exist(parms.fname_DTI,'file')
    fprintf('%s: ERROR: DTI file %s not found\n',mfilename,parms.fname_DTI);
    errcode = 1;
    return;
  end;
  % set name of regT1 mat file
  [tpath_DTI,tstem_DTI,text_DTI] = fileparts(parms.fname_DTI);
  fname_reg = [ContainerPath '/' tstem_DTI '_regT1.mat'];
else
  fname_reg = parms.fname_reg;
end;

if ~exist(fname_reg)
  fprintf('%s: ERROR: regT1 file %s not found\n',mfilename,fname_reg);
  errcode = 1;
  return;
end;

load(fname_reg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Force fname_T1 and fname_T2 to be on the same path. which is almost always true;
%%Feng Xue 2017/11/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(RegInfo.fname_T1,'file')
  [procpath,~] = fileparts(ContainerPath);
  [projpath,procname] = fileparts(procpath);
  [ContainerPath_T1,fstem_T1,ext_T1] = fileparts(RegInfo.fname_T1);
  [procpath_T1,ContainerName_T1,ContainerName_ext_T1] = fileparts(ContainerPath_T1);
  [~,procname_T1] = fileparts(procpath_T1);
  RegInfo.fname_T1 = strcat([projpath '/' procname_T1 '/' ContainerName_T1 ContainerName_ext_T1 '/' fstem_T1 ext_T1]);
  
  [ContainerPath_T2,fstem_T2,ext_T2] = fileparts(RegInfo.fname_T2);
  [procpath_T2,ContainerName_T2,ContainerName_ext_T2] = fileparts(ContainerPath_T2);
  [~,procname_T2] = fileparts(procpath_T2);
  RegInfo.fname_T2 = strcat([projpath '/' procname_T2 '/' ContainerName_T2 ContainerName_ext_T2 '/' fstem_T2 ext_T2]);
end

return

