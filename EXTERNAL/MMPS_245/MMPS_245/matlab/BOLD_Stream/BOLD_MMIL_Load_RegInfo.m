function [RegInfo,fname_reg,errcode] = BOLD_MMIL_Load_RegInfo(ContainerPath,varargin)
%function [RegInfo,fname_reg,errcode] = BOLD_MMIL_Load_RegInfo(ContainerPath,[options])
%
% Purpose: Load precalculated RegInfo from BOLD to T1 registration
%
% Usage:
%  [RegInfo,fname_reg,errcode] = BOLD_MMIL_Load_RegInfo(ContainerPath,'key1', value1,...);
%
% Required Input
%  ContainerPath: full path of directory containing processed BOLD data
%
% Optional Input:
%   'fname_reg': name of mat file with RegInfo struct
%     If empty, will generate based on fname_DTI
%     {default = []}
%   'fname_BOLD': full path name of 4d BOLD data registered to T1
%     If empty, will choose reference BOLD file based on the following
%     parameters using BOLD_MMIL_Get_ScanInfo's regT1_ref
%     {default = []}
%   'snums': vector of scan numbers that were processed (if empty, assume all)
%     {default = []}
%   'refsnum': reference scan number registered to T1 volume
%     if empty, will use reference number determined by BOLD_MMIL_Get_ScanInfo
%      {default = []}
%   'infix': BOLD file name infix
%      e.g. '', 'corr', 'corr_resBOLD', 'corr_regT1'
%      {default = []}
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
% Created:  08/16/10 by Don Hagler
% Last Mod: 03/07/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output
RegInfo = [];
fname_reg = [];
errcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_reg',[],[],...
  'fname_BOLD',[],[],...
  'snums',[],[],...
  'refsnum',[],[],...
  'infix',[],[],...
...
  'fnamestem','BOLD',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fname_reg)
  % choose BOLD file
  ref_snum = [];
  if isempty(parms.fname_BOLD)
    [ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
      'snums',parms.snums,'fnamestem',parms.fnamestem);
    if errcode || isempty(ScanInfo)
      errcode = 1;
      return;
    end;

    % choose reference BOLD scan
    if isempty(parms.refsnum) % use SessInfo ref
      parms.refsnum = SessInfo.regT1_ref;
    end;
    % check for bad scan num
    if ~ismember(parms.refsnum,SessInfo.snums_valid)
      fprintf('%s: ERROR: bad BOLD scan num (%d)\n',mfilename,parms.refsnum);
      errcode = 1;
      return;
    end;
    fstem_ref = [ContainerPath '/' ScanInfo(parms.refsnum).fstem];
    % check for first frame copy
    if ~isempty(parms.infix)
      parms.fname_BOLD = [fstem_ref '_f0_' parms.infix parms.ext];
    else
      parms.fname_BOLD = [fstem_ref '_f0' parms.ext];
    end;
    if ~exist(parms.fname_BOLD,'file')
      % look for full file
      if ~isempty(parms.infix)
        parms.fname_BOLD = [fstem_ref '_' parms.infix parms.ext];
      else
        parms.fname_BOLD = [fstem_ref parms.ext];
     end;
    end;
  end;
  if ~exist(parms.fname_BOLD,'file')
    fprintf('%s: ERROR: file %s not found\n',mfilename,parms.fname_BOLD);
    errcode = 1;
    return;
  end;
  [tpath_BOLD,tstem_BOLD,text_BOLD] = fileparts(parms.fname_BOLD);
  fname_reg = [ContainerPath '/' tstem_BOLD '_regT1.mat'];
else
  fname_reg = parms.fname_reg;
end;

if ~exist(fname_reg)
  fprintf('%s: ERROR: regT1 file %s not found\n',mfilename,fname_reg);
  errcode = 1;
  return;
end;

load(fname_reg);


return

