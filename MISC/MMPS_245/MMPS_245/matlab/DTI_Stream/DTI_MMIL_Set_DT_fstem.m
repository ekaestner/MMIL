function [fstem,snums] = DTI_MMIL_Set_DT_fstem(ContainerPath,varargin)
%function [fstem,snums] = DTI_MMIL_Set_DT_fstem(ContainerPath,[options])
%
% Usage:
%  [fstem,snums] = ...
%     DTI_MMIL_Set_DT_fstem(ContainerPath,'key1', value1,...);
%
% Required Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%
% Optional Parameters:
%   'outdir': output directory
%     absolute or relative to ContainerPath
%     {default = 'DTcalc'}
%   'outfix': string attached to output file names
%     {default = []}
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr_regT1'
%     {default = []}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'min_ndirs': require at least this many diffusion directions to be valid
%     {default = 6}
%   'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%   'max_bval': maximum b-value used in tensor fit
%     {default = Inf}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%
% Output:
%   fstem: full path file stem for DT calculation output files
%   snums: vector of scan numbers used for DT calculations
%
% Created:  09/04/07 by Don Hagler
% Last Mod: 06/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options, set defaults

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir','DTcalc',[],...
  'outfix',[],[],...
  'snums',[],[],...
  'infix',[],[],...
  'revflag',0,[0:2],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'max_bval',Inf,[100,Inf],...
  'flex_flag',false,[false true],...
  'min_nb0',1,[],...
  'nob0_flag',false,[false true],...
...
  'fnamestem','DTI',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'info_tags',{'snums','revflag','min_nb0','min_ndirs',...
               'min_bval','flex_flag'},[],...
});

fstem = []; snums = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get valid DTI scan numbers
args = mmil_parms2args(parms,parms.info_tags);
[ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
if errcode || isempty(ScanInfo)
  fprintf('%s: ERROR: no DTI scans\n',mfilename);
  return;
end;

% exclude snums for files that do not exist (e.g. processing errors)
for s=SessInfo.snums_DT
  revflag_list = set_revflag_list(ScanInfo(s),parms.revflag);
  for revflag=revflag_list
    fstem = sprintf('%s/%s%d',ContainerPath,parms.fnamestem,s);
    if revflag, fstem = [fstem '_rev']; end;
    if ~isempty(parms.infix), fstem = [fstem '_' parms.infix]; end;
    fname = [fstem parms.ext];
    if exist(fname,'file'), snums = [snums,s]; end;
  end;
end;
if isempty(snums)
  fprintf('%s: ERROR: no valid DTI scans\n',mfilename);
  return;
end;

% set full output directory
if mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;

% set file stem
fstem = [parms.outdir '/' parms.fnamestem];
scaninfix = sprintf('%d',snums(1));
if length(snums)>1
  for i=2:length(snums)
    scaninfix = sprintf('%s_%d',scaninfix,snums(i));
  end;
  scaninfix = sprintf('_scans_%s',scaninfix);
end;
switch parms.revflag
  case 1
    scaninfix = sprintf('%s_rev',scaninfix);  
  case 2
    scaninfix = sprintf('%s_crev',scaninfix);
end;
fstem = sprintf('%s%s',fstem,scaninfix);
if ~isempty(parms.infix), fstem = sprintf('%s_%s',fstem,parms.infix); end;
fstem = sprintf('%s_DT',fstem);
if parms.nob0_flag
  fstem = [fstem 'nb'];
end;
if parms.min_bval > 1
  fstem = sprintf('%s_minb%d',fstem,parms.min_bval);
end
if parms.max_bval < Inf
  fstem = sprintf('%s_maxb%d',fstem,parms.max_bval);
end
if parms.flex_flag
  fstem = [fstem '_flex'];
end;
if ~isempty(parms.outfix)
  fstem = [fstem '_' parms.outfix];
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revflag_list = set_revflag_list(ScanInfo,revflag)
  revflag_list = [];
  if ismember(ScanInfo.pepolar,[0,1])
    revflag_list = ScanInfo.pepolar;
  elseif revflag==2
    if ismember(ScanInfo.DTI_Sequence_Type,[4,5,6])
      if ScanInfo.pepolar==2 % main scan is "forward"
        revflag_list = 0;
      else                   % main scan is "reverse"
        revflag_list = 1;
      end;
    else
      revflag_list = [0,1];
    end;
  else
    revflag_list = revflag;
  end;
return;

