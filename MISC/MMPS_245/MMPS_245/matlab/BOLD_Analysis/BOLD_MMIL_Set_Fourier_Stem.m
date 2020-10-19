function [stem,errcode] = BOLD_MMIL_Set_Fourier_Stem(ContainerPath,varargin)
%function [stem,errcode] = BOLD_MMIL_Set_Fourier_Stem(ContainerPath,[options])
%
% Purpose: set full path file stems of Fourier analysis output
%
% Required Input:
%   ContainerPath: full path of processed MRI container
%
% Optional Input ('key', value):
%  'snums': vector of BOLD scan numbers used in Fourier analysis
%    If empty, will assume all valid 'for' or 'rev' scans
%     (depending on SessInfo.revflag from BOLD_MMIL_Get_ScanInfo)
%    {default = []}
%  'snums_valid': vector of BOLD scan numbers used in processing
%    If empty, will assume all scans
%    {default = []}
%  'adir': analysis dir relative to ContainerPath
%    If empty, will construct from snums
%    {default = []}
%  'stem': file stem of Fourier analysis output
%    If empty, will construct from snums
%    {default = []}
%  'infix': string inside BOLD file names (e.g. 'corr_resBOLD')
%    {default = 'corr_resBOLD'}
%  'fstats_infix': string inside Fourier analysis file names
%    {default = 'fstats_pval'}
%
% Output:
%   stem: full path file stem for Fourier analysis files
%   errcode: returns 1 if unable to read ContainerInfo or no BOLD scans found
%
% Created:  03/12/11 by Don Hagler
% Last Mod: 03/23/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'snums_valid',[],[],...
  'adir',[],[],...
  'stem',[],[],...
  'infix','corr_resBOLD',[],...
  'fstats_infix','fstats_pval',[],...
...
  'fnamestem','BOLD',[],...
});
adir = [];
stem = [];
errcode = 0;

if isempty(parms.adir) || isempty(parms.stem)
  if isempty(parms.snums) || length(parms.snums)==1
    % get info about BOLD MRI scans
    [ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
      'snums',parms.snums_valid,'fnamestem',parms.fnamestem);
    if isempty(ScanInfo), errcode = 1; end;
    if errcode, return; end;
    if isempty(parms.snums)
      if SessInfo.revflag
        parms.snums = SessInfo.snums_rev;
      else
        parms.snums = SessInfo.snums_for;
      end;
    end;
  end;

  if length(parms.snums) > 1
    % construct scaninfix for combined output
    fstem = sprintf('%s_scans%s_%s',...
      parms.fnamestem,sprintf('_%d',parms.snums),parms.infix);
    dstem = [fstem '_analysis'];
    fstem = [fstem '_avg'];
  else
    fstem = sprintf('%s%d',parms.fnamestem,parms.snums);
    dstem = [fstem '_' parms.infix '_analysis'];
    if strcmp(ScanInfo(parms.snums).ScanType,'BOLD_ape')
      fstem = [fstem '_ape'];
    elseif SessInfo.revflag
      fstem = [fstem '_rev'];
    end;
    fstem = [fstem '_' parms.infix];
  end;

  if isempty(parms.adir), parms.adir = dstem; end;
  if isempty(parms.stem), parms.stem = fstem; end;
end;
parms.stem = [parms.stem '_' parms.fstats_infix];

stem = sprintf('%s/%s/%s',ContainerPath,parms.adir,parms.stem);
