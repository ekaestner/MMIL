function [stem,errcode] = BOLD_MMIL_Set_GLM_Stem(ContainerPath,varargin)
%function [stem,errcode] = BOLD_MMIL_Set_GLM_Stem(ContainerPath,[options])
%
% Purpose: set full path file stems of BOLD GLM analysis output
%
% Required Input:
%   ContainerPath: full path of processed MRI container
%
% Optional Input ('key', value):
%  'snums': vector of BOLD scan numbers used in GLM analysis
%    If empty, will assume all valid 'for' or 'rev' scans
%     (depending on SessInfo.revflag from BOLD_MMIL_Get_ScanInfo)
%    {default = []}
%  'snums_valid': vector of BOLD scan numbers used in processing
%    If empty, will assume all scans
%    {default = []}
%  'adir': analysis dir relative to ContainerPath
%    If empty, will construct from snums
%    {default = []}
%  'stem': file stem of GLM analysis output
%    If empty, will construct from snums
%    {default = []}
%  'infix': string inside BOLD file names (e.g. 'corr_resBOLD')
%    {default = 'corr_resBOLD'}
%  'results_infix': string inside GLM analysis file names
%    {default = '3dDeconv'}
%  'pthresh': probability threshold applied to f-stats for each condition
%    {default: 0}
%  'resamp_flag': [0|1] whether resampled to 1x1x1 before painting
%    {default: 0}
%  'smoothsteps': smoothing steps on surface after painting
%    {default: 0}
%  'sphere_flag': [0|1] whether sampled to icosohedral sphere
%    {default: 0}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default: 0}
%
% Output:
%   stem: full path file stem for GLM analysis files (except '-lh.mgh')
%   errcode: returns 1 if unable to read ContainerInfo or no BOLD scans found
%
% Created:  03/22/11 by Don Hagler
% Prev Mod: 03/23/11 by Don Hagler
% Last Mod: 05/15/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'snums_valid',[],[],...
  'adir',[],[],...
  'stem',[],[],...
  'infix','corr_resBOLD',[],...
  'results_infix','3dDeconv',[],...
  'resamp_flag',false,[false true],...
  'pthresh',0,[0,1],...
  'smooth',0,[0,50],...
  'smoothsteps',0,[0,1000],...
  'sphere_flag',false,[false true],...
  'sphsmoothsteps',0,[0,1000],...
...
  'fnamestem','BOLD',[],...
  'dirstem','BOLD',[],...
  'analysis_outfix','analysis',[],...
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
      parms.dirstem,sprintf('_%d',parms.snums),parms.infix);
    dstem = sprintf('%s_%s',fstem,parms.analysis_outfix);
  else
    fstem = sprintf('%s_scan_%d',parms.dirstem,parms.snums);
    dstem = sprintf('%s_%s_%s',fstem,parms.infix,parms.analysis_outfix);
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

parms.stem = [parms.stem '_' parms.results_infix];
if parms.resamp_flag
  parms.stem = [parms.stem '_resT1'];
end;
if parms.smoothsteps
  parms.stem = sprintf('%s-sm%d',parms.stem,parms.smoothsteps);
end;
if parms.pthresh
  parms.stem = sprintf('%s_pthresh%0.1e',parms.stem,parms.pthresh);
end;
if parms.sphere_flag
  parms.stem = [parms.stem '-sphere'];
  if parms.sphsmoothsteps
    parms.stem = sprintf('%s-sm%d',parms.stem,parms.sphsmoothsteps);
  end;
end;

stem = sprintf('%s/%s/%s',ContainerPath,parms.adir,parms.stem);
