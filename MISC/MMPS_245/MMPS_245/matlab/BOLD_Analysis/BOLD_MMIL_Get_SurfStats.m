function [fname_surfstats,SubjIDs] = BOLD_MMIL_Get_SurfStats(RootDirs,varargin)
%function [fname_surfstats,SubjIDs] = BOLD_MMIL_Get_SurfStats(RootDirs,[options])
%
% Required Input:
%  RootDirs:
%    a struct which must contain the following fields:
%         proc_bold, fsurf
%    and may contain the following fields:
%         orig, raw, proc, proc_dti, proc_bold
%         fsurf, fsico, fsclean, long
%    these specify the locations of data
%
% Optional Input:
%  'SubjIDlist': Can be one subject ID or a cell array of IDs
%    If empty, use all subjects in StudyInfo
%     (or in RootDirs.fsurf if StudyInfo not supplied)
%    {defult = []}
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%      must contain these fields:
%        SubjID
%        StudyDate
%      may contain these fields
%        proc_bold
%        fsurf
%      if proc_bold and fsurf are unspecified, will look for Containers
%        with SubjID and StudyDate
%       (will choose first one if more than one)
%      if empty, use all subjects found in RootDirs.proc_bold and RootDirs.fsurf
%     {default = []}
%  'infix': if empty, will look for files like 'BOLD1.mgz'
%     otherwise, input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%     example infix = 'corr_resBOLD'
%     {default = []}
%  'hemi': cortical hemisphere
%     'lh' or 'rh'
%     {default = 'lh'}
%
% Output:
%   fname_surfstats: cell array of surface stats file names
%     Size will be number of subjects X 1
%     If more than one GLM_snums for a subject, will be nested cell array
%
% Created:  06/23/09 by Don Hagler
% Prev Mod: 03/23/12 by Don Hagler
% Prev Mod: 09/10/12 by Don Hagler
% Last Mod: 04/07/17 by Don Hagler
%

%% todo: use BOLD_MMIL_Set_GLM_Stem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_surfstats = [];
SubjIDs = [];

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'SubjIDlist',[],[],...
  'StudyInfo',[],[],...
  'infix',[],[],...
  'hemi','lh',{'lh' 'rh'},...
...
  'smoothsteps',0,[0,Inf],...
  'sphsmoothsteps',10,[0,Inf],...
  'snums',[],[],...
  'snums_tag','GLM_snums',[],...
  'fnamestem','BOLD',[],...
  'resT1flag',true,[false true],...
...
  'verbose',true,[false true],...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
});

% filter StudyInfo or generate StudyInfo struct if none supplied
args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(parms.StudyInfo,RootDirs,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.SubjIDlist)
  SubjIDs = {StudyInfo.SubjID};
elseif ~iscell(parms.SubjIDlist)
  SubjIDs = {parms.SubjIDlist};
else
  SubjIDs = parms.SubjIDlist;
end;
nsubs = length(SubjIDs);

if nsubs==0 | (nsubs==1 & isempty(SubjIDs{1}))
  error('no valid subjects given');
end;

fname_surfstats = cell(nsubs,1);
skipflags = zeros(nsubs,1);
for s=1:nsubs
  SubjID = SubjIDs{s};
  S_ind = find(strcmp(SubjIDs,SubjID));
  if isempty(S_ind)
    if parms.verbose
      fprintf('%s: WARNING: no studies with SubjID %s\n',...
        mfilename,SubjID);
    end;
    skipflags(s) = 1;
    continue;
  elseif length(S_ind)>1
    if parms.verbose
      fprintf('%s: WARNING: %d studies with SubjID %s\n',...
        mfilename,length(S_ind),SubjID);
    end;
    S_ind = S_ind(1);
  end;
  tmpinfo = StudyInfo(S_ind);
  StudyDate = tmpinfo.StudyDate;

  ContainerPath = sprintf('%s/%s',RootDirs.proc_bold,tmpinfo.proc_bold);
  FSContainerPath = sprintf('%s/%s',RootDirs.fsurf,tmpinfo.fsurf);

  % get number of scans, etc.
  if isfield(tmpinfo,'BOLDScaNums') & isempty(parms.snums)
    parms.snums = tmpinfo.BOLDScanNums;
  end;
  [ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
    'snums',parms.snums,'fnamestem',parms.fnamestem);
  if (SessInfo.nscans==0)
    if parms.verbose
      fprintf('%s: WARNING: no BOLD scans found in %s\n',...
        mfilename,ContainerPath);
    end;
    skipflags(s) = 1;
    continue;
  end;

  if ~isfield(tmpinfo,parms.snums_tag)
    if isempty(parms.snums)
      snums = [1:SessInfo.nscans];
    else
      snums = parms.snums;
    end;
  else
    snums = getfield(tmpinfo,parms.snums_tag);
    if isempty(snums)
      if parms.verbose
        fprintf('%s: WARNING: no %s for %s\n',...
          mfilename,parms.snums_tag,SubjID);
      end;
      skipflags(s) = 1;
      continue;
    end;
  end;

  % construct analysis file name(s)
  fnames = [];
  for i=1:length(snums)
    analdir = sprintf('%s/%s%d',...
      ContainerPath,parms.fnamestem,snums(i));
    if ~isempty(parms.infix)
      analdir = [analdir '_' parms.infix];
    end;
    analdir = [analdir '_analysis'];
    if ~exist(analdir,'dir')
      if parms.verbose
        fprintf('%s: WARNING: %s not found\n',mfilename,analdir);
      end;
      skipflags(s) = 1;
      break;
    end;
    tmp_fname = sprintf('%s/%s',analdir,ScanInfo(snums(i)).fstem);
    if ~isempty(parms.infix)
      tmp_fname = [tmp_fname '_' parms.infix];
    end;
    tmp_fname = [tmp_fname '_3dDeconv'];
    if parms.resT1flag
      tmp_fname = [tmp_fname '_resT1'];
    end;
    if parms.smoothsteps
      tmp_fname = sprintf('%s-sm%d',tmp_fname,parms.smoothsteps);
    end;
    tmp_fname = [tmp_fname '-sphere'];
    if parms.sphsmoothsteps
      tmp_fname = sprintf('%s-sm%d',tmp_fname,parms.sphsmoothsteps);
    end;
    tmp_fname = [tmp_fname '-' parms.hemi '.mgh'];
    if ~exist(tmp_fname,'file')
      if parms.verbose
        fprintf('%s: WARNING: %s not found\n',mfilename,tmp_fname);
      end;
      skipflags(s) = 1;
      break;
    end;
    fname_surfstats{s,i} = tmp_fname;
  end;
  if skipflags(s), continue; end;
end;

SubjIDs = SubjIDs(~skipflags);
fname_surfstats = fname_surfstats(~skipflags,:);

