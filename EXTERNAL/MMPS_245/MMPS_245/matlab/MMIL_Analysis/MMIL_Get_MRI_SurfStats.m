function [fname_surfstats,SubjIDs,SessIDs,fspaths] =...
   MMIL_Get_MRI_SurfStats(RootDirs,varargin)
%function [fname_surfstats,SubjIDs,SessIDs,fspaths] =...
%  MMIL_Get_MRI_SurfStats(RootDirs,[options])
%
% Required Input:
%  RootDirs:
%    a struct which must contain the following fields:
%         proc, fsurf
%    and may contain the following fields:
%         orig, raw, proc
%         fsurf, fsclean, long
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
%        proc, fsurf
%      if proc and fsurf are unspecified, will look for Containers
%        with SubjID and StudyDate
%       (will choose first one if more than one)
%      if empty, use all subjects found in RootDirs.proc and RootDirs.fsurf
%     {default = []}
%  'datatype': Type of surface stats data
%     Allowed Data Types: {'thick','area','T1surf','T1cont','dv'}
%     {Default = 'thick'}
%  'hemi': cortical hemisphere
%     'lh' or 'rh'
%     {default = 'lh'}
%  'smoothing': surface smoothing steps (on sphere)
%     slope of FWHM vs. sqrt(N) is ~1.13 for fsaverage (v3)
%     (FWHM = full-width-half-max smoothing kernel
%         N = number of smoothing steps)
%      with [176,705,2819], approx FWHM (mm) = 15,30,60
%      with [100,300,3000], approx FWHM (mm) = 11.3, 19.6, 61.9
%     {default = 176}
%  'mask_midbrain_flag', [0|1] whether mid brain and other
%     cortical regions marked "unknown" were masked out
%     {default = 1}
%  'analysis_outdir': where MMIL_Analyze_MRI_Exam puts output files
%     relative to ContainerDirs
%     {default = 'analysis'}
%  'projdist_list': vector of projection distances for painting T1 onto surface
%     {default = [-0.2,0.2]}
%  'visitA': for dv, baseline visit number
%     {default = 1}
%  'visitB': for dv, followup visit number
%     {default = 2}
%
% Output:
%   fname_surfstats: cell array of surface stats file names
%     Size will be number of subjects X 1
%     If datatype = 'T1surf', size will be nsubs X 2
%       (for inside and outside white/gray boundary)
%   SubjIDs: cell array of subject IDs included in fname_surfstats
%   SessIDs: cell array of session IDs included in fname_surfstats
%   fspaths: cell array of full paths of freesurfer containers   
%
% Notes:
%   'thick' = cortical thickness
%   'area' = cortical area corrected for expansion relative to atlas
%   'T1surf' = T1-weighted image sampled on both sides of white/gray boundary
%   'T1cont' = T1-weighted contrast image (white vs. gray)
%
% Created:  05/12/09 by Don Hagler
% Last Mod: 01/22/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_surfstats = [];
SubjIDs = [];

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'SubjIDlist',[],[],...
  'StudyInfo',[],[],...
  'datatype','thick',{'thick','area','T1surf','T1cont','dv'},...
  'smoothing',176,[0,Inf],...
  'mask_midbrain_flag',true,[false true],...
  'analysis_outdir','analysis',[],...
  'hemi','lh',{'lh' 'rh'},...
  'projdist_list',[-0.2,0.2],[-5,5],...
  'visitA',1,[],...
  'visitB',2,[],...
...
  'ext','.mgz',{'.mgh','.mgz'},...
  'verbose',true,[false true],...
  'modality','MRI',[],...
});

switch parms.datatype
  case 'dv'
    parms.required_containers = {'long'};
  otherwise
    parms.required_containers = {'fsurf'};
end;

% filter StudyInfo or generate StudyInfo struct if none supplied
args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(parms.StudyInfo,RootDirs,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

% set generic file names
switch upper(parms.datatype)
  case 'THICK'
    fstem = 'thickness';
    if parms.mask_midbrain_flag, fstem = [fstem '-mbmask']; end;
    if parms.smoothing==0
      fnames = {sprintf('%s-sphere-%s%s',fstem,parms.hemi,parms.ext)};
    else
      fnames = {sprintf('%s-sphere-sm%d-%s%s',...
       fstem,parms.smoothing,parms.hemi,parms.ext)};
    end;
  case 'AREA'
    fstem = 'area';
    if parms.mask_midbrain_flag, fstem = [fstem '-mbmask']; end;
    if parms.smoothing==0
      fnames = {sprintf('%s-sphere-%s%s',fstem,parms.hemi,parms.ext)};
    else
      fnames = {sprintf('%s-sphere-sm%d-%s%s',...
       fstem,parms.smoothing,parms.hemi,parms.ext)};
    end;
  case 'T1SURF'
    fnames = [];
    for i=1:length(parms.projdist_list)
      fstem = sprintf('T1w_pdist%0.1f',parms.projdist_list(i));
      if parms.mask_midbrain_flag, fstem = [fstem '-mbmask']; end;
      if parms.smoothing==0
        fnames{i} = sprintf('%s-sphere-%ss%s',...
          fstem,parms.hemi,parms.ext);
      else
        fnames{i} = sprintf('%s-sphere-sm%d-%ss%s',...
          fstem,parms.smoothing,parms.hemi,parms.ext);
      end;
    end;
  case 'T1CONT'
    fnames = [];
    ii=0;
    for projdist=parms.projdist_list,
      if projdist <=0, continue; end
      ii=ii+1;
      fstem = sprintf('T1w_contrast_pdist%0.1f',projdist);
      if parms.mask_midbrain_flag, fstem = [fstem '-mbmask']; end;
      if parms.smoothing==0
        fnames{ii} = sprintf('%s-sphere-%ss%s',...
          fstem,parms.hemi,parms.ext);
      else
        fnames{ii} = sprintf('%s-sphere-sm%d-%ss%s',...
          fstem,parms.smoothing,parms.hemi,parms.ext);
      end;
    end;
  case 'DV'
    fstem = 'dv_Fine_paint';
    if parms.mask_midbrain_flag, fstem = [fstem '-mbmask']; end;
    if parms.smoothing==0
      fnames = {sprintf('%s-sphere-%ss%s',fstem,parms.hemi,parms.ext)};
    else
      fnames = {sprintf('%s-sphere-sm%d-%ss%s',...
       fstem,parms.smoothing,parms.hemi,parms.ext)};
    end;
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.SubjIDlist)
  SubjIDs = {StudyInfo.SubjID};
elseif ~iscell(parms.SubjIDlist)
  SubjIDs = {parms.SubjIDlist};
else
  SubjIDs = parms.SubjIDlist;
end;

if strcmp(upper(parms.datatype),'DV'),
  SubjIDs = unique(SubjIDs);
end

nsubs = length(SubjIDs);

if nsubs==0 | (nsubs==1 & isempty(SubjIDs{1}))
  error('no valid subjects given');
end;

fname_surfstats = cell(nsubs,length(fnames));
skipflags = zeros(nsubs,1);
if strcmp(upper(parms.datatype),'DV')
  SessIDs = cell(length(SubjIDs),2);
else
  SessIDs = cell(length(SubjIDs),1);
end
fspaths = cell(length(SubjIDs),1);
for s=1:nsubs
  SubjID = SubjIDs{s};
  S_ind = find(strcmp(SubjIDs,SubjID));
  if isempty(S_ind)
    if parms.verbose
      fprintf('%s: WARNING: no studies with SubjID %s\n',...
        mfilename,SubjID);
    end;
    continue;
  elseif ~strcmp(upper(parms.datatype),'DV') & length(S_ind)>1
    if parms.verbose
      fprintf('%s: WARNING: %d studies with SubjID %s\n',...
        mfilename,length(S_ind),SubjID);
    end;
    S_ind = S_ind(1);
  end;
  if strcmp(upper(parms.datatype),'DV')
    ind = strmatch(SubjID,{StudyInfo.SubjID},'exact');
    tt = find([StudyInfo(ind).VisitNumber]==parms.visitA);
    tmpinfoA = StudyInfo(ind(tt));
    tt = find([StudyInfo(ind).VisitNumber]==parms.visitB);
    tmpinfoB = StudyInfo(ind(tt));
    ContainerPath = sprintf('%s/%s',RootDirs.long,StudyInfo(s).long);
    if ~exist(ContainerPath)
      if parms.verbose
        fprintf('%s: WARNING: longitudinal container %s not found\n',...
          mfilename,ContainerPath);
      end;
      continue;
    end;
  else
    tmpinfo = StudyInfo(S_ind);
    StudyDate = tmpinfo.StudyDate;
    if isempty(tmpinfo.fsurf)
      if parms.verbose
        fprintf('%s: WARNING: no freesurfer recon for %s (%d)\n',...
          mfilename,SubjID,StudyDate);
      end;
      continue;
    end;
    ContainerPath = sprintf('%s/%s',RootDirs.fsurf,StudyInfo(s).fsurf);
    if ~exist(ContainerPath)
      if parms.verbose
        fprintf('%s: WARNING: freesurfer recon %s not found\n',...
          mfilename,ContainerPath);
      end;
      continue;
    end;
  end
  for i=1:length(fnames)
    if strcmp(upper(parms.datatype),'DV')
      analysis_outdir = sprintf('analysis/visit_%d_VS_%d',parms.visitA,parms.visitB);
      full_fname = sprintf('%s/%s/%s',ContainerPath,analysis_outdir,fnames{i});
    else
      full_fname = sprintf('%s/%s/%s',ContainerPath,parms.analysis_outdir,fnames{i});
    end;
    if exist(full_fname,'file')
      fname_surfstats{s,i} = full_fname;
    else
      if parms.verbose
        fprintf('%s: WARNING: %s not found\n',mfilename,full_fname);
      end;
      skipflags(s) = 1;
      break;
    end;
  end;
 
  if strcmp(upper(parms.datatype),'DV')
    fspaths{s} = sprintf('%s/%s',ContainerPath,analysis_outdir);
    nA = regexp(tmpinfoA.fsurf,'FSURF_(?<sess>.+)','names');
    nB = regexp(tmpinfoB.fsurf,'FSURF_(?<sess>.+)','names');
    if ~isempty(nA) & ~isempty(nB),
       SessIDs(s,1) = {nA.sess};
       SessIDs(s,2) = {nB.sess};
    end
  else
    fspaths{s} = ContainerPath;
    n = regexp(tmpinfo.fsurf,'FSURF_(?<sess>.+)','names');
    if ~isempty(n), SessIDs{s} = n.sess; end;
  end
end;
fname_surfstats = fname_surfstats(~skipflags,:);
SubjIDs = SubjIDs(~skipflags);
SessIDs = SessIDs(~skipflags,:);
fspaths = fspaths(~skipflags);
