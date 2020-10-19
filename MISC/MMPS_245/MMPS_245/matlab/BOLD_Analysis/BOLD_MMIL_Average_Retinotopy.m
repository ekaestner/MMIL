function BOLD_MMIL_Average_Retinotopy(ProjID,varargin)
%function BOLD_MMIL_Average_Retinotopy(ProjID,[options])
%
% Usage:
%  BOLD_MMIL_Average_Retinotopy(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string, used to to load ProjInfo and StudyInfo
%     from user's home directory
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    If empty, will use ProjID to get StudyInfo
%    May control which subjects are included in average with
%     fields 'GroupAvg_Polar' and 'GroupAvg_Ecc'
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    If empty, will use ProjID to get RootDirs
%    {default = []}
%  'datatypes': cell array of types of retinotopy data
%    {default = {'polar','eccen'}}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'multi_flag': [0|1|2] indicates source of input analyses
%    0: use results from individual scans
%    1: use results from multiple scans within session
%    2: use results from multiple sessions
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0 }
%
% Optional Parameters that specify how analysis was done
%  'infix': string attached to processed BOLD file names
%    {default = 'corr_resBOLD'}
%  'fstats_type': [0|1|2] how output Fourier components should be scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default = 0}
%  'resamp_flag': [0|1] whether fourier results were resampled to 1x1x1
%     before painting
%    { default = 0 }
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 10}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default = 3}
%
% Created:  04/10/10 by Don Hagler
% Last Mod: 08/01/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
if isempty(ProjID), error('empty ProjID'); end;

parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'outdir','Average_Retinotopy',[],...
  'hemilist',{'lh','rh'},[],...
  'fstats_type',0,[0:2],...
  'resamp_flag',false,[false true],...
  'sphere_flag',true,[false true],...
  'smoothsteps',10,[0,1000],...
  'sphsmoothsteps',3,[0,1000],...
  'infix','corr_resBOLD',[],...
  'datatypes',{'polar','eccen'},{'polar','eccen'},...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'suffixlist',{'_r','_i'},{'_r','_i'},...
  'qcflag',true,[false true],...
  'multi_flag',1,[0:2],...
  'forceflag',false,[false true],...
...
  'fnamestem','BOLD',[],...
  'required_containers',{'proc_bold','fsurf'},[],...
  'cxfstatsflag',true,[false true],...
...
  'snums_tags',{'FP_snums','FE_snums'},[],...
  'avg_tags',{'GroupAvg_Polar','GroupAvg_Ecc'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if  mmil_isrelative(parms.outdir)
  parms.outdir = ...
    sprintf('%s/MetaData/%s/%s',...
  RootDirs.home,ProjID,parms.outdir);
end;
parms.OutContainerPath = parms.outdir;
parms = rmfield(parms,'outdir');

% which to pass to MMIL_Combine_Retinotopy_Exam
tags_excl = {'RootDirs','StudyInfo','required_containers','qcflag',....
  'snums_tags','avg_tags','fnamestem','hemilist','suffixlist'};
tags = setdiff(fieldnames(parms),tags_excl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d=1:length(parms.datatypes)
  snums_tag = parms.snums_tags{d};
  avg_tag = parms.avg_tags{d};

  if ~isfield(StudyInfo,snums_tag)
    fprintf('%s: WARNING: StudyInfo is missing field %s\n',...
      mfilename,snums_tag);
    continue;
  end;
  ind = find(~cellfun(@isempty,{StudyInfo.(snums_tag)}));
  if isempty(ind)
    fprintf('%s: WARNING: no sessions with %s retintotopy data\n',...
      mfilename,parms.datatypes{d});
    continue;
  end;

  if isfield(StudyInfo,avg_tag)
    tmp = {StudyInfo.(avg_tag)};
    ind_avg = find(~cellfun(@isempty,tmp));
    tmp2 = zeros(length(StudyInfo));
    tmp2(ind_avg) = tmp{ind_avg};
    ind_avg = find(tmp2);
    if isempty(ind_avg)
      fprintf('%s: WARNING: no sessions with %s = 1\n',...
        mfilename,avg_tag);
      continue;
    end;
    ind = intersect(ind,ind_avg);
  end;

  tmp_parms = parms;
  tmp_parms.datatypes = parms.datatypes(d);
  args = mmil_parms2args(tmp_parms,tags);

  MMIL_Combine_Retinotopy_Exam(RootDirs,StudyInfo(ind),args{:});
end;

