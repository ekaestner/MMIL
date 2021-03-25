function MMIL_Combine_Retinotopy_Exam(RootDirs,StudyInfo,varargin)
%function MMIL_Combine_Retinotopy_Exam(RootDirs,StudyInfo,[options])
%
% Purpose: Combine Fourier analysis results from multiple sessions
%   for both polar angle and eccentricity retinotopy
%
% Required Input:
%   RootDirs: RootDirs struct containing full paths of project data directories
%   StudyInfo: StudyInfo struct containing all visits for a single subject
%
% Optional Paramters:
%  'infix': string attached to processed BOLD data files
%    e.g. 'corr_resBOLD'
%    {default = []}
%  'multi_flag': [0|1|2] indicates source of input analyses
%    0: use results from individual scans
%    1: use results from multiple scans within session
%    2: use results from multiple sessions
%    {default = 0}
%
% Optional Parameters specifying how Fourier analysis was done:
%  'fstats_type' - [0|1|2] how output Fourier components were scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default: 2}
%  'resamp_flag' - [0|1] whether Fourier results were resampled to 1x1x1
%     before painting
%    {default: 0}
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 0}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default = 0}
%   sphere_flag': [0|1] if data were resampled to spherical atlas space
%    {default = 0}
%
% Optional Parameters for viewing results:
%  'tksmooth' - number of surface smoothing iterations for tksurfer script
%    (relevant only if paint_flag = 1)
%    {default: 2}
%  'fthresh' - threshold for scaling amplitudes in tksurfer
%    {default: 0}
%  'fmid' - mid point for scaling amplitudes in tksurfer
%    {default: 1.5}
%  'fslope' - slope for scaling amplitudes in tksurfer
%    {default: 1.0}
%  'view' - view of cortical surface view in tksurfer
%    ('med','ven','lat','pos','dor')
%    {default: 'med'}
%  'surf' - surface to display in tksurfer
%    {default: 'inflated'}
%
% Other optional parameters:
%  'r_min_factor': value multiplied by r_max to get r_min (if empty)
%    {default = 0.02}
%  'cxfstatsflag' [0|1] whether to calculate cross-scan complex f-stats
%    (in addition to amplitude f-stats)
%    { default: 1 }
%  'OutContainerPath' - output ContainerPath
%    if empty, will use first proc_bold Container with polar retinotopy data
%    or first proc_bold Container with eccen retinotopy data if no polar data
%    {default = []}
%  'forceflag' - [0|1] whether to run calculations even if output files exist
%    {default: 0}
%
% Created:  03/03/11 by Don Hagler
% Last Mod: 03/11/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'infix',[],[],...
  'multi_flag',0,[0:2],...
  'forceflag',false,[false true],...
...
  'fstats_type',2,[0:2],...
  'resamp_flag',false,[false true],...
  'smoothsteps',0,[0,1000],...
  'sphsmoothsteps',0,[0,1000],...
  'sphere_flag',false,[false true],...
... % parameters for viewing results
  'tksmooth',2,[0,1000],...
  'fthresh',0,[0,100],...
  'fmid',1.5,[0.1,100],...
  'fslope',1.0,[0.1,100],...
  'view','med',{'lat','med','ven','pos','dor'},...
  'surf','inflated',[],...
... % other
  'r_min_factor',0.02,[],...
  'cxfstatsflag',true,[false true],...
  'OutContainerPath',[],[],...
...
  'datatypes',{'polar','eccen'},[],...
  'fsaverage','fsaverage',[],... % used if sphere_flag = 1
  'fsavg_rootdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
});

% parameters to come from StudyInfo
SI_tags = {...
  'snums'...
  'revflags'...
  'phase_offset'...
  'r_max'...
  'r_min'...
  'logtrans_flag'...
};
% parameters to come from StudyInfo that are cell arrays
SI_cell_tags = {...
  'snums'...
  'revflags'...
};

% parameters to not  pass to BOLD_MMIL_Combine_Fourier_Sessions
tags_excl = {'datatypes'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that study info is for single subject
SubjID = unique({StudyInfo.SubjID});
if ~parms.sphere_flag && length(SubjID) > 1
  error('multiple subjects not allowed in StudyInfo unless sphere_flag = 1');
end;

% find pol and ecc sessions
ind_pol = find(~cellfun(@isempty,{StudyInfo.FP_snums}));
ind_ecc = find(~cellfun(@isempty,{StudyInfo.FE_snums}));

if isempty(ind_pol) && isempty(ind_ecc)
  fprintf('%s: WARNING: no retintotopy data (FP_snums and FE_snums) in StudyInfo\n',...
    mfilename);
  return;
end;

if ~parms.sphere_flag
  nsess = length(StudyInfo);
  FSContainerPaths = cell(nsess,1);
  for i=1:nsess
    FSContainerPaths{i} = [RootDirs.fsurf '/' StudyInfo(i).fsurf];
  end;
  [tmp,ind_FS] = unique(FSContainerPaths,'first');
  if ~parms.sphere_flag && length(ind_FS)>1
    fprintf('%s: ERROR: multiple FreeSurfer recons not allowed unless sphere_flag = 1\n',...
      mfilename);
    return;
  end;
  FSContainerPath = FSContainerPaths{1};
else
  % will be set to fsaverage by BOLD_MMIL_Combine_Fourier_Sessions
  FSContainerPath = []; 
end;

% combine sessions
for d=1:length(parms.datatypes)
  datatype = parms.datatypes{d};
  switch datatype
    case 'polar'
      ind = ind_pol;
      tag_prefix = 'FP';
    case 'eccen'
      ind = ind_ecc;
      tag_prefix = 'FE';
  end;
  if isempty(ind)
    fprintf('%s: WARNING: no %s retintotopy data\n',...
      mfilename,datatype);
    continue;
  end;
  ContainerDirs = {StudyInfo(ind).proc_bold};
  nsess = length(ContainerDirs);
  ContainerPaths = cell(nsess,1);
  for i=1:nsess
    ContainerPaths{i} = [RootDirs.proc_bold '/' ContainerDirs{i}];
  end;
  if isempty(parms.OutContainerPath)
    parms.OutContainerPath = ContainerPaths{1};
  end;
  for f=1:length(SI_tags)
    out_tag = SI_tags{f};
    in_tag = sprintf('%s_%s',tag_prefix,out_tag);
    if isfield(StudyInfo,in_tag)
      if ismember(out_tag,SI_cell_tags)
        parms.(out_tag) = {StudyInfo(ind).(in_tag)};
      else
        parms.(out_tag) = [StudyInfo(ind).(in_tag)];
      end;
    else
      parms.(out_tag) = [];
    end;
  end;
  parms.datatype = datatype;
  %% todo: option to specify outstem/outdir
  parms.outdir = ['BOLD_multisess_' datatype '_analysis'];
  parms.outstem = ['BOLD_multisess_' datatype];

  tags = setdiff(fieldnames(parms),tags_excl);
  args = mmil_parms2args(parms,tags);
  BOLD_MMIL_Combine_Fourier_Sessions(ContainerPaths,FSContainerPath,args{:})
end;

