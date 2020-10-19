function MMIL_Analyze_Long_Exam(RootDirs,StudyInfo,varargin)
%function MMIL_Analyze_Long_Exam(RootDirs,StudyInfo,[options])
%
% Usage:
%   MMIL_Analyze_Long_Exam(RootDirs,StudyInfo,'key',value,...)
%
% Required Input:
%   RootDirs: RootDirs struct containing full paths of project data directories
%   StudyInfo: StudyInfo struct containing all visits for a single subject
%
% Optional Parameters:
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1[ whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere (can be vector)
%     {default = 0}
%   'mask_midbrain_flag', [0|1] whether to mask out mid brain and other
%     cortical regions marked "unknown" (masking done before smoothing)
%     {default = 0}
%   'outdir': where to place output files
%     can be full path, otherwise will be relative to longitudinal Container
%     {default = 'analysis'}
%   'aseg_flag': use aseg ROIs (non-cortical, volume segmentation)
%     {default = 1}
%   'aseg_roigroups_flag': use groups of aseg roi codes defined by
%      fs_aseg_roigroups
%       (includes 'WholeBrain', 'LatVentricles', and 'AllVentricles')
%     {default = 1}
%   'subhippo_flag': use subdivided hippocampal ROIs
%     {default = 1}
%   'aseg_roigroups': struct array containing the following fields:
%      roiname: name of new ROI
%      roicode: new ROI code number
%      roicodes: vector of aseg ROI code numbers
%     {default: aseg_roigroups = fs_define_aseg_roigroups}
%       (includes 'WholeBrain', 'LatVentricles', and 'AllVentricles')
%   'subhippo_flag': subdivided hippocampal ROIs
%     {default = 1}
%   'aparc_flag': use aparc ROIs (cortical surface parcellation)
%     {default = 1}
%   'erode_flag': create and use eroded ROIs (aseg, groups, hippo)
%     {default = 1}
%   'baseflag': [0|1] assume all visits were registered to visit 1
%      otherwise, registered all visits to all earlier visits
%     {default: 1}
%   'nobiasflag':[0|1] register T1 using Quarc with no bias (forward and reverse)
%      {default: 1}
%   'forceflag': overwrite existing output
%     {default: 0}
%
% Created:  03/01/10 by Don Hagler
% Last Mod: 03/13/14 by Don Hagler
%

%% todo: aparc_roigroups for lobar analysis
%% todo: use mmil_cortsurf_analysis
%% todo: use mmil_aseg_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'projdist',1,[-5,5],...
  'projfrac',[],[],...
  'sphere_flag',false,[false true],...
  'smoothsteps',0,[0,Inf],...
  'sphsmoothsteps',0,[0,Inf],...
  'mask_midbrain_flag',false,[false true],...
  'outdir','analysis',[],...
  'aseg_flag',true,[false true],...
  'aparc_flag',true,[false true],...
  'aseg_roigroups_flag',true,[false true],...
  'aseg_roigroups',[],[],...
  'subhippo_flag',true,[false true],...
  'erode_flag',true,[false true],...
  'baseflag',true,[false true],...
  'nobiasflag',true,[false true],...
  'forceflag',false,[false true],...
... % undocumented:
  'baseline_outdir','analysis',[],... % in FSContainerPath
  'dirprefix','LONG',[],...
  'erode_nvoxels',1,[1:100],...
  'hemilist',{'lh','rh'},{'lh' 'rh'},...
  'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,10001:10003,10011:10016,20001:20003,20009,20010],[1,Inf],...
  'hippo_code_base',10010,[],...
... % regdirs
  'vol_regtypes',{'Fine','ROI'},[],...
  'surf_regtypes',{'Fine'},[],...
  'hippo_regtypes',{'ROI'},[],...
  'group_regtypes',{'Fine'},[],...
...
  'analyze_mri_tags',{'smoothsteps' 'sphsmoothsteps' 'mask_midbrain_flag'...
                      'aseg_roigroups' 'forceflag' 'aseg_roigroups_flag'...
                      'subhippo_flag' 'erode_flag' 'sphere_flag'...
                      'erode_nvoxels' 'hemilist' 'aseg_roilist' 'outdir'},[],...
});

if parms.erode_flag
  parms.erode_outfix = 'erode';
  if parms.erode_nvoxels>1
    parms.erode_outfix = sprintf('%s%d',parms.erode_outfix,parms.erode_nvoxels);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that study info is for single subject
SubjID = unique({StudyInfo.SubjID});
if length(SubjID) > 1
  error('multiple subjects not allowed in StudyInfo');
else
  SubjID = SubjID{1};
end;

% sort studies by VisitNumber
ind = find(cell2mat({StudyInfo.VisitNumber})>=1);
if isempty(ind)
  error('no visits with valid VisitNumber');
end;
StudyInfo = StudyInfo(ind);
VisitNumbers = cell2mat({StudyInfo.VisitNumber});
[VisitNumbers,ind] = sort(VisitNumbers);
StudyInfo = StudyInfo(ind);
if min(VisitNumbers)>1
  error('no VisitNumber=1');
end;

% set longitudinal container path
ContainerPath = sprintf('%s/%s_%s',RootDirs.long,parms.dirprefix,SubjID);
if ~exist(ContainerPath,'file')
  error('Container %s not found',ContainerPath);
end;

% set FreeSurfer container path
ind_baseline = find(cell2mat({StudyInfo.VisitNumber})==1);
subjdir  = RootDirs.fsurf;
baseline_subj = StudyInfo(ind_baseline).fsurf;
FSContainerPath = [subjdir '/' baseline_subj];
if ~exist(FSContainerPath,'file')
  error('FSContainerPath %s not found',FSContainerPath);
end;

% check status of FreeSurfer recon
[status,message] = MMIL_Get_FSReconStatus(FSContainerPath);
if ~ismember(status,[2,5,6])
  fprintf('%s: WARNING: incomplete recon for %s\n',mfilename,baseline_subj);
  return;
end;
if status==6
  fprintf('%s: WARNING: only volume recon is complete for %s\n',...
    mfilename,baseline_subj);
  parms.aparc_flag = 0;
end;

if mmil_isrelative(parms.outdir)  
  root_outdir = [ContainerPath '/' parms.outdir];
else
  root_outdir = parms.outdir;
end;
mmil_mkdir(root_outdir);

if isempty(parms.aseg_roigroups)
  parms.aseg_roigroups = fs_define_aseg_roigroups;
end;

% check files necessary for surface analyses
for h=1:length(parms.hemilist)
  if parms.sphere_flag
    spherefile = sprintf('%s/surf/%s.sphere.reg',...
      FSContainerPath,parms.hemilist{h});
    if ~exist(spherefile,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,spherefile);
      parms.aparc_flag = 0;
      continue;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze baseline recon
%   includes eroding aseg, creating hippocampal ROIs, ROI groups
args = mmil_parms2args(parms,parms.analyze_mri_tags);
MMIL_Analyze_MRI_Exam(FSContainerPath,args{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over pairs of sessions
if parms.baseflag
  nbase = 1;
else
  nbase = length(StudyInfo);
end;
for i=1:nbase
  vA = StudyInfo(i).VisitNumber;
  sessA = StudyInfo(i).SessID;
  dirA = sprintf('%s/visit%d',ContainerPath,vA);
  if ~exist(dirA,'dir')
    fprintf('%s: WARNING: dir %s not found',mfilename,dirA);
    continue;
  end;
  for j=i+1:length(StudyInfo)
    vB = StudyInfo(j).VisitNumber;
    sessB = StudyInfo(j).SessID;
    dirB = sprintf('%s/visit%d',ContainerPath,vB);
    if ~exist(dirB,'dir')      
      fprintf('%s: WARNING: dir %s not found\n',mfilename,dirB);
      continue;
    end;

    % create outdir for analysis
    outdir = sprintf('%s/visit_%d_VS_%d',root_outdir,vA,vB);
    mmil_mkdir(outdir);

    % set nonlinreg dir
    if parms.nobiasflag
      root_regdir = sprintf('%s/nonlinreg_visit%d/f2b',dirB,vA);
    else
      root_regdir = sprintf('%s/nonlinreg_visit%d',dirB,vA);
    end;
    if ~exist(root_regdir,'dir')
      fprintf('%s: WARNING: dir %s not found\n',mfilename,root_regdir);
      continue;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parms.aseg_flag
      regdirs = [];
      fstems = [];
      for r=1:length(parms.vol_regtypes)
        regdirs{r} = sprintf('%s/nonlinReg%s',root_regdir,parms.vol_regtypes{r});
        fstems{r} = sprintf('%s/dv_%s_aseg',outdir,parms.vol_regtypes{r});
        if ~exist(regdirs{r},'dir')
          fprintf('%s: WARNING: dir %s not found\n',mfilename,regdirs{r});
          continue;
        end;
      end;

      % warp dv to volume atlas
      if parms.nobiasflag
        fname_in = sprintf('%s/dvAvg.mgz',regdirs{1});
      else
        fname_in = sprintf('%s/dv.mgz',regdirs{1});
      end;
      fname_out = sprintf('%s/dv_warp.mgz',regdirs{1});
      if exist(fname_in,'file') && (~exist(fname_out,'file') || parms.forceflag)
        fs_warp_vol2atlas(baseline_subj,fname_in,fname_out,...
          'subjdir',subjdir,'overwrite_flag',parms.forceflag);
      end

      fprintf('%s: extracting dv for volume ROIs...\n',mfilename);
      for r=1:length(regdirs)
        regdir = regdirs{r};
        fstem = fstems{r};
        if parms.nobiasflag
          fname_dv = sprintf('%s/dvAvg.mgz',regdir);
        else
          fname_dv = sprintf('%s/dv.mgz',regdir);
        end;
        if ~exist(regdir,'dir'), continue; end;
        if ~exist(fname_dv,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname_dv);
          continue;
        end;

        % extract dv measures for ROIs
        if ~parms.erode_flag
          fname_out = [fstem '.mat'];
          fname_aseg = sprintf('%s/mri/aseg.mgz',FSContainerPath);
        else
          fname_out = sprintf('%s_%s.mat',fstem,parms.erode_outfix);
          fname_aseg = sprintf('%s/%s/aseg_%s.mgz',...
            FSContainerPath,parms.baseline_outdir,parms.erode_outfix);
        end;            
        if ~exist(fname_aseg,'file')
          fprintf('%s: WARNING: aseg file %s not found\n',...
            mfilename,fname_aseg);
          continue;
        end;
        if ~exist(fname_out,'file') || parms.forceflag
          roi_data = mmil_aseg_roi(fname_dv,fname_aseg,'minval',0);
          if ~isempty(roi_data)
            save(fname_out,'roi_data','baseline_subj');
          else
            fprintf('%s: WARNING: failed to extract dv ROI averages\n',...
              mfilename);
          end;
        end;
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract dv measures for ROI groups
    if parms.aseg_roigroups_flag
      regdirs = [];
      fstems = [];
      for r=1:length(parms.group_regtypes)
        regdirs{r} = sprintf('%s/nonlinReg%s',root_regdir,parms.group_regtypes{r});
        fstems{r} = sprintf('%s/dv_%s_roigroups',outdir,parms.group_regtypes{r});
        if ~exist(regdirs{r},'dir')
          fprintf('%s: WARNING: dir %s not found\n',mfilename,regdirs{r});
          continue;
        end;
      end;

      fprintf('%s: extracting dv for ROI groups...\n',mfilename);
      for r=1:length(regdirs)
        regdir = regdirs{r};
        fstem = fstems{r};
        if parms.nobiasflag
          fname_dv = sprintf('%s/dvAvg.mgz',regdir);
        else
          fname_dv = sprintf('%s/dv.mgz',regdir);
        end;
        if ~exist(regdir,'dir'), continue; end;
        if ~exist(fname_dv,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname_dv);
          continue;
        end;
        if ~parms.erode_flag
          fname_out = [fstem '.mat'];
        else
          fname_out = sprintf('%s_%s.mat',fstem,parms.erode_outfix);
        end;
        rois_exist = 1;
        for k=1:length(parms.aseg_roigroups)
          fname_aseg = sprintf('%s/%s/aseg_roigroup_%d',...
            FSContainerPath,parms.baseline_outdir,...
            parms.aseg_roigroups(k).roicode);
          if parms.erode_flag
            fname_aseg = sprintf('%s_%s',fname_aseg,parms.erode_outfix);
          end;
          fname_aseg = [fname_aseg '.mgz'];
          if ~exist(fname_aseg,'file')
            fprintf('%s: WARNING: aseg ROI file %s not found\n',...
              mfilename,fname_aseg);
            rois_exist = 0;
            break;
          end;
        end;
        if rois_exist && (~exist(fname_out,'file') || parms.forceflag)
          roi_data = [];
          for k=1:length(parms.aseg_roigroups)
            fname_aseg = sprintf('%s/%s/aseg_roigroup_%d',...
              FSContainerPath,parms.baseline_outdir,...
              parms.aseg_roigroups(k).roicode);
            if parms.erode_flag
              fname_aseg = sprintf('%s_%s',fname_aseg,parms.erode_outfix);
            end;
            fname_aseg = [fname_aseg '.mgz'];
            tmp_roi_data = mmil_aseg_roi(fname_dv,fname_aseg,'minval',0,...
              'roigroups',parms.aseg_roigroups(k)); % calculate change for this ROI group
            roi_data = [roi_data tmp_roi_data]; % append it to the list of individual ROI's
          end
          if ~isempty(roi_data)
            save(fname_out,'roi_data','baseline_subj');
          else
            fprintf('%s: WARNING: failed to extract dv ROI averages\n',...
              mfilename);
          end;
        end;
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract dv for sub-sections of hippocampus
    if parms.subhippo_flag
      flist_hippo_masks = fs_create_subhippo_masks(...
        baseline_subj,RootDirs.fsurf,parms.baseline_outdir,1);
      if isempty(flist_hippo_masks)
        fprintf('%s: WARNING: fs_create_subhippo_masks failed\n',mfilename);
      end;
    else
      flist_hippo_masks = [];
    end;
    if ~isempty(flist_hippo_masks)
      regdirs = [];
      fstems = [];
      pstems = [];
      if parms.erode_flag
        for k=1:length(flist_hippo_masks)
          fname_in = flist_hippo_masks{k};
          [tmp_fpath,tmp_fstem,tmp_fext]=fileparts(flist_hippo_masks{k});
          fname_out = sprintf('%s/%s/%s_%s%s',...
            FSContainerPath,parms.baseline_outdir,...
            tmp_fstem,parms.erode_outfix,tmp_fext);
          flist_hippo_masks{k} = fname_out;
        end;
      end;
      for r=1:length(parms.hippo_regtypes)
        regdirs{r} = sprintf('%s/nonlinReg%s',root_regdir,parms.hippo_regtypes{r});
        fstems{r} = sprintf('%s/dv_%s_hippo',outdir,parms.hippo_regtypes{r});
        if ~exist(regdirs{r},'dir')
          fprintf('%s: WARNING: dir %s not found\n',mfilename,regdirs{r});
          continue;
        end;
      end;

      fprintf('%s: extracting dv for hippocampal ROIs...\n',mfilename);
      for r=1:length(regdirs)
        regdir = regdirs{r};
        fstem = fstems{r};
        if parms.erode_flag
          fname_out = sprintf('%s_%s',fstem,parms.erode_outfix);
        else
          fname_out = fstem;
        end;
        fname_out = [fname_out '.mat'];
        if parms.nobiasflag
          fname_dv = sprintf('%s/dvAvg.mgz',regdir);
        else
          fname_dv = sprintf('%s/dv.mgz',regdir);
        end;
        if ~exist(regdir,'dir'), continue; end;
        if ~exist(fname_dv,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname_dv);
          continue;
        end;
        % extract ROI averages for each hippocampus mask
        if ~exist(fname_out,'file') || parms.forceflag
          fprintf('%s: getting hippo ROI data from %s...\n',mfilename,fname_dv);
          roi_data = mmil_multi_roi(fname_dv,flist_hippo_masks);
          if isempty(roi_data)
            fprintf('%s: WARNING: hippo ROI analysis failed\n',mfilename);
          else
            % add roi names
            for k=1:length(flist_hippo_masks)
              [tmp_fpath,tmp_fstem,tmp_fext]=fileparts(flist_hippo_masks{k});
              roi_data(k).roiname = tmp_fstem;
              roi_data(k).roicode = parms.hippo_code_base + k;
            end;
            % save results
            save(fname_out,'roi_data');
          end;
        end;
      end
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % longitudinal dv on (baseline) surface
    if parms.aparc_flag
      clear regdirs fstems;
      regdirs = [];
      fstems = [];
      pstems = [];
      for r=1:length(parms.surf_regtypes)
        regdirs{r} = sprintf('%s/nonlinReg%s',root_regdir,parms.surf_regtypes{r});
        fstems{r} = sprintf('%s/dv_%s_aparc',outdir,parms.surf_regtypes{r});
        pstems{r} = sprintf('%s/dv_%s_paint',outdir,parms.surf_regtypes{r});
        if ~exist(regdirs{r},'dir')
          fprintf('%s: WARNING: dir %s not found\n',mfilename,regdirs{r});
          continue;
        end;
      end;

      fprintf('%s: sampling dv onto baseline surface...\n',mfilename);
      for r=1:length(regdirs)
        regdir = regdirs{r};
        fstem = fstems{r};
        pstem = pstems{r};
        if parms.nobiasflag
          fname_dv = sprintf('%s/dvAvg.mgz',regdir);
        else
          fname_dv = sprintf('%s/dv.mgz',regdir);
        end;
        if ~exist(regdir,'dir'), continue; end;
        if ~exist(fname_dv,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname_dv);
          continue;
        end;
        % paint onto surface
        if ~isempty(parms.projfrac)
          projfrac_flag = 1;
        else
          projfrac_flag = 0;
        end;
        fs_paint(baseline_subj,fname_dv,...
          'outtype','mgz','sphere_flag',0,...
          'outstem',pstem,...
          'subjdir',subjdir,...
          'projdist',parms.projdist,...
          'projfrac',parms.projfrac,...
          'projfrac_flag',projfrac_flag,...
          'mask_midbrain_flag',parms.mask_midbrain_flag,...
          'smoothsteps',0,...
          'overwrite_flag',parms.forceflag);
        if parms.sphere_flag
          for i=1:length(parms.sphsmoothsteps)
            sphsmoothsteps = parms.sphsmoothsteps(i);
            % paint onto ico sphere surface, smooth on ico
            fs_paint(baseline_subj,fname_dv,...
              'outtype','mgz','sphere_flag',1,...
              'outstem',pstem,...
              'subjdir',subjdir,...
              'projdist',parms.projdist,...
              'projfrac',parms.projfrac,...
              'projfrac_flag',projfrac_flag,...
              'mask_midbrain_flag',parms.mask_midbrain_flag,...
              'sphsmoothsteps',sphsmoothsteps,...
              'overwrite_flag',parms.forceflag);
          end;
        end;
      end;

      fprintf('%s: extracting dv for surface ROIs...\n',mfilename);
      for r=1:length(regdirs)
        regdir = regdirs{r};
        fstem = fstems{r};
        pstem = pstems{r};
        % get stats from surface ROIs
        fname_out = [fstem '.mat'];
        if ~exist(regdir,'dir'), continue; end;
        if ~exist(fname_out,'file') || parms.forceflag
          roi_data = [];
          for h=1:length(parms.hemilist)
            hemi = parms.hemilist{h};
            fname_aparc = sprintf('%s/label/%s.aparc.annot',...
              FSContainerPath,hemi);
            if ~exist(fname_aparc,'file')
              error('%s not found',fname_aparc);
            end;
            surffuncname = [pstem '-' hemi '.mgz'];
            fprintf('%s: getting longitudinal aparc ROI data...\n',mfilename);
            tmp_roi_data = mmil_surf_roi(surffuncname,...
              'fname_aparc',fname_aparc,'minval',0);
            if isempty(tmp_roi_data)
              error('failed to get aparc ROI data');
            end;
            roi_data = [roi_data,tmp_roi_data];
          end;
          save(fname_out,'roi_data','baseline_subj');
        end;
      end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: finished\n',mfilename);
