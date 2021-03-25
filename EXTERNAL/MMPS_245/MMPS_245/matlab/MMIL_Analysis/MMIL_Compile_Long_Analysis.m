function results = MMIL_Compile_Long_Analysis(ProjID,varargin)
% function results = MMIL_Compile_Long_Analysis(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
% 
% Optional Parameters that specify project specific information
%  'RootDirs':
%    struct containing the following fields: proc, fsurf, long
%    these specify the locations of data
%  'StudyInfo': struct array of study information
%    (e.g. read from csv file with MMIL_Read_StudyInfo)
%    must contain these fields: SubjID, StudyDate, VisitNumber
%    may contain these fields: proc, fsurf
%    if proc and fsurf are unspecified, will look for Containers
%      with SubjID and StudyDate (will choose first one if more than one)
%    if empty, use all subjects found in RootDirs.proc and RootDirs.fsurf
%    {default = []}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%    {default = 1}
%
% Optional Parameters Controlling which Stats to Compile:
%  'outdir': where analyzed mat files can be found relative to ContainerDir
%     {default = 'analysis'}
%  'regtype': nonlinear registration type (e.g. 'Fine', 'ROI')
%     {default = 'Fine'}
%  'aseg_flag': use aseg ROIs (non-cortical, volume segmentation)
%     {default = 1}
%  'aseg_roigroups_flag': [0|1] use masks for groups of aseg roi codes
%     includes 'WholeBrain', 'LatVentricles', and 'AllVentricles'
%     {default = 1}
%  'subhippo_flag': use subdivided hippocampal ROIs
%     {default = 1}
%  'aparc_flag': use aparc ROIs (cortical surface parcellation)
%     {default = 1}
%
% Output:
%   results: struct containing these fields:
%     SubjIDs: cell array of subject IDs
%     Cort: matrix of cortical ROI volume change data
%     Cort_labels: cell array of row labels for cortical ROI data
%     Vol: matrix of volume ROI volume change data
%     Vol_labels: cell array of row labels for volume ROI data
%
%
%   NOTE: subjects with incomplete FreeSurfer recons will be excluded
%    from output.  NaN's will be inserted in the case of missing ROI data.
%
% Created:  03/02/10 by Don Hagler
% Last Mod: 03/26/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(ProjID,varargin);

% initialize output struct
results = init_results(parms);

% compile aseg results
if parms.aseg_flag
  results = compile_aseg_results(parms,results);
end;

% compile aparc results
if parms.aparc_flag
  results = compile_aparc_results(parms,results);
end;

% compile hippocampus sub-section results
if parms.subhippo_flag
  results = compile_subhippo_results(parms,results);
end;

% concatenate aseg, aparc, and subhippo results
results = concat_results(parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'RootDirs',[],[],...
    'StudyInfo',[],[],...
    'qcflag',true,[false true],...
  ...
    'outdir','analysis',[],...
    'baseflag',true,[false true],...
    'regtype','Fine',[],...
    'aseg_flag',true,[false true],...
    'subhippo_flag',true,[false true],...
    'aseg_roigroups_flag',true,[false true],...
    'aparc_flag',true,[false true],...
  ...
    'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,10001:10003,20001:20003,20009,20010],[1,Inf],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'hippo_roilist',[10011:10016],[1,Inf],...
    'fname_fscolorlut',[],[],...
    'dirprefix','LONG',[],...
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
  ...
    'required_containers',{'proc'},[],...
    'required_rootdirs',{'proc','fsurf','long'},[],...
    'modality','MRI',[],...
  };
  parms = mmil_args2parms(options,parms_filter);

  % filter StudyInfo or generate StudyInfo struct if none supplied
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(StudyInfo), error('empty StudyInfo'); end;

  if ~isempty(ProjInfo)
    % For arg names present in both varargin and ProjInfo
    % the varargin values will appear in merged_args
    ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
    if isfield(ProjInfo,'LONG_baseflag')
      parms.baseflag = ProjInfo.LONG_baseflag;
    end;
  end;

  parms.StudyInfo = StudyInfo;
  parms.RootDirs = RootDirs;

  % set erode_outfix
  if parms.erode_flag
    parms.erode_outfix = 'erode';
    if parms.erode_nvoxels>1
      parms.erode_outfix = sprintf('%s%d',parms.erode_outfix,parms.erode_nvoxels);
    end;
  end;

  % get ROI codes and names
  parms = get_roi_codes(parms);

  % check visit information
  parms = check_visits(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = get_roi_codes(parms)
  % read roi codes and names from FreeSurfer color lookup table
  if isempty(parms.fname_fscolorlut)
    MMPS_parms = getenv('MMPS_PARMS');
    parms.fname_fscolorlut = [MMPS_parms '/MMIL_FSColorLUT.txt'];
  end;
  if ~exist(parms.fname_fscolorlut,'file')
    error('FreeSurfer color lookup table file %s not found',...
      parms.fname_fscolorlut);
  end;

  [all_roicodes,all_roinames] = fs_colorlut(parms.fname_fscolorlut);

  % get aseg roi codes and names
  [tmp,ind_aseg]=intersect(all_roicodes,parms.aseg_roilist);
  parms.aseg_roicodes = all_roicodes(ind_aseg);
  parms.aseg_roinames = all_roinames(ind_aseg);
  parms.aseg_nrois = length(ind_aseg);

  % get aparc roi codes and names
  [tmp,ind_aparc]=intersect(all_roicodes,parms.aparc_roilist);
  parms.aparc_roicodes = all_roicodes(ind_aparc);
  parms.aparc_roinames = all_roinames(ind_aparc);
  parms.aparc_nrois = length(ind_aparc);

  % get hippo roi codes and names
  [tmp,ind_hippo]=intersect(all_roicodes,parms.hippo_roilist);
  parms.hippo_roicodes = all_roicodes(ind_hippo);
  parms.hippo_roinames = all_roinames(ind_hippo);
  parms.hippo_nrois = length(ind_hippo);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_visits(parms)
  % identify subjects with multiple visits
  if ~isfield(parms.StudyInfo,'VisitNumber')
    [parms.StudyInfo(:).VisitNumber] = deal(0);
  end
  visitnums = cell2mat({parms.StudyInfo.VisitNumber});
  ind = find(visitnums>1);
  parms.SubjIDs = unique(({parms.StudyInfo(ind).SubjID}));
  parms.nsubs = length(parms.SubjIDs);
  if parms.nsubs==0
    error('no subjects with followup scans found -- check parms.StudyInfo VisitNumber');
  end;

  % exclude empty SubjIDs
  ind = find(~cellfun(@isempty,parms.SubjIDs));
  if isempty(ind)
    error('SubjIDs are all empty');
  end;
  parms.SubjIDs = parms.SubjIDs(ind);

  % find studies for valid SubjIDs
  ind = find(ismember({parms.StudyInfo.SubjID},parms.SubjIDs));
  parms.StudyInfo = parms.StudyInfo(ind);

  % check for baseline scan
  excl_SubjIDs = [];
  for s=1:parms.nsubs
    ind = find(strcmp({parms.StudyInfo.SubjID},parms.SubjIDs{s}));
    if ~any([parms.StudyInfo(ind).VisitNumber]==1)
      fprintf('%s: WARNING: subject %s has no visit 1 - skipping\n',...
        mfilename,parms.SubjIDs{s});
      excl_SubjIDs{end+1} = parms.SubjIDs{s};
    end;
  end;

  % check baseline parms.StudyInfo.fsurf
  for i=1:length(parms.StudyInfo)
    SubjID = parms.StudyInfo(i).SubjID;
    if parms.StudyInfo(i).VisitNumber==1 
      FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(i).fsurf];
      % check exists
      if isempty(parms.StudyInfo(i).fsurf) || ~exist(FSContainerPath,'dir')
        fprintf('%s: WARNING: subj %s missing FreeSurfer recon for visit 1\n',...
          mfilename,parms.StudyInfo(i).SubjID);
        excl_SubjIDs{end+1} = parms.StudyInfo(i).SubjID;
        continue;
      end;
      % check status
      if parms.aparc_flag || parms.aseg_flag
        [status,message] = MMIL_Get_FSReconStatus(FSContainerPath);
        if (parms.aparc_flag & ~ismember(status,[2,5])) | ...
           (~parms.aseg_flag & ~ismember(status,[2,5,6]))
          fprintf('%s: WARNING: incomplete recon for %s\n',mfilename,parms.StudyInfo(i).SubjID);
            excl_SubjIDs{end+1} = parms.StudyInfo(i).SubjID;
        end;
      end;
    end;
  end;

  % exclude subjects
  parms.SubjIDs = setdiff(parms.SubjIDs,excl_SubjIDs);
  parms.nsubs = length(parms.SubjIDs);
  if parms.nsubs==0
    error('no subjects with baseline recons found');
  end;
  ind = find(ismember({parms.StudyInfo.SubjID},parms.SubjIDs));
  parms.StudyInfo = parms.StudyInfo(ind);

  % check longitudinal containers
  excl_SubjIDs = [];
  for i=1:length(parms.SubjIDs)
    ContainerPath = sprintf('%s/%s_%s',...
      parms.RootDirs.long,parms.dirprefix,parms.SubjIDs{i});
    if ~exist(ContainerPath,'file')
      fprintf('%s: WARNING: container %s not found\n',mfilename,ContainerPath);
      excl_SubjIDs{end+1} = SubjID;
    end;
  end;
  parms.SubjIDs = setdiff(parms.SubjIDs,excl_SubjIDs);
  parms.nsubs = length(parms.SubjIDs);
  if parms.nsubs==0
    error('no subjects with longitudinal containers');
  end;
  ind = find(ismember({parms.StudyInfo.SubjID},parms.SubjIDs));
  parms.StudyInfo = parms.StudyInfo(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  % initialize output struct
  results.StudyInfo = parms.StudyInfo;
  results.unique_SubjIDs = parms.SubjIDs; % unique subject IDs
  results.nsubs = parms.nsubs; % unique subject IDs
  results.maxVisitNumber = max(cell2mat({parms.StudyInfo.VisitNumber}));

  % determine size of output matrices
  %   based on number of pairs of visits for each subject
  r = 0;
  results.SubjIDs = {};
  results.VisitNumbers_Baseline = {};
  results.VisitNumbers_Followup = {};
  results.StudyDates_Baseline = {};
  results.StudyDates_Followup = {};
  for s=1:results.nsubs
    SubjID = parms.SubjIDs{s};
    % get StudyInfo for this subject, sorted by VisitNumber
    [SubjInfo,nvisits,nbase] = get_SubjInfo(parms,SubjID);
    % loop over baseline-followup pairs
    for i=1:nbase
      vA = SubjInfo(i).VisitNumber;
      StudyDate_Baseline = SubjInfo(i).StudyDate;
      for j=i+1:nvisits
        vB = SubjInfo(j).VisitNumber;
        StudyDate_Followup = SubjInfo(j).StudyDate;
        r = r + 1;
        results.SubjIDs{r,1} = SubjID;
        results.VisitNumbers_Baseline{r,1} = vA;
        results.VisitNumbers_Followup{r,1} = vB;
        results.StudyDates_Baseline{r,1} = StudyDate_Baseline;
        results.StudyDates_Followup{r,1} = StudyDate_Followup;
      end;
    end;
  end;
  results.nrows = r;

  if parms.aseg_flag
    results.SubCort.nrois = parms.aseg_nrois;
    results.SubCort.data_baseline = nan(results.nrows,parms.aseg_nrois);
    results.SubCort.data_change = nan(results.nrows,parms.aseg_nrois);
    results.SubCort.roicodes = parms.aseg_roicodes;
    results.SubCort.roinames = parms.aseg_roinames;
  end;

  if parms.aparc_flag
    results.Cort.nrois = parms.aparc_nrois;
    results.Cort.data_baseline = nan(results.nrows,parms.aparc_nrois);
    results.Cort.data_change = nan(results.nrows,parms.aparc_nrois);
    results.Cort.roicodes = parms.aparc_roicodes;
    results.Cort.roinames = parms.aparc_roinames;
  end;

  if parms.subhippo_flag
    results.Hippo.nrois = parms.hippo_nrois;
    results.Hippo.data_baseline = nan(results.nrows,parms.hippo_nrois);
    results.Hippo.data_change = nan(results.nrows,parms.hippo_nrois);
    results.Hippo.roicodes = parms.hippo_roicodes;
    results.Hippo.roinames = parms.hippo_roinames;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SubjInfo,nvisits,nbase] = get_SubjInfo(parms,SubjID)
  ind = find(strcmp(SubjID,{parms.StudyInfo.SubjID}));
  nvisits = length(ind);
  SubjInfo = parms.StudyInfo(ind);
  % sort by visit number
  [visitnums,ind] = sort(cell2mat({SubjInfo.VisitNumber}));
  SubjInfo = SubjInfo(ind);
  % select number of baselines
  if parms.baseflag
    nbase = 1;
  else
    nbase = length(visitnums)-1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_aseg_results(parms,results)
  r = 0;
  for s=1:results.nsubs
    SubjID = parms.SubjIDs{s};
    % set longitudinal container path
    ContainerPath = sprintf('%s/%s_%s',...
      parms.RootDirs.long,parms.dirprefix,SubjID);
    % get StudyInfo for this subject, sorted by VisitNumber
    [SubjInfo,nvisits,nbase] = get_SubjInfo(parms,SubjID);

    % get data for each baseline-followup pair
    for i=1:nbase
      vA = SubjInfo(i).VisitNumber;
      for j=i+1:nvisits
        % set analysis dir
        vB = SubjInfo(j).VisitNumber;
        r = r + 1;

        % get baseline data from FreeSurfer recon container
        baseline_subj = SubjInfo(i).fsurf;
        FSContainerPath = [parms.RootDirs.fsurf '/' baseline_subj];
        % load baseline aseg_stats
        matfile = sprintf('%s/%s/aseg_stats.mat',...
          FSContainerPath,parms.outdir);
        aseg_stats = [];
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: missing file %s\n',mfilename,matfile);
        else
          load(matfile);
          if isempty(aseg_stats)
            fprintf('%s: WARNING: empty aseg_stats in %s\n',...
              mfilename,matfile);
          else
            roicodes = cell2mat({aseg_stats.roicode});
            [tmp,ind_all,ind_subj]=intersect(parms.aseg_roicodes,roicodes);
            tmp_data = cell2mat({aseg_stats(ind_subj).volume});
            results.SubCort.data_baseline(r,ind_all) = tmp_data;
          end;
        end;

        % get change aseg roi_data from longitudinal container
        indir = sprintf('%s/%s/visit_%d_VS_%d',ContainerPath,parms.outdir,vA,vB);
        matfile = [indir '/dv_' parms.regtype '_aseg'];
        if parms.erode_flag
          matfile = sprintf('%s_%s',matfile,parms.erode_outfix);
        end;
        matfile = [matfile '.mat'];
        roi_data = [];
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: missing file %s\n',mfilename,matfile);
        else
          load(matfile);
          if isempty(roi_data)
            fprintf('%s: WARNING: empty roi_data in %s\n',...
              mfilename,matfile);
          else
            roicodes = cell2mat({roi_data.roicode});
            [tmp,ind_all,ind_subj]=intersect(parms.aseg_roicodes,roicodes);
            tmp_data = cell2mat({roi_data(ind_subj).avg});
            results.SubCort.data_change(r,ind_all) = tmp_data;
          end;
        end;

        % load aseg groups roi_data
        if parms.aseg_roigroups_flag
          matfile = [indir '/dv_' parms.regtype '_roigroups'];
          if parms.erode_flag
            matfile = sprintf('%s_%s',matfile,parms.erode_outfix);
          end;
          matfile = [matfile '.mat'];
          roi_data = [];
          if ~exist(matfile,'file')
            fprintf('%s: WARNING: missing file %s\n',mfilename,matfile);
          else
            load(matfile);
            if isempty(roi_data)
              fprintf('%s: WARNING: empty roi_data in %s\n',...
                mfilename,matfile);
            else
              roicodes = cell2mat({roi_data.roicode});
              [tmp,ind_all,ind_subj]=intersect(parms.aseg_roicodes,roicodes);
              tmp_data = cell2mat({roi_data(ind_subj).avg});
              results.SubCort.data_change(r,ind_all) = tmp_data;
            end;
          end;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_aparc_results(parms,results)
  r = 0;
  for s=1:results.nsubs
    SubjID = parms.SubjIDs{s};
    % set longitudinal container path
    ContainerPath = sprintf('%s/%s_%s',...
      parms.RootDirs.long,parms.dirprefix,SubjID);
    % get StudyInfo for this subject, sorted by VisitNumber
    [SubjInfo,nvisits,nbase] = get_SubjInfo(parms,SubjID);

    % get data for each baseline-followup pair
    for i=1:nbase
      vA = SubjInfo(i).VisitNumber;
      for j=i+1:nvisits
        % set analysis dir
        vB = SubjInfo(j).VisitNumber;
        r = r + 1;

        % get baseline aparc_stats from FreeSurfer recon container
        baseline_subj = SubjInfo(i).fsurf;
        FSContainerPath = [parms.RootDirs.fsurf '/' baseline_subj];
        matfile = sprintf('%s/%s/aparc_stats.mat',FSContainerPath,parms.outdir);
        aparc_stats = [];
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: missing file %s\n',mfilename,matfile);
        else
          load(matfile);
          if isempty(aparc_stats)
            fprintf('%s: WARNING: empty aparc_stats in %s\n',...
              mfilename,matfile);
          else
            roicodes = cell2mat({aparc_stats.roicode});
            [tmp,ind_all,ind_subj]=intersect(parms.aparc_roicodes,roicodes);
            tmp_data_thick = cell2mat({aparc_stats(ind_subj).thickavg});
            tmp_data_area = cell2mat({aparc_stats(ind_subj).surfarea});
            tmp_data = tmp_data_thick.*tmp_data_area;
            results.Cort.data_baseline(r,ind_all) = tmp_data;
          end;
        end;

        % load change aparc roi data from longitudinal container
        indir = sprintf('%s/%s/visit_%d_VS_%d',ContainerPath,parms.outdir,vA,vB);
        matfile = [indir '/dv_' parms.regtype '_aparc.mat'];
        roi_data = [];
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: missing file %s\n',mfilename,matfile);
        else
          load(matfile);
          if isempty(roi_data)
            fprintf('%s: WARNING: empty roi_data in %s\n',...
              mfilename,matfile);
          else
            roicodes = cell2mat({roi_data.roicode});
            [tmp,ind_all,ind_subj]=intersect(parms.aparc_roicodes,roicodes);
            tmp_data = cell2mat({roi_data(ind_subj).avg});
            results.Cort.data_change(r,ind_all) = tmp_data;
          end;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_subhippo_results(parms,results)
  r = 0;
  for s=1:results.nsubs
    SubjID = parms.SubjIDs{s};
    % set longitudinal container path
    ContainerPath = sprintf('%s/%s_%s',...
      parms.RootDirs.long,parms.dirprefix,SubjID);
    % get StudyInfo for this subject, sorted by VisitNumber
    [SubjInfo,nvisits,nbase] = get_SubjInfo(parms,SubjID);

    % get data for each baseline-followup pair
    for i=1:nbase
      vA = SubjInfo(i).VisitNumber;
      for j=i+1:nvisits
        % set analysis dir
        vB = SubjInfo(j).VisitNumber;
        r = r + 1;

        % load change subhippo roi data from longitudinal container
        indir = sprintf('%s/%s/visit_%d_VS_%d',ContainerPath,parms.outdir,vA,vB);
        matfile = [indir '/dv_' parms.regtype '_hippo'];
        if parms.erode_flag
          matfile = sprintf('%s_%s',matfile,parms.erode_outfix);
        end;
        matfile = [matfile '.mat'];
        roi_data = [];
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: missing file %s\n',mfilename,matfile);
        else
          load(matfile);
          if isempty(roi_data)
            fprintf('%s: WARNING: empty roi_data in %s\n',...
              mfilename,matfile);
          else
            roicodes = cell2mat({roi_data.roicode});
            [tmp,ind_all,ind_subj]=intersect(parms.hippo_roicodes,roicodes);
            tmp_data = cell2mat({roi_data(ind_subj).avg});
            results.Hippo.data_change(r,ind_all) = tmp_data;
            % use nvoxels as baseline volume          
            tmp_data = cell2mat({roi_data(ind_subj).nvals});
            results.Hippo.data_baseline(r,ind_all) = tmp_data;
          end;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = concat_results(parms,results)
  results.All = [];
  results.All.nrois = 0;
  results.All.data_baseline = [];
  results.All.data_change = [];
  results.All.roicodes = [];
  results.All.roinames = [];
  if parms.aseg_flag
    results.All.roicodes = ...
      cat(1,results.All.roicodes,results.SubCort.roicodes);
    results.All.roinames = ...
      cat(1,results.All.roinames,results.SubCort.roinames);
    results.All.data_baseline = ...
      cat(2,results.All.data_baseline,results.SubCort.data_baseline);
    results.All.data_change = ...
      cat(2,results.All.data_change,results.SubCort.data_change);
  end;
  if parms.subhippo_flag
    results.All.roicodes = ...
      cat(1,results.All.roicodes,results.Hippo.roicodes);
    results.All.roinames = ...
      cat(1,results.All.roinames,results.Hippo.roinames);
    results.All.data_baseline = ...
      cat(2,results.All.data_baseline,results.Hippo.data_baseline);
    results.All.data_change = ...
      cat(2,results.All.data_change,results.Hippo.data_change);
  end;
  if parms.aparc_flag
    results.All.roicodes = ...
      cat(1,results.All.roicodes,results.Cort.roicodes);
    results.All.roinames = ...
      cat(1,results.All.roinames,results.Cort.roinames);
    results.All.data_baseline = ...
      cat(2,results.All.data_baseline,results.Cort.data_baseline);
    results.All.data_change = ...
      cat(2,results.All.data_change,results.Cort.data_change);
  end;
  results.All.nrois = length(results.All.roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

