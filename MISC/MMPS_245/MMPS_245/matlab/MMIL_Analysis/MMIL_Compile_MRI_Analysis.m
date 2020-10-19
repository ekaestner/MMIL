function results = MMIL_Compile_MRI_Analysis(ProjID,varargin)
% function results = MMIL_Compile_MRI_Analysis(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%      fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'QC_raw' : [0|1] use good raw QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RawQC.mat
%     {default = 1}
%  'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%
% Optional Parameters Controlling which Stats to Compile:
%  'thick_flag': cortical thickness from aparc stats or fuzzy cluster ROIs
%    {default = 1}
%  'sulc_flag': sulcal depth from aparc stats or fuzzy cluster ROIs
%    {default = 1}
%  'area_flag': cortical area from aparc stats or fuzzy cluster ROIs
%    {default = 1}
%  'cortvol_flag': cortical volume from aparc stats or fuzzy cluster ROIs
%    {default = 1}
%  'cortT1w_flag': average T1-weighted values for aparc ROIs
%      or fuzzy cluster ROIs
%    and calculate gray-white contrast
%    {default = 1}
%  'cortproc_flag': average values for aparc or fuzzy cluster ROIs
%    also calculate gray-white contrast
%    with input determined by 'proc_outputlist'
%    {default = 1}
%  'subvol_flag': subcortical volume from aseg stats
%    {default = 1}
%  'subT1w_flag': average T1-weighted values for aseg ROIs 
%    {default = 1}
%  'subproc_flag': average values for aseg ROIs 
%    with input determined by 'proc_outputlist'
%    {default = 1}
%  'aparc_flag': [0|1] use aparc analysis results
%    for thickness, sulc, area, cortvol, and T1w
%    {default = 1}
%  'aparc_infix': string in aparc annot and stats file
%     e.g. 'aparc', 'aparc.a2009s'
%    {default = 'aparc'}
%  'fuzzy_flag': [0|1] use fuzzy cluster weighted surface ROIs
%      for thickness, sulc, area, cortvol, and T1w
%    applies if one or more of thick_flag, sulc_flag, area_flag,
%      cortvol_flag, cortT1w_flag, or cortproc_flag = 1
%    {default = 1}
%  'fuzzy_dir': input directory for fuzzy cluster ROIs
%     if empty, will use [getenv('MMPS_dir') '/atlases/fuzzy_clusters']
%     {default = []}
%  'fuzzy_fstem': input file stem for fuzzy cluster ROIs
%     with expected names like {fstem}{order}-{hemi}.mgz
%     {default = 'fuzzy'}
%  'fuzzy_order': [0|2|4|12|18] number of fuzzy cluster ROIs
%     note: set of 18 includes combined sets of 2, 4, and 12
%     if order=0, names are like {fstem}-{hemi}.mgz
%     {default = 18}
%  'concat_flag': [0|1] concatenate results for all analysis types
%    {default = 0}
%
% Other Optional Parameters:
%  'proc_outputlist': list of output file stems; e.g. {'T1w','T2w'}
%       if proc_outputlist contains 'T1w'
%         cortT1w_flag and subT1w_flag are irrelevant
%     {default = []}
%  'continfo_flag': [0|1] include container-related information in StudyInfo
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%  'subjinfo_flag': [0|1] include subject information in StudyInfo
%     e.g. Sex, Site, DOB, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%  'analysis_outdir': where analyzed mat files can be found
%     relative to FreeSurfer Container
%     {default = 'analysis'}
%  'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-0.2,0.2]}
%   'check_complete_flag': [0|1] whether to require that recon is complete
%     {default = 1}
%   'FS_version': which version of Freesurfer used (e.g. 305, 450, 510, 530)
%     for checking whether recon is complete
%     if empty, will use FREESURER_VER environment variable
%       or get from ContainerInfo
%     {default = []}
%   'baseflag': [0|1] only include baseline visits (VisitNumber = 1)
%     ignored if VisitNumber not specified in VisitInfo
%     {default = 1}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%
% Output:
%   results: struct containing these fields:
%     nsubs: number of subjects
%     SubjIDs: cell array of subject IDs
%     VisitIDs: cell array of visit IDs
%     STRUCT_VisitIDs: structural (i.e. FreeSurfer) visit IDs
%     StudyInfo: struct array containing info for each subject
%     cort_thick: struct containing cortical thickness results
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     cort_sulc: struct containing sulcal depth results
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     cort_area: struct containing cortical area results
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     cort_T1w: struct containing cortical T1w results
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI x 2)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     cort_contrast: struct containing cortical gray-white contrast
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     subcort_vol: struct containing subcortical ROI volumes
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     subcort_T1w: struct containing subcortical T1w results
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%     all: all results combined (if concat_flag = 1)
%       nrois: number of ROIs
%       data: matrix of ROI averages (subject X ROI)
%       roicodes: ROI code numbers
%       roinames: ROI names
%
%   NOTE: NaN's will be inserted for subjects with incomplete FreeSurfer recons
%         and in the case of missing data
%
% Created:  05/22/09 by Don Hagler
% Prev Mod: 01/12/17 by Don Hagler
% Last Mod: 07/18/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

results = init_results(parms);

if parms.subvol_flag
  results = compile_aseg_stats(parms,results);
end;

if parms.thick_flag || parms.sulc_flag || parms.area_flag || parms.cortvol_flag
  % compile thickness, area, and cortvol from aparc stats
  if parms.aparc_flag
    results = compile_aparc_stats(parms,results);
  end;
  % compile results for sulc values
  if parms.sulc_flag && parms.aparc_flag
    results = compile_sulc_aparc_data(parms,results);
  end;
  % compile results for fuzzy clusters
  if parms.fuzzy_flag
    if parms.thick_flag
      results = compile_fuzzy_data(parms,results,'thick');
    end;
    if parms.sulc_flag
      results = compile_fuzzy_data(parms,results,'sulc');
    end;
    if parms.area_flag
      results = compile_fuzzy_data(parms,results,'area');
    end;
    if parms.cortvol_flag
      results = compile_fuzzy_data(parms,results,'vol');
    end;
  end;
  % calculate mean cort thickness weighted by area for aparc_roilist
  if parms.thick_flag && parms.aparc_flag
    results = calc_mean_cort_vals(parms,results,'thick');
  end;
  % calculate mean sulcal depth weighted by area for aparc_roilist
  if parms.sulc_flag && parms.aparc_flag
    results = calc_mean_cort_vals(parms,results,'sulc');
  end;
  % calculate total cort area for aparc_roilist
  if parms.area_flag && parms.aparc_flag
    results = calc_mean_cort_vals(parms,results,'area');
  end;  
  % calculate total cort volume for aparc_roilist
  if parms.cortvol_flag && parms.aparc_flag
    results = calc_mean_cort_vals(parms,results,'vol');
  end;  
end;

% compile results for image values in aseg ROIs
if parms.subproc_flag
  for n=1:parms.ninputs
    meas = parms.proc_outputlist{n};
    results = compile_proc_aseg_data(parms,results,meas);
  end;
end;

% compile results for cort white and gray image values
if parms.cortproc_flag
  for n=1:parms.ninputs
    meas = parms.proc_outputlist{n};
    results = compile_proc_aparc_data(parms,results,meas);
    % compile proc results for fuzzy clusters
    if parms.fuzzy_flag
      results = compile_proc_fuzzy_data(parms,results,meas);  
    end;
    if parms.contrast_flag
      % calculate gray-white contrast for each ROI
      results = calc_proc_contrast(parms,results,meas);
      if parms.aparc_flag
        % calculate mean cort contrast for aparc_roilist
        results = calc_mean_cort_vals(parms,results,...
                                      sprintf('contrast_%s',meas));
      end;
    end;
    if parms.aparc_flag
      % calculate mean cort proc white and gray for aparc_roilist
      results = calc_mean_cort_vals(parms,results,meas);
    end;
  end;
end;

if parms.concat_flag
  results = concat_results(parms,results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'QC_raw',true,[false true],...
    'QC_recon',true,[false true],...
    'qcflag',true,[false true],...
  ...
    'thick_flag',true,[false true],...
    'sulc_flag',true,[false true],...
    'area_flag',true,[false true],...
    'cortvol_flag',true,[false true],...
    'cortT1w_flag',true,[false true],...
    'cortproc_flag',true,[false true],...
    'subvol_flag',true,[false true],...
    'subT1w_flag',true,[false true],...
    'subproc_flag',true,[false true],...
    'aparc_flag',true,[false true],...
    'aparc_infix','aparc',[],...
    'fuzzy_flag',true,[false true],...
    'fuzzy_dir',[],[],...
    'fuzzy_fstem','fuzzy',[],...
    'fuzzy_order',18,[0,2,4,12,18],...
    'concat_flag',false,[false true],...
  ...
    'proc_outputlist',[],[],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'analysis_outdir','analysis',[],...
    'projdist_list',[-0.2,0.2],[-5,5],...
    'check_complete_flag',true,[false true],...
    'FS_version',[],[],...
    'baseflag',true,[false true],...
    'verbose',1,[0:2],...
  ...
    'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,251:255],[1,Inf],...
    'extra_roilist',[10001:10003,20001:20003,20009,20010],[1,Inf],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'fname_fscolorlut',[],[],...
  ...
    'required_containers',{'fsurf'},[],...
    'modality','MRI',[],...
    'erode_nvoxels',1,[1:100],...
    'erode_flag',false,[false true],...  
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'fuzzy_name_tags',{'fuzzy_fstem','fuzzy_order'},[],...
  };

  parms = mmil_args2parms(options,parms_filter);
  
  parms.nhemi = length(parms.hemilist);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(parms.StudyInfo), error('empty StudyInfo'); end;

  % exclude non-baseline subjects
  if parms.baseflag && isfield(parms.StudyInfo,'VisitNumber')
    visitnums = cell2mat({parms.StudyInfo.VisitNumber});
    ind = find(visitnums==1);
    parms.StudyInfo = parms.StudyInfo(ind);
    if isempty(parms.StudyInfo)
      error('no baseline subjects -- check StudyInfo VisitNumber');
    end;
  end;

  % include only unique STRUCT_VisitIDs
  STRUCT_VisitIDs = unique({parms.StudyInfo.STRUCT_VisitID});
  VisitIDs = {parms.StudyInfo.VisitID};
  [tmp,ind_uniq] = intersect(VisitIDs,STRUCT_VisitIDs);
  parms.StudyInfo = parms.StudyInfo(ind_uniq);
  parms.nsubs = length(parms.StudyInfo);

  % determine whether surfaces are needed for each subject
  if (parms.aparc_flag || parms.fuzzy_flag) && ...
     (parms.thick_flag || parms.sulc_flag || parms.area_flag ||...
     parms.cortvol_flag || parms.cortT1w_flag || parms.cortproc_flag)
    parms.surf_flag = 1;
  else
    parms.surf_flag = 0;
  end;

  % exclude subjects with incomplete FreeSurfer recon
  if parms.check_complete_flag
    keep_flags = boolean(zeros(parms.nsubs,1));
    for s=1:parms.nsubs
      subj = parms.StudyInfo(s).fsurf;
      if isempty(subj), continue; end;
      FSContainerPath = [parms.RootDirs.fsurf '/' subj];
      [status,message] = MMIL_Get_FSReconStatus(FSContainerPath,parms.FS_version);
      if (parms.surf_flag && ~ismember(status,[2,5])) || ...
         (~parms.subvol_flag && ~ismember(status,[2,5,6]))
        if parms.verbose
          fprintf('%s: WARNING: incomplete recon for %s\n',mfilename,subj);
        end;
        continue;
      end;
      keep_flags(s) = 1;
    end;
    if ~all(keep_flags==1)
      parms.StudyInfo = parms.StudyInfo(keep_flags);
      parms.nsubs = length(parms.StudyInfo);
      if parms.nsubs==0, error('no valid subjects'); end;
    end;
  end;

  % get container-related info
  if parms.continfo_flag
    parms.StudyInfo = MMIL_Get_ContInfo(parms.StudyInfo,parms.RootDirs);
  end;

  % get subject-related info
  if parms.subjinfo_flag
    parms.StudyInfo = ...
      MMIL_Get_SubjInfo(parms.StudyInfo,parms.RootDirs,parms.ProjID);
  end;

  % set string included in eroded aseg roi_data files
  if parms.erode_flag
    parms.erode_outfix = '_erode';
    if parms.erode_nvoxels>1
      parms.erode_outfix = sprintf('%s%d',...
        parms.erode_outfix,parms.erode_nvoxels);
    end;
  else
    parms.erode_outfix = [];
  end;

  % read roi codes and names from FreeSurfer color lookup table
  if isempty(parms.fname_fscolorlut)
    parms.fname_fscolorlut = [getenv('MMPS_PARMS') '/MMIL_FSColorLUT.txt'];
  end;
  if ~exist(parms.fname_fscolorlut,'file')
    error('FreeSurfer color lookup table file %s not found',...
      parms.fname_fscolorlut);
  end;
  [parms.all_roicodes,parms.all_roinames] = fs_colorlut(parms.fname_fscolorlut);

  % get aseg roi codes and names
  [tmp,ind_aseg]=intersect(parms.all_roicodes,parms.aseg_roilist);
  parms.aseg_roicodes = parms.all_roicodes(ind_aseg);
  parms.aseg_roinames = parms.all_roinames(ind_aseg);
  parms.aseg_nrois = length(ind_aseg);

  % get subcort roi codes and names
  subcort_roilist = [parms.aseg_roilist,parms.extra_roilist];
  [tmp,ind_subcort]=intersect(parms.all_roicodes,subcort_roilist);
  parms.subcort_roicodes = parms.all_roicodes(ind_subcort);
  parms.subcort_roinames = parms.all_roinames(ind_subcort);
  parms.subcort_nrois = length(ind_subcort);

  if parms.surf_flag
    % initialize roicodes and names
    parms.surf_roicodes = [];
    parms.surf_roinames = [];
    parms.surf_nrois = 0;

    % create roi codes and names for aparc ROIs
    if parms.aparc_flag
      parms = set_aparc_parms(parms);
    end;

    % create roi codes and names for fuzzy cluster ROIs
    if parms.fuzzy_flag
      parms = set_fuzzy_parms(parms);
    end;

    % whether to calculate gray-white contrast depends on if 2 projdist
    parms.nprojdist = length(parms.projdist_list);
    if parms.nprojdist==2
      parms.contrast_flag = 1;
      parms.projdist_list = sort(parms.projdist_list); % white matter first
    else
      parms.contrast_flag = 0;
    end;
  end;

  % check proc_outputlist
  if isempty(parms.proc_outputlist)
    parms.cortproc_flag = 0;
    parms.subproc_flag = 0;
    parms.ninputs = 0;
  elseif ~iscell(parms.proc_outputlist)
    parms.proc_outputlist = {parms.proc_outputlist};
    parms.ninputs = length(parms.proc_outputlist);
  end;

  % check T1w flag
  if parms.cortT1w_flag || parms.subT1w_flag
    if ~parms.cortproc_flag && ~parms.subproc_flag
      parms.proc_outputlist = {'T1w'};
    else
      parms.proc_outputlist = union(parms.proc_outputlist,{'T1w'});
    end;
    if parms.cortT1w_flag, parms.cortproc_flag = 1; end;
    if parms.subT1w_flag, parms.subproc_flag = 1; end;
    parms.ninputs = length(parms.proc_outputlist);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_aparc_parms(parms)
  % get aparc roi codes and names
  [tmp,ind_aparc]=intersect(parms.all_roicodes,parms.aparc_roilist);
  parms.aparc_roicodes = parms.all_roicodes(ind_aparc);
  parms.aparc_roinames = parms.all_roinames(ind_aparc);
  parms.aparc_nrois = length(ind_aparc);
  % add to surf roi codes and names
  parms.surf_roicodes = cat(1,parms.surf_roicodes,parms.aparc_roicodes);
  parms.surf_roinames = cat(1,parms.surf_roinames,parms.aparc_roinames);
  parms.surf_nrois = parms.surf_nrois + parms.aparc_nrois;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_fuzzy_parms(parms)
  if isempty(parms.fuzzy_dir)
    parms.fuzzy_dir = [getenv('MMPS_DIR') '/atlases/fuzzy_clusters'];
  end;
  if ~exist(parms.fuzzy_dir,'dir')
    error('fuzzy cluster dir %s not found',parms.fuzzy_dir);
  end;
  if parms.fuzzy_order~=0
    parms.fuzzy_full_fstem = ...
      sprintf('%s%d',parms.fuzzy_fstem,parms.fuzzy_order);
    parms.fuzzy_nrois = 2*parms.fuzzy_order;
    args = mmil_parms2args(parms,parms.fuzzy_name_tags);
    fuzzy_names = mmil_fuzzy_names(args{:});
    parms.fuzzy_roicodes = [];
    parms.fuzzy_roinames = [];
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      tmp_codes = 20000 + (h-1)*5000 + ...
                  parms.fuzzy_order*100 + [1:parms.fuzzy_order]';
      tmp_names = cell(parms.fuzzy_order,1);
      for f=1:parms.fuzzy_order
        tmp_names{f} = sprintf('ctx-%s-%s',hemi,fuzzy_names{f});
      end;
      parms.fuzzy_roicodes = cat(1,parms.fuzzy_roicodes,tmp_codes);
      parms.fuzzy_roinames = cat(1,parms.fuzzy_roinames,tmp_names);
    end;
  else
    parms.fuzzy_full_fstem = parms.fuzzy_fstem;
    parms.fuzzy_nrois = 0;
    parms.fuzzy_roicodes = [];
    parms.fuzzy_roinames = [];
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fname_weights = sprintf('%s/%s-%s.mgz',...
        parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
      if ~exist(fname_weights,'file')
        error('fuzzy cluster file %s not found',fname_weights);
      end;
      [M,volsz] = fs_read_header(fname_weights);
      nrois = volsz(4);
      parms.fuzzy_nrois = parms.fuzzy_nrois + nrois;
      tmp_codes = 20000 + (h-1)*5000 + ...
                  nrois*100 + [1:nrois]';
      % load roinames
      fname_txt = sprintf('%s/%s_%s_roinames.txt',...
        parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
      if ~exist(fname_txt,'file')
        fprintf('%s: WARNING: roinames file %s not found',...
          fname_txt);
        tmp_names = cell(nrois,1);
        for f=1:nrois
          tmp_names{f} = sprintf('ctx-%s_%d-%s',hemi,parms.fuzzy_fstem,f);
        end;
      else
        tmp_names = mmil_readtext(fname_txt);
      end;
      parms.fuzzy_roicodes = cat(1,parms.fuzzy_roicodes,tmp_codes);
      parms.fuzzy_roinames = cat(1,parms.fuzzy_roinames,tmp_names);
    end;
  end;
  parms.surf_roicodes = cat(1,parms.surf_roicodes,parms.fuzzy_roicodes);
  parms.surf_roinames = cat(1,parms.surf_roinames,parms.fuzzy_roinames);
  parms.surf_nrois = parms.surf_nrois + parms.fuzzy_nrois;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results = [];
  results.StudyInfo = parms.StudyInfo;
  results.RootDirs = parms.RootDirs;
  results.SubjIDs = {parms.StudyInfo.SubjID};
  results.VisitIDs = {parms.StudyInfo.VisitID};
  results.STRUCT_VisitIDs = {parms.StudyInfo.STRUCT_VisitID};
  if ~isfield(parms.StudyInfo,'VisitNumber')
    results.VisitNumbers = zeros(parms.nsubs,1);
  else
    results.VisitNumbers = cell2mat({parms.StudyInfo.VisitNumber});
  end;
  results.nsubs = parms.nsubs;
  % cortical thickness
  if parms.thick_flag
    results.cort_thick.nrois = parms.surf_nrois;
    results.cort_thick.data = nan(results.nsubs,parms.surf_nrois);
    results.cort_thick.roicodes = parms.surf_roicodes;
    results.cort_thick.roinames = parms.surf_roinames;
  end;
  % cortical sulcal depth
  if parms.sulc_flag
    results.cort_sulc.nrois = parms.surf_nrois;
    results.cort_sulc.data = nan(results.nsubs,parms.surf_nrois);
    results.cort_sulc.roicodes = parms.surf_roicodes;
    results.cort_sulc.roinames = parms.surf_roinames;
  end;
  % cortical area
  if parms.area_flag
    results.cort_area.nrois = parms.surf_nrois;
    results.cort_area.data = nan(results.nsubs,parms.surf_nrois);
    results.cort_area.roicodes = parms.surf_roicodes;
    results.cort_area.roinames = parms.surf_roinames;
  end;
  % cortical volume
  if parms.cortvol_flag
    results.cort_vol.nrois = parms.surf_nrois;
    results.cort_vol.data = nan(results.nsubs,parms.surf_nrois);
    results.cort_vol.roicodes = parms.surf_roicodes;
    results.cort_vol.roinames = parms.surf_roinames;
  end;
  % cortical image values
  if parms.cortproc_flag
    for n=1:parms.ninputs
      fname = sprintf('cort_%s',parms.proc_outputlist{n});
      results.(fname).nrois = parms.surf_nrois;
      results.(fname).data = nan(results.nsubs,...
        parms.surf_nrois,length(parms.projdist_list));
      results.(fname).roicodes = parms.surf_roicodes;
      results.(fname).roinames = parms.surf_roinames;
      % contrast ratio
      if parms.contrast_flag
        fname = sprintf('cort_contrast_%s',parms.proc_outputlist{n});
        results.(fname).nrois = parms.surf_nrois;
        results.(fname).data = nan(results.nsubs,parms.surf_nrois);
        results.(fname).roicodes = parms.surf_roicodes;
        results.(fname).roinames = parms.surf_roinames;
      end;
    end;
  end;
  % subcortical (aseg) volume
  if parms.subvol_flag
    results.subcort_vol.nrois = parms.subcort_nrois;
    results.subcort_vol.data = nan(results.nsubs,parms.subcort_nrois);
    results.subcort_vol.roicodes = parms.subcort_roicodes;
    results.subcort_vol.roinames = parms.subcort_roinames;
  end;
  % subcortical (aseg) image values
  if parms.subproc_flag
    for n=1:parms.ninputs
      fname = sprintf('subcort_%s',parms.proc_outputlist{n});
      results.(fname).nrois = parms.aseg_nrois;
      results.(fname).data = nan(results.nsubs,parms.aseg_nrois);
      results.(fname).roicodes = parms.aseg_roicodes;
      results.(fname).roinames = parms.aseg_roinames;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_aseg_stats(parms,results)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load aseg_stats
    fname_mat = sprintf('%s/%s/aseg_stats.mat',...
      FSContainerPath,parms.analysis_outdir);
    aseg_stats = [];
    if ~exist(fname_mat,'file')
      if parms.verbose
        fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
      end;
    else
      if parms.verbose==2
        fprintf('%s: loading aseg stats for %s...\n',...
          mfilename,parms.StudyInfo(s).VisitID);
      end;
      load(fname_mat);
      if isempty(aseg_stats)
        if parms.verbose
          fprintf('%s: WARNING: empty aseg_stats in %s\n',...
            mfilename,fname_mat);
        end;
      else
        roicodes = cell2mat({aseg_stats.roicode});
        [tmp,ind_all,ind_subj]=intersect(parms.subcort_roicodes,roicodes);
        tmp_data = cell2mat({aseg_stats(ind_subj).volume});
        results.subcort_vol.data(s,ind_all) = tmp_data;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_aparc_stats(parms,results)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load aparc_stats
    fname_mat = sprintf('%s/%s/%s_stats.mat',...
      FSContainerPath,parms.analysis_outdir,parms.aparc_infix);
    aparc_stats = [];
    if ~exist(fname_mat,'file')
      if parms.verbose
        fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
      end;
    else
      aparc_stats = [];
      if parms.verbose==2
        fprintf('%s: loading aparc stats for %s...\n',...
          mfilename,parms.StudyInfo(s).VisitID);
      end;
      load(fname_mat);
      if isempty(aparc_stats)
        if parms.verbose
          fprintf('%s: WARNING: empty aparc_stats in %s\n',...
            mfilename,fname_mat);
        end;
      else
        roicodes = cell2mat({aparc_stats.roicode});
        [tmp,ind_all,ind_subj]=intersect(parms.aparc_roicodes,roicodes);
        % thickness
        if parms.thick_flag
          tmp_data = cell2mat({aparc_stats(ind_subj).thickavg});
          results.cort_thick.data(s,ind_all) = tmp_data;
        end;
        % area
        if parms.area_flag
          tmp_data = cell2mat({aparc_stats(ind_subj).surfarea});
          results.cort_area.data(s,ind_all) = tmp_data;
        end;
        % volume
        if parms.cortvol_flag
          tmp_data = cell2mat({aparc_stats(ind_subj).grayvol});
          results.cort_vol.data(s,ind_all) = tmp_data;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_mean_cort_vals(parms,results,meas)
  % results struct field name for meas
  tag = ['cort_' meas];
  % get roicodes from results
  roicodes = results.(tag).roicodes;
  roinames = results.(tag).roinames;
  codes_all = parms.aparc_roilist;
  codes_lh = codes_all(codes_all<2000);
  codes_rh = codes_all(codes_all>=2000);
  % identify ROIs from aparc (excluding fuzzy clusters)
  [tmp,ind_lh]=intersect(roicodes,codes_lh);
  [tmp,ind_rh]=intersect(roicodes,codes_rh);
  [tmp,ind_all]=intersect(roicodes,codes_all);
  % get data
  data_lh = results.(tag).data(:,ind_lh,:);
  data_rh = results.(tag).data(:,ind_rh,:);
  data_all = results.(tag).data(:,ind_all,:);
  % calculate mean or sum of ROI measures
  switch meas
    case {'area','vol'}
      mean_lh = sum(data_lh,2);
      mean_rh = sum(data_rh,2);
      mean_all = sum(data_all,2);
      new_roinames = {'ctx-lh-total';'ctx-rh-total';'ctx-total'};
    otherwise
      if parms.area_flag
        np = size(data_all,3);
        % get area data
        area_lh = repmat(results.cort_area.data(:,ind_lh),[1,1,np]);
        area_rh = repmat(results.cort_area.data(:,ind_rh),[1,1,np]);
        area_all = repmat(results.cort_area.data(:,ind_all),[1,1,np]);
      else
        area_lh = ones(size(data_lh));
        area_rh = ones(size(data_rh));
        area_all = ones(size(data_all));
      end;
      % calculate mean cort val weighted by area for aparc_roilist
      mean_lh = mmil_wtd_mean(data_lh,area_lh,2);
      mean_rh = mmil_wtd_mean(data_rh,area_rh,2);
      mean_all = mmil_wtd_mean(data_all,area_all,2);
      new_roinames = {'ctx-lh-mean';'ctx-rh-mean';'ctx-mean'};
  end;
  % set roicodes
  new_roicodes = [1040;2040;2050];
  % concatenate to data matrix, roicodes, and roinames
  results.(tag).data = cat(2,results.(tag).data,mean_lh,mean_rh,mean_all);
  results.(tag).roicodes = cat(1,results.(tag).roicodes,new_roicodes);
  results.(tag).roinames = cat(1,results.(tag).roinames,new_roinames);
  results.(tag).nrois = length(results.(tag).roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_proc_aseg_data(parms,results,meas)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load average aseg ROI data
    fname_mat = sprintf('%s/%s/%s_aseg%s_roi_data.mat',...
      FSContainerPath,parms.analysis_outdir,meas,parms.erode_outfix);
    roi_data = [];
    if ~exist(fname_mat,'file')
      if parms.verbose
        fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
      end;
    else
      roi_data = [];
      if parms.verbose==2
        fprintf('%s: loading %s aseg ROI data for %s...\n',...
          mfilename,meas,parms.StudyInfo(s).VisitID);
      end;
      load(fname_mat);
      if isempty(roi_data)
        if parms.verbose
          fprintf('%s: WARNING: missing %ss aseg roi_data in %s\n',...
            mfilename,meas,fname_mat);
        end;
      else
        roicodes = cell2mat({roi_data.roicode});
        [tmp,ind_all,ind_subj]=intersect(parms.aseg_roicodes,roicodes);
        tmp_data = cell2mat({roi_data(ind_subj).avg});
        results.(['subcort_' meas]).data(s,ind_all) = tmp_data;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_sulc_aparc_data(parms,results)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load sulc data
    fname_mat = sprintf('%s/%s/sulc_%s_roi_data.mat',...
      FSContainerPath,parms.analysis_outdir,parms.aparc_infix);
    if ~exist(fname_mat,'file')
      if parms.verbose
        fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
      end;
    else
      roi_data = [];
      if parms.verbose==2
        fprintf('%s: loading sulc aparc ROI data for %s...\n',...
          mfilename,parms.StudyInfo(s).VisitID);
      end;
      load(fname_mat);
      if isempty(roi_data)
        if parms.verbose
          fprintf('%s: WARNING: missing roi_data in %s\n',...
            mfilename,fname_mat);
        end;
      else
        roicodes = cell2mat({roi_data.roicode});
        [tmp,ind_all,ind_subj]=intersect(parms.surf_roicodes,roicodes);
        tmp_data = cell2mat({roi_data(ind_subj).avg});
        results.cort_sulc.data(s,ind_all) = tmp_data;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_proc_aparc_data(parms,results,meas)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load contrast data
    for p=1:parms.nprojdist
      projdist = parms.projdist_list(p);
      fname_mat = sprintf('%s/%s/%s_%s_pdist%0.1f_roi_data.mat',...
        FSContainerPath,parms.analysis_outdir,meas,parms.aparc_infix,projdist);
      if ~exist(fname_mat,'file')
        if parms.verbose
          fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
        end;
      else
        roi_data = [];
        if parms.verbose==2
          fprintf('%s: loading %s aparc ROI data for %s...\n',...
            mfilename,meas,parms.StudyInfo(s).VisitID);
        end;
        load(fname_mat);
        if isempty(roi_data)
          if parms.verbose
            fprintf('%s: WARNING: missing roi_data in %s\n',...
              mfilename,fname_mat);
          end;
        else
          roicodes = cell2mat({roi_data.roicode});
          [tmp,ind_all,ind_subj]=intersect(parms.surf_roicodes,roicodes);
          tmp_data = cell2mat({roi_data(ind_subj).avg});
          results.(['cort_' meas]).data(s,ind_all,p) = tmp_data;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_fuzzy_data(parms,results,meas)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load fuzzy ROI data
    switch meas
      case 'thick'
        instem = 'thickness';
      case 'sulc'
        instem = 'sulc';
      case 'area'
        instem = 'area';
      case 'vol'
        instem = 'cortvol';
    end;
    fname_mat = sprintf('%s/%s/%s_%s_roi_data.mat',...
      FSContainerPath,parms.analysis_outdir,...
      instem,parms.fuzzy_full_fstem);
    if ~exist(fname_mat,'file')
      if parms.verbose
        fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
      end;
    else
      roi_data = [];
      if parms.verbose==2
        fprintf('%s: loading fuzzy ROI %s data for %s...\n',...
          mfilename,meas,parms.StudyInfo(s).VisitID);
      end;
      load(fname_mat);
      if isempty(roi_data)
        if parms.verbose
          fprintf('%s: WARNING: missing roi_data in %s\n',...
            mfilename,fname_mat);
        end;
      else
        % match fuzzy ROIs by name
        roinames = {roi_data.roiname};
        [tmp,ind_all,ind_subj]=intersect(parms.surf_roinames,roinames);
        tmp_data = cell2mat({roi_data(ind_subj).avg});
        results.(['cort_' meas]).data(s,ind_all) = tmp_data;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_proc_fuzzy_data(parms,results,meas)
  for s=1:results.nsubs
    FSContainerPath = [parms.RootDirs.fsurf '/' parms.StudyInfo(s).fsurf];
    % load contrast data
    for p=1:parms.nprojdist
      projdist = parms.projdist_list(p);
      fname_mat = sprintf('%s/%s/%s_%s_pdist%0.1f_roi_data.mat',...
        FSContainerPath,parms.analysis_outdir,meas,...
        parms.fuzzy_full_fstem,projdist);
      if ~exist(fname_mat,'file')
        if parms.verbose
          fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
        end;
      else
        roi_data = [];
        if parms.verbose==2
          fprintf('%s: loading fuzzy ROI %s data for %s...\n',...
            mfilename,meas,parms.StudyInfo(s).VisitID);
        end;
        load(fname_mat);
        if isempty(roi_data)
          if parms.verbose
            fprintf('%s: WARNING: missing roi_data in %s\n',...
              mfilename,fname_mat);
          end;
        else
          % match fuzzy ROIs by name
          roinames = {roi_data.roiname};
          [tmp,ind_all,ind_subj]=intersect(parms.surf_roinames,roinames);
          tmp_data = cell2mat({roi_data(ind_subj).avg});
          % add to results.cort_meas.data
          results.(['cort_' meas]).data(s,ind_all,p) = tmp_data;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_proc_contrast(parms,results,meas)
  itag = sprintf('cort_%s',meas);
  otag = sprintf('cort_contrast_%s',meas);
  % calculate contrast ratio
  results.(otag).data = ...
    -diff(results.(itag).data,[],3) ./ mean(results.(itag).data,3);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = concat_results(parms,results)
  % concatenate cort_thick, cort_sulc, cort_vol, cort_area,
  %    cort_T1w, cort_contrast, subcort_vol, subcort_T1w
  results.all = [];
  results.all.nrois = 0;
  results.all.data = [];
  results.all.roicodes = [];
  results.all.roinames = [];
  if parms.thick_flag
    results.all.roicodes = ...
      cat(1,results.all.roicodes,results.cort_thick.roicodes);
    tmp_roinames = results.cort_thick.roinames;
    stem = 'cort_thick';
    for i=1:length(tmp_roinames)
      tmp_roinames{i} = [stem '-' tmp_roinames{i}];
    end;
    results.all.roinames = ...
      cat(1,results.all.roinames,tmp_roinames);
    results.all.data = ...
      cat(2,results.all.data,results.cort_thick.data);
  end;
  if parms.sulc_flag
    results.all.roicodes = ...
      cat(1,results.all.roicodes,results.cort_sulc.roicodes);
    tmp_roinames = results.cort_sulc.roinames;
    stem = 'cort_sulc';
    for i=1:length(tmp_roinames)
      tmp_roinames{i} = [stem '-' tmp_roinames{i}];
    end;
    results.all.roinames = ...
      cat(1,results.all.roinames,tmp_roinames);
    results.all.data = ...
      cat(2,results.all.data,results.cort_sulc.data);
  end;
  if parms.area_flag
    results.all.roicodes = ...
      cat(1,results.all.roicodes,results.cort_area.roicodes);
    tmp_roinames = results.cort_area.roinames;
    stem = 'cort_area';
    for i=1:length(tmp_roinames)
      tmp_roinames{i} = [stem '-' tmp_roinames{i}];
    end;
    results.all.roinames = ...
      cat(1,results.all.roinames,tmp_roinames);
    results.all.data = ...
      cat(2,results.all.data,results.cort_area.data);
  end;
  if parms.cortvol_flag
    results.all.roicodes = ...
      cat(1,results.all.roicodes,results.cort_vol.roicodes);
    tmp_roinames = results.cort_vol.roinames;
    stem = 'cort_vol';
    for i=1:length(tmp_roinames)
      tmp_roinames{i} = [stem '-' tmp_roinames{i}];
    end;
    results.all.roinames = ...
      cat(1,results.all.roinames,tmp_roinames);
    results.all.data = ...
      cat(2,results.all.data,results.cort_vol.data);
  end;
  if parms.cortproc_flag
    for n=1:parms.ninputs
      meas = parms.proc_outputlist{n};
      tag = ['cort_' meas];
      for p=1:parms.nprojdist
        projdist = parms.projdist_list(p);
        if projdist<0
          stem = sprintf('cort_%s_white%0.1f',meas,projdist);
        else
          stem = sprintf('cort_%s_gray+%0.1f',meas,projdist);
        end;
        results.all.roicodes = ...
          cat(1,results.all.roicodes,results.(tag).roicodes);
        tmp_roinames = results.(tag).roinames;
        for i=1:length(tmp_roinames)
          tmp_roinames{i} = [stem '-' tmp_roinames{i}];
        end;
        results.all.roinames = ...
          cat(1,results.all.roinames,tmp_roinames);
        results.all.data = ...
          cat(2,results.all.data,results.(tag).data(:,:,p));
      end;
      % contrast ratio
      if length(parms.projdist_list)==2
        tag = ['cort_contrast_' meas];
        results.all.roicodes = ...
          cat(1,results.all.roicodes,results.(tag).roicodes);
        tmp_roinames = results.cort_T1w.roinames;
        for i=1:length(tmp_roinames)
          tmp_roinames{i} = [tag '-' tmp_roinames{i}];
        end;
        results.all.roinames = ...
          cat(1,results.all.roinames,tmp_roinames);
        results.all.data = ...
          cat(2,results.all.data,results.(tag).data);
      end;
    end;
  end
  if parms.subvol_flag
    results.all.roicodes = ...
      cat(1,results.all.roicodes,results.subcort_vol.roicodes);
    tmp_roinames = results.subcort_vol.roinames;
    stem = 'subcort_vol';
    for i=1:length(tmp_roinames)
      tmp_roinames{i} = [stem '-' tmp_roinames{i}];
    end;
    results.all.roinames = ...
      cat(1,results.all.roinames,tmp_roinames);
    results.all.data = ...
      cat(2,results.all.data,results.subcort_vol.data);
  end;
  if parms.subproc_flag
    for n=1:parms.ninputs
      meas = parms.proc_outputlist{n};
      tag = ['subcort_' meas];
      results.all.roicodes = ...
        cat(1,results.all.roicodes,results.(tag).roicodes);
      tmp_roinames = results.(tag).roinames;
      for i=1:length(tmp_roinames)
        tmp_roinames{i} = [tag '-' tmp_roinames{i}];
      end;
      results.all.roinames = ...
        cat(1,results.all.roinames,tmp_roinames);
      results.all.data = ...
        cat(2,results.all.data,results.(tag).data);
    end;
  end;
  results.all.nrois = length(results.all.roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

