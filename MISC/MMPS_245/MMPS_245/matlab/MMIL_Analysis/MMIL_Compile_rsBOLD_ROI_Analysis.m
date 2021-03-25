function results = MMIL_Compile_rsBOLD_ROI_Analysis(ProjID,varargin)
%function results = MMIL_Compile_rsBOLD_ROI_Analysis(ProjID,[options])
%
% Usage:
%  results = MMIL_Compile_rsBOLD_ROI_Analysis(ProjID,'key1', value1,...);
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
%      proc_bold, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%
% Optional Parameters that specify which analysis results to compile:
%   'indir': rsBOLD analysis subdirectory inside each BOLDPROC dir
%     {default = 'rsBOLD_analysis'}
%   'instem': rsBOLD analysis file stem
%     {default = 'rsBOLD'}
%   'roi_infix': additional string added to ROI output files
%     {default = []}
%   'corr_infix': string attached to rsBOLD correlation output files
%     {default = []}
%
% Other Optional Parameters:
%   'qcflag': [0|1] only include Visits from StudyInfo with QC = 1
%     {default = 1}
%   'continfo_flag': [0|1] include container-related information in StudyInfo
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%   'subjinfo_flag': [0|1] include subject information in StudyInfo
%     e.g. Sex, Site, DOB, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%
% Output:
%   results: struct containing compiled rsBOLD analysis results with fields:
%     StudyInfo: struct array containing info for each subject
%     RootDirs: struct containing paths of root directories
%     SubjIDs: cell array of subject IDs
%     VisitIDs: cell array of visit IDs
%     STRUCT_VisitIDs: cell array of structural (i.e. FreeSurfer) visit IDs
%     VisitNumbers: vector of visit numbers
%
% Created:  10/29/12 by Don Hagler
% Prev Mod: 08/29/13 by Don Hagler
% Last Mod: 04/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);
if parms.nsubs==0, return; end;

results = init_results(parms);

results = compile_roi_data(parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
... % specify which results to compile
    'indir','rsBOLD_analysis',[],...
    'instem','rsBOLD',[],...
    'roi_infix',[],[],...
    'corr_infix',[],[],...
    'qcflag',true,[false true],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'verbose',1,[0:2],...
... % hidden
    'required_containers',{'proc_bold','fsurf'},[],...
    'QC_raw',true,[false true],... % only applies if manual raw QC exists
    'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
    'QC_recon',true,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(parms.StudyInfo), error('empty StudyInfo'); end;

  parms = check_roi_data(parms);
  if parms.nsubs==0
    fprintf('%s: ERROR: no valid subjects\n',mfilename);
    return;
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
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  
  % set roinames based on first subject's data
  subj_data = load_subj_data(parms,1);
  results.roinames = subj_data.roinames;
  results.seed_roinames = subj_data.seed_roinames;
  results.seed_roinames_full = subj_data.seed_roinames_full;

  % store relevant numbers
  results.nseeds = length(results.seed_roinames);
  results.nrois = length(results.roinames);
  results.nsubs = parms.nsubs;

  % initialize TR, nreps, numTRs, mean_motion, etc
  results.TR = nan(1,results.nsubs);
  results.nreps = nan(1,results.nsubs);
  results.numTRs = nan(1,results.nsubs);
  results.max_rz = nan(1,results.nsubs);
  results.max_rx = nan(1,results.nsubs);
  results.max_ry = nan(1,results.nsubs);
  results.max_dz = nan(1,results.nsubs);
  results.max_dx = nan(1,results.nsubs);
  results.max_dy = nan(1,results.nsubs);
  results.mean_trans = nan(1,results.nsubs);
  results.max_trans = nan(1,results.nsubs);
  results.mean_rot = nan(1,results.nsubs);
  results.max_rot = nan(1,results.nsubs);
  results.mean_motion = nan(1,results.nsubs);
  results.mode_motion = nan(1,results.nsubs);
  results.med_motion = nan(1,results.nsubs);
  results.min_motion = nan(1,results.nsubs);
  results.max_motion = nan(1,results.nsubs);
  results.mean_accel = nan(1,results.nsubs);
  
  % initialize roi_var, R, Z
  results.roi_var = nan(results.nrois,results.nsubs);
  results.R = nan(results.nrois,results.nseeds,results.nsubs);
  results.Z = nan(results.nrois,results.nseeds,results.nsubs);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnamelist = set_fnames(parms,s)
  fnamelist = [];
  instem = [parms.RootDirs.proc_bold '/' parms.StudyInfo(s).proc_bold '/' ...
            parms.indir '/' parms.instem];
  % info
  fnamelist{end+1} = [instem '_info.mat'];
  % motion
  fnamelist{end+1} = [instem '_motion.mat'];
  % variance
  fstem = instem;
  if ~isempty(parms.roi_infix)
    fstem = [fstem '_' parms.roi_infix];
  end;
  fnamelist{end+1} = [fstem '_roi_data_var.mat'];
  % correlations
  fstem = instem;
  if ~isempty(parms.roi_infix)
    fstem = [fstem '_' parms.roi_infix];
  end;
  if ~isempty(parms.corr_infix)
    fstem = [fstem '_' parms.corr_infix];
  end;
  fnamelist{end+1} = [fstem '_roi_corr.mat'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_roi_data(parms,results)
  parms.nsubs = length(parms.StudyInfo);
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    fnames = set_fnames(parms,s);
    for f=1:length(fnames)
      if ~exist(fnames{f},'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fnames{f});
        valid_flags(s) = 0;
        break;
      end;
    end;
  end;
  ind_valid = find(valid_flags);
  if length(ind_valid)<parms.nsubs
    parms.StudyInfo = parms.StudyInfo(ind_valid);
    parms.nsubs = length(parms.StudyInfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_roi_data(parms,results)
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).VisitID;
    subj_data = load_subj_data(parms,s);
    if length(subj_data.roinames) ~= results.nrois
      fprintf('%s: WARNING: wrong number of ROIs for %s (%d instead of %d)\n',...
        mfilename,VisitID,length(subj_data.roinames),results.nrois);
      continue;
    end;
    if ~all(strcmp(subj_data.roinames,results.roinames))
      fprintf('%s: WARNING: mismatched ROI names for %s\n',...
        mfilename,VisitID);
      continue;
    end;
    if length(subj_data.seed_roinames) ~= results.nseeds
      fprintf('%s: WARNING: wrong number of seed ROIs for %s (%d instead of %d)\n',...
        mfilename,VisitID,length(subj_data.seed_roinames),results.nseeds);
      continue;
    end;
    if ~all(strcmp(subj_data.seed_roinames,results.seed_roinames))
      fprintf('%s: WARNING: mismatched seed ROI names for %s\n',...
        mfilename,VisitID);
      continue;
    end;
    results.TR(s) = mean(subj_data.TRs);
    results.nreps(s) = sum(subj_data.nreps);
    results.numTRs(s) = sum(subj_data.numTRs);
    results.max_rz(s) = subj_data.max_rz;
    results.max_rx(s) = subj_data.max_rx;
    results.max_ry(s) = subj_data.max_ry;
    results.max_dz(s) = subj_data.max_dz;
    results.max_dx(s) = subj_data.max_dx;
    results.max_dy(s) = subj_data.max_dy;
    results.mean_trans(s) = subj_data.mean_trans;
    results.max_trans(s) = subj_data.max_trans;
    results.mean_rot(s) = subj_data.mean_rot;
    results.max_rot(s) = subj_data.max_rot;
    results.mean_motion(s) = subj_data.mean_motion;
    results.mode_motion(s) = subj_data.mode_motion;
    results.med_motion(s) = subj_data.med_motion;
    results.min_motion(s) = subj_data.min_motion;
    results.max_motion(s) = subj_data.max_motion;
    results.mean_accel(s) = subj_data.mean_accel;
    results.roi_var(:,s) = subj_data.roi_var;
    results.R(:,:,s) = subj_data.R;
    results.Z(:,:,s) = subj_data.Z;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subj_data = load_subj_data(parms,s)
  subj_data = [];
  fnamelist = set_fnames(parms,s);
  for f=1:length(fnamelist)
    load(fnamelist{f});
  end;
  subj_data.TRs = info.TRs;
  subj_data.nreps = info.nreps;
  subj_data.numTRs = info.numTRs;
  subj_data.max_rz = max_rz;
  subj_data.max_rx = max_rx;
  subj_data.max_ry = max_ry;
  subj_data.max_dz = max_dz;
  subj_data.max_dx = max_dx;
  subj_data.max_dy = max_dy;
  subj_data.mean_trans = mean_trans;
  subj_data.max_trans = max_trans;
  subj_data.mean_rot = mean_rot;
  subj_data.max_rot = max_rot;
  subj_data.mean_motion = mean_motion;
  subj_data.mode_motion = mode_motion;
  subj_data.med_motion = med_motion;
  subj_data.min_motion = min_motion;
  subj_data.max_motion = max_motion;
  subj_data.mean_accel = mean_accel;
  subj_data.roi_var = roi_var;
  subj_data.roinames = roinames;
  subj_data.seed_roinames = seed_roinames;
  subj_data.seed_roinames_full = seed_roinames_full;
  subj_data.R = R;
  subj_data.Z = Z;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

