function results = ABCD_Summarize_rsBOLD_Analysis(ProjID,varargin)
%function results = ABCD_Summarize_rsBOLD_Analysis(ProjID,[options])
%
% Usage:
%  results = ABCD_Summarize_rsBOLD_Analysis(ProjID,'key1', value1,...);
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
% Optional Parameters that specify which values to include in output
%   'var_flag': [0|1] time series variance for each ROI
%     {default = 0}
%   'roi_flag': [0|1] correlation for each seed with each ROI
%     {default = 0}
%   'network_flag': [0|1] average correlation within networks of ROIs
%     0: no network correlations
%     1: mean correlation between each network and each ROI
%     2: mean correlation within each network
%     3: mean correlation within and between each network
%     {default = 2}
%   'continfo_flag': [0|1] include container-related information in StudyInfo
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%   'subjinfo_flag': [0|1] include subject information in StudyInfo
%     e.g. Sex, Site, DOB, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%
% Other Optional Parameters:
%   'qcflag': [0|1] only include Visits from StudyInfo with QC = 1
%     {default = 1}
%   'max_motion': exclude subjects with mean relative motion greater than this
%     {default = 0.5}
%   'roi_type': type of surface ROI to be used ('fparc', 'gordon', 'talrc', or 'aparc')
%     {default = 'fparc'}
%   'roinames': cell array of ROIs to be included in output
%     applies if var_flag = 1, roi_flag = 1, or network_flag = 1
%     if empty, use all ROIs
%     {default = []}
%   'networks': struct containing network definitions
%     fields of networks struct are names of networks with values
%       that are cell arrays of ROI names
%     if empty, will generate automatically based on roi_type
%     {default = []}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'outdir': output directory for summary csv files
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ROI_Summaries'}
%   'outstem': output file stem
%     {default = 'rsBOLD'}
%   'outfix': optional output file suffix attached to csv file
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
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
% Created:  10/31/12 by Don Hagler
% Prev Mod: 08/11/17 by Don Hagler
% Last Mod: 10/09/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

% check input parameters, get StudyInfo, RootDirs
parms = check_input(ProjID,varargin);

% check if output file exists
if exist(parms.fname_csv,'file') && ~parms.forceflag
  return;
end;

% compile rsBOLD ROI analysis results for each subject
results = compile_results(parms);
if isempty(results), return; end;

% identify subjects with excessive motion and exclude from results
results = exclude_motion(parms,results);

% select specific ROIs, exclude others
if ~isempty(parms.roinames)
  results = select_rois(parms,results);
end;

% calculate mean correlation for each network
if parms.network_flag
  results = calc_network_corr(parms,results);
end;

% create summary csv file
summarize_analysis(parms,results);

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
... % specify with values to include in output
    'var_flag',false,[false true],...
    'roi_flag',false,[false true],...
    'network_flag',2,[0:3],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
... % other
    'qcflag',true,[false true],...
    'max_motion',0.5,[],...
    'networks',[],[],...
    'roinames',[],[],...
    'roi_type','fparc',{'fparc','gordon','talrc','aparc'},...
    'verbose',1,[0:2],...
    'outdir','ROI_Summaries',[],...
    'outstem','rsBOLD',[],...
    'outfix',[],[],...
    'forceflag',false,[false true],...
...
    'QC_raw',true,[false true],... % only applies if manual raw QC exists
    'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
    'QC_recon',true,[false true],...
    'csv_outfix',[],[],...
    'mat_outfix',[],[],...
    'results_tags',{'StudyInfo','SubjIDs','VisitIDs','STRUCT_VisitIDs',...
                   'VisitNumbers','TR','numTRs'},[],... % 'nreps'
    'motion_tags',{'mean_motion','max_motion'},[],... % 'mean_trans','max_trans','mean_rot','max_rot'
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'excl_tags',{'networks','max_motion','networks','roi_type','roinames',...
                 'var_flag','roi_flag','network_flag',...
                 'excl_tags','ProjID','outdir','outstem','outfix',...
                 'forceflag','csv_outfix','mat_outfix',...
                 'results_tags','motion_tags','info_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  parms.compile_tags = setdiff(fieldnames(parms),parms.excl_tags);

  if mmil_isrelative(parms.outdir)
    parms.outdir = [getenv('HOME') '/MetaData/' parms.ProjID '/' parms.outdir];
  end;

  if isempty(parms.networks)
    parms.networks = abcd_define_rsBOLD_networks(parms.roi_type);
  end;
  
  if isempty(parms.csv_outfix)
    parms.csv_outfix = set_csv_outfix(parms);
  end;  
  if isempty(parms.mat_outfix)
    parms.mat_outfix = set_mat_outfix(parms);
  end;  

  parms.fname_csv = sprintf('%s/%s_%s.csv',...
    parms.outdir,parms.outstem,parms.csv_outfix);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(parms)
  fname_out = sprintf('%s/%s_%s.mat',...
    parms.outdir,parms.outstem,parms.mat_outfix);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose==2
      fprintf('%s: compiling rsBOLD ROI analysis results...\n',mfilename);
    end;
    mmil_mkdir(parms.outdir);
    args = mmil_parms2args(parms,parms.compile_tags);
    results = MMIL_Compile_rsBOLD_ROI_Analysis(parms.ProjID,args{:});
    if isempty(results), return; end;
    if ~isempty(fname_out), save(fname_out,'results'); end;
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = exclude_motion(parms,results)
  % select subjects with motion less than threshold
  ind_valid = find(results.mean_motion < parms.max_motion);
  results.nsubs = length(ind_valid);
  if parms.verbose==2
    fprintf('%s: n = %d with motion < %0.2f\n',...
      mfilename,results.nsubs,parms.max_motion);
  end;
  % keep results only for valid subjects
  tags = cat(2,parms.results_tags,parms.motion_tags);
  for i=1:length(tags)
    tag = tags{i};
    results.(tag) = results.(tag)(ind_valid);
  end;
  results.roi_var = results.roi_var(:,ind_valid);
  results.R = results.R(:,:,ind_valid);
  results.Z = results.Z(:,:,ind_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = select_rois(parms,results)
  % select values for ROIs in parms.roinames
  [tmp,ind_parms,ind_results] = intersect(parms.roinames,results.roinames);
  if isempty(ind_results)
    error('input roinames are not valid');
  end;
  % maintain order in parms.roinames
  [tmp,ind_sort] = sort(ind_parms);
  ind_results = ind_results(ind_sort);
  results.nrois = length(ind_results);
  if parms.verbose==2
    fprintf('%s: %d ROIs selected\n',mfilename,results.nrois);
  end;
  % keep results only for selected ROIs
  results.roinames = results.roinames(ind_results);
  results.roi_var = results.roi_var(ind_results,:);
  results.R = results.R(ind_results,:,:);
  results.Z = results.Z(ind_results,:,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_network_corr(parms,results)
  if parms.verbose==2
    fprintf('%s: calculating average network correlations...\n',mfilename);
  end;
  results.networks = parms.networks;
  results.network_names = fieldnames(parms.networks);
  results.nn = length(results.network_names);
  switch parms.network_flag
    case 1
      results.nt = results.nrois;
    case 2
      results.nt = 1;
    case 3
      results.nt = results.nn;
  end;
  % calculate mean correlation for each network
  results.network_corr = zeros(results.nn,results.nt,results.nsubs);
  for n=1:results.nn
    network_name = results.network_names{n};
    network_seed_roinames = results.networks.(network_name);
    ns = length(network_seed_roinames);
    switch parms.network_flag
      case 1
        targ_names = results.roinames;
      case 2
        targ_names = {network_name};
      case 3
        targ_names = results.network_names;
    end;
    for t=1:results.nt
      targ_name = targ_names{t};
      if parms.network_flag>1
        % get list of ROIs to average, excluding seed
        targ_roinames = results.networks.(targ_name);
      else
        targ_roinames = {targ_name};
      end;
      % find all ROIs with this name (e.g. lh and rh)
      ind_rois = [];
      for r=1:length(targ_roinames)
        ind_rois = [ind_rois,...
          find(~cellfun(@isempty,regexp(results.roinames,targ_roinames{r})))];
      end;
      nrois = length(ind_rois);
      for r=1:nrois
        ind_roi = ind_rois(r);
        roiname = results.roinames{ind_roi};
        % exclude seed ROI from targets
        ind_tmp = find(cellfun(@isempty,regexp(roiname,network_seed_roinames)));
        seed_roinames = network_seed_roinames(ind_tmp);
        [tmp,ind_seeds] = intersect(results.seed_roinames,seed_roinames);
        if isempty(ind_seeds)
          if parms.verbose
            fprintf('%s: WARNING: ind_seeds is empty for %s and %s\n',...
              mfilename,network_name,roiname);
          end;
          continue;
        end;
        % sum Z-transformed correlations for this ROI
        tmp_Z = mean(results.Z(ind_roi,ind_seeds,:),2);
        results.network_corr(n,t,:) = results.network_corr(n,t,:) + tmp_Z;
      end;
      % divide by the number of ROIs for this network
      results.network_corr(n,t,:) = results.network_corr(n,t,:) / nrois;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_analysis(parms,results)
  if ~exist(parms.fname_csv,'file') || parms.forceflag
    if parms.verbose==2
      fprintf('%s: writing file %s...\n',mfilename,parms.fname_csv);
    end;
    % initialize data cell array
    data = cell(results.nsubs,0);
    col_labels = cell(1,0);
    % add SubjIDs as first column    
    data = cat(2,data,results.SubjIDs');
    col_labels = cat(2,col_labels,'SubjID');
    % add VisitIDs as second column
    data = cat(2,data,results.VisitIDs');
    col_labels = cat(2,col_labels,'VisitID');
    % add TR
    data = cat(2,data,num2cell(results.TR'));
    col_labels = cat(2,col_labels,'TR');
    % add nreps
    data = cat(2,data,num2cell(results.nreps'));
    col_labels = cat(2,col_labels,'nreps');
    % add numTRs
    data = cat(2,data,num2cell(results.numTRs'));
    col_labels = cat(2,col_labels,'numTRs');
    % add motion stats
    for i=1:length(parms.motion_tags)
      tag = parms.motion_tags{i};
      data = cat(2,data,num2cell(results.(tag)'));
      col_labels = cat(2,col_labels,tag);
    end;
    % add info from StudyInfo
    if parms.continfo_flag || parms.subjinfo_flag
      for t=1:length(parms.info_tags)
        tag = parms.info_tags{t};
        if isfield(results.StudyInfo,tag)
          info = {results.StudyInfo.(tag)}';
          data = cat(2,data,info);
          col_labels = cat(2,col_labels,tag);
        end;
      end;
    end;
    % add var values
    if parms.var_flag
      data = cat(2,data,num2cell(results.roi_var'));
      tmp_labels = cell(1,results.nrois);
      for r=1:results.nrois
        tmp_labels{r} = ['var-' results.roinames{r}];
      end;
      col_labels = cat(2,col_labels,tmp_labels);   
    end;
    % add ROI values
    if parms.roi_flag
      % reshape corr matrix and add to data cell array
      data = cat(2,data,...
        num2cell(reshape(results.Z,...
                        [results.nrois*results.nseeds,results.nsubs])'));
      tmp_labels = cell(1,results.nrois*results.nseeds);
      k = 1;
      for s=1:results.nseeds
        seedname = results.seed_roinames{s};
        for r=1:results.nrois
          roiname = results.roinames{r};
          tmp_labels{k} = [seedname '->' roiname];
          k = k + 1;
        end;
      end;
      col_labels = cat(2,col_labels,tmp_labels);
    end;
    % add network values
    if parms.network_flag
      % reshape corr matrix and add to data cell array
      data = cat(2,data,...
        num2cell(reshape(results.network_corr,...
                        [results.nn*results.nt,results.nsubs])'));
      tmp_labels = cell(1,results.nn*results.nt);
      k = 1;
      for n=1:results.nn
        network_name = results.network_names{n};
        switch parms.network_flag
          case 1
            targ_names = results.roinames;
          case 2
            targ_names = {network_name};
          case 3
            targ_names = results.network_names;
        end;
        for t=1:results.nt
          targ_name = targ_names{t};
          tmp_labels{k} = [network_name '->' targ_name];
          k = k + 1;
        end;
      end;
      col_labels = cat(2,col_labels,tmp_labels);
    end;
    % add column labels
    data = cat(1,col_labels,data);
    % write file
    mmil_write_csv(parms.fname_csv,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function csv_outfix = set_csv_outfix(parms)
  csv_outfix = parms.roi_type;
  if ~isempty(parms.outfix)
    csv_outfix = [csv_outfix '_' parms.outfix];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat_outfix = set_mat_outfix(parms)
  mat_outfix = [parms.roi_infix '_roi_data'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


