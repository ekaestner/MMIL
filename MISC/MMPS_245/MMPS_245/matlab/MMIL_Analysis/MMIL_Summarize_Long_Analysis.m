function MMIL_Summarize_Long_Analysis(ProjID,varargin)
% function MMIL_Summarize_Long_Analysis(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
% 
% Optional parameters that specify project specific information
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
% Optional parameters:
%  'fstem': file stem of output spreadsheet files (csv format)
%     {default = 'long_dv'}
%  'outdir': output directory
%     Created in user's home directory if not full path
%     if empty, will have name like: 'MetaData/ROI_Summaries' or
%                                    'MetaData/ProjID/ROI_Summaries'
%     {default = []}
%  'analysis_outdir': where analyzed mat files can be found placed relative
%     to Container dirs
%     {default = 'analysis'}
%  'baseflag': [0|1] register all visits to visit 1
%      otherwise, register all visits to all earlier visits
%     {default: 1}
%  'regtypes': nonlinear registration types (e.g. 'Fine', 'ROI')
%     {default = {'Fine','ROI'}}
%  'aseg_flag': use aseg ROIs (non-cortical, volume segmentation)
%     {default = 1}
%  'aparc_flag': use aparc ROIs (cortical surface parcellation)
%     {default = 1}
%  'aseg_roigroups_flag': [0|1] use masks for groups of aseg roi codes
%     includes 'WholeBrain', 'LatVentricles', and 'AllVentricles'
%     {default = 1}
%  'subhippo_flag': use subdivided hippocampal ROIs
%     {default = 1}
%  'indiv_flag': create individual spreadsheets for aseg, aparc, and hippo
%     Otherwise, only create one spreadsheet with all ROIs
%     {default = 0}
%  'baseline_data_flag': [0|1] include interleaved baseline data columns
%     {default = 1}
%  'followup_data_flag': [0|1] calculate followup values based on
%      baseline and change. if 0, output fractional change
%     {default = 0}
%  'forceflag': overwrite existing output
%     {default = 0}
%
%   NOTE: subjects with incomplete baseline FreeSurfer recons will be excluded
%
% Created:  03/02/10 by Don Hagler
% Last Mod: 03/17/17 by Don Hagler
%

%% todo: selectively combine results from Fine and ROI regtypes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'qcflag',true,[false true],...
...
  'fstem','long_dv',[],...
  'outdir',[],[],...
  'analysis_outdir','analysis',[],...
  'baseflag',true,[false true],...
  'regtypes',{'Fine','ROI'},[],...
  'aseg_flag',true,[false true],...
  'aparc_flag',true,[false true],...
  'subhippo_flag',true,[false true],...
  'aseg_roigroups_flag',true,[false true],...
  'indiv_flag',false,[false true],...
  'baseline_data_flag',true,[false true],...
  'followup_data_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,10001:10003,20001:20003,20009,20010],[1,Inf],...
  'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
  'hippo_roilist',[10011:10016],[1,Inf],...
  'fname_fscolorlut',[],[],...
  'dirprefix','LONG',[],...
  'erode_flag',true,[false true],...
  'erode_nvoxels',1,[1:100],...
  'SubjID_label','SubjID',[],...
...
  'required_containers',{'proc'},[],...
  'required_rootdirs',{'proc','fsurf','long'},[],...
  'modality','MRI',[],...
...
  'compile_tags',{'StudyInfo','qcflag','outdir','regtype',...
                  'aseg_flag','aparc_flag','subhippo_flag',...
                  'aseg_roigroups_flag','baseflag',...
                  'aseg_roilist','aparc_roilist',...
                  'hippo_roilist','fname_fscolorlut',...
                  'dirprefix','erode_flag','erode_nvoxels'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~parms.aseg_flag && ~parms.subhippo_flag && ~parms.aparc_flag
  fprintf('%s: WARNING: nothing to do\n',mfilename);
  return;
end;

% check outdir
if isempty(parms.outdir)
  parms.outdir = [getenv('HOME') '/MetaData'];
  if ~isempty(ProjID)
    parms.outdir = [parms.outdir '/' ProjID];
  end;
  parms.outdir = [parms.outdir '/ROI_Summaries'];
elseif mmil_isrelative(parms.outdir)
  parms.outdir = [getenv('HOME') '/' parms.outdir];
end;
% create output directory
mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r=1:length(parms.regtypes)
  regtype = parms.regtypes{r};

  if ~parms.forceflag
    fname_list = [];
    if parms.indiv_flag
      if parms.aseg_flag
        fname_list{end+1} = ...
          sprintf('%s/%s_%s_SubCort.csv',parms.outdir,parms.fstem,regtype);
      end;
      if parms.aseg_flag
        fname_list{end+1} = ...
          sprintf('%s/%s_%s_Cort.csv',parms.outdir,parms.fstem,regtype);
      end;
      if parms.aseg_flag
        fname_list{end+1} = ...
          sprintf('%s/%s_%s_Hippo.csv',parms.outdir,parms.fstem,regtype);
      end;
    else
      fname_list{end+1} = ...
        sprintf('%s/%s_%s_All.csv',parms.outdir,parms.fstem,regtype);
    end;
    runflag = 0;
    for f=1:length(fname_list)
      if ~exist(fname_list{f},'file'), runflag=1; break; end;
    end;
  else
    runflag = 1;
  end;
  if ~runflag, continue; end;  

  fprintf('%s: Compiling longitudinal analysis for regtype %s...\n',...
    mfilename,regtype);
  tmp_parms = parms;
  tmp_parms.outdir = parms.analysis_outdir;
  tmp_parms.regtype = parms.regtypes{r};
  args = mmil_parms2args(tmp_parms,parms.compile_tags);
  results = MMIL_Compile_Long_Analysis(ProjID,args{:});

  SubjIDs = results.SubjIDs;
  VisitNumbers = [results.StudyInfo.VisitNumber];
  StudyInfo_SubjIDs = {results.StudyInfo.SubjID};
  nrows = length(SubjIDs);
  stind = zeros(nrows,2);
  for i=1:nrows
    ind_tmp = find(VisitNumbers==results.VisitNumbers_Baseline{i} &...
                      strcmp(StudyInfo_SubjIDs,SubjIDs{i}));
    if length(ind_tmp)>1
      fprintf('%s: WARNING: duplicate baseline visit for %s; using %s not %s\n',...
        mfilename,results.StudyInfo(ind_tmp(1)).SubjID,...
        results.StudyInfo(ind_tmp(1)).VisitID,...
        results.StudyInfo(ind_tmp(2)).VisitID);
      ind_tmp = ind_tmp(1);
    end;
    stind(i,1) = ind_tmp;
    ind_tmp = find(VisitNumbers==results.VisitNumbers_Followup{i} &...
                      strcmp(StudyInfo_SubjIDs,SubjIDs{i}));
    if length(ind_tmp)>1
      fprintf('%s: WARNING: duplicate followup visit %d for %s; using %s not %s\n',...
        mfilename,results.VisitNumbers_Followup{i},...
        results.StudyInfo(ind_tmp(1)).SubjID,...
        results.StudyInfo(ind_tmp(1)).VisitID,...
        results.StudyInfo(ind_tmp(2)).VisitID);
      ind_tmp = ind_tmp(1);
    end;
    stind(i,2) = ind_tmp;
  end
  visit_data = cat(2,results.VisitNumbers_Baseline,...
                     results.VisitNumbers_Followup,...
                     results.StudyDates_Baseline,...
                     results.StudyDates_Followup,...
                     {results.StudyInfo(stind(:,1)).VisitID}',...
                     {results.StudyInfo(stind(:,2)).VisitID}');
  visit_labels = {'VisitNumber_Baseline', 'VisitNumber_Followup',...
                  'StudyDate_Baseline', 'StudyDate_Followup',...
                  'VisitID_Baseline', 'VisitID_Followup'};

  if parms.indiv_flag
    % Volume ROIs
    if parms.aseg_flag
      write_data(parms,results,visit_data,visit_labels,regtype,'SubCort')
    end;
    % Cortical Surface ROIs
    if parms.aparc_flag
      write_data(parms,results,visit_data,visit_labels,regtype,'Cort')
    end;
    % Subdivided Hippocampal ROIs
    if parms.subhippo_flag
      write_data(parms,results,visit_data,visit_labels,regtype,'Hippo')
    end;
  else
    % All ROIs
    write_data(parms,results,visit_data,visit_labels,regtype,'All')
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_data(parms,results,visit_data,visit_labels,regtype,roitype)
  
  fname = sprintf('%s/%s_%s_%s.csv',parms.outdir,parms.fstem,regtype,roitype);
  if ~exist(fname,'file') || parms.forceflag
    fprintf('%s: writing %s...\n',mfilename,fname);
    tmp_labels = results.(roitype).roinames';
    if parms.baseline_data_flag
      tmp_data1 = results.(roitype).data_baseline;
      if parms.followup_data_flag
        tmp_data2 = tmp_data1.*(1 + results.(roitype).data_change);
        tmp1 = cellfun(@(x) sprintf('%s-baseline',x),tmp_labels,'UniformOutput',false);
        tmp2 = cellfun(@(x) sprintf('%s-followup',x),tmp_labels,'UniformOutput',false);
      else
        tmp_data2 = results.(roitype).data_change;
        tmp1 = cellfun(@(x) sprintf('%s-baseline',x),tmp_labels,'UniformOutput',false);
        tmp2 = cellfun(@(x) sprintf('%s-change',x),tmp_labels,'UniformOutput',false);
      end;
      tmp = intlv_data(tmp1,tmp2);
      tmp_labels = reshape(tmp',[1,numel(tmp)]);
      tmp_data = num2cell(intlv_data(tmp_data1,tmp_data2));
    else
      if parms.followup_data_flag
        tmp_data = results.(roitype).data_baseline.*...
                    (1 + results.(roitype).data_change);
        tmp_labels = mmil_splitstr(sprintf('%s-followup\n',tmp_labels{:}),'\n')';
      else	
        tmp_data = results.(roitype).data_change;
        tmp_labels = mmil_splitstr(sprintf('%s-change\n',tmp_labels{:}),'\n');
      end;
    end;
    tmp_data = cat(2,visit_data,tmp_data);
    tmp_labels = cat(2,visit_labels,tmp_labels);
    mmil_write_csv(fname,tmp_data,...
        'row_labels', results.SubjIDs,...
        'col_labels', tmp_labels,...
        'firstcol_label', parms.SubjID_label);
  end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = intlv_data(data1,data2)
  [nr,nc] = size(data1);
  if iscell(data1)
    data = cell(nr,2*nc);
  else
    data = zeros(nr,2*nc);
  end;
  data(:,1:2:2*nc) = data1;
  data(:,2:2:2*nc) = data2;
return;


