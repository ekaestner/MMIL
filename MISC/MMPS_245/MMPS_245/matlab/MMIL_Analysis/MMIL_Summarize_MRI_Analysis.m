function results = MMIL_Summarize_MRI_Analysis(ProjID,varargin)
% function results = MMIL_Summarize_MRI_Analysis(ProjID,[options])
%
% Purpose: create summary spreadsheets with MRI ROI measures
%
% Usage:
%  results = MMIL_Summarize_MRI_Analysis(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information:
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
%  'cortT1w_flag': average T1-weighted values for aparc or fuzzy cluster ROIs
%    also calculate gray-white contrast
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
%  'aparc_flag': [0|1] use aparc stats
%      for thickness, sulc, area, cortvol, and T1w
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
%  'concat_flag': [0|1|2] summarize concatenated results for all analysis types
%     if 0, create separate summaries for each analysis type
%     if 1, create a single summary for all analysis types
%     if 2, create both individual and concatenated summaries
%    {default = 1}
%
% Other Optional Parameters:
%  'proc_outputlist': list of output file stems; e.g. {'T1w','T2w'}
%       if proc_outputlist contains 'T1w'
%         cortT1w_flag and subT1w_flag are irrelevant
%     {default = []}
%  'continfo_flag': [0|1] include container-related information
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%  'subjinfo_flag': [0|1] include subject information
%     e.g. Age, Sex, Site, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%  'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-0.2,0.2]}
%  'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%  'outdir': output directory for summary csv files
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ROI_Summaries'}
%  'outstem': output file stem
%     {default = 'MRI'}
%  'save_mat_flag': [0|1] save compiled ROI data in a mat file
%     {default = 0}
%  'check_complete_flag': [0|1] whether to require that recon is complete
%     {default = 1}
%  'FS_version': which version of Freesurfer used (e.g. 305, 450, 510, 530)
%    for checking whether recon is complete
%    if empty, will use FREESURER_VER environment variable
%       or get from ContainerInfo
%    {default = []}
%  'baseflag': [0|1] only include baseline visits (VisitNumber = 1)
%     ignored if VisitNumber not specified in VisitInfo
%     {default = 1}
%  'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
%   NOTE: subjects with incomplete FreeSurfer recons will be excluded
%
% Created:  05/23/09 by Don Hagler
% Prev Mod: 02/01/16 by Don Hagler
% Last Mod: 07/17/17 by Don Hagler
%

%% todo: requires testing! and changes to MMIL_Compile_MRI_Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

% check input parameters
parms = check_input(ProjID,varargin);

% check whether anything to do
if ~parms.thick_flag ~parms.sulc_flag && ~parms.area_flag &&...
   ~parms.cortvol_flag && ~parms.cortT1w_flag && ~parms.cortproc_flag &&...
   ~parms.subvol_flag && ~parms.subT1w_flag && ~parms.subproc_flag
  fprintf('%s: nothing to do\n',mfilename);
  return;
end;

% check whether output files exist
if ~parms.forceflag
  runflag = check_output(parms);
  if ~runflag, return; end;
end;

% compile MRI analysis results for each subject
results = compile_results(parms);

% create summary csv files for individual roitypes and measuress
if parms.concat_flag~=1
  if parms.thick_flag
    summarize_analysis(parms,results,'thick');
  end;
  if parms.sulc_flag
    summarize_analysis(parms,results,'sulc');
  end;
  if parms.area_flag
    summarize_analysis(parms,results,'area');
  end;
  if parms.cortvol_flag
    summarize_analysis(parms,results,'cortvol');
  end;
  if parms.cortproc_flag
    for n=1:parms.ninputs
      fname = sprintf('cort%s',parms.proc_outputlist{n});
      for p=1:parms.nprojdist
        summarize_analysis(parms,results,fname,p);
      end;
      if parms.contrast_flag
        fname = sprintf('contrast%s',parms.proc_outputlist{n});
        summarize_analysis(parms,results,fname);
      end;
    end;
  end;
  if parms.subvol_flag
    summarize_analysis(parms,results,'subvol');
  end;
  if parms.subproc_flag
    for n=1:parms.ninputs
      fname = sprintf('sub%s',parms.proc_outputlist{n});
      summarize_analysis(parms,results,fname);
    end;
  end;
end;

% create summary csv files for all analysis types
if parms.concat_flag
  summarize_analysis(parms,results,'all');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'QC_raw',true,[false true],...
    'QC_recon',true,[false true],...
    'qcflag',true,[false true],...
... % specify which results to compile
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
    'concat_flag',1,[0:2],...
...
    'proc_outputlist',[],[],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'projdist_list',[-0.2,0.2],[-5,5],...
    'verbose',1,[0:2],...
    'outdir','ROI_Summaries',[],...
    'outstem','MRI',[],...
    'save_mat_flag',false,[false true],...
    'check_complete_flag',true,[false true],...
    'FS_version',[],[],...
    'baseflag',true,[false true],...
    'forceflag',false,[false true],...
...
    'erode_flag',false,[false true],...
    'erode_nvoxels',1,[1:100],...
    'analysis_outdir','analysis',[],...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'required_containers',{'fsurf'},[],...
    'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,251:255],[1,Inf],...
    'extra_roilist',[10001:10003,20001:20003,20009,20010],[1,Inf],...
    'aparc_roilist',[1001:1003,1005:1035,2001:2003,2005:2035],[1,Inf],...
    'fname_fscolorlut',[],[],...
...
    'excl_tags',{'excl_tags','ProjID',...
                 'outdir','outstem','save_mat_flag',...
                 'forceflag','info_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  parms.compile_tags = setdiff(fieldnames(parms),parms.excl_tags);

  % check outdir
  if mmil_isrelative(parms.outdir)
    parms.outdir = [getenv('HOME') '/MetaData/' parms.ProjID '/' parms.outdir];
  end;

  parms.nprojdist = length(parms.projdist_list);
  if parms.nprojdist==2
    parms.contrast_flag = 1;
    parms.projdist_list = sort(parms.projdist_list); % white matter first
  else
    parms.contrast_flag = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runflag = check_output(parms)
  runflag = 0;
  if parms.concat_flag~=1
    if parms.thick_flag
      outfix = set_outfix(parms,'thick');
      fname_out = sprintf('%s/%s_%s.csv',...
        parms.outdir,parms.outstem,outfix);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
    end;
    if parms.sulc_flag
      outfix = set_outfix(parms,'sulc');
      fname_out = sprintf('%s/%s_%s.csv',...
        parms.outdir,parms.outstem,outfix);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
    end;
    if parms.area_flag
      outfix = set_outfix(parms,'area');
      fname_out = sprintf('%s/%s_%s.csv',...
        parms.outdir,parms.outstem,outfix);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
    end;
    if parms.cortvol_flag
      outfix = set_outfix(parms,'cortvol');
      fname_out = sprintf('%s/%s_%s.csv',...
        parms.outdir,parms.outstem,outfix);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
    end;
    if parms.cortT1w_flag
      for p=1:parms.nprojdist
        projdist = parms.projdist_list(p);
        outfix = set_outfix(parms,'cortT1w',p);
        fname_out = sprintf('%s/%s_%s.csv',...
          parms.outdir,parms.outstem,outfix);
        if ~exist(fname_out,'file'), runflag = 1; return; end;
      end;
      if parms.contrast_flag
        outfix = set_outfix(parms,'contrast');
        fname_out = sprintf('%s/%s_%s.csv',...
          parms.outdir,parms.outstem,outfix);
        if ~exist(fname_out,'file'), runflag = 1; return; end;
      end;
    end;
    if parms.cortproc_flag
      for n=1:parms.ninputs
        fname = sprintf('cort%s',parms.proc_outputlist{n});
        for p=1:parms.nprojdist
          projdist = parms.projdist_list(p);
          outfix = set_outfix(parms,fname,p);
          fname_out = sprintf('%s/%s_%s.csv',...
            parms.outdir,parms.outstem,outfix);
          if ~exist(fname_out,'file'), runflag = 1; return; end;
        end;
        if parms.contrast_flag
          fname = sprintf('contrast%s',parms.proc_outputlist{n});
          outfix = set_outfix(parms,fname);
          fname_out = sprintf('%s/%s_%s.csv',...
            parms.outdir,parms.outstem,outfix);
          if ~exist(fname_out,'file'), runflag = 1; return; end;
        end;
      end;
    end;
    if parms.subvol_flag
      outfix = set_outfix(parms,'subvol');
      fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
    end;
    if parms.subT1w_flag
      outfix = set_outfix(parms,'subT1w');
      fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
    end;
    if parms.subproc_flag
      for n=1:parms.ninputs
        fname = sprintf('sub%s',parms.proc_outputlist{n});
        outfix = set_outfix(parms,fname);
        fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
        if ~exist(fname_out,'file'), runflag = 1; return; end;
      end;
    end;
  end;
  if parms.concat_flag
    outfix = set_outfix(parms,'all');
    fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
    if ~exist(fname_out,'file'), runflag = 1; return; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(parms)
  fname_out = sprintf('%s/%s_roi_data.mat',parms.outdir,parms.outstem);
  if ~exist(fname_out,'file') || parms.forceflag || ~parms.save_mat_flag
    mmil_mkdir(parms.outdir);
    parms.concat_flag = (parms.concat_flag>0);
    args = mmil_parms2args(parms,parms.compile_tags);
    results = MMIL_Compile_MRI_Analysis(parms.ProjID,args{:});
    if parms.save_mat_flag, save(fname_out,'results'); end;
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(fname,results,parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  fieldname = set_fieldname(parms,roitype);
  % initialize data cell array
  data = cell(results.nsubs,0);
  col_labels = cell(1,0);
  % add SubjIDs as first column    
  data = cat(2,data,results.SubjIDs');
  col_labels = cat(2,col_labels,'SubjID');
  % add VisitIDs as second column
  data = cat(2,data,results.VisitIDs');
  col_labels = cat(2,col_labels,'VisitID');
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
  % add ROI data
  tmp_data = results.(fieldname).data(:,:,p);
  data = cat(2,data,num2cell(tmp_data));
  col_labels = cat(2,col_labels,mmil_rowvec(results.(fieldname).roinames));
  % add column labels
  data = cat(1,col_labels,data);
  % write file
  mmil_write_csv(fname,data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_analysis(parms,results,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  outfix = set_outfix(parms,roitype,p);
  fname = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
  if ~exist(fname,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    if parms.verbose==2
      fprintf('%s: writing %s...\n',mfilename,fname);
    end;
    write_csv(fname,results,parms,roitype,p);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_outfix(parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  if parms.ninputs>0
    proc_types = cellfun(@(x) sprintf('cort%s',x),parms.proc_outputlist,...
      'UniformOutput',false);
    contrast_types = cellfun(@(x) sprintf('contrast%s',x),parms.proc_outputlist,...
      'UniformOutput',false);
    sub_types = cellfun(@(x) sprintf('sub%s',x),parms.proc_outputlist,...
      'UniformOutput',false);
  else
    proc_types = {'cortT1w'};
    contrast_types = {'contrast'};
    sub_types = {'subT1w'};
  end;
  sub_types{end+1} = 'subvol';
  switch roitype
    case 'thick'
      outfix = 'cort_thick';
    case 'sulc'
      outfix = 'cort_sulc';
    case 'area'
      outfix = 'cort_area';
    case 'cortvol'
      outfix = 'cort_vol';
    case proc_types
      projdist = parms.projdist_list(p);
      prefix = regexprep(roitype,'cort','cort_');
      if projdist<0
        outfix = sprintf('%s_white%0.1f',prefix,projdist);
      else
        outfix = sprintf('%s_gray+%0.1f',prefix,projdist);
      end;
    case contrast_types
      prefix = regexprep(roitype,'contrast','cort_contrast');
    case sub_types
      prefix = regexprep(roitype,'sub','subcort_');
    case 'all'
      outfix = 'all';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fieldname,p] = set_fieldname(parms,roitype)
  if parms.ninputs>0
    proc_types = cellfun(@(x) sprintf('cort%s',x),parms.proc_outputlist,...
      'UniformOutput',false);
    contrast_types = cellfun(@(x) sprintf('contrast%s',x),parms.proc_outputlist,...
      'UniformOutput',false);
    sub_types = cellfun(@(x) sprintf('sub%s',x),parms.proc_outputlist,...
      'UniformOutput',false);
  else
    proc_types = {'cortT1w'};
    contrast_types = {'contrast'};
    sub_types = {'subT1w'};
  end;
  sub_types{end+1} = 'subvol';
  switch roitype
    case 'thick'
      fieldname = 'cort_thick';
    case 'sulc'
      fieldname = 'cort_sulc';
    case 'area'
      fieldname = 'cort_area';
    case 'cortvol'
      fieldname = 'cort_vol';
    case proc_types
      fieldname = regexprep(roitype,'cort','cort_');
    case contrast_types
      fieldname = regexprep(roitype,'contrast','cort_contrast_');
    case sub_types
      fieldname = regexprep(roitype,'sub','subcort_');
    case 'all'
      fieldname = 'all';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
