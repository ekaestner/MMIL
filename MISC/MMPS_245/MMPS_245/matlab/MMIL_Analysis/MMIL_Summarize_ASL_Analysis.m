function results = MMIL_Summarize_ASL_Analysis(ProjID,varargin)
% function results = MMIL_Summarize_ASL_Analysis(ProjID,[options])
%
% Purpose: create summary spreadsheets with ASL ROI measures
%
% Usage:
%  results = MMIL_Summarize_ASL_Analysis(ProjID,'key1', value1,...);
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
%      proc_ASL, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%
% Optional Parameters that specify which analysis results to compile:
%   'measlist': list of ASL "measures" to extract for each ROI
%      e.g. 'CBF'
%      'CBF' is cerebral blood flow
%     {default = {'CBF'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     if empty, will use 1 for each
%     if not empty, length must match length of measlist
%     {default = []}
%   'aseg_flag': [0|1] compile aseg ROI results
%     {default = 1}
%   'cortsurf_flag': [0|1] compile cortical surface ROI results
%     {default = 1}
%   'area_flag': [0|1] when calculating mean values, weight by cortical area
%     requires that MRI analysis have been run
%     {default = 1}
%   'concat_flag': [0|1|2] summarize concatenated results for all analysis types
%     if 0, create separate summaries for each analysis type and measure
%     if 1, create a single summary for all analysis types and measures
%     if 2, create both individual and concatenated summaries
%     {default = 1}
%   'concat_outfix': string attached to concatenated summary file
%     {default = 'all'}
%
% Optional Parameters specific to aseg analysis:
%   'erode_flag': [0|1] whether aseg ROIs were "eroded" by
%     smoothing and thresholding to reduce partial voluming
%     {default = 1}
%   'erode_nvoxels': number of voxels eroded
%     {default = 1}
%   'aseg_aparc_flag': [0|1|2] whether cortical parcellation volume ROIs were used
%     0: aseg only
%     1: aparc only
%     2: aparc+aseg
%     {default = 0}
%
% Optional Parameters specific to cortical surface analysis:
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     {default = []}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast maps
%     {default = 1}
%   'gwnorm_flag': [0|1] when calculating gray-white contrast, normalize by mean
%     {default = 1}
%
% Other Optional Parameters:
%   'analdir': name of ASL analysis folder relative to ContainerPath
%     {default = 'analysis'}
%   'continfo_flag': [0|1] include container-related information
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%   'subjinfo_flag': [0|1] include subject information
%     e.g. Age, Sex, Site, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'outdir': output directory for summary csv files
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ROI_Summaries'}
%   'outstem': output file stem
%     {default = 'ASL'}
%   'save_mat_flag': [0|1] save compiled ROI data in a mat file
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/28/13 by Don Hagler
% Last Mod: 08/29/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

% check input parameters
parms = check_input(ProjID,varargin);

% check whether anything to do
if ~parms.aseg_flag && ~parms.cortsurf_flag
  fprintf('%s: nothing to do\n',mfilename);
  return;
end;

% check whether output files exist
if ~parms.forceflag
  runflag = check_output(parms);
  if ~runflag, return; end;
end;

% compile ASL analysis results for each subject
results = compile_results(parms);

% create summary csv files for individual roitypes and measuress
if parms.concat_flag~=1
  if parms.aseg_flag
    summarize_roitype_analysis(parms,results,'aseg');
  end;
  if parms.cortsurf_flag
    for p=1:parms.nprojdist
      summarize_roitype_analysis(parms,results,'cortsurf',p);
    end;
    if parms.contrast_flag
      summarize_roitype_analysis(parms,results,'contrast');
    end;
  end;
end;

% create summary csv files for all roitypes and measuress
if parms.concat_flag
  summarize_concat_analysis(parms,results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'QC_recon',true,[false true],...
    'qcflag',true,[false true],...
... % specify which results to compile
    'measlist',{'CBF'},[],...
    'scalefacts',[],[],...
    'aseg_flag',true,[false true],...
    'cortsurf_flag',true,[false true],...
    'area_flag',true,[false true],...
    'concat_flag',1,[0:2],...
    'concat_outfix','all',[],...
... % aseg analysis
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
... % cortsurf analysis
    'fnames_fparc',[],[],...
    'projdist_list',[1],[-5,5],...
    'gwnorm_flag',true,[false true],...
... % other
    'analdir','analysis',[],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'verbose',1,[0:2],...
    'outdir','ROI_Summaries',[],...
    'outstem','ASL',[],...
    'save_mat_flag',false,[false true],...
    'forceflag',false,[false true],...
... % hidden
    'required_containers',{'proc_asl','fsurf'},[],...
    'fname_fscolorlut',[],[],...
    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[],...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'excl_tags',{'excl_tags','ProjID',...
                 'outdir','outstem',...
                 'concat_outfix','save_mat_flag',...
                 'forceflag','combined_measlist','info_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  ProjInfo = MMIL_Get_ProjInfo(ProjID);
  
  parms.compile_tags = setdiff(fieldnames(parms),parms.excl_tags);

  if mmil_isrelative(parms.outdir)
    parms.outdir = [getenv('HOME') '/MetaData/' parms.ProjID '/' parms.outdir];
  end;

  if parms.cortsurf_flag
    % whether to calculate gray-white contrast depends on if 2 projdist
    parms.nprojdist = length(parms.projdist_list);
    if parms.nprojdist==2
      parms.contrast_flag = 1;
      parms.projdist_list = sort(parms.projdist_list); % white matter first
    else
      parms.contrast_flag = 0;
    end;
  else
    parms.contrast_flag = 0;
  end;

  parms.nmeas = length(parms.measlist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runflag = check_output(parms)
  runflag = 0;
  if parms.concat_flag~=1
    if parms.aseg_flag
      all_exist = check_output_files(parms,'aseg');
      if ~all_exist
        runflag = 1;
        return;
      end;
    end;
    if parms.cortsurf_flag
      for p=1:parms.nprojdist
        all_exist = check_output_files(parms,'cortsurf',p);
        if ~all_exist
          runflag = 1;
          return;
        end;
      end;
      if parms.contrast_flag
        all_exist = check_output_files(parms,'contrast');
        if ~all_exist
          runflag = 1;
          return;
        end;
      end;
    end;
  end;
  if parms.concat_flag
    outfix = set_concat_outfix(parms);
    fname = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
    if ~exist(fname,'file')
      runflag = 1;
      return;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_exist = check_output_files(parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  all_exist = 1;
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    switch roitype
      case 'aseg'
        outfix = set_aseg_outfix(parms);
      case 'cortsurf'  
        outfix = set_cortsurf_outfix(parms,p);
      case 'contrast'  
        outfix = set_cortsurf_outfix(parms,0);
    end;
    fname = sprintf('%s/%s_%s_%s.csv',parms.outdir,parms.outstem,meas,outfix);
    if ~exist(fname,'file')
      all_exist = 0;
      return;
    end;
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(parms)
  fname_out = sprintf('%s/%s_roi_data.mat',parms.outdir,parms.outstem);
  if ~exist(fname_out,'file') || parms.forceflag || ~parms.save_mat_flag
    mmil_mkdir(parms.outdir);
    parms.concat_flag = (parms.concat_flag>0);
    args = mmil_parms2args(parms,parms.compile_tags);
    results = MMIL_Compile_ASL_Analysis(parms.ProjID,args{:});
    if parms.save_mat_flag, save(fname_out,'results'); end;
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(fname,results,parms,roitype,meas,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  if ~exist('meas','var'), meas = []; end;
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
  if strcmp(roitype,'all')
    tmp_data = results.(roitype).data;
  else
    tmp_data = results.(roitype).(meas).data(:,:,p);
  end;
  data = cat(2,data,num2cell(tmp_data));
  col_labels = cat(2,col_labels,mmil_rowvec(results.(roitype).roinames));
  % add column labels
  data = cat(1,col_labels,data);
  % write file
  mmil_write_csv(fname,data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_roitype_analysis(parms,results,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    switch roitype
      case 'aseg'
        outfix = set_aseg_outfix(parms);
      case 'cortsurf'  
        outfix = set_cortsurf_outfix(parms,p);
      case 'contrast'  
        outfix = set_cortsurf_outfix(parms,0);
    end;
    fname = sprintf('%s/%s_%s_%s.csv',parms.outdir,parms.outstem,meas,outfix);
    if ~exist(fname,'file') || parms.forceflag
      mmil_mkdir(parms.outdir);
      if parms.verbose==2
        fprintf('%s: writing %s...\n',mfilename,fname);
      end;
      write_csv(fname,results,parms,roitype,meas,p);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_concat_analysis(parms,results)
  outfix = set_concat_outfix(parms);
  fname = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
  if ~exist(fname,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    if parms.verbose==2
      fprintf('%s: writing %s...\n',mfilename,fname);
    end;
    write_csv(fname,results,parms,'all');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_aseg_outfix(parms)
  switch parms.aseg_aparc_flag
    case 0
      outfix = 'aseg';
    case 1
      outfix = 'aparc';
    case 2
      outfix = 'aparc+aseg';
  end;
  if parms.erode_flag
    outfix = [outfix '_erode'];
  end
  if parms.aseg_disp_flag
    outfix = sprintf('%s_%s%0.1f',...
      outfix,parms.aseg_disp_suffix,parms.aseg_dispfact);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_cortsurf_outfix(parms,p)
  outfix = 'cortsurf';
  if p==0
    outfix = [outfix '_contrast'];
  else
    projdist = parms.projdist_list(p);
    if projdist<0
      outfix = sprintf('%s_white%0.1f',outfix,projdist);
    else
      outfix = sprintf('%s_gray+%0.1f',outfix,projdist);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_concat_outfix(parms)
  outfix = parms.concat_outfix;
return;

