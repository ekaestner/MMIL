function results = MMIL_Compile_ASL_Analysis(ProjID,varargin)
%function results = MMIL_Compile_ASL_Analysis(ProjID,[options])
%
% Usage:
%  results = MMIL_Compile_ASL_Analysis(ProjID,'key1', value1,...);
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
%   'QC_raw' : [0|1] use good raw QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RawQC.mat
%     {default = 1}
%   'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%   'QC_ASL' : [0|1] use ASLreg QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_ASLQC.csv
%     {default = 1}
%   'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
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
%     {default = 0}
%   'cortsurf_flag': [0|1] compile cortical surface ROI results
%     {default = 0}
%   'area_flag': [0|1] when calculating mean values, weight by cortical area
%     requires that MRI analysis have been run
%     {default = 1}
%   'concat_flag': [0|1] concatenate results for all analysis types
%     {default = 0}
%
% Optional Parameters specific to how aseg analysis was done:
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
% Optional Parameters specific to how cortical surface analysis was done:
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     {default = []}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-1,1]}
%   'gwnorm_flag': [0|1] when calculating gray-white contrast, normalize by mean
%     {default = 1}
%
% Other Optional Parameters:
%   'analdir': name of ASL analysis folder relative to ContainerPath
%     {default = 'analysis'}
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
%   results: struct containing compiled ASL analysis results with fields:
%     StudyInfo: struct array containing info for each subject
%     RootDirs: struct containing paths of root directories
%     SubjIDs: cell array of subject IDs
%     VisitIDs: cell array of visit IDs
%     STRUCT_VisitIDs: cell array of structural (i.e. FreeSurfer) visit IDs
%     VisitNumbers: vector of visit numbers
%     nsubs: number of subjects/visits in StudyInfo
%     nmeas: number of measures in measlist
%     aseg       (if aseg_flag = 1)
%     cortsurf   (if cortsurf_flag = 1)
%     all   (if concat_flag = 1)
%       nrois: number of ROIs of all types
%       roicodes: ROI codes
%       roinames: ROI names with meas prepended
%       data: matrix of ROI averages for all measures and ROI types
%       std: matrix of standard deviations
%       nvals: matrix of number of voxels for each subject and ROI
%       nvals_valid: number of values included in ROI averages
%       nvals_invalid: number of values excluded from ROI averages
%
% Created:  02/28/13 by Don Hagler
% Last Mod: 08/29/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

results = init_results(parms);

if parms.aseg_flag
  results = compile_roi_data(parms,results,'aseg');
end;

if parms.cortsurf_flag
  results = compile_roi_data(parms,results,'cortsurf');
  % calculate difference between gray and white matter
  if parms.contrast_flag
    results = calc_contrast(parms,results);
  end;
  % compile cortical area for calculating weighted mean values
  if parms.area_flag
    results = compile_cort_area(parms,results);
  end;
  % calculate mean cortical values
  results = calc_mean_cort_vals(parms,results);
  if parms.contrast_flag
    results = calc_mean_cort_vals(parms,results,1);
  end;
end;

if parms.concat_flag
  results = concat_results(parms,results);
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
    'inputlist',[],[],...
    'aseg_flag',false,[false true],...
    'cortsurf_flag',false,[false true],...
    'area_flag',true,[false true],...
    'concat_flag',false,[false true],...
... % aseg analysis
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
... % cortsurf analysis
    'fnames_fparc',[],[],...
    'projdist_list',[-1,1],[-5,5],...
    'gwnorm_flag',true,[false true],...
... % other
    'analdir','analysis',[],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'verbose',1,[0:2],...
... % hidden
    'required_containers',{'proc_asl','fsurf'},[],...
    'fname_fscolorlut',[],[],...
    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[],...
...
    'info_tags',{'snums','revflag','min_nb0','min_ndirs',...
                 'min_bval','flex_flag'},[],...
    'code_tags',{'aseg_aparc_flag','aseg_roilist','aparc_roilist',...
                 'exclude_roilist'},[],...
    'compile_area_tags',{'StudyInfo','RootDirs','QC_raw','QC_recon',...
      'qcflag','thick_flag','area_flag','cortvol_flag','cortT1w_flag',...
      'subvol_flag','subT1w_flag','aparc_flag','fuzzy_flag','fuzzy_dir',...
      'fuzzy_fstem','fuzzy_order','concat_flag','continfo_flag',...
      'subjinfo_flag','analysis_outdir','projdist_list','check_complete_flag',...
      'verbose','baseflag','aseg_roilist','extra_roilist',...
      'aparc_roilist','fname_fscolorlut','required_containers','modality',...
      'erode_nvoxels','erode_flag','hemilist','fuzzy_name_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(StudyInfo), error('empty StudyInfo'); end;

  parms.StudyInfo = StudyInfo;
  parms.RootDirs = RootDirs;
  parms.nsubs = length(parms.StudyInfo);
  
  % get container-related info
  if parms.continfo_flag
    parms.StudyInfo = MMIL_Get_ContInfo(parms.StudyInfo,parms.RootDirs);
  end;

  % get subject-related info
  if parms.subjinfo_flag
    parms.StudyInfo = ...
      MMIL_Get_SubjInfo(parms.StudyInfo,parms.RootDirs,parms.ProjID);
  end;

  if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;
  parms.nmeas = length(parms.measlist);
  if isempty(parms.scalefacts)
    parms.scalefacts = ones(1,parms.nmeas);
  elseif numel(parms.scalefacts)~=parms.nmeas
    error('number of scalefacts (%d) does not match measlist (%d)',...
      numel(parms.scalefacts),parms.nmeas);
  end;

  if parms.aseg_flag || parms.cortsurf_flag
    % read roi codes and names from FreeSurfer color lookup table
    if isempty(parms.fname_fscolorlut)
      MMPS_parms = getenv('MMPS_PARMS');
      parms.fname_fscolorlut = [MMPS_parms '/MMIL_FSColorLUT.txt'];
    end;
    if ~exist(parms.fname_fscolorlut,'file')
      error('FreeSurfer color lookup table file %s not found',...
        parms.fname_fscolorlut);
    end;
    [fs_roicodes,fs_roinames] = fs_colorlut(parms.fname_fscolorlut);
  end;

  % set parameters for aseg analysis
  if parms.aseg_flag
    args = mmil_parms2args(parms,parms.code_tags);
    parms.aseg_roicodes = mmil_aseg_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aseg_roicodes);
    parms.aseg_roicodes = fs_roicodes(ind_fs);
    parms.aseg_roinames = fs_roinames(ind_fs);
    parms.aseg_nrois = length(parms.aseg_roicodes);
  end;

  % set parameters for cortsurf analysis
  if parms.cortsurf_flag
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aparc_roilist);
    parms.cortsurf_roicodes = fs_roicodes(ind_fs);
    parms.cortsurf_roinames = fs_roinames(ind_fs);
    if ~isempty(parms.fnames_fparc)
      [fparc_roicodes,fparc_roinames] = read_fparc_roinames(parms);
      parms.cortsurf_roicodes = cat(1,parms.cortsurf_roicodes,fparc_roicodes);
      parms.cortsurf_roinames = cat(1,parms.cortsurf_roinames,fparc_roinames);
    end;
    parms.cortsurf_nrois = length(parms.cortsurf_roicodes);
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
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roicodes,roinames] = read_fparc_roinames(parms)
  roicodes = []; roinames = [];
  nrois = 0;
  for f=1:length(parms.fnames_fparc)
    fname = parms.fnames_fparc{f};
    if parms.verbose==2
      fprintf('%s: loading annotation from %s...\n',mfilename,fname);
    end;
    n = regexp(fname,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
    if isempty(n)
      error('unexpected naming convention for annot file %s\n',fname);
    end;
    hemi = n.hemi;
    % read annotation file
    [annot_nums,annot_names] = fs_read_annotation(fname);
    roicode = 0;
    for i=1:length(annot_names)
      roicode = roicode + 1;    
      annot_name = annot_names{i};
      if strcmp(annot_name,'unknown'), continue; end;
      nrois = nrois + 1;  
      roinames{nrois} = ['ctx-' hemi '-' annot_name];
      roicodes(nrois) = roicode;
    end;
  end;
  roinames = mmil_colvec(roinames);
  roicodes = mmil_colvec(roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function infix = set_infix(parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  switch roitype
    case 'aseg'
      switch parms.aseg_aparc_flag
        case 0
          infix = 'aseg';
        case 1
          infix = 'aparc';
        case 2
          infix = 'aparc+aseg';
      end;
      if parms.erode_flag
        infix = sprintf('%s_erode',infix);
        if parms.erode_nvoxels>1,
          infix = sprintf('%s%d',infix,parms.erode_nvoxels);
        end;
      end;
      infix = [infix '_roi_data'];
    case 'cortsurf'
      infix = sprintf('cortsurf_pdist%0.1f_roi_data',...
        parms.projdist_list(p));
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
  results.nsubs = parms.nsubs;
  results.nmeas = parms.nmeas;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_roi_results(parms,results,roitype)
  results.(roitype) = [];
  nrois = parms.([roitype '_nrois']);
  if strcmp(roitype,'cortsurf') && parms.nprojdist>1
    np = parms.nprojdist;
  else
    np = 1;
  end;
  results.(roitype).nrois = nrois;
  results.(roitype).roicodes = parms.([roitype '_roicodes']);  
  results.(roitype).roinames = parms.([roitype '_roinames']);
  init_vals = nan(parms.nsubs,nrois,np);
  for m=1:parms.nmeas
    results.(roitype).(parms.measlist{m}).data = init_vals;
    results.(roitype).(parms.measlist{m}).std = init_vals;
    results.(roitype).(parms.measlist{m}).nvals = init_vals;
    results.(roitype).(parms.measlist{m}).nvals_valid = init_vals;
    results.(roitype).(parms.measlist{m}).nvals_invalid = init_vals;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_all_results(parms,results)
  results.all = [];
  results.all.nrois = 0;
  results.all.roicodes = [];
  results.all.roinames = {};
  results.all.data  = [];
  results.all.std   = [];
  results.all.nvals = [];
  results.all.nvals_valid = [];
  results.all.nvals_invalid = [];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_roi_data(parms,results,roitype)
  results = init_roi_results(parms,results,roitype);
  if strcmp(roitype,'cortsurf')
    pvec = [1:parms.nprojdist];
  else
    pvec = 1;
  end;
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    sf = parms.scalefacts(m);
    if parms.verbose==2
      fprintf('%s: compiling %s %s results for %d subjects...\n',...
        mfilename,meas,roitype,parms.nsubs);
    end;
    for s=1:parms.nsubs
      for p=pvec
        roi_data = load_roi_data(parms,s,meas,roitype,p);
        if isempty(roi_data), continue; end;
        [tmp,ind_data,ind_sel] = ...
          intersect({roi_data.roiname},parms.([roitype '_roinames']));
%        [tmp,ind_data,ind_sel] = ...
%          intersect([roi_data.roicode],parms.([roitype '_roicodes']));
        if isempty(ind_data)
          if parms.verbose
            fprintf('%s: WARNING: no valid %s ROIs for %s\n',...
              mfilename,roitype,parms.StudyInfo(s).VisitID);
          end;
          continue;
        end;
        results.(roitype).(meas).data(s,ind_sel,p)  = [roi_data(ind_data).avg]*sf;
        results.(roitype).(meas).std(s,ind_sel,p)   = [roi_data(ind_data).stdv]*sf;
        results.(roitype).(meas).nvals(s,ind_sel,p)  = [roi_data(ind_data).nvals];
        results.(roitype).(meas).nvals_valid(s,ind_sel,p) = ...
          [roi_data(ind_data).nvals_valid];
        results.(roitype).(meas).nvals_invalid(s,ind_sel,p) = ...
          [roi_data(ind_data).nvals_invalid];
      end;
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_contrast(parms,results)
  % calculate contrast ratio
  results.contrast = [];
  results.contrast.nrois = results.cortsurf.nrois;
  results.contrast.roicodes = results.cortsurf.roicodes;
  results.contrast.roinames = results.cortsurf.roinames;
  for m=1:parms.nmeas
    if parms.gwnorm_flag
      % normalize by mean
      mean_vals = mean(results.cortsurf.(parms.measlist{m}).data,3);
    else
      mean_vals = ones(size(diff_vals));
    end;
    % calculate contrast ratio
    results.contrast.(parms.measlist{m}).data = ...
      -diff(results.cortsurf.(parms.measlist{m}).data,[],3) ./ mean_vals;
    results.contrast.(parms.measlist{m}).std = ...
      sqrt(sum(results.cortsurf.(parms.measlist{m}).std.^2,3)) ./ mean_vals;
    results.contrast.(parms.measlist{m}).nvals = ...
      mean(results.cortsurf.(parms.measlist{m}).nvals,3); % should be identical
    results.contrast.(parms.measlist{m}).nvals_valid = ...
      min(results.cortsurf.(parms.measlist{m}).nvals_valid,[],3);
    results.contrast.(parms.measlist{m}).nvals_invalid = ...
      max(results.cortsurf.(parms.measlist{m}).nvals_invalid,[],3);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function results = compile_cort_area(parms,results)
  parms.thick_flag = 0;
  parms.area_flag = 1;
  parms.cortvol_flag = 0;
  parms.cortT1w_flag = 0;
  parms.subvol_flag = 0;
  parms.subT1w_flag = 0;
  parms.aparc_flag = 1;
  parms.fuzzy_flag = 0;
  parms.concat_flag = 0;
  parms.subjinfo_flag = 0;
  args = mmil_parms2args(parms,parms.compile_area_tags);
  tmp_results = MMIL_Compile_MRI_Analysis(parms.ProjID,args{:});
  % handle potential for subjects wtih missing MRI data (e.g. incomplete recon)
  DTI_VisitIDs = results.STRUCT_VisitIDs;
  MRI_VisitIDs = tmp_results.VisitIDs;
  [tmp,ind_MRI,ind_DTI] = intersect(MRI_VisitIDs,DTI_VisitIDs);
  results.cort_area = tmp_results.cort_area;
  results.cort_area.data = nan(results.nsubs,results.cort_area.nrois);
  results.cort_area.data(ind_DTI,:) = tmp_results.cort_area.data(ind_MRI,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_mean_cort_vals(parms,results,contrast_flag)
  if ~exist('contrast_flag','var') || isempty(contrast_flag)
    contrast_flag = 0;
  end;
  % results struct field name
  if contrast_flag
    tag = 'contrast';
  else
    tag = 'cortsurf';
  end;  
  % get roicodes from results
  roicodes = results.(tag).roicodes;
  roinames = results.(tag).roinames;
  codes_all = parms.aparc_roilist;
  codes_lh = codes_all(codes_all<2000);
  codes_rh = codes_all(codes_all>=2000);
  % identify ROIs from aparc
  [tmp,ind_lh]=intersect(roicodes,codes_lh);
  [tmp,ind_rh]=intersect(roicodes,codes_rh);
  [tmp,ind_all]=intersect(roicodes,codes_all);
  % calculate mean values for each measure
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    % prepare cortical area values for calculating weighted averages
    if parms.area_flag
      np = size(results.(tag).(meas).data,3);
      % get area data
      area_lh = repmat(results.cort_area.data(:,ind_lh),[1,1,np]);
      area_rh = repmat(results.cort_area.data(:,ind_rh),[1,1,np]);
      area_all = repmat(results.cort_area.data(:,ind_all),[1,1,np]);
    else
      area_lh = ones(size(results.(tag).(meas).data(:,ind_lh,:)));
      area_rh = ones(size(results.(tag).(meas).data(:,ind_rh,:)));
      area_all = ones(size(results.(tag).(meas).data(:,ind_all,:)));
    end;
    % calculate mean of ROI measures
    statlist = {'data','std'};
    for s=1:length(statlist)
      stat = statlist{s};
      % get data
      tmp_lh = results.(tag).(meas).(stat)(:,ind_lh,:);
      tmp_rh = results.(tag).(meas).(stat)(:,ind_rh,:);
      tmp_all = results.(tag).(meas).(stat)(:,ind_all,:);
      % calculate mean value weighted by area
      tmp_lh = mmil_wtd_mean(tmp_lh,area_lh,2);
      tmp_rh = mmil_wtd_mean(tmp_rh,area_rh,2);
      tmp_all = mmil_wtd_mean(tmp_all,area_all,2);
      % concatenate to matrix
      results.(tag).(meas).(stat) = ...
        cat(2,results.(tag).(meas).(stat),tmp_lh,tmp_rh,tmp_all);
    end;
    % calculate sum of nvals
    statlist = {'nvals','nvals_valid','nvals_invalid'};
    for s=1:length(statlist)
      stat = statlist{s};
      % get data
      tmp_lh = results.(tag).(meas).(stat)(:,ind_lh,:);
      tmp_rh = results.(tag).(meas).(stat)(:,ind_rh,:);
      tmp_all = results.(tag).(meas).(stat)(:,ind_all,:);
      % calculate sum
      tmp_lh = sum(tmp_lh,2);
      tmp_rh = sum(tmp_rh,2);
      tmp_all = sum(tmp_all,2);
      % concatenate to matrix
      results.(tag).(meas).(stat) = ...
        cat(2,results.(tag).(meas).(stat),tmp_lh,tmp_rh,tmp_all);
    end;
  end;
  % set roicodes
  new_roicodes = [1040;2040;2050];
  new_roinames = {'ctx-lh-mean';'ctx-rh-mean';'ctx-mean'};
  % concatenate to roicodes, and roinames
  results.(tag).roicodes = cat(1,results.(tag).roicodes,new_roicodes);
  results.(tag).roinames = cat(1,results.(tag).roinames,new_roinames);
  results.(tag).nrois = length(results.(tag).roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = concat_results(parms,results)
  results = init_all_results(parms,results);
  roitypes = {};
  if parms.aseg_flag,     roitypes{end+1} = 'aseg';     end;
  if parms.cortsurf_flag, roitypes{end+1} = 'cortsurf'; end;
  if parms.contrast_flag,   roitypes{end+1} = 'contrast'; end;
  for t=1:length(roitypes)
    roitype = roitypes{t};
    roicodes = results.(roitype).roicodes;
    if strcmp(roitype,'cortsurf')
      pvec = [1:parms.nprojdist];
    else
      pvec = 1;
    end;
    measlist = parms.measlist;
    nmeas = length(measlist);
    for p=pvec
      for m=1:nmeas
        meas = measlist{m};
        results.all.roicodes = cat(1,results.all.roicodes,roicodes);
        roinames = results.(roitype).roinames;
        for i=1:length(roinames)
          stem = [roitype '_' meas];
          if strcmp(roitype,'cortsurf')
            projdist = parms.projdist_list(p);
            if projdist < 0
              pstr = sprintf('white%0.1f',projdist);
            else
              pstr = sprintf('gray+%0.1f',projdist);
            end;
            stem = [stem '_' pstr];
          end;
          roinames{i} = [stem '-' roinames{i}];
        end;
        results.all.roinames = cat(1,results.all.roinames,roinames);
        results.all.data = ...
          cat(2,results.all.data,results.(roitype).(meas).data(:,:,p));
        results.all.std = ...
          cat(2,results.all.std,results.(roitype).(meas).std(:,:,p));
        results.all.nvals = ...
          cat(2,results.all.nvals,results.(roitype).(meas).nvals(:,:,p));
        results.all.nvals_valid = ...
          cat(2,results.all.nvals_valid,results.(roitype).(meas).nvals_valid(:,:,p));
        results.all.nvals_invalid = ...
          cat(2,results.all.nvals_invalid,results.(roitype).(meas).nvals_invalid(:,:,p));
      end;
    end;
  end;
  results.all.nrois = length(results.all.roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = load_roi_data(parms,s,meas,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  roi_data = [];
  % set input file name
  indir = parms.StudyInfo(s).proc_asl;
  if isempty(indir), return; end;
  indir = [parms.RootDirs.proc_asl '/' indir '/' parms.analdir];
  if ~exist(indir,'dir')
    if parms.verbose
      fprintf('%s: WARNING: input dir %s not found for %s\n',...
        mfilename,indir,parms.StudyInfo(s).VisitID);
    end;
    return;
  end;
  infix = set_infix(parms,roitype,p);
  fname_in = sprintf('%s/%s_%s.mat',indir,meas,infix);
  if ~exist(fname_in,'file')
    if parms.verbose
      fprintf('%s: WARNING: missing %s for %s\n',...
        mfilename,fname_in,parms.StudyInfo(s).VisitID);
    end;
    return;
  end;
  % load ROI data
  if parms.verbose==2
    fprintf('%s: loading %s results for %s...\n',...
      mfilename,meas,parms.StudyInfo(s).VisitID);
  end;
  load(fname_in);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

