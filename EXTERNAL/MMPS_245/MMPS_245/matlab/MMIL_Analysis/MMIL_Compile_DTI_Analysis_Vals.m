function results = MMIL_Compile_DTI_Analysis_Vals(ProjID,varargin)
%function results = MMIL_Compile_DTI_Analysis_Vals(ProjID,[options])
%
% Usage:
%  results = MMIL_Compile_DTI_Analysis_Vals(ProjID,'key1', value1,...);
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
%      proc_dti, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%
% Optional Parameters that specify which analysis results to compile:
%   'measlist': list of DT "measures" to extract for each fiber
%      e.g. 'FA', 'MD', 'LD', 'TD', 'b0', 'b0N', 'T2w', 'T1w'
%      'b0' is the b=0 volume collected with diffusion scan
%      'b0N' is b=0 volume normalized by mean signal in ventricles
%      'T2w' is b=0 volume normalized by linear fit to MD
%      'T1w' is the T1-weighted FreeSurfer nu.mgz volume
%       use 'none' to extract none of these measures
%     {default = {'MD'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     if empty, will use 1 for each
%     if not empty, length must match length of measlist
%     {default = []}
%   'inputlist': list of input files in addition to or in place of measlist
%     relative to each ContainerPath
%     {default = []}
%   'fiber_flag': [0|1] compile white matter fiber ROI results
%     {default = 1}
%   'aseg_flag': [0|1] compile aseg ROI results
%     {default = 0}
%   'wmparc_flag': [0|1] extract wmparc ROI results
%     {default = 0}
%   'cortsurf_flag': [0|1] compile cortical surface ROI results
%     {default = 0}
%
% Optional Parameters that determine input diffusion data:
%  'snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%    {default = 1}
%  'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2) to use in spread sheet (must have run DT
%     calculations separately for each scan)
%     If empty, use DT measures calculated from all DTIScanNums
%     {default = []}
%  'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix: 'corr_resDTI'
%     {default = []}
%  'auto_infix_flag': [0|1] set infix automatically based on typical
%     processing and settings in ProjInfo
%     ignored if infix is not empty
%     {default = 1}
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%  'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'min_ndirs': minimum number of diffusion directions to be included
%     {default = 6}
%  'full_fstem_flag': [0|1] full stem included in DT measures analysis output
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     if 0, options above (snums_flag, snum_index, infix, etc.) are ignored
%     {default = 0}
%
% Optional Parameters specific to fiber analysis:
%   'fiber_indir': name of DTI fiber analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%   'fibers': vector of fiber numbers
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so,
%      what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'weighted_avg_flag': [0|1] use weighted average fiber ROI measures
%     {default = 1}
%   'thresh_FA': FA threshold used to calculate fiber ROI measures
%     {default = 0}
%   'thresh_prob': fiber probability threshold
%     {default = 0}
%   'fiber_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'fiber_dispfact': multiplicative factor applied to fiber dispersion values
%     {default = 4}
%   'xcg_flag': [0|1] exclude CSF and gray-mattter from fiber ROIs
%     {default = 0}
%   'masksf_flag': [0|1] voxels with multiple fibers excluded
%     {default = 0}
%   'fiber_atlasname': name of the atlas used in AtlasTrack
%     and attached to output analysis mat files (empty = default atlas)
%     {default = []}
%   'fname_fiber_legend': name of fiber legend csv file
%     if empty, use default ($MMPS_PARMS/DTI_Fiber/DTI_Fiber_Legend.csv)
%     {default = []}
%
% Optional Parameters specific to aseg analysis:
%   'aseg_indir': name of DTI aseg analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%   'aseg_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'aseg_dispfact': multiplicative factor applied to aseg dispersion values
%     {default = 4}
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
% Optional Parameters specific to wmparc ROIs:
%   'wmparc_aparc_flag': [0|1] whether to use cortical parcellation ROIs
%     0: wmparc only
%     1: wmparc and aparc
%     {default = 0}
%
% Optional Parameters specific to cortical surface analysis:
%   'cortsurf_indir': name of DTI cortsurf analysis folder
%     relative to ContainerPath
%     {default = 'DTanalysis'}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     {default = [-1,1]}
%
% Other Optional Parameters:
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%
% Output:
%   results: struct containing compiled DTI analysis results with fields:
%     StudyInfo: struct array containing info for each subject
%     RootDirs: struct containing paths of root directories
%     SubjIDs: cell array of subject IDs
%     VisitIDs: cell array of visit IDs
%     STRUCT_VisitIDs: cell array of structural (i.e. FreeSurfer) visit IDs
%     VisitNumbers: vector of visit numbers
%     measlist: cell array of measure names
%     nsubs: number of subjects/visits in StudyInfo
%     nmeas: number of measures in measlist
%     nrois: number of ROIs compiled
%     StudyData: struct array containing ROI data for each subject
%       roi_data: struct array containing data for each ROI
%         roicode: numeric code of ROI
%         roiname: name of ROI
%         nvals: number of voxels or vertices in ROI
%         nvals_valid: number of valid values in ROI
%         nvals_invalid: number of valid values in ROI
%         nmeas: number of measures in measlist
%         vals: matrix of values with size = [nvals,nmeas]
%         weights: matrix of values with size = [nvals,nmeas] (optional)
%
% Created:  10/21/12 by Don Hagler
% Last Mod: 07/31/15 by Don Hagler
%

% created based on MMIL_Compile_DTI_Analysis, created 10/06/12 by Don Hagler

%% todo: wmparc_disp_flag, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

[results,parms] = init_results(parms);

if parms.fiber_flag
  [results,parms] = compile_roi_vals(parms,results,'fiber');
end;

if parms.aseg_flag
  [results,parms] = compile_roi_vals(parms,results,'aseg');
end;

if parms.wmparc_flag
  [results,parms] = compile_roi_vals(parms,results,'wmparc');
end;

if parms.cortsurf_flag
  [results,parms] = compile_roi_vals(parms,results,'cortsurf');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
... % specify which results to compile
    'measlist',{'MD'},[],...
    'scalefacts',[],[],...
    'inputlist',[],[],...
    'fiber_flag',true,[false true],...
    'aseg_flag',false,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',false,[false true],...
... % specify DTI data used for tensor calculations and analysis
    'snums_flag',1,[0:3],...
    'snum_index',[],[],...
    'infix',[],[],...
    'auto_infix_flag',true,[false true],...
    'revflag',0,[0:2],...
    'nob0_flag',false,[false true],...
    'min_bval',1,[],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'min_ndirs',6,[],...
    'full_fstem_flag',false,[false true],...
... % fiber analysis
    'fiber_indir','DTanalysis',[],...
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'atlas_flag',2,[0:4],...
    'weighted_avg_flag',true,[false true],...
    'thresh_FA',0,[],...
    'thresh_prob',0,[],...
    'fiber_disp_flag',false,[false true],...
    'fiber_dispfact',4,[1e-6,1e6],...
    'xcg_flag',false,[false true],...
    'masksf_flag',false,[false true],...
    'fiber_atlasname',[],[],...
    'fname_fiber_legend',[],[],...
... % aseg analysis
    'aseg_indir','DTanalysis',[],...
    'aseg_disp_flag',false,[false true],...
    'aseg_dispfact',4,[1e-6,1e6],...
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
... % wmparc parameters
    'wmparc_aparc_flag',false,[false true],...    
... % cortsurf analysis
    'cortsurf_indir','DTanalysis',[],...
    'projdist_list',[-1,1],[-5,5],...
...
    'verbose',1,[0:2],...
... % hidden
    'qcflag',true,[false true],...
    'required_rootdirs',{'proc_dti','fsurf'},[],...
    'QC_raw',true,[false true],...
    'QC_recon',true,[false true],...
    'QC_DTI',true,[false true],...
    'resample_flag',true,[false true],...
    'regT1flag',1,[0:2],...
    'aseg_disp_suffix','dwtd',[],...
    'fiber_disp_suffix','dwtd',[],...
    'xcg_suffix','xcg',[],...
    'masksf_suffix','masksf',[],...
    'fname_fscolorlut',[],[],...
    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[],...
    'wmparc_roilist',[3001:3003,3005:3034,4001:4003,4005:4034],[1,Inf],...
    'DT_outdir','DTcalc',[],...
    'DT_analdir','DTanalysis',[],...
    'DT_outfix',[],[],...
    'RSI_outdir','RSIcalc',[],...
    'RSI_analdir','RSIanalysis',[],...
    'RSI_outfix',[],[],...
    'DT_measlist',{'FA','MD','LD','TD','b0','b0N','T2w'},[],...
    'RSI_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                    'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
...
    'ProjInfo_tags',{'min_bval','flex_flag',...
                     'xcg_flag','masksf_flag','fseg_flag',...
                     'cortsurf_flag','atlas_flag','revflag','snums_flag',...
                     'fiber_atlasname','resample_flag','regT1flag'},[],...
    'info_tags',{'snums','revflag','min_nb0','min_ndirs',...
                 'min_bval','flex_flag'},[],...
    'fstem_tags',{'snums','infix','revflag','min_bval','flex_flag',...
                  'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
    'aseg_code_tags',{'aseg_aparc_flag','aseg_roilist','aparc_roilist',...
                      'exclude_roilist'},[],...
    'wmparc_code_tags',{'wmparc_aparc_flag','wmparc_roilist','aparc_roilist',...
                        'exclude_roilist'},[],...
    'suffix_tags',{'xcg_flag','xcg_suffix',...
                   'masksf_flag','masksf_suffix',...
                   'disp_flag','disp_suffix','dispfact',...
                   'thresh_prob','thresh_FA'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(StudyInfo), error('empty StudyInfo'); end;

  % for arg names present in both varargin and ProjInfo,
  %   the options values will appear in merged_args
  %   i.e. command-line takes precedence, then ProjInfo, then defaults
  if ~isempty(ProjInfo)
    % get all 'DTI_' parms, strip 'DTI_'
    ProjInfo_args = MMIL_Args(ProjInfo,'DTI');
    % convert from args cell array to parms struct
    ProjInfo_parms = mmil_args2parms(ProjInfo_args,[],0);
    % keep only those in ProjInfo_tags
    ProjInfo_args = mmil_parms2args(ProjInfo_parms,parms.ProjInfo_tags);
    % merge arguments, giving varargin precedence over ProjInfo_args
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
  end;

  parms.StudyInfo = StudyInfo;
  parms.RootDirs = RootDirs;
  
  % set infix
  if parms.auto_infix_flag && isempty(parms.infix)
    parms.infix = 'corr';
    if parms.resample_flag
      if parms.regT1flag==2
        parms.infix = [parms.infix '_regT1'];
      else
        parms.infix = [parms.infix '_resDTI'];
      end;
    end;
  end;

  % set DTI-specific fields in StudyInfo
  [parms.StudyInfo,parms.SIflags] = DTI_MMIL_Check_StudyInfo(...
    parms.StudyInfo,parms.RootDirs,...
    'snums_flag',parms.snums_flag,... % DTIScanNums field will get changed
    'qcflag',parms.qcflag,...
    'DTI_flag',1);
  if isempty(parms.StudyInfo), error('empty StudyInfo after check for DTI'); end;
  parms.nsubs = length(parms.StudyInfo);

  if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;
  parms.nmeas = length(parms.measlist);
  if isempty(parms.scalefacts)
    parms.scalefacts = ones(1,parms.nmeas);
  elseif numel(parms.scalefacts)~=parms.nmeas
    error('number of scalefacts (%d) does not match measlist (%d)',...
      numel(parms.scalefacts),parms.nmeas);
  end;
  parms = combine_measlist(parms);

  % set parameters for fiber analysis
  if parms.fiber_flag
    FiberLegend = DTI_MMIL_Read_Fiber_Legend(parms.fname_fiber_legend);
    legend_fibernums = cell2mat({FiberLegend.FiberNumber});
    legend_fibernames = {FiberLegend.FiberName};    
    [tmp,ind_fibers,ind_legend] = intersect(parms.fibers,legend_fibernums);
    parms.fibers = parms.fibers(ind_fibers);
    parms.fiber_roicodes = mmil_colvec(parms.fibers);
    parms.fiber_roinames = mmil_colvec(legend_fibernames(ind_legend));
    parms.fiber_nrois = length(parms.fibers);
    if parms.fiber_nrois==0
      error('no valid fiber numbers');
    end;      
  else
    parms.fiber_nrois = 0;
  end;

  if parms.aseg_flag || parms.wmparc_flag || parms.cortsurf_flag
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
    args = mmil_parms2args(parms,parms.aseg_code_tags);
    parms.aseg_roicodes = mmil_aseg_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aseg_roicodes);
    parms.aseg_roicodes = fs_roicodes(ind_fs);
    parms.aseg_roinames = fs_roinames(ind_fs);
    parms.aseg_nrois = length(parms.aseg_roicodes);
  else
    parms.aseg_nrois = 0;
  end;

  % set parameters for wmparc analysis
  if parms.wmparc_flag
    args = mmil_parms2args(parms,parms.wmparc_code_tags);
    parms.wmparc_roicodes = mmil_wmparc_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.wmparc_roicodes);
    parms.wmparc_roicodes = fs_roicodes(ind_fs);
    parms.wmparc_roinames = fs_roinames(ind_fs);
    parms.wmparc_nrois = length(parms.wmparc_roicodes);
  else
    parms.wmparc_nrois = 0;
  end;

  % set parameters for cortsurf analysis
  if parms.cortsurf_flag
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aparc_roilist);
    parms.cortsurf_roicodes = fs_roicodes(ind_fs);
    parms.cortsurf_roinames = fs_roinames(ind_fs);
    parms.cortsurf_nrois = length(parms.cortsurf_roicodes);
    parms.nprojdist = length(parms.projdist_list);
  else
    parms.cortsurf_nrois = 0;
  end;

  parms.nrois = parms.fiber_nrois + parms.aseg_nrois +...
                parms.wmparc_nrois + parms.cortsurf_nrois;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = combine_measlist(parms)
  for m=1:parms.nmeas
    meas = parms.measlist{m};
  end;
  m = parms.nmeas + 1;
  for i=1:length(parms.inputlist)
    [tmp,meas] = fileparts(parms.inputlist{i});
    parms.measlist{m} = meas;
    m = m + 1;
  end;
  parms.nmeas = length(parms.measlist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = set_fstem(parms,s,meas)
  % set input file stem
  fstem = [];
  ContainerPath = sprintf('%s/%s',...
    parms.RootDirs.proc_dti,parms.StudyInfo(s).proc_dti);
  parms.snums = parms.StudyInfo(s).DTIScanNums;
  args = mmil_parms2args(parms,parms.info_tags);
  [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
  if errcode || isempty(ScanInfo)
    fprintf('%s: WARNING: no DTI scans for %s\n',mfilename,ContainerPath);
    return;
  end;
  switch meas
    case parms.DT_measlist
      if parms.full_fstem_flag
        parms.outdir = parms.DT_outdir;
        parms.outfix = parms.DT_outfix;
        args = mmil_parms2args(parms,parms.fstem_tags);
        fstem = DTI_MMIL_Set_DT_fstem(ContainerPath,args{:});
        if isempty(fstem), return; end;
        [tmp_path,fstem,tmp_ext] = fileparts(fstem);
      else
        fstem = parms.DT_outfix;
        if parms.flex_flag
          if isempty(fstem)
            fstem = 'flex';
          else
            fstem = ['flex_' fstem];
          end;
        end;
      end;
      indir = parms.DT_analdir;
    case parms.RSI_measlist
      if parms.full_fstem_flag
        parms.outdir = parms.RSI_outdir;
        parms.outfix = parms.RSI_outfix;
        args = mmil_parms2args(parms,parms.fstem_tags);
        fstem = DTI_MMIL_Set_RSI_fstem(ContainerPath,args{:});
        if isempty(fstem), return; end;
        [tmp_path,fstem,tmp_ext] = fileparts(fstem);
      else
        fstem = parms.RSI_outfix;
        if parms.flex_flag
          if isempty(fstem)
            fstem = 'flex';
          else
            fstem = ['flex_' fstem];
          end;
        end;
      end;
      indir = parms.RSI_analdir;
    otherwise
      fstem = [];
      indir = parms.DT_analdir;
  end;
  if ~isempty(fstem), fstem = [fstem '_']; end;
  fstem = sprintf('%s/%s/%s',ContainerPath,indir,fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function infix = set_infix(parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  switch roitype
    case 'fiber'
      % set file name stems (must match those set by DTI_MMIL_Analyze_Fibers_Exam)
      infix = 'fiber';
      switch parms.atlas_flag
        case 1 % loc only, count atlas
          infix = [infix '_loc_countatlas'];
        case 2 % loc+dir, count atlas
          infix = [infix '_countatlas'];
        case 3 % loc only, mask atlas
          infix = [infix '_loc_atlas'];
        case 4 % loc+dir, mask atlas
          infix = [infix '_atlas'];
      end;
      if parms.weighted_avg_flag
        infix = [infix '_wtd'];
      end;
      parms.disp_flag = parms.fiber_disp_flag;
      parms.dispfact = parms.fiber_dispfact;
      parms.disp_suffix = parms.fiber_disp_suffix;
      if parms.atlas_flag==0
        parms.thresh_prob = 0;
      end;
      args = mmil_parms2args(parms,parms.suffix_tags);
      suffix = dti_set_fiber_suffix(args{:});
      if ~isempty(suffix)
        infix = [infix '_' suffix];
      end;
      if ~isempty(parms.fiber_atlasname)
        infix = [infix '_' parms.fiber_atlasname];
      end
      infix = [infix '_roi_data'];
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
      if parms.aseg_disp_flag
        infix = sprintf('%s_%s%0.1f',...
          infix,parms.aseg_disp_suffix,parms.aseg_dispfact);
      end;
      infix = [infix '_roi_data'];
    case 'wmparc'
      infix = 'wmparc_roi_data';
    case 'cortsurf'
      infix = sprintf('pdist%0.1f_roi_data',...
        parms.projdist_list(p));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results,parms] = init_results(parms)
  results = [];
  results.RootDirs = parms.RootDirs;
  results.StudyInfo = parms.StudyInfo;
  results.StudyData = [];
  for s=1:parms.nsubs
    results.StudyData(s).roi_data = [];
    for r=1:parms.nrois
      results.StudyData(s).roi_data(r).roicode = 0;
      results.StudyData(s).roi_data(r).roiname = [];
      results.StudyData(s).roi_data(r).roitype = [];
      results.StudyData(s).roi_data(r).nvals = 0;
      results.StudyData(s).roi_data(r).nvals_valid = 0;
      results.StudyData(s).roi_data(r).nvals_invalid = 0;
      results.StudyData(s).roi_data(r).nmeas = parms.nmeas;
      results.StudyData(s).roi_data(r).vals = [];
      results.StudyData(s).roi_data(r).weights = [];
    end;    
    r = 1;
    if parms.fiber_flag
      for i=1:parms.fiber_nrois
        results.StudyData(s).roi_data(r).roicode = parms.fiber_roicodes(i) + 1e4;
        results.StudyData(s).roi_data(r).roiname = parms.fiber_roinames{i};
        results.StudyData(s).roi_data(r).roitype = 'fiber';
        r = r + 1;
      end;
    end;
    if parms.aseg_flag
      for i=1:parms.aseg_nrois
        results.StudyData(s).roi_data(r).roicode = parms.aseg_roicodes(i);
        results.StudyData(s).roi_data(r).roiname = parms.aseg_roinames{i};
        results.StudyData(s).roi_data(r).roitype = 'aseg';
        r = r + 1;
      end;
    end;
    if parms.wmparc_flag
      for i=1:parms.wmparc_nrois
        results.StudyData(s).roi_data(r).roicode = parms.wmparc_roicodes(i);
        results.StudyData(s).roi_data(r).roiname = parms.wmparc_roinames{i};
        results.StudyData(s).roi_data(r).roitype = 'wmparc';
        r = r + 1;
      end;
    end;
    if parms.cortsurf_flag
      for i=1:parms.cortsurf_nrois
        results.StudyData(s).roi_data(r).roicode = parms.cortsurf_roicodes(i);
        results.StudyData(s).roi_data(r).roiname = parms.cortsurf_roinames{i};
        results.StudyData(s).roi_data(r).roitype = 'cortsurf';
        r = r + 1;
      end;
    end;
  end;
  results.roicodes = [results.StudyData(1).roi_data.roicode]';
  results.roinames = {results.StudyData(1).roi_data.roiname}';
  results.roitypes = {results.StudyData(1).roi_data.roitype}';
  results.measlist = parms.measlist;
  results.nsubs = parms.nsubs;
  results.nrois = parms.nrois;
  results.nmeas = parms.nmeas;
  parms.nrois_compiled = 0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results,parms] = compile_roi_vals(parms,results,roitype)
  if strcmp(roitype,'cortsurf')
    pvec = [1:parms.nprojdist];
  else
    pvec = 1;
  end;
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).VisitID;
    if parms.verbose==2
      fprintf('%s: compiling %s results for %s...\n',...
        mfilename,roitype,VisitID);
    end;
    for m=1:parms.nmeas
      meas = parms.measlist{m};
      sf = parms.scalefacts(m);
      for p=pvec
        roi_data = load_roi_data(parms,s,meas,roitype,p);
        if isempty(roi_data)
          continue;
        end;
        [tmp,ind_data,ind_sel] = ...
          intersect([roi_data.roicode],parms.([roitype '_roicodes']));
        if isempty(ind_data)
          if parms.verbose
            fprintf('%s: WARNING: no valid %s ROIs for %s\n',...
              mfilename,roitype,VisitID);
          end;
          continue;
        end;
        for i=1:length(ind_sel)
          j = ind_sel(i);
          k = ind_data(i);
          r = parms.nrois_compiled + j;
          nvals = roi_data(k).nvals;
          if m==1 && p==1
            results.StudyData(s).roi_data(r).nvals = nvals;
            results.StudyData(s).roi_data(r).nvals_valid = ...
              roi_data(k).nvals_valid;
            results.StudyData(s).roi_data(r).nvals_invalid = ...
              roi_data(k).nvals_invalid;
            if nvals
              results.StudyData(s).roi_data(r).vals = ...
                zeros(nvals,parms.nmeas,length(pvec));
            end;
            if isfield(roi_data,'weights')
              results.StudyData(s).roi_data(r).weights = roi_data(k).weights;
            end;
          elseif nvals ~= results.StudyData(s).roi_data(r).nvals
            error('number of values for meas %s does not match %s',...
              parms.measlist{m},parms.measlist{1});
          end;
          if nvals
            results.StudyData(s).roi_data(r).vals(:,m,p) = roi_data(k).vals*sf;
          end;
        end;
      end;
    end
  end;
  parms.nrois_compiled = parms.nrois_compiled + parms.([roitype '_nrois']);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = load_roi_data(parms,s,meas,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  roi_data = [];
  % set input file name
  fstem = set_fstem(parms,s,meas);
  if isempty(fstem), return; end;
  infix = set_infix(parms,roitype,p);
  fname_in = sprintf('%s%s_%s.mat',fstem,meas,infix);
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

