function outfiles=DTI_MMIL_Compile_ROIs(RootDirs,varargin)
%function outfiles=DTI_MMIL_Compile_ROIs(RootDirs,[options])
%
% Purpose: Compile DTI ROI data across subjects
%   and save struct array in one matfile for each measure
%
% Usage:
%  outfiles=DTI_MMIL_Compile_ROIs(RootDirs,'key1', value1,...);
%
% Required Parameters:
%  RootDirs: struct that must contain the following fields:
%       proc_dti
%    these specify the full paths of root directories containing data containers
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with EPDTI_MMIL_Read_StudyInfo)
%     If empty, use all subjects in RootDirs.proc_dti
%     {default = []}
%
% Optional Parameters specifying diffusion data used to calculate ROI averages:
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'gruw_reg_mc_ecc_B0uw_iso'
%     {default = []}
%   'snums_flag': [0|1|2|3] which set of "snums" to use
%      0: use all available scans
%      1: use scan numbers in StudyInfo.DTIScanNums
%      2: use scan numbers in StudyInfo.DTIScanNums2
%      3: use scan numbers in StudyInfo.DTIScanNums3
%     {default = 0}
%   'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2) to use in spread sheet
%     (must have run DT calculations separately for each scan)
%     If empty, use DT measures calculated from all DTIScanNums
%     {default = []}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     (DT measures only)
%     {default = 0}
%   'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%
% Optional Parameters that specify how fiber ROI averages were calculated:
%   'measlist': cell array of DT-derived measures
%     {default = {'FA','ADC'}}
%   'atlas_flag': whether to use atlas fibers and if so,
%      what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'weighted_avg_flag': [0|1] whether weighted averages were calculated for
%      each ROI using fiber counts or fiber probabilities to weight
%      contribution from each voxel
%     {default = 1}
%   'thresh_FA': whether fiber ROIs was masked by applying threshold to FA image
%     {default = 0}
%   'thresh_prob': whether atlas-derived fiber ROIs were masked by applying
%      probability threshold
%     {default = 0}
%   'xcg_flag': [0|1] CSF and gray-matter voxels excluded
%     using FreeSurfer-derived aseg segmentation
%     {default = 1}
%   'masksf_flag': [0|1] exclude voxels with multiple fibers
%     {default = 0}
%  'fiber_atlasname': Name of atlas used in AtlasTrack {default=[]} Uses the default atlas
%
% Optional Parameters (misc)
%   'outdir': output directory
%     {default = pwd}
%   'forceflag': overwrite existing output files
%     {default = 0}
%
% Output:
%   outfiles: cell array of output mat file names
%
% Created:  11/07/08 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%% todo: add aseg and aparc ROIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULT_MIN_BVAL = 1;
DEFAULT_MIN_NDIRS = 6;

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'infix',[],[],...
  'snums_flag',0,[0:3],...
  'snum_index',[],[],...
  'revflag',0,{0,1,2},...
  'nob0_flag',false,[false true],...
  'min_bval',DEFAULT_MIN_BVAL,[],...
  'flex_flag',false,[false true],...
  'atlas_flag',2,[0:4],...
  'measlist',{'FA','ADC'},[],...
  'weighted_avg_flag',true,[false true],...
  'thresh_FA',0,[0,1],...
  'thresh_prob',0,[0,1],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
  'outdir',pwd,[],...
  'forceflag',false,[false true],...
... % hidden parameters:
  'min_ndirs',DEFAULT_MIN_NDIRS,[],...
  'min_nb0',1,[],...
  'xcg_suffix','xcg',[],...
  'masksf_suffix','masksf',[],...
  'fiber_atlasname',[],[],...
});

outfiles = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure RootDirs has required fields 'proc_dti'
RootDirs = MMIL_Check_RootDirs(RootDirs,{'proc_dti'});

if isempty(parms.measlist)
  error('measlist is empty')
end;
if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set file name infixes (must match those set by DTI_MMIL_Analyze_Fibers_Exam)
switch parms.atlas_flag
  case 0 % manual
    instr = [];
  case 1 % loc only, count atlas
    instr = '_loc_countatlas';
  case 2 % loc+dir, count atlas
    instr = '_countatlas';
  case 3 % loc only, mask atlas
    instr = '_loc_atlas';
  case 4 % loc+dir, mask atlas
    instr = '_atlas';
end;    
if parms.weighted_avg_flag
  instr = [instr '_wtd'];
end;
if parms.thresh_prob>0 && parms.atlas_flag>0
  instr = sprintf('%s_pthresh%0.2f',instr,parms.thresh_prob);
end;
if parms.thresh_FA>0
  instr = sprintf('%s_FAthresh%0.2f',instr,parms.thresh_FA);
end;
if parms.xcg_flag
  instr = [instr '_' parms.xcg_suffix];
end;
if parms.masksf_flag
  instr = [instr '_' parms.masksf_suffix];
end;
if ~isempty(parms.fiber_atlasname)
  instr = [instr '_' parms.fiber_atlasname];
end
outstr = ['fibers' instr];
instr = ['fiber_data' instr];
if parms.snums_flag
  outstr = sprintf('%s_snums%d',outstr,parms.snums_flag);
end;
if ~isempty(parms.snum_index)
  outstr = sprintf('%s_snumidx%d',outstr,parms.snum_index);
end;
outstr = outstr;
if parms.nob0_flag
  outstr = [outstr '_nob0'];
end;
if parms.min_bval ~= DEFAULT_MIN_BVAL
  outstr = sprintf('%s_minb%d',outstr,parms.min_bval);
end
if parms.flex_flag
  outstr = sprintf('%s_flex',outstr);
end;
if parms.min_ndirs ~= DEFAULT_MIN_NDIRS
  outstr = sprintf('%s_mind%d',outstr,parms.min_ndirs);
end
outstr = [outstr '_groupdata'];

% create output directory
if isempty(parms.outdir), parms.outdir = pwd; end;
mmil_mkdir(parms.outdir);

% generate StudyInfo struct if none supplied
[StudyInfo,SIflags] = DTI_MMIL_Check_StudyInfo(...
  parms.StudyInfo,RootDirs,...
  'snums_flag',parms.snums_flag);
nsubs = length(StudyInfo);
if nsubs==0
  error('no valid subjects found for summary');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=1:length(parms.measlist)
  nrois = [];
  meas = parms.measlist{m};
  outfiles{m} = sprintf('%s/%s',parms.outdir,meas);
  if ~exist(outfiles{m},'file') || parms.forceflag
    fprintf('%s: compiling %s results for %d subjects...\n',...
      mfilename,meas,nsubs);
    clear group_data; d = 1;
    for i = 1:length(StudyInfo)
      group_data(d).SubjID = StudyInfo(i).SubjID;
      group_data(d).StudyDate = num2str(StudyInfo(i).StudyDate);
      snums = StudyInfo(i).DTIScanNums;
      tmp_parms = parms;
      tmp_parms.snums = StudyInfo(i).DTIScanNums;
      if SIflags.group
        group_data(d).Group = StudyInfo(i).Group;
      else
        group_data(d).Group = 'all';
      end;
      if SIflags.visitcode
        group_data(d).VisitNumber = StudyInfo(i).VisitNumber;
      else
        group_data(d).VisitNumber = 'none';
      end;
      if SIflags.age
        group_data(d).Age = StudyInfo(i).Age;
      else
        group_data(d).Age = 0;
      end;
      if SIflags.stats
        group_data(d).StatsFlag = StudyInfo(i).StatsFlag;
      else
        group_data(d).StatsFlag = 1;
      end;
      ContainerPath = sprintf('%s/%s',RootDirs.proc_dti,StudyInfo(i).proc_dti);
      % load ROI results
      fprintf('%s: Loading %s results for %s %s\n',...
        mfilename,meas,group_data(d).SubjID,group_data(d).StudyDate);
      tags = {'snums','infix','revflag','min_bval','flex_flag',...
              'min_ndirs','min_nb0','nob0_flag'};
      args = mmil_parms2args(tmp_parms,tags);
      fstem = DTI_MMIL_Set_DT_fstem(ContainerPath,args{:});
      [tmp_path,tmp_stem,tmp_ext] = fileparts(fstem);
      %% todo: make DTanalysis a parameter
      out_fstem = sprintf('%s/DTanalysis/%s%s',ContainerPath,tmp_stem,tmp_ext);
      fname_in_mat = sprintf('%s_%s_%s.mat',out_fstem,meas,instr);
      if ~exist(fname_in_mat,'file')
        fprintf('%s: WARNING: missing %s for %s %s\n',...
          mfilename,fname_in_mat,group_data(d).SubjID,group_data(d).StudyDate);
        continue;
      end;
      fiber_data = [];
      load(fname_in_mat);
      if isempty(fiber_data), continue; end;
      tmp_nrois = length(fiber_data);
      if isempty(nrois)
        nrois = tmp_nrois;
      elseif tmp_nrois~=nrois
        % make sure nrois is same for each subject
        error('number of fiber ROIs for %s (%d) does not match first subject (%d)',...
          group_data(d).SubjID,tmp_nrois,nrois);
      end;
      for f=1:nrois
        group_data(d).fiber_data(f).avg = fiber_data(f).avg;
        group_data(d).fiber_data(f).stdv = fiber_data(f).stdv;
        group_data(d).fiber_data(f).nvox = fiber_data(f).nvox;
        group_data(d).fiber_data(f).nvals = fiber_data(f).nvals;
        group_data(d).fiber_data(f).fibernum = fiber_data(f).fibernum;
      end;
      d=d+1;
    end
    fprintf('%s: saving group data to %s...\n',mfilename,outfiles{m});
    save(outfiles{m},'group_data');
  end;
end;

