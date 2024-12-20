function MMIL_RetFit_BOLD_Exams(ProjID,varargin)
% function MMIL_RetFit_BOLD_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    including these fields: SubjID, VisitID, STRUCT_VisitID,
%    This option can be used to specify subject specific values for the
%     "RF" parameters below (e.g. RF_pol_dir, RF_pol_stem, RF_pol_snums, etc.)
%    Must have 'VisitID' to indicate name of the orig data directory
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: proc_bold, fsurf
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied, MMIL_ProjInfo.csv is not required
%    {default = []}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'batchname': name of ougtput batchdir
%    {default = 'MMIL_RetFit_BOLD_Exams'}
%  'verbose': output messages about VisitIDs skipped because of missing data
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0 }
%
% Optional Parameters for Input Data:
%  'multisess_flag': [0|1] indicates source of input analyses
%    0: use results from one or more scans within session
%    1: use results from multiple sessions
%    if 1, RF_pol_dir, RF_ecc_dir, etc. will be ignored
%    {default = 0}
%  'RF_pol_dir': polar angle BOLD analysis subdirectory of ContainerPath
%    if empty, will attempt to find suitable value
%    {default = []}
%  'RF_ecc_dir': eccentricity BOLD analysis subdirectory of ContainerPath
%    if empty, will attempt to find suitable value
%    {default = []}
%  'RF_pol_stem': polar angle file stem
%    if empty, will attempt to find suitable value
%    {default = []}
%  'RF_ecc_stem' eccentricity file stem
%    if empty, will attempt to find suitable value
%    {default = []}
%  'RF_pol_snums': polar angle BOLD scan numbers
%    if empty, will be set to odd scans
%      ('for' or 'rev' depending on SessInfo.revflag)
%    {default = []}
%  'RF_ecc_snums': eccentricity BOLD scan numbers
%    if empty, will be set to even scans
%      ('for' or 'rev' depending on SessInfo.revflag)
%    {default = []}
%  'RF_BOLD_infix': string inside BOLD file names (e.g. 'corr_resBOLD')
%    {default = 'corr_resBOLD'}
%
% Optional Parameters for Retinotopy Fitting:
%  'RF_outdir': output directory
%    full path or relative to ContainerPath
%    {default = 'retfit'}
%  'RF_outstem': output file stem
%    {default = 'retfit'}
%  'RF_roi_name': file stem of ROI file (e.g. lh.roi_name.label)
%    {default = 'v123'}
%  'RF_r_max': maximum eccentricity (for stimulus presentation)
%    {default = 15}
%  'RF_multisess_r_max': if multisess_flag=1, will use this value
%     instead of RF_r_max (if not empty)
%    {default = []}
%  'RF_r_min': minimum eccentricity (for stimulus presentation)
%    If empty, will use RF_r_min_factor*RF_r_max
%    {default = []}
%  'RF_r_min_factor': multiplier of RF_r_max to obtain RF_r_min, if empty
%    {default = 0.02}
%  'RF_logtrans_flag': [0|1] whether log transform was used for ecc stimulus
%    {default = 0}
%  'RF_map_v123_flag': [0|1] whether to model V1-V2-V3 complex or
%     a single area mapping entire hemifield
%     {default = 1}
%  'RF_map_poly_flag': [0|1] use polynomial function to deform template
%    {default = 1}
%  'RF_map_poly_order': order of polynomial function (n+1 additional parameters)
%     used to deform template
%    {default = 4}
%  'RF_map_model_type': [0|1|2] model used for initial estimates of u and v
%    0: rectangle
%    1: wedge
%    2: radial wedge
%    {default = 2}
%  'RF_map_area_name': area label if map_v123_flag=0
%     {default = 'v'}
%  'RF_map_rev_polar_flag': [0|1] whether to reverse direction of polar angle
%     {default = 0}
%
% Created:  02/21/11 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

%% todo: get FP_snums and FE_snums to override RF_pol_snums and RF_ecc_snums?
%%       use FP_r_max or FE_r_max instead of RF_r_max?

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'qcflag',true,[false true],...
  'batchname','MMIL_RetFit_BOLD_Exams',[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
...
  'multisess_flag',false,[false true],...
  'RF_pol_dir',[],[],...
  'RF_ecc_dir',[],[],...
  'RF_pol_stem',[],[],...
  'RF_ecc_stem',[],[],...
  'RF_pol_snums',[],[],...
  'RF_ecc_snums',[],[],...
  'RF_BOLD_infix','corr_resBOLD',[],...
  'RF_Fourier_infix','fstats_pval',[],...
  'RF_outdir','retfit',[],...
  'RF_outstem','retfit',[],...
  'RF_roi_name','v123',[],...
  'RF_r_max',15,[1,100],...
  'RF_multisess_r_max',[],[1,100],...
  'RF_r_min',[],[0,100],...
  'RF_r_min_factor',0.02,[0,1],...
  'RF_logtrans_flag',false,[false true],...
  'RF_fnamestem','BOLD',[],...
... % retinotopy data / tksurfer
  'RF_hemilist',{'lh','rh'},{'lh','rh'},...
  'RF_suffixlist',{'_r','_i'},{'_r','_i'},...
  'RF_outstemlist',{'pol','ecc'},[],...
  'RF_surf','sphere',{'white','pial','inflated','sphere'},...
  'RF_smooth',10,[0,100],...
  'RF_fthresh',0,[],...
  'RF_fmid',1.5,[],...
  'RF_fslope',3,[],...
  'RF_revflag',false,[false true],...
  'RF_sph_rot',{[45 0 90],[45 -20 -90]},[],...
... % retfit parameters
  'RF_roi_dilate_niters',0,[0,1000],...
  'RF_roi_rotation',65,[-180,180],...
  'RF_roi_shift_u',0,[],...
  'RF_roi_shift_v',0,[],...
  'RF_roi_scale_u',0.7,[],...
  'RF_roi_scale_v',0.7,[],...
  'RF_prereg_nruns_quick',2,[],...
  'RF_prereg_niter_quick',100,[],...
  'RF_prereg_step_size_quick',[0.1,0.05,0.01],[],...
  'RF_prereg_nruns',200,[],...
  'RF_prereg_niter',100,[],...
  'RF_map_poly_coef_range',[-5,5],[-100,100],...
  'RF_map_radial_wedge_fact_range',[0.05 0.4],[],...
  'RF_map_radial_offset_range',[1 4],[],...
  'RF_map_scale_u_range',[0.6 0.8],[],...
  'RF_map_scale_v_range',[0.3 0.6],[],...
  'RF_map_rotation_range',[-10 10],[],...
  'RF_map_shift_u_range',[-0.1,0.1],[],...
  'RF_map_shift_v_range',[-0.1,0.1],[],...
  'RF_map_wedge_fact_range',[0.7 1.3],[],...
  'RF_map_r_min_range',[1 1],[],...
  'RF_map_r_max_range',[12 12],[],...
  'RF_map_v1_width_range',[1 1],[],...
  'RF_map_v2_width_range',[0.6 1],[],...
  'RF_map_v3_width_range',[0.5 1],[],...
  'RF_map_v1_length_range',[0.8 1.2],[],...
  'RF_map_v2_length_range',[0.8 1.2],[],...
  'RF_map_v3_length_range',[0.8 1.2],[],...
  'RF_nruns',1,[],...
  'RF_niter',2000,[],...
  'RF_ecc_fact',1,[0,1],...
  'RF_smooth_fact',1,[0,Inf],...
  'RF_fold_fact',15,[0,Inf],...
  'RF_vacancy_fact',0,[0,Inf],...
  'RF_max_outbound_penalty',20,[0 10000],...
  'RF_data_smooth_sigma',0.1,[0,1],...
  'RF_cost_include_percentile',100,[50,100],...
  'RF_map_v123_flag',true,[false true],...
  'RF_map_poly_flag',true,[false true],...
  'RF_map_poly_order',4,[1,10],...
  'RF_map_model_type',2,[0,1,2],...
  'RF_map_logtrans_flag',true,[false true],...
  'RF_map_area_name','v',[],...
  'RF_map_rev_polar_flag',false,[false true],...
...
  'required_containers',{'proc_bold','fsurf'},[],...
  'QC_raw',true,[false true],... % only applies if manual raw QC exists
  'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
  'QC_recon',true,[false true],...
};
%% todo: set numvec_tags to include relevant RF options (e.g. RF_map_v2_width_range)
%%       so that range options may be set from VisitInfo (but what about LH and RH?)

parms = mmil_args2parms(varargin,parms_filter);

% excl_tags are fields that should not be passed to MMIL_RetFit_BOLD_Exam
excl_tags = {'StudyInfo','RootDirs','qcflag','batchname',...
  'RF_multisess_r_max','required_containers','verbose',...
  'QC_raw','QC_BOLD','QC_recon'};
tags = setdiff(fieldnames(parms),excl_tags);

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid<0
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

% create jobs for each subject in StudyInfo
j = 1;
for i=1:length(StudyInfo)
  SubjID = StudyInfo(i).SubjID;
  VisitID = StudyInfo(i).VisitID;
  ContainerPath = sprintf('%s/%s',RootDirs.proc_bold,StudyInfo(i).proc_bold);
  FSContainerPath = sprintf('%s/%s',RootDirs.fsurf,StudyInfo(i).fsurf);

  % replace values in parms with (non-empty) values from StudyInfo
  tmp_parms = parms;
  for t=1:length(tags)
    tmp_val = mmil_getfield(StudyInfo(i),tags{t});
    if ~isempty(tmp_val), tmp_parms.(tags{t}) = tmp_val; end;
  end;
  % replace r_max for multisess
  if parms.multisess_flag
    tmp_val = mmil_getfield(StudyInfo(i),'RF_multisess_r_max',...
      parms.RF_multisess_r_max);
    if ~isempty(tmp_val), tmp_parms.RF_r_max = tmp_val; end;
  end;

  % check that we have BOLD analysis files (polar and eccen)
  % NOTE: this may take ~1 sec because of loading ContainerInfo
  [pol_stem,ecc_stem,errcode] = BOLD_MMIL_Check_RetFit_Data(ContainerPath,...
   'infix',tmp_parms.RF_BOLD_infix,'fstats_infix',tmp_parms.RF_Fourier_infix,...
   'pol_dir',tmp_parms.RF_pol_dir,'ecc_dir',tmp_parms.RF_ecc_dir,...
   'pol_stem',tmp_parms.RF_pol_stem,'ecc_stem',tmp_parms.RF_ecc_stem,...
   'pol_snums',tmp_parms.RF_pol_snums,'ecc_snums',tmp_parms.RF_ecc_snums,...
   'fnamestem',tmp_parms.RF_fnamestem,'multisess_flag',tmp_parms.multisess_flag);
  if errcode
    if parms.verbose
      fprintf('%s: WARNING: polar and eccentricity data not found for %s\n',...
        mfilename,VisitID);
    end;
    continue;
  elseif parms.verbose
    fprintf('%s: SUCCESS: polar and eccentricity data found for %s\n',...
      mfilename,VisitID);
  end;

  for h=1:length(parms.RF_hemilist)
    tmp_parms.RF_hemilist = parms.RF_hemilist(h);  
    jstem = regexprep(VisitID,'\^','_');
    jstem = jstem(1:min(20,length(jstem)));
    jobID = sprintf('job_%03d_%s_%s',j,jstem,parms.RF_hemilist{h}); j = j+1;
    jobfname = sprintf('%s/%s.m',batchdir,jobID);
    mmil_write_script(jobfname,'MMIL_RetFit_BOLD_Exam',...
      {ContainerPath,FSContainerPath},tags,tmp_parms);
    fid = fopen(scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
  end;
end

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

