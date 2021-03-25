function MMIL_Long_Setup_Exam(RootDirs,StudyInfo,OutContainerPath,varargin)
%function MMIL_Long_Setup_Exam(RootDirs,StudyInfo,OutContainerPath,[options])
%
% Purpose: Setup file structure, copy files, create masks for 
%          longitudinal structural MRI and DTI analysis
%
% Required Input:
%   RootDirs: RootDirs struct containing full paths of project data directories
%   StudyInfo: StudyInfo struct containing all visits for a single subject
%   OutContainerPath: full path of output directory (will create if does not exist)
%
% Optional Input:
%  'T1type': which type of T1 series ('MPR' or 'hiFA') to use
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%  'LongDTI_flag': [0|1] setup for longitudinal DTI analysis
%     {default: 0}
%  'forceflag': if output files exist, delete them and then run setup
%     {default: 0}
%
% Created:  05/19/11 by Vijay Venkatraman
% Rcnt Mod: 08/29/13 by Don Hagler
% Last Mod: 09/25/13 Chun C. Fan
%

% based on code by Don Hagler

if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin, { ...
  'LongDTI_flag',false,[false true],...
  'forceflag',false,[false true],...
  'T1type',2,[0 1 2 3]...
...
  'snums_flag',0,[0:3],...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'nob0_flag',false,[false true],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
  'min_nb0',1,[],...
  'nonlin_flag',false,[false true],...
  'infix','corr_regT1',[],...
  'revflag',0,[0:2],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
  'outext','.mgz',{'.mgh','.mgz'},...
  'fiberdir_resT1','AtlasTrack/fiber_maps_resT1',[],...
  'atlas_flag',2,[0,1,2,3],...
  'atl_tensor_smooth_sigma',5,[],...
  'first_only_flag',true,[false true],...
  'count_flag',true,[false true],...
  'atlasname',[],[],...
  'resT1flag',true,[false true],...
  'DT_regT1flag',2,[0 1 2],...
  'lesion_flag',false,[true false]...
  'atlasname',[],[],...
...
  'track_tags',{'snums','infix','revflag','resT1flag',...
                'xcg_flag','masksf_flag'},[],...
  'infix_tags',{'atlas_flag','xcg_flag','masksf_flag','resT1flag'},[],...
  'info_tags',{'snums','revflag','min_nb0','min_ndirs',...
               'min_bval','flex_flag'},[],...
  'fstem_tags',{'snums','infix','revflag','min_bval','flex_flag','min_ndirs',...
                'min_nb0','nob0_flag'},[],...
  'calcDT_tags',{'snums','infix','revflag','nob0_flag','min_bval','flex_flag',...
                 'min_ndirs','min_nb0','nonlin_flag','regT1flag','forceflag'},[],...
});

parms.regT1flag = parms.DT_regT1flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(parms.atlasname)
  parms.fiberdir_resT1 = [parms.fiberdir_resT1 '_' parms.atlasname];
end;

% set fiber infix
args = mmil_parms2args(parms,parms.infix_tags);
parms.fiber_infix_resT1 = dti_set_fiber_infix(args{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that study info is for single subject
SubjID = unique({StudyInfo.SubjID});
if length(SubjID) > 1
  error('multiple subjects not allowed in StudyInfo');
else
  SubjID = SubjID{1};
end;

% exclude visits with visitnum<1
ind = find(cell2mat({StudyInfo.VisitNumber})>=1);
StudyInfo = StudyInfo(ind);
if length(StudyInfo)==0
  error('no visits with valid VisitNumber');
end;

all_input_fnames = {};
all_output_fnames = {};

% loop over visits
for i=1:length(StudyInfo)
  % check DTI proc container exists
  if parms.LongDTI_flag
    DTIContainerPath = [RootDirs.proc_dti '/' StudyInfo(i).proc_dti];
    if isempty(StudyInfo(i).proc_dti) | ~exist(DTIContainerPath,'dir')
      if StudyInfo(i).VisitNumber == 1
        error('baseline proc container %s not found',DTIContainerPath);
      else
        fprintf('%s: WARNING: followup container %s not found\n',...
          mfilename,DTIContainerPath);
        continue;
      end;
    end;
  end;
  % check T1 proc container exists
  T1ContainerPath = ...
    MMIL_Get_Container(RootDirs,StudyInfo(i).STRUCT_VisitID,'proc');
  if ~exist(T1ContainerPath,'dir')
     if StudyInfo(i).VisitNumber == 1
      error('baseline proc container %s not found',T1ContainerPath);
    else
      fprintf('%s: WARNING: followup container %s not found\n',...
        mfilename,T1ContainerPath);
      continue;
    end;
  end;
  % check freesurfer container exists
  FSContainerPath = [RootDirs.fsurf '/' StudyInfo(i).fsurf];
  if isempty(StudyInfo(i).fsurf) | ~exist(FSContainerPath,'dir')
    if StudyInfo(i).VisitNumber == 1
      error('baseline fsurf container %s not found',FSContainerPath);
    else
      if parms.lesion_flag | (parms.LongDTI_flag & parms.xcg_flag)
        error('followup fsurf container for %s not found in %s',...
          StudyInfo(i).STRUCT_VisitID,FSContainerPath);
      elseif ~isempty(StudyInfo(i).fsurf)
        fprintf('%s: WARNING: followup fsurf container %s not found\n',...
          mfilename,FSContainerPath);
      end;
    end;
  end;
  
  if parms.LongDTI_flag
    % get snums used for DT calculations
    parms.snums = [];
    if parms.snums_flag
      switch parms.snums_flag
        case 1
          parms.snums = StudyInfo(i).DTIScanNums;
        case 2
          parms.snums = StudyInfo(i).DTIScanNums2;
        case 3
          parms.snums = StudyInfo(i).DTIScanNums3;
      end;
      if isempty(parms.snums), continue; end;
    end;
    args = mmil_parms2args(parms,parms.info_tags);
    [ScanInfo,SessInfo,errcode] = ...
      DTI_MMIL_Get_ScanInfo(DTIContainerPath,args{:});
    parms.snums = SessInfo.snums_DT;
    if isempty(parms.snums)
      if StudyInfo(i).VisitNumber == 1
        error('no valid DTI scan numbers for baseline container %s',DTIContainerPath);
      else
        fprintf('%s: WARNING: no valid DTI scan numbers for followup container %s\n',...
          DTIContainerPath);
        continue;
      end;
    end;
    % for convenience, run DTI_MMIL_AtlasTrack_Exam
    %   to make sure fibers exist (i.e. with resT1, xcg, masksf)
    args = mmil_parms2args(parms,parms.track_tags);
    DTI_MMIL_AtlasTrack_Exam(DTIContainerPath,...
      'FSContainerPath',FSContainerPath,...
      'forceflag',0,args{:});
  end;

  % check if output exists
  outdir = sprintf('%s/visit%d',OutContainerPath,StudyInfo(i).VisitNumber);
  if exist(outdir) && parms.forceflag
    fprintf('%s: forcing complete redo by removing existing directory %s...\n',mfilename,outdir);
    cmd = ['rm -r ' outdir];
    [s,r] = unix(cmd);
    if s, error('cmd %s failed:\n%s',cmd,r); end;
  end;

  % create output directories
  fprintf('%s: creating directories for %s...\n',mfilename,StudyInfo(i).SessID);
  outdir_images_resT1 = [outdir '/images_resT1'];
  mmil_mkdir(outdir_images_resT1);

  % choose T1 file
  [fname_in,errcode] = MMIL_Choose_T1(T1ContainerPath,'T1type',parms.T1type);
  if errcode, return; end;

  input_fnames = {fname_in};
  output_fnames = {[outdir_images_resT1 '/T1' parms.outext]};

  if StudyInfo(i).VisitNumber == 1 || ~isempty(StudyInfo(i).fsurf)
    input_fnames = cat(2,input_fnames,...
                         [FSContainerPath '/mri/aseg.mgz']);
    output_fnames = cat(2,output_fnames,...
                          [outdir_images_resT1 '/aseg' parms.outext]);
  end;
  if parms.LongDTI_flag
    outdir_fibers_resT1 = [outdir '/fibers_resT1'];
    mmil_mkdir(outdir_fibers_resT1);
    % set fstem
    args = mmil_parms2args(parms,parms.fstem_tags);
    fstem = DTI_MMIL_Set_DT_fstem(DTIContainerPath,args{:});
    imgtypes = {'b0','FA','V0','MD','LD','TD'};
    for j=1:length(imgtypes)
      if ~exist([fstem '_' imgtypes{j} '_resT1.mgz'],'file')
        args = mmil_parms2args(parms,parms.calcDT_tags);
        DTI_MMIL_Calc_DT(DTIContainerPath,args{:});
      end;
    end;
    [RegInfo,fname_reg,errcode]= DTI_MMIL_Load_RegInfo(DTIContainerPath,...
      'revflag',parms.revflag,'infix',parms.infix);
    % set input/output file names
    input_fnames = cat(2,input_fnames,{...
      [fstem '_FA_resT1.mgz'] ...
      [fstem '_MD_resT1.mgz'] ...
      [fstem '_LD_resT1.mgz'] ...
      [fstem '_TD_resT1.mgz'] ...
      [fstem '_V0_resT1.mgz'] ...
      [fstem '_b0_resT1.mgz'] ...
      [fstem '_b0N_resT1.mgz'] ...
      [fname_reg] ...
    });
    output_fnames = cat(2,output_fnames,{...
      [outdir_images_resT1 '/FA' parms.outext] ...
      [outdir_images_resT1 '/MD' parms.outext] ...
      [outdir_images_resT1 '/LD' parms.outext] ...
      [outdir_images_resT1 '/TD' parms.outext] ...
      [outdir_images_resT1 '/V0' parms.outext] ...
      [outdir_images_resT1 '/b0' parms.outext] ...
      [outdir_images_resT1 '/b0N' parms.outext] ...
      [outdir '/reg_DTItoT1.mat'] ...
    });
    if parms.xcg_flag
      input_fnames = cat(2,input_fnames,{...
        [DTIContainerPath '/fseg_resT1_xcg.mgz'] ...
        [DTIContainerPath '/aseg_xcg_mask.mgz'] ...
      });
      output_fnames = cat(2,output_fnames,{...
        [outdir_images_resT1 '/fseg' parms.outext] ...
        [outdir_images_resT1 '/xcg_mask' parms.outext] ...
      });
    else
      input_fnames = cat(2,input_fnames,{...
        [DTIContainerPath '/fseg_resT1.mgz'] ...
      });
      output_fnames = cat(2,output_fnames,{...
        [outdir_images_resT1 '/fseg' parms.outext] ...
      });
    end;
  end;
  % copy files
  fprintf('%s: copying files for %s...\n',mfilename,StudyInfo(i).SessID);
  for f=1:length(input_fnames)
    if exist(output_fnames{f},'file'), continue; end;
    if ~exist(input_fnames{f},'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,input_fnames{f});
      continue;
    end;
    copy_file(input_fnames{f},output_fnames{f});
  end;

  if parms.LongDTI_flag
    % copy fibers
    fprintf('%s: copying fibers for %s...\n',mfilename,StudyInfo(i).SessID);
    for f=parms.fibers
      % copy resT1 fiber, The _pthresh0.08 part is the one generated in previous proc step
      fname_in = sprintf('%s/%s/fiber_%03d_%s_pthresh0.08.mat',... 
        DTIContainerPath,parms.fiberdir_resT1,f,parms.fiber_infix_resT1);
      suffix = [];
      if parms.xcg_flag
        suffix = '_xcg';
      end;
      if parms.masksf_flag
        suffix = [suffix '_masksf'];
      end;
      suffix = [suffix parms.outext];      
      fname_out = sprintf('%s/fiber_%03d%s',outdir_fibers_resT1,f,suffix);
      if exist(fname_out,'file'), continue; end;
      if ~exist(fname_in,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
        continue;
      end;
      copy_file(fname_in,fname_out);
    end;
    % create masks for resT1 files
    fprintf('%s: creating masks for %s...\n',mfilename,StudyInfo(i).SessID);
    fname_aseg = [outdir_images_resT1 '/aseg' parms.outext];
    fname_T1_mask = [outdir_images_resT1 '/T1_mask' parms.outext];
    fname_FA =  [outdir_images_resT1 '/FA' parms.outext];
    fname_FA_mask =  [outdir_images_resT1 '/FA_mask' parms.outext];
    fname_MD =  [outdir_images_resT1 '/MD' parms.outext];
    fname_MD_mask =  [outdir_images_resT1 '/MD_mask' parms.outext];
    if exist(fname_aseg,'file')
      MMIL_Long_Mask_T1(fname_aseg,fname_T1_mask,parms.forceflag);
    else
      sprintf('%s: WARNING: unable to create T1 mask for resT1 (missing ASEG!)\n',mfilename);
    end;
    if exist(fname_FA,'file') && exist(fname_T1_mask,'file')
      MMIL_Long_Mask_DTI(fname_FA,fname_T1_mask,fname_FA_mask,parms.forceflag)
    else
      sprintf('%s: WARNING: unable to create FA mask for resT1\n',mfilename);
    end;
    if exist(fname_MD,'file') && exist(fname_T1_mask,'file')
      MMIL_Long_Mask_DTI(fname_MD,fname_T1_mask,fname_MD_mask,parms.forceflag)
    else
      sprintf('%s: WARNING: unable to create MD mask for resT1\n',mfilename);
    end;
  elseif 0 % not needed for QUARC, but may be useful in future
    % create masks
    fprintf('%s: creating masks for %s...\n',mfilename,StudyInfo(i).SessID);
    fname_aseg = [outdir_images '/aseg' parms.outext];
    fname_T1_mask = [outdir_images '/T1_mask' parms.outext];
    if exist(fname_aseg,'file')
      MMIL_Long_Mask_T1(fname_aseg,fname_T1_mask,parms.forceflag);
    else
      sprintf('%s: WARNING: unable to create mask (missing aseg!)\n',mfilename);
    end;
  end;

   all_input_fnames = cat(2,all_input_fnames,input_fnames);
   all_output_fnames = cat(2,all_output_fnames,output_fnames);
end;


% create ContainerInfo
fname_info = [OutContainerPath '/ContainerInfo.mat'];
ContainerInfo = [];
ContainerInfo.SubjID = StudyInfo(i).SubjID;
ContainerInfo.StudyInfo = StudyInfo;
ContainerInfo.input_fnames = all_input_fnames;
ContainerInfo.output_fnames = all_output_fnames;
save(fname_info,'ContainerInfo');


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_file(fname_in,fname_out,outext)
  [tmpA,tmpB,ext_in] = fileparts(fname_in);
  [tmpA,tmpB,ext_out] = fileparts(fname_out);
  if ~strcmp(ext_in,ext_out) & ismember(ext_in,{'.mgz','.mgh'})
    [vol,M,mr_parms] = fs_load_mgh(fname_in);
    fs_save_mgh(vol,fname_out,M,mr_parms);
  elseif ~strcmp(ext_in,ext_out) & ismember(ext_in,{'.mat'})
    [vol,M] = mmil_load_sparse(fname_in);
    fs_save_mgh(vol,fname_out,M);
  else
    cmd = sprintf('cp %s %s',fname_in,fname_out);
    [s,r]=unix(cmd);
    if s, error('cmd %s failed:\n%s',cmd,r); end;
  end;
return;
