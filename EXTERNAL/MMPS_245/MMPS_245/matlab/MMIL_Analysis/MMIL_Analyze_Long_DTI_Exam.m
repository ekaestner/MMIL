function MMIL_Analyze_Long_DTI_Exam(StudyInfo,ContainerPath,varargin)
%function MMIL_Analyze_Long_DTI_Exam(StudyInfo,ContainerPath,[options])
%
% Purpose: Analysis of Longitudinal DTI registration
%
% Required Input:
%   StudyInfo: StudyInfo struct containing all visits for a single subject
%   ContainerPath: full path of subject's longitudinal DTI directory
%
% Optional Input:
%   'FSPath': full path of FreeSurfer recon container (for painting to surface)
%     {default = []}
%   'resDTI_flag': [0|1] whether registration was done on DTI resolution images
%     or T1 resolution images
%     {default = 0}
%   'fibers': vector of fiber numbers to include in analysis
%     {default: [101:110,115:123,133:138,141:150,1011,1021]}
%   'surf_flag':[0/1] whether to paint DT measure and diff on cortical surface 
%     {default = 0}
%   'lesion_flag': [0/1] whether to create the lesion mask or not 
%     {default = 0}
%   'lesion_T1_mask_flag': [0|1] use T1 masks to create lesion mask (otherwise FA)
%     {default = 1}
%   'diff_smooth': smoothing sigma (voxels) applied to difference volumes
%    {default = 10}
%   'diff_scale': scaling factor applied to difference volumes
%    {default = 100}
%   'xcg_flag': [0|1] exclude CSF and gray-matter from fiber ROIs
%     {default = 1}
%   'masksf_flag': [0|1] exclude voxels with two or more fibers from fiber ROIs
%     {default = 0}
%   'forceflag': [0|1] run registrations even if output files exist
%     {default = 0}
%
% Created:  05/10/11 by Vijay Venkatraman
% Last Mod: 01/15/13 by Don Hagler
%

%% todo: xcg_suffix as input option

% based on EPDTI_MMIL_PrePost_Analyze_Exam, created 02/12/10 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'FSPath',[],[],...
  'resDTI_flag',false,[false true],...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'lesion_T1_mask_flag',true,[false true],...
  'lesion_flag',false,[false true],...
  'diff_smooth',10,[0,100],...
  'diff_scale',100,[-10000,10000],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'measlist',{'FA','MD','LD','TD','b0N','T2w'},[],...
  'scalefacts',[1,1000,1000,1000,1],[],...
  'thresh',0.5,[0,Inf],...
  'weighted_avg_flag',true,[false true],...
  'outext','.mgz',{'.mgh','.mgz'},...
  'lesion_thresh0',0.99,[0,100],...
  'lesion_smooth1',10,[0,100],...
  'lesion_thresh1',0.01,[0,100],...
  'lesion_smooth2',10,[0,100],...
  'lesion_thresh2',0.01,[0,100],...
  'lesion_smooth3',0,[0,100],...
  'combo_thresh',0.1,[0,100],...
...
  'projdist_list',[-1,1],[-5,5],...
  'dv_projdist',1,[-5,5],...
  'surf_flag',false,[false true],...
  'sphere_flag',false,[false true],...
  'surf_smoothsteps',0,[0,Inf],...
  'surf_sphsmoothsteps',0,[0,Inf],...
  'surf_mbmask_flag',false,[false true],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'fiber_outfix','fiber_roi_data',[],...
  'aparc_outfix','aparc_roi_data',[],...
});

if parms.xcg_flag
  parms.fiber_outfix = [parms.fiber_outfix '_xcg'];
end;
if parms.masksf_flag
  parms.fiber_outfix = [parms.fiber_outfix '_masksf'];
end;

if parms.resDTI_flag
  dir_infix = 'resDTI';
else
  dir_infix = 'resT1';
end;

if ~isempty(parms.FSPath)
  if parms.resDTI_flag
    fprintf('%s: WARNING: painting to surface with resDTI_flag not currently implemented\n',...
      mfilename);
    % NOTE: to paint to surface with resDTI_flag=1,
    %   need to create register.dat file from M_T1_to_DTI registration matrix
    %   using fs_write_regdat with ras2tk_flag=1
  elseif ~exist(parms.FSPath,'dir')
    fprintf('%s: WARNING: FreeSurfer recon %s not found\n',...
      mfilename,parms.FSPath);
  else
    [FSRootDir,FSContainerDir,tmp_ext] = fileparts(parms.FSPath);
    FSContainerDir = [FSContainerDir tmp_ext];
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort studies by VisitNumber
ind = find(cell2mat({StudyInfo.VisitNumber})>=1);
if length(ind)==0
  error('no visits with valid VisitNumber');
end;
StudyInfo = StudyInfo(ind);
VisitNumbers = cell2mat({StudyInfo.VisitNumber});
[VisitNumbers,ind] = sort(VisitNumbers);
StudyInfo = StudyInfo(ind);
if min(VisitNumbers)>1
  error('no VisitNumber=1');
end;

% loop over pairs of sessions
for i=1:length(StudyInfo)
  vA = StudyInfo(i).VisitNumber;
  sessA = StudyInfo(i).SessID;
  dirA = sprintf('%s/visit%d',ContainerPath,vA);
  if ~exist(dirA,'dir')
    fprintf('%s: WARNING: dir %s not found\n',mfilename,dirA);
    continue;
  end;
  for j=i+1:length(StudyInfo)
    vB = StudyInfo(j).VisitNumber;
    sessB = StudyInfo(j).SessID;
    dirB = sprintf('%s/visit%d',ContainerPath,vB);
    if ~exist(dirB,'dir')      
      fprintf('%s: WARNING: dir %s not found\n',mfilename,dirB);
      continue;
    end;

    % set nonlinreg dirs
    regdir_FA = sprintf('%s/nonlinreg_FA_%s_visit%d',...
      dirB,dir_infix,vA);
    if ~exist(regdir_FA,'dir')
      fprintf('%s: WARNING: dir %s not found\n',mfilename,regdir_FA);
      continue;
    end;
    regdir_T1 = sprintf('%s/nonlinreg_T1_%s_visit%d',...
      dirB,dir_infix,vA);
    if ~exist(regdir_T1,'dir')
      fprintf('%s: WARNING: dir %s not found\n',mfilename,regdir_T1);
      continue;
    end;

    % create outdir for analysis
    outdir = sprintf('%s/analysis_%s_visit%d',...
      dirB,dir_infix,vA);
    if ~exist(outdir,'dir')
      [succ,msg] = mkdir(outdir);
      if ~succ
        error('failed to create dir %s:\n%s',outdir,msg);
      end;      
    end;

    % apply 3D displacement to FA mask
    fname_in = sprintf('%s/images_%s/FA_mask%s',dirB,dir_infix,parms.outext);
    fname_out = sprintf('%s/FA_mask_TP2_NonLinReg%s',regdir_FA,parms.outext);
    fname_dx = sprintf('%s/dx%s',regdir_FA,parms.outext);
    fname_dy = sprintf('%s/dy%s',regdir_FA,parms.outext);
    fname_dz = sprintf('%s/dz%s',regdir_FA,parms.outext);
    if ~exist(fname_out,'file') || parms.forceflag
      fprintf('%s: warping FA mask to baseline...\n',mfilename);
      mmil_mriApply3DDispField(fname_in,fname_out,fname_dx,fname_dy,fname_dz);
    end;

    % to do: make mriRegister work for T1 non-linear registration
    % apply 3D displacement to T1 mask
    %fname_in = sprintf('%s/images_%s/T1_mask%s',dirB,dir_infix,parms.outext);
    %fname_out = sprintf('%s/T1_mask_TP2_NonLinReg%s',regdir_T1,parms.outext);
    %fname_dx = sprintf('%s/dx%s',regdir_T1,parms.outext);
    %fname_dy = sprintf('%s/dy%s',regdir_T1,parms.outext);
    %fname_dz = sprintf('%s/dz%s',regdir_T1,parms.outext);
    %if ~exist(fname_out,'file') || parms.forceflag
    %  fprintf('%s: warping T1 mask to baseline...\n',mfilename);
    %  mmil_mriApply3DDispField(fname_in,fname_out,fname_dx,fname_dy,fname_dz);
    %end;
    
    if parms.lesion_flag
    % create lesion mask
      if vA~=1
        % combine lesion masks for vA VS v1 and vB VS v1
        v1 = 1;
        sess1 = StudyInfo(1).SessID;
        outdirA = sprintf('%s/analysis_%s_visit%d',...
         dirA,dir_infix,v1);
        outdirB = sprintf('%s/analysis_%s_visit%d',...
         dirB,dir_infix,v1);
        if parms.lesion_T1_mask_flag
         fnameA = sprintf('%s/T1_lesion_mask%s',outdirA,parms.outext);
         fnameB = sprintf('%s/T1_lesion_mask%s',outdirB,parms.outext);
         fname_lesion = sprintf('%s/T1_lesion_mask%s',outdir,parms.outext);
        else
         fnameA = sprintf('%s/FA_lesion_mask%s',outdirA,parms.outext);
         fnameB = sprintf('%s/FA_lesion_mask%s',outdirB,parms.outext);
         fname_lesion = sprintf('%s/FA_lesion_mask%s',outdir,parms.outext);
        end;
        if ~exist(fnameA,'file')
         error('file %s not found',fnameA);
        end;
        if ~exist(fnameB,'file')
         error('file %s not found',fnameB);
        end;
      
        if ~exist(fname_lesion,'file') || parms.forceflag
          % resample lesion masks from v1 to vA
          if parms.lesion_T1_mask_flag
            fnameA_tmp = sprintf('%s/T1_lesion_mask_v%d%s',...
              outdir,vA,parms.outext);
            fnameB_tmp = sprintf('%s/T1_lesion_mask_v%d%s',...
              outdir,vB,parms.outext);
          else
           fnameA_tmp = sprintf('%s/FA_lesion_mask_v%d%s',...
              outdir,vA,parms.outext);
           fnameB_tmp = sprintf('%s/FA_lesion_mask_v%d%s',...
              outdir,vB,parms.outext);
          end;

          % find affine reg matrix for vA to v1
          v1 = StudyInfo(1).VisitNumber;
          sess1 = StudyInfo(1).SessID;
          dir1 = sprintf('%s/visit%d',ContainerPath,v1);
          [dir1_path,dir1_stem,dir1_ext] = fileparts(dir1);
          T1_outdir = [dirA '/nonlinreg_T1_' dir_infix '_'...
              dir1_stem,dir1_ext];
          fname_affine_mat = [T1_outdir '/affineRegMatrixS2T.txt'];
          if ~exist(fname_affine_mat,'file')
            error('affine registration matrix file %s not found',...
              fname_affine_mat);
          end;

          % resample masks into vA space
          mmil_mriApplyTransforms(fnameA,...
            'fname_out',fnameA_tmp,...
            'fname_mat',fname_affine_mat,...
            'forceflag',parms.forceflag);
          mmil_mriApplyTransforms(fnameB,...
            'fname_out',fnameB_tmp,...
            'fname_mat',fname_affine_mat,...
            'forceflag',parms.forceflag);

          % combine masks
          if parms.lesion_T1_mask_flag
           fprintf('%s: combining T1 lesion masks...\n',mfilename);
          else
           fprintf('%s: combining FA lesion masks...\n',mfilename);
          end;
          [volA,MA] = fs_load_mgh(fnameA_tmp);
          [volB,MB] = fs_load_mgh(fnameB_tmp);
          vol = sqrt(volA.*volB);
          fs_save_mgh(vol,fname_lesion,MA);
        end;
      else
        if parms.lesion_T1_mask_flag
         % calculate difference of T1 masks
          fnameA = sprintf('%s/images_%s/T1_mask%s',dirA,dir_infix,parms.outext);
          fnameB = sprintf('%s/T1_mask_TP2_AffReg%s',regdir_T1,parms.outext);
          fname_out = sprintf('%s/diff_T1_mask%s',outdir,parms.outext);
          fname_lesion = sprintf('%s/T1_lesion_mask%s',outdir,parms.outext);
        else
         % calculate difference of FA masks
         fnameA = sprintf('%s/images_%s/FA_mask%s',dirA,dir_infix,parms.outext);
         fnameB = sprintf('%s/FA_mask_TP2_NonLinReg%s',regdir_FA,parms.outext);
         fname_out = sprintf('%s/diff_FA_mask%s',outdir,parms.outext);
         fname_lesion = sprintf('%s/FA_lesion_mask%s',outdir,parms.outext);
        end;
        if ~exist(fnameA,'file')
          error('file %s not found',fnameA);
        end;
        if ~exist(fnameB,'file')
          error('file %s not found',fnameB);
        end;
        if ~exist(fname_out,'file') || parms.forceflag
          if parms.lesion_T1_mask_flag
            fprintf('%s: calculating T1 mask difference...\n',mfilename);
          else
            fprintf('%s: calculating FA mask difference...\n',mfilename);
          end;
          [volA,MA] = fs_load_mgh(fnameA);
          [volB,MB] = fs_load_mgh(fnameB);
          vol_diff = volA - volB;
          fs_save_mgh(vol_diff,fname_out,MA);
        end;
        % dilate lesion mask
        vol = ctx_load_mgh(fname_out);
        vol.imgs = abs(vol.imgs);
        mmil_dilate_mask(vol,...
          'fname_out',fname_lesion,...
          'thresh0',parms.lesion_thresh0,...
          'smooth1',parms.lesion_smooth1,...
          'thresh1',parms.lesion_thresh1,...
          'smooth2',parms.lesion_smooth2,...
          'thresh2',parms.lesion_thresh2,...
          'smooth3',parms.lesion_smooth3);
     end;
   end;
   
   % combine masks: FA first visit, FA second visit, T1 lesion
   fnameA = sprintf('%s/images_%s/FA_mask%s',dirA,dir_infix,parms.outext);
   fnameB = sprintf('%s/FA_mask_TP2_NonLinReg%s',regdir_FA,parms.outext); 
   fname_mask = sprintf('%s/combo_mask%s',outdir,parms.outext);
   [volA,MA] = fs_load_mgh(fnameA);
   [volB,MB] = fs_load_mgh(fnameB);
   if parms.lesion_flag 
     fnameC = fname_lesion;
     if ~exist(fname_mask,'file') || parms.forceflag
       fprintf('%s: combining masks...\n',mfilename);
       [volC,MB] = fs_load_mgh(fnameC);
       vol_mask = (volA>parms.combo_thresh) & ...
                  (volB>parms.combo_thresh) & ...
                  (volC<1-parms.combo_thresh);
       fs_save_mgh(vol_mask,fname_mask,MA);  
     end;
   else
     if ~exist(fname_mask,'file') || parms.forceflag
       fprintf('%s: combining masks...\n',mfilename);
       vol_mask = (volA>parms.combo_thresh) & ...
                   (volB>parms.combo_thresh);
       fs_save_mgh(vol_mask,fname_mask,MA);  
     end;
   end;

    
    % mask first visit fibers for each follow-up session
    fprintf('%s: masking visit A fibers...\n',mfilename);
    if ~exist(fname_mask,'file')
      error('file %s not found',fname_mask);
    end;
    indir = sprintf('%s/fibers_%s',dirA,dir_infix);
    outdir = sprintf('%s/fibers_%s_visit%d',dirB,dir_infix,vA);
    if ~exist(outdir,'dir')
      [succ,msg] = mkdir(outdir);
      if ~succ
        error('failed to create dir %s:\n%s',outdir,msg);
      end;      
    end;
    for f = 1:length(parms.fibers)
      suffix = [];
      if parms.xcg_flag
        suffix = '_xcg';
      end;
      if parms.masksf_flag
        suffix = [suffix '_masksf'];
      end;
      suffix = [suffix parms.outext];
      fname_in = sprintf('%s/fiber_%02d%s',indir,parms.fibers(f),suffix);
      fname_out = sprintf('%s/fiber_%02d%s',outdir,parms.fibers(f),suffix);
      if ~exist(fname_in,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
        continue;
      end;
      fs_maskvol(fname_in,fname_mask,fname_out,...
        'thresh',parms.thresh,'forceflag',parms.forceflag);
    end

    % calculate differences in DT measures
    outdir = sprintf('%s/analysis_%s_visit%d',dirB,dir_infix,vA);
    for m=1:length(parms.measlist)
      meas = parms.measlist{m};
      scalefact = parms.scalefacts(m);
      
      % apply displacements if not FA
      if ~strcmp(meas,'FA')
        fname_in = sprintf('%s/images_%s/%s%s',dirB,dir_infix,meas,parms.outext);
        fname_out = sprintf('%s/%s_TP2_NonLinReg%s',regdir_FA,meas,parms.outext);
        fname_dx = sprintf('%s/dx%s',regdir_FA,parms.outext);
        fname_dy = sprintf('%s/dy%s',regdir_FA,parms.outext);
        fname_dz = sprintf('%s/dz%s',regdir_FA,parms.outext);
        if ~exist(fname_in,'file') || ~exist(fname_dx,'file') || ~exist(fname_dy,'file') || ~exist(fname_dz,'file') 
          fprintf('%s: %s File not found',mfilename,fname_in);
          continue;
        end; 
        if ~exist(fname_out,'file') || parms.forceflag
          fprintf('%s: warping %s to baseline...\n',mfilename,meas);
          mmil_mriApply3DDispField(fname_in,fname_out,fname_dx,fname_dy,fname_dz);
        end;
      end;

      % calculate difference
      fnameA = sprintf('%s/images_%s/%s%s',dirA,dir_infix,meas,parms.outext);
      fnameB = sprintf('%s/%s_TP2_NonLinReg%s',regdir_FA,meas,parms.outext);
      fname_out = sprintf('%s/diff%s%s',outdir,meas,parms.outext);
      if ~exist(fnameA,'file')
        error('file %s not found',fnameA);
      end;
      if ~exist(fnameB,'file')
        error('file %s not found',fnameB);
      end;
      if ~exist(fname_out,'file') || parms.forceflag
        fprintf('%s: calculating %s differences...\n',mfilename,meas);
        [volA,MA] = fs_load_mgh(fnameA);
        [volB,MB] = fs_load_mgh(fnameB);
        vol_diff = volA - volB;
        fs_save_mgh(vol_diff,fname_out,MA);
      end;
      fname_diff = fname_out;

      % apply mask to meas for v1
      fname_out = sprintf('%s/%s_vA_masked%s',outdir,meas,parms.outext);
      fs_maskvol(fnameA,fname_mask,fname_out,...
        'thresh',parms.thresh,'forceflag',parms.forceflag);
      fname_vA = fname_out;

      % apply mask to meas for v2
      fname_out = sprintf('%s/%s_vB_masked%s',outdir,meas,parms.outext);
      fs_maskvol(fnameB,fname_mask,fname_out,...
        'thresh',parms.thresh,'forceflag',parms.forceflag);
      fname_vB = fname_out;

      % apply mask to diff
      fname_out = sprintf('%s/diff%s_masked%s',outdir,meas,parms.outext);
      fs_maskvol(fname_diff,fname_mask,fname_out,...
        'thresh',parms.thresh,'forceflag',parms.forceflag);
      fname_diff = fname_out;
      
      if parms.surf_flag  
        % paint meas and diff to cortical surface for white and gray matter
        if ~isempty(parms.FSPath) && ~parms.resDTI_flag && exist(parms.FSPath,'dir')
          fprintf('%s: painting %s and %s difference to cortical surface...\n',...
           mfilename,meas,meas);
          fname_list = {fname_vA,fname_vB,fname_diff};
          for f=1:length(fname_list)
            fname_in = fname_list{f};
            [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_in);
            for projdist = parms.projdist_list
              outstem = sprintf('%s/%s_pdist%0.1f',...
                tmp_path,tmp_stem,projdist);
              fs_paint(FSContainerDir,fname_in,...
                'outstem',outstem,...
                'sphere_flag',parms.sphere_flag,...
                'projdist',projdist,...
                'smoothsteps',parms.surf_smoothsteps,...
                'sphsmoothsteps',parms.surf_sphsmoothsteps,...
                'mask_midbrain_flag',parms.surf_mbmask_flag,...
                'subjdir',FSRootDir,...
                'overwrite_flag',parms.forceflag);

              % extract values for cortical ROIs
              fname_out = [outstem '_' parms.aparc_outfix '.mat'];
              if ~exist(fname_out,'file') || parms.forceflag
               for h=1:length(parms.hemilist)
                  hemi = parms.hemilist{h};
                  fname_aparc = sprintf('%s/label/%s.aparc.annot',parms.FSPath,hemi);
                  if ~exist(fname_aparc,'file'), error('%s not found',fname_aparc); end;
                  fname_in = sprintf('%s-%s.mgh',outstem,hemi);
                  aparc_roi_data = mmil_surf_roi(fname_in,...
                    'fname_aparc',fname_aparc,'minval',0);
                  if isempty(aparc_roi_data)
                    error('failed to get aparc ROI data');
                  end;
                  switch hemi
                    case 'lh'
                      lh_aparc_roi_data = aparc_roi_data;
                    case 'rh'
                      rh_aparc_roi_data = aparc_roi_data;
                  end;
                end;
                save(fname_out,'lh_aparc_roi_data','rh_aparc_roi_data');
              end;
            end;
         end;          
        end;
      end;

      % apply xcg mask
      if parms.xcg_flag
        fname_xcgmask = sprintf('%s/images_%s/xcg_mask%s',dirA,dir_infix,parms.outext);
        if ~exist(fname_xcgmask,'file')
          error('xcg mask file %s not found',fname_xcgmask);
        end;
        [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_diff);
        fname_out = sprintf('%s/%s_xcg%s',tmp_path,tmp_stem,tmp_ext);
        fs_maskvol(fname_diff,fname_xcgmask,fname_out,...
          'thresh',parms.thresh,'forceflag',parms.forceflag);
        fname_diff = fname_out;
      end;

      % extract mean meas and meas difference for each fiber
      fiberdir = sprintf('%s/fibers_%s_visit%d',dirB,dir_infix,vA);
      flist_fibers = {};
      suffix = [];
      if parms.xcg_flag
        suffix = '_xcg';
      end;
      if parms.masksf_flag
        suffix = [suffix '_masksf'];
      end;
      suffix = [suffix parms.outext];
      for f=1:length(parms.fibers)
        flist_fibers{f} = sprintf('%s/fiber_%02d%s%s',fiberdir,parms.fibers(f),suffix);
      end;
      fname_list = {fname_vA,fname_vB,fname_diff};
      outfix_list = {[meas '_vA'],[meas '_vB'],['diff' meas]};
      for f=1:length(fname_list)
        fname_in = fname_list{f};
        fname_out = sprintf('%s/%s_%s.mat',...
          outdir,outfix_list{f},parms.fiber_outfix);
        if ~exist(fname_out,'file') || parms.forceflag
          fprintf('%s: getting fiber data from %s...\n',mfilename,fname_in);
          fiber_data = mmil_multi_roi(fname_in,flist_fibers,...
            'weighted_avg_flag',parms.weighted_avg_flag);
          if isempty(fiber_data)
            error('fiber_data is empty');
          end;
          for f=1:length(parms.fibers)
            fiber_data(f).fibernum = parms.fibers(f);
          end;
          fprintf('%s: saving results to %s...\n',mfilename,fname_out);
          save(fname_out,'fiber_data');
        end;
      end;

      % smooth diff (for viewing with tractoview or tkmedit)
      if parms.diff_smooth>0
        sm = parms.diff_smooth;
        [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_diff);
        fname_out = sprintf('%s/%s_sm%d%s',...
          tmp_path,tmp_stem,parms.diff_smooth,tmp_ext);
        if ~exist(fname_out,'file') || parms.forceflag
          [vol,M] = fs_load_mgh(fname_diff);
          vol = mmil_smooth3d(vol,sm,sm,sm);
          fs_save_mgh(vol,fname_out,M);
        end;
        fname_diff = fname_out;    
      end;

      % scale diff (for viewing with tractoview or tkmedit)
      if parms.diff_scale~=1
        [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_diff);
        fname_out = sprintf('%s/%s_scale%d%s',...
          tmp_path,tmp_stem,parms.diff_scale,tmp_ext);
        if ~exist(fname_out,'file') || parms.forceflag
          [vol,M] = fs_load_mgh(fname_diff);
          vol = scalefact * parms.diff_scale * vol;
          fs_save_mgh(vol,fname_out,M);
        end;
        fname_diff = fname_out;    
      end;

      % save as mgh (for tractoview)
      if strcmp(parms.outext,'.mgz')
        [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_diff);
        fname_out = sprintf('%s/%s.mgh',tmp_path,tmp_stem);
        if ~exist(fname_out,'file') || parms.forceflag
          [vol,M] = fs_load_mgh(fname_diff);
          fs_save_mgh(vol,fname_out,M);
        end;
        fname_diff = fname_out;
      end;

      % resample to DTI resolution (for tractoview)
      if ~parms.resDTI_flag
        [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_diff);
        fname_out = sprintf('%s/%s_resDTI.mgh',tmp_path,tmp_stem);
        fname_DTI = sprintf('%s/images_resT1/%s%s',dirA,meas,parms.outext);
        fname_reg = sprintf('%s/reg_DTItoT1.mat',dirA);
        load(fname_reg);
        if ~exist('M_T1_to_EPI','var')
            error('M_T1_to_EPI not found in %s',fname_reg);
        end;    
        if ~exist(fname_out,'file') || parms.forceflag
          fprintf('%s: resampling diff back to DTI resolution...\n',mfilename);
          vol = ctx_load_mgh(fname_diff);
          vol_DTI = ctx_load_mgh(fname_DTI);
          vol_res = vol_resample_pad(vol,vol_DTI,inv(M_T1_to_EPI),1,0);
          ctx_save_mgh(vol_res,fname_out);
        end;
        fname_diff = fname_out;
      end;
      %% todo: create csh script to run tractoview?
      %% todo: create csh script to run tkmedit to compare diff and fseg
      %% todo: create csh script to run tksurfer to view diff +- pdist
    end; % loop over measlist
    
    if parms.surf_flag
     % paint gray matter dv to cortical surface
     if ~isempty(parms.FSPath) && ~parms.resDTI_flag && exist(parms.FSPath,'dir')
        fprintf('%s: painting dv to cortical surface...\n',mfilename);
        [FSRootDir,FSContainerDir,tmp_ext] = fileparts(parms.FSPath);
        FSContainerDir = [FSContainerDir tmp_ext];
        fname_list = {...
          sprintf('%s/dv%s',regdir_FA,parms.outext) ...
          sprintf('%s/dv%s',regdir_T1,parms.outext) ...
        };
        outstem_list = {'dvFA','dvT1'};
        for f=1:length(fname_list)
          fname_in = fname_list{f};
          outstem = sprintf('%s/%s',outdir,outstem_list{f});
          fs_paint(FSContainerDir,fname_in,...
            'outstem',outstem,...
            'sphere_flag',parms.sphere_flag,...
            'projdist',parms.dv_projdist,...
            'smoothsteps',parms.surf_smoothsteps,...
            'sphsmoothsteps',parms.surf_sphsmoothsteps,...
            'mask_midbrain_flag',parms.surf_mbmask_flag,...
            'subjdir',FSRootDir,...
            'overwrite_flag',parms.forceflag);

          % extract values for cortical ROIs
          fname_out = [outstem '_' parms.aparc_outfix '.mat'];
          if ~exist(fname_out,'file') || parms.forceflag
            for h=1:length(parms.hemilist)
              hemi = parms.hemilist{h};
              fname_aparc = sprintf('%s/label/%s.aparc.annot',parms.FSPath,hemi);
              if ~exist(fname_aparc,'file'), error('%s not found',fname_aparc); end;
              fname_in = sprintf('%s-%s.mgh',outstem,hemi);
              aparc_roi_data = mmil_surf_roi(fname_in,...
                'fname_aparc',fname_aparc,'minval',0);
              if isempty(aparc_roi_data)
                error('failed to get aparc ROI data');
              end;
              switch hemi
                case 'lh'
                 lh_aparc_roi_data = aparc_roi_data;
                case 'rh'
                 rh_aparc_roi_data = aparc_roi_data;
              end;
            end;
            save(fname_out,'lh_aparc_roi_data','rh_aparc_roi_data');
          end;
        end;
      end; % dv
     end;  
  end;
end;

