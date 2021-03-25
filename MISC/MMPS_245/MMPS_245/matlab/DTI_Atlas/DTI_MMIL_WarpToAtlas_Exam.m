function DTI_MMIL_WarpToAtlas_Exam(ContainerPath,varargin);
%function DTI_MMIL_WarpToAtlas_Exam(ContainerPath,varargin);
%
% Required Input:
%  ContainerPath: processed container path
%
% Optional Parameters:
%  'outdir': output directory
%     full path or relative to ContainerPath
%     {default = 'WarpToAtlas'}
%  'fname_T1': file name of T1-weighted image volume
%     if empty, will look for MPR_res.mgz or hiFA_res.mgz in T1ContainerPath
%       or ContainerPath if T1ContainerPath not supplied
%     (preference depends on STRUCT_T1type)
%     {default = []}
%  'STRUCT_T1type': which T1 series to reg the DTI to
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default = 2}
%  'T1ContainerPath': full path of directory containing T1-weighted images
%     (ignored if fname_T1 is also supplied)
%     {default = []}
%  'FSContainerPath': full path of directory containing freesurfer recon
%     will use mri/nu.mgz for fname_T1
%     (ignored if fname_T1 is also supplied)
%     {default = []}
%  'DTI_snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%  'DTI_min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%  'DTI_nob0_flag': [0|1] whether to exclude b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'DTI_min_ndirs': minimum number of directions to be included in tensor fit 
%     {default = 6}
%  'DTI_min_nb0': minimum number of B0 images
%     {default = 1}
%  'DTI_fiber_infix': infix of fibers used to warp to atlas
%     {default = 'corr_regT1'}
%  'DTI_revflag': whether rev scans were used for tensor calculations
%     {default = 0}
%  'DTI_fibers': Fibers to used warped to atlas
%     {default = [101:110,115:123,133:138,141:150]}
%  'DTI_measlist': DTI measures to be warped to atlas
%     {default = {'FA','b0'}}
%  'fiber_countflag': [0|1] whether fiber was imported as count map or mask map
%     {default = 1}
%  'sparse_flag': [0|1] save volumes as sparse mat files
%     {default = 1}
%
% Optional Parameters controlling processing steps:
%  'warp_masks_flag': [0|1] whether to warp brain masks to atlas
%     {default = 1}
%  'warp_tensors_flag': [0|1] whether to warp tensors, eigen vectors,
%                              and eigen values to atlas
%     {default = 1}
%  'warp_DTmeas_flag': [0|1] whether to warp DTI measures (FA, etc.) to atlas
%     {default = 1}
%  'warp_fibers_flag': [0|1] whether to warp fibers to atlas
%     {default = 1}
%  'warp_fiber_tensors_flag': [0|1] whether to smooth tensors within fibers
%                              and warp to atlas
%     {default = 1}
%  'first_only_flag': [0|1] whether to remove 2nd and 3rd Eigenvectors
%     {default = 1}
%  'atlasfibers_flag': [0|1] whether to use original, manual fibers (0)
%                             or new fibers from atlas (1) {default = 0}
%  'smfiber_sigma',[0,Inf]: sigma for smoothing the fiber (mm)
%     {default = 10}
%  'smfiber_win',[0,Inf]: window for smoothing the fiber (voxels)
%     {default = 3}
% 
% parameters for smoothing aseg mask
%  'smooth1': smoothing parameter 1 for aseg mask
%     {default = 20}
%  'thresh1': threshold parameter 1 for aseg mask
%     {default = 0.5}
%  'smooth2': smoothing parameter 2 for aseg mask
%     {default = 40}
%  'thresh2': threshold parameter 2 for aseg mask
%     {default = 0.2}
%  'smooth3': smoothing parameter 3 for aseg mask
%     {default = 10}
%  'thresh3': threshold parameter 3 for aseg mask
%     {default = 0.01}
%  'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/28/11 by Vijay Venkatraman
% Last Mod: 06/07/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
...% parameters for selecting data
  'DTI_min_ndirs',6,[],...
  'DTI_snums',[],[],...
  'DTI_min_nb0',1,[],...
  'DTI_fiber_infix','corr_regT1',[],...
  'DTI_min_bval',1,[],...
  'DTI_flex_flag',false,[false true],...
  'DTI_nob0_flag',false,[false true],...
  'STRUCT_T1type',2,[],... 
  'DTI_revflag',0,[0,1,2],...
  'DTI_fibers',[101:110,115:123,133:138,141:150],[],...
  'DTI_measlist',{'FA','b0'},[],...
  'ext','.mgz',{'.mgh','.mgz'},...
...  
  'outdir','WarpToAtlas',[],...
  'fname_T1',[],[],...
  'T1ContainerPath',ContainerPath,[],...
  'FSContainerPath',[],[],...
...
  'fiber_countflag',true,[false true],...
  'sparse_flag',true,[false true],...
  'warp_masks_flag',true,[false true],...
  'warp_tensors_flag',true,[false true],...
  'warp_DTmeas_flag',true,[false true],...
  'warp_fibers_flag',true,[false true],...
  'warp_fiber_tensors_flag',true,[false true],...
  'first_only_flag',true,[false true],...
  'atlasfibers_flag',false,[false true],...
... % parameters for smoothing tensors within fibers
  'smfiber_sigma',10,[0,Inf],... % mm
  'smfiber_win',3,[0,Inf],... % voxels
... % parameters for smoothing aseg mask
  'smooth1',20,[0,100],...
  'thresh1',0.5,[0,1],...
  'smooth2',40,[0,100],...
  'thresh2',0.2,[0,1],...
  'smooth3',10,[0,100],...
  'thresh3',0.01,[0,1],...
...
  'forceflag',false,[false,true],...
... % hidden parameters
  'fiber_maps_infix','sm5.00_probV0_countatlas',[],... %% todo: should be constructed from other input parms
  'atlasname',[],[],...
  'tensor_components',[1,5,9,2,3,6],[1:9],...  % do not change
});

parms.T1type = parms.STRUCT_T1type;
parms.num_tensor_components = length(parms.tensor_components);

if mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;

if parms.sparse_flag
  parms.outext = '.mat';
else
  parms.outext = '.mgz';
end;

% output directories
dir_tensor_atlas = sprintf('%s/tensor_atlas',parms.outdir);
dir_DTmeas_atlas = sprintf('%s/DTmeas_atlas',parms.outdir);
dir_masks_atlas = sprintf('%s/masks_atlas',parms.outdir);

if parms.atlasfibers_flag
  % input fiber directory
  if isempty(parms.atlasname)
    dir_fiber = sprintf('%s/AtlasTrack/fiber_maps',ContainerPath);
  else
    dir_fiber = sprintf('%s/AtlasTrack_%s/fiber_maps',ContainerPath,parms.atlasname);
  end;
  % output directories
  if isempty(parms.atlasname)
    dir_fiber_resT1 = sprintf('%s/AtlasTrack/fiber_maps_resT1',ContainerPath);
    dir_fiber_atlas = sprintf('%s/fiber_fromatlas_atlas',parms.outdir);
  else
    dir_fiber_resT1 = sprintf('%s/AtlasTrack_%s/fiber_maps_resT1',ContainerPath,parms.atlasname);
    dir_fiber_atlas = sprintf('%s/fiber_fromatlas_atlas_%s',parms.outdir,parms.atlasname);
  end; 
else
  dir_fiber = sprintf('%s/DTIStudio_fiber_masks',ContainerPath); 
  dir_fiber_resT1 = sprintf('%s/fiber_resT1',ContainerPath);
  dir_fiber_atlas = sprintf('%s/fiber_atlas',parms.outdir);
end;

if parms.warp_fibers_flag && ~exist(dir_fiber,'dir')
  fprintf('%s: WARNING: fiber dir %s not found\n',mfilename,dir_fiber);
  parms.warp_fibers_flag = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preliminaries

% check that processed data container exists
if ~exist(ContainerPath,'dir')
  error('Processed container path %s not found',ContainerPath);
end;
mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose T1 file
tags = {'fname_T1','T1type','FSContainerPath'};
args = mmil_parms2args(parms,tags);
[parms.fname_T1,errcode] = MMIL_Choose_T1(parms.T1ContainerPath,args{:});
if errcode, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that DT calculations exist
[DT_fstem,parms.DTI_snums] = DTI_MMIL_Set_DT_fstem(ContainerPath,...
  'snums',parms.DTI_snums,'infix',parms.DTI_fiber_infix,...
  'min_ndirs',parms.DTI_min_ndirs,'min_bval',parms.DTI_min_bval,...
  'flex_flag',parms.DTI_flex_flag,...
  'revflag',parms.DTI_revflag,...
  'min_nb0',parms.DTI_min_nb0);
if isempty(DT_fstem)
  fprintf('%s: ERROR: no valid DTI scans\n',mfilename);
  return;
end;
fname_DTout = [DT_fstem '_meas.mat'];
fname_Fitout = [DT_fstem '_fit.mat'];
if ~exist(fname_DTout,'file')
  error('file %s not found, run DTI_MMIL_CalcDT first',fname_DTout);
end;
if ~exist(fname_Fitout,'file')
  error('file %s not found, run DTI_MMIL_CalcDT first',fname_Fitout);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load registration info
clear RegInfo;
[RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,'infix',...
    parms.DTI_fiber_infix,'revflag',parms.DTI_revflag);
M_T1_to_EPI= RegInfo.M_T1_to_T2;
if ~exist('RegInfo','var'),
  error('RegInfo not found in %s',fname_reg);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear M M_DTI DTmeas DTfit;
% local copy is not necessary, but maybe useful for something
fname_tmp = sprintf('%s/T1%s',parms.outdir,parms.ext);
[vol,M,mr_parms] = fs_load_mgh(parms.fname_T1);
fs_save_mgh(vol,fname_tmp,M,mr_parms);
parms.fname_T1 = fname_tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% register T1-weighted MRI to atlas and warp T1 to atlas
[fname_atl,fname_atlreg] = mmil_warp_to_atlas(parms.fname_T1,'fname_T1',parms.fname_T1,...
  'outdir',parms.outdir,'forceflag',parms.forceflag);
load(fname_atlreg);
M_Atlas_to_DTI = M_T1_to_EPI*M_Atlas_to_Subj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% warp masks to atlas
if parms.warp_masks_flag
  mmil_mkdir(dir_masks_atlas);

  % warp MPR_res_mask to atlas
  fname_mask = [parms.T1ContainerPath '/MPR_res_mask' parms.ext];
  if exist(fname_mask,'file')
      mmil_warp_to_atlas(fname_mask,'fname_T1',parms.fname_T1,...
        'outdir',dir_masks_atlas,...
        'fname_reg',fname_atlreg,...
        'forceflag',parms.forceflag);
  end;

  % warp hiFA_res_mask to atlasfiber_maps_infix
  fname_mask = [parms.T1ContainerPath '/hiFA_res_mask' parms.ext];
  if exist(fname_mask,'file')
      mmil_warp_to_atlas(fname_mask,'fname_T1',parms.fname_T1,...
        'outdir',dir_masks_atlas,...
        'fname_reg',fname_atlreg,...
        'forceflag',parms.forceflag);
  end;

  % warp filled aseg mask to atlas
  fname_aseg = [parms.FSContainerPath '/mri/aseg.mgz'];
  if ~isempty(parms.FSContainerPath) & exist(fname_aseg,'file')
    fname_mask = [dir_masks_atlas '/aseg_mask.mgz'];
    mmil_dilate_mask([],'fname_in',fname_aseg,'fname_out',fname_mask,...
      'smooth1',parms.smooth1,'thresh1',parms.thresh1,...
      'smooth2',parms.smooth2,'thresh2',parms.thresh2,...
      'smooth3',parms.smooth3,'thresh3',parms.thresh3,...
      'forceflag',parms.forceflag);
    mmil_warp_to_atlas(fname_mask,'fname_T1',parms.fname_T1,...
      'outdir',dir_masks_atlas,...
      'fname_reg',fname_atlreg,...
      'forceflag',parms.forceflag);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load DTI measures, tensor data and T1 data
if ~exist('DTmeas','var') | ~exist('M_DTI','var') | ~exist('DTfit','var')
  load(fname_DTout); load(fname_Fitout);
  if ~exist('DTmeas','var')
   error('missing DTmeas from %s',fname_DTout);
  end;
  if ~exist('DTfit','var')
   error('missing DTfit from %s',fname_Fitout);
  end;
  if ~exist('M','var')
   error('missing M from %s',fname_DTout);
  end;
  M_DTI = M;
end;
if ~exist('vol_T1','var') | isempty(vol_T1)
  % load freesurfer T1-weighted T1 volume or fname_T1
  fprintf('%s: loading T1 volume %s...\n',mfilename,parms.fname_T1);
  [vol_T1,M_T1,mr_parms,volsz_T1] = fs_load_mgh(parms.fname_T1);
  vol_T1 = ctx_mgh2ctx(vol_T1,M_T1);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotate tensors and eigen vectors and warp to atlas
if parms.warp_tensors_flag
  mmil_mkdir(dir_tensor_atlas);
 [tmp_path,tmp_DT_fstem] = fileparts(DT_fstem);
  fname_tensor_atl = sprintf('%s/%s_tensor_atlas%s',...
    dir_tensor_atlas,tmp_DT_fstem,parms.outext);
  fname_eigvect_atl = sprintf('%s/%s_eigvect_atlas%s',...
    dir_tensor_atlas,tmp_DT_fstem,parms.outext);
  fname_eigval_atl = sprintf('%s/%s_eigval_atlas%s',...
    dir_tensor_atlas,tmp_DT_fstem,parms.outext);
  if ~exist(fname_tensor_atl,'file') |...
     ~exist(fname_eigvect_atl,'file') | ...
     ~exist(fname_eigval_atl,'file') | parms.forceflag
     if ~exist(fname_tensor_atl,'file') | parms.forceflag
      % rotate tensors
      fprintf('%s: rotating tensors to atlas...\n',mfilename);
      %% todo: verify that this is correct:
      volT = dti_rotate_tensors(DTmeas.volT,inv(M_Atlas_to_DTI));
      if isempty(volT)
        error('failed to rotate tensors');
      end;
      volT = reshape(volT,...
        [size(volT,1),size(volT,2),size(volT,3),9]);
      % warp tensors to atlas
      fprintf('%s: warping tensors to atlas...\n',mfilename);
      volT_atl = zeros([volsz_T1(1:3),parms.num_tensor_components]);
      for i=1:parms.num_tensor_components
        t = parms.tensor_components(i);
        vol = squeeze(volT(:,:,:,t));
        [vol,volres,M_atl] = applywarp(vol,M_DTI,vol_T1,M_T1_to_EPI,1,0,1,regStruct); 
        %interpm = 1;bclamp = 0;padding = 1
        volT_atl(:,:,:,i) = vol;
      end;
      if parms.sparse_flag
        mmil_save_sparse(volT_atl,fname_tensor_atl,M_atl);
      else
        fs_save_mgh(volT_atl,fname_tensor_atl,M_atl);
      end;
      clear volT volT_atl vol volres;
    end;

    if ~exist(fname_eigvect_atl,'file') | parms.forceflag
      % rotate eigen vectors
      fprintf('%s: rotating eigen vectors to atlas...\n',mfilename);
      volV = zeros([DTfit.volsz(1:3),9]);
      for i=1:3
        %% todo: verify that this is correct:
        tmp_volV = dti_rotate_vectors(squeeze(DTmeas.volV(:,:,:,i,:)),...
          inv(M_Atlas_to_DTI));
        if isempty(tmp_volV)
          error('failed to rotate eigen vectors');
        end;
        volV(:,:,:,1+3*(i-1):3+3*(i-1)) = tmp_volV;
        clear tmp_volV;
      end;
      % warp eigen vectors to atlas
      fprintf('%s: warping eigen vectors to atlas...\n',mfilename);
      volV_atl = zeros([volsz_T1(1:3),9]);
      for i=1:9
        vol = squeeze(volV(:,:,:,i));
        [vol,volres,M_atl] = applywarp(vol,M_DTI,vol_T1,M_T1_to_EPI,0,0,0,regStruct); 
        %interpm = 0; bclamp = 0;padding = 0
        volV_atl(:,:,:,i) = vol;
      end;
      if parms.sparse_flag
        mmil_save_sparse(volV_atl,fname_eigvect_atl,M_atl);
      else
        fs_save_mgh(volV_atl,fname_eigvect_atl,M_atl);
      end;
      clear volV_atl volV vol volres;
    end;

    if ~exist(fname_eigval_atl,'file') | parms.forceflag
      % warp eigen values to atlas
      fprintf('%s: warping eigen values to atlas...\n',mfilename);
      volE_atl = zeros([volsz_T1(1:3),3]);
      for i=1:3
        vol = squeeze(DTmeas.volE(:,:,:,i));
        [vol,volres,M_atl] = applywarp(vol,M_DTI,vol_T1,M_T1_to_EPI,1,0,1,regStruct); 
        %interpm = 1; bclamp = 0;padding = 1
        volE_atl(:,:,:,i) = vol;
      end;
      if parms.sparse_flag
        mmil_save_sparse(volE_atl,fname_eigval_atl,M_atl);
      else
        fs_save_mgh(volE_atl,fname_eigval_atl,M_atl);
      end;
      clear volE_atl vol volres;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% warp DT measures to atlas
if parms.warp_DTmeas_flag
  mmil_mkdir(dir_DTmeas_atlas);
  [tmp_path,tmp_DT_fstem] = fileparts(DT_fstem);
  nmeas = length(parms.DTI_measlist);
  outfiles = cell(nmeas,1);needfiles = 0;
  for i=1:nmeas
    meas = parms.DTI_measlist{i};
    outfiles{i} = sprintf('%s/%s_%s_atlas%s',...
      dir_DTmeas_atlas,tmp_DT_fstem,meas,parms.outext);
    if ~exist(outfiles{i},'file'), needfiles = 1; end;
  end;
  if needfiles || parms.forceflag
    for i=1:nmeas
      meas = parms.DTI_measlist{i};
      fname_out = outfiles{i};
      if ~exist(fname_out,'file') || parms.forceflag
        if exist('DTmeas','var') && ismember(meas,{'FA'})
          vol = getfield(DTmeas,sprintf('vol%s',meas));
        else % load from indiv volumes
          fname_in = sprintf('%s_%s%s',DT_fstem,meas,parms.ext);
          vol = fs_load_mgh(fname_in);          
        end;
        % warp meas values to atlas
        fprintf('%s: warping %s to atlas...\n',mfilename,meas);
        [vol_atl,volres,M_atl] = applywarp(vol,M_DTI,vol_T1,M_T1_to_EPI,1,1,1,regStruct); 
        %interpm = 1; bclamp = 1;padding = 1
        if parms.sparse_flag
          mmil_save_sparse(vol_atl,fname_out,M_atl);
        else
          fs_save_mgh(vol_atl,fname_out,M_atl);
        end;
        clear vol vol_atl volres;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% warp fiber masks or count masks and tensors to atlas
if parms.warp_fibers_flag
  mmil_mkdir(dir_fiber_atlas);
  mmil_mkdir(dir_fiber_resT1);
  if parms.atlasfibers_flag
    ext = '.mat';
    fiber_infix = parms.fiber_maps_infix;
  else
    if parms.fiber_countflag
      ext = '.mat';
      fiber_infix = 'count';
    else
      ext = '.mat';
      fiber_infix = 'mask';
    end;
  end;
  fiber_list = dir(sprintf('%s/fiber_*_%s%s',dir_fiber,fiber_infix,ext));

  if isempty(fiber_list)
    fprintf('%s: WARNING: no fiber files with infix %s found in %s\n',...
      mfilename,fiber_infix,dir_fiber);
  end;

  if exist(dir_fiber,'dir') & ~isempty(fiber_list)
    for f=parms.DTI_fibers
      vol_fiber = [];
      fname_stem_list = {...
          sprintf('fiber_%03d_%s',f,fiber_infix),...
        };
      for i=1:length(fname_stem_list)
        fname_stem = fname_stem_list{i};
        fname_fiber = sprintf('%s/%s%s',dir_fiber,fname_stem,ext);
        fname_res = sprintf('%s/%s_resT1%s',...
          dir_fiber_resT1,fname_stem,parms.outext);
        fname_atl = sprintf('%s/%s_atlas%s',...
          dir_fiber_atlas,fname_stem,parms.outext);
        if ~exist(fname_fiber,'file'), continue; end;
        if ~exist(fname_res,'file') | ~exist(fname_atl,'file') | parms.forceflag
          [vol_fiber,M_fiber] = mmil_load_sparse(fname_fiber);
          [vol_fiber,volres,M_atl] = applywarp(vol_fiber,M_fiber,vol_T1,M_T1_to_EPI,1,1,1,regStruct); 
          %interpm = 1; bclamp = 1;padding = 1
          if ~exist (fname_res,'file') | parms.forceflag
            fprintf('%s: resampling %s to T1...\n',mfilename,fname_fiber);
            [volres,Mres] = ctx_ctx2mgh(volres);
            if parms.sparse_flag
              mmil_save_sparse(volres,fname_res,Mres);
            else
              fs_save_mgh(volres,fname_res,Mres);
            end;
          end;
          if ~exist(fname_atl,'file') | parms.forceflag
            fprintf('%s: warping %s to atlas...\n',mfilename,fname_res);
            if parms.sparse_flag
              mmil_save_sparse(vol_fiber,fname_atl,M_atl);
            else
              fs_save_mgh(vol_fiber,fname_atl,M_atl);
            end;
          end;
          clear volres vol_fiber M_fiber;
        end;
      end;

      if parms.warp_fiber_tensors_flag
        fname_stem = sprintf('fiber_%03d',f);
        fname_fiber = sprintf('%s/%s_%s%s',...
          dir_fiber,fname_stem,fiber_infix,parms.outext);
        if ~exist(fname_fiber,'file'), continue; end;
        if parms.first_only_flag
          fname_tensor_sm = sprintf('%s/%s_tensor_smoothV0%s',...
            dir_fiber_atlas,fname_stem,parms.outext);
        else
          fname_tensor_sm = sprintf('%s/%s_tensor_smooth%s',...
            dir_fiber_atlas,fname_stem,parms.outext);
        end;
        % smooth tensors within fiber
        clear volT;
        if ~exist(fname_tensor_sm,'file') | parms.forceflag
          fprintf('%s: smoothing fiber %d tensors...\n',mfilename,f);
          % weights for tensor smoothing
          numer = (DTmeas.volE(:,:,:,1)-DTmeas.volE(:,:,:,2));
          denom = DTmeas.volE(:,:,:,1)+eps;
          % calculate DR
          volwt = numer./denom;
          % load fiber
          vol_fiber = mmil_load_sparse(fname_fiber);
          % smooth tensors within fiber
          ps = sqrt(sum(M_DTI(1:3,1:3).^2,1));
          sigmavec = parms.smfiber_sigma*ones(1,3)./ps;
          [vol_fiber,volT] = dti_smfiber_tensor(vol_fiber,DTmeas.volT,...
            volwt,sigmavec,parms.smfiber_win,[],[],parms.first_only_flag);
          % note: DTI_smfiber_tensor will reshape output into
          %  [nx,ny,nz,9] if input is [nx,ny,nz,3,3]
          if isempty(volT)
            error('fiber tensor smoothing failed for %s',fname_fiber);
          end;
          if parms.sparse_flag
            mmil_save_sparse(volT,fname_tensor_sm,M_DTI);
          else
            fs_save_mgh(volT,fname_tensor_sm,M_DTI);
          end;
          clear vol_fiber volT;
        end;
      
        if parms.first_only_flag
          fname_tensor_atl = sprintf('%s/%s_tensorV0_atlas%s',...
            dir_fiber_atlas,fname_stem,parms.outext);
        else
          fname_tensor_atl = sprintf('%s/%s_tensor_atlas%s',...
            dir_fiber_atlas,fname_stem,parms.outext);
        end;
        if ~exist(fname_tensor_atl,'file') | parms.forceflag
          if ~exist('volT','var') | isempty(volT)
            volT = mmil_load_sparse(fname_tensor_sm);
          end;
          % rotate tensors
          fprintf('%s: rotating fiber %d tensors to atlas...\n',mfilename,f);
          %% todo: verify that this is correct:
          volT = dti_rotate_tensors(volT,inv(M_Atlas_to_DTI));
          if isempty(volT)
            error('failed to rotate tensors');
          end;
          % warp tensors to atlas
          fprintf('%s: warping fiber %d tensors to atlas...\n',mfilename,f);
          volT_atl = zeros([volsz_T1(1:3),parms.num_tensor_components]);
          for i=1:parms.num_tensor_components
            t = parms.tensor_components(i);
            vol = squeeze(volT(:,:,:,t));
           [vol,volres,M_atl] = applywarp(vol,M_DTI,vol_T1,M_T1_to_EPI,1,0,1,regStruct);
            %interpm = 1;bclamp = 0; padding = 1
            volT_atl(:,:,:,i) = vol;
          end;
          if parms.sparse_flag
            mmil_save_sparse(volT_atl,fname_tensor_atl,M_atl);
          else
            fs_save_mgh(volT_atl,fname_tensor_atl,M_atl);
          end;
          clear volT volT_atl vol volres;
        end;
      end;
    end;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volwarp,volres,Mwarp] = applywarp(volin,M_in,vol_T1in,M_T12DTI,interpm,bclamp,padding,regStructin)
  volin = ctx_mgh2ctx(volin,M_in);
  volres = vol_resample_pad(volin,vol_T1in,M_T12DTI,interpm,bclamp);
  volwarp = volMorph(regStructin.volm, volres,...
                regStructin.VL, regStructin.VP, regStructin.VH,...
                interpm, padding, bclamp);
  volwarp.imgs(find(isnan(volwarp.imgs))) = 0;
  [volwarp,Mwarp] = ctx_ctx2mgh(volwarp);
return;
