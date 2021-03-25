function errcode=DTI_MMIL_Calc_DT(ContainerPath,varargin)
%function errcode=DTI_MMIL_Calc_DT(ContainerPath,[options])
%
% Usage:
%  errcode = DTI_MMIL_Calc_DT(ContainerPath,'key1', value1,...);
%
% Required Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%
% Optional Parameters:
%   'outdir': output directory
%     absolute or relative to ContainerPath
%     {default = 'DTcalc'}
%   'outfix': string attached to output file names
%     {default = []}
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'regT1flag': [0|1|2] whether to resample DT calculations into T1-alignment
%     0: do not resample DT results to T1 resolution
%     1: apply DTI to T1 registration but keep DTI resolution,
%        (padded with extra slices to match number of in-plane voxels)
%     2: apply DTI to T1 registration and resample to T1 resolution
%     {default = 0}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr_resDTI'
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
%     {default = 0}
%   'min_ndirs': minimum number of gradient directions allowed
%     for tensor calculations
%     {default = 6}
%   'min_bval': minimum b-value allowed for tensor calculations
%     {default = 0}
%   'max_bval': maximum b-value used in tensor fit
%     {default = Inf}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': minimum number of b=0 images required for tensor calculations
%     {default = 1}
%   'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit
%     {default = 0}
%   'err_flag': calculate tensor fit error
%     {default = 1}
%   'DTmeas_flag': calculate DT measures (Eigenvalues, Eigvenvectors, FA, etc.)
%     {default = 1}
%   'mask_DTmeas_flag': whether or not to mask DTmeas output volumes
%     {default = 1}
%   'dat_flag': convert FA and V0 to dat format (for DTI Studio)
%     {default = 0}
%   'nii_flag': convert DT measures to nii format (for FSL)
%     {default = 0}
%   'nii_out_orient': [0 1] output orientation for nii format
%     if empty, maintain original orientation
%     {default = []}
%   'scalefacts_flag': [0|1] calculate scaling factors from b=0 images
%     and apply them to all subsequent frames (for multiple acquisitions)
%     {default = 0}
%   'forceflag': run calculations even if output files exist
%     {default = 0}
%
%
% Created:  12/11/06 by Don Hagler
% Last Mod: 06/23/17 by Don Hagler
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
errcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set constants
smf = 10^-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir','DTcalc',[],...
  'outfix',[],[],...
  'snums',[],[],...
  'regT1flag',0,[0 1 2],...
  'infix',[],[],...
  'revflag',0,{0,1,2},...
  'nob0_flag',false,[false true],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'max_bval',Inf,[100,Inf],...
  'flex_flag',false,[false true],...
  'min_nb0',1,[],...
  'nonlin_flag',false,[false true],...
  'err_flag',true,[false true],...
  'DTmeas_flag',true,[false true],...
  'mask_DTmeas_flag',true,[false true],...
  'dat_flag',false,[false true],...
  'nii_flag',false,[false true],...
  'nii_out_orient',[],[],...
  'scalefacts_flag',false,[false true],...
  'forceflag',false,[false true],...
... % hidden parameters
  'measlist',{'b0','FA','V0','MD','LD','TD'},[],...
  'dat_measlist',{'FA','V0'},[],...
  'nii_measlist',{'b0','FA','MD','LD','TD'},[],...
  'fnamestem','DTI',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms, set dependent parms

tags = {'snums','infix','revflag','min_bval','max_bval','flex_flag',...
  'min_ndirs','min_nb0','nob0_flag','outdir','outfix'};
args = mmil_parms2args(parms,tags);
[fstem,parms.snums] = DTI_MMIL_Set_DT_fstem(ContainerPath,args{:});
if isempty(fstem)
  errcode = 1;
  return;
end;
mmil_mkdir(fileparts(fstem));

if parms.nob0_flag
  % check that we have multiple b values
  tags = {'snums','revflag','min_nb0','min_ndirs',...
    'min_bval','flex_flag','fnamestem'};
  args = mmil_parms2args(parms,tags);
  [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:})
  if errcode
    errcode = 1;
    return;
  end;
  tmp_bvals = [];
  for i=1:length(parms.snums)
    snum = parms.snums(i);
    bval = ScanInfo(snum).bval;
    tmp_bvals = [tmp_bvals;bval];
  end;
  if length(unique(tmp_bvals))<2
    fprintf('%s: ERROR: multiple b values are required if nob0_flag=1\n',mfilename);
    errcode = 1;
    return;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate tensor fit

DTfit = [];
vol = [];
fname_fit = [fstem '_fit.mat'];
if ~exist(fname_fit,'file') || parms.forceflag
  tags = {'snums','infix','revflag','min_bval','flex_flag',...
    'min_ndirs','min_nb0'};
  args = mmil_parms2args(parms,tags);
  [vol,M,qmat,bvals]=DTI_MMIL_Load_Data(ContainerPath,args{:});
  if isempty(vol), return; end;

  nframes = size(qmat,1);
  fprintf('%s: running tensor calculations with %d frames...\n',...
    mfilename,nframes);
  tic
  DTfit = dti_fit_tensor(vol,qmat,...
    'bvals',bvals,'M',M,...
    'nob0_flag',parms.nob0_flag,...
    'max_bval',parms.max_bval,...
    'nonlin_flag',parms.nonlin_flag,...
    'scalefacts_flag',parms.scalefacts_flag);
  toc
  if isempty(DTfit)
    fprintf('%s: ERROR: tensor calculations were unsuccesful\n',mfilename);
    errcode = 1;
    return;
  end;

  % save output structure in matfile
  save(fname_fit,'DTfit','M'); % include M for saving mgh files
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate error

if parms.err_flag
  fname_err = [fstem '_err.mgz'];
  if ~exist(fname_err,'file') || parms.forceflag
    if isempty(DTfit)
      fprintf('%s: loading %s...\n',mfilename,fname_fit);
      load(fname_fit);
    end;
    if isempty(vol)
      tags = {'snums','infix','revflag','min_bval','flex_flag',...
        'min_ndirs','min_nb0'};
      args = mmil_parms2args(parms,tags);
      [vol,M,qmat,bvals]=DTI_MMIL_Load_Data(ContainerPath,args{:});
      if isempty(vol), return; end;
    end;
    % select frames used in tensor fit
    vol = vol(:,:,:,DTfit.i_fit);

    fprintf('%s: calculating DT error...\n',mfilename);
    % calculate synthesized volume
    tic
    vol_synth = dti_synth_vol(DTfit); % save memory by reusing variable later
    
    % calculate sum of squared errors and normalize
    % one slice at a time to save memory
    vol_ss = zeros(DTfit.volsz(1:3));
    vol_err = zeros(DTfit.volsz(1:3));
    for z = 1:size(vol_err,3)
      slc = vol(:,:,z,:);
      slc_synth = vol_synth(:,:,z,:);
      vol_err(:,:,z) = sum((slc_synth-slc).^2,4);
      vol_ss(:,:,z) = sum(slc.^2,4);
    end;
    vol_err = vol_err/max(vol_ss(:));

    % apply mask
    vol_err = vol_err.*DTfit.volmask_dilated;

    % save err in mgz file
    fs_save_mgh(vol_err,fname_err,M);
    toc
  end;
end;
clear vol vol_synth vol_err vol_ss slc slc_synth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate Eigenvectors, FA, etc.

if ~parms.DTmeas_flag, return; end;
DTmeas = [];
fname_meas = [fstem '_meas.mat'];
if ~exist(fname_meas,'file') || parms.forceflag
  if isempty(DTfit)
    fprintf('%s: loading %s...\n',mfilename,fname_fit);
    load(fname_fit);
  end;

  fprintf('%s: calculating DT measures...\n',mfilename);
  tic
  DTmeas = dti_calc_DTmeas(DTfit,parms.mask_DTmeas_flag);
  toc

  % save output structure in matfile
  save(fname_meas,'DTmeas','M'); % include M for saving mgh files
end;
clear DTfit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert DT measures to mgz, dat, or nii

if ~parms.forceflag
  flist = []; j=1;
  for i=1:length(parms.measlist)
    flist{j} = [fstem '_' parms.measlist{i} '.mgz'];
    j=j+1;
  end;
  if parms.dat_flag
    for i=1:length(parms.dat_measlist)
      flist{j} = [fstem '_' parms.dat_measlist{i} '.dat'];
      j=j+1;
    end;
  end;
  if parms.nii_flag
    for i=1:length(parms.nii_measlist)
      flist{j} = [fstem '_' parms.nii_measlist{i} '.nii'];
      j=j+1;
    end;
  end;
  convert_flag = 0;
  for j=1:length(flist)
    if ~exist(flist{j},'file'), convert_flag = 1; end;
  end;

  if parms.regT1flag && ~convert_flag
    if parms.regT1flag==2
      suffix = 'resT1';
    else
      suffix = 'regT1';
    end;

    flist = []; j=1;
    for i=1:length(parms.measlist)
      flist{j} = [fstem '_' parms.measlist{i} '_' suffix '.mgz'];
      j=j+1;
    end;
    if parms.dat_flag
      for i=1:length(parms.dat_measlist)
        flist{j} = [fstem '_' parms.dat_measlist{i} '_' suffix '.dat'];
        j=j+1;
      end;
    end;
    if parms.nii_flag
      for i=1:length(parms.nii_measlist)
        flist{j} = [fstem '_' parms.nii_measlist{i} '_' suffix '.nii'];
        j=j+1;
      end;
    end;
    convert_flag = 0;
    for j=1:length(flist)
      if ~exist(flist{j},'file'), convert_flag = 1; end;
    end;
  end;
else
  convert_flag = 1;
end;

if convert_flag
  if isempty(DTmeas)
    fprintf('%s: loading %s...\n',mfilename,fname_meas);
    load(fname_meas);
  end;

  if parms.regT1flag
    tags = {'infix','revflag'};
    args = mmil_parms2args(parms,tags);
    [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,args{:});
    if errcode
      fprintf('%s: WARNING: RegInfo not found in %s, cannot resample DTmeas to T1\n',...
        mfilename,ContainerPath);
      parms.regT1flag = 0;
    end;
  else
    RegInfo = [];
  end;

  dti_convert_DTmeas(DTmeas,fstem,'M',M,'RegInfo',RegInfo,...
    'regT1flag',parms.regT1flag,...
    'measlist',parms.measlist,...
    'dat_flag',parms.dat_flag,...
    'dat_measlist',parms.dat_measlist,...
    'nii_flag',parms.nii_flag,...
    'nii_measlist',parms.nii_measlist,...
    'nii_out_orient',parms.nii_out_orient,...
    'forceflag',parms.forceflag);
end;

tags = {'infix','revflag'};
args = mmil_parms2args(parms,tags);
[RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,args{:});
if ~errcode
  fname_FA = [fstem '_FA.mgz'];
  if exist(fname_FA,'file')
    mmil_write_regscripts(fname_reg,'fname_T2',fname_FA);
  else
    fprintf('%s: WARNING: file %s not found\n',mfilename,fname_FA);
  end;
end;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


