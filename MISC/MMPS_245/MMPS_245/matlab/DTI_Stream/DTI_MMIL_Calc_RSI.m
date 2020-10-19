function errcode=DTI_MMIL_Calc_RSI(ContainerPath,varargin)
%function errcode=DTI_MMIL_Calc_RSI(ContainerPath,[options])
%
% Usage:
%  errcode = DTI_MMIL_Calc_RSI(ContainerPath,'key1', value1,...);
%
% Required Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%
% Optional Parameters:
%   'outdir': output directory
%     absolute or relative to ContainerPath
%     {default = 'RSIcalc'}
%   'outfix': string attached to output file names
%     {default = []}
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'regT1flag': [0|1|2] whether to resample RSI calculations into T1-alignment
%     0: do not resample RSI results to T1 resolution
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
%     for RSI calculations
%     {default = 6}
%   'min_bval': minimum b-value allowed for RSI calculations
%     {default = 0}
%   'flex_flag': [0|1] DTI_flex scans included in RSI calculations
%     {default = 0}
%   'min_nb0': minimum number of b=0 images required for RSI calculations
%     {default = 1}
%   'err_flag': calculate RSI fit error
%     {default = 1}
%   'RSImeas_flag': calculate RSI measures
%     {default = 1}
%   'measlist': cell array of DT measures to convert to mgh/mgz
%     {default = {'T','F0','N0','F2','N2','F4','N4','FD','ND',...
%                 'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'}}
%   'mask_RSImeas_flag': whether or not to mask RSImeas output volumes
%     {default = 1}
%   'dat_flag': convert RSI measures to dat format (for DTI Studio)
%     {default = 0}
%   'dat_measlist': cell array of DT measures to convert to dat
%     Must be included in 'measlist'
%     {default = {'F2','V0'}
%   'nii_flag': convert RSI measures to nii format (for FSL)
%     {default = 0}
%   'nii_measlist': cell array of DT measures to convert to nii
%     must be included in 'measlist'
%     {default = {'T','F0','N0','F2','N2','F4','N4','FD','ND',...
%                 'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'}}
%   'nii_out_orient': [0 1] output orientation for nii format
%     if empty, maintain original orientation
%    {default = []}
%   'scalefacts_flag': [0|1] calculate scaling factors from b=0 images
%     and apply them to all subsequent frames (for multiple acquisitions)
%     {default = 0}
%   'forceflag': run calculations even if output files exist
%     {default = 0}
%
% Optional RSI Input Parameters:
%   'lambda': regularization constant
%     {default = 0.1}
%   'iso_restricted_flag': [0|1] model isotropic diffusion of restricted water
%     {default = 1}
%   'iso_hindered_flag': [0|1] model isotropic diffusion of hindered water
%     {default = 1}
%   'iso_free_flag': [0|1] model isotropic diffusion of free water
%     {default = 1}
%   'ADC_hindered': ADC of isotropic hindered water (e.g. edema)
%     {default = 1.5e-3}
%   'ADC_free': apparent diffusion coefficient (ADC) of
%               isotropic free water (e.g. CSF)
%     {default = 3e-3}
%   'ADC_long': longitudinal ADC
%     {default = 1e-3}
%   'ADC_trans_min': minimum transverse ADC
%     {default = 0}
%   'ADC_trans_max': maximum transverse ADC
%     {default = 0.9e-3}
%   'num_ADC_trans': number of transverse ADC size scales
%     {default = 5}
%   'SH_order': spherical harmonic order -- must be even
%     {default = 4}
%   'norm_flag': normalize data to b=0 image
%     {default = 0}
%
% Created:  02/05/13 by Don Hagler
% Prev Mod: 07/31/15 by Don Hagler
% Last Mod: 11/03/17 by Don Hagler
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
parms = mmil_args2parms(varargin,{...
  'outdir','RSIcalc',[],...
  'outfix',[],[],...
  'snums',[],[],...
  'regT1flag',0,[0 1 2],...
  'infix',[],[],...
  'revflag',0,{0,1,2},...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
  'min_nb0',1,[],...
  'nob0_flag',false,[false true],...
  'err_flag',true,[false true],...
  'RSImeas_flag',true,[false true],...
  'measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
              'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
  'mask_RSImeas_flag',true,[false true],...
  'dat_flag',false,[false true],...
  'dat_measlist',{'F2','V0'},[],...
  'nii_flag',false,[false true],...
  'nii_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                  'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
  'nii_out_orient',[],[],...
  'scalefacts_flag',false,[false true],...
  'forceflag',false,[false true],...
... % RSI parameters
  'lambda',0.1,[],...
  'iso_free_flag',true,[false true],...
  'iso_hindered_flag',true,[false true],...
  'iso_restricted_flag',true,[false true],...
  'ADC_free',3e-3,[],...
  'ADC_hindered',1.5e-3,[],...
  'ADC_long',1e-3,[],...
  'ADC_trans_min',0,[],...
  'ADC_trans_max',0.9e-3,[],...
  'num_ADC_trans',5,[],...
  'SH_order',4,[2:2:10],...
  'norm_flag',false,[false true],...
... % hidden parameters
  'fnamestem','DTI',[],...
...
  'fstem_tags',{'snums','infix','revflag','min_bval','flex_flag',...
                'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
  'info_tags',{'snums','revflag','min_nb0','min_ndirs',...
               'min_bval','flex_flag','fnamestem'},[],...
  'data_tags',{'snums','infix','revflag','min_bval','flex_flag',...
               'min_ndirs','min_nb0'},[],...
  'reg_tags',{'infix','revflag'},[],...
  'rsi_tags',{'lambda','iso_free_flag','iso_hindered_flag',...
              'iso_restricted_flag','ADC_free','ADC_hindered',...
              'ADC_long','ADC_trans_min','ADC_trans_max','num_ADC_trans',...
              'SH_order','norm_flag','nob0_flag','scalefacts_flag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms, set dependent parms

args = mmil_parms2args(parms,parms.fstem_tags);
[fstem,parms.snums] = DTI_MMIL_Set_RSI_fstem(ContainerPath,args{:});
if isempty(fstem)
  errcode = 1;
  return;
end;
mmil_mkdir(fileparts(fstem));

if parms.nob0_flag
  % check that we have multiple b values
  args = mmil_parms2args(parms,parms.info_tags);
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

% remove Ir, Ih, and If from measlist
%   depending on iso_restricted_flag, iso_hindered_flag, and iso_free_flag
if ~parms.iso_restricted_flag
  parms.measlist = setdiff(parms.measlist,{'Ir'});
end;
if ~parms.iso_hindered_flag
  parms.measlist = setdiff(parms.measlist,{'Ih'});
end;
if ~parms.iso_free_flag
  parms.measlist = setdiff(parms.measlist,{'If'});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate RSI fit

RSIfit = [];
vol = [];
fname_fit = [fstem '_fit.mat'];
if ~exist(fname_fit,'file') || parms.forceflag
  args = mmil_parms2args(parms,parms.data_tags);
  [vol,M,qmat,bvals]=DTI_MMIL_Load_Data(ContainerPath,args{:});
  if isempty(vol), return; end;

  nframes = size(qmat,1);
  fprintf('%s: running RSI calculations with %d frames...\n',...
    mfilename,nframes);
  tic
  args = mmil_parms2args(parms,parms.rsi_tags);
  RSIfit = rsi_fit_MFOD(vol,qmat,'bvals',bvals,'M',M,args{:});
  toc
  if isempty(RSIfit)
    fprintf('%s: ERROR: tensor calculations were unsuccesful\n',mfilename);
    errcode = 1;
    return;
  end;

  % save output structure in matfile
  save(fname_fit,'RSIfit','M'); % include M for saving mgh files
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate error

if parms.err_flag
  fname_err = [fstem '_err.mgz'];
  if ~exist(fname_err,'file') || parms.forceflag
    if isempty(vol)
      args = mmil_parms2args(parms,parms.data_tags);
      [vol,M,qmat,bvals]=DTI_MMIL_Load_Data(ContainerPath,args{:});
      if isempty(vol), return; end;
    end;
    if isempty(RSIfit)
      fprintf('%s: loading %s...\n',mfilename,fname_fit);
      load(fname_fit);
    end;
    % select frames used in RSI fit
    vol = vol(:,:,:,RSIfit.i_fit);

    fprintf('%s: calculating RSI error...\n',mfilename);
    % calculate synthesized volume
    tic
    vol_synth = rsi_synth_vol(RSIfit); % save memory by reusing variable later
    
    % calculate sum of squared errors and normalize
    % one slice at a time to save memory
    vol_ss = zeros(RSIfit.volsz(1:3));
    vol_err = zeros(RSIfit.volsz(1:3));
    for z = 1:size(vol_err,3)
      slc = vol(:,:,z,:);
      slc_synth = vol_synth(:,:,z,:);
      vol_err(:,:,z) = sum((slc_synth-slc).^2,4);
      vol_ss(:,:,z) = sum(slc.^2,4);
    end;
    vol_err = vol_err/max(vol_ss(:));

    % apply mask
    vol_err = vol_err.*RSIfit.volmask_dilated;

    % save err in mgz file
    fs_save_mgh(vol_err,fname_err,M);
    toc
  end;
end;
clear vol vol_synth vol_err vol_ss slc slc_synth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate F0, F2, etc.

if ~parms.RSImeas_flag, return; end;
RSImeas = [];
fname_meas = [fstem '_meas.mat'];
if ~exist(fname_meas,'file') || parms.forceflag
  if isempty(RSIfit)
    fprintf('%s: loading %s...\n',mfilename,fname_fit);
    load(fname_fit);
  end;

  fprintf('%s: calculating RSI measures...\n',mfilename);
  tic
  RSImeas = rsi_calc_MFmeas(RSIfit,'mask_flag',parms.mask_RSImeas_flag);
  toc

  % save output structure in matfile
  save(fname_meas,'RSImeas','M'); % include M for saving mgh files
end;
clear RSIfit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert RSI measures to mgz, dat, or nii

if ~parms.forceflag
  flist = []; j=1;
  for i=1:length(parms.measlist)
    flist{j} = [fstem '_' parms.measlist{i} '.mgz'];
    if strcmp(parms.measlist{i}, 'ND')
      fname_ND = flist{j}; 
    end;
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
  if isempty(RSImeas)
    fprintf('%s: loading %s...\n',mfilename,fname_meas);
    load(fname_meas);
  end;

  if parms.regT1flag
    args = mmil_parms2args(parms,parms.reg_tags);
    [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,args{:});
    if errcode
      fprintf('%s: WARNING: RegInfo not found in %s, cannot resample RSImeas to T1\n',...
        mfilename,ContainerPath);
      parms.regT1flag = 0;
    end;
  else
    RegInfo = [];
  end;

  rsi_convert_MFmeas(RSImeas,fstem,'M',M,'RegInfo',RegInfo,...
    'regT1flag',parms.regT1flag,...
    'measlist',parms.measlist,...
    'dat_flag',parms.dat_flag,...
    'dat_measlist',parms.dat_measlist,...
    'nii_flag',parms.nii_flag,...
    'nii_measlist',parms.nii_measlist,...
    'nii_out_orient',parms.nii_out_orient,...
    'forceflag',parms.forceflag);
end;

args = mmil_parms2args(parms,parms.reg_tags);
[RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(ContainerPath,args{:});
if ~errcode
  mmil_write_regscripts(fname_reg,'fname_T2',fname_ND);
end;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


