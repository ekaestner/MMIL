function dti_convert_DTmeas(DTmeas,fstem,varargin)
%function dti_convert_DTmeas(DTmeas,fstem,[options])
%
% Usage:
%  dti_convert_DTmeas(DTmeas,fstem,'key1', value1,...);
%
% Required Parameters:
%   DTmeas: output structure with fields:
%     volT : tensor matrix                                size = [nx,ny,nz,3,3]
%     volV : eigen vectors                                size = [nx,ny,nz,3,3]
%     volE : eigen values                                 size = [nx,ny,nz,3]
%     volFA : fractional anisotropy                       size = [nx,ny,nz]
%   fstem: output file stem (include full path)
%
% Optional Input Parameters:
%   'regT1flag': [0|1|2] whether to resample DT calculations into T1-alignment
%     This is done in addition to converting in DTI resolution
%     0: do not resample DT results to T1 resolution
%     1: apply DTI to T1 registration but keep DTI resolution,
%        (padded with extra slices to match number of in-plane voxels)
%     2: apply DTI to T1 registration and resample to T1 resolution
%     {default = 0}
%   'zpadflag': [0|1] pad slices to match number of in-plane voxels
%     avoids clipping of image after registration to T1
%     (only applicable if regT1flag = 1)
%     {default = 1}
%   'RegInfo': structure with fields (Required if regT1flag > 0):
%     M_T1_to_T2 : transformation matrix
%     M_T2 : DTI vox2ras matrix
%     M_T1 : T1 vox2ras matrix
%     volsz_T2 : vector of DTI voxel sizes [x,y,z]
%     volsz_T1 : vector of T1 voxel sizes [x,y,z]
%     {default = []}
%   'M': vox2ras matrix of input diffusion data
%     ignored if RegInfo is supplied
%     {default = eye(4)}
%   'measlist': cell array of DT measures to convert to mgh/mgz
%     {default = {'b0','FA','V0','MD','LD','TD'}
%   'mgz_flag': [0|1] use compressed mgz format (otherwise mgh)
%     {default = 1}
%   'dat_flag': convert DT measures to dat format (for DTI Studio)
%     {default = 0}
%   'dat_measlist': cell array of DT measures to convert to dat
%     Must be included in 'measlist'
%     {default = {'FA','V0'}
%   'nii_flag': convert DT measures to nii format (for FSL)
%     {default = 0}
%   'nii_measlist': cell array of DT measures to convert to nii
%     Must be included in 'measlist'
%     {default = {'b0','FA','MD','LD','TD'}
%   'nii_out_orient': output orientation for nii files
%     e.g. 'RAS', 'LPI', etc.
%     If empty, keep original orientation
%     {default = []}
%   'verbose_flag': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/16/10 by Don Hagler
% Last Mod: 11/01/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'regT1flag',0,[0:2],...
  'zpadflag',true,[false true],...
  'M',eye(4),[],...
  'RegInfo',[],[],...
  'measlist',{'b0','FA','V0','MD','LD','TD'},[],...
  'mgz_flag',true,[false true],...
  'dat_flag',false,[false true],...
  'dat_measlist',{'FA','V0'},[],...
  'nii_flag',false,[false true],...
  'nii_measlist',{'b0','FA','MD','LD','TD'},[],...
  'nii_out_orient',[],[],...
  'verbose_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'dat_revsliceflag',true,[false true],... %% todo: should depend on orient of M
  'dat_permvec',[1,2,3],[],...  
});

if parms.mgz_flag
  parms.mgh_ext = '.mgz';
else
  parms.mgh_ext = '.mgh';
end;
if ~parms.dat_flag, parms.dat_measlist = []; end;
if ~parms.nii_flag, parms.nii_measlist = []; end;

if parms.regT1flag && isempty(parms.RegInfo)
  error('missing RegInfo with regT1flag = %d',parms.regT1flag);
end;

switch parms.regT1flag
  case 1
    parms.res_dat_permvec = parms.dat_permvec;
    outfix = ['regT1'];
    nx_res = size(DTmeas.volb0,1);
    ny_res = size(DTmeas.volb0,2);
    if parms.zpadflag
      % pad with extra slices to avoid clipping
      nz_res = size(DTmeas.volb0,1); %% todo: number of output slices as option?
    else
      nz_res = size(DTmeas.volb0,3);
    end;
    % adjust offset of output M so that it has same offset as T1
    M_out = parms.RegInfo.M_T2;
    nvox = [nx_res,ny_res,nz_res];
    M_out(1:3,4) = M_out(1:3,4) - parms.M(1:3,:)*[nvox/2+1 1]' +...
      parms.RegInfo.M_T1(1:3,:)*[parms.RegInfo.volsz_T1/2+1 1]';
    % create reference volume with correct output resolution, orientation, and offset
    vol_ref = ctx_mgh2ctx(zeros(nx_res,ny_res,nz_res),M_out);
  case 2
    parms.res_dat_permvec = [1,3,2];
    outfix = ['resT1'];
    vol_ref = ctx_mgh2ctx(zeros(parms.RegInfo.volsz_T1),parms.RegInfo.M_T1);
    nx_res = parms.RegInfo.volsz_T1(1);
    ny_res = parms.RegInfo.volsz_T1(2);
    nz_res = parms.RegInfo.volsz_T1(3);
end;
if ~isempty(parms.RegInfo)
  parms.M = parms.RegInfo.M_T2;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save measures in mgh,dat, nii formats
for m=1:length(parms.measlist)
  meas = parms.measlist{m};
  convert_flag = 0;
  fname_out = [fstem '_' meas parms.mgh_ext];
  fname_dat = [fstem '_' meas '.dat'];
  fname_nii = [fstem '_' meas '.nii'];
  if ~exist(fname_out,'file') ||...
    (~exist(fname_dat,'file') && any(ismember(meas,parms.dat_measlist))) ||...
    (~exist(fname_nii,'file') && any(ismember(meas,parms.nii_measlist))) ||...
    parms.forceflag
    convert_flag = 1;
  end;
  if parms.regT1flag
    fname_out_res = [fstem '_' meas '_' outfix parms.mgh_ext];
    fname_dat_res = [fstem '_' meas '_' outfix '.dat'];
    fname_nii_res = [fstem '_' meas '_' outfix '.nii'];
    if ~exist(fname_out_res,'file') ||...
      (~exist(fname_dat_res,'file') && any(ismember(meas,parms.dat_measlist))) ||...
      (~exist(fname_nii_res,'file') && any(ismember(meas,parms.nii_measlist))) ||...
      parms.forceflag
      convert_flag = 1;
    end;
  end;

  if convert_flag
    if parms.verbose_flag
      fprintf('%s: converting DT %s...\n',mfilename,meas);
    end;
    switch meas
      case {'b0','FA'}
        vol = getfield(DTmeas,['vol' meas]);
      case 'V0'
        vol = squeeze(DTmeas.volV(:,:,:,1,:));
      case 'MD'
        vol = mean(DTmeas.volE,4);
      case 'LD'
        vol = squeeze(DTmeas.volE(:,:,:,1));
      case 'TD'
        vol = mean(DTmeas.volE(:,:,:,2:3),4);
      case {'DR','DR1'} % dominance ratio 1
        numer = DTmeas.volE(:,:,:,1) - DTmeas.volE(:,:,:,2);
        denom = DTmeas.volE(:,:,:,1) + eps;
        vol = numer./denom;
      case 'DR2' % dominance ratio 2
        numer = DTmeas.volE(:,:,:,2) - DTmeas.volE(:,:,:,3);
        denom = DTmeas.volE(:,:,:,2) + eps;
        vol = numer./denom;
      otherwise
        fprintf('%s: WARNING: DT measure %s not supported\n',mfilename,meas);
        continue;
    end;
    if ~exist(fname_out,'file') || parms.forceflag
      fs_save_mgh(vol,fname_out,parms.M);
    end;

    if any(ismember(meas,parms.dat_measlist)) &&...
      (~exist(fname_dat,'file') || parms.forceflag)
      if strcmp(meas,'V0')
        dti_mgh2dat(fname_out,[4,parms.dat_permvec],parms.dat_revsliceflag);
      else
        dti_mgh2dat(fname_out,parms.dat_permvec,parms.dat_revsliceflag);
      end;
    end;

    if any(ismember(meas,parms.nii_measlist))
      fs_mri_convert(fname_out,fname_nii,...
        'out_orient',parms.nii_out_orient,'forceflag',parms.forceflag);
    end;

    % save measures in mgh,dat, nii formats after resampling to T1 resolution
    if parms.regT1flag
      if strcmp(meas,'V0')
        % rotate vectors
        M_T2_to_T1 = inv(parms.RegInfo.M_T1_to_T2);
        vol = dti_rotate_vectors(vol,M_T2_to_T1);  
        vol_res = zeros(nx_res,ny_res,nz_res,3,'single');
        % resample each component
        for i=1:3
          vol_DTI = ctx_mgh2ctx(squeeze(vol(:,:,:,i)),parms.RegInfo.M_T2);
          vol_DTI_res = vol_resample_pad(vol_DTI,vol_ref,...
            parms.RegInfo.M_T1_to_T2,0,0);
          vol_res(:,:,:,i) = single(vol_DTI_res.imgs);
        end;
        M_out = M_LPH_TO_RAS*vol_ref.Mvxl2lph;
        if ~exist(fname_out_res,'file') || parms.forceflag
          fs_save_mgh(vol_res,fname_out_res,M_out);
        end;
        if ismember(meas,parms.dat_measlist)
          dti_mgh2dat(fname_out_res,[4,parms.res_dat_permvec],parms.dat_revsliceflag);
        end;
      else
        vol_DTI = ctx_mgh2ctx(vol,parms.RegInfo.M_T2);
        vol_DTI_res = vol_resample_pad(vol_DTI,vol_ref,parms.RegInfo.M_T1_to_T2,2);
        if ~exist(fname_out_res,'file') || parms.forceflag
          ctx_save_mgh(vol_DTI_res,fname_out_res);
        end;
        if ismember(meas,parms.dat_measlist)
          dti_mgh2dat(fname_out_res,parms.res_dat_permvec,parms.dat_revsliceflag);
        end;
      end;

      if ismember(meas,parms.nii_measlist)
        fs_mri_convert(fname_out_res,fname_nii_res,...
          'out_orient',parms.nii_out_orient,'forceflag',parms.forceflag);
      end;
    end;
  end;
end;

