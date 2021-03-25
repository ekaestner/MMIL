function dti_mask_fibers(fiber_dir,varargin)
%function dti_mask_fibers(fiber_dir,[options])
%
% Required Input:
%   fiber_dir: full path of directory containing fiber ROI files
%
% Optional Parameters:
%   'outdir': output directory
%     {default = fiber_dir}
%   'fibers': fiber numbers to mask
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'resT1flag': [0|1] whether to use fibers already resampled to T1 resolution
%     {default = 0}
%   'atlas_flag': whether to use atlas fibers and if so, what type
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'mask_suffix': suffix attached to output fiber file names
%     {default = 'masked'}
%   'save_mgz_flag': [0|1] save fibers in mgz format in addition to sparse
%     {default = 0}
%   'verbose': display status messages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Optional Parameters for aseg masking:
%   'fname_aseg': full path name of aseg file
%     aseg volume should have dimensions matching DTI data,
%       unless M_T1_to_DTI is supplied
%     if supplied, will be used to exclude voxels based on aseg_codes
%     {default = []}
%   'M_T1_to_DTI': registration matrix between DTI data and high-res T1 data
%     if not empty, will be used to resample aseg to DTI resolution
%     {default = []}
%   'aseg_codes': aseg ROI code numbers to be excluded or included
%     {default = [0,24,4,5,14,15,43,44,72,75,76,3,8,42,47] (CSF and gray-matter)
%   'exclude_flag': exclude voxels with aseg_codes
%     otherwise only include voxels with aseg_codes
%     {default = 1}
%
% Optional Parameters for dispersion weighting:
%  'fname_vals': full path of value file (with dimensions matching fibers)
%     required for dispersion weighting
%     {default = []}
%  'disp_scalefact': scaling factor applied to values in fname_vals
%    {default = 1}
%  'dispvec': vector of dispersion values (MAD estimates)
%    may be a single value or vector with size matching number of fibers
%    {default = 0.1}
%  'dispfact': multiplicative factor applied to dispersion values
%    {default = 4}
%  'disp_mask_flag': [0|1] interaction between disp weighting and aseg mask
%     0: calculate dispersion weighting everywhere, but mask out aseg_codes
%     1: calculate dispersion weighting only in aseg_codes voxels
%     {default = 0}
%
% Optional Parameters for thresholding:
%   'thresh_FA': fractional anisotropy threshold applied to fiber ROIs
%     {default = 0}
%   'fname_FA': full path name of FA (fractional anisotropy) image file
%     Required if thresh_FA>0
%     {default = []}
%   'thresh_prob': fiber probability threshold applied to fiber ROIs
%     {default = 0}
%
% Created:  10/10/12 by Don Hagler
% Last Mod: 11/21/13 by Don Hagler
%

% based on DTI_MMIL_AsegMask_Fibers, created 06/08/09 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters and check for problems
parms = check_input(fiber_dir,varargin);

% check if output files already exist
if ~parms.forceflag
  all_exist = check_output(parms);
  if all_exist, return; end;
end;

% resample aseg to DTI resolution if necessary
if parms.aseg_flag
  parms = resamp_aseg(parms);
end;

% load file containing values for dispersion weighting
if ~isempty(parms.fname_vals)
  parms = load_values(parms);
end;

% apply aseg mask, dispersion, thresh_prob, and/or thresh_FA to fibers
mask_fibers(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fiber_dir,options)
  parms_filter = {...
    'fiber_dir',fiber_dir,[],...
  ...
    'outdir',fiber_dir,[],...
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'resT1flag',false,[false true],...
    'atlas_flag',2,[0:4],...
    'mask_suffix','masked',[],...
    'save_mgz_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'fname_aseg',[],[],...
    'M_T1_to_DTI',[],[],...
    'aseg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
    'exclude_flag',true,[false true],...
  ...
    'fname_vals',[],[],...
    'dispvec',0.1,[],...
    'disp_scalefact',1,[1e-10,1e10],...
    'dispfact',4,[1e-6,1e6],...
    'disp_mask_flag',false,[false true],...
  ...
    'thresh_FA',0,[0 1],...
    'fname_FA',[],[],...
    'thresh_prob',0,[0 1],...
  ... % hidden parameters
    'count_flag',true,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % check if input exists
  if ~exist(fiber_dir,'dir')
    error('fiber_dir %s not found',fiber_dir);
  end;
  parms.nfibers = length(parms.fibers);
  % set fiber_infix and fiber_ext
  [parms.fiber_infix,parms.fiber_ext] = dti_set_fiber_infix(...
    'resT1flag',parms.resT1flag,'count_flag',parms.count_flag,...
    'atlas_flag',parms.atlas_flag);
  % check fiber files
  flist = dir(sprintf('%s/fiber_*_%s%s',...
    parms.fiber_dir,parms.fiber_infix,parms.fiber_ext));
  if isempty(flist)
    error('no %s fiber files with fiber infix %s found in %s\n',...
      parms.fiber_ext,parms.fiber_infix,parms.fiber_dir);
  end;
  if isempty(parms.fname_aseg) &&...
     isempty(parms.fname_vals) &&...
     parms.thresh_prob==0 && parms.thresh_FA==0
    fprintf('%s: WARNING: output fibers will be identical to input fibers\n',...
      mfilename);
  end;
  % check aseg file
  if ~isempty(parms.fname_aseg)
    parms.aseg_flag = 1;
    if ~exist(parms.fname_aseg,'file')
      error('file %s not found',parms.fname_aseg);
    end;
  else
    parms.aseg_flag = 0;
  end;
  % check dispersion parameters
  if ~isempty(parms.fname_vals)
    parms.disp_flag = 1;
    if ~exist(parms.fname_vals)
      error('file %s not found',parms.fname_vals);
    end;
    if length(parms.dispvec)==1
      parms.dispvec = parms.dispvec*ones(1,parms.nfibers);
    elseif length(parms.dispvec)~=parms.nfibers
      error('# of dispersion values (%d) does not match number of fibers (%d)',...
        length(parms.dispvec),parms.nfibers);    
    end;
  end;
  parms.dispvec = parms.dispfact * parms.dispvec;
  % check FA threshold
  if parms.thresh_FA>0
    if isempty(parms.fname_FA)
      error('must specify fname_FA if thresh_FA>0');
    end;
    if ~exist(parms.fname_FA,'file')
      error('file %s not found',parms.fname_FA);
    end;
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_exist = check_output(parms)
  all_exist = 1;
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_out = sprintf('%s/fiber_%d_%s_%s%s',...
      parms.outdir,fnum,parms.fiber_infix,parms.mask_suffix,parms.fiber_ext);
    if ~exist(fname_out,'file')
      all_exist = 0;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_aseg(parms)
  [volsz_DTI,M_DTI] = read_DTI_volsz(parms);
  [volsz_aseg,M_aseg] = read_volsz(parms.fname_aseg,parms);
  if all(volsz_aseg==volsz_DTI) && all(M_aseg(:)==M_DTI(:))
    % no resampling needed
    return;
  end;
  if isempty(parms.M_T1_to_DTI)
    error('aseg must be resampled to DTI resolution but M_T1_to_DTI not supplied');
  end;
  [tmp,fstem,fext] = fileparts(parms.fname_aseg);
  fname_aseg_resDTI = sprintf('%s/%s_resDTI%s',...
    parms.outdir,fstem,fext);
  if ~exist(fname_aseg_resDTI,'file') || parms.forceflag
    [vol,M] = fs_load_mgh(parms.fname_aseg);
    [vol,M] = mmil_resample_vol(vol,M,...
      'nvox_ref',volsz_DTI,'M_ref',M_DTI,...
      'interpm',0,...
      'M_reg',inv(parms.M_T1_to_DTI));
    fs_save_mgh(vol,fname_aseg_resDTI,M);
  end;
  parms.fname_aseg = fname_aseg_resDTI;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volsz_DTI,M_DTI] = read_DTI_volsz(parms)
  volsz_DTI = [];
  M_DTI = [];
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_fiber = sprintf('%s/fiber_%d_%s%s',...
      parms.fiber_dir,fnum,parms.fiber_infix,parms.fiber_ext);
    if ~exist(fname_fiber,'file')
      continue;
    end;  
    [volsz_DTI,M_DTI] = read_volsz(fname_fiber,parms);
    break;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_values(parms)
  [parms.vol_vals,M_vals,volsz_vals] = load_vol(parms.fname_vals);
  [volsz_DTI,M_DTI] = read_DTI_volsz(parms);
  if any(volsz_vals~=volsz_DTI) || any(M_vals(:)~=M_DTI(:))
    error('volume in fname_vals does not match fibers');
  end;
  if size(parms.vol_vals,4)>1
    error('multi-frame vals volume not currently supported');
  end;
  if parms.disp_scalefact ~= 1
    parms.vol_vals = parms.disp_scalefact * parms.vol_vals;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mask_fibers(parms)
  if parms.aseg_flag
    vol_aseg = load_vol(parms.fname_aseg);
    if parms.exclude_flag
      ind_exclude = find(ismember(vol_aseg(:),parms.aseg_codes));
    else
      ind_exclude = find(~ismember(vol_aseg(:),parms.aseg_codes));
    end;
  else
    ind_exclude = [];
  end;
  if parms.thresh_FA>0
    vol_FA = load_vol(parms.fname_FA);
  end;
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_in = sprintf('%s/fiber_%d_%s%s',...
      parms.fiber_dir,fnum,parms.fiber_infix,parms.fiber_ext);
    fname_out = sprintf('%s/fiber_%d_%s_%s%s',...
      parms.outdir,fnum,parms.fiber_infix,parms.mask_suffix,parms.fiber_ext);
    if exist(fname_in,'file')
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: masking fiber file %s...\n',mfilename,fname_in);
        end;
        [vol_fiber,M] = load_vol(fname_in);
        if ~isempty(parms.fname_vals)
          % dispersion weighting
          if parms.disp_mask_flag && parms.aseg_flag
            ind_mask = ind_exclude;
          else
            ind_mask = [];
          end;
          vol_fiber = calc_weights(vol_fiber,parms.vol_vals,...
            parms.dispvec(f),ind_mask);
        end;
        if parms.aseg_flag && (~parms.disp_mask_flag || ~parms.disp_flag)
          % aseg masking
          vol_fiber(ind_exclude) = 0;
        end;
        if parms.thresh_prob>0
          % probability threshold
          vol_fiber(vol_fiber < parms.thresh_prob) = 0;
        end;
        if parms.thresh_FA>0
          % FA threshold
          vol_fiber(vol_FA < parms.thresh_FA) = 0;
        end;
        save_vol(vol_fiber,fname_out,M);
      end;
      if parms.save_mgz_flag && strcmp(parms.fiber_ext,'.mat')
        fname_mgz = regexprep(fname_out,'.mat','.mgz');
        mmil_sparse2mgh(fname_out,fname_mgz,parms.forceflag);
      end;
    else      
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_fiber);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_fiber = calc_weights(vol_fiber,vol_vals,dispersion,ind_mask)
  % define Tukey's bisquare function
  weightfun = @(x,mu,nu) max(0,1-((x-mu)/nu).^2);
  % get values and weights within fiber ROI
  ind_roi = find(vol_fiber);
  x = vol_vals(ind_roi);
  w = vol_fiber(ind_roi);
  % calculate weighted median  
  mu = mmil_wtd_median(x,w,1);
  % get values in "excluded" voxels
  if ~isempty(ind_mask)
    ind_roi = intersect(ind_roi,ind_mask);
    x = vol_vals(ind_roi);
  end;
  % calulate weights
  nu = dispersion;
  w = weightfun(x,mu,nu); 
  vol_fiber(ind_roi) = w .* vol_fiber(ind_roi);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volsz,M] = read_volsz(fname,parms)
  volsz = []; M = [];
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [tmp,M,volsz] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);
  end;
  volsz = volsz(1:3);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,M,volsz] = load_vol(fname)
  vol = []; M = []; volsz = [];
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [vol,M,volsz] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [vol,M,tmp,volsz] = fs_load_mgh(fname);
  end;
  volsz = volsz(1:3);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_vol(vol,fname,M)
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      mmil_save_sparse(vol,fname,M);
    case {'.mgh','.mgz'}
      fs_save_mgh(vol,fname,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

