function errcode=mmil_average_volumes(fname_inlist,fname_mask,varargin)
%function errcode=mmil_average_volumes(fname_inlist,fname_mask,[options])
%
%  Purpose: Registering and averaging two or more MRI volumes
%             for a single subject, with the same contrast properties
%           e.g. averaging two T1-weighted images together
%
%  Required Input:
%    fname_inlist: cell array of input file names to be registered and averaged
%      input files must be mgh or mgz format
%    fname_mask: name of mgh file containing 3D binary mask registered to 
%      volume specified by first member of fname_inlist (or fname_ref)
%
%  Optional Input:
%    'fname_ref': name of mgh file containing reference volume
%      If empty, uses first file name in fname_inlist
%      If not empty and ref_scantype = input_scantype,
%        assumes this is a resampled version of the first file in fname_inlist
%      { default: [] }
%    'fname_out': output file name
%      should have mgh or mgz extension
%      { default: 'average_vol.mgh' }
%    'ref_scantype': scan type of fname_ref
%       Allowed: 'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo'
%      { default: 'MPR' }
%    'input_scantype': scan type of files in fname_inlist
%       Allowed: 'MPR','XetaT2','FLASHhi','FLASHlo','MEDIChi','MEDIClo'
%       If ref and input scantypes are different, will use jpdf registration
%         for first input file, then reg for subsequent scans
%      { default: 'MPR' }
%    'M_init': initial estimate of registration matrix for jpdf
%      e.g. registration matrix for reference (M_res_to_orig)
%      { default: identity matrix }
%    'regflag': [0|1] whether to register first scan to reference
%       or just assume they are already registered (after applying M_init)
%      { default: 1 }
%    'image_range': range of values to clip output image
%      if empty, no clipping
%      {default = []}
%    'cleanupflag': [0|1] whether to remove temporary files created by reg
%      { default: 1 }
%    'forceflag': [0|1] whether to run calculations even if output file exists
%      { default: 0 }
%
%  Output:
%   errcode: 0 if successful, 1 if error
%
% Created:  09/21/10 by Don Hagler
% Prev Mod: 05/11/17 by Don Hagler
% Last Mod: 08/07/17 by Don Hagler
%

%% todo: use mmil_reg instead of mmil_rbreg_vol2vol_jpdf

errcode = 0;
if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_ref',[],[],...
  'fname_out','./average_vol.mgh',[],...
  'ref_scantype','MPR',{'MPR','FLASHhi','MEDIChi'},...
  'input_scantype','MPR',{'MPR','XetaT2','FLASHhi','FLASHlo','MEDIChi','MEDIClo'},...
  'M_init',eye(4),[],...
  'regflag',true,[false true],...
  'image_range',[],[],...
  'cleanupflag',true,[false true],...
  'forceflag',false,[false true],...
...
  'binfile','reg',[],...
  'interpm',2,[0:5],...
  'smoothmask_flag',true,[false,true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input list
if isempty(fname_inlist), error('fname_inlist is empty'); end;
if ~iscell(fname_inlist), fname_inlist = {fname_inlist}; end;

% check input files
for i=1:length(fname_inlist)
  fname = fname_inlist{i};
  if ~exist(fname,'file')
    fprintf('%s: ERROR: input file %s not found\n',mfilename,fname);
    errcode = 1;
    return;
  end;
  [tmp,tmp_stem,tmp_ext]=fileparts(fname);
  if ~ismember(tmp_ext,{'.mgh','.mgz'})
    error('input files must be mgh or mgz format (%s)',fname);
  end;
end;

% check ref file
if isempty(parms.fname_ref)
  parms.fname_ref = fname_inlist{1};
else
  if ~exist(parms.fname_ref,'file')
    fprintf('%s: ERROR: ref file %s not found\n',mfilename,parms.fname_ref);
    errcode = 1;
    return;
  end;
  [tmp,tmp_stem,tmp_ext]=fileparts(parms.fname_ref);
  if ~ismember(tmp_ext,{'.mgh','.mgz'})
    error('ref file must be mgh or mgz format (%s)',parms.fname_ref);
  end;
end;

% check output file
[outdir,tmp_stem,tmp_ext]=fileparts(parms.fname_out);
if isempty(outdir), outdir = pwd; end;
if ~ismember(tmp_ext,{'.mgh','.mgz'})
  fname = parms.fname_out;
  parms.fname_out = [outdir '/' tmp_stem '.mgh'];
  fprintf('%s: WARNING: changing parms.fname_out %s to %s...',...
    mfilename,fname,parms.fname_out);
else
  parms.fname_out = [outdir '/' tmp_stem tmp_ext];
end;
if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% check mask file
[tmp,tmp_stem,tmp_ext]=fileparts(fname_mask);
if ~ismember(tmp_ext,{'.mgh','.mgz'})
  error('input files must be mgh or mgz format (%s)',fname_mask);
end;
if ~exist(fname_mask,'file')
  fprintf('%s: ERROR: mask file %s not found\n',mfilename,fname_mask);
  errcode = 1;
  return;
end;
[M,volsz_mask,mr_parms] = mmil_load_mgh_info(fname_mask,parms.forceflag,outdir);
[M,volsz,mr_parms] = mmil_load_mgh_info(parms.fname_ref,parms.forceflag,outdir);
if any(volsz~=volsz_mask)
  fprintf('%s: ERROR: mask file %s image dimensions do not match reference file %s\n',...
    mfilename,fname_mask,parms.fname_ref);
  errcode = 1;
  return;
end;

% create temporary directory for reg output
tmpdir = [outdir '/tmp_reg'];
mmil_mkdir(tmpdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(parms.ref_scantype,parms.input_scantype)
  fname_in = fname_inlist{1};
  [tmp,fstem,fext] = fileparts(fname_in);
  fname_res = [outdir '/' fstem '_res' fext];
  [tmp,fstem_ref] = fileparts(parms.fname_ref);
  if parms.regflag
    % register first input file to reference using jpdf
    fname_reg = [outdir '/' fstem '_reg_' fstem_ref '.mat'];
    if ~exist(fname_res,'file') || ~exist(fname_reg,'file') || parms.forceflag
      % load vol and ref
      fprintf('%s: loading %s...\n',mfilename,fname_in);
      [vol,mr_parms] = ctx_load_mgh(fname_in);
      fprintf('%s: loading %s...\n',mfilename,parms.fname_ref);
      vol_ref = ctx_load_mgh(parms.fname_ref);
      % use different registration tool for XetaT2      
      if strcmp(parms.input_scantype,'XetaT2')
        fname_T1 = parms.fname_ref;
        fname_T2 = fname_in;
        M_res_to_ref = eye(4);
        fprintf('%s: registering %s to %s with jpdf...\n',...
          mfilename,parms.input_scantype,parms.ref_scantype);
        M_ref_to_orig = mmil_jpdfreg_T1T2(fname_T1,fname_T2,...
          'outdir',tmpdir,'fname_T1_mask',fname_mask,...
          'smoothmask_flag',parms.smoothmask_flag,...
          'cleanup_flag',parms.cleanupflag,'forceflag',parms.forceflag);
      else
        fprintf('%s: loading %s...\n',mfilename,fname_mask);
        vol_mask = ctx_load_mgh(fname_mask);
        M_res_to_ref = parms.M_init;
        fprintf('%s: registering %s to %s with jpdf...\n',...
          mfilename,parms.input_scantype,parms.ref_scantype);
        vol.Mvxl2lph = inv(M_res_to_ref)*vol.Mvxl2lph;
        M_ref_to_orig = mmil_rbreg_vol2vol_jpdf(vol_ref,vol,...
          'type1',parms.ref_scantype,...
          'type2',parms.input_scantype,...
          'volmask',vol_mask);
      end;
      fprintf('%s: resampling...\n',mfilename);
      vol_res = vol_resample_pad(vol,vol_ref,M_ref_to_orig,parms.interpm);
      M_res_to_orig = M_ref_to_orig*M_res_to_ref;
      ctx_save_mgh(vol_res,fname_res,mr_parms);
      save(fname_reg,'M_res_to_orig','M_ref_to_orig','M_res_to_ref');
      %% todo: create scripts for checking registration?
    end;
  else
    % resample first input file to reference using M_init
    if ~exist(fname_res,'file') || parms.forceflag
      fprintf('%s: loading %s...\n',mfilename,parms.fname_ref);
      vol_ref = ctx_load_mgh(parms.fname_ref);
      fprintf('%s: loading %s...\n',mfilename,fname_in);
      [vol,mr_parms] = ctx_load_mgh(fname_in);
      M_res_to_ref = parms.M_init;
      fprintf('%s: resampling...\n',mfilename);
      vol.Mvxl2lph = inv(M_res_to_ref)*vol.Mvxl2lph;
      M_ref_to_orig = eye(4);
      vol_res = vol_resample_pad(vol,vol_ref,M_ref_to_orig,parms.interpm);
      ctx_save_mgh(vol_res,fname_res,mr_parms);
    end;
  end;
  parms.fname_ref = fname_res;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% register volumes and average
[tmpM,tmpvolsz,mr_parms] = mmil_load_mgh_info(fname_inlist{1},parms.forceflag,outdir);
nvols = 0;
vol_sum = [];

fprintf('%s: loading %s...\n',mfilename,parms.fname_ref);
[tmp,ref_stem,tmp_ext] = fileparts(parms.fname_ref);
vol_sum = ctx_load_mgh(parms.fname_ref);
nvols = 1;
for i=2:length(fname_inlist)
  fname = fname_inlist{i};
  [tmp,tmp_stem,tmp_ext] = fileparts(fname);
  fname_reg = [outdir '/' tmp_stem '_reg_' ref_stem '.mat'];
  if ~exist(fname_reg,'file') || parms.forceflag
    tmpdir2 = sprintf('%s/scan%d',tmpdir,i);
    mmil_mkdir(tmpdir2);
    fprintf('%s: registering %s to %s...\n',mfilename,fname,parms.fname_ref);
    M_volA_to_volB = mmil_reg(parms.fname_ref,fname,...
      'fname_maskA',fname_mask,...
      'binfile',parms.binfile,...
      'outdir',tmpdir2);
    % save output file
    save(fname_reg,'M_volA_to_volB');
  else
    load(fname_reg);
  end;
  fprintf('%s: loading %s...\n',mfilename,fname);
  vol = ctx_load_mgh(fname);
  fprintf('%s: resampling %s...\n',mfilename,fname);
  vol_res = vol_resample_pad(vol,vol_sum,M_volA_to_volB,parms.interpm);
  vol_sum.imgs = vol_sum.imgs + vol_res.imgs;
  nvols = nvols + 1;
end;
vol_sum.imgs = vol_sum.imgs / nvols;
% clip values
if ~isempty(parms.image_range)
  vol_sum.imgs(vol_sum.imgs<parms.image_range(1)) = parms.image_range(1);
  vol_sum.imgs(vol_sum.imgs>parms.image_range(2)) = parms.image_range(2);
end;
fprintf('%s: saving %s...\n',mfilename,parms.fname_out);
ctx_save_mgh(vol_sum,parms.fname_out,mr_parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove temporary files
if parms.cleanupflag
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result] = unix(cmd);
  if status, error('cmd %s failed:\n%s',cmd,status); end;
end;

