function errcode=mmil_average_volumes_old(fname_inlist,fname_mask,varargin)
%function errcode=mmil_average_volumes_old(fname_inlist,fname_mask,[options])
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
%       Allowed: 'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo'
%       If ref and input scantypes are different, will use jpdf registration
%         for first input file, then mriRegister for subsequent scans
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
%    'cleanupflag': [0|1] whether to remove temporary files created by mriRegister
%      { default: 1 }
%    'paramfile' - full or relative path of parameter file
%      If relative, relative to $MMPS_PARMS/MRIREG
%      {default = 'inputParamsRigid.txt'}
%    'forceflag': [0|1] whether to run calculations even if output file exists
%      { default: 0 }
%
%  Output:
%   errcode: 0 if successful, 1 if error
%
% Created:  07/09/08 by Don Hagler
% Last Mod: 06/24/16 by Don Hagler
%

errcode = 0;
if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_ref',[],[],...
  'fname_out','./average_vol.mgh',[],...
  'ref_scantype','MPR',{'MPR','FLASHhi','MEDIChi'},...
  'input_scantype','MPR',{'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo'},...
  'M_init',eye(4),[],...
  'regflag',true,[false true],...
  'image_range',[],[],...
  'cleanupflag',true,[false true],...
  'paramfile','inputParamsRigid.txt',[],...
  'forceflag',false,[false true],...
...
  'parmsdir',[],[],...
  'interpm',2,[0:5],...
});

if isempty(parms.parmsdir)
  parms.parmsdir = [getenv('MMPS_PARMS') '/MRIREG'];
end;
if isempty(regexp(parms.paramfile,'^/')) % relative path
  parms.paramfile = [parms.parmsdir '/' parms.paramfile];
end;

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
  [tmp_path,tmp_stem,tmp_ext]=fileparts(fname);
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
  [tmp_path,tmp_stem,tmp_ext]=fileparts(parms.fname_ref);
  if ~ismember(tmp_ext,{'.mgh','.mgz'})
    error('ref file must be mgh or mgz format (%s)',parms.fname_ref);
  end;
end;

% check mask file
[tmp_path,tmp_stem,tmp_ext]=fileparts(fname_mask);
if ~ismember(tmp_ext,{'.mgh','.mgz'})
  error('input files must be mgh or mgz format (%s)',fname_mask);
end;
if ~exist(fname_mask,'file')
  fprintf('%s: ERROR: mask file %s not found\n',mfilename,fname_mask);
  errcode = 1;
  return;
end;
[vol,M,mr_parms,volsz_mask]=fs_load_mgh(fname_mask,[],[],1);
[vol,M,mr_parms,volsz]=fs_load_mgh(parms.fname_ref,[],[],1);
if any(volsz~=volsz_mask)
  fprintf('%s: ERROR: mask file %s image dimensions do not match reference file %s\n',...
    mfilename,fname_mask,parms.fname_ref);
  errcode = 1;
  return;
end;

% check output file
[tmp_path,tmp_stem,tmp_ext]=fileparts(parms.fname_out);
if isempty(tmp_path), tmp_path = pwd; end;
if ~ismember(tmp_ext,{'.mgh','.mgz'})
  fname = parms.fname_out;
  parms.fname_out = [tmp_path '/' tmp_stem '.mgh'];
  fprintf('%s: WARNING: changing parms.fname_out %s to %s...',...
    mfilename,fname,parms.fname_out);
else
  parms.fname_out = [tmp_path '/' tmp_stem tmp_ext];
end;

if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(parms.ref_scantype,parms.input_scantype)
  fname_in = fname_inlist{1};
  [tmp_path,tmp_stem,tmp_ext]=fileparts(fname);
  if isempty(tmp_path), tmp_path = pwd; end;
  fname_res = [tmp_path '/' tmp_stem '_res' tmp_ext];
  if parms.regflag
    % register first input file to reference using jpdf
    fname_reg = [tmp_path '/' parms.input_scantype 'reg2' parms.ref_scantype '.mat'];
    if ~exist(fname_res,'file') || ~exist(fname_reg,'file') || parms.forceflag
      fprintf('%s: loading %s...\n',mfilename,parms.fname_ref);
      vol_ref = ctx_load_mgh(parms.fname_ref);
      fprintf('%s: loading %s...\n',mfilename,fname_mask);
      vol_mask = ctx_load_mgh(fname_mask);
      fprintf('%s: loading %s...\n',mfilename,fname_in);
      [vol,mr_parms] = ctx_load_mgh(fname_in);
      M_res_to_orig = [];
      M_res_to_ref = parms.M_init;
      fprintf('%s: registering %s to %s with jpdf...\n',...
        mfilename,parms.input_scantype,parms.ref_scantype);
      vol.Mvxl2lph = inv(M_res_to_ref)*vol.Mvxl2lph;
      M_ref_to_orig = mmil_rbreg_vol2vol_jpdf(vol_ref,vol,...
        'type1',parms.ref_scantype,...
        'type2',parms.input_scantype,...
        'volmask',vol_mask);
      vol_res = vol_resample_pad(vol,vol_ref,M_ref_to_orig,parms.interpm);
      M_res_to_orig = M_ref_to_orig*M_res_to_ref;
      ctx_save_mgh(vol_res,fname_res,mr_parms);
      save(fname_reg,'M_res_to_orig','M_ref_to_orig','M_res_to_ref');
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

% create temporary directory for mriRegister output
if isempty(tmp_path), tmp_path = pwd; end;
tmpdir = [tmp_path '/tmp_mriregister'];
if ~exist(tmpdir,'dir')
  [success,msg] = mkdir(tmpdir);
  if ~success, error('failed to create tmpdir %s:\n%s',tmpdir,msg); end;
end;
% register volumes and average
[tmpvol,tmpM,mr_parms] = fs_load_mgh(fname_inlist{1},[],[],1); % get mr_parms for first file
nvols = 0;
vol_sum = [];

fprintf('%s: loading %s...\n',mfilename,parms.fname_ref);
vol_sum = ctx_load_mgh(parms.fname_ref);
nvols = 1;
for i=2:length(fname_inlist)
  fname = fname_inlist{i};
  [tmp_path,tmp_stem,tmp_ext] = fileparts(fname);
  orient = fs_read_orient(fname);
  if ~strcmp(orient,'LIA')
    fname_LIA = sprintf('%s/%s_LIA%s',tmpdir,tmp_stem,tmp_ext);
    fname_tmp = sprintf('%s/%s_LIA_TP2_RigidReg%s',tmpdir,tmp_stem,tmp_ext);
  else
    fname_LIA = [];
    fname_tmp = sprintf('%s/%s_TP2_RigidReg%s',tmpdir,tmp_stem,tmp_ext);
  end;
  if ~exist(fname_tmp,'file') || parms.forceflag
    % reorient to LIA if needed
    if ~strcmp(orient,'LIA')
      mmil_reorient_to_LIA(fname,fname_LIA);
      fname = fname_LIA;
    end;
    fprintf('%s: registering %s to %s...\n',mfilename,fname,parms.fname_ref);
    tic
    cmd = sprintf('mriRegister -ip %s -s %s -sm %s -t %s -od %s',...
      parms.paramfile,parms.fname_ref,fname_mask,fname,tmpdir);
    disp(cmd);
    [status,result] = unix(cmd);
    toc
    disp(result);
    if status
      error('cmd %s failed:\n%s',cmd,result);
    end;
    if ~exist(fname_tmp)
      error('file %s not created as expected',fname_tmp);
    end;
  end;
  fprintf('%s: loading %s...\n',mfilename,fname_tmp);
  vol = ctx_load_mgh(fname_tmp);
  vol_sum.imgs = vol_sum.imgs + vol.imgs;
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
% cleanup

if parms.cleanupflag
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result] = unix(cmd);
  if status, error('cmd %s failed:\n%s',cmd,status); end;
end;
