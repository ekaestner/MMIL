function epi_resample_data(fname,varargin)
%function epi_resample_data(fname,[options])
%
% Purpose: Resample diffusion data to have the desired resolution.
%          Optionally, can register to T1 volume, rotate, and translate
%
% Required Parameters:
%   fname: full path name of mgh/mgz file containing 4d EPI volume
%
% Optional Parameters:
%   'fname_out': output file fname for corrected data
%     If empty, will append '_res.mgz' to file stem of fname
%     {default=[]}
%   'save_reg_flag': [0|1] whether to save new registration info in fname_reg
%     {default = 1}
%   'fname_reg': output file name for adjusted registration info
%     If empty, will append '_reg.mat' to file stem of fname
%     {default = []}
%   'qmat': matrix of diffusion directions
%     If supplied, will adjust directions for rotation related to
%     between scan registration, registration to T1, and additional rotation
%     and save adjusted qmat in fname_qmat
%     {default = []}
%   'fname_qmat': output file name for adjusted qmat and registration info
%     If empty, will append '_res_qmat.mat' to file stem of fname
%     {default = []}
%   'native_flag': resample but keep original nvoxels and resolution
%     {default = 0}
%   'nvoxels': vector of [nx,ny,nz] (in numbers of voxels)
%     ignored if native_flag = 1
%     {default = [120 120 70]}
%   'resolution': vector of [resx,resy,resz] (in mm)
%     ignored if native_flag = 1
%     {default = [2 2 2]}
%   'std_orient' : specify the resampled output orientation (e.g. 'LPI', 'RAS')
%     ignored if native_flag = 1 or empty
%      {default = []}   
%   'deoblique_flag': [0|1] whether to resample oblique slices to on-axis
%     ignored if native_flag = 1
%     {default = 1}
%   'M_reg': 4x4 registration matrix from reference volume to fname
%     (i.e. for between scan motion)
%     {default = []}
%   'M_ref': 4x4 vox2ras matrix for inter-scan reference
%     if empty, use M from fname
%     {default=[]}
%   'volsz_ref': vector of number of voxels for inter-scan reference
%     if empty, assume same as vol from fname
%     {default=[]}
%   'regT1flag': [0|1] whether to align diffusion data with T1 scan
%     requires fname_T1 and M_T1_to_EPI
%     {default = 0}
%   'fname_T1': full path name of T1-weighted image
%     required if regT1flag = 1
%     {default = []}
%   'M_T1_to_EPI': 4x4 registration matrix between T1 volume and
%     reference EPI volume
%     {default = []}
%   'smooth': isotropic smoothing kernel sigma (in voxels)
%     {default = 0}
%   'rot': three value vector of x, y, z rotation (in degrees)
%     order of rotations: 'x','y','z'
%     {default: [0,0,0]}
%   'trans': three value vector of x, y, z translation (in mm)
%     {default: [0,0,0]}
%   'interpm': [0|1|2|3|4|5] interpolation method
%      0:nearest  1:linear  2:cubic  3:Key's spline  4:cubic spline  5: hamming sinc
%     {default = 2}
%   'bclamp' : [0|1] set negative values to zero
%     {default = 1}
%   'EPI_type': type of EPI scan (e.g. 'EPI','BOLD','DTI')
%      used for output T1_resEPI.mgz file name
%     {default = 'EPI'}
%   'ext': extension used for output T1_resEPI file (i.e. '.mgh' or '.mgz')
%     {default = '.mgz'}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Created:  03/08/10 by Don Hagler
% Last Mod: 04/08/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
errcode = 0;
outfix = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options, set defaults

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_out',[],[],...
  'save_reg_flag',true,[false true],...
  'fname_reg',[],[],...  
  'qmat',[],[],...
  'fname_qmat',[],[],...  
  'native_flag',false,[false true],...
  'nvoxels',[120,120,70],[],...
  'resolution',[2 2 2],[],...
  'deoblique_flag',true,[false true],...
  'std_orient',[],[],...  
  'M_reg',[],[],...
  'M_ref',[],[],...
  'volsz_ref',[],[],...
  'regT1flag',false,[false true],...
  'fname_T1',[],[],...
  'M_T1_to_EPI',[],[],...
  'smooth',0,[0,10],...
  'rot',[0 0 0],[],...
  'trans',[0 0 0],[],...
  'interpm',2,[1:5],...
  'bclamp',true,[false,true],...
  'EPI_type','EPI',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'resample_tags',{'nvoxels','resolution','orient','M_ref','nvox_ref',...
                   'M_reg','deoblique_flag','smooth','interpm','bclamp'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fname_out)
  [tpath,tstem,text] = fileparts(fname);
  parms.fname_out = [tpath '/' tstem '_res' text];
else
  [tpath,tstem,text] = fileparts(parms.fname_out);  
end;
parms.outdir = tpath;
if isempty(parms.fname_qmat)
  [tpath,tstem,text] = fileparts(fname);
  parms.fname_qmat = [tpath '/' tstem '_res_qmat.mat'];
end;
if isempty(parms.fname_reg)
  [tpath,tstem,text] = fileparts(fname);
  parms.fname_reg = [tpath '/' tstem '_reg.mat'];
end;

if exist(parms.fname_out,'file') &&...
    (exist(parms.fname_qmat,'file') || isempty(parms.qmat)) &&...
    (~parms.save_reg_flag || exist(parms.fname_reg,'file') || isempty(parms.M_T1_to_EPI)) &&...
   ~parms.forceflag
  return;
end;

if parms.native_flag, parms.deoblique_flag = false; end;

parms.orient = parms.std_orient;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load header from EPI
[M_EPI,volsz_EPI] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);

if parms.native_flag
  parms.nvoxels = volsz_EPI(1:3);
  parms.resolution = sqrt(sum(M_EPI(1:3,1:3).^2,1));
end;

% load header from T1
if ~isempty(parms.M_T1_to_EPI)
  if isempty(parms.fname_T1)
    error('M_T1_to_EPI not empty but missing fname_T1');
  end;
  if ~exist(parms.fname_T1)
    error('file %s not found',parms.fname_T1);
  end;
  [M_T1,volsz_T1] = mmil_load_mgh_info(parms.fname_T1,parms.forceflag);
else  
  M_T1 = [];
  volsz_T1 = [];
end;

if isempty(parms.volsz_ref)
  parms.volsz_ref = volsz_EPI;
end;
parms.volsz_ref = parms.volsz_ref(1:3);

% apply registration to T1
if parms.regT1flag
  parms.M_ref = M_T1;
  parms.nvox_ref = volsz_T1(1:3);
  M_ref_to_EPI = parms.M_T1_to_EPI;
  M_T1_to_EPI = eye(4); % new T1 to EPI registration (identity)
else
  if isempty(parms.M_ref)
    parms.M_ref = M_EPI;
  end;
  parms.nvox_ref = parms.volsz_ref;
  M_ref_to_EPI = eye(4);
  M_T1_to_EPI = parms.M_T1_to_EPI;
end;

% apply registration to between-scan reference
if ~isempty(parms.M_reg)
  M_ref_to_EPI = parms.M_reg*M_ref_to_EPI;
end;

% apply rotation and translation to M_ref_to_EPI
orig_M_ref_to_EPI = M_ref_to_EPI;
if any(parms.rot~=0)
  M_rot = mmil_construct_M('rot',parms.rot,'trans',parms.trans);
  M_ref_to_EPI = M_ref_to_EPI*inv(M_rot);
  if ~isempty(M_T1_to_EPI)
    M_T1_to_EPI = M_rot*M_T1_to_EPI; % new T1 to EPI registration
  end;
end;

parms.M_reg = M_ref_to_EPI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample volume
M_res = [];
if ~exist(parms.fname_out,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: resampling %s...\n',mfilename,fname);
  end;
  [vol,M,mr_parms] = fs_load_mgh(fname);
  args = mmil_parms2args(parms,parms.resample_tags);
  [vol_res,M_res] = mmil_resample(vol,M,args{:});
  volsz_res = size(vol_res);
  if parms.verbose
    fprintf('%s: saving output as %s...\n',mfilename,parms.fname_out);
  end;
  fs_save_mgh(vol_res,parms.fname_out,M_res,mr_parms);
  clear vol vol_res;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(parms.rot~=0) || any(parms.trans~=0)
  rot = parms.rot;
  trans = parms.trans;
else
  rot = 0;
  trans = 0;
end;

RegInfo = [];
if ~isempty(parms.M_T1_to_EPI) && (parms.save_reg_flag || ~isempty(parms.qmat))
  [M_res,volsz_res] = mmil_load_mgh_info(parms.fname_out,parms.forceflag);
  RegInfo.fname_T1 = parms.fname_T1;
  RegInfo.fname_T2 = parms.fname_out;
  RegInfo.M_T1_to_T2 = M_T1_to_EPI;
  RegInfo.M_T1 = M_T1;
  RegInfo.M_T2 = M_res;
  if ~isempty(volsz_T1)
    RegInfo.volsz_T1 = volsz_T1(1:3);
  else
    RegInfo.volsz_T1 = parms.nvox_ref;
  end;
  RegInfo.volsz_T2 = volsz_res(1:3);
end;

if parms.save_reg_flag & ~isempty(RegInfo)
  % resample T1 to new EPI space
  if parms.regT1flag
    fname_T1_lowres = sprintf('%s/T1_res%s_regT1%s',...
      parms.outdir,parms.EPI_type,parms.ext);
  else
    fname_T1_lowres = sprintf('%s/T1_res%s%s',...
      parms.outdir,parms.EPI_type,parms.ext);
  end;
  if ~exist(fname_T1_lowres,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: resampling T1 volume to new EPI space...\n',mfilename);
    end;
    [vol_T1,M_T1] = fs_load_mgh(parms.fname_T1);
    [vol_T1_lowres,M_T1_lowres] = mmil_resample_vol(vol_T1,M_T1,...
      'M_ref',M_res,'nvox_ref',parms.nvoxels,...
      'interpm',parms.interpm,'bclamp',1,...
      'M_reg',inv(M_T1_to_EPI));
    fs_save_mgh(vol_T1_lowres,fname_T1_lowres,M_T1_lowres);
  end;
  save(parms.fname_reg,...
    'RegInfo','M_ref_to_EPI','orig_M_ref_to_EPI','M_T1_to_EPI','rot','trans');
  mmil_write_regscripts(parms.fname_reg);
end;

% adjust qmat if needed
if ~isempty(parms.qmat)
  qmat = parms.qmat;
  orig_qmat = qmat;
  if parms.regT1flag || any(parms.rot~=0)
    % adjust qmat for rotation component of registration
    qmat = dti_rotate_qmat(qmat,inv(M_ref_to_EPI));
  end;
  save(parms.fname_qmat,'qmat','orig_qmat',...
    'RegInfo','M_ref_to_EPI','orig_M_ref_to_EPI','M_T1_to_EPI','rot','trans');
end;


