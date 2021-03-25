function M_atl_to_subj = mmil_rbreg_atlas(fname_subj,fname_atl,varargin)
%function M_atl_to_subj = mmil_rbreg_atlas(fname_subj,fname_atl,[options])
%
% Purpose: Rigidly register a subject brain image to an atlas brain image
%   Does a two-step registration, with mutual information for second iteration
%
% Required Input:
%   fname_subj: file name of subject brain volume (mgh/mgz format)
%   fname_atl: file name of atlas brain volume (mgh/mgz format)
%
% Optional Parameters:
%   'fname_mask': name of atlas brain mask for first iteration
%     if empty, will generate a brain mask using mmil_quick_brainmask
%     {default = []}
%   'fname_mask_broad': name of dilated atlas brain mask for second iteration
%     if empty, will generate from fname_mask using mmil_dilate_mask
%     {default = []}
%   'tmpdir': name of directory for temporary files
%     full path or relative to directory containing fname_subj
%     {default = 'rbreg_atlas'}
%   'cleanup_flag': [0|1] remove temporary files
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output (in tmpdir)
%     {default = 0}
%
% Created:  05/28/12 by Don Hagler (based on code by Vijay Venkatraman)
% Last Mod: 05/28/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
M_atl_to_subj = [];
parms = mmil_args2parms(varargin, { ...
  'fname_mask',[],[],...
  'fname_mask_broad',[],[],...
  'tmpdir','rbreg_atlas',[],...
  'cleanup_flag',true,[false,true],...
  'forceflag',false,[false,true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_subj,'file'), error('file %s not found',fname_subj); end;
if ~exist(fname_atl,'file'), error('file %s not found',fname_atl); end;

if ~isempty(parms.fname_mask) && ~exist(parms.fname_mask,'file')
  error('file %s not found',parms.fname_mask);
end;
if ~isempty(parms.fname_mask_broad) && ~exist(parms.fname_mask_broad,'file')
  error('file %s not found',parms.fname_mask_broad);
end;

[tpath,fstem_subj,text] = fileparts(fname_subj);

if mmil_isrelative(parms.tmpdir)
  parms.tmpdir = [tpath '/' parms.tmpdir];
end;

mmil_mkdir(parms.tmpdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create brain mask
if isempty(parms.fname_mask)
  parms.fname_mask = sprintf('%s/mask_atl.mgz',parms.tmpdir);
  mmil_quick_brainmask([],...
    'fname_in',fname_atl,...
    'fname_out',parms.fname_mask,...
    'forceflag',parms.forceflag);
end;

% dilate brain mask
if isempty(parms.fname_mask_broad)
  parms.fname_mask_broad = sprintf('%s/mask_atl_broad.mgz',parms.tmpdir);
  mmil_dilate_mask([],...
    'fname_in',parms.fname_mask,...
    'fname_out',parms.fname_mask_broad,...
    'forceflag',parms.forceflag);
end;

% register atlas to subject iteration 1
tmp_outdir_itr1 = sprintf('%s/rr',parms.tmpdir);
fname_itr1 = sprintf('%s/%s_reg2atl_itr1.mgz',tmp_outdir_itr1,fstem_subj);
M_atl_to_subj_itr1 = mmil_reg(fname_atl,fname_subj,...
  'fname_maskA',parms.fname_mask, ...
  'rigid_flag',1,'outdir',tmp_outdir_itr1,...
  'cleanup_flag',parms.cleanup_flag,'forceflag',parms.forceflag);

% resample subject to atlas
vol_subj = ctx_load_mgh(fname_subj);
[vol_atl,mr_parms] = ctx_load_mgh(fname_atl);  
vol_itr1 = vol_resample_pad(vol_subj,vol_atl,M_atl_to_subj_itr1,2,1);  
clear vol_subj vol_atl;
ctx_save_mgh(vol_itr1,fname_itr1,mr_parms);
clear vol_itr1;

% register atlas to subject iteration 2
tmp_outdir_itr2 = sprintf('%s/rrjpdf',parms.tmpdir);	
fname_itr2 = sprintf('%s/%s_reg2atl_itr2.mgz',tmp_outdir_itr2,fstem_subj);
M_atl_to_subj_itr2 = mmil_reg(fname_atl,fname_itr1,...
  'fname_maskA',parms.fname_mask_broad,...
  'rigid_flag',1,'options','-jpr -jpbrrMultScale',...
  'outdir',tmp_outdir_itr2,...
  'cleanup_flag',parms.cleanup_flag,'forceflag',parms.forceflag);

% remove temporary directory
if parms.cleanup_flag
  cmd = sprintf('rm -r %s',parms.tmpdir);
  [status,result] = unix(cmd);
  if status
    error('failed to delete the temporary directory %s:\n%s\n',parms.tmpdir,result);
  end;
end;

% combine transformation matrices
M_atl_to_subj = M_atl_to_subj_itr1 * M_atl_to_subj_itr2;

