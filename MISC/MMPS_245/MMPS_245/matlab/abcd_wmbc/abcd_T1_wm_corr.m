function fname_out = abcd_T1_wm_corr(fname_in,varargin)
%function fname_out = abcd_T1_wm_corr(fname_in,[options])
%
% Required Input:
%   fname_in: full path of input file (mgh/mgz format)
%
% Optional Parameters:
%  'fname_out': output file name for bias corrected image (mgh/mgz format)
%    {default = 'corr.mgz'}
%  'fname_wm': output file for white matter segmentation
%    {default = 'wmseg.mgz'}
%  'fname_bm': output file for brain mask
%    {default = 'brainmask.mgz'}
%  'image_range': range of values to clip output image
%    if empty, no clipping
%    {default = []}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  11/06/17 by Don Hagler
% Last Mod: 11/06/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output file
if exist(parms.fname_out) && ~parms.forceflag
  fname_out = parms.fname_out;
  return;
end;

% segment white matter and perform bias correction
bias_corr(parms);

fname_out = parms.fname_out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_out','corr.mgz',[],...
    'fname_wm','wmseg.mgz',[],...
    'fname_bm','brainmask.mgz',[],...
    'image_range',[],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'ext','.mgz',[],...
  });
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = bias_corr(parms)
  % perform bias correction using white matter segmentation
  if parms.verbose
    fprintf('%s: performing wm bias correction...\n',mfilename);
  end;
  vol_orig = ctx_load_mgh(parms.fname_in);
  [vol_corr,vol_wm,vol_bm,vol_gm] = ABCD_wmbc(vol_orig);
  if ~isempty(parms.image_range)
    vol_corr.imgs(vol_corr.imgs<parms.image_range(1)) = parms.image_range(1);
    vol_corr.imgs(vol_corr.imgs>parms.image_range(2)) = parms.image_range(2);
  end;
  ctx_save_mgh(vol_corr,parms.fname_out);
  ctx_save_mgh(vol_wm,parms.fname_wm);
  ctx_save_mgh(vol_bm,parms.fname_bm);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

