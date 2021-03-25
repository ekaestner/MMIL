function fname_out = mmil_wm_corr(fname_in,varargin)
%function fname_out = mmil_wm_corr(fname_in,[options])
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
%  'n3flag': [0|1] apply N3 bias correction as pre-normalization step
%    {default = 0}
%  'fname_n3': output file for N3 bias corrected image
%    {default = 'n3.mgz'}
%  'vol_target': target volume for white matter intensity
%    may be a scalar value
%    if empty, will use median white matter intensity value
%    {default = []}
%  'ratio_range': range of values to constrain ratio
%    if empty, no clipping
%    {default = []}
%  'image_range': range of values to clip output image
%    if empty, no clipping
%    {default = []}
%  'lambda0': weighting factor for difference volume
%    {default = 1}
%  'lambda1': weighting factor for smooth volume
%    {default = 0}
%  'lambda2': weighting factor for Laplacian of smooth volume
%    higher values mean stronger smoothness constraint
%    {default = 5}
%  'niters_spsm': number of iterations of sparse smoothing
%    {default = 3}
%  'tmpdir': temporary directory containing intermediate output
%    created if fname_wm or fname_bm are empty
%    {default = 'tmp_wm_corr'}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 0}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  11/27/15 by Don Hagler
% Last Mod: 08/04/17 by Don Hagler
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

% white matter segmentation
parms = create_wmseg(parms);

% perform bias correction
parms = bias_corr(parms);

% optionally remove temporary directory
cleanup_tmpdir(parms);

fname_out = parms.fname_out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_out','corr.mgz',[],...
    'fname_wm',[],[],...
    'fname_bm',[],[],...
    'n3flag',true,[false true],...
    'fname_n3',[],[],...
    'vol_target',[],[],...
    'ratio_range',[],[],...
    'image_range',[],[],...
    'lambda0',1,[0,1e10],...
    'lambda1',0,[0,1e10],...
    'lambda2',5,[0,1e10],...
    'niters_spsm',3,[1,10],...
    'tmpdir','tmp_wm_corr',[],...
    'cleanup_flag',false,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % for N3 bias correction
    'estop',0.001,[],...
    'maxiter',100,[],...
  ... % for dilation of brain mask
    'smooth1',30,[],...
    'thresh1',0.05,[],...
    'smooth2',30,[],...
    'thresh2',0.05,[],...
    'smooth3',0,[],...
    'thresh3',0,[],...  
  ...
    'ext','.mgz',[],...
  ...
    'wmseg_tags',{'fname_wm','fname_bm','erode_wm_flag','dilate_bm_flag',...
                  'verbose','forceflag','estop','maxiter','nvoxels',...
                  'sigma','thresh_init','thresh','smooth1','thresh1',...
                  'smooth2','thresh2','smooth3','thresh3'},[],...
    'bias_tags',{'vol_target','ratio_range','image_range',...
                 'lambda0','lambda1','lambda2','niters_spsm'},[],...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;

  % create temporary output directory
  if isempty(parms.fname_wm) || isempty(parms.fname_bm) ||...
     (isempty(parms.fname_n3) && parms.n3flag)
    if mmil_isrelative(parms.tmpdir)
      parms.tmpdir = [pwd '/' parms.tmpdir];
    end;
    mmil_mkdir(parms.tmpdir);
  else
    parms.tmpdir = [];
  end;

  % set names of intermediate files
  if isempty(parms.fname_wm)
    parms.fname_wm = sprintf('%s/wmseg%s',parms.tmpdir,parms.ext);
  end;
  if isempty(parms.fname_bm)
    parms.fname_bm = sprintf('%s/brain%s',parms.tmpdir,parms.ext);
  end;
  if isempty(parms.fname_n3) && parms.n3flag
    parms.fname_n3 = sprintf('%s/n3%s',parms.tmpdir,parms.ext);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_wmseg(parms)
  % create white matter segmentation
  args = mmil_parms2args(parms,parms.wmseg_tags);
  mmil_wmseg(parms.fname_in,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = bias_corr(parms)
  if parms.n3flag
    if parms.verbose
      fprintf('%s: performing N3 bias correction...\n',mfilename);
    end;
    if ~exist(parms.fname_n3,'file') || parms.forceflag
      vol_brain = ctx_load_mgh(parms.fname_in);
      vol_brain = vol_correct_bias_field_n3(vol_brain,parms.estop,parms.maxiter);
      ctx_save_mgh(vol_brain,parms.fname_n3);
    end;
    parms.fname_in = parms.fname_n3;
  end;

  % perform bias correction using white matter segmentation
  if parms.verbose
    fprintf('%s: performing wm bias correction...\n',mfilename);
  end;
  [vol_in,M] = fs_load_mgh(parms.fname_in);
  vol_wm = fs_load_mgh(parms.fname_wm);
  vol_bm = fs_load_mgh(parms.fname_bm);
  vol_bm = 1.0 * (vol_bm > 0);
  args = mmil_parms2args(parms,parms.bias_tags);
  [vol_out,vol_ratio,vol_pred] = ...
    mmil_corr_wm_bias_spsm(vol_in,vol_wm,vol_bm,args{:});
  fs_save_mgh(vol_out,parms.fname_out,M);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanup_flag && ~isempty(parms.tmpdir) && exist(parms.tmpdir,'dir')
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,msg] = unix(cmd);
    if status
      error('tmpdir cleanup failed:\n%s\n%s\n',cmd,msg);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

