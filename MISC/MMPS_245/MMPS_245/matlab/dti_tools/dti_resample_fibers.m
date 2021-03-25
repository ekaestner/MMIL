function dti_resample_fibers(fiber_dir,varargin)
%function dti_resample_fibers(fiber_dir,[options])
%
% Required Input:
%   fiber_dir: full path of directory containing fiber ROI files
%
% Optional Parameters:
%   'outdir': output directory
%     {default = fiber_dir}
%   'fibers': fiber numbers to resample
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so, what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%       (DTIStudio_fiber_masks directory must exist
%        -- imported with Import_DTIStudio_FiberMasks)
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'M_T1_to_DTI': registration matrix between DTI data and high-res T1 data
%     {default = eye(4)}
%   'volsz_T1': dimensions of high-res T1 volume
%     {default = [256 256 256]}
%   'M_T1': vox2ras matrix for high-res T1 volume
%     {default = [-1,0,0,129;0,0,1,-129;0,-1,0,129;0,0,0,1]}
%   'save_mgz_flag': [0|1] save fibers in mgz format in addition to sparse
%     {default = 0}
%   'verbose': display status messages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  10/10/12 by Don Hagler
% Last Mod: 01/15/13 by Don Hagler
%

% based on DTI_MMIL_Resample_Fibers, created 08/12/07 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters and check for problems
parms = check_input(fiber_dir,varargin);

% load fibers, resample to T1 resolution
resample_fibers(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fiber_dir,options)
  parms_filter = {...
    'fiber_dir',fiber_dir,[],...
  ...
    'outdir',fiber_dir,[],...
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'M_T1_to_DTI',eye(4),[],...
    'volsz_T1',[256,256,256],[],...
    'M_T1',[-1,0,0,129;0,0,1,-129;0,-1,0,129;0,0,0,1],[],...
    'atlas_flag',2,[0:4],...
    'save_mgz_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden parameters
    'fiber_outext',[],[],...
    'count_flag',true,[false true],...
    'interpm',1,[0:5],...
    'bclamp',true,[false,true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % check if input exists
  if ~exist(fiber_dir,'dir')
    error('fiber_dir %s not found',fiber_dir);
  end;
  parms.nfibers = length(parms.fibers);
  % set fiber_infix and fiber_ext
  [parms.fiber_infix,parms.fiber_ext] = dti_set_fiber_infix(...
    'count_flag',parms.count_flag,'atlas_flag',parms.atlas_flag);
  % check fiber files
  flist = dir(sprintf('%s/fiber_*_%s%s',...
    parms.fiber_dir,parms.fiber_infix,parms.fiber_ext));
  if isempty(flist)
    error('no %s fiber files with fiber infix %s found in %s\n',...
      parms.fiber_ext,parms.fiber_infix,parms.fiber_dir);
  end;
  % set fiber_outfix
  parms.fiber_outfix = dti_set_fiber_infix(...
    'resT1flag',1,...
    'count_flag',parms.count_flag,'atlas_flag',parms.atlas_flag);
  if isempty(parms.fiber_outext)
    parms.fiber_outext = parms.fiber_ext;
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resample_fibers(parms)
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_in = sprintf('%s/fiber_%02d_%s%s',...
      parms.fiber_dir,fnum,parms.fiber_infix,parms.fiber_ext);
    fname_out = sprintf('%s/fiber_%02d_%s%s',...
      parms.outdir,fnum,parms.fiber_outfix,parms.fiber_outext);
    if exist(fname_in,'file')
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: loading fiber file %s...\n',mfilename,fname_in);
        end;
        [vol,M] = load_vol(fname_in);
        if parms.verbose
          fprintf('%s: resampling volume to T1 resolution...\n',mfilename);
        end;
        [vol,M] = resample_vol(vol,M,parms);
        if parms.verbose
          fprintf('%s: saving output fiber file %s...\n',mfilename,fname_out);
        end;
        save_vol(vol,fname_out,M);
      end;
      if parms.save_mgz_flag && strcmp(parms.fiber_ext,'.mat')
        mmil_sparse2mgh(fname_out,[],parms.forceflag)
      end;
    else
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,M] = resample_vol(vol,M,parms)
  % resample vol using M_T1_to_DTI
  [vol,M] = mmil_resample_vol(vol,M,...
    'nvox_ref',parms.volsz_T1,'M_ref',parms.M_T1,...
    'interpm',parms.interpm,'bclamp',parms.bclamp,...
    'M_reg',parms.M_T1_to_DTI);
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

