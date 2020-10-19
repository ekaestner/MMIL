function mmil_fiber_analysis(fname,fiber_dir,varargin)
%function mmil_fiber_analysis(fname,fiber_dir,[options])
%
% Required Parameters:
%   fname: full path name of file containing image data
%   fiber_dir: full path of directory containing fiber ROI files
%
% Optional Parameters:
%   'outdir': output directory
%     if empty, will attempt to use path of fname
%     {default = []}
%   'outstem': output file stem
%     if empty, will use file stem of fname
%     {default = []}
%   'csv_flag': [0|1] summarize results in csv files
%     {default = 1}
%   'fibers': fiber numbers to inclue in analysis
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so, what type
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     ignored if fiber_infix not empty
%     {default = 2}
%   'resT1flag': [0|1] use fibers resampled to T1 resolution
%     ignored if fiber_infix not empty
%     {default = 0}
%   'disp_flag': [0|1] generate fiber ROIs excluding CSF and gray matter
%     as defined by FreeSurfer aseg
%     {default = 0}
%   'disp_suffix':  suffix attached to output fiber file names
%       after calculating dispersion weighting
%     {default = 'dwtd'}
%   'dispfact': multiplicative factor applied to dispersion values
%     {default = 4}
%   'xcg_flag': [0|1] use fibers with CSF and gray matter excluded
%     ignored if fiber_infix not empty
%     {default = 0}
%   'xcg_suffix':  suffix attached to output fiber file names
%       after excluding CSF and gray matter
%     {default = 'xcg'}
%   'masksf_flag': [0|1] use fibers with multiple fiber voxels excluded
%     ignored if fiber_infix not empty
%     {default = 0}
%   'masksf_suffix': suffix attached to output fiber file names
%       after excluding voxels with multiple fibers
%     {default = 'masksf'}
%   'fiber_atlasname': name of fiber atlas used
%     if empty, default atlas used, no extra string included
%     {default=[])
%   'fiber_infix': string added to fiber ROI file name
%     e.g. fiber_101_<infix>.mat
%     if empty, set based on atlas_flag, resT1flag, xcg_flag, and masksf_flag
%     {default = []}
%   'fiber_ext': file extension for fiber ROI files
%     must be '.mgh' or '.mgz' or '.mat'
%     {default = '.mat'}
%   'M_reg': 4x4 transformation matrix between fname and fibers
%     should be like M_fibers_to_fname
%     if supplied, fname will be resampled to space of fibers
%     {default = []}
%   'res_outfix': output string added to resampled fname
%     ignored if M_reg is empty
%     {default = 'res'}
%   'fname_FA': full path of FA volume (required if thresh_FA not 0)
%     {default = []}
%   'thresh_FA': mask fiber ROIs by applying threshold to FA image
%     {default = 0}
%   'thresh_prob': mask atlas-derived fiber ROIs with probability threshold
%     {default = 0}
%   'weighted_avg_flag': [0|1] whether to calculate weighted averages for
%     each ROI using fiber counts (manual) or fiber probabilities (atlas)
%     to weight contribution from each voxel
%     {default = 1}
%   'scalefact': scaling factor applied to average values
%     {default = 1}
%   'minval': minimum value; if NaN, include all voxels
%     {default = 1e-6}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] run calculations even if output files exist
%     {default = 0}
%
% Created:  10/12/12 by Don Hagler
% Last Mod: 02/05/13 by Don Hagler
%

% based on DTI_MMIL_Analyze_Fibers, created 11/10/09 by Don Hagler
%      and DTI_MMIL_Analyze_Fibers_Exam, created 02/12/07 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% check parameters, check input files and directories
parms = check_input(fname,fiber_dir,varargin);

% get list of fiber files, apply thresh_FA and thresh_prob
parms = prep_fibers(parms);

% resample fname to fibers
if ~isempty(parms.M_reg)
  parms = resamp_input(parms);
end;

% calculate ROI averages
extract_values(parms);

% summarize results in csv file
if parms.csv_flag, write_csv(parms); end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,fiber_dir,options)
  % parse input arguments
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
    'fiber_dir',fiber_dir,[],...
  ...
    'outdir',[],[],...
    'outstem',[],[],...
    'csv_flag',true,[false true],...
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'atlas_flag',2,[0:4],...
    'resT1flag',false,[false true],...
    'disp_flag',false,[false true],...
    'disp_suffix','dwtd',[],...
    'dispfact',4,[1e-6,1e6],...
    'xcg_flag',false,[false true],...
    'xcg_suffix','xcg',[],...
    'masksf_flag',false,[false true],...
    'masksf_suffix','masksf',[],...
    'fiber_atlasname',[],[],...
    'fiber_infix',[],[],...
    'fiber_ext','.mat',{'.mgz','.mgh','.mat'},...
    'M_reg',[],[],...
    'res_outfix','res',[],...
    'fname_FA',[],[],...
    'thresh_FA',0,[0,1],...
    'thresh_prob',0,[0,1],...
    'weighted_avg_flag',true,[false true],...
    'scalefact',1,[-Inf,Inf],...
    'minval',1e-6,[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden parameters
    'outfix','fiber',[],...
    'frames',1,[],...
  ...
    'infix_tags',{'atlas_flag','count_flag','resT1flag',...
                  'xcg_flag','xcg_suffix',...
                  'masksf_flag','masksf_suffix',...
                  'disp_flag','disp_suffix','dispfact',...
                  'thresh_prob','thresh_FA'},[],...
    'roi_tags',{'frames','minval','weighted_avg_flag','scalefact','verbose'},[],...
    'suffix_tags',{'xcg_flag','xcg_suffix',...
                   'masksf_flag','masksf_suffix',...
                   'disp_flag','disp_suffix','dispfact',...
                   'thresh_prob','thresh_FA'},[],...
  });

  % check input file
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;

  % check fiber_dir
  if ~exist(fiber_dir,'dir')
    error('fiber_dir %s not found',fiber_dir);
  end;

  % check fname_FA if thresh_FA > 0
  if parms.thresh_FA>0
    if isempty(parms.fname_FA)
      error('must specify fname_FA if thresh_FA>0');
    end;
    if ~exist(parms.fname_FA,'file')
      error('file %s not found',parms.fname_FA);
    end;
  end;

  % set fiber_infix
  if isempty(parms.fiber_infix) || isempty(parms.fiber_ext)
    args = mmil_parms2args(parms,parms.infix_tags);
    [fiber_infix,fiber_ext] = dti_set_fiber_infix(args{:});
    if isempty(parms.fiber_infix), parms.fiber_infix = fiber_infix; end;
    if isempty(parms.fiber_ext), parms.fiber_ext = fiber_ext; end;
  end;

  % check fiber files
  flist = dir(sprintf('%s/fiber_*_%s%s',...
    parms.fiber_dir,parms.fiber_infix,parms.fiber_ext));
  if isempty(flist)
    error('no %s fiber files with fiber infix %s found in %s\n',...
      parms.fiber_ext,parms.fiber_infix,parms.fiber_dir);
  end;

  % set output directory and file stem
  if isempty(parms.outstem)
    [tmp,parms.outstem] = fileparts(parms.fname);
  end;
  % set output directory and create if necessary
  if ~mmil_isrelative(parms.outstem)
    parms.outdir = fileparts(parms.outstem);
  else
    % set output directory from outstem
    if isempty(parms.outdir)
      parms.outdir = fileparts(parms.outstem);
      if isempty(parms.outdir)
        parms.outdir = pwd;
      end;
    else
      parms.outstem = [parms.outdir '/' parms.outstem];
    end;
  end;

  % set outfix
  switch parms.atlas_flag
    case 1 % loc only, count atlas
      parms.outfix = [parms.outfix '_loc_countatlas'];
    case 2 % loc+dir, count atlas
      parms.outfix = [parms.outfix '_countatlas'];
    case 3 % loc only, mask atlas
      parms.outfix = [parms.outfix '_loc_atlas'];
    case 4 % loc+dir, mask atlas
      parms.outfix = [parms.outfix '_atlas'];
  end;
  if parms.weighted_avg_flag
    parms.outfix = [parms.outfix '_wtd'];
  end;
  args = mmil_parms2args(parms,parms.suffix_tags);
  suffix = dti_set_fiber_suffix(args{:});
  if ~isempty(suffix)
    parms.outfix = [parms.outfix '_' suffix];
  end;
  if ~isempty(parms.fiber_atlasname)
    parms.outfix = [parms.outfix '_' parms.fiber_atlasname];
  end
  parms.outfix = [parms.outfix '_roi_data'];

  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_fibers(parms)
  % create list of fiber files and get dimensions from first
  parms.volsz = []; parms.M = [];
  parms.nfibers = length(parms.fibers);
  parms.roi_files = cell(parms.nfibers,1);
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_fiber = sprintf('%s/fiber_%d_%s%s',...
      parms.fiber_dir,fnum,parms.fiber_infix,parms.fiber_ext);
    if ~exist(fname_fiber,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_fiber);
      continue;
    end;
    if isempty(parms.volsz)
      [parms.volsz,parms.M] = read_volsz(fname_fiber,parms);
      if parms.thresh_FA>0
        [vol_FA,M_FA,volsz_FA] = load_vol(parms.fname_FA);
        if any(volsz_FA~=parms.volsz) || any(M_FA(:)~=parms.M)
          error('dimensions of volume in %s does not match %s',...
            parms.fname_FA,fname_fiber);
        end;
      end;
    end;
    parms.roi_files{f} = fname_fiber;
  end;
return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_input(parms)
  % resample fname to fiber space
  fname_res = [parms.outstem '_' parms.res_outfix  '.mgz'];
  if ~exist(fname_res,'file') || parms.forceflag
    [vol,M_in] = fs_load_mgh(parms.fname);
    [vol,M_in] = mmil_resample_vol(vol,M_in,...
      'nvox_ref',parms.volsz,'M_ref',parms.M,'M_reg',parms.M_reg);
    fs_save_mgh(vol,fname_res,M_in);
  end
  parms.fname = fname_res;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extract_values(parms)
  % extract averages for fiber tract ROIs
  fname_mat = [parms.outstem '_' parms.outfix '.mat'];
  if ~exist(fname_mat,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: extracting fiber data from %s...\n',mfilename,parms.fname);
    end;
    args = mmil_parms2args(parms,parms.roi_tags);
    roi_data = mmil_multi_roi(parms.fname,parms.roi_files,args{:});
    if isempty(roi_data), error('fiber analysis failed'); end;
    for f=1:parms.nfibers
      roi_data(f).roicode = parms.fibers(f);
      %% todo: set roiname from fiber legend
    end;
    save(fname_mat,'roi_data');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(parms)
  % summarize results in csv file
  fname_mat = [parms.outstem '_' parms.outfix '.mat'];
  fname_csv = [parms.outstem '_' parms.outfix '.csv'];
  if ~exist(fname_csv,'file') || parms.forceflag
    roi_data = [];
    load(fname_mat);
    fid = fopen(fname_csv,'wt');
    if fid==-1, error('failed to open %s for writing',fname_csv); end;
    fprintf(fid,'"fiber #","mean","stdv","nvals","nvals valid"\n');
    for f=1:parms.nfibers
      fprintf(fid,'%d,%0.6f,%0.6f,%d,%d\n',...
        roi_data(f).roicode,roi_data(f).avg,...
        roi_data(f).stdv,roi_data(f).nvals,roi_data(f).nvals_valid);
    end;
    %% todo: use roiname instead of roicode
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
