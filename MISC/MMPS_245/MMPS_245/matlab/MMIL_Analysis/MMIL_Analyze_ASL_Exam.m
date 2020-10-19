function errcode = MMIL_Analyze_ASL_Exam(InPath,OutPath,FSContainerPath,varargin)
%function errcode = MMIL_Analyze_ASL_Exam(InPath,OutPath,FSContainerPath,[options])
%
% Purpose:
%  Extract ROI averages of arterial spin labeling (ASL) measure CBF
%  ROIs include those derived from freesurfer's subcortical segmentation (aseg)
%    and cortical surface parcellation (aparc)
%
% Usage:
%  MMIL_Analyze_ASL_Exam(InPath,OutPath,FSContainerPath,'key1', value1,...)
%
% Required Parameters:
%   InPath: full path of directory containing processed ASL data (BRIK format)
%   OutPath: full path of output directory to be created
%   FSContainerPath: full path of directory containing freesurfer recon
%
% Optional Parameters:
%   'outdir': analysis output directory
%     full path or relative to OutPath
%     {default = 'analysis'}
%   'regdir': output directory for registration files
%     full path or relative to OutPath
%     {default = 'reg2T1'}
%   'aseg_flag': [0|1] extract aseg ROI results
%     {default = 1}
%   'cortsurf_flag': [0|1] paint to cortical surface
%     and extract aparc ROI results
%     {default = 1}
%   'measlist': list of ASL "measures" to extract for each ROI
%      e.g. 'CBF'
%      'CBF' is cerebral blood flow
%     {default = {'CBF'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     {default = [1]}
%   'minval': minimum value; if NaN, include all voxels
%     {default = 1e-6}
%   'csv_flag': [0|1] output results in csv files
%     {default = 1}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Optional Parameters specific to aseg ROIs:
%   'erode_flag': [0|1] whether to "erode" ROIs
%     by smoothing and thresholding (to reduce edge effects)
%     {default = 1}
%   'erode_nvoxels': number of voxels to erode (integer)
%     {default = 1}
%   'aseg_aparc_flag': [0|1|2] whether to use cortical parcellation ROIs
%     0: aseg only
%     1: aparc only   NOTE: for best results with aparc, use erode_flag=0
%     2: aparc+aseg
%     {default = 0}
%
% Optional Parameters specific to cortical surface:
%   'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     if empty, will use ?h.aparc.annot files in fspath/label
%     if not full path, assumed to be relative to fspath/label
%     {default = []}
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     will be resampled to to individual subject space before use
%     {default = []}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast maps
%     {default = 1}
%   'gwnorm_flag': [0|1] for gray/white contrast, normalize difference by mean
%     {default = 0}
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1[ whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere
%     {default = 0}
%
% Optional Parameters for resampling volumes to atlas space
%   'atlas_warp_flag': [0|1] nonlinearly resample input volumes to atlas space
%     {default = 0}
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%
% Created:  02/07/13 by Don Hagler
% Last Mod: 09/23/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
errcode = 0;

[parms,errcode] = check_input(InPath,OutPath,FSContainerPath,varargin);
if errcode, return; end;
if ~parms.aseg_flag && ~parms.cortsurf_flag &&...
   ~parms.atlas_warp_flag
  fprintf('%s: WARNING: nothing to do\n',mfilename);
  return;
end;

mmil_mkdir(parms.outdir);
[parms,errcode] = prep_data(parms);
if errcode, return; end;

if ~isempty(parms.fnames_fparc)
  parms = resample_fparc_annot(parms);
end;

errcode = analyze_data(parms);

if parms.atlas_warp_flag
  errcode = atlas_warp_data(parms);
end;

if errcode
  fprintf('\n%s: WARNING: one or more analysis steps failed\n\n',mfilename);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input(InPath,OutPath,FSContainerPath,options)
  errcode = 0;
  parms = mmil_args2parms(options,{...
    'inpath',InPath,[],...
    'outpath',OutPath,[],...
    'fspath',FSContainerPath,[],...
  ...
    'outdir','analysis',[],...
    'regdir','reg2T1',[],...
    'aseg_flag',true,[false true],...
    'cortsurf_flag',false,[false true],...
    'measlist',{'CBF'},[],...
    'scalefacts',[1],[],...
    'csv_flag',true,[false true],...
    'minval',1e-6,[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % aseg parameters
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
  ... % cortsurf parameters
    'fnames_aparc',[],[],...
    'fnames_fparc',[],[],...
    'projdist_list',[1],[-5,5],...
    'gwnorm_flag',false,[false true],...
    'smoothsteps',0,[0,Inf],...
    'sphere_flag',false,[false true],...
    'sphsmoothsteps',0,[0,Inf],...
  ... % resample atlas parameters
    'atlas_warp_flag',false,[false true],...
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas',[],...    
  ... % hidden parameters
    'input_files',{'intermediate/anat+reg+rs+orig.BRIK','CBF+orig.BRIK'},[],...
    'output_files',{'T1_lowres.mgz','CBF.mgz'},[],...
    'resT1_outfix','resT1',[],...
    'fnum',1,[1,Inf],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'fname_colorlut',[],[],...
  ... % sets of parameter names to pass to other functions
    'aseg_tags',{'outdir','outstem','fname_out','csv_flag','fname_aseg',...
                 'aseg_aparc_flag','fname_vals','dispvec',...
                 'disp_roicodes','disp_scalefact','dispfact',...
                 'erode_flag','erode_nvoxels',...
                 'scalefact','frames','minval','M_reg','res_outfix','fname_colorlut',...
                 'verbose','forceflag'},[],...
    'cortsurf_tags',{'outdir','outstem','fnames_aparc','csv_flag','M_reg',...
                     'resT1flag','res_outfix','projdist_list',...
                     'gwnorm_flag','smoothsteps','sphere_flag',...
                     'sphsmoothsteps',...
                     'mask_midbrain_flag',...
                     'scalefact','frames','minval',...
                     'fname_colorlut','verbose','forceflag','hemilist'},[],...
    'warp_tags',{'atlasdir','atlasname'...
                 'smoothflag','sampling','nK','tstep','astep',...
                 'scales','ns','sf','thresh','stdflag','maskflag',...
                 'stdbgval','fname_reg','fname_T1',...
                 'outdir','outstem','outstem_T1',...
                 'interpm','bclamp','padding','vxlmap_flag','forceflag'},[],...
    'mask_tags',{'fname_mask','brain_flag','force_flag'},[],...
    'reg_tags',{'outdir','fname_maskA','fname_maskB','mask_flag','affine_flag'},[],...
  });

  parms.nhemi = length(parms.hemilist);

  % check that number of measures matches number of scalefacts
  if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;
  parms.nmeas = length(parms.measlist);
  if strcmp(parms.measlist{1},'none')
    parms.measlist = [];
  elseif length(parms.scalefacts)==1
    parms.scalefacts = parms.scalefacts*ones(parms.nmeas,1);
  elseif length(parms.measlist)~=length(parms.scalefacts)
    error('number of measures in measlist does not match number of scalefacts');
  end;

  % check that input data container exists
  if ~exist(parms.inpath,'dir')
    fprintf('%s: ERROR: %s not found\n',mfilename,parms.inpath);
    errcode = 1; return;
  end;

  % check that Freesurfer recon exists
  if ~exist(parms.fspath,'dir')
    fprintf('%s: ERROR: %s not found\n',mfilename,parms.fspath);
    errcode = 1; return;
  end
  parms.fname_T1 = [parms.fspath '/mri/orig.mgz'];
  if ~exist(parms.fname_T1,'file')
    fprintf('%s: ERROR: %s not found\n',mfilename,parms.fname_T1);
    errcode = 1; return;
  end;

  % set subj name
  [parms.subjdir,parms.subj,tmp] = fileparts(parms.fspath);
  parms.subj = [parms.subj tmp];

  % check input files
  errcode = check_input_files(parms);
  if errcode, return; end;
  
  % check that aseg exists
  if parms.aseg_flag || parms.brainmask_flag ||...
     (parms.fiber_flag && parms.xcg_flag)
    parms = check_aseg(parms);
  end;

  % check cortsurf
  if parms.cortsurf_flag
    parms = check_cortsurf(parms);
  end;

  % set regdir
  if mmil_isrelative(parms.regdir)
    parms.regdir = [parms.outpath '/' parms.regdir];
  end;

  % set outdir
  if mmil_isrelative(parms.outdir)
    parms.outdir = [parms.outpath '/' parms.outdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = check_input_files(parms)
  errcode = 0;
  % check that all necessary files are found in parms.inpath
  for i=1:length(parms.input_files)
    fname = [parms.inpath '/' parms.input_files{i}];
    if ~exist(fname,'file')
      fprintf('%s: ERROR: missing %s\n',mfilename,fname);
      errcode = 1;
      return;
    end;      
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_aseg(parms)
  if parms.verbose
    fprintf('%s: checking aseg...\n',mfilename);
  end;
  % set name of aseg
  if parms.aseg_aparc_flag
    parms.aseg_name = 'aparc+aseg';
  else
    parms.aseg_name = 'aseg';
  end;
  parms.fname_aseg = sprintf('%s/mri/%s.mgz',parms.fspath,parms.aseg_name);
  if ~exist(parms.fname_aseg,'file')
    fprintf('%s: WARNING: aseg file %s not found\n',...
      mfilename,parms.fname_aseg);
    parms.aseg_flag = 0;
    parms.brainmask_flag = 0;
    if parms.xcg_flag, parms.fiber_flag = 0; end;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_cortsurf(parms)
  if parms.verbose
    fprintf('%s: checking cortical surface...\n',mfilename);
  end;
  parms.aparc_hemis = [];
  if isempty(parms.fnames_aparc)
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      parms.aparc_hemis{h} = hemi;
      parms.aparc_names{h} = 'aparc';
      parms.fnames_aparc{h} = sprintf('%s/label/%s.aparc.annot',...
        parms.fspath,hemi);
    end;
  else
    if ~iscell(parms.fnames_aparc), parms.fnames_aparc = {parms.fnames_aparc}; end;
    for f=1:length(parms.fnames_aparc)
      if mmil_isrelative(parms.fnames_aparc{h})
        parms.fnames_aparc{h} = [parms.fspath '/label/' parms.fnames_aparc{h}];
      end;
      fname = parms.fnames_aparc{f};
      n = regexp(fname,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
      if isempty(n)
        error('unexpected naming convention for aparc file %s\n',fname);
      end;
      parms.aparc_hemis{f} = n.hemi;
      parms.aparc_names{f} = n.name;
    end;
  end;
  parms.num_aparcs = length(parms.fnames_aparc);
  for f=1:parms.num_aparcs
    if ~exist(parms.fnames_aparc{f},'file')
      fprintf('%s: WARNING: aparc annot file %s not found\n',...
        mfilename,parms.fnames_aparc{f});
      parms.cortsurf_flag = 0;
      return;
    end;
  end;
  if ~isempty(parms.fnames_fparc)
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end;
    for f=1:length(parms.fnames_fparc)
      if ~exist(parms.fnames_fparc{f},'file')
        fprintf('%s: ERROR: fparc annot file %s not found\n',...
          mfilename,parms.fnames_fparc{f});
        errcode = 1;
        return;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = prep_data(parms)
  errcode = 0;
  % make local copy of parms.fname_T1
  fname_T1 = [parms.outpath '/T1_hires.mgz'];
  if ~exist(fname_T1,'file') || parms.forceflag
    fs_copy_mgh(parms.fname_T1,fname_T1);
  end;
  parms.fname_T1 = fname_T1;
  % convert BRIK files in inpath to mgz files in outpath
  for i=1:length(parms.input_files)
    fname_in = [parms.inpath '/' parms.input_files{i}];
    fname_out = [parms.outpath '/' parms.output_files{i}];
    mmil_BRIK2mgh(fname_in,fname_out,parms.forceflag);
  end;
  fname_T1_CBF = [parms.outpath '/T1_lowres.mgz'];
  parms.fname_T1_CBF = fname_T1_CBF;
  % register T1 to ASL
  parms = register_T1_to_ASL(parms);
  % resample each ASL measure to T1 space
  valid_flag = ones(1,parms.nmeas);
  parms.fnamelist = cell(1,parms.nmeas);
  for i=1:parms.nmeas
    meas = parms.measlist{i};
    fname = [parms.outpath '/' meas '.mgz'];
    fname = resample_data(fname,meas,parms);
    parms.fnamelist{i} = fname;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_fparc_annot(parms)
  i = parms.num_aparcs;
  for f=1:length(parms.fnames_fparc)
    fname_in = parms.fnames_fparc{f};
    if parms.verbose
      fprintf('%s: resampling annotation file %s from fsaverage to %s...\n',...
        mfilename,fname_in,parms.subj);
    end;
    % call fs_annot2annot
    tmp_parms = [];
    tmp_parms.outdir = [parms.outdir '/fparc_annot'];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj = parms.subj;
    tmp_parms.subjdir = parms.subjdir;
    tmp_parms.verbose = (parms.verbose==2);
    tmp_parms.forceflag = parms.forceflag;
    args = mmil_parms2args(tmp_parms);
    fname_out = fs_annot2annot(fname_in,args{:});
    n = regexp(fname_out,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
    if isempty(n)
      error('unexpected naming convention for aparc file %s\n',fname_out);
    end;
    i = i + 1;
    parms.fnames_aparc{i} = fname_out;
    parms.aparc_hemis{i} = n.hemi;
    parms.aparc_names{i} = n.name;
  end;
  parms.num_aparcs = length(parms.fnames_aparc);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = analyze_data(parms)
  errcode = 0;
  for i=1:parms.nmeas
    meas = parms.measlist{i};
    scalefact = parms.scalefacts(i);
    fname = parms.fnamelist{i};
    if parms.aseg_flag
      errcode = aseg_analysis(fname,meas,scalefact,parms);
      if errcode, parms.aseg_flag = 0; end;
    end;
    if parms.cortsurf_flag
      errcode = cortsurf_analysis(fname,meas,scalefact,parms);
      if errcode, parms.cortsurf_flag = 0; end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = atlas_warp_data(parms)
  errcode = 0;
  % register T1-weighted MRI to atlas and warp T1 to atlas
  if ~exist(parms.fname_T1,'file')
    fprintf('%s: WARNING: T1 file %s not found, skipping atlas warp\n',...
      mfilename,parms.fname_T1);
    return;
    errcode = 1;
  end;
  args = mmil_parms2args(parms,parms.warp_tags);
  [fname_atl,parms.fname_reg] = mmil_warp_to_atlas(parms.fname_T1,args{:});
  % warp each input volume to atlas
  for i=1:parms.nmeas
    meas = parms.measlist{i};
    fname = parms.fnamelist{i};
    parms.outstem = meas;
    args = mmil_parms2args(parms,parms.warp_tags);
    mmil_warp_to_atlas(fname,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = register_T1_to_ASL(parms)
  if parms.verbose
    fprintf('%s: registering T1 to ASL...\n',mfilename);
  end;
  
  % get info about T1 volume
  [parms.M_T1,parms.volsz_T1] = fs_read_header(parms.fname_T1);

  % create output directory for registration files
  mmil_mkdir(parms.regdir);

  % create brain mask from aseg using mmil_aseg_brainmask
  fname_mask = [parms.regdir '/aseg_mask.mgz'];
  tmp_parms = parms;
  tmp_parms.fname_mask = fname_mask;
  tmp_parms.brain_flag = 0;
  args = mmil_parms2args(tmp_parms,parms.mask_tags);
  mmil_aseg_brainmask(parms.fspath,args{:});

  % use mmil_reg to register lowres T1 to highres T1
  tmp_parms = parms;
  tmp_parms.fname_maskA = fname_mask;
  tmp_parms.outdir = parms.regdir;
  % tmp_parms.fname_T1_orig = parms.fname_T1;
  args = mmil_parms2args(tmp_parms,parms.reg_tags);
  parms.M_T1_to_ASL = mmil_reg(parms.fname_T1,parms.fname_T1_CBF,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = resample_data(fname,meas,parms)
  % resample data to T1 resolution using M_T1_to_ASL
  if ~exist('resT1_flag','var'), resT1_flag = 0; end;
  nvox_ref = parms.volsz_T1;
  M_ref = parms.M_T1;
  M_reg = parms.M_T1_to_ASL;
  fname_out = [parms.outdir '/' meas '_' parms.resT1_outfix '.mgz'];
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: resampling %s to T1 resolution...\n',...
        mfilename,fname);
    end;
    [vol,M] = fs_load_mgh(fname,[],1);
    [vol,M] = mmil_resample_vol(vol,M,...
      'nvox_ref',nvox_ref,'M_ref',M_ref,'M_reg',M_reg);
    fs_save_mgh(vol,fname_out,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = aseg_analysis(fname,meas,scalefact,parms)
  errcode = 0;
  % calculate averages for aseg ROIs
  if parms.verbose
    fprintf('%s: aseg analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' meas];
  parms.scalefact = scalefact;
  args = mmil_parms2args(parms,parms.aseg_tags);
  try
    mmil_aseg_analysis(fname,parms.fspath,args{:});
  catch
    fprintf('\n%s: WARNING: aseg analysis failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = cortsurf_analysis(fname,meas,scalefact,parms)
  errcode = 0;
  % calculate averages for cortsurf ROIs
  if parms.verbose
    fprintf('%s: cortsurf analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' meas '_cortsurf'];
  parms.outfix = 'cortsurf';
  parms.scalefact = scalefact;
  parms.hemilist = parms.aparc_hemis;
  args = mmil_parms2args(parms,parms.cortsurf_tags);
  try
    mmil_cortsurf_analysis(fname,parms.fspath,args{:});
  catch
    fprintf('\n%s: WARNING: cortsurf analysis failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

