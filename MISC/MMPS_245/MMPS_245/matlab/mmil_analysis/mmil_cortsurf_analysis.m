function mmil_cortsurf_analysis(fname,fspath,varargin)
%function mmil_cortsurf_analysis(fname,fspath,[options])
%
% Required Parameters:
%   fname: full path name of file containing image data
%   fspath: full path of freesurfer recon
%
% Optional Parameters:
%   'outdir': output directory
%     if empty, will attempt to use path of fname
%     {default = []}
%   'outstem': output file stem
%     if empty, will use file stem of fname
%     {default = []}
%   'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     if empty, will use ?h.aparc.annot files in fspath/label
%     if not full path, assumed to be relative to fspath/label
%     {default = []}
%   'fnames_weights': cell array of weighted ROI mgh files (one each hemi)
%     may be multi-frame for multiple ROIs
%     {default = []}
%   'weights_thresh': threshold applied to weights file
%     {default = 0}
%   'csv_flag': [0|1] summarize results in csv files
%     {default = 1}
%   'M_reg': 4x4 transformation matrix between fname and fspath/mri/orig.mgz
%     should be like M_orig_to_fname
%     if empty, will assume fname volume is registered to orig.mgz
%     {default = []}
%   'resT1flag': [0|1] resample volume to T1 resolution before painting
%     {default = 0}
%   'res_outfix': output string added to resampled fname
%     ignored if resT1flag = 0
%     {default = 'resT1'}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast maps
%     {default = [-1,1]}
%   'gwnorm_flag': [0|1] for gray/white contrast, normalize difference by mean
%     {default = 1}
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1] whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere (can be vector)
%     {default = 0}
%   'mask_midbrain_flag', [0|1] whether to mask out mid brain and other
%     cortical regions marked "unknown" (masking done before smoothing)
%     {default = 0}
%   'outtype': output file type ('mgh', or 'mgz')
%    {default = 'mgz'}
%   'scalefact': scaling factor applied to extracted ROI values
%     {default = 1}
%   'minval': minimum value; if 0, include all voxels except NaNs
%     {default = 1e-6}
%   'fname_colorlut': full path of color look up table file (to get ROI names)
%     if empty, use FreeSurferColorLUT.txt in $FREESURFER_HOME
%    {default = []}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% NOTE: surface maps resulting from smoothing and sampling to sphere
%   are intended for vertex-wise group-analysis and are
%   not used for extracting ROI measures
%
% Created:  10/12/12 by Don Hagler
% Last Mod: 09/22/15 by Don Hagler
%

% based on DTI_MMIL_Analyze_CortSurf_Exam, created 10/29/08 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% check parameters, check input files and directories
parms = check_input(fname,fspath,varargin);

% sample volume data onto cortical surface mesh
parms = paint_surf(parms);

% calculate gray-white contrast maps
if parms.contrast_flag, parms = calc_contrast_maps(parms); end;

% calculate average values for cortical parcellation ROIs
extract_all_values(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,fspath,options)
  % parse input arguments
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
    'fspath',fspath,[],...
...
    'outdir',[],[],...
    'outstem',[],[],...
    'fnames_aparc',[],[],...
    'fnames_weights',[],[],...
    'weights_thresh',0,[0,Inf],...
    'csv_flag',true,[false true],...
    'M_reg',[],[],...
    'resT1flag',false,[false true],...
    'res_outfix','resT1',[],...
    'projdist_list',[-1,1],[-5,5],...
    'gwnorm_flag',true,[false true],...
    'smoothsteps',0,[0,Inf],...
    'sphere_flag',false,[false true],...
    'sphsmoothsteps',0,[0,Inf],...
    'mask_midbrain_flag',false,[false true],...
    'outtype','mgz',{'mgh','mgz'},...
    'scalefact',1,[1e-10,1e10],...
    'minval',1e-6,[],...
    'fname_colorlut',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
... % hidden parameters
    'mbmask_flag',false,[false true],...
    'outfix',[],[],...
    'contrast_outfix','contrast',[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'frames',1,[],...
    'interpm',2,[0:5],...
    'bclamp',false,[false true],...
...
    'paint_tags',{'subjdir','sphere_flag','smoothsteps',...
                  'sphsmoothsteps','mask_midbrain_flag',...
                  'outstem','outtype','regfile',...
                  'projdist','mask_midbrain_flag',...
                  'verbose','forceflag'},[],...
    'roi_tags',{'fname_aparc','fname_label','fname_weights','frames',...
                'minval','scalefact','fname_colorlut','hemi',...
                'weights_thresh','verbose'},[],...
  });

  % check input file
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;

  % check FreeSurfer recon exists
  if ~exist(parms.fspath,'dir')
    error('FreeSurfer recon dir %s not found',parms.fspath);
  end;
  [parms.subjdir,parms.subj,text] = fileparts(parms.fspath);
  parms.subj = [parms.subj text];

  % check fnames_aparc
  parms.nhemi = length(parms.hemilist);
  if isempty(parms.fnames_aparc)
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      parms.fnames_aparc{h} = sprintf('%s/label/%s.aparc.annot',...
        parms.fspath,hemi);
    end;
  else
    if ~iscell(parms.fnames_aparc)
      parms.fnames_aparc = {parms.fnames_aparc};
    end;
    if length(parms.fnames_aparc) ~= parms.nhemi
      error('must have %d elements in fnames_aparc (have %d)',...
        parms.nhemi,length(parms.fnames_aparc));
    end;
  end;
  for h=1:parms.nhemi
    if mmil_isrelative(parms.fnames_aparc{h})
      parms.fnames_aparc{h} = [parms.fspath '/label/' parms.fnames_aparc{h}];
    end;
    if ~exist(parms.fnames_aparc{h},'file')
      error('file %s not found',parms.fnames_aparc{h});
    end;
  end;

  % check fnames_weights
  if ~isempty(parms.fnames_weights)
    if ~isempty(parms.fnames_weights{h}) &&...
       ~exist(parms.fnames_weights{h},'file')
      error('file %s not found',parms.fnames_weights{h});
    end;
  else
    parms.fnames_weights = cell(parms.nhemi,1);
  end;

  % set contrast_flag depending on number of elements in projdist_list
  parms.nprojdist = length(parms.projdist_list);
  if parms.nprojdist==2
    parms.contrast_flag = 1;
    parms.projdist_list = sort(parms.projdist_list); % white matter first
  else
    parms.contrast_flag = 0;
  end;
  if ~parms.sphere_flag, parms.sphsmoothsteps = 0; end;

  % reset parameter name for fs_paint
  parms.mask_midbrain_flag = parms.mbmask_flag;

  % set M_reg if empty and resT1flag = 1
  if parms.resT1flag && isempty(parms.M_reg)
    parms.M_reg = eye(4);
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
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = paint_surf(parms)
  fname_regdat = [];
  if ~isempty(parms.M_reg)
    fname_orig = sprintf('%s/mri/orig.mgz',parms.fspath);
    [M_ref,nvox_ref] = mmil_load_mgh_info(fname_orig,...
                                          parms.forceflag,parms.outdir);
    if parms.resT1flag
      % resample to T1 resolution
      fname_res = [parms.outstem '_' parms.res_outfix '.mgz'];
      if ~exist(fname_res,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: resampling %s to T1 resolution...\n',...
            mfilename,parms.fname);
        end;
        [vol,M] = fs_load_mgh(parms.fname,[],parms.frames);
        [vol,M] = mmil_resample_vol(vol,M,...
          'nvox_ref',nvox_ref,'M_ref',M_ref,...
          'interpm',parms.interpm,'bclamp',parms.bclamp,...
          'M_reg',parms.M_reg);
        fs_save_mgh(vol,fname_res,M);
      end;
      parms.fname = fname_res;    
      parms.frames = [];
    else
      %  create register.dat file from M_reg
      [M_reg,nvox_reg] = mmil_load_mgh_info(parms.fname,...
                                            parms.forceflag,parms.outdir);
      fname_regdat =  sprintf('%s_register.dat',parms.outstem);
      fs_write_regdat(fname_regdat,'M',parms.M_reg,...
        'M_ref',M_ref,'nvox_ref',nvox_ref,...
        'M_reg',M_reg,'nvox_reg',nvox_reg,...
        'ras2tk_flag',1,'forceflag',parms.forceflag);
    end;
  end;
  for projdist = parms.projdist_list
    tmp_parms = parms;
    tmp_parms.projdist = projdist;
    tmp_parms.outstem = sprintf('%s_pdist%0.1f',parms.outstem,projdist);
    tmp_parms.regfile = fname_regdat;
    for sphsmoothsteps = parms.sphsmoothsteps
      tmp_parms.sphsmoothsteps = sphsmoothsteps;
      args = mmil_parms2args(tmp_parms,parms.paint_tags);
      fs_paint(parms.subj,parms.fname,args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maps of gray-white contrast
function parms = calc_contrast_maps(parms)
  % calculate difference between 1st and 2nd element of projdist_list
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    for sphsmoothsteps = parms.sphsmoothsteps
      infix  = [];
      if parms.mask_midbrain_flag
        infix = sprintf('%s-mbmask',infix);
      end;
      if parms.smoothsteps
        infix = sprintf('%s-sm%d',infix,parms.smoothsteps);
      end;
      if parms.sphere_flag
        infix = [infix '-sphere'];
      end;
      if sphsmoothsteps
        infix = sprintf('%s-sm%d',infix,sphsmoothsteps);
      end;
      fname_gray = sprintf('%s_pdist%0.1f%s-%s.%s',...
        parms.outstem,parms.projdist_list(2),infix,hemi,parms.outtype);
      fname_white = sprintf('%s_pdist%0.1f%s-%s.%s',...
        parms.outstem,parms.projdist_list(1),infix,hemi,parms.outtype);
      fname_contrast = sprintf('%s_%s_pdist%0.1f%s-%s.%s',...
        parms.outstem,parms.contrast_outfix,parms.projdist_list(2),...
        infix,hemi,parms.outtype);
      calc_contrast_map(fname_gray,fname_white,fname_contrast,parms);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate contrast for two input files
function calc_contrast_map(fname_gray,fname_white,fname_contrast,parms)
  if ~exist(fname_contrast,'file') || parms.forceflag
    % load values
    vals_gray = fs_load_mgh(fname_gray);
    vals_white = fs_load_mgh(fname_white);
    % calculate difference between gray and white matter values
    vals_contrast = vals_white - vals_gray;
    % normalize by mean
    if parms.gwnorm_flag
      vals_contrast = 2*vals_contrast./(vals_white + vals_gray);
    end;
    fs_save_mgh(vals_contrast,fname_contrast);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extract_all_values(parms)
  % extract ROI values for each element of projdist_list
  if ~isempty(parms.outfix)
    parms.outfix = [parms.outfix '_'];
  end;
  fnamelist = cell(parms.nprojdist,1);
  for p=1:parms.nprojdist
    projdist = parms.projdist_list(p);
    fstem = sprintf('%s_pdist%0.1f',parms.outstem,projdist);
    outfix = sprintf('%spdist%0.1f_roi_data',parms.outfix,projdist);
    fname_mat = [parms.outstem '_' outfix '.mat'];
    fname_csv = [parms.outstem '_' outfix '.csv'];
    extract_values(fstem,fname_mat,parms.scalefact,parms);
    if parms.csv_flag
      write_csv(fname_mat,fname_csv,parms);
    end;
    fnamelist{p} = fname_mat;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extract_values(fstem,fname_mat,scalefact,parms)
  if ~exist(fname_mat,'file') || parms.forceflag
    roi_data = [];
    for h=1:parms.nhemi
      tmp_parms = parms;
      hemi = parms.hemilist{h};
      fname_vals = [fstem '-' hemi '.' parms.outtype];
      tmp_parms.fname_aparc = parms.fnames_aparc{h};
      tmp_parms.fname_weights = parms.fnames_weights{h};
      if parms.verbose
        fprintf('%s: extracting aparc ROI data from %s...\n',...
          mfilename,fname_vals);
      end;
      tmp_parms.scalefact = scalefact;
      args = mmil_parms2args(tmp_parms,parms.roi_tags);
      tmp_data = mmil_surf_roi(fname_vals,args{:});
      if isempty(tmp_data)
        error('failed to get %s cortsurf ROI data',hemi);
      end;
      roi_data = [roi_data,tmp_data];
    end;
    save(fname_mat,'roi_data');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(fname_mat,fname_csv,parms)
  if ~exist(fname_csv,'file') || parms.forceflag
    roi_data = [];
    load(fname_mat);
    fid = fopen(fname_csv,'wt');
    if fid==-1, error('failed to open %s for writing',fname_csv); end;
    fprintf(fid,'"ROI","mean","median","stdv","nvals","nvals valid"\n');
    nrois = length(roi_data);
    for i=1:nrois
      fprintf(fid,'"%s",',roi_data(i).roiname);
      fprintf(fid,'%0.6f,%0.6f,%0.6f,%d,%d\n',...
        roi_data(i).avg,roi_data(i).median,roi_data(i).stdv,...
        roi_data(i).nvals_valid,roi_data(i).nvals_valid);
    end
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

