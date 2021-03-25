function mmil_gwcsurf_analysis(fname,fspath,varargin)
%function mmil_gwcsurf_analysis(fname,fspath,[options])
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
%   'pmin': minimum distance (mm) to project along surface vertex normal
%     {default = 0}
%   'pmax': maximum distance (mm) to project along surface vertex normal
%     {default = 5}
%   'pstep': step distance (mm) along surface vertex normal
%     {default = 1}
%   'interpmethod': interpolation method ('nearest','linear','spline','cubic')
%     {default = 'linear'}
%   'avg_flag': [0|1] average values across projection steps
%     otherwise create multi-frame output volume
%     {default = 1}
%   'gwmask_flag': [0|1] use cortical ribbon volume for gray/white matter masks
%     to calculate weighted average
%     {default = 0}
%   'tukey_flag': use Tukey's bisquare function to transform weights
%     {default = 1}
%   'tukey_fact': scaling factor used to transform weights with Tukey's bisqaure
%     {default = 0.5}
%   'gwnorm_flag': [0|1] for gray/white contrast, normalize difference by mean
%     {default = 1}
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1] whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere (can be vector)
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
% Created:  09/03/15 by Don Hagler
% Last Mod: 09/04/15 by Don Hagler
%

%% todo: option to analyze only gm, or only wm?

% based on mmil_cortsurf_analysis, created 10/12/12 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% check parameters, check input files and directories
parms = check_input(fname,fspath,varargin);

% create wm and gm masks
if parms.gwmask_flag
  parms = create_gwmasks(parms);
end;

% resample to T1 resolution
if parms.resT1flag && ~isempty(parms.M_reg)
  parms = resample_data(parms);
end;

% sample volume data onto cortical surface mesh
paint_gw_surf(parms);

% calculate gray-white contrast maps
calc_contrast_maps(parms);

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
    'pmin',0,[0,10],...
    'pmax',5,[0,10],...
    'pstep',1,[0.001,10],...
    'interpmethod','linear',{'nearest','linear','spline','cubic'},...
    'avg_flag',true,[false true],...
    'gwmask_flag',false,[false true],...
    'tukey_flag',true,[false true],...
    'tukey_fact',0.5,[1e-10,1e10],...
    'gwnorm_flag',true,[false true],...
    'smoothsteps',0,[0,Inf],...
    'sphere_flag',false,[false true],...
    'sphsmoothsteps',0,[0,Inf],...
    'outtype','mgz',{'mgh','mgz'},...
    'scalefact',1,[1e-10,1e10],...
    'minval',1e-6,[],...
    'fname_colorlut',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
... % hidden parameters
    'outfix',[],[],...
    'contrast_outfix','contrast',[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'frames',1,[],...
    'interpm',2,[0:5],...
    'bclamp',false,[false true],...
    'wm_codes',[2,41],[],...
    'gm_codes',[3,42],[],...
    'wm_outstem','wmmask',[],...
    'gm_outstem','gmmask',[],...
    'wm_outfix','wm',[],...
    'gm_outfix','gm',[],...
    'gwc_outfix','gwc',[],...
    'cortex_flag',true,[false true],...
    'icolevel',7,[0 7],...
    'surfname','white',[],...
...
    'paint_tags',{'hemi','outstem','outdir','outfix','outext',...
                  'regmat','pmin','pmax','pstep','interpmethod',...
                  'avg_flag','fname_weights','tukey_flag','tukey_fact',...
                  'subjdir','surfname','frame','forceflag','verbose','hemilist'},[],...
    'surf2surf_tags',{'trgsubj','fname_out','outdir','hemi','intype',...
                      'outtype','subjdir','smooth_in','smooth_out',...
                      'cortex_flag','icolevel','options','forceflag',...
                      'verbose'},[],...
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

  % check cortical ribbon file
  if parms.gwmask_flag
    parms.fname_ribbon = sprintf('%s/mri/ribbon.mgz',parms.fspath);
    if ~exist(parms.fname_ribbon,'file')
      error('file %s not found',parms.fname_ribbon);
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

  if ~parms.sphere_flag, parms.sphsmoothsteps = 0; end;

  % set M_reg if empty and resT1flag = 1
  if parms.resT1flag && isempty(parms.M_reg)
    parms.M_reg = eye(4);
  end;
  parms.regmat = parms.M_reg;

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

function parms = create_gwmasks(parms)
  % create gray matter mask from ribbon file
  parms.fname_gmmask = ...
    create_mask_from_ribbon(parms,parms.gm_codes,parms.gm_outstem);
  % create white matter mask from ribbon file
  parms.fname_wmmask = ...
    create_mask_from_ribbon(parms,parms.wm_codes,parms.wm_outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_mask = create_mask_from_ribbon(parms,codes,outstem)
  fname_mask = [];
  % create mask in FS resolution
  fname_mask = sprintf('%s/%s.mgz',parms.outdir,outstem);
  if ~exist(fname_mask,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: creating %s from ribbon file...\n',mfilename,outstem);
    end;
  end;
  fs_aseg_mask(parms.fname_ribbon,'fname_mask',fname_mask,...
    'aseg_codes',codes,'forceflag',parms.forceflag);
  % resample mask to input data resolution
  if ~isempty(parms.M_reg) && ~parms.resT1flag
    fname_res = sprintf('%s/%s_res.mgz',parms.outdir,outstem);
    if ~exist(fname_res,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: resampling %s to input data resolution...\n',...
          mfilename,outstem);
      end;
      % load M and volsz from input file
      [M_ref,nvox_ref] = mmil_load_mgh_info(parms.fname,...
                                            parms.forceflag,parms.outdir);
      % load mask
      [vol_mask,M] = fs_load_mgh(fname_mask);
      % resample mask to input data resolution
      [vol_res,M_res] = mmil_resample_vol(vol_mask,M,...
        'nvox_ref',nvox_ref,'M_ref',M_ref,...
        'interpm',parms.interpm,'bclamp',parms.bclamp,...
        'M_reg',inv(parms.M_reg));
      % clip resampled mask values to be between 0 and 1
      vol_res = max(min(vol_res,1),0);
      % save resampled mask
      fs_save_mgh(vol_res,fname_res,M_res);
    end;
    fname_mask = fname_res;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_data(parms)
  % resample input data to FS resolution
  fname_res = [parms.outstem '_' parms.res_outfix '.mgz'];
  if ~exist(fname_res,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: resampling %s to T1 resolution...\n',...
        mfilename,parms.fname);
    end;
    fname_orig = sprintf('%s/mri/orig.mgz',parms.fspath);
    [M_ref,nvox_ref] = mmil_load_mgh_info(fname_orig,...
                                          parms.forceflag,parms.outdir);
    [vol,M] = fs_load_mgh(parms.fname,[],parms.frames);
    [vol,M] = mmil_resample_vol(vol,M,...
      'nvox_ref',nvox_ref,'M_ref',M_ref,...
      'interpm',parms.interpm,'bclamp',parms.bclamp,...
      'M_reg',parms.M_reg);
    fs_save_mgh(vol,fname_res,M);
  end;
  parms.fname = fname_res;    
  parms.frames = [];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paint_gw_surf(parms)
  % sample values for gray matter onto surface
  paint_surf(parms,parms.gm_outfix);
  % sample values for white matter onto surface
  paint_surf(parms,parms.wm_outfix);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paint_surf(parms,outfix)
  % sample values onto surface
  tmp_parms = parms;
  tmp_parms.outstem = parms.outstem;
  tmp_parms.outfix = ['_' outfix];
  tmp_parms.outext = ['.' parms.outtype];
  switch outfix
    case 'gm'
      if parms.gwmask_flag
        tmp_parms.fname_weights = parms.fname_gmmask;
      end;      
    case 'wm'
      if parms.gwmask_flag
        tmp_parms.fname_weights = parms.fname_wmmask;
      end;      
      tmp_parms.pmax = -parms.pmin;
      tmp_parms.pmin = -parms.pmax;
  end;
  if isempty(parms.frames)
    tmp_parms.frame = 1;
  else
    tmp_parms.frame = parms.frames(1);
  end;
  args = mmil_parms2args(tmp_parms,parms.paint_tags);
  mmil_paint(parms.subj,parms.fname,args{:});

  % smooth values on surface, optionally resample to sphere
  if parms.smoothsteps || parms.sphere_flag
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      smooth_vals(outfix,hemi,parms);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% smooth values on surface and/or resample to sphere
function fname_out = smooth_vals(outfix,hemi,parms)
  for sphsmoothsteps = parms.sphsmoothsteps
    fname_in = set_fname(parms,outfix,0,0,0,hemi);
    tmp_parms = [];
    tmp_parms.fname_out = set_fname(parms,outfix,...
      parms.smoothsteps,parms.sphere_flag,sphsmoothsteps,hemi);
    tmp_parms.smooth_in = parms.smoothsteps;
    if parms.sphere_flag
      tmp_parms.trgsubj = 'ico';
      tmp_parms.smooth_out = sphsmoothsteps;
    end;
    args = mmil_parms2args(tmp_parms,parms.surf2surf_tags);
    fnames_out = fs_surf2surf(fname_in,parms.subj,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_fname(parms,outfix,smoothsteps,sphere_flag,sphsmoothsteps,hemi)
  fname = sprintf('%s_%s',parms.outstem,outfix);
  if smoothsteps
    fname = sprintf('%s-sm%d',fname,smoothsteps);
  end;
  if sphere_flag
    fname = [fname '-sphere'];
    if sphsmoothsteps
      fname = sprintf('%s-sm%d',fname,sphsmoothsteps);
    end;
  end;
  fname = sprintf('%s-%s.%s',fname,hemi,parms.outtype);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maps of gray-white contrast
function calc_contrast_maps(parms)
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    for sphsmoothsteps = parms.sphsmoothsteps
      fname_gm = set_fname(parms,'gm',...
        parms.smoothsteps,parms.sphere_flag,sphsmoothsteps,hemi);
      fname_wm = set_fname(parms,'wm',...
        parms.smoothsteps,parms.sphere_flag,sphsmoothsteps,hemi);
      fname_gwc = set_fname(parms,'gwc',...
        parms.smoothsteps,parms.sphere_flag,sphsmoothsteps,hemi);
      calc_contrast_map(fname_gm,fname_wm,fname_gwc,parms);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate contrast for two input files
function calc_contrast_map(fname_gm,fname_wm,fname_gwc,parms)
  if ~exist(fname_gwc,'file') || parms.forceflag
    % load values
    vals_gray = fs_load_mgh(fname_gm);
    vals_white = fs_load_mgh(fname_wm);
    % calculate difference between gray and white matter values
    vals_contrast = vals_white - vals_gray;
    % normalize by mean
    if parms.gwnorm_flag
      vals_contrast = 2*vals_contrast./(vals_white + vals_gray);
    end;
    fs_save_mgh(vals_contrast,fname_gwc);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extract_all_values(parms)
  % extract ROI values for gm and wm
  if ~isempty(parms.outfix)
    parms.outfix = ['_' parms.outfix];
  end;
  infix_list = {'gm','wm'};
  for i=1:length(infix_list)
    infix = infix_list{i};
    instem = sprintf('%s_%s',parms.outstem,infix);
    outstem = sprintf('%s%s_roi_data',instem,parms.outfix);
    fname_mat = sprintf('%s.mat',outstem);
    extract_values(instem,fname_mat,parms.scalefact,parms);
    if parms.csv_flag
      fname_csv = sprintf('%s.csv',outstem);
      write_csv(fname_mat,fname_csv,parms);
    end;
    switch infix
      case 'gm'
        fname_gm = fname_mat;
      case 'wm'
        fname_wm = fname_mat;      
    end;
  end;

  % calculate gray-white contrast
  outstem = sprintf('%s_%s%s_roi_data',...
    parms.outstem,parms.gwc_outfix,parms.outfix);
  fname_gwc = sprintf('%s.mat',outstem);
  calculate_contrast_roi(fname_gm,fname_wm,fname_gwc,parms)
  if parms.csv_flag
    fname_csv = sprintf('%s.csv',outstem);
    write_csv(fname_gwc,fname_csv,parms);
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
      tmp_parms.frames = [];
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

function calculate_contrast_roi(fname_gm,fname_wm,fname_gwc,parms)
  if ~exist(fname_gwc,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating gray-white contrast for ROIs...\n',mfilename);
    end;
    % load values
    gm = load(fname_gm);
    wm = load(fname_wm);
    % calculate difference between gray and white matter values
    roi_data = wm.roi_data;
    for i=1:length(roi_data)
      roi_data(i).vals = wm.roi_data(i).vals - gm.roi_data(i).vals;
      roi_data(i).median = wm.roi_data(i).median - gm.roi_data(i).median;
      roi_data(i).avg = wm.roi_data(i).avg - gm.roi_data(i).avg;
      roi_data(i).stdv = sqrt(wm.roi_data(i).stdv.^2 + gm.roi_data(i).stdv.^2);
      % normalize by mean
      if parms.gwnorm_flag
        mean_vals = (wm.roi_data(i).vals + gm.roi_data(i).vals)/2;
        roi_data(i).vals = roi_data(i).vals ./ mean_vals;
        mean_val = (wm.roi_data(i).avg + gm.roi_data(i).avg)/2;
        roi_data(i).avg = roi_data(i).avg ./ mean_val;
        roi_data(i).stdv = roi_data(i).stdv ./ mean_val;
        mean_val = (wm.roi_data(i).median + gm.roi_data(i).median)/2;
        roi_data(i).median = roi_data(i).median ./ mean_val;
      end;
    end;
    save(fname_gwc,'roi_data');
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

