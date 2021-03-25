function mmil_fparc2fibers(vol,M,qmat,bvals,subj,varargin)
%function mmil_fparc2fibers(vol,M,qmat,bvals,subj,[options])
%
% Required parameters:
%   vol: 4D volume containing multiple diffusion weighted volumes
%   M: 4x4 vox2ras matrix for vol
%   qmat: matrix of diffusion direction vectors
%   bvals: vector of b values
%     one for all, or one for each diffusion direction
%   subj: FreeSurfer subject name
%
% Optional parameters:
%  'M_reg': registration matrix between FreeSurfer recon and data
%    if empty, assumes data have been resampled to FreeSurfer space
%    {default = []}
%  'annotname': name of annotation
%    {default = 'fparc'}
%  'wmparcname' : wmparc file stem
%    {default = 'wmparc'}
%  'asegname': aseg file stem
%    {default = 'aseg'}
%  'outdir': output directory; if empty, will use pwd
%    {default = []}
%  'fnames_fparc': cell array of annotation files in fsaverage space
%    will be resampled to individual subject space before use
%    if supplied, fnames_aparc will be ignored
%    {default = []}
%  'fnames_aparc': cell array of annotation files (one for each hemisphere)
%    if empty, will use ?h.aparc.annot files in fspath/label
%    if not full path, assumed to be relative to fspath/label
%    {default = []}
%  'fname_aseg': name of aseg file
%    if empty, will use asegname
%    if not full path, assumed to be relative to fspath/mri
%    {default = []}
%  'fname_aseg_clut': name of color LUT file for aseg
%    if empty, will use $FREESURFER_HOME/FreeSurferColorLUT.txt
%    {default = []}
%  'subjdir': directory containing FreeSurfer subject
%    {default = [getenv('FREESURFER_HOME') '/subjects']}
%  'fname_roi_info': full path of csv file containing ROI codes and names
%    if supplied, will determine which ROIs are used for fiber selection
%    {default = []}
%  'fname_fiber_rois': full path of csv file containing fiber names
%    and matching ROI names used to select the fiber
%    if not supplied, will generate all combinations of ipsilateral
%     connections and all connections between contralateral homologs
%    {default = []}
%  'contra_homo_flag': [0|1] whether to select fibers for connections
%     between contralateral homologs
%    {default = 1}
%  'not_contra_wm_flag': [0|1] use aseg to exclude streamlines that pass
%     through contralateral white matter (except for contra_homo fibers)
%    {default = 1}
%  'cutflag': [0|1] trim streamlines as they exit terminal ROIs
%    {default = 0}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%   
%
% Created:  09/10/15 by Don Hagler
% Last Mod: 10/19/15 by Don Hagler
%

%% todo: allow use of AtlasTrack ROIs
%          e.g. AND AtlasTrack fiber
%   see: rt_AtlasTrack_select_fibers.m, run_FACT_aseg_wmparc.m

%% todo: create ROIs through combinations of ROIs
%%       to allow an implicit OR
%%       while keeping things simple (e.g. streamlines between two ROIs)

%% todo: allow fname_wmparc as input
%%       also require matching colorLUT file

%% todo: fiber outstem?

if ~mmil_check_nargs(nargin,5), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(vol,M,qmat,bvals,subj,varargin);

% create wmparc from fparc
parms = create_wmparc(parms);

% resample wmparc and aseg to DTI space
if ~isempty(parms.M_reg)
  parms = resample_wmparc(parms);
  parms = resample_aseg(parms);
end;

% create individual ROI mgz from wmparc and aseg
parms = create_rois(parms);

% generate spreadsheets with ROI names for each fiber
parms = create_fiber_rois(parms);

% generate fiber tracts
parms = generate_streamlines(parms);

% select fiber tracts based on spreadsheet
parms = select_fibers(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
function parms = check_input(vol,M,qmat,bvals,subj,options)
  parms = mmil_args2parms(options,{...
    'vol',vol,[],...
    'M',M,[],...
    'qmat',qmat,[],...
    'bvals',bvals,[],...
    'subj',subj,[],...
  ...
    'M_reg',[],[],...
    'annotname','fparc',[],...
    'wmparcname','wmparc',[],...
    'asegname','aseg',[],...
    'outdir',[],[],...
    'fnames_fparc',[],[],...
    'fnames_aparc',[],[],...
    'fname_aseg',[],[],...
    'fname_aseg_clut',[getenv('FREESURFER_HOME') '/FreeSurferColorLUT.txt'],[],...
    'subjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
    'fname_roi_info',[],[],...
    'fname_fiber_rois',[],[],...
    'contra_homo_flag',true,[false true],...
    'not_contra_wm_flag',true,[false true],...
    'cutflag',false,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'wmparc_tags',{'annotname' 'wmparcname' 'outdir' 'fnames_fparc'...
                   'fnames_aparc' 'subjdir' 'verbose' 'forceflag'...
                   'tmpdir' 'cleanup_flag' 'verbose' 'hemilist'...
                   'gm_codes' 'wm_codes' 'extra_lut_lines'},[],...
    'tracto_tags',{'outdir' 'outstem' 'orient' 'forceflag'...
                   'seed_point_sampling' 'step_size' 'FOD_thresh'...
                   'angle_thresh' 'fiber_length_range' 'max_order'...
                   'thresh_FA' 'smf' 'orient_ref'},[],...
    'select_tags',{'fiber' 'fname_out' 'fname_NOT' 'cutflag' 'revsliceflag'...
                   'thresh' 'verbose' 'forceflag'},[],...
  ...
    'verbose',false,[false true],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'gm_codes',[1000,2000],[],...
    'wm_codes',[3000,4000],[],...
    'sc_codes',[10,11,12,49,50,51],[],... % thalamus, caudate, putamen
    'sc_flag',true,[false true],... % whether to include sub-cortical fibers
    'gm_flag',2,[0:2],... % 0=wm, 1=gm, 2=wm+gm
    'aseg_wm_codes',[2,41],[],...
    'aseg_codes',[1:255],[],...
    'extra_lut_lines',[],[],...
    'roi_outfix',[],[],...
    'roi_outtype','mat',{'mat','mgh','mgz'},...
    'roi_dilate_flag',false,[false true],...
  });

  parms.nhemi = length(parms.hemilist);

  if isempty(parms.subjdir)
    parms.subjdir = getenv('SUBJECTS_DIR');
    if isempty(parms.subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  else
    setenv('SUBJECTS_DIR',parms.subjdir);
  end;

  % set fspath
  parms.fspath = sprintf('%s/%s',parms.subjdir,parms.subj);
  if ~exist(parms.fspath,'file')
    error('FS recon dir %s not found',parms.fspath);
  end;

  % check fparc file
  if ~isempty(parms.fnames_fparc)
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end;
    for f=1:length(parms.fnames_fparc)
      if ~exist(parms.fnames_fparc{f},'file')
        error('fparc annot file %s not found',parms.fnames_fparc{f});
      end;
    end;
  end;

  % check aparc file if no fparc file
  if isempty(parms.fnames_fparc)
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
      hemi = parms.hemilist{h};
      if mmil_isrelative(parms.fnames_aparc{h})
        parms.fnames_aparc{h} = [parms.fspath '/label/' parms.fnames_aparc{h}];
      end;
      if ~exist(parms.fnames_aparc{h},'file')
        error('file %s not found',parms.fnames_aparc{h});
      end;
      parms.aparc_hemis{h} = hemi;
      parms.aparc_names{h} = 'aparc';
    end;
  else
    parms.fnames_aparc = [];
  end;

  % check aseg file
  if isempty(parms.fname_aseg)
    parms.fname_aseg = [parms.asegname '.mgz'];
  end;
  if mmil_isrelative(parms.fname_aseg)
    parms.fname_aseg = sprintf('%s/mri/%s',parms.fspath,parms.fname_aseg);
  end;
  if ~exist(parms.fname_aseg,'file')
    error('aseg file %s not found',parms.fname_aseg);
  end;

  if ~isempty(parms.M_reg)
    parms.roi_outfix = '_resDTI';
  else
    fprintf('%s: NOTE: no M_reg supplied... assuming data are resampled to FreeSurfer recon\n',...
      mfilename);
  end;

  % set subdirectory for ROI files
  parms.roi_outdir = [parms.outdir '/rois'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create wmparc from fparc
function parms = create_wmparc(parms)
  parms.wmparc_outdir = [parms.outdir '/wmparc'];
  tparms = parms;
  tparms.outdir = parms.wmparc_outdir;
  tparms.tmpdir = 'fparc2wmparc';
  tparms.cleanup_flag = 0;
  args = mmil_parms2args(tparms,parms.wmparc_tags);
  [parms.fname_wmparc,parms.fname_clut] = mmil_fparc2wmparc(parms.subj,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_wmparc(parms)
  % resample wmparc to DTI space
  fname_wmparc_res = sprintf('%s/%s%s.mgz',...
    parms.wmparc_outdir,parms.wmparcname,parms.roi_outfix);
  if ~exist(fname_wmparc_res,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: resampling wmparc...\n',mfilename);
    end;
    [vol_wmparc,M_wmparc] = fs_load_mgh(parms.fname_wmparc);
    [vol_wmparc_res,M_res] = mmil_resample_vol(vol_wmparc,M_wmparc,...
      'nvox_ref',size(parms.vol),'M_ref',parms.M,...
      'interpm',0,'M_reg',inv(parms.M_reg));
    fs_save_mgh(vol_wmparc_res,fname_wmparc_res,M_res);
  end;
  parms.fname_wmparc = fname_wmparc_res;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_aseg(parms)
  % resample aseg to DTI space
  parms.aseg_outdir = [parms.outdir '/aseg'];
  mmil_mkdir(parms.aseg_outdir);
  fname_aseg_res = sprintf('%s/%s%s.mgz',...
    parms.aseg_outdir,parms.asegname,parms.roi_outfix);
  if ~exist(fname_aseg_res,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: resampling aseg...\n',mfilename);
    end;
    [vol_aseg,M_aseg] = fs_load_mgh(parms.fname_aseg);
    [vol_aseg_res,M_res] = mmil_resample_vol(vol_aseg,M_aseg,...
      'nvox_ref',size(parms.vol),'M_ref',parms.M,...
      'interpm',0,'M_reg',inv(parms.M_reg));
    fs_save_mgh(vol_aseg_res,fname_aseg_res,M_res);
  end;
  parms.fname_aseg = fname_aseg_res;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create individual ROI mgz from wmparc and aseg
function parms = create_rois(parms)
  % create csv file containing ROI codes and names
  if ~isempty(parms.fname_roi_info) && ~exist(parms.fname_roi_info,'file')
    fname_roi_info = parms.fname_roi_info;
  else
    fname_roi_info = sprintf('%s/%s%s_roi_info.csv',...
      parms.outdir,parms.wmparcname,parms.roi_outfix);
  end;
  if ~isempty(parms.fname_roi_info) && exist(parms.fname_roi_info,'file')
    roi_info = mmil_csv2struct(parms.fname_roi_info);
    valid_roicodes = [roi_info.code];
  else
    valid_roicodes = [];
  end;  
  if ~exist(fname_roi_info,'file') ||...
     ~exist(parms.roi_outdir) || parms.forceflag
    % save each wmparc ROI
    if parms.verbose
      fprintf('%s: creating wmparc ROIs...\n',mfilename);
    end;
    parms = create_wmparc_rois(parms,valid_roicodes);
    % save each aseg ROI
    if parms.verbose
      fprintf('%s: creating aseg ROIs...\n',mfilename);
    end;
    parms = create_aseg_rois(parms,valid_roicodes);
    % combine wmparc_roicodes and aseg_roicodes
    roicodes = cat(1,parms.wmparc_roicodes,parms.aseg_roicodes);
    % combine wmparc_roinames and aseg_roinames
    roinames = cat(1,parms.wmparc_roinames,parms.aseg_roinames);
    % quit if no roicodes were selected
    if isempty(roicodes)
      error('roicodes is empty');
    end;
    % write csv file containing roicodes and roinames
    info = cat(2,num2cell(roicodes),roinames);
    info = cat(1,{'code','name'},info);
    mmil_write_csv(fname_roi_info,info);
  end;
  if isempty(parms.fname_roi_info)
    parms.fname_roi_info = fname_roi_info;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create individual ROI mgz from wmparc and aseg
function parms = create_wmparc_rois(parms,valid_roicodes)
  % load color lookup table
  [roicodes, roinames, rgbv] = fs_colorlut(parms.fname_clut);
  % exclude aseg roi codes
  [roicodes,ind] = setdiff(roicodes,parms.aseg_codes');
  roinames = roinames(ind);
  nroi = length(roicodes);
  % save each ROI as sparse mat file
  mmil_mkdir(parms.roi_outdir);
  vol_wmparc = [];
  i_keep = [];
  for r=1:nroi
    roicode = roicodes(r);
    roiname = roinames{r};
    if roicode==0, continue; end;
    if ~isempty(valid_roicodes) && ~ismember(roicode,valid_roicodes)
      continue;
    end;
    fname_roi = sprintf('%s/%s%s',...
      parms.roi_outdir,roiname,parms.roi_outfix);
    if parms.roi_dilate_flag
      fname_roi = [fname_roi '_dilate'];
    end;
    fname_roi = [fname_roi '.' parms.roi_outtype];
    if ~exist(fname_roi,'file') || parms.forceflag
      % load wmparc file
      if isempty(vol_wmparc)
        if parms.verbose
          fprintf('%s: loading %s...\n',mfilename,parms.fname_wmparc);
        end;
        [vol_wmparc,M_wmparc] = fs_load_mgh(parms.fname_wmparc);
        codes = unique(vol_wmparc(vol_wmparc>0));
        ncodes = length(codes);
      end;
      i_roi = find(ismember(vol_wmparc,roicode));
      if isempty(i_roi), continue; end;
      if parms.verbose
        fprintf('%s: saving %s...\n',mfilename,fname_roi);
      end;
      vol_roi = zeros(size(vol_wmparc));
      vol_roi(i_roi) = 1;
      if parms.roi_dilate_flag
        vol_roi_ctx = ctx_mgh2ctx(vol_roi,M_wmparc);
        %% todo: options for degree of dilation
        vol_roi_ctx = mmil_dilate_mask(vol_roi_ctx,'smooth3',0);
        vol_roi = vol_roi_ctx.imgs;
      end;
      if strcmp(parms.roi_outtype,'mat')
        mmil_save_sparse(vol_roi,fname_roi,M_wmparc);
      else
        fs_save_mgh(vol_roi,fname_roi,M_wmparc);
      end;
    end;
    if exist(fname_roi,'file'), i_keep = [i_keep,r]; end;
  end;
  % save selected roicodes, roinames
  parms.wmparc_roicodes = roicodes(i_keep);
  parms.wmparc_roinames = roinames(i_keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create individual ROI mgz from aseg
function parms = create_aseg_rois(parms,valid_roicodes)
  % load color lookup table
  [roicodes, roinames, rgbv] = fs_colorlut(parms.fname_aseg_clut);
  % exclude all but aseg roi codes
  [roicodes,ind] = intersect(roicodes,parms.aseg_codes');
  roinames = roinames(ind);
  nroi = length(roicodes);
  % save each ROI as sparse mat file
  mmil_mkdir(parms.roi_outdir);
  vol_aseg = [];
  i_keep = [];
  for r=1:nroi
    roicode = roicodes(r);
    roiname = roinames{r};
    if roicode==0, continue; end;
    if ~isempty(valid_roicodes) && ~ismember(roicode,valid_roicodes)
      continue;
    end;
    fname_roi = sprintf('%s/%s%s',...
      parms.roi_outdir,roiname,parms.roi_outfix);
    if parms.roi_dilate_flag
      fname_roi = [fname_roi '_dilate'];
    end;
    fname_roi = [fname_roi '.' parms.roi_outtype];
    if ~exist(fname_roi,'file') || parms.forceflag
      % load aseg file
      if isempty(vol_aseg)
        if parms.verbose
          fprintf('%s: loading %s...\n',mfilename,parms.fname_aseg);
        end;
        [vol_aseg,M_aseg] = fs_load_mgh(parms.fname_aseg);
        codes = unique(vol_aseg(vol_aseg>0));
        ncodes = length(codes);
      end;
      i_roi = find(ismember(vol_aseg,roicode));
      if isempty(i_roi), continue; end;
      if parms.verbose
        fprintf('%s: saving %s...\n',mfilename,fname_roi);
      end;
      vol_roi = zeros(size(vol_aseg));
      vol_roi(i_roi) = 1;
      if parms.roi_dilate_flag
        vol_roi_ctx = ctx_mgh2ctx(vol_roi,M_aseg);
        %% todo: options for degree of dilation
        vol_roi_ctx = mmil_dilate_mask(vol_roi_ctx,'smooth3',0);
        vol_roi = vol_roi_ctx.imgs;
      end;
      if strcmp(parms.roi_outtype,'mat')
        mmil_save_sparse(vol_roi,fname_roi,M_aseg);
      else
        fs_save_mgh(vol_roi,fname_roi,M_aseg);
      end;
    end;
    if exist(fname_roi,'file'), i_keep = [i_keep,r]; end;
  end;
  % save selected roicodes, roinames
  parms.aseg_roicodes = roicodes(i_keep);
  parms.aseg_roinames = roinames(i_keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate spreadsheets with ROI names for each fiber
function parms = create_fiber_rois(parms)
  % create csv file containing ROI codes and names
  if isempty(parms.fname_fiber_rois)
    parms.fname_fiber_rois = sprintf('%s/%s%s_fiber_rois.csv',...
      parms.outdir,parms.wmparcname,parms.roi_outfix);
  end;
  if ~exist(parms.fname_fiber_rois,'file') || parms.forceflag
    % get info about ROIs available
    roi_info = mmil_csv2struct(parms.fname_roi_info);
    roinames = {roi_info.name};
    roicodes = [roi_info.code];

    % identify sub-cortical ROIs
    if parms.sc_flag
      [tmp,i_sc] = intersect(roicodes,parms.sc_codes);
    else
      i_sc = [];
    end;
    n_sc = length(i_sc);
    i_left = find(~cellfun(@isempty,regexp(roinames(i_sc),'Left')));
    i_right = find(~cellfun(@isempty,regexp(roinames(i_sc),'Right')));

    % identify cortical gray and white matter ROIs
    all_gm_codes = [parms.gm_codes(1)+[1:999],parms.gm_codes(2)+[1:999]];
    all_wm_codes = [parms.wm_codes(1)+[1:999],parms.wm_codes(2)+[1:999]];
    if parms.gm_flag==2
      [tmp,i_gm,i_gm_all] = intersect(roicodes,all_gm_codes);
      [tmp,i_wm,i_wm_all] = intersect(roicodes,all_wm_codes);
      % only use matching gm and wm ROIs
      i_all = intersect(i_gm_all,i_wm_all);
    else
      i_all = [1:length(all_gm_codes)];
    end;
    [tmp,i_gm,i_gm_all] = intersect(roicodes,all_gm_codes(i_all));
    [tmp,i_wm,i_wm_all] = intersect(roicodes,all_wm_codes(i_all));
    if ~parms.gm_flag
      i_gm = i_wm;
    end;
    n_gm = length(i_gm);
    i_lh = find(~cellfun(@isempty,regexp(roinames(i_gm),'lh')));
    i_rh = find(~cellfun(@isempty,regexp(roinames(i_gm),'rh')));
    
    % initialize output cell arrays
    fiber_names = cell(n_gm,1);
    fiber_roinames = cell(n_gm,1);
    fiber_not_roinames = cell(n_gm,1);
    f = 1;

    % create entry for every combination of subcortical and cortical ROIs
    if parms.sc_flag
      for i=1:n_sc
        switch i
          case num2cell(i_left)
            hemi1 = 1;
          case num2cell(i_right)
            hemi1 = 2;
          otherwise
            hemi1 = 0;
        end;
        for j=1:n_gm
          switch j
            case num2cell(i_lh)
              hemi2 = 1;
            case num2cell(i_rh)
              hemi2 = 2;
            otherwise
              hemi2 = 0;
          end;
          % eliminate anything between hemisphers
          if hemi1 && hemi2 && hemi1~=hemi2
            continue;
          end;
          switch parms.gm_flag
            case 0
              i_rois = [i_sc(i),i_wm(j)];
            case 1
              i_rois = [i_sc(i),i_gm(j)];
            case 2
              i_rois = [i_sc(i),i_wm(j),i_gm(j)]; % order matters for cutflag=1
          end;
          fiber_roinames{f} = roinames(i_rois);
          % create fiber name based on combination
          fiber_names{f} = sprintf('fiber%s',sprintf('-%s',fiber_roinames{f}{:}));
          % specify NOT ROI with contra aseg wm
          fiber_not_roinames{f} = [];
          if parms.not_contra_wm_flag
            switch hemi2
              case 1
                h_contra = 2;
              case 2
                h_contra = 1;
              otherwise
                h_contra = 0;
            end;
            if h_contra
              wm_code = parms.aseg_wm_codes(h_contra);
              i_not = find(roicodes == wm_code);
              fiber_not_roinames{f} = roinames(i_not);
            end;
          end;
          f = f + 1;
        end;
      end;
    end;
    
    % create entry for every combination of cortical ROIs
    for i=1:n_gm
      switch i
        case num2cell(i_lh)
          hemi1 = 1;
        case num2cell(i_rh)
          hemi1 = 2;
        otherwise
          hemi1 = 0;
      end;
      for j=i+1:n_gm
        switch j
          case num2cell(i_lh)
            hemi2 = 1;
          case num2cell(i_rh)
            hemi2 = 2;
          otherwise
            hemi2 = 0;
        end;
        % eliminate anything between hemisphers
        %% todo: make this optional?
        if hemi1 && hemi2 && hemi1~=hemi2
          continue;
        end;
        switch parms.gm_flag
          case 0
            i_rois = [i_wm(i),i_wm(j)];
          case 1
            i_rois = [i_gm(i),i_gm(j)];
          case 2
            i_rois = [i_gm(i),i_wm(i),i_wm(j),i_gm(j)]; % order matters for cutflag=1
        end;
        fiber_roinames{f} = roinames(i_rois);
        % create fiber name based on combination
        fiber_names{f} = sprintf('fiber%s',sprintf('-%s',fiber_roinames{f}{:}));
        % specify NOT ROI with contra aseg wm
        fiber_not_roinames{f} = [];
        if parms.not_contra_wm_flag
          switch hemi2
            case 1
              h_contra = 2;
            case 2
              h_contra = 1;
            otherwise
              h_contra = 0;
          end;
          if h_contra
            wm_code = parms.aseg_wm_codes(h_contra);
            i_not = find(roicodes == wm_code);
            fiber_not_roinames{f} = roinames(i_not);
          end;
        end;
        f = f + 1;
      end;
    end;

    % create entry for ever combination of cortical ROI with contrlateral
    if parms.contra_homo_flag
      for i=1:length(i_lh)
        roiname = roinames{i_gm(i_lh(i))};
        n = regexp(roiname,'\w+-?h-(?<name>[\w\-]+)','names');
        name = n.name;
        i_bi = find(~cellfun(@isempty,regexp(roinames(i_gm),name)));
        % skip if there is contralateral connection (no exact homolog)
        if ~any(ismember(i_bi,i_lh)) || ~any(ismember(i_bi,i_rh))
          continue;
        end;
        % pick wm, gm, or both
        switch parms.gm_flag
          case 0
            i_rois = [i_wm(i_bi)];
          case 1
            i_rois = [i_gm(i_bi)];
          case 2
            i_rois = [i_gm(i_bi),i_wm(i_bi)];
        end;
        fiber_roinames{f} = roinames(i_rois);
        % create fiber name based on combination
        fiber_names{f} = sprintf('fiber%s',sprintf('-%s',fiber_roinames{f}{:}));
        % specify NOT ROI
        fiber_not_roinames{f} = [];
        f = f + 1;
      end;
    end;

    % write csv file containing fiber names and ROI names
    info = cat(2,fiber_names,fiber_roinames,fiber_not_roinames);
    info = cat(1,{'fiber','roinames','not_roinames'},info);
    mmil_write_csv(parms.fname_fiber_rois,info);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate fiber tracts
function parms = generate_streamlines(parms)
  % run whole brain CSD tractography
  if parms.verbose
    fprintf('%s: running whole brain CSD tractography...\n',mfilename);
  end;
  parms.csd_outdir = [parms.outdir '/csd_tracto'];
  tparms = parms;
  tparms.outdir = parms.csd_outdir;
  %% todo: fiber outstem?
  args = mmil_parms2args(tparms,parms.tracto_tags);
  dti_CSD_tracto(parms.vol,parms.M,parms.qmat,parms.bvals,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select fiber tracts based on spreadsheet
function parms = select_fibers(parms)
  if parms.verbose
    fprintf('%s: selecting fibers...\n',mfilename);
  end;
  % read file defining ROIs for each fiber
  fiber_rois = mmil_csv2struct(parms.fname_fiber_rois);
  % set fname_fiber containing whole-brain streamlines
  %% todo: fiber outstem?
  fname_fiber = sprintf('%s/fiber_path.grp',parms.csd_outdir);
  % set fiber output dir
  parms.fiber_outdir = [parms.outdir '/fibers'];
  mmil_mkdir(parms.fiber_outdir);
  % loop over fibers
  nfibers = length(fiber_rois);
  tparms = parms;
  tparms.fiber = [];
  for f=1:nfibers
    if parms.verbose
      fprintf('%s: selecting streamlines for %s...\n',...
        mfilename,fiber_rois(f).fiber);
    end;
    % set output file name
    %% todo: fiber outstem?
    tparms.fname_out = sprintf('%s/%s.grp',...
      parms.fiber_outdir,fiber_rois(f).fiber);
    % load fiber file
    if isempty(tparms.fiber)
      if ~exist(tparms.fname_out,'file') || parms.forceflag
        if parms.verbose
          fprintf('%s: loading fiber file %s...\n',mfilename,fname_fiber);
        end;
        tic;
        tparms.fiber = dti_read_DTIStudio_fiber(fname_fiber);
        toc;
      end;
    end;
    % set fname_rois
    roinames = fiber_rois(f).roinames;
    nroi = length(roinames);
    fname_rois = cell(nroi,1);
    for r=1:nroi
      fname_rois{r} = sprintf('%s/%s%s.%s',...
        parms.roi_outdir,roinames{r},parms.roi_outfix,parms.roi_outtype);
    end;
    % set fname_NOT
    roinames = fiber_rois(f).not_roinames;
    nroi = length(roinames);
    if nroi
      tparms.fname_NOT = cell(nroi,1);
      for r=1:nroi
        tparms.fname_NOT{r} = sprintf('%s/%s%s.%s',...
          parms.roi_outdir,roinames{r},parms.roi_outfix,parms.roi_outtype);
      end;
    else
      tparms.fname_NOT = [];
    end;
    % select fibers
    args = mmil_parms2args(tparms,parms.select_tags);
    dti_select_fibers(fname_fiber,fname_rois,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

