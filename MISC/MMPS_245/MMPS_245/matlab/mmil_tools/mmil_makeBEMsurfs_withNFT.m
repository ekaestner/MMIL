function mmil_makeBEMsurfs_withNFT(cpath,fspath,varargin)
%function mmil_makeBEMsurfs_withNFT(cpath,fspath,[options])
%
% Purpose: create brain and skull surface meshses
%   using NFT (Neuroelectromagnetic Forward Head Modeling Toolbox)
%   based on FreeSurfer aseg and existing BEM surfaces
%
% Required Input:
%   cpath: full path of processing MRI container
%   fspath: full path of FreeSurfer recon
%
% Optional Input:
%   'fstem_aseg': file stem of FreeSurfer aseg
%     {default = 'aseg'}
%   'fname_brain': full path name of brain image
%     used in auto-generated scripts for displaying surfaces
%     if empty, will use brain mask from aseg
%     {default = []}
%   'input_bem_dir': input subdirectory in fspath containing BEM surfaces
%     {default = 'bem_PD'}
%   'output_bem_dir': output subdirectory in fspath
%     may be full path instead
%     {default = 'bem_PD_NFT'}
%   'outstem': output file stem
%     {default = 'PD_NFT'}
%   'surf_names': names for input mesh files
%     {default = {'none','inner_skull4','outer_skull4','outer_scalp4'}}
%   'mask_names': names for output mask files
%     {default = {'brainmask','innerskullmask','outerskullmask','scalpmask'}
%   'mesh_names': names for output mesh files
%     {default = {'brain','inner_skull','outer_skull','scalp'}}
%   'smooth1': first round smoothing kernel for creating brain mask from aseg
%     {default = 40
%   'thresh1': first round threshold for creating brain mask from aseg
%     {default = 0.5}
%   'smooth2': second round smoothing kernel for creating brain mask from aseg
%     {default = 12}
%   'thresh2': second round threshold for creating brain mask from aseg
%     {default = 0.2}
%   'smooth3': third round smoothing kernel for creating brain mask from aseg
%     {default = 0}
%   'check_flag': [0|1] check outer masks include inner mask voxels
%     if dilate_flag, also includes dilated brain masks
%     {default = 1}
%   'dilate_flag': [0|1|2] interatively dilate brain mask
%     0: no dilation
%     1: dilation for inner skull, outer skull, and outer scalp
%     2: dilation for brain, inner skull, outer skull, and outer scalp
%     {default = 2}
%   'dilate_smooth': smoothing kernel for brain mask dilation
%     {default = 25}
%   'dilate_thresh': threshold for brain mask dilation
%     {default = 0.2}
%   'nmesh': number of vertices for initial mesh creation
%     {default = 10000}
%   'downsample_flag': [0|1] downsample meshes
%     {default = 1}
%   'dmesh': approximate number of vertices for downsampled meshes
%     {default = 2000}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  06/16/14 by Don Hagler
% Last Mod: 01/14/15 by Don Hagler
%

%% todo: output mesh fnames with simple names - create links?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% check input parameters
parms = check_input(cpath,fspath,varargin);

% create brain mask from aseg
parms = create_brain_mask(parms);

% create skull and scalp masks from meshes
parms = create_masks_from_meshes(parms);

% iteratively dilate brain mask
parms = dilate_brain_mask(parms);

% make sure outer masks include inner masks
parms = check_masks(parms);

% create nft meshes
parms = create_meshes(parms);

% create scripts to view surfaces with tkmedit
parms = write_scripts(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(cpath,fspath,options)
  parms = mmil_args2parms(options,{...
    'cpath',cpath,[],...
    'fspath',fspath,[],...
  ...
    'fstem_aseg','aseg',[],...
    'fname_brain',[],[],...
    'input_bem_dir','bem_PD',[],...
    'output_bem_dir','bem_PD_NFT',[],...
    'outstem','PD_NFT',[],...
    'surf_names',{'none','inner_skull4','outer_skull4','outer_scalp4'},[],...
    'mask_names',{'brainmask','innerskullmask','outerskullmask','scalpmask'},[],...
    'mesh_names',{'brain','inner_skull','outer_skull','outer_scalp'},[],...
    'smooth1',40,[],...
    'thresh1',0.5,[],...
    'smooth2',12,[],...
    'thresh2',0.2,[],...
    'smooth3',0,[],...
    'check_flag',true,[false true],...
    'dilate_flag',2,[0:2],...
    'dilate_smooth',25,[],...
    'dilate_thresh',0.2,[],...
    'nmesh',10000,[],...
    'downsample_flag',true,[false true],...
    'dmesh',2000,[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'aseg_tags',{'smooth1','thresh1','smooth2','thresh2','smooth3',...
                 'forceflag','fname_in','fname_out'},[],...
  });

  % check input directories
  if ~exist(parms.cpath,'dir')
    error('directory %s not found',parms.cpath);
  end;
  if ~exist(parms.fspath,'dir')
    error('directory %s not found',parms.fspath);
  end;
  parms.indir = [parms.cpath '/' parms.input_bem_dir];
  if ~exist(parms.indir,'dir')
    error('input directory %s not found',parms.indir);
  end;
  
  % check aseg
  parms.fname_aseg = sprintf('%s/mri/%s.mgz',parms.fspath,parms.fstem_aseg);
  if ~exist(parms.fname_aseg,'file')
    error('file %s not found',parms.fname_aseg);
  end;

  % check brain image
  if ~isempty(parms.fname_brain) && ~exist(parms.fname_brain,'file')
    error('file %s not found',parms.fname_brain);
  end;

  parms.nmasks = length(parms.mask_names);

  % create output directory
  if mmil_isrelative(parms.output_bem_dir)
    parms.outdir = [parms.cpath '/' parms.output_bem_dir];
  else
    parms.outdir = parms.output_bem_dir;
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_masks_from_meshes(parms)
  parms.vol = [];
  for i=2:parms.nmasks
    surf_name = parms.surf_names{i};
    mask_name = parms.mask_names{i};
    fname_mesh = sprintf('%s/%s.tri',parms.indir,surf_name);
    if ~exist(fname_mesh,'file')
      error('file %s not found',fname_mesh);
    end;
    fname_out = sprintf('%s/%s.mgz',parms.outdir,mask_name);
    if ~exist(fname_out,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: creating mask from %s...\n',mfilename,fname_mesh);
      end;
      if isempty(parms.vol)
        parms.vol = ctx_load_mgh(parms.fname_brain);
      end;
      surf = fs_read_trisurf(fname_mesh);
      surf.vertices(:,1) = -surf.vertices(:,1) + parms.vol.lphcent(1);
      surf.vertices(:,2) = -surf.vertices(:,2) + parms.vol.lphcent(2);
      surf.vertices(:,3) = surf.vertices(:,3) + parms.vol.lphcent(3);
      msurf = preprocessQ(surf);
      [tmp,vol_mask] = getmaskvol(parms.vol,msurf,eye(4));
      ctx_save_mgh(vol_mask,fname_out);
    end;
    parms.fname_masks{i} = fname_out;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_brain_mask(parms)
  parms.fname_masks = cell(1,parms.nmasks);
  parms.fname_brainmask = sprintf('%s/%s.mgz',parms.outdir,parms.mask_names{1});
  if ~exist(parms.fname_brainmask,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: creating brain mask from aseg...\n',mfilename);
    end;
    tparms = parms;
    tparms.fname_in = parms.fname_aseg;
    tparms.fname_out = parms.fname_brainmask;
    args = mmil_parms2args(tparms,parms.aseg_tags);
    mmil_dilate_mask([],args{:});
  end;
  parms.fname_masks{1} = parms.fname_brainmask;
  if isempty(parms.fname_brain)
    parms.fname_brain = parms.fname_brainmask;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = dilate_brain_mask(parms)
  for i=1:parms.nmasks
    parms.fname_brainmasks{i} = parms.fname_brainmask;
  end;
  if parms.dilate_flag
    if parms.dilate_flag==2
      nmasks = 4;
    else
      nmasks = 3;
    end;
    tparms = [];
    tparms.smooth1 = parms.dilate_smooth;
    tparms.thresh1 = parms.dilate_thresh;
    tparms.smooth2 = 0;
    tparms.smooth3 = 0;
    tparms.forceflag = parms.forceflag;
    tparms.fname_in = parms.fname_brainmask;
    j = 1;
    for i=1:parms.nmasks
      if i>1 || parms.dilate_flag==2
        % replace brain mask with dilated mask
        parms.fname_brainmasks{i} = sprintf('%s/%s_diln%d.mgz',...
          parms.outdir,parms.mask_names{1},j);
        tparms.fname_out = parms.fname_brainmasks{i};
        if ~exist(tparms.fname_out,'file') || parms.forceflag
          if parms.verbose
            fprintf('%s: dilating brain mask (%d)...\n',mfilename,j);
          end;
          args = mmil_parms2args(tparms);
          mmil_dilate_mask([],args{:});
        end;
        tparms.fname_in = tparms.fname_out;
        j = j + 1;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_masks(parms)
  if parms.check_flag
    for i=1:parms.nmasks
      mask_name = parms.mask_names{i};
      fname_mask = parms.fname_masks{i};
      fname_brainmask = parms.fname_brainmasks{i};
      if i>1
        fname_inner = parms.fname_masks{i-1};
      else
        fname_inner = parms.fname_masks{i};
      end;
      fname_masks = unique({fname_mask,fname_brainmask,fname_inner});
      if length(fname_masks)==1
        fname_out = fname_mask;
      else
        if parms.dilate_flag
          fname_out = sprintf('%s/%s_check_dil%d.mgz',...
            parms.outdir,mask_name,parms.dilate_flag);
        else
          fname_out = sprintf('%s/%s_check.mgz',...
            parms.outdir,mask_name);
        end;
        if ~exist(fname_out,'file') || parms.forceflag
          if parms.verbose
            fprintf('%s: checking %s...\n',mfilename,mask_name);
          end;
          % compare masks
          vol1 = [];
          for j=1:length(fname_masks)
            if isempty(vol1)
              [vol1,M] = fs_load_mgh(fname_masks{j});
              continue;
            end;
            [vol2,M] = fs_load_mgh(fname_masks{j});
            vol1 = 1.0*(vol1 | vol2);
          end;
          fs_save_mgh(vol1,fname_out,M);
        end;
      end;
      parms.fname_masks{i} = fname_out;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_meshes(parms)
  % reverse order of masks for NFT
  mask_names = fliplr(parms.mask_names);
  mesh_names = fliplr(parms.mesh_names);
  fname_masks = fliplr(parms.fname_masks);

  % create segmentation file in format expected by NFT
  fname_seg = sprintf('%s/%s_segments.mat',...
    parms.outdir,parms.outstem);
  if ~exist(fname_seg,'file') || parms.forceflag
    % initialize seg struct
    Segm = [];
    Segm.parameters.LRflip = 0;
    % load masks and add to struct
    for i=1:parms.nmasks
      mask_name = mask_names{i};
      fname_mask = fname_masks{i};
      vol_mask = mmil_reorient_for_nft(fname_mask,parms.outdir,parms.forceflag);
      Segm.(mask_name) = (vol_mask>0);
    end;
    save(fname_seg,'Segm');
  end;

  % create meshes
  mesh_suffix = sprintf('_n%d',parms.nmesh);
  mesh_outstem = sprintf('%s_mesh%s',parms.outstem,mesh_suffix);
  mmil_create_nft_mesh(fname_seg,...
    'outdir',parms.outdir,...
    'outstem',mesh_outstem,...
    'nmesh',parms.nmesh,...
    'verbose',parms.verbose,...
    'forceflag',parms.forceflag);

  % convert to tri files
  fname_mesh = [parms.outdir '/' mesh_outstem '.mat'];
  mmil_convert_nft_mesh(fname_mesh,...
    'outdir',parms.outdir,...
    'mesh_names',mesh_names,...
    'mesh_suffix',mesh_suffix,...
    'forceflag',parms.forceflag);
  parms.fname_meshes = cell(1,parms.nmasks);
  for i=1:parms.nmasks
    parms.fname_meshes{i} = sprintf('%s/%s%s.tri',...
      parms.outdir,mesh_names{i},mesh_suffix);
  end;

  % downsample tri files
  if parms.downsample_flag
    dmesh_suffix = sprintf('%s_d%d',mesh_suffix,parms.dmesh);
    for i=1:parms.nmasks
      fname_in = sprintf('%s/%s%s.tri',...
        parms.outdir,mesh_names{i},mesh_suffix);
      fname_out = sprintf('%s/%s%s.tri',...
        parms.outdir,mesh_names{i},dmesh_suffix);
      mmil_downsample_mesh(fname_in,fname_out,...
        'dmesh',parms.dmesh,...
        'forceflag',parms.forceflag);
      parms.fname_meshes{i} = fname_out;
    end;
  end;

  % reverse order of meshes
  parms.fname_meshes = fliplr(parms.fname_meshes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = write_scripts(parms)
  for i=1:parms.nmasks
    mesh_name = parms.mesh_names{i};
    fname_mesh = parms.fname_meshes{i};
    fname_csh = sprintf('%s/view_%s.csh',parms.outdir,mesh_name);
    if ~exist(fname_csh,'file') || parms.forceflag
      fid = fopen(fname_csh,'wt');
      if fid<0
        error('unable to create %s',fname_csh);
      end;
      fprintf(fid,'#!/bin/csh -f\n');
      fprintf(fid,'# view_%s.csh  %s\n',date,mesh_name);
      fprintf(fid,'tkmedit -f "%s" -surface "%s"\n',...
        parms.fname_brain,fname_mesh);
      fclose(fid);
      cmd = sprintf('chmod +x %s',fname_csh);
      [status,result] = unix(cmd);
      if status
        fprintf('%s: WARNING: chmod exited with errors:\n%s\n',mfilename,result);
      end;
    end;
  end;
return;


