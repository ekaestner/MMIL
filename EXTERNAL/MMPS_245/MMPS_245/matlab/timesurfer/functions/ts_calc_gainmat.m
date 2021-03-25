function G_xyz = ts_calc_gainmat(avg_data,varargin)
%function G_xyz = ts_calc_gainmat(avg_data,[options])
%
% Purpose: calculate gain matrix (forward solution)
%   using boundary element method or spherical shell model
%
% Required input:
%  avg_data: structure containing averaged MEG or EEG data
%    e.g. output of ts_avg_fif_data
%
% Optional parameters for specifying dipoles (two alternatives):
%  'fname_lh_surf': name of FreeSurfer binary format surface file
%     {default = []}
%  'fname_rh_surf': name of FreeSurfer binary format surface file
%     {default = []}
%  'lh_dip_info': matrix with size [6,nverts]
%    containing locations and orientations for each dipole in left hemisphere
%    generated with ts_dip_info or ts_read_dip_file
%    NOTE: required if fname_lh_surf not specified 
%    {default = []}
%  'rh_dip_info':  matrix with size [6,nverts]
%    containing locations and orientations for each dipole in left hemisphere
%    NOTE: required if fname_rh_surf not specified 
%    {default = []}
%
% Optional input:
%  'trans': 4x4 transformation matrix specifying registration between
%    mri (freesurfer brain space) and head (MEG/EEG sensor space)
%    can be generated with loadtrans function (fiff access) or with ts_pointreg
%    if not specified, will use avg_data.coor_trans.mri2head
%    {default = []}
%  'bem_surf_files': cell array of file names for MRI-derived boundary surfaces
%    in FreeSurfer tri file text format; order should be inside-out
%    NOTE: required if bem_flag = 1
%    if four or more layers, openmeeg_flag must = 1
%    {default = []}
%  'bem_matfile': name of matfile to be created for storage of
%    BEM transfer functions; if empty, will not save
%    {default = []}
%  'lh_dec_dips': vector containing zeros and ones specifying which left
%    hemisphere dipoles to include; can be generated with ts_read_dec_file
%    if not specified, will use all dipoles
%    {default = []}
%  'rh_dec_dips': vector containing zeros and ones specifying which right
%    hemisphere dipoles to include
%    if not specified, will use all dipoles
%    {default = []}
%  'bem_flag': [0|1] use boundary element method
%    otherwise use spherical shell model
%    {default = 1}
%  'openmeeg_flag': [0|1] use OpenMEEG to calculate BEM solution
%    otherwise use M.X. Huang's BEM solution function
%    {default = 0}
%  'useMEG_flag': [0|1] calculate MEG gain matrix if there are MEG sensors
%    {default = 1}
%  'useEEG_flag': [0|1] calculate EEG gain matrix if there are EEG sensors
%    NOTE: either useMEG_flag or useEEG_flag must = 1 (both cannot be 0)
%    {default = 1}
%  'conductivities': vector of conductivity values of each tissue layer:
%    e.g. brain, skull, and scalp
%    {default = [0.3 0.01 0.3]}
%  'cen_sph': coordinates of center of sphere for spherical head model
%    supply as [x y z] (in mm)
%    {default = [0 0 0]}
%  'radii': vector of three values specifying radius (in mm) of each shell:
%    NOTE: required if bem_flag = 0 and useEEG_flag = 1
%    {default: []}
%  'rootoutdir': root output directory (temporary subdirectory may be created)
%    {default = pwd} (current working directory)
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output (bem_matfile)
%    {default = 0}
%
% Output:
%  G_xyz: gain matrix forward solution with three rows (x,y,z) for each dipole
%         units for gradiometers are fT/cm per nA source
%         units for magnetometers are fT per nA source
%         units for EEG sensors are uV per nA source
%
%         size of G_xyz is 3*ndips x nsensors
%         if useEEG_flag=0, EEG columns of G_xyz will contain zeros
%         if useMEG_flag=0, MEG columns of G_xyz will contain zeros
%         for non-MEG/EEG columns, G_xyz will contain zeros
%
% Created:  11/17/05 by Don Hagler
% Last Mod: 10/28/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_xyz = [];
if ~mmil_check_nargs(nargin,2), return; end;

% check input parameters
parms = check_input(avg_data,varargin);

% extract necessary info from sensor info
parms = prep_sensor_info(parms,avg_data);

% transform dipole coordinates from to "head" space
[source_locs,parms] = transform_dipoles(parms);

% calculate gain matrices for MEG and/or EEG
if parms.bem_flag
  % load boundary surfaces
  [bem_verts,bem_faces] = load_bem_surfs(parms);
  if parms.openmeeg_flag
    % calculate BEM gain matrix with fieldtrip/openmeeg code
    G_xyz = calc_bem_gainmat_openmeeg(parms,source_locs,bem_verts,bem_faces);
  else
    % calculate BEM gain matrix with mxhuang code
    G_xyz = calc_bem_gainmat_mxhuang(parms,source_locs,bem_verts,bem_faces);
  end;
else
  % calculate spherical model gain matrix
  G_xyz = calc_sph_gainmat(parms,source_locs);
end;

cleanup_tempfiles(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(avg_data,options)
  parms_filter = {...
    'avg_data',avg_data,[],...
  ... % dipoles
    'fname_lh_surf',[],[],...
    'fname_rh_surf',[],[],...
    'lh_dip_info',[],[],...
    'rh_dip_info',[],[],...
  ... % optional
    'trans',[],[],...
    'bem_surf_files',[],[],...
    'bem_matfile',[],[],...
    'lh_dec_dips',[],[],...
    'rh_dec_dips',[],[],...
    'bem_flag',true,[false true],...
    'openmeeg_flag',false,[false true],...
    'useMEG_flag',true,[false true],...
    'useEEG_flag',true,[false true],...
    'conductivities',[0.3 0.01 0.3],[],...
    'cen_sph',[0 0 0],[],...
    'radii',[],[],...
    'rootoutdir',pwd,[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'batchsize',3000,[1,1e10],...
  };
  parms = mmil_args2parms(options,parms_filter);

  % check trans
  if isempty(parms.trans)
    parms.trans = avg_data.coor_trans.mri2head;
    if isempty(parms.trans)
      error('trans not specified and not found in avg_data.coor_trans');
    end;
  end;

  % check that either useMEG_flag or useEEG_flag is 1
  if ~parms.useMEG_flag & ~parms.useEEG_flag
    error('both useMEG_flag and useEEG_flag cannot be zero');
  end;

  % check that dipoles are specified
  if isempty(parms.fname_lh_surf) && isempty(parms.lh_dip_info)
    error('either fname_lh_surf or lh_dip_info must be specified');
  end;
  if isempty(parms.fname_rh_surf) && isempty(parms.rh_dip_info)
    error('either fname_rh_surf or rh_dip_info must be specified');
  end;

  % check dipole surf files exist
  if ~isempty(parms.fname_lh_surf) && ~exist(parms.fname_lh_surf,'file')
    error('file %s not found',parms.fname_lh_surf);
  end;
  if ~isempty(parms.fname_rh_surf) && ~exist(parms.fname_rh_surf,'file')
    error('file %s not found',parms.fname_rh_surf);
  end;

  % check inputs for BEM or sphere
  if parms.bem_flag
    % check that bem_surf_files is not empty
    if isempty(parms.bem_surf_files)
      error('bem_surf_files required because bem_flag = 1');
    end;
    % set number of layers
    if ~iscell(parms.bem_surf_files)
      parms.bem_surf_files = {parms.bem_surf_files};
    end;
    parms.nlayers = length(parms.bem_surf_files);
    % check that nlayers is 1 or 3 if openmeeg ~= 1
    if ~parms.openmeeg_flag && ~ismember(parms.nlayers,[1,3])
      error('openmeeg_flag=1 required if number of layers not 1 or 3');    
    end;
    % check bem_surf_files exist
    for i=1:parms.nlayers
      fname = parms.bem_surf_files{i};
      if isempty(fname), error('bem_surf_files(%d) is empty',i); end;
      if ~exist(fname,'file')
        error('file %s not found',fname);
      end;
    end;
    % check if bem_matfile specified
    if isempty(parms.bem_matfile) || parms.openmeeg_flag
      parms = create_tempdir(parms);
    end;
    if isempty(parms.bem_matfile)
      parms.bem_matfile = [parms.tempdir '/bem_transfer.mat'];
      parms.cleanup_bem_matfile_flag = 1;
    else
      parms.cleanup_bem_matfile_flag = 0;    
    end;
    % set scaling factors
    if parms.openmeeg_flag
      parms.MEG_sf = 1e6;
      parms.EEG_sf = 1e-3;
    else
      parms.MEG_sf = 1;
      parms.EEG_sf = 1;
    end; 
  else
    parms.nlayers = length(parms.radii);
    % set scaling factors
    parms.MEG_sf = 1;
    parms.EEG_sf = 1e-3;
  end;

  % check that conductivities vector has enough values
  if length(parms.conductivities)<parms.nlayers
    msg = sprintf('conductivities must have at least %d value',parms.nlayers);
    if parms.nlayers>1, msg = [msg 's']; end;
    error(msg);
  elseif length(parms.conductivities)>parms.nlayers
    parms.conductivities = parms.conductivities(1:parms.nlayers);
  end;

  % complain if number of layers is too low for EEG
  if parms.useEEG_flag & parms.nlayers < 3 && parms.verbose
    fprintf('%s: WARNING: for EEG, should use at least 3 layers (using %d)\n',...
      mfilename,parms.nlayers);
  end;

  % add new fieldtrip functions to path for OpenMEEG
  if parms.openmeeg_flag
    mmil_addpath_new_fieldtrip;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_tempdir(parms)
  % create a temporary directory with a randomly generated number
  rand('twister',sum(100*clock));
  dirname = sprintf('%s/temp_%d',parms.rootoutdir,round(rand(1)*1e7));
  mmil_mkdir(dirname)
  % create a function to override the built-in tempdir function
  % NOTE: this is so that multiple instances of openmeeg will not fill up /tmp
  fname_script = sprintf('%s/tempdir.m',dirname);
  fid = fopen(fname_script,'wt');
  fprintf(fid,'function dirname = tempdir()\n');
  fprintf(fid,'  dirname = ''%s'';\n',dirname);
  fprintf(fid,'return\n');
  fclose(fid);
  % add dirname to path so that tempdir function will be used
  addpath(dirname);
  parms.tempdir = dirname;
  parms.tempdir_script = fname_script;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tempfiles(parms)
  if isfield(parms,'tempdir') && exist(parms.tempdir,'dir')
    % if bem_matfile was temporary, delete now
    if parms.cleanup_bem_matfile_flag && exist(parms.bem_matfile,'file')
      delete(parms.bem_matfile);
    end;
    if exist(parms.tempdir_script,'file')
      delete(parms.tempdir_script);
    end;
    rmpath(parms.tempdir);
    [s,m] = rmdir(parms.tempdir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_sensor_info(parms,avg_data)
  if parms.verbose
    fprintf('%s: preparing sensor info...\n',mfilename);
  end;
  parms.nsensors = length(avg_data.sensor_info);
  chan_types = lower({avg_data.sensor_info.typestring});
  % get integration points from sensor info
  if parms.useMEG_flag
    parms.MEG_info = ts_prep_MEG_info(avg_data);
    grad_chans = find(strncmp('grad',chan_types,length('grad')));
    mag_chans  = find(strcmp('mag',chan_types));
    parms.MEG_chans = union(grad_chans,mag_chans);
  else
    parms.MEG_chans = [];
  end;
  if parms.useEEG_flag
    parms.EEG_info = ts_prep_EEG_info(avg_data);
    parms.EEG_chans = find(strcmp('eeg',chan_types));
  else
    parms.EEG_chans = [];
  end;
  if parms.openmeeg_flag
    % prepare elec struct required by fieldtrip
    if parms.useEEG_flag
      elec = [];
      elec.pnt = parms.EEG_info.intpnt;
      parms.npnts_eeg = size(elec.pnt,1);
      elec.label = cell(1,parms.npnts_eeg);
      for i=1:parms.npnts_eeg
        elec.label{i} = sprintf('vertex%03d',i);
      end
      parms.elec = elec;
    end;
    % prepare grad struct required by fieldtrip
    if parms.useMEG_flag
      grad = [];
      grad.pnt = parms.MEG_info.intpnt_loc;
      grad.ori = parms.MEG_info.intpnt_ori;
      parms.npnts_meg = size(grad.pnt,1);
      grad.label = cell(1,parms.npnts_meg);
      for i=1:parms.npnts_meg
        grad.label{i} = sprintf('vertex%03d',i);
      end
      parms.grad = grad;
    end;
  end;
return;

s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [source_locs,parms] = transform_dipoles(parms)
  source_locs = [];

  if parms.verbose
    fprintf('%s: transforming dipole locations...\n',mfilename);
  end;

  % load dipoles from surface files or use locations from dip_info
  if ~isempty(parms.fname_lh_surf)
    surf_lh = fs_read_surf(parms.fname_lh_surf);
    lh_dip_locs = surf_lh.vertices';
  else
    lh_dip_locs = parms.lh_dip_info(1:3,:);
  end;
  if ~isempty(parms.fname_rh_surf)
    surf_rh = fs_read_surf(parms.fname_rh_surf);
    rh_dip_locs = surf_rh.vertices';
  else
    rh_dip_locs = parms.rh_dip_info(1:3,:);
  end;

  % use dec dips if supplied
  if ~isempty(parms.lh_dec_dips)
    ind = find(parms.lh_dec_dips==1);
    lh_dip_locs = lh_dip_locs(:,ind);
  end;
  if ~isempty(parms.rh_dec_dips)
    ind = find(parms.rh_dec_dips==1);
    rh_dip_locs = rh_dip_locs(:,ind);
  end;

  % combine info for both hemispheres
  grid_mri = [lh_dip_locs';rh_dip_locs'];
  parms.nsources=size(grid_mri,1);
  if parms.verbose
    fprintf('%s: total number of dipoles in forward = %d\n',mfilename,parms.nsources);
  end;

  % apply mri2head coordinate transformation to get to "head"-space
  grid_head = parms.trans*[grid_mri'/1000;ones(1,parms.nsources)]; % mm to meters
  grid_head = grid_head(1:3,:)'; % discard extra row of ones and permute back

  % set source locs
  source_locs = grid_head;
  if ~parms.bem_flag
    source_locs = bsxfun(@minus,source_locs,parms.cen_sph/1000);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [verts,faces] = load_bem_surfs(parms)
  [verts,faces]=ts_load_bem_surfs(parms.bem_surf_files);
  for i=1:parms.nlayers
    tmp_verts=verts{i};
    % rescale units as necessary
    maxdist=max(tmp_verts(:,1))-min(tmp_verts(:,1));
    if maxdist > 100 % original unit in mm
      tmp_verts=tmp_verts(:,1:3)/1000; % mm to m
    elseif maxdist > 10 & maxdist < 30 % original unit in cm
      tmp_verts=tmp_verts(:,1:3)/100; % cm to m
    else
      tmp_verts=tmp_verts(:,1:3); % in m
    end
    nvert=size(tmp_verts,1);
    % apply coordinate transformation to get to "head"-space
    tmp_verts=(parms.trans*[tmp_verts';ones(1,nvert)])';
    verts{i}=tmp_verts(:,1:3);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G_xyz = calc_bem_gainmat_openmeeg(parms,source_locs,bem_verts,bem_faces)
  % create vol struct from bem_verts and bem_faces
  vol = [];
  vol.bnd = [];
  vol.cond = parms.conductivities;
  for i=1:parms.nlayers
    vol.bnd(i).pnt = bem_verts{i};
    vol.bnd(i).tri = bem_faces{i};
  end;

  % calculate BEM transfer matrix
  if ~exist(parms.bem_matfile,'file') || parms.forceflag
    % set options
    cfg = [];
    cfg.method = 'openmeeg';
    % call ft_prepare_bemmodel
    if parms.verbose
      fprintf('%s: calculating BEM transfer matrix...\n',mfilename);
      tic;
    end;
    vol_trans = ft_prepare_bemmodel(cfg,vol);
    if parms.verbose, toc; end;
    if ~isfield(vol_trans,'mat')
      error('ft_prepare_bemmodel failed');
    else
      save(parms.bem_matfile,'vol_trans','-v7.3');
    end;
  else
    if parms.verbose
      fprintf('%s: loading pre-calculated BEM transfer matrix...\n',mfilename);
      tic;
    end;
    load(parms.bem_matfile);
    if parms.verbose, toc; end;
  end;

  % calculate gain matrix for MEG
  if parms.useMEG_flag
    cfg = [];
    cfg.grad = parms.grad;
    cfg.vol = vol_trans;
    cfg.method = 'openmeeg';
    cfg.reducerank = 'no';
    G_xyz_inp.meg = zeros(parms.npnts_meg,3,parms.nsources);
    % loop over batches of sources (seg fault if too many)
    nbatch = ceil(parms.nsources/parms.batchsize);
    if parms.verbose
      fprintf('%s: calculating MEG lead fields...\n',mfilename);
    end;
    for b=1:nbatch
      j = 1 + (b-1)*parms.batchsize;
      k = min(j + parms.batchsize - 1,parms.nsources);
      cfg.grid.pos = source_locs(j:k,:);
      if nbatch>1 && parms.verbose
        fprintf('%s: batch %d of %d (dipoles %d to %d)...\n',...
          mfilename,b,nbatch,j,k);
      end;
      if parms.verbose, tic; end;
      tmp = ft_prepare_leadfield(cfg);
      if parms.verbose, toc; end;
      % create G_xyz from leadfield cell array
      for i=1:length(tmp.leadfield)
        G_xyz_inp.meg(:,:,j+i-1) = tmp.leadfield{i};
      end;
    end;
    % reshape 3D matrix to 2D
    G_xyz_inp.meg = reshape(G_xyz_inp.meg,[parms.npnts_meg,3*parms.nsources]);
  end;
  
  % calculate gain matrix for EEG
  if parms.useEEG_flag
    cfg = [];
    cfg.elec = parms.elec;
    cfg.vol = vol_trans;
    cfg.method = 'openmeeg';
    cfg.reducerank = 'no';
    G_xyz_inp.eeg = zeros(parms.npnts_eeg,3,parms.nsources);
    % loop over batches of sources (seg fault if too many)
    nbatch = ceil(parms.nsources/parms.batchsize);
    if parms.verbose
      fprintf('%s: calculating EEG lead fields...\n',mfilename);
    end;
    for b=1:nbatch
      j = 1 + (b-1)*parms.batchsize;
      k = min(j + parms.batchsize - 1,parms.nsources);
      cfg.grid.pos = source_locs(j:k,:);
      if nbatch>1 && parms.verbose
        fprintf('%s: batch %d of %d (dipoles %d to %d)...\n',...
          mfilename,b,nbatch,j,k);
      end;
      if parms.verbose, tic; end;
      tmp = ft_prepare_leadfield(cfg);
      if parms.verbose, toc; end;
      % create G_xyz from leadfield cell array
      for i=1:length(tmp.leadfield)
        G_xyz_inp.eeg(:,:,j+i-1) = tmp.leadfield{i};
      end;
    end;
    % reshape 3D matrix to 2D
    G_xyz_inp.eeg = reshape(G_xyz_inp.eeg,[parms.npnts_eeg,3*parms.nsources]);
  end;

  % convert integration points into channels
  if parms.verbose
    fprintf('%s: converting integration points to channels...\n',mfilename);
  end;
  if parms.useMEG_flag
    G_xyz_MEG = ts_gain_intpnt2chan(parms.MEG_info.Coil,G_xyz_inp.meg);
  else
    G_xyz_MEG = [];
  end;
  if parms.useEEG_flag
    G_xyz_EEG = ts_gain_intpnt2chan(parms.EEG_info.sensor,G_xyz_inp.eeg);
  else
    G_xyz_EEG = [];
  end;

  % combine MEG and EEG gain matrices
  G_xyz = combine_gain_matrices(parms,G_xyz_MEG,G_xyz_EEG);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G_xyz = calc_bem_gainmat_mxhuang(parms,source_locs,bem_verts,bem_faces)
  % set bem_input
  if parms.useMEG_flag
    bem_input.R_meg = parms.MEG_info.intpnt_loc; % location of MEG integration points  
    bem_input.O_meg = parms.MEG_info.intpnt_ori; % orientation of MEG integration points
  else
    bem_input.R_meg = [];
    bem_input.O_meg = [];
  end;
  if parms.useEEG_flag
    bem_input.R_eeg = parms.EEG_info.intpnt; % location of EEG integration points
  else
    bem_input.R_eeg = [];
  end;
  if parms.useMEG_flag & parms.useEEG_flag
    bem_input.mode = 3; % MEG and EEG
    if parms.verbose
      fprintf('%s: will calculate BEM gain matrix for MEG and EEG\n',mfilename);
    end;
  elseif parms.useMEG_flag
    bem_input.mode = 2; % MEG only
    if parms.verbose
      fprintf('%s: will calculate BEM gain matrix for MEG only\n',mfilename);
    end;
  elseif parms.useEEG_flag
    bem_input.mode = 1; % EEG only
    if parms.verbose
      fprintf('%s: will calculate BEM gain matrix for EEG only\n',mfilename);
    end;
  end;
  bem_input.vertices = bem_verts; % location for BEM vertices
  bem_input.faces = bem_faces; % face connection matrix for BEM mesh
  bem_input.sigma = parms.conductivities; % conductivity SI unit
  bem_input.basis_opt = 1; % linear potential function
  bem_input.test_opt = 0; % collocation
  bem_input.ISA = 0; % Inhibit Isolated Skull Approach
  bem_input.fn_eeg = 'bem_eeg_transfer';
  bem_input.fn_meg = 'bem_meg_transfer';

  % calculate BEM transfer matrix
  if parms.verbose
    fprintf('%s: calculating BEM transfer matrix...\n',mfilename);
    tic
  end;
  if exist(parms.bem_matfile,'file') && parms.forceflag
    delete(parms.bem_matfile);
  end;
  if parms.nlayers==1
    bem_transfer = bem_transfer_1shell(bem_input,[],parms.bem_matfile,1);
  else
    bem_transfer = bem_transfer_3shell(bem_input,[],parms.bem_matfile,1);
  end;
  if parms.verbose, toc; end;

  % calculate BEM gain matrix (leadfield) for integration points
  if parms.verbose
    fprintf('%s: calculating BEM gain matrix...\n',mfilename);
    tic
  end;
  G_xyz_inp=bem_gainmat(source_locs,bem_transfer,0);
  if parms.verbose, toc; end;

  % convert integration points into channels
  if parms.verbose
    fprintf('%s: converting integration points to channels...\n',mfilename);
  end;
  if parms.useMEG_flag
    G_xyz_MEG = ts_gain_intpnt2chan(parms.MEG_info.Coil,G_xyz_inp.meg);
  else
    G_xyz_MEG = [];
  end;
  if parms.useEEG_flag
    G_xyz_EEG = ts_gain_intpnt2chan(parms.EEG_info.sensor,G_xyz_inp.eeg);
  else
    G_xyz_EEG = [];
  end;

  % combine MEG and EEG gain matrices
  G_xyz = combine_gain_matrices(parms,G_xyz_MEG,G_xyz_EEG);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G_xyz = calc_sph_gainmat(parms,source_locs)
  % calculate gain matrix for MEG
  if parms.useMEG_flag
    if parms.verbose
      fprintf('%s: calculating spherical shell gain matrix for MEG...\n',...
        mfilename);
    end;
    psensor = [];
    psensor.sensor = ...
      bsxfun(@minus,parms.MEG_info.intpnt_loc*100,parms.cen_sph/10)'; % in cm
    psensor.orient = parms.MEG_info.intpnt_ori';
    psensor.weight = 1000;
    G_xyz_inp.meg = sarvas5(source_locs'*100,psensor,-1); % fT and fT/cm
  end;

  % calculate gain matrix for EEG
  if parms.useEEG_flag
    if parms.verbose
      fprintf('%s: calculating spherical shell gain matrix for EEG...\n',...
        mfilename);
    end;
    G_xyz_inp.eeg = gain3(source_locs,parms.EEG_info.intpnt,...
      fliplr(parms.radii/1000),fliplr(parms.conductivities),1);
  end;

  % convert integration points into channels
  if parms.verbose
    fprintf('%s: converting integration points to channels...\n',mfilename);
  end;
  if parms.useMEG_flag
    G_xyz_MEG = ts_gain_intpnt2chan(parms.MEG_info.Coil,G_xyz_inp.meg);
  else
    G_xyz_MEG = [];
  end;
  if parms.useEEG_flag
    G_xyz_EEG = ts_gain_intpnt2chan(parms.EEG_info.sensor,G_xyz_inp.eeg);
  else
    G_xyz_EEG = [];
  end;

  % combine MEG and EEG gain matrices
  G_xyz = combine_gain_matrices(parms,G_xyz_MEG,G_xyz_EEG);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G_xyz = combine_gain_matrices(parms,G_xyz_MEG,G_xyz_EEG)
  G_xyz=zeros(parms.nsensors,3*parms.nsources);
  if parms.useMEG_flag
    G_xyz(parms.MEG_chans,:)=G_xyz_MEG*parms.MEG_sf;
  end;
  if parms.useEEG_flag
    G_xyz(parms.EEG_chans,:)=G_xyz_EEG*parms.EEG_sf;
  end;
  if parms.verbose
    fprintf('%s: finished calculating gain matrix\n',mfilename);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

