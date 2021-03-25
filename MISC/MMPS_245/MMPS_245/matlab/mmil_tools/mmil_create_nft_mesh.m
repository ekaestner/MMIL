function mmil_create_nft_mesh(fname_in,varargin)
%function mmil_create_nft_mesh(fname_in,[options])
%
% Purpose: create brain and skull surface meshses
%   using NFT (Neuroelectromagnetic Forward Head Modeling Toolbox)
%
% Required Input:
%   fname_in: name of input mat file containing segmentation volumes
%
% Optional Input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mesh'}
%   'nmesh': number of vertices for initial mesh creation
%     {default = 10000}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Acknowledgement: based on NFT's nft_mesh_generation
%
% Created:  09/28/13 by Don Hagler
% Last Mod: 06/16/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output files
if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% save segmentation volumes in raw format
parms = save_raw_volumes(parms);

% create initial meshes from volumes
init_mesh(parms);

% coarsen and correct meshes
correct_mesh(parms);

% final mesh correction
final_correct_mesh(parms);

% delete temporary files
cleanup_files(parms);

% save output file
save_output(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  parms = mmil_args2parms(options,{,...
    'fname_in',fname_in,[],...
  ...
    'outdir',pwd,[],...
    'outstem','mesh',[],...
    'nmesh',10000,[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'lmr_flag',true,[false true],...
    'quad_flag',false,[false true],...
    'ratio_lmr',2.1,[],...
    'segvol_names',{'scalpmask','brainmask',...
                    'innerskullmask','outerskullmask'},[],...
    'layer_names',{'Scalp','Brain','Csf','Skull'},[],...
  });
  % set number of layers
  parms.nlayers = length(parms.layer_names);
  % check input file exists
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  % load mesh configuration for NFT path names
  parms.conf = nft_get_config;
  % add '/' to end of outdir because NFT functions expect it
  lof = length(parms.outdir);
  if parms.outdir(lof) ~= filesep
    parms.outdir(lof+1) = filesep;
  end
  % set output file name
  parms.fname_out = [parms.outdir parms.outstem '.mat'];
  % create output directory
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = save_raw_volumes(parms)
  % load input mat file containing segmentation volumes
  if parms.verbose
    fprintf('%s: loading input segmentation volumes...\n',mfilename);
  end;
  Segm = [];
  load(parms.fname_in);
  if isempty(Segm), error('missing Segm struct in %s',parms.fname_in); end;

  % save volume dimensions
  [parms.K,parms.L,parms.M] = size(Segm.(parms.segvol_names{1}));
  transform = size(Segm.(parms.segvol_names{1}))/2;
  save([parms.outdir 'transform'],'transform','-ascii'); 

  % save volumes in raw format
  if parms.verbose
    fprintf('%s: saving volumes in raw...\n',mfilename);
  end;
  for k=1:parms.nlayers
    fname_raw = [parms.outdir parms.layer_names{k} '.raw'];
    write_raw_vol(Segm.(parms.segvol_names{k}),fname_raw);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_raw_vol(vol,fname_out);
  vol = single(vol);
  maxval = max(vol(:));
  vol = vol * 20 / maxval;
  fid=fopen(fname_out, 'w+');
  fwrite(fid, vol, 'uint8');
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function init_mesh(parms)
  if parms.verbose
    fprintf('%s: triangulating volumes...\n',mfilename);
  end;
  for k=1:parms.nlayers
    tt = parms.layer_names{k};
    a = sprintf('"%s" -t 10 -dr1 "%s%s.raw" %d %d %d -f "%s%s.asc" -ot',...
      parms.conf.asc,parms.outdir,tt,parms.K,parms.L,parms.M,parms.outdir,tt);
    [status, result] = unix(a);
    if status
      error('Mesh_Generation:system','Failed to execute: %s',result);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function correct_mesh(parms)
  % generate a file for StepSc.txt for coarsening and smoothing
  f=fopen(sprintf('%sStepSc.txt',parms.outdir), 'w');
  fprintf(f, 'correct 5\n');
  fprintf(f, 'smooth 1\n');
  fprintf(f, 'correct 5\n');
  fprintf(f, 'improve 2 0.1 0.05\n');
  fprintf(f, 'improve 2 0.3 0.2\n');
  fprintf(f, 'correct 5\n');
  fprintf(f, 'prune all\n');
  fprintf(f, 'save %sScS.smf\n',parms.outdir);
  fprintf(f, 'quit\n');
  fclose(f);

  % generate a file for StepSc2.txt for final improvement
  f=fopen(sprintf('%sStepSc2.txt',parms.outdir), 'w');
  fprintf(f, 'correct 2\n');
  fprintf(f, 'split intersect\n');
  fprintf(f, 'fill holes\n');
  fprintf(f, 'correct 5\n');
  fprintf(f, 'improve 2 0.1 0.05\n');
  fprintf(f, 'improve 2 0.3 0.2\n');
  fprintf(f, 'correct 5\n');
  fprintf(f, 'split intersect\n');
  fprintf(f, 'correct 2\n');
  fprintf(f, 'prune all\n');
  fprintf(f, 'split intersect\n');
  fprintf(f, 'fill holes\n');
  fprintf(f, 'correct 2\n');
  fprintf(f, 'prune all\n');
  fprintf(f, 'improve 2 0.1 0.05\n');
  fprintf(f, 'correct 2\n');
  fprintf(f, 'save %sScS.smf\n',parms.outdir);
  fprintf(f, 'quit\n');
  fclose(f);

  csi = 300000; i = 1;
  while parms.nmesh < csi
    i = i+1;
    csi(i) = round(csi(i-1) / 1.5);
  end
  csi(i) = parms.nmesh;
  nsteps = length(csi); % number coarsening steps

  % coarsening and correcting
  if parms.verbose
    fprintf('%s: coarsening and correcting...\n',mfilename);
  end;
  for k=1:parms.nlayers
    tt = parms.layer_names{k};
    a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.asc"',...
      parms.conf.showmesh,parms.outdir,parms.outdir,tt);
    [status, result] = unix(a);
    if status
      error('Mesh_Generation:system','Failed to execute: %s',result);
    end
    for iter = 1:nsteps
      a = sprintf('"%s" -c 0.5 -m 5 -o "%s%s.smf" -t %d "%sScS.smf"',...
        parms.conf.qslim,parms.outdir,tt,csi(iter),parms.outdir);
      [status, result] = unix(a);
      if status
        error('Mesh_Generation:system','Failed to execute: %s',result);
      end
      a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.smf"',...
        parms.conf.showmesh,parms.outdir,parms.outdir,tt);
      [status, result] = unix(a);
      if status
        error('Mesh_Generation:system','Failed to execute: %s',result);
      end
      if parms.quad_flag & iter == 6
        copyfile([parms.outdir 'ScS.smf'],[parms.outdir tt 'f.smf'])
      end
    end
    a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"',...
      parms.conf.showmesh,parms.outdir,parms.outdir,tt);
    [status, result] = unix(a);
    if status
      error('Mesh_Generation:system','Failed to execute: %s',result);
    end
    movefile([parms.outdir 'ScS.smf'],[parms.outdir tt '.smf']);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function final_correct_mesh(parms)
  if parms.verbose
    fprintf('%s: final mesh correction...\n',mfilename);
  end;
  % final correction with Matlab functions
  mesh_final_correction(parms.outdir,parms.nlayers);
  % local mesh refinement
  if parms.lmr_flag
    if parms.verbose
      fprintf('%s: local mesh refinement...\n',mfilename);
    end;
    mesh_local_refinement(parms.outdir,parms.nlayers,parms.ratio_lmr);
    handles.param.lmr = parms.ratio_lmr;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_files(parms)
  % delete unnecessary files (.raw, .asc, Scs.smf, and StepSc.txt)
  delete([parms.outdir 'StepSc.txt']);
  delete([parms.outdir 'StepSc2.txt']);
  for k=1:parms.nlayers
    tt = parms.layer_names{k};
    delete([parms.outdir tt '.raw']);
    delete([parms.outdir tt '.asc']);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_output(parms)
  if parms.verbose
    fprintf('%s: mesh generation complete. saving to %s...\n',...
      mfilename,parms.fname_out);
  end;
  % generate nl-layer head model
  mesh_read_write(parms.outdir,parms.outstem,parms.nlayers,parms.quad_flag);
  % record mesh parameters
  parameters = [];
  parameters.no_layers = parms.nlayers;
  parameters.linear = 1;
  parameters.mesh_name = parms.outstem;
  if parms.lmr_flag
    parameters.lmr = parms.ratio_lmr;
  end;
  save(parms.fname_out,'-STRUCT','parameters')
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

