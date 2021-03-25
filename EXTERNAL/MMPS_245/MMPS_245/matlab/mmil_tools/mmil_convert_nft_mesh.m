function mmil_convert_nft_mesh(fname_in,varargin)
%function mmil_convert_nft_mesh(fname_in,[options])
%
% Purpose: convert NFT meshes to FreeSurfer tri format
%
% Required Input:
%   fname_in: name of input mat file containing brain and skull meshes
%     created by NFT (Neuroelectromagnetic Forward Head Modeling Toolbox)
%
% Optional Input:
%   'outdir': output directory
%     {default = pwd}
%   'mesh_names': cell array of names for each mesh
%     {default = {'scalp','outer_skull','inner_skull', 'brain'}}
%   'mesh_suffix': output file suffix
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  09/28/13 by Don Hagler
% Last Mod: 06/16/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output files
all_exist = check_output(parms);
if all_exist && ~parms.forceflag, return; end;

% load input file
[mesh,shells] = load_mesh_file(parms);

% convert meshes
save_mesh_files(mesh,shells,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  parms = mmil_args2parms(options,{,...
    'fname_in',fname_in,[],...
  ...
    'outdir',pwd,[],...
    'mesh_names',{'outer_scalp','outer_skull','inner_skull', 'brain'},[],...
    'mesh_suffix',[],[],...
    'forceflag',false,[false true],...
  });
  % set number of layers
  parms.nlayers = length(parms.mesh_names);
  % check input file exists
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_exist = check_output(parms)
  all_exist = 1;
  for i=1:parms.nlayers
    fname_out = sprintf('%s/%s%s.tri',...
      parms.outdir,parms.mesh_names{i},parms.mesh_suffix);
    if ~exist(fname_out,'file'), all_exist = 0; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mesh,shells] = load_mesh_file(parms)
  [tpath,tstem] = fileparts(parms.fname_in);
  mesh = bem_load_mesh([tpath '/' tstem]);
  if parms.nlayers > mesh.num_boundaries
    error('mesh contains only %d layers (expected %d)',...
      mesh.num_boundaries,parms.nlayers);
  end;
  % determine which nodes are in which surface
  shells = zeros(mesh.num_nodes,1);
  k = 1;
  for i=1:parms.nlayers
    for j=1:mesh.bnd(i,1)
      shells(mesh.elem(k,:)) = i;
      k = k + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_mesh_files(mesh,shells,parms)
  % separate the surfaces
  surfs = [];
  j = 1;
  for i=1:parms.nlayers
    surfs(i).nfaces = mesh.bnd(i,1);
    k = j + surfs(i).nfaces - 1;
    surfs(i).faces = mesh.elem(j:k,:);
    surfs(i).faces = surfs(i).faces - min(surfs(i).faces(:)) + 1;
    j = k+1;
    surfs(i).vertices = mesh.coord(find(shells==i),:);
    surfs(i).nverts = length(surfs(i).vertices);
  end;
  % transform coordinates
  for i=1:parms.nlayers
    surfs(i).vertices = surfs(i).vertices(:,[1,2,3])-128;
  end;
  % save as tri files
  for i=1:parms.nlayers
    fname_out = sprintf('%s/%s%s.tri',...
      parms.outdir,parms.mesh_names{i},parms.mesh_suffix);
    fs_write_trisurf(surfs(i),fname_out);
  end;
return;

