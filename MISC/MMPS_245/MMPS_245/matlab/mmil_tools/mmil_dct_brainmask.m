function [vol_out,M_Atl_to_Subj] = mmil_dct_brainmask(vol_in,varargin)
%function [vol_out,M_Atl_to_Subj] = mmil_dct_brainmask(vol_in,[options])
%
% Purpose: create mask from T1-weighted image volume using dct morph to atlas
%
% Usage: [vol_out,M_Atl_to_Subj] = mmil_dct_brainmask(vol_in,'key1', value1,...);
%
% Required Input:
%   vol_in: input mask volume (ctx format)
%     Can be [] if 'fname_in' is specified
%
% Optional Parameters:
%  'fname_in': full path of input image volume file  (mgh/mgz format)
%    If supplied, vol_in is ignored (can be empty []) and fname_in is loaded instead
%    {default: []}
%  'fname_reg': output/input file name of mat file containing DCT regStruct
%    If exists, will load registration and create brain mask
%    If not exists, will calculate registration and save to this file name
%    If empty, will calcualte registration and not save
%    {default: []}
%  'fname_mask': output file name for brain mask (mgh/mgz format)
%    If empty or ommitted, mask will not be saved
%    {default: []}
%  'fname_brain': output file name for masked brain volume (mgh/mgz format)
%    If empty or ommitted, mask brain volume will not be saved
%    {default: []}
%  'fname_surf': output file name for brain mesh (tri file format)
%    If empty or ommitted, brain mesh will not be saved
%    {default: []}
%  'smooth': smoothing sigma (mm) applied to mask
%    {default: 5}
%  'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%  'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%  'forceflag': [0|1] whether to overwrite output files if they already exists
%    {default: 0}
%
% Output:
%   vol_out: output mask volume (ctx format)
%   M_Atl_to_Subj: rigid body registration matrix to atlas
%
% See Also:
%   mmil_aseg_brainmask: to create mask from FreeSurfer aseg
%   mmil_dilate_mask: to create mask from directly specified aseg volume
%   mmil_masksurf: to create surface mesh from mask
%
% Created:  02/02/07 by Don Hagler
% Last Mod: 02/28/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_in',[],[],...
  'fname_reg',[],[],...
  'fname_mask',[],[],...
  'fname_brain',[],[],...
  'fname_surf',[],[],...
  'smooth',5,[0,100],...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...
  'forceflag',false,[false true],...
});
vol_out = [];
M_Atl_to_Subj = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.atlasdir)
  parms.atlasdir = [getenv('MMPS_DIR') '/atlases'];
end;
if mmil_isrelative(parms.atlasname)
  parms.atlasname = [parms.atlasdir '/' parms.atlasname];
end;

runflag = 0;
if isempty(parms.fname_mask) & isempty(parms.fname_brain) & isempty(parms.fname_surf)
  fprintf('%s: WARNING: no output names (fname_mask, fname_brain, fname_surf) specified \n',...
    mfilename);
  runflag = 1;
else 
  if ~isempty(parms.fname_mask) & (~exist(parms.fname_mask,'file') | parms.forceflag)
    runflag = 1;
  end;
  if ~isempty(parms.fname_brain) & (~exist(parms.fname_brain,'file') | parms.forceflag)
    runflag = 1;
  end;
  if ~isempty(parms.fname_surf) & (~exist(parms.fname_surf,'file') | parms.forceflag)
    runflag = 1;
  end;
end;

if ~isempty(parms.fname_mask) && exist(parms.fname_mask,'file') &&...
   ~parms.forceflag && nargout>0
  vol_out = ctx_load_mgh(parms.fname_mask);
end;

if ~runflag, return; end;


if ~isempty(parms.fname_in)
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [vol_in,mr_parms] = ctx_load_mgh(parms.fname_in);
else
  mr_parms = [];
end;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parms.fname_reg) || ~exist(parms.fname_reg,'file') || parms.forceflag
  % get mask from atlas morph
  [M_Atl_to_Subj,regStruct] = mmil_dct_reg_atlas(vol_in,...
    'atlasdir',parms.atlasdir,'atlasname',parms.atlasname);
  if ~isempty(parms.fname_reg)
    save(parms.fname_reg,'M_Atl_to_Subj','regStruct');
  end;
elseif ~isempty(parms.fname_reg) && exist(parms.fname_reg,'file')
  load(parms.fname_reg);
end;

% deform atlas brain mesh to subject, create mask volume
brainmesh = regStruct.bmesh;
lph = brainmesh.vertices;
lph(:,4)=1;
dL = vol_getvxlsval(lph, regStruct.VL, eye(4,4));
dP = vol_getvxlsval(lph, regStruct.VP, eye(4,4));
dH = vol_getvxlsval(lph, regStruct.VH, eye(4,4));
brainmesh.vertices = brainmesh.vertices+[dL dP dH];
brainmesh.color='r';
[vol_brain, vol_mask] = getmaskvol(vol_in,brainmesh,eye(4));

% smooth mask
if parms.smooth>0
  ps = fs_voxsize(vol_mask.Mvxl2lph); % voxel dimensions in mm
  smoothvec = parms.smooth./ps;
  vol_mask.imgs = mmil_smooth3d(vol_mask.imgs,...
    smoothvec(1),smoothvec(2),smoothvec(3));
  ind = find(vol_mask.imgs<=0);
  vol_mask.imgs(ind)=0;  
end;

% save masked brain vol
if ~isempty(parms.fname_brain) & (~exist(parms.fname_brain) | parms.forceflag)
  vol_brain.imgs = vol_in.imgs.*vol_mask.imgs;
  ctx_save_mgh(vol_brain,parms.fname_brain,mr_parms);
end;

% save mask vol
if ~isempty(parms.fname_mask) & (~exist(parms.fname_mask) | parms.forceflag)
  ctx_save_mgh(vol_mask,parms.fname_mask,mr_parms);
end;

% save surface mesh
if ~isempty(parms.fname_surf) & (~exist(parms.fname_surf) | parms.forceflag)
  skull_surf = [];
  skull_surf.nverts = brainmesh.V;
  skull_surf.nfaces = brainmesh.F;
  skull_surf.faces = brainmesh.faces;
  skull_surf.vertices = brainmesh.vertices;
  skull_surf.vertices(:,1) = -skull_surf.vertices(:,1) + vol.lphcent(1);
  skull_surf.vertices(:,2) = -skull_surf.vertices(:,2) + vol.lphcent(2);
  skull_surf.vertices(:,3) = skull_surf.vertices(:,3) + vol.lphcent(3);
  fs_write_trisurf(skull_surf,parms.fname_surf);
end;

vol_out = vol_mask;

