function [M_Atl_to_Subj,regStruct] = mmil_dct_reg_atlas(vol,varargin)
%function [M_Atl_to_Subj,regStruct] = mmil_dct_reg_atlas(vol,[options])
%
% Purpose: register volume to atlas using Discrete Cosine Transform morph
%
% Required Input:
%   vol: input brain volume (ctx format)
%
% Optional Parameters to specify atlas:
%  'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%  'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%
% Optional Parameters to control registration:
%  'smoothflag': [0|1] apply Gaussian smoothing to input volume
%     {default = 1}
%  'sampling': sampling rate (num voxels) in each dimension
%     {default = [4 4 4]}
%  'nK': number of DCT basis functions in each dimension
%     {default = [5 5 5]}
%  'tstep': size of translation step (mm)
%     {default = 0.5}
%  'astep': size of angle step (degrees)
%     {default = 0.25}
%  'scales': vector of scales for multi-scale search
%     {default = [0 83 49 27 16 9 5 3 2 1]}
%  'ns': number of samples
%     {default = 64}
%  'sf': scaling factor
%     {default = 1}
%  'thresh': threshold applied to vol
%     to prevent zero values affecting registration
%     {default = 20}
%
% Optional Parameters:
%  'verbose': [0|1] whether to display progress messages
%     {default = 1}
%  'stdflag': [0|1] whether to use atlas std volume
%     {default = 1}
%  'maskflag': [0|1] whether to use atlas mask volume
%     {default = 1}
%  'stdbgval': value of of std outside brain mask
%     {default = 75}
%
% Output:
%   M_Atl_to_Subj: rigid body registration matrix (atlas to subject)
%   regStruct: structure containing fields:
%     M_atl_to_vol_rb: rigid body registration matrix
%     M_atl_to_vol_af: affine registration matrix
%     VL: volume containing displacements in L-R axis
%     VP: volume containing displacements in A-P axis
%     VH: volume containing displacements in I-S axis
%     Tr: output from nonlinear registration
%     sf_rb: scaling factor for rigid body registration
%     bconverged
%     bcombined
%     min_cost_rb
%     min_cost_af
%     min_cost_m
%     volm: atlas mean volume struct
%     bmesh: outer brain surface struct for atlas (not subj)
%
% See Also:
%   mmil_dct_brainmask: runs this function and create subject brain mask
%   mmil_warp_to_atlas: runs this function and calculates voxel mapping
%
% Created:  11/20/09 by Don Hagler
% Last Mod: 09/08/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...
...
  'smoothflag',true,[false true],...
  'sampling',[4 4 4],[],...
  'nK',[5 5 5],[],...
  'tstep',0.5,[],...
  'astep',0.25,[],...
  'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  'ns',64,[],...
  'sf',1,[],...
  'thresh',20,[0,Inf],...
...
  'verbose',true,[false true],...
  'stdflag',true,[false true],...
  'maskflag',true,[false true],...
  'stdbgval',75,[],...
...
  'morph_tags',{'volstd','volmask','stdbgval','smoothflag','sampling','nK',...
                'tstep','astep','scales','ns','sf','thresh','verbose'},[],...
});
M_Atl_to_Subj = [];
regStruct = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load atlas
if isempty(parms.atlasdir)
  parms.atlasdir = [getenv('MMPS_DIR') '/atlases'];
end;
if mmil_isrelative(parms.atlasname)
  parms.atlasname = [parms.atlasdir '/' parms.atlasname];
end;
parms.atlasname = sprintf('%s.mat',parms.atlasname);
if ~exist(parms.atlasname,'file')
  error('file %s not found\n',parms.atlasname);
end;
atlas = load(parms.atlasname);
volm = atlas.volm;
bm = mmil_getfield(atlas,'bm');
volmask = mmil_getfield(atlas,'volmask');
if parms.stdflag
  parms.volstd = mmil_getfield(atlas,'volstd');
else
  parms.volstd = [];
end;
clear atlas;

% get mask from atlas brain mesh
if parms.maskflag
  if ~isempty(volmask)
    if parms.verbose
      fprintf('%s: using atlas brain mask...\n',mfilename);
    end
    parms.volmask = volmask;
  elseif ~isempty(bm)
    if parms.verbose
      fprintf('%s: getting atlas brain mask from mesh...\n',mfilename);
    end
    [tmp,parms.volmask] = getmaskvol(volm,bm,eye(4,4));
    clear tmp;
  else
    if parms.verbose
      fprintf('%s: no atlas brain mask to use...\n',mfilename);
    end
    parms.volmask = [];
  end;
else
  parms.volmask = [];
end;

args = mmil_parms2args(parms,parms.morph_tags);
regStruct = mmil_dct_morph(vol,volm,args{:});

M_Atl_to_Subj = regStruct.M_atl_to_vol_rb;
regStruct.volm = volm;
regStruct.bmesh = bm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

