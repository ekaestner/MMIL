function [fname_atl,fname_reg,fname_vxl]=mmil_warp_to_atlas(fname,varargin)
%function [fname_atl,fname_reg,fname_vxl]=mmil_warp_to_atlas(fname,[options])
%
% Usage:
%  [fname_atl,fname_reg]=mmil_warp_to_atlas(fname,'key1', value1,...);
%
% Required Input:
%   fname: full path name of input volume file
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
%  'verbose': [0|1] whether to display progress messages
%     {default = 1}
%  'stdflag': [0|1] whether to use atlas std volume
%     {default = 1}
%  'maskflag': [0|1] whether to use atlas mask volume
%     {default = 1}
%  'stdbgval': value of of std outside brain mask
%     {default = 75}
%
% Optional Parameters:
%  'fname_reg': full path name of matlab file containing registration info
%                 (e.g. output of this function)
%     {default = []}
%  'fname_T1': file name of T1-weighted image to use to register to atlas
%                 (instead of fname; ignored if fname_reg is not empty)
%     {default = fname}
%  'outdir': output directory
%     {default = directory of fname (must have permission!)}
%  'outstem': output file stem
%     if empty, will use file stem of input fname
%     {default = []}
%  'outstem_T1': output file stem used in place of stem of fname_T1
%     {default = []}
%  'interpm': interpolation method
%      0 = nearest neighbor, 1 = linear, 2 = cubic
%      3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%     {default = 2}
%  'bclamp': [0|1] whether to set negative values to zero
%     { default = 1 }
%  'vxlmap_flag':  [0|1] whether to generate voxel mapping
%     {default = 0}
%  'forceflag': [0|1] whether to run calculations even if output files exist
%     {default = 0}
%
%  Output:
%   fname_atl: full path of output atlas-space volume file
%   fname_reg: full path of matlab file containing registration info
%   fname_vxl: full path of matlab file containing voxel mapping from
%                 atlas to subject space (generated if vxlmap_flag=1)
%
% Created:  09/02/07 by Don Hagler
% Last Mod: 02/07/13 by Don Hagler
%

%% todo: allow fname_atl, fname_reg, and fname_vxl to be specified
%%       as optional parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
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
  'verbose',true,[false true],...
  'stdflag',true,[false true],...
  'maskflag',true,[false true],...
  'stdbgval',75,[],...
...
  'fname_reg',[],[],...
  'fname_T1',[],[],...
  'outdir',[],[],...
  'outstem',[],[],...
  'outstem_T1',[],[],...
  'interpm',2,[0,5],...
  'bclamp',true,[false true],...
  'padding',1,[1,5],...
  'vxlmap_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'reg_tags',{'atlasdir','atlasname','smoothflag','sampling','nK',...
              'tstep','astep','scales','ns','sf','thresh','verbose',...
              'stdflag','maskflag','stdbgval'},[],...
});

fname_atl = [];
fname_reg = [];
fname_vxl = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_Atlas_to_Subj=[];
regStruct=[];

if ~exist(fname,'file'), error('file %s not found',fname); end;
[fpath,fstem,fext] = fileparts(fname);
if isempty(parms.outdir), parms.outdir = fpath; end;
if isempty(parms.outstem), parms.outstem = fstem; end;

if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to create outdir %s:\n%s',parms.outdir,msg);
   end;
end;

if isempty(parms.fname_T1)
  parms.fname_T1 = fname;
elseif ~exist(parms.fname_T1,'file')
  error('file %s not found',parms.fname_T1);
end;

if ~isempty(parms.fname_reg) && ~exist(parms.fname_reg,'file')
  error('file %s not found',parms.fname_reg);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% register parms.fname_T1 to atlas

[fpath_T1,fstem_T1,fext_T1] = fileparts(parms.fname_T1);
if isempty(parms.outstem_T1)
  parms.outstem_T1 = fstem_T1;
end;
if ~isempty(parms.fname_reg)
  fname_reg = parms.fname_reg;
else
  fname_reg = sprintf('%s/%s_dctReg2Atlas.mat',...
    parms.outdir,parms.outstem_T1);
  if ~exist(fname_reg,'file') || parms.forceflag
    fprintf('%s: loading subj volume %s...\n',mfilename,parms.fname_T1);
    vol = ctx_load_mgh(parms.fname_T1);
    M_Subj = M_LPH_TO_RAS * vol.Mvxl2lph;
    volsz_Subj = size(vol.imgs);
    fprintf('%s: registering volume %s to atlas %s...\n',...
      mfilename,parms.fname_T1,parms.atlasname);
    args = mmil_parms2args(parms,parms.reg_tags);
    [M_Atlas_to_Subj,regStruct] = mmil_dct_reg_atlas(vol,args{:});
    save(fname_reg,'M_Atlas_to_Subj','regStruct','M_Subj','volsz_Subj');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warp fname to atlas

fname_atl = sprintf('%s/%s_atlas%s',parms.outdir,parms.outstem,fext);
if ~exist(fname_atl,'file') || parms.forceflag
  fprintf('%s: loading volume %s...\n',mfilename,fname);
  [vol,mr_parms] = ctx_load_mgh(fname);
  if ~exist('regStruct','var') | isempty(regStruct)
    load(fname_reg);
  end;
  fprintf('%s: warping volume %s to atlas...\n',mfilename,fname);
  vol_atl = volMorph(regStruct.volm, vol,...
    regStruct.VL, regStruct.VP, regStruct.VH,...
    parms.interpm, parms.padding, parms.bclamp);
  ctx_save_mgh(vol_atl,fname_atl,mr_parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate voxel map for morphing from atlas to subj

if parms.vxlmap_flag
  fname_vxl = sprintf('%s/%s_dctReg2Atlas_vxlmap.mat',...
    parms.outdir,parms.outstem_T1);
  clear vxlmap vol_atl;
  if ~exist(fname_vxl,'file') || parms.forceflag
    if ~exist('vol_atl','var') | isempty(vol_atl)
      fprintf('%s: loading volume %s...\n',mfilename,fname_atl);
      [vol_atl,M,mr_parms,volsz_atl] = fs_load_mgh(fname_atl,[],1);
      vol_atl = ctx_mgh2ctx(vol_atl,M);
    end;
    if ~exist('regStruct','var') | isempty(regStruct)
      load(fname_reg);
    end;
    vol_atl.maxI = max(vol_atl.imgs(:));
    vol_atl.minI = min(vol_atl.imgs(:));
    regStruct.range = [1 size(vol_atl.imgs,1);...
                       1 size(vol_atl.imgs,2);...
                       1 size(vol_atl.imgs,3)];
    fprintf('%s: generating voxel map for morphing from atlas to subject...\n',...
      mfilename);
    vxlmap = getDCTVxlMapping_amd(vol_atl, regStruct);
    save(fname_vxl,'vxlmap');
  end;
end;
