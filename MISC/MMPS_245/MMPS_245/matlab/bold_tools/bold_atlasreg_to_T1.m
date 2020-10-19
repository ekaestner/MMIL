function [M_T1_to_BOLD,M_atlas_to_BOLD,M_atlas_to_T1] = ...
  bold_atlasreg_to_T1(vol_BOLD,vol_T1)
%function [M_T1_to_BOLD,M_atlas_to_BOLD,M_atlas_to_T1] = ...
%  bold_atlasreg_to_T1(vol_BOLD,[vol_T1])
%
% Required Input:
%   vol_BOLD: BOLD b=0 3D volume
%
% Optional Input:
%   vol_T1: T1-weighted 3d volume (structural scan)
%     if empty or omitted, will register to BOLD atlas
%       (multi-subject average of BOLD images aligned to T1 atlas)
%     if vol_T1 is supplied, will register BOLD to T1 via BOLD and T1 atlases
%
% Output:
%   M_T1_to_BOLD: registration matrix from T1 to BOLD
%
% Note: vol_BOLD and vol_T1 should be in ctx format
%   use vol_ctx=ctx_mgh2ctx(vol,M);
%
% Created : 06/13/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

M_T1_to_BOLD = [];
M_atlas_to_BOLD = [];
M_atlas_to_T1 = [];


M_atlasT1_to_atlasBOLD = [     1     0     0     -2; ...
                               0     1     0     0; ...
                               0     0     1     0; ...
                               0     0     0     1    ];

if ~mmil_check_narg(nargin,1), return; end;
if ~exist('vol_T1','var'), vol_T1 = []; end;

dir_atlas = [getenv('MMPS_DIR') '/atlases'];
fname_T1_atlas = sprintf('%s/BOLD_reg_targ_T1.mgh',dir_atlas);
fname_BOLD_atlas = sprintf('%s/BOLD_reg_targ_BOLD.mgh',dir_atlas);
fname_T1_mask_atlas = sprintf('%s/BOLD_reg_targ_mask.mgh',dir_atlas);
fname_BOLD_mask_atlas = sprintf('%s/BOLD_reg_targ_BOLD_mask.mgh',dir_atlas);

if ~exist(fname_BOLD_atlas,'file')
  error('BOLD atlas file %s not found',fname_BOLD_atlas);
end;
if ~isempty(vol_T1) && ~exist(fname_T1_atlas,'file')
  error('T1 atlas file %s not found',fname_T1_atlas);
end;

fprintf('%s: registering BOLD to atlas...\n',mfilename);
[vol_BOLD_atlas,mr_parms,volsz] = ctx_load_mgh(fname_BOLD_atlas);
vol_BOLD_mask = ctx_load_mgh(fname_BOLD_mask_atlas);
M_atlasBOLD_to_BOLD = ...
  rbreg_vol2vol_icorr_mask(vol_BOLD_atlas,vol_BOLD,vol_BOLD_mask,0,4);


if ~isempty(vol_T1)
  fprintf('%s: registering T1 to atlas...\n',mfilename);
  [vol_T1_atlas,mr_parms,volsz] = ctx_load_mgh(fname_T1_atlas);
  vol_T1_mask = ctx_load_mgh(fname_T1_mask_atlas);
  M_atlasT1_to_T1 = ...
    rbreg_vol2vol_icorr_mask(vol_T1_atlas,vol_T1,vol_T1_mask,0,4);
  M_T1_to_BOLD = M_atlasBOLD_to_BOLD*M_atlasT1_to_atlasBOLD*inv(M_atlasT1_to_T1);
else
  M_T1_to_BOLD = M_atlasBOLD_to_BOLD*M_atlasT1_to_atlasBOLD;
end;

