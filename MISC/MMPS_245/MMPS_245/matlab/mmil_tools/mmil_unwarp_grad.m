function vol_uw = mmil_unwarp_grad(vol,varargin)
%function vol_uw = mmil_unwarp_grad(vol,varargin)
%
% Usage:
%   vol_uw = mmil_unwarp_grad(vol,'key1', value1,...);
%
% Required Input:
%   vol: volume in ctx structure to be grad unwarped
% 
% Optional Input:
%   'gwtype': gradient type number
%     0:  Siemens Sonata
%     1:  Siemens Allegra
%     2:  GE BRM
%     3:  GE CRM
%     4:  Siemens Avanto
%     5:  Siemens AXXess/Espree
%     6:  Siemens Quantum/Symphony
%     7:  GE Twin Speed Whole Body
%     8:  GE Twin Speed Zoom
%     9:  GE mr450 or mr750
%     10: GE MR750W
%     11: Siemens Skyra
%     12: Siemens Connectome Skyra
%     13: Siemens Prisma
%     {default = 0}
%   'jacobian_flag': [0|1] whether to apply jacobian
%     {default = 1}
%   'jacobian_dist': distance for jacobian estimate
%     {default = 1}
%   'unwarpflag': [0|1|2]
%     0: unwarp 3D
%     1: unwarp through plan only
%     2: unwarp inplane only
%     {default = 0}
%   'interpm': interpolation method number
%     0: Nearest Neighbor 1: Linear  2: Cubic
%     3: Key's spline 4: Cubic spline 5: Hamming Sinc
%     {default = 2}
%   'bclamp': [0|1] whether to set negative value to zero
%     {default = 1}
%   'isoctrflag': [0|1] whether to adjust for isocenter coordinates
%     {default = 1}
%
%  See: ctx_load_mgh
%
% Created:  07/19/08 by Don Hagler
% Prev Mod: 11/16/10 by Don Hagler
% Last Mod: 01/08/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_uw = [];
if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms( varargin, {...
  'gwtype',0,[0:13],...
  'jacobian_flag',true,[false true],...
  'jacobian_dist',1,[1 10],...
  'unwarpflag',0,[0:2],...
  'interpm',2,[1:5],...
  'bclamp',true,[false true],...
  'isoctrflag',true,[false true],...
  'npad',5,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input volume (must be ctx format)
if ~isfield(vol,'imgs')
  error('input vol must be ctx structure')
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pad volume with extra slices, copying edge slices

if parms.interpm==0
  tmp_vol = vol;
else
  % interpolation results in loss of edge slices
  %   unless we first pad the volume with extra slices
  nx = size(vol.imgs,1);
  ny = size(vol.imgs,2);
  nz = size(vol.imgs,3);
  % allocate volume with extra slices
  tmp_vol = vol;
  tmp_vol.imgs = zeros(nx,ny,nz+2*parms.npad);
  tmp_vol.imgs(:,:,parms.npad+1:nz+parms.npad) = vol.imgs;
  for i=1:parms.npad
    tmp_vol.imgs(:,:,i) = vol.imgs(:,:,1); % pad with copies of first slice
  end;
  for i=nz+parms.npad+1:nz+2*parms.npad
    tmp_vol.imgs(:,:,i) = vol.imgs(:,:,nz); % pad with copies of last slice
  end;
  tmp_vol.Mvxl2lph = shift_M(tmp_vol.Mvxl2lph,[nx,ny,nz],[nx,ny,nz+2*parms.npad]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply grad warp

if parms.isoctrflag
  Mvxl2lph = tmp_vol.Mvxl2lph;
  firstLPH = Mvxl2lph*[1; 1; 1; 1];
  lastLPH = Mvxl2lph*[tmp_vol.dimr; tmp_vol.dimc; tmp_vol.dimd; 1];
  tmp_vol.Mvxl2lph(3,4) = tmp_vol.Mvxl2lph(3,4)-(firstLPH(3)+lastLPH(3))/2;
end

vol_uw = vol_unwarp_grad(tmp_vol,parms.gwtype,parms.jacobian_flag,...
  parms.jacobian_dist,parms.unwarpflag,parms.interpm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove extra slices

vol_uw.Mvxl2lph = vol.Mvxl2lph;
vol_uw.imgs = squeeze(vol_uw.imgs(:,:,parms.npad+1:nz+parms.npad));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

function M_shift = shift_M(M,nvox,nvox_shift);
  M_shift = M;
  M_shift(1:3,4) = M_shift(1:3,4) + M(1:3,:)*[nvox/2+1 1]' - M_shift(1:3,:)*[nvox_shift/2+1 1]';
return;

