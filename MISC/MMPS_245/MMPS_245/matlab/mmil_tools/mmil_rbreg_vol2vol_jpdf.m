function [M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_jpdf(vol1,vol2,varargin)
%function [M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_jpdf(vol1,vol2,[options])
%
% Required Input:
%   vol1: 3D volume of reference scan (ctx struct)
%   vol2: 3D volume of scan to be registered
%
% Optional Parameters:
%   'type1': image type for vol1 (reference image)
%     (e.g. 'MPR','HIFA', or 'FSNU' (i.e. Freesurer nu.mgz))
%     {default = 'MPR'} 
%   'type2': image type for vol2
%     (e.g.'BOLD','DTI','GECAL','PET','LOFA', or 'HIFA')
%     {default = 'LOFA'}
%   'volmask': mask for vol1 (e.g. brain mask)
%     if not supplied, will automatically generate one by registering to atlas
%     NOTE: mask generation will only work if vol1 is T1-weighted
%     {default = []}
%   'interpm':  0:Nearest 1:Linear (default) 2:Cubic
%             3:Key's spline 4:Cubic spline 5:Hamming Sinc
%   'Mreg0': 4x4 matrix defining initial estimate of registration matrix
%     {default = identity}
%   'scales': vector of scales (i.e. displacement sizes) for multi-scale search
%     {default:[0 83 49 27 16 9 5 3 2 1]}
%   'jpdfdir': directory for jdpf file (input/output)
%     {default = $MMPS_DIR/atlases/jpdfs}
%   'initflag': [0|1|2] indicate desired behavior for initalizing jpdf
%     0: if jpdf does not already exist, exit with error message
%     1: if jpdf does not already exist (or if forceflag = 0), initialize it
%        (this requires that vol1 and vol2 are already well-aligned)
%     2: use existing jpdf if it exists (otherwise initialize)
%        and reinitialize with values from vol1 and vol2 after registration
%     {default = 0}
%   'forceflag': [0|1] overwrite existing jpdf file
%     ignored if initflag = 0
%     {default = 0}
%
% Output:
%   M_v1_to_v2: registration matrix from vol1 to vol2
%   volmask: brain mask for vol1
%
% NOTE: vol2, vol1, and volmask must be in ctx format
%   use vol_ctx=ctx_mgh2ctx(vol,M);
%
% Based on code by Anders Dale
% Created : 12/16/06 by Don Hagler
% Last Mod: 06/09/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

M_v1_to_v2 = []; volmask = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'type1','MPR',[],...
  'type2','LOFA',[],...
  'volmask',[],[],...
  'interpm',1,[0:5],...
  'Mreg0',eye(4),[],...
  'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  'jpdfdir',[],[],...
  'initflag',false,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (parms.scales(1)~=0) parms.scales=[0,parms.scales]; end;

parms.type1 = upper(parms.type1);
parms.type2 = upper(parms.type2);

if isempty(parms.jpdfdir)
  parms.jpdfdir = sprintf('%s/atlases/jpdfs',getenv('MMPS_DIR'));
end;

parms.jpdffname = sprintf('%s/%s_%s_jpdf.mat',...
  parms.jpdfdir,parms.type1,parms.type2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove existing jpdf file if it exists
if exist(parms.jpdffname,'file') && parms.initflag==1 && parms.forceflag
  delete(parms.jpdffname);
end;

% load jpdf file if it exists
if exist(parms.jpdffname,'file')
  fprintf('%s: using jpdf file %s\n',mfilename,parms.jpdffname);
  load(parms.jpdffname);
  if parms.initflag==1, parms.initflag=0; end;
elseif ~parms.initflag
  error('jpdf file %s not found',parms.jpdffname);
else
  fprintf('%s: initializing jpdf file %s\n',mfilename,parms.jpdffname);
  mmil_mkdir(parms.jpdfdir);
end;

% compute scalefactor for vol1
[hc,bv] = hist(vol1.imgs(:),1000);
if ~parms.initflag
  sf1 = ComputeHistScalefactor(hc,bv,v1_hc,v1_bv,v1_bn0);
  vol1.imgs = sf1*vol1.imgs;
  fprintf('%s: sf1 = %f\n',mfilename,sf1);
else
  sf1 = 1;
  v1_hc = hc;
  v1_bv = bv;
  v1_bn0 = 50;
end

% compute scalefactor for vol2
[hc,bv] = hist(vol2.imgs(:),1000);
if ~parms.initflag
  sf2 = ComputeHistScalefactor(hc,bv,v2_hc,v2_bv,v2_bn0);
  vol2.imgs = sf2*vol2.imgs;
  fprintf('%s: sf2 = %f\n',mfilename,sf2);
else
  sf2 = 1;
  v2_hc = hc;
  v2_bv = bv;
  v2_bn0 = 50;
end

% make brain mask for v1
if isempty(parms.volmask)
  parms.volmask = mmil_dct_brainmask(vol1);
end;
fprintf('\n');

% register and resample vol2 to vol1
if parms.initflag
  vol2_res = vol_resample_pad(vol2,vol1,parms.Mreg0,1);
  vol2_res.imgs(find(isnan(vol2_res.imgs))) = 0;
  [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=...
    vols_jhist_mask_amd(vol1, vol2_res, parms.volmask, 1, eye(4), eye(4));
  jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
end
[M_v1_to_v2, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol1, vol2, parms.volmask,...
  parms.Mreg0,0,4,parms.scales,[],[],[],jpdf_mask,bins1_mask,bins2_mask,parms.interpm);
vol2_res = vol_resample_pad(vol2, vol1, M_v1_to_v2, 1);
vol2_res.imgs(find(isnan(vol2_res.imgs))) = 0;

if parms.initflag
  % calculate jpdf from registered volumes and reregister
  [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]= ...
    vols_jhist_mask_amd(vol1,vol2_res,parms.volmask,1,eye(4),eye(4));
  jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
  parms.Mreg0 = M_v1_to_v2;
  [M_v1_to_v2, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol1,vol2,parms.volmask,...
    parms.Mreg0,0,4,parms.scales,[],[],[],jpdf_mask,bins1_mask,bins2_mask,parms.interpm);
  vol2_res = vol_resample_pad(vol2, vol1, M_v1_to_v2, 1);
  vol2_res.imgs(find(isnan(vol2_res.imgs))) = 0;

  fprintf('%s: saving jpdf file %s\n',mfilename,parms.jpdffname);
  save(parms.jpdffname,'jpdf_mask','bins1_mask','bins2_mask','v2_hc','v2_bv','v2_bn0','v1_hc','v1_bv','v1_bn0');
end

fprintf('%s: min_cost = %f\n',mfilename,min_cost);

if nargout>1
  volmask = parms.volmask;
end;

