function [M_S_to_F,RegInfo,volmask]=mmil_atlas_jpdfreg_func2struct(vol_F,vol_S,varargin)
%function [M_S_to_F,RegInfo,volmask]=mmil_atlas_jpdfreg_func2struct(vol_F,vol_S,[options]);
%
% Required Input:
%   vol_F: 3D volume (i.e. first frame) of functional scan
%   vol_S: 3D volume of structural scan
%
% Optional Input:
%  'ftype': either 'BOLD', 'DTI', 'PET', or 'GECAL'
%    {default = 'BOLD'}
%  'stype': either 'MPR', 'HIFA', 'LOFA', or 'FSNU'
%    {default = 'FSNU'} (i.e. Freesurer nu.mgz)
%  'volmask': mask for S volume (e.g. brain mask)
%    {default = []}
%  'reinitflag': [0|1] toggle whether to replace existing jpdf
%     with new one calculated after registering data
%    {default = 0}
%
% Note: vol_F, vol_S, and volmask should be in ctx format
% Note: if jpdf does not yet exist, must initialize by
%         running this program with volumes that are already aligned
%
% Based on code by Anders Dale
%
% Created : 12/16/06 by Don Hagler
% Last Mod: 04/17/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_S_to_F = []; RegInfo = []; volmask = [];

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'ftype','BOLD',{'BOLD','DTI','PET','GECAL'},...
  'stype','FSNU',{'MPR','HIFA','FLASHhi','LOFA','FLASHlo','FSNU'},...
  'volmask',[],[],...
  'reinitflag',false,[false true],...
  ...
  'atlasdir',[getenv('MMPS_DIR') '/atlases/jpdfs'],[],...
});

parms.ftype = upper(parms.ftype);
parms.stype = upper(parms.stype);

if strcmp(parms.stype,'FLASHHI')
  parms.stype = 'HIFA';
end;
if strcmp(parms.stype,'FLASHLO')
  parms.stype = 'LOFA';
end;

parms.jpdffname = sprintf('%s/%s_%s_jpdf.mat',...
  parms.atlasdir,parms.stype,parms.ftype);

% load jpdf file if it exists
if exist(parms.jpdffname,'file')
  fprintf('%s: using %s\n',mfilename,parms.jpdffname);
  load(parms.jpdffname);
  parms.initflag = 0;
else
  fprintf('%s: initializing %s\n',mfilename,parms.jpdffname);
  parms.initflag = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute F scalefactor
[hc,bv] = hist(vol_F.imgs(:),1000);
if ~parms.initflag
  sf_F = ComputeHistScalefactor(hc,bv,Fhc,Fbv,Fbn0);
  vol_F.imgs = sf_F*vol_F.imgs;
  fprintf(1,'%s: sf_F = %f\n',mfilename,sf_F);
else
  sf_F = 1;
  Fhc = hc;
  Fbv = bv;
  Fbn0 = 50;
%  figure(1); clf; plot(cumsum(Fhc(Fbn0:end))/sum(Fhc(Fbn0:end))); axis tight;
end

% compute S scalefactor
[hc,bv] = hist(vol_S.imgs(:),1000);
if ~parms.initflag
  sf_S = ComputeHistScalefactor(hc,bv,Shc,Sbv,Sbn0);
  vol_S.imgs = sf_S*vol_S.imgs;
  fprintf(1,'%s: sf_S = %f\n',mfilename,sf_S);
else
  sf_S = 1;
  Shc = hc;
  Sbv = bv;
  Sbn0 = 50;
%  figure(2); clf; plot(cumsum(Shc(Sbn0:end))/sum(Shc(Sbn0:end))); axis tight;
end

% make brain mask for S
if isempty(parms.volmask)
  parms.volmask = mmil_dct_brainmask(vol_S);
end;

% register and resample F to S
Mreg0 = eye(4);
if parms.initflag
  vol_F_res = vol_resample_pad(vol_F, vol_S, Mreg0, 1);
  vol_F_res.imgs(find(isnan(vol_F_res.imgs))) = 0;
  [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=...
    vols_jhist_mask_amd(vol_S,vol_F_res,parms.volmask,1,eye(4),eye(4));
  jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
end
[M_S_to_F, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol_S,vol_F,...
  parms.volmask,Mreg0,1,4,[],[],[],[],jpdf_mask,bins1_mask,bins2_mask);
vol_F_res = vol_resample_pad(vol_F, vol_S, M_S_to_F, 1);
vol_F_res.imgs(find(isnan(vol_F_res.imgs))) = 0;

if parms.initflag | parms.reinitflag
  % calculate jpdf from registered S & F volumes and reregister
  for iter = 1:2
    [sum_log10_val,jpdf_mask,bins1_mask,bins2_mask,jentropy] = ...
      vols_jhist_mask_amd(vol_S, vol_F_res,parms.volmask,1,eye(4),eye(4));
    jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
    Mreg0 = M_S_to_F;
    [M_S_to_F, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol_S,vol_F,...
      parms.volmask,Mreg0,1,4,[],[],[],[],jpdf_mask,bins1_mask, bins2_mask);
    vol_F_res = vol_resample_pad(vol_F, vol_S, M_S_to_F, 1);
    vol_F_res.imgs(find(isnan(vol_F_res.imgs))) = 0;
  end
end;

% save jdpf file
if parms.initflag | parms.reinitflag
  fprintf('%s: saving jpdf file %s\n',mfilename,parms.jpdffname);
  save(parms.jpdffname,'jpdf_mask','bins1_mask','bins2_mask','Fhc','Fbv','Fbn0','Shc','Sbv','Sbn0');
end

% fill output S with registration info
RegInfo.M_S_to_F = M_S_to_F;
RegInfo.min_cost = min_cost;
RegInfo.jpdf_mask = jpdf_mask;
RegInfo.bins1_mask = bins1_mask;
RegInfo.bins2_mask = bins2_mask;
RegInfo.sf_F = sf_F;
RegInfo.sf_S = sf_S;

volmask = parms.volmask;
