function MFmeas = rsi_calc_MFmeas(MFfit,varargin)
%function MFmeas = rsi_calc_MFmeas(MFfit,[options])
%
% Purpose: calculate various measures based on multi-compartment
%   fiber orientation density (FOD) function fit
%   of multi-shell diffusion MRI data 
%
% Required Parameters:
%   MFfit: structure containing FOD calculations with fields:
%     volMF : multi-FOD fit parameters size = [nx,ny,nz,nb]
%     volb0 : average b=0 image        size = [nx,ny,nz]
%     volmask : mask volume            size = [nx,ny,nz]
%     volmask_dilated : dilated mask   size = [nx,ny,nz]
%     M : vox2ras matrix
%     volsz : size of input 4D vol          = [nx,ny,nz,nf]
%     qmat_orig : input q matrix       size = [nf,3]
%     qmat : unit vector q matrix      size = [nf,3]
%     bvals : input b values           size = [nf,1]
%     b0_scalefacts : scale factors calculated from b=0 images
%     i_b0 : indices of b=0 images
%     ADC_long : longitudinal ADC
%     ADC_free : free water ADC (e.g. CSF)
%     ADC_hindered : hindered water ADC (e.g. edema)
%     ADC_trans_vals : transverse ADC values
%     num_ADC_trans : number of transverse ADC values (size scales)
%     ADC_iso_vals : isotropic ADC values
%     num_ADC_iso : number of isotropic ADC values (size scales)
%     SH_order : spherical harmonic order
%     lambda : regularization parameter
%     icoverts : icosahedral sphere vertex coordinates
%     icostruct : icosahedral sphere surface struct
%     beta2ico : matrix mapping from FOD to ico3 surface
%     nb : number of parameters (betas)
%     A : RSI forward matrix           size = [nf,nb]
%     Ainv : RSI inverse matrix        size = [nb,nf]
%     B : tensor forward matrix        size = [nf,7]
%     Binv : tensor inverse matrix     size = [7,nf]
%
% Optional Parameters:
%   'mask_flag': [0|1] apply dilated brain mask to output volumes
%     {default = 1}
%
% Output:
%   MFmeas: structure containing multi-FOD measures
%     volmask : dilated brain mask         size = [nx,ny,nz]  
%     volsz   : size of original data           = [nx,ny,nz,nf]
%     ns      : number of size scales (num_ADC_trans)
%     volT    :  norm of all parameters    size = [nx,ny,nz]
%     volF0   : 0th order FOD              size = [nx,ny,nz,ns]
%     volN0   : volF0 normalized by volT   size = [nx,ny,nz,ns]
%     volF2   : norm of 2nd order FOD      size = [nx,ny,nz,ns]
%     volN2   : volF2 normalized by volT   size = [nx,ny,nz,ns]
%     volF4   : norm of 4th order FOD      size = [nx,ny,nz,ns]
%     volN4   : volF4 normalized by volT   size = [nx,ny,nz,ns]
%     volFD   : norm of directional FOD    size = [nx,ny,nz,ns]
%               (F2 and F4 combined)
%     volND   : volFD normalized by volT   size = [nx,ny,nz,ns]
%     volFT   : norm of total FOD          size = [nx,ny,nz,ns]
%               (F0, F2, and F4 combined)
%     volNT   : volFT normalized by volT   size = [nx,ny,nz,ns]
%     volI    : isotropic vol              size = [nx,ny,nz,ni]
%     volNI   : volI normalized by volT    size = [nx,ny,nz,ns]
%     volIr   : restricted isotropic vol   size = [nx,ny,nz]
%     volNIr  : volIr normalized by volT   size = [nx,ny,nz,ns]
%     volIh   : hindered isotropic vol     size = [nx,ny,nz]
%     volNIh  : volIh normalized by volT   size = [nx,ny,nz,ns]
%     volIf   : free isotropic vol         size = [nx,ny,nz]
%     volNIf  : volIf normalized by volT   size = [nx,ny,nz,ns]
%     volV0   : orientation of maximum FOD size = [nx,ny,nz,3]
%     volAU   : angular uncertainty        size = [nx,ny,nz]
%
% Created:  05/08/12 by Don Hagler
% Last Mod: 07/31/15 by Don Hagler
%

% based on code by Nate White

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
MFmeas = [];
parms = mmil_args2parms( varargin, { ...
  'mask_flag',true,[false true],...
... % constants
  'ind_F0',1,[],...
  'ind_F2',[2:6],[],...
  'ind_F4',[7:15],[],...
  'ind_FD',[],[],...
  'ind_FT',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MFmeas.volsz = MFfit.volsz;
if parms.mask_flag
  MFmeas.volmask = MFfit.volmask_dilated;
else
  MFmeas.volmask = 1.0*(MFfit.volb0~=0);
end;
MFmeas.ns = MFfit.num_ADC_trans;
MFmeas.ni = MFfit.num_ADC_iso;
volsz = MFmeas.volsz(1:3);
npf = (MFfit.SH_order+1)*(MFfit.SH_order+2)/2;

if isempty(parms.ind_FD)
  parms.ind_FD = [2:npf];
end;
if isempty(parms.ind_FT)
  parms.ind_FT = [1:npf];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all parameters combined
MFmeas.volT = sqrt(sum(MFfit.volMF.^2,4)).*MFmeas.volmask;

% spectrum of FOD-derived measures
MFmeas.volF0 = zeros([volsz,MFmeas.ns]);
MFmeas.volN0 = zeros([volsz,MFmeas.ns]);
MFmeas.volF2 = zeros([volsz,MFmeas.ns]);
MFmeas.volN2 = zeros([volsz,MFmeas.ns]);
if MFfit.SH_order>=4
  MFmeas.volF4 = zeros([volsz,MFmeas.ns]);
  MFmeas.volN4 = zeros([volsz,MFmeas.ns]);
else
  MFmeas.volF4 = [];
  MFmeas.volN4 = [];
end;
MFmeas.volFD = zeros([volsz,MFmeas.ns]);
MFmeas.volND = zeros([volsz,MFmeas.ns]);
MFmeas.volFT = zeros([volsz,MFmeas.ns]);
MFmeas.volNT = zeros([volsz,MFmeas.ns]);

for f=1:MFmeas.ns
  % F0: isotropic component
  ind = parms.ind_F0+(f-1)*npf;
  MFmeas.volF0(:,:,:,f) = ...
    squeeze(MFfit.volMF(:,:,:,ind)).*MFmeas.volmask;
  MFmeas.volN0(:,:,:,f) = ...
    (MFmeas.volF0(:,:,:,f) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
  % F2: lower-order directional component
  ind = parms.ind_F2+(f-1)*npf;
  MFmeas.volF2(:,:,:,f) = ...
    sqrt(sum(MFfit.volMF(:,:,:,ind).^2,4)).*MFmeas.volmask;
  MFmeas.volN2(:,:,:,f) = ...
    (MFmeas.volF2(:,:,:,f) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
  % F4: higher-order directional component
  if MFfit.SH_order>=4
    ind = parms.ind_F4+(f-1)*npf;
    MFmeas.volF4(:,:,:,f) = ...
      sqrt(sum(MFfit.volMF(:,:,:,ind).^2,4)).*MFmeas.volmask;
    MFmeas.volN4(:,:,:,f) = ...
      (MFmeas.volF4(:,:,:,f) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
  end;
  % FD: combined directional component
  ind_FD = intersect(parms.ind_FD,[1:npf]);
  ind  = ind_FD+(f-1)*npf;
  MFmeas.volFD(:,:,:,f) = ...
    sqrt(sum(MFfit.volMF(:,:,:,ind).^2,4)).*MFmeas.volmask;
  MFmeas.volND(:,:,:,f) = ...
    (MFmeas.volFD(:,:,:,f) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
  % FT: combined total
  ind_FT = intersect(parms.ind_FT,[1:npf]);
  ind  = ind_FT+(f-1)*npf;
  MFmeas.volFT(:,:,:,f) = ...
    sqrt(sum(MFfit.volMF(:,:,:,ind).^2,4)).*MFmeas.volmask;
  MFmeas.volNT(:,:,:,f) = ...
    (MFmeas.volFT(:,:,:,f) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
end;
ind = MFmeas.ns*npf+1;

% isotropic volume fractions
if MFmeas.ni>0
  MFmeas.volI = zeros([volsz,MFmeas.ni]);
  MFmeas.volNI = zeros([volsz,MFmeas.ni]);
  for f=1:MFmeas.ni
    MFmeas.volI(:,:,:,f) = ...
      squeeze(MFfit.volMF(:,:,:,ind)).*MFmeas.volmask;
    MFmeas.volNI(:,:,:,f) = ...
      (MFmeas.volI(:,:,:,f) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
    ind = ind + 1;
  end;
else
  MFmeas.volI = [];
  MFmeas.volNI = [];
end;

% specific isotropic volume fractions
iso_flags = {'iso_restricted_flag','iso_hindered_flag','iso_free_flag'};
iso_vols = {'volIr','volIh','volIf'};
iso_norm_vols = {'volNIr','volNIh','volNIf'};
for i=1:length(iso_flags)
  if MFfit.(iso_flags{i})
    MFmeas.(iso_vols{i}) = squeeze(MFfit.volMF(:,:,:,ind)).*MFmeas.volmask;
    MFmeas.(iso_norm_vols{i}) = ...
      (MFmeas.(iso_vols{i}) ./ max(MFmeas.volT,eps)).*MFmeas.volmask;
    ind = ind + 1;
  else
    MFmeas.(iso_vols{i}) = [];
    MFmeas.(iso_norm_vols{i}) = [];
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orientation measures from restricted FOD
MFmeas.volAU = zeros(volsz);
MFmeas.volV0 = zeros([volsz,3]);
if MFfit.num_ADC_trans
  vox = find(MFmeas.volmask);
  nvox = length(vox);
  for m=1:nvox
    [i j k] = ind2sub(volsz,vox(m));
    % maximum orientation of restricted FOD
    beta = squeeze(MFfit.volMF(i,j,k,:));
    fod = max(0,MFfit.beta2ico*beta(1:npf));
    [val ind] = max(fod); %% todo: implement more robust max finding procedure
    V0 = MFfit.icoverts(ind,:)';
    MFmeas.volV0(i,j,k,:) = V0;
    % angular uncertainty
    dprod = dot(MFfit.icoverts,repmat(V0',size(MFfit.icoverts,1),1),2);
    angular_uncertainty = sqrt(max(0,1-dprod.^2)).*fod;
    angular_uncertainty = asin(sum(angular_uncertainty)/max(eps,sum(fod)));
    MFmeas.volAU(i,j,k) = angular_uncertainty;
  end
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

