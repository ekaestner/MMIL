function FODmeas = rsi_calc_FODmeas(FODfit,varargin)
%function FODmeas = rsi_calc_FODmeas(FODfit,[options])
%
% Purpose: calculate various measures based on
%   fiber orientation density (FOD) function fit of diffusion MRI data 
%
% Required Parameters:
%   FODfit: structure containing FOD calculations with fields:
%     volFOD : FOD fit parameters      size = [nx,ny,nz,nb]
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
%     ADC_trans : transverse ADC
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
%   FODmeas: structure containing FOD measures
%     volmask : dilated brain mask         size = [nx,ny,nz]  
%     volsz   : size of original data           = [nx,ny,nz,nf]
%     volT    : norm of all parameters     size = [nx,ny,nz]
%     volF0   : 0th order FOD              size = [nx,ny,nz]
%     volN0   : volF0 normalized by volT   size = [nx,ny,nz]
%     volF2   : norm of 2nd order FOD      size = [nx,ny,nz]
%     volN2   : volF2 normalized by volT   size = [nx,ny,nz]
%     volF4   : norm of 4th order FOD      size = [nx,ny,nz]
%     volN4   : volF4 normalized by volT   size = [nx,ny,nz]
%     volFD   : norm of directional FOD    size = [nx,ny,nz,ns]
%               (F2 and F4 combined)
%     volND   : volFD normalized by volT   size = [nx,ny,nz,ns]
%     volFT   : norm of total FOD          size = [nx,ny,nz,ns]
%               (F0, F2, and F4 combined)
%     volNT   : volFT normalized by volT   size = [nx,ny,nz,ns]
%     volV0   : orientation of maximum FOD size = [nx,ny,nz,3]
%     volAU   : angular uncertainty        size = [nx,ny,nz]
%
% Created:  05/18/12 by Don Hagler
% Last Mod: 07/31/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
FODmeas = []; DTmeas = [];
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

FODmeas.volsz = FODfit.volsz;
if parms.mask_flag
  FODmeas.volmask = FODfit.volmask_dilated;
else
  FODmeas.volmask = 1.0*(FODfit.volb0~=0);
end;
volsz = FODmeas.volsz(1:3);
npf = (FODfit.SH_order+1)*(FODfit.SH_order+2)/2;

if isempty(parms.ind_FD)
  parms.ind_FD = [2:npf];
end;
if isempty(parms.ind_FT)
  parms.ind_FT = [1:npf];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all parameters combined
FODmeas.volT = sqrt(sum(FODfit.volFOD.^2,4)).*FODmeas.volmask;

% spectrum of FODs
FODmeas.volF0 = zeros(volsz);
FODmeas.volN0 = zeros(volsz);
FODmeas.volF2 = zeros(volsz);
FODmeas.volN2 = zeros(volsz);
if FODfit.SH_order>=4
  FODmeas.volF4 = zeros(volsz);
  FODmeas.volN4 = zeros(volsz);
else
  FODmeas.volF4 = [];
  FODmeas.volN4 = [];
end;
FODmeas.volFD = zeros(volsz);
FODmeas.volND = zeros(volsz);
FODmeas.volFT = zeros(volsz);
FODmeas.volNT = zeros(volsz);

% F0: isotropic component
FODmeas.volF0 = ...
  squeeze(FODfit.volFOD(:,:,:,parms.ind_F0)).*FODmeas.volmask;
FODmeas.volN0 = ...
  (FODmeas.volF0 ./ max(FODmeas.volT,eps)).*FODmeas.volmask;
% F2: lower-order directional component
FODmeas.volF2 = ...
  sqrt(sum(FODfit.volFOD(:,:,:,parms.ind_F2).^2,4).*FODmeas.volmask);
FODmeas.volN2 = ...
  (FODmeas.volF2 ./ max(FODmeas.volT,eps)).*FODmeas.volmask;
% F4: higher-order directional component
if FODfit.SH_order>=4
  FODmeas.volF4 = ...
    sqrt(sum(FODfit.volFOD(:,:,:,parms.ind_F4).^2,4).*FODmeas.volmask);
  FODmeas.volN4 = ...
    (FODmeas.volF4 ./ max(FODmeas.volT,eps)).*FODmeas.volmask;
end;

% FD: combined directional component
ind_FD = intersect(parms.ind_FD,[1:npf]);
FODmeas.volFD = ...
  sqrt(sum(FODfit.volFOD(:,:,:,ind_FD).^2,4)).*FODmeas.volmask;
FODmeas.volND = ...
  (FODmeas.volFD ./ max(FODmeas.volT,eps)).*FODmeas.volmask;
% FT: combined total
ind_FT = intersect(parms.ind_FT,[1:npf]);
FODmeas.volFT = ...
  sqrt(sum(FODfit.volFOD(:,:,:,ind_FT).^2,4)).*FODmeas.volmask;
FODmeas.volNT = ...
  (FODmeas.volFT ./ max(FODmeas.volT,eps)).*FODmeas.volmask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orientation measures from FOD
FODmeas.volAU = zeros(volsz);
FODmeas.volV0 = zeros([volsz,3]);
vox = find(FODmeas.volmask);
nvox = length(vox);
for m=1:nvox
  [i j k] = ind2sub(volsz,vox(m));
  % maximum orientation of FOD
  beta = squeeze(FODfit.volFOD(i,j,k,:));
  fod = max(0,FODfit.beta2ico*beta);
  [val ind] = max(fod); %% todo: implement more robust max finding procedure
  V0 = FODfit.icoverts(ind,:)';
  FODmeas.volV0(i,j,k,:) = V0;
  % angular uncertainty
  dprod = dot(FODfit.icoverts,repmat(V0',size(FODfit.icoverts,1),1),2);
  angular_uncertainty = sqrt(max(0,1-dprod.^2)).*fod;
  angular_uncertainty = asin(sum(angular_uncertainty)/max(eps,sum(fod)));
  FODmeas.volAU(i,j,k) = angular_uncertainty;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

