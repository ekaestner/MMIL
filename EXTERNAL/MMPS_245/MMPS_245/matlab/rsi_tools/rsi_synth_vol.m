function vol = rsi_synth_vol(MFfit,varargin)
%function vol = rsi_synth_vol(MFfit,[options])
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
%   'rFOD_flag': [0|1] synthesize data using only most restricted FOD
%     {default = 0}
%   'rFOD_nscales': number of size scales to include in synthesized data
%     {default = 1}
%   'iso_flag': [0|1] whether to add isotropic signal to synthesized data
%     {default = 0}
%   'iso_ADC': ADC of isotropic signal added to synthesized data
%     to make b=0 images match; may be a scalar or volume
%     {default = 1e-3}
%   'iso_b0': b=0 intensity used to add isotropic signal to synthesized data
%     may be a scalar or volume
%     if not supplied, will use MFfit.volb0
%     {default = []}
%
% Created:  05/10/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse arguments
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms( varargin, { ...
  'rFOD_flag',false,[false true],...
  'rFOD_nscales',1,[1:100],...
  'iso_flag',false,[false true],...
  'iso_ADC',1e-3,[],...
  'iso_b0',[],[],...
});

nx = MFfit.volsz(1);
ny = MFfit.volsz(2);
nz = MFfit.volsz(3);
nf = MFfit.volsz(4);
nb = MFfit.nb;

% number of parameters for FOD at one size scale
npf = (MFfit.SH_order+1)*(MFfit.SH_order+2)/2;

% number of scales may not be greater MFfit.num_ADC_trans
nscales = min(parms.rFOD_nscales,MFfit.num_ADC_trans);

% indices of FOD parameters at 1 or more size scales
npfs = npf*nscales;
ind_rFOD = [1:npf*nscales];

% index of first b=0 image (usually first frame)
i_b0 = MFfit.i_b0(1);

if isempty(parms.iso_b0)
  parms.iso_b0 = MFfit.volb0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('%s: synthesizing volume from multi-FOD fit...\n',mfilename);
vol = zeros(MFfit.volsz,'single');
for z = 1:nz
  if ~parms.rFOD_flag
    betas = reshape(MFfit.volMF(:,:,z,:),[nx*ny,nb])';
    ysynth = max(0,MFfit.A*betas);
  else
    betas = reshape(MFfit.volMF(:,:,z,ind_rFOD),[nx*ny,npfs])';
    ysynth = max(0,MFfit.A(:,ind_rFOD)*betas);
  end;
  if MFfit.norm_flag
    yb0 = repmat(reshape(MFfit.volb0(:,:,z),[nx*ny,1])',[nf,1]);    
    ysynth = ysynth .* yb0;
  end;
  vol(:,:,z,:) = single(reshape(ysynth',[nx ny nf]));
end;

synth_b0 = vol(:,:,:,i_b0);
if numel(parms.iso_b0)==1
  synth_b0 = mean(synth_b0(MFfit.volmask>0))*MFfit.volmask_dilated;
  parms.iso_b0 = parms.iso_b0*MFfit.volmask_dilated;
end;
iso_b0 = parms.iso_b0 - synth_b0;
if numel(parms.iso_ADC)==1
  iso_ADC = parms.iso_ADC*MFfit.volmask_dilated;
else
  iso_ADC = parms.iso_ADC;
end;

if parms.iso_flag
  for f=1:nf
    vol(:,:,:,f) = vol(:,:,:,f) + iso_b0 .* exp(-iso_ADC*MFfit.bvals(f));
  end;
end;

