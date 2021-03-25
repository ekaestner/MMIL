function vol = rsi_synth_iso_vol(MIfit)
%function vol = rsi_synth_iso_vol(MIfit)
%
% Required Parameters:
%   MIfit: structure containing FOD calculations with fields:
%     volMI : multi-iso fit parameters size = [nx,ny,nz,nb]
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
%     ADC_iso_vals : isotropic ADC values
%     num_ADC_iso : number of isotropic ADC values (size scales)
%     lambda : regularization parameter
%     nb : number of parameters (betas)
%     A : RSI forward matrix           size = [nf,nb]
%     Ainv : RSI inverse matrix        size = [nb,nf]
%
% Created:  07/29/13 by Don Hagler
% Last Mod: 07/29/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse arguments
if ~mmil_check_nargs(nargin, 1), return; end;

nx = MIfit.volsz(1);
ny = MIfit.volsz(2);
nz = MIfit.volsz(3);
nf = MIfit.volsz(4);
nb = MIfit.nb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('%s: synthesizing volume from multi-iso fit...\n',mfilename);
vol = zeros(MIfit.volsz,'single');
for z = 1:nz
  betas = reshape(MIfit.volMI(:,:,z,:),[nx*ny,nb])';
  ysynth = max(0,MIfit.A*betas);
  if MIfit.norm_flag
    yb0 = repmat(reshape(MIfit.volb0(:,:,z),[nx*ny,1])',[nf,1]);    
    ysynth = ysynth .* yb0;
  end;
  vol(:,:,z,:) = single(reshape(ysynth',[nx ny nf]));
end;

