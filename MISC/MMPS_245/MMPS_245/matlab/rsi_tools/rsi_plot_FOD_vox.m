function rsi_plot_FOD_vox(FODfit,xyz,varargin)
%function rsi_plot_FOD_vox(FODfit,xyz,[options])
%
% Purpose: plot FOD for a single voxel
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
%   xyz: vector of voxel indices [row,column,slice]
%
%  NOTE: may supply MFfit structure instead of FODfit
%
% Optional Parameters:
%   'scale': index of size scale (must be less than or equal to num_ADC_trans)
%     {default = 1}
%   'clim': vector of lower and upper bounds for color scale used by imagesc
%     {default = []}
%
% see also: rsi_fit_FOD, rsi_fit_MFOD, rsi_plot_FOD
%
% Created:  02/28/14 by Don Hagler
% Last Mod: 02/28/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms( varargin, { ...
  'scale',1,[],...
  'clim',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(FODfit,'volMF')
  volFOD = FODfit.volMF;
else
  volFOD = FODfit.volFOD;
end;

volsz = FODfit.volsz(1:3);

if any(xyz < 1) | any(volsz-xyz < 0)
  error('xyz out of bounds');
end;

if parms.scale>FODfit.num_ADC_trans
  error('scale must be <= %d',FODfit.num_ADC_trans);
end;

npf = (FODfit.SH_order+1)*(FODfit.SH_order+2)/2;

[TH,PHI,RHO] = ...
  cart2sph(FODfit.icoverts(:,1),FODfit.icoverts(:,2),FODfit.icoverts(:,3));
FOD.faces = FODfit.icostruct.faces;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = squeeze(volFOD(xyz(1),xyz(2),xyz(3),:));
ind_scale = (parms.scale - 1)*npf + [1:npf];
icoF = max(0,FODfit.beta2ico*beta(ind_scale));
[val ind] = max(icoF);
V0 = FODfit.icoverts(ind,:)';
[X,Y,Z] = sph2cart(TH,PHI,0.5*icoF/max(icoF));
FOD.vertices = [Y,X,Z];
p = patch(FOD);
tcolor = repmat(max(icoF/max(icoF),eps),1,3).*abs(FODfit.icoverts);
set(p,'FaceVertexCData',tcolor,'LineStyle','none','FaceColor','interp','Facelighting','phong');

axis ij; axis image;
axis off;
