function FODfit = rsi_fit_FOD(vol,qmat,varargin)
%function FODfit = rsi_fit_FOD(vol,qmat,[options])
%
% Purpose: fit diffusion MRI data with fiber orientation density (FOD) function
% 
% Required Parameters:
%   vol: 4D volume containing multiple diffusion weighted volumes
%   qmat: matrix of diffusion direction vectors
%
% Optional Parameters:
%   'bvals': vector of b values
%     one for all, or one for each diffusion direction
%     {default = 1000}
%   'volmask': mask volume to restrict analysis to portion of image
%     if empty, will automatically generated from b=0 image
%     {default = []}
%   'M': 4x4 vox2ras matrix for vol (recommended for mask creation)
%     {default = identity matrix}
%   'lambda': regularization constant
%     {default = 0.1}
%   'ADC_long': longitudinal ADC
%     {default = 1e-3}
%   'ADC_trans': minimum transverse ADC
%     {default = 1e-4}
%   'SH_order': spherical harmonic order -- must be even
%     {default = 4}
%   'norm_flag': normalize data to b=0 image
%     {default = 0}
%
% Output:
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
% Created:  05/18/12 by Don Hagler
% Last Mod: 10/27/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
FODfit = [];
parms = mmil_args2parms( varargin, { ...
  'bvals',1000,[],...
  'volmask',[],[],...
  'M',eye(4),[],...
  'lambda',0.1,[],...
  'ADC_long',1e-3,[],...
  'ADC_trans',1e-4,[],...
  'SH_order',4,[2:2:10],...
  'norm_flag',false,[false true],...
... % parameters for brain mask creation
  'thresh',0.5,[],...
  'fill1_smooth1',15,[],...
  'fill1_thresh1',0.8,[],...
  'fill1_smooth2',0,[],...
  'fill1_thresh2',0,[],...
  'fill1_smooth3',0,[],...
  'fill1_erode_flag',true,[false true],...
  'clip_edges_flag',true,[false true],...
  'clip_edges_width',2,[1 10],...
  'fill2_smooth1',10,[],...
  'fill2_thresh1',0.1,[],...
  'fill2_smooth2',30,[],...
  'fill2_thresh2',0.7,[],...
  'fill2_smooth3',0,[],...
  'fill2_thresh3',0,[],...
  'fill2_erode_flag',false,[false true],...
  'binary_flag',true,[false true],...
... % parameters for brain  mask dilation
  'mask_smooth1',25,[],...
  'mask_thresh1',0.1,[],...
  'mask_smooth2',10,[],...
  'mask_thresh2',0.5,[],...
  'mask_smooth3',5,[],...
  'mask_thresh3',0.1,[],...  
... % misc parameters
  'nob0_flag',false,[false true],...
  'smf',10^-5,[10^-100,10^-1],...
... % parameters to be passed
  'mask_tags',{'log_flag' 'thresh' 'fill1_smooth1' 'fill1_thresh1'...
               'fill1_smooth2' 'fill1_thresh2' 'fill1_smooth3'...
               'fill1_erode_flag' 'clip_edges_flag' 'clip_edges_width'...
               'fill2_smooth1' 'fill2_thresh1' 'fill2_smooth2' 'fill2_thresh2'...
               'fill2_smooth3' 'fill2_thresh3' 'fill2_erode_flag' 'binary_flag'...
               'forceflag'},[],...
});

excl_tags = {'ADC_trans'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parms for running multi-shell FOD with single shell, no isotropic
parms.num_ADC_trans = 1;
parms.ADC_trans_min = parms.ADC_trans;
parms.ADC_trans_max = parms.ADC_trans;
parms.iso_restricted_flag = 0;
parms.iso_hindered_flag = 0;
parms.iso_free_flag = 0;

% run FOD calculations
tags = setdiff(fieldnames(parms),excl_tags);
args = mmil_parms2args(parms,tags);
FODfit = rsi_fit_MFOD(vol,qmat,args{:});

% change output
FODfit.volFOD = FODfit.volMF;
FODfit.ADC_trans = FODfit.ADC_trans_vals;

% remove inappropriate fields
FODfit = rmfield(FODfit,'volMF');
FODfit = rmfield(FODfit,'iso_restricted_flag');
FODfit = rmfield(FODfit,'iso_hindered_flag');
FODfit = rmfield(FODfit,'iso_free_flag');
FODfit = rmfield(FODfit,'ADC_free');
FODfit = rmfield(FODfit,'ADC_hindered');
FODfit = rmfield(FODfit,'ADC_trans_vals');
FODfit = rmfield(FODfit,'num_ADC_trans');


