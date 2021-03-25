function MFfit = rsi_fit_MFOD(vol,qmat,varargin)
%function MFfit = rsi_fit_MFOD(vol,qmat,[options])
%
% Purpose: fit multi-shell diffusion MRI data with multi-compartment
%   fiber orientation density (FOD) function
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
%   'ADC_trans_min': minimum transverse ADC
%     {default = 0}
%   'ADC_trans_max': maximum transverse ADC
%     {default = 0.9e-3}
%   'num_ADC_trans': number of transverse ADC size scales
%     {default = 5}
%   'SH_order': spherical harmonic order -- must be even
%     {default = 4}
%   'iso_restricted_flag': [0|1] model isotropic diffusion of restricted water
%     {default = 1}
%   'iso_hindered_flag': [0|1] model isotropic diffusion of hindered water
%     {default = 1}
%   'iso_free_flag': [0|1] model isotropic diffusion of free water
%     {default = 1}
%   'ADC_hindered': ADC of isotropic hindered water (e.g. edema)
%     {default = 1.5e-3}
%   'ADC_free': apparent diffusion coefficient (ADC) of
%               isotropic free water (e.g. CSF)
%     {default = 3e-3}
%   'ADC_iso_min': minimum isotropic ADC
%     {default = 0}
%   'ADC_iso_max': maximum isotropic ADC
%     {default = 3e-3}
%   'num_ADC_iso': number of isotropic ADC size scales
%     {default = 0}
%   'ADC_iso_vals': vector of isotropic ADC values
%     if empty, will be set according to
%       ADC_iso_min, ADC_iso_max, and num_ADC_iso
%     {default = []}
%   'norm_flag': normalize data to b=0 image
%     {default = 0}
%   'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit
%     {default = 0}
%   'b0_thresh': threshold used for considering a b-value to be 0
%     {default = 10}
%   'scalefacts_flag': [0|1] calculate scaling factors from b=0 images
%     and apply them to all subsequent frames (for multiple acquisitions)
%     {default = 0}
%
% Output:
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
%     iso_free_flag
%     iso_hindered_flag
%     iso_restricted_flag
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
% Created:  05/08/12 by Don Hagler based on code by Nate White
% Prev Mod: 10/27/15 by Don Hagler
% Last Mod: 11/03/17 by Don Hagler
%

%% todo: remove iso_restricted, iso_hindered, and iso_free (use ADC_iso_vals)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MFfit = [];

% check input arguments
if ~mmil_check_nargs(nargin, 2), return; end;
parms = check_input(varargin);

% initialize output struct
[MFfit,parms] = init_MFfit(vol,qmat,parms);

% linear RSI fit
MFfit = fit_rsi(vol,MFfit,parms);

% nonlinear RSI fit
if parms.nonlin_flag
  MFfit = nlfit_rsi(vol,MFfit,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'bvals',1000,[],...
    'volmask',[],[],...
    'M',eye(4),[],...
    'lambda',0.1,[],...
    'ADC_trans_min',0,[],...
    'ADC_trans_max',0.9e-3,[],...
    'num_ADC_trans',5,[],...
    'SH_order',4,[2:2:10],...
    'iso_free_flag',true,[false true],...
    'iso_hindered_flag',true,[false true],...
    'iso_restricted_flag',true,[false true],...
    'ADC_free',3e-3,[],...
    'ADC_hindered',1.5e-3,[],...
    'ADC_long',1e-3,[],...
    'ADC_iso_min',0,[],...
    'ADC_iso_max',3e-3,[],...
    'num_ADC_iso',0,[],...
    'ADC_iso_vals',[],[],...
    'norm_flag',false,[false true],...
    'nonlin_flag',false,[false true],...
    'b0_thresh',10,[0,500],...
    'scalefacts_flag',false,[false true],...
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
    'nlfit_method','fmincon',{'lsqnonlin','fmincon'},...
    'nlfit_display','none',[],...
    'nlfit_F0_bounds',[0,100],[-Inf,Inf],...
    'nlfit_FX_bounds',[-15,15],[-Inf,Inf],...
    'nlfit_FX_flag',false,[false true],...
  ... % parameters to be passed
    'mask_tags',{'log_flag' 'thresh' 'fill1_smooth1' 'fill1_thresh1'...
                 'fill1_smooth2' 'fill1_thresh2' 'fill1_smooth3'...
                 'fill1_erode_flag' 'clip_edges_flag' 'clip_edges_width'...
                 'fill2_smooth1' 'fill2_thresh1' 'fill2_smooth2' 'fill2_thresh2'...
                 'fill2_smooth3' 'fill2_thresh3' 'fill2_erode_flag' 'binary_flag'...
                 'forceflag'},[],...
  });
  
  % set ADC iso vals
  if isempty(parms.ADC_iso_vals) && parms.num_ADC_iso>0
    parms.ADC_iso_vals = ...
      linspace(parms.ADC_iso_min,parms.ADC_iso_max,parms.num_ADC_iso);
  else
    parms.num_ADC_iso = length(parms.ADC_iso_vals);
  end;
  
  % initialize nonlinear fit parameters
  if parms.nonlin_flag
    switch parms.nlfit_method
      case 'lsqnonlin'
        opts = {'Display',parms.nlfit_display,...
                'MaxIter',100,...
                'TolFun',1e-10};
        try % for R2009b
          parms.nlfit_options = optimset(opts{:},...
            'Algorithm','levenberg-marquardt');
        catch % for R2007a
          parms.nlfit_options = optimset(opts{:},...
            'LargeScale','off');
        end;
      case 'fmincon'
        opts = {'Display',parms.nlfit_display,...
                'MaxFunEvals',Inf,...
                'MaxIter',100,...
                'TolFun',1e-6,...
                'TolX',1e-7};
        try % for R2009b
          parms.nlfit_options = optimset(opts{:},...
            'Algorithm','active-set');
        catch % for R2007a
          parms.nlfit_options = optimset(opts{:},...
            'LargeScale','off');
        end;
    end;

    % set ind_L0 to identify volume fraction parameters
    parms.FOD_nb = (parms.SH_order^2 + parms.SH_order + 2)/2 + parms.SH_order;
    nb = parms.FOD_nb*parms.num_ADC_trans;
    parms.ind_L0 = [1:parms.FOD_nb:nb];
    if parms.num_ADC_iso>0
      parms.ind_L0 = cat(2,parms.ind_L0,[nb+1:nb+parms.num_ADC_iso]);
      nb = nb+parms.num_ADC_iso;
    end;
    % isotropic restricted
    if parms.iso_restricted_flag
      nb = nb+1;
      parms.ind_L0 = cat(2,parms.ind_L0,nb);
    end;
    % isotropic hindered water fractions
    if parms.iso_hindered_flag
      nb = nb+1;
      parms.ind_L0 = cat(2,parms.ind_L0,nb);
    end;
    % isotropic free and hindered water fractions
    if parms.iso_free_flag
      nb = nb+1;
      parms.ind_L0 = cat(2,parms.ind_L0,nb);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MFfit,parms] = init_MFfit(vol,qmat,parms)
  MFfit = [];
  MFfit.volsz = size(vol);
  [parms.nx,parms.ny,parms.nz,parms.nf] = size(vol);

  % extend bvals to nf elements if only 1 supplied
  if numel(parms.bvals)==1
    MFfit.bvals = parms.bvals*ones(parms.nf,1);
  elseif numel(parms.bvals)==parms.nf
    MFfit.bvals = reshape(parms.bvals,[parms.nf,1]);
  else
    error('number of bvals must match frames in vol');
  end;

  % check qmat is correct size
  if size(qmat,1)~=parms.nf
    error('number of directions in qmat must match frames in vol (%d)',parms.nf);
  end;
  if size(qmat,2)~=3
    error('qmat must have 3 columns');
  end;

  % ensure that qmat contains unit vectors
  qlength = max(sqrt(sum(qmat.^2,2)),eps);
  MFfit.qmat_orig = qmat;
  MFfit.qmat = qmat ./ repmat(qlength,1,3);

  % set ADC values
  MFfit.ADC_long = parms.ADC_long;
  MFfit.ADC_trans_vals =...
    linspace(parms.ADC_trans_min,parms.ADC_trans_max,parms.num_ADC_trans);
  MFfit.num_ADC_trans = parms.num_ADC_trans;

  % set ADC_iso values
  if parms.iso_free_flag
    MFfit.ADC_free = parms.ADC_free;
  else
    MFfit.ADC_free = 0;
  end;  
  if parms.iso_hindered_flag
    MFfit.ADC_hindered = parms.ADC_hindered;
  else
    MFfit.ADC_hindered = 0;
  end;  
  MFfit.ADC_iso_vals = parms.ADC_iso_vals;
  MFfit.num_ADC_iso = parms.num_ADC_iso;
  
  % identify b=0 frames
  qlength = sqrt(sum(MFfit.qmat.^2,2));
  i_b0 = find(qlength<parms.smf | MFfit.bvals<=parms.b0_thresh);
  i_bx = setdiff([1:parms.nf],i_b0);
  parms.nb0 = length(i_b0);
  ndirs = length(i_bx);
  volb0 = vol(:,:,:,i_b0);
  MFfit.bvals(i_b0) = 0;
  MFfit.i_b0 = i_b0;

  % calculate scaling factors and apply to each frame
  if parms.scalefacts_flag && length(i_b0)>1
    MFfit.b0_scalefacts = dti_calc_b0_scale(volb0);
    fprintf('%s: scaling factors from %d b=0 images: %s\n',...
      mfilename,parms.nb0,sprintf('%0.2f  ',MFfit.b0_scalefacts));
    vol = dti_apply_b0_scale(vol,i_b0,MFfit.b0_scalefacts);
  else
    MFfit.b0_scalefacts = 1;
  end;

  % average b=0 images
  if parms.nb0>0
    MFfit.volb0 = mean(volb0,4);
    clear volb0;
  else
    MFfit.volb0 = [];
    if parms.norm_flag
      fprintf('%s: WARNING: setting norm_flag to 0 because no b=0 images were found\n',...
        mfilename);
      parms.norm_flag = 0;
    end;
  end;

  % create brain mask
  if isempty(parms.volmask)
    if ~isempty(MFfit.volb0)
      % create approximate brain mask from average b=0 image
      MFfit.volmask = create_brain_mask(MFfit.volb0,parms,1);
    else
      % exclude voxels with zeros
      MFfit.volmask = ones(parms.nx,parms.ny,parms.nz);
      tmp = min(vol,[],4);
      MFfit.volmask(tmp==0) = 0;
    end;
  else
    MFfit.volmask = parms.volmask;
  end;

  MFfit.i_fit = [1:MFfit.volsz(4)];

  % exclude b=0 images
  if parms.nob0_flag
    i_keep = find(MFfit.bvals > parms.b0_thresh);
    parms.nf = length(i_keep);
    vol = vol(:,:,:,i_keep);
    MFfit.volsz = size(vol);
    MFfit.bvals = MFfit.bvals(i_keep);
    MFfit.qmat = MFfit.qmat(i_keep,:);
    MFfit.i_fit = MFfit.i_fit(i_keep);
    parms.nb0 = 0;
  end;

  fprintf('%s: including %d directions in multi-FOD fit...\n',...
  mfilename,length(MFfit.i_fit));
  
  MFfit.i_b0 = find(MFfit.bvals <= parms.b0_thresh);
  MFfit.i_bx = find(MFfit.bvals > parms.b0_thresh);
  MFfit.nb0 = length(MFfit.i_b0);

  % store info about data volume
  MFfit.M = parms.M;
  MFfit.volsz = [parms.nx,parms.ny,parms.nz,parms.nf];

  % store additional values
  MFfit.SH_order = parms.SH_order;
  MFfit.lambda = parms.lambda;
  MFfit.norm_flag = parms.norm_flag;
  MFfit.iso_restricted_flag = parms.iso_restricted_flag;
  MFfit.iso_hindered_flag = parms.iso_hindered_flag;
  MFfit.iso_free_flag = parms.iso_free_flag;

  % create forward matrix for fitting tensor
  % construct tensor forward matrix
  [MFfit.B,MFfit.Binv] = create_tensor_matrix(MFfit);

  % construct RSI multi-FOD forward matrix
  [MFfit.A,MFfit.Ainv,MFfit.icoverts,MFfit.icostruct,MFfit.beta2ico]=...
    create_rsi_matrix(MFfit);
  MFfit.nb = size(MFfit.A,2);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,Binv] = create_tensor_matrix(MFfit)
  ndirs = size(MFfit.qmat,1);
  B = zeros(ndirs,7);
  B(:,7) = 1; % S0
  for i = 1:ndirs
    outerprod = -MFfit.bvals(i).*MFfit.qmat(i,:)'*MFfit.qmat(i,:);
    B(i,1:3) = diag(outerprod)';
    B(i,4) = 2*outerprod(1,2);
    B(i,5) = 2*outerprod(1,3);
    B(i,6) = 2*outerprod(2,3);
  end
  Binv = pinv(B);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,Ainv,V,VV,F] = create_rsi_matrix(MFfit)
  VV = icosaSphere(3);
  V = VV.vertices;
  Q = MFfit.qmat.*repmat(sqrt(MFfit.bvals),1,3);
  F = rsi_SH_matrix(V,MFfit.SH_order);
  F0 = rsi_SH_matrix(V,0);
  A = [];

  % series of FODs for varying size scales
  for t=1:MFfit.num_ADC_trans
    R = rsi_FOD_matrix(Q,V,MFfit.ADC_long,MFfit.ADC_trans_vals(t));
    A = [A R*F];
  end

  % series of isotropic FODs for varying size scales
  for t=1:MFfit.num_ADC_iso
    R = rsi_FOD_matrix(Q,V,MFfit.ADC_iso_vals(t),MFfit.ADC_iso_vals(t));
    A = [A R*F0];
  end

  % isotropic restricted
  if MFfit.iso_restricted_flag
    R0 = rsi_FOD_matrix(Q,V,0,0);
    A = [A R0*F0];
  end;

  % isotropic hindered water fractions
  if MFfit.iso_hindered_flag
    R_hindered = rsi_FOD_matrix(Q,V,MFfit.ADC_hindered,MFfit.ADC_hindered);
    A = [A R_hindered*F0];
  end;

  % isotropic free and hindered water fractions
  if MFfit.iso_free_flag
    R_free = rsi_FOD_matrix(Q,V,MFfit.ADC_free,MFfit.ADC_free);
    A = [A R_free*F0];
  end;

  % compute regularized inverse
  AtA = A'*A;
  Ainv = inv(AtA+MFfit.lambda*mean(diag(AtA))*eye(size(AtA)))*A';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lb,ub] = set_bounds(MFfit,parms)
  % create lower and upper bounds vectors for nonlinear fitting
  lb = [];
  ub = [];
  % use all parameters or volume fractions only
  if parms.nlfit_FX_flag
    lb = parms.nlfit_FX_bounds(1)*ones(MFfit.nb,1);
    ub = parms.nlfit_FX_bounds(2)*ones(MFfit.nb,1);
    lb(parms.ind_L0) = parms.nlfit_FX_bounds(1);
    ub(parms.ind_L0) = parms.nlfit_FX_bounds(2);
  else
    lb = parms.nlfit_F0_bounds(1)*ones(length(parms.ind_L0),1);
    ub = parms.nlfit_F0_bounds(2)*ones(length(parms.ind_L0),1);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MFfit = fit_rsi(vol,MFfit,parms)
  % initialize multi-FOD parameter volumes
  MFfit.volMF = zeros(parms.nx,parms.ny,parms.nz,MFfit.nb);
  % loop over slices
  for z = 1:parms.nz
    y = reshape(vol(:,:,z,:),[parms.nx*parms.ny,parms.nf])';
    if parms.norm_flag % normalize data to b=0
      y0 = repmat(reshape(MFfit.volb0(:,:,z),[parms.nx*parms.ny,1])',[parms.nf 1]);
      y = y./y0;
    end;
    betas = MFfit.Ainv*y;
    betas = reshape(betas',[parms.nx parms.ny MFfit.nb]);
    betas(isnan(betas)) = 0;
    MFfit.volMF(:,:,z,:) = betas;
  end
  % dilate brain mask
  MFfit.volmask_dilated = dilate_brain_mask(MFfit.volmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MFfit = nlfit_rsi(vol,MFfit,parms)
  % set lower and upper bounds
  [lb,ub] = set_bounds(MFfit,parms);
  % loop over all voxels inside mask
  ind_mask = find(MFfit.volmask_dilated);
  for m=1:length(ind_mask)
    [i,j,k]=ind2sub([parms.nx,parms.ny,parms.nz],ind_mask(m));
    % get data vector for this voxel
    data_vec = double(squeeze(vol(i,j,k,:)));
    if parms.norm_flag
      y0 = double(MFfit.volb0(i,j,k));
    else
      y0 = 1;
    end;
    % linear fit as initial estimate
    beta_lin = squeeze(MFfit.volMF(i,j,k,:));
    % use all parameters or volume fractions only
    if parms.nlfit_FX_flag
      beta_cur = beta_lin;
    else
      beta_cur = beta_lin(parms.ind_L0);
    end;
    % nonlinear fit
    switch parms.nlfit_method
      case 'lsqnonlin'
        beta_nonlin = lsqnonlin(@(beta) costfunc(beta,...
                        beta_lin,MFfit.A,data_vec,y0,parms),...
                        beta_cur,lb,ub,parms.nlfit_options);
      case 'fmincon'
        beta_nonlin = fmincon(@(beta) sum(costfunc(beta,...
                        beta_lin,MFfit.A,data_vec,y0,parms).^2),...
                        beta_cur,[],[],[],[],lb,ub,[],parms.nlfit_options);
    end;
    % use volume fraction parameters only
    if ~parms.nlfit_FX_flag
      beta_cur = beta_nonlin;
      beta_nonlin = beta_lin;
      beta_nonlin(parms.ind_L0) = beta_cur;
    end;    
    MFfit.volMF(i,j,k,:) = beta_nonlin;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = costfunc(beta,beta_full,A,y,y0,parms)
  cost = [];
  % use all parameters or volume fractions only
  if parms.nlfit_FX_flag
    beta_full = beta;
  else
    beta_full(parms.ind_L0) = beta;
  end;
  yh = y0*A*beta_full;
  cost = y-yh;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function volmask = create_brain_mask(vol,parms,log_flag)
  parms.log_flag = log_flag;
  vol = ctx_mgh2ctx(vol,parms.M);
  args = mmil_parms2args(parms,parms.mask_tags);
  vol = mmil_quick_brainmask(vol,args{:});
  volmask = vol.imgs;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function volmask = dilate_brain_mask(volmask,parms)
  if all(volmask==1), return; end;
  volmask = ctx_mgh2ctx(volmask,parms.M);
  volmask = mmil_dilate_mask(volmask,...
    'smooth1',parms.mask_smooth1,'thresh1',parms.mask_thresh1,...
    'smooth2',parms.mask_smooth2,'thresh2',parms.mask_thresh2,...
    'smooth3',parms.mask_smooth3,'thresh3',parms.mask_thresh3);
  volmask = 1.0*(volmask.imgs>0);  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
