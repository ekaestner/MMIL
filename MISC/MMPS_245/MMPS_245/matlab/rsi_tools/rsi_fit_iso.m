function MIfit = rsi_fit_iso(vol,qmat,varargin)
%function MIfit = rsi_fit_iso(vol,qmat,[options])
%
% Purpose: fit multi-shell diffusion MRI data
%   with multi-compartment isotropic difusion function
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
%   'ADC_iso_min': minimum isotropic ADC
%     {default = 0}
%   'ADC_iso_max': maximum isotropic ADC
%     {default = 3e-3}
%   'num_ADC_iso': number of isotropic ADC size scales
%     {default = 10}
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
%   MIfit: structure containing iso calculations with fields:
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
%     A : iso forward matrix           size = [nf,nb]
%     Ainv : iso inverse matrix        size = [nb,nf]
%
% Created:  07/29/13 by Don Hagler
% Prev Mod: 10/27/15 by Don Hagler
% Last Mod: 11/03/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIfit = [];

% check input arguments
if ~mmil_check_nargs(nargin, 2), return; end;
parms = check_input(varargin);

% initialize output struct
[MIfit,parms] = init_MIfit(vol,qmat,parms);

% linear iso fit
MIfit = fit_iso(vol,MIfit,parms);

% nonlinear iso fit
if parms.nonlin_flag
  MIfit = nlfit_iso(vol,MIfit,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'bvals',1000,[],...
    'volmask',[],[],...
    'M',eye(4),[],...
    'lambda',0.1,[],...
    'ADC_iso_min',0,[],...
    'ADC_iso_max',3e-3,[],...
    'num_ADC_iso',10,[],...
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
    'nlfit_bounds',[0,100],[-Inf,Inf],...
  ... % parameters to be passed
    'mask_tags',{'log_flag' 'thresh' 'fill1_smooth1' 'fill1_thresh1'...
                 'fill1_smooth2' 'fill1_thresh2' 'fill1_smooth3'...
                 'fill1_erode_flag' 'clip_edges_flag' 'clip_edges_width'...
                 'fill2_smooth1' 'fill2_thresh1' 'fill2_smooth2' 'fill2_thresh2'...
                 'fill2_smooth3' 'fill2_thresh3' 'fill2_erode_flag' 'binary_flag'...
                 'forceflag'},[],...
  });
  
  % set ADC iso vals
  if isempty(parms.ADC_iso_vals)
    parms.ADC_iso_vals = ...
      linspace(parms.ADC_iso_min,parms.ADC_iso_max,parms.num_ADC_iso);
  else
    parms.num_ADC_iso = length(parms.ADC_iso_vals);
  end;
  
  % initialize nonlinear fit output
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
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MIfit,parms] = init_MIfit(vol,qmat,parms)
  MIfit = [];
  MIfit.volsz = size(vol);
  [parms.nx,parms.ny,parms.nz,parms.nf] = size(vol);

  % extend bvals to nf elements if only 1 supplied
  if numel(parms.bvals)==1
    MIfit.bvals = parms.bvals*ones(parms.nf,1);
  elseif numel(parms.bvals)==parms.nf
    MIfit.bvals = reshape(parms.bvals,[parms.nf,1]);
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
  MIfit.qmat_orig = qmat;
  MIfit.qmat = qmat ./ repmat(qlength,1,3);

  % set ADC_iso values
  MIfit.ADC_iso_vals = parms.ADC_iso_vals;
  MIfit.num_ADC_iso = parms.num_ADC_iso;
  
  % identify b=0 frames
  qlength = sqrt(sum(MIfit.qmat.^2,2));
  i_b0 = find(qlength<parms.smf | MIfit.bvals<=parms.b0_thresh);
  i_bx = setdiff([1:parms.nf],i_b0);
  parms.nb0 = length(i_b0);
  ndirs = length(i_bx);
  volb0 = vol(:,:,:,i_b0);
  MIfit.bvals(i_b0) = 0;
  MIfit.i_b0 = i_b0;

  % calculate scaling factors and apply to each frame
  if parms.scalefacts_flag && length(i_b0)>1
    MIfit.b0_scalefacts = dti_calc_b0_scale(volb0);
    fprintf('%s: scaling factors from %d b=0 images: %s\n',...
      mfilename,parms.nb0,sprintf('%0.2f  ',MIfit.b0_scalefacts));
    vol = dti_apply_b0_scale(vol,i_b0,MIfit.b0_scalefacts);
  else
    MIfit.b0_scalefacts = 1;
  end;

  % average b=0 images
  if parms.nb0>0
    MIfit.volb0 = mean(volb0,4);
    clear volb0;
  else
    MIfit.volb0 = [];
    if parms.norm_flag
      fprintf('%s: WARNING: setting norm_flag to 0 because no b=0 images were found\n',...
        mfilename);
      parms.norm_flag = 0;
    end;
  end;

  % create brain mask
  if isempty(parms.volmask)
    if ~isempty(MIfit.volb0)
      % create approximate brain mask from average b=0 image
      MIfit.volmask = create_brain_mask(MIfit.volb0,parms,1);
    else
      % exclude voxels with zeros
      MIfit.volmask = ones(parms.nx,parms.ny,parms.nz);
      tmp = min(vol,[],4);
      MIfit.volmask(tmp==0) = 0;
    end;
  else
    MIfit.volmask = parms.volmask;
  end;

  MIfit.i_fit = [1:MIfit.volsz(4)];

  % exclude b=0 images
  if parms.nob0_flag
    i_keep = find(MIfit.bvals > parms.b0_thresh);
    parms.nf = length(i_keep);
    vol = vol(:,:,:,i_keep);
    MIfit.volsz = size(vol);
    MIfit.bvals = MIfit.bvals(i_keep);
    MIfit.qmat = MIfit.qmat(i_keep,:);
    MIfit.i_fit = MIfit.i_fit(i_keep);
    parms.nb0 = 0;
  end;

  fprintf('%s: including %d directions in multi-iso fit...\n',...
  mfilename,length(MIfit.i_fit));
  
  MIfit.i_b0 = find(MIfit.bvals <= parms.b0_thresh);
  MIfit.i_bx = find(MIfit.bvals > parms.b0_thresh);
  MIfit.nb0 = length(MIfit.i_b0);

  % store info about data volume
  MIfit.M = parms.M;
  MIfit.volsz = [parms.nx,parms.ny,parms.nz,parms.nf];

  % store additional values
  MIfit.lambda = parms.lambda;
  MIfit.norm_flag = parms.norm_flag;

  % construct iso multi-FOD forward matrix
  [MIfit.A,MIfit.Ainv] = create_iso_matrix(MIfit);
  MIfit.nb = size(MIfit.A,2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,Ainv] = create_iso_matrix(MIfit)
  VV = icosaSphere(3);
  V = VV.vertices;
  Q = MIfit.qmat.*repmat(sqrt(MIfit.bvals),1,3);
  F0 = rsi_SH_matrix(V,0);
  A = [];

  % series of isotropic FODs for varying size scales
  for t=1:MIfit.num_ADC_iso
    R = rsi_FOD_matrix(Q,V,MIfit.ADC_iso_vals(t),MIfit.ADC_iso_vals(t));
    A = [A R*F0];
  end

  % compute regularized inverse
  AtA = A'*A;
  Ainv = inv(AtA+MIfit.lambda*mean(diag(AtA))*eye(size(AtA)))*A';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lb,ub] = set_bounds(MIfit,parms)
  % create lower and upper bounds vectors for nonlinear fitting
  lb = [];
  ub = [];
  lb = parms.nlfit_bounds(1)*ones(MIfit.nb,1);
  ub = parms.nlfit_bounds(2)*ones(MIfit.nb,1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MIfit = fit_iso(vol,MIfit,parms)
  % initialize multi-FOD parameter volumes
  MIfit.volMI = zeros(parms.nx,parms.ny,parms.nz,MIfit.nb);
  % loop over slices
  for z = 1:parms.nz
    y = reshape(vol(:,:,z,:),[parms.nx*parms.ny,parms.nf])';
    if parms.norm_flag % normalize data to b=0
      y0 = repmat(reshape(MIfit.volb0(:,:,z),[parms.nx*parms.ny,1])',[parms.nf 1]);
      y = y./max(y0,eps);
    end;
    betas = MIfit.Ainv*y;
    betas = reshape(betas',[parms.nx parms.ny MIfit.nb]);
    betas(isnan(betas)) = 0;
    MIfit.volMI(:,:,z,:) = betas;
  end
  % dilate brain mask
  MIfit.volmask_dilated = dilate_brain_mask(MIfit.volmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MIfit = nlfit_iso(vol,MIfit,parms)
  % set lower and upper bounds
  [lb,ub] = set_bounds(MIfit,parms);
  % loop over all voxels inside mask
  ind_mask = find(MIfit.volmask_dilated);
  for m=1:length(ind_mask)
    [i,j,k]=ind2sub([parms.nx,parms.ny,parms.nz],ind_mask(m));
    % get data vector for this voxel
    data_vec = double(squeeze(vol(i,j,k,:)));
    if parms.norm_flag
      y0 = double(MIfit.volb0(i,j,k));
    else
      y0 = 1;
    end;
    % linear fit as initial estimate
    beta_lin = squeeze(MIfit.volMI(i,j,k,:));
    % nonlinear fit
    switch parms.nlfit_method
      case 'lsqnonlin'
        beta_nonlin = lsqnonlin(@(beta) costfunc(beta,...
                        MIfit.A,data_vec,y0,parms),...
                        beta_lin,lb,ub,parms.nlfit_options);
      case 'fmincon'
        beta_nonlin = fmincon(@(beta) sum(costfunc(beta,...
                        MIfit.A,data_vec,y0,parms).^2),...
                        beta_lin,[],[],[],[],lb,ub,[],parms.nlfit_options);
    end;
    MIfit.volMI(i,j,k,:) = beta_nonlin;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = costfunc(beta,A,y,y0,parms)
  cost = [];
  yh = y0*A*beta;
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
