function DTfit = dti_fit_tensor(vol,qmat,varargin)
%function DTfit = dti_fit_tensor(vol,qmat,[options])
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
%     if empty, will create an approximate brain mask from b=0 image
%     {default = []}
%   'M': 4x4 vox2ras matrix for vol (recommended for mask creation)
%     {default = identity matrix}
%   'nob0_flag': [0|1] whether to exclude b=0 images from fits
%     {default = 0}
%   'max_bval': maximum b-value used in tensor fit
%     {default = Inf}
%   'censor_niter': number of iterations of censoring to perform
%     {default = 0}
%   'censor_mat': matrix of slices by frame to exclude from fit
%     ignored if censor_niter = 0
%     {default = []}
%   'censor_thresh': error threshold for censoring bad frames
%     normalized to median error for each slice
%     higher values mean less censoring
%     {default = 3.2}
%   'censor_min_ndirs': minimum number of diffusion directions (not including
%     b=0 images) required for tensor fit after censoring
%     will not do censoring if it means reducing number of directions below min
%     {default = 6}
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
%   DTfit: structure containing tensor calculations with fields:
%     volDT : tensor fit parameters    size = [nx,ny,nz,7]
%     volb0 : average b=0 image        size = [nx,ny,nz]
%     volmask : mask volume            size = [nx,ny,nz]
%     volmask_dilated : dilated mask   size = [nx,ny,nz]
%     volsz : size of input 4D vol          = [nx,ny,nz,nf]
%     M : input vox2ras matrix         size = [4,4]
%     qmat : input q matrix            size = [nf,3]
%     bvals : input b values           size = [nf,1]
%     Q : tensor fit forward matrix    size = [nf,7]
%     censor_mat : censored slices per frame       size = [nz,nf]
%     censor_err : fit error per slices per frame  size = [nz,nf] 
%     censor_niter : number of censoring iterations
%
% Created:  11/02/06 by Don Hagler
% Last Mod: 06/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DTfit = [];

% check input arguments
if ~mmil_check_nargs(nargin, 2), return; end;
parms = check_input(varargin);

% initialize output struct
[DTfit,vol,parms] = init_DTfit(vol,qmat,parms);

% linear tensor fit
if parms.censor_niter==0
  DTfit = fit_tensor(vol,DTfit,parms);
else
  DTfit = fit_tensor_censor(vol,DTfit,parms);
end;

% nonlinear tensor fit
if parms.nonlin_flag
  DTfit = nlfit_tensor(vol,DTfit,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'bvals',1000,[],...
    'volmask',[],[],...
    'M',eye(4),[],...
    'nob0_flag',false,[false true],...
    'max_bval',Inf,[100,Inf],...
    'censor_niter',0,[0,100],...
    'censor_mat',[],[],...
    'censor_min_ndirs',6,[],...
    'censor_thresh',3.2,[],...
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
    'smf',10^-5,[10^-100,10^-1],...
  ... % parameters to be passed
    'mask_tags',{'log_flag' 'thresh' 'fill1_smooth1' 'fill1_thresh1'...
                 'fill1_smooth2' 'fill1_thresh2' 'fill1_smooth3'...
                 'fill1_erode_flag' 'clip_edges_flag' 'clip_edges_width'...
                 'fill2_smooth1' 'fill2_thresh1' 'fill2_smooth2' 'fill2_thresh2'...
                 'fill2_smooth3' 'fill2_thresh3' 'fill2_erode_flag' 'binary_flag'...
                 'forceflag'},[],...
  });

  % initialize nonlinear fit output
  if parms.nonlin_flag
    try % for R2009b
      parms.nlfit_options = optimset(...
        'Display','none','Algorithm','levenberg-marquardt',...
        'MaxIter',100,'TolFun',1e-10);
    catch % for R2007a
      parms.nlfit_options = optimset(...
        'Display','none','LargeScale','off',...
        'MaxIter',100,'TolFun',1e-10);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DTfit,vol,parms] = init_DTfit(vol,qmat,parms)
  DTfit = [];
  DTfit.volsz = size(vol);
  [parms.nx,parms.ny,parms.nz,parms.nf] = size(vol);

  % extend bvals to nf elements if only 1 supplied
  if numel(parms.bvals)==1
    DTfit.bvals = parms.bvals*ones(parms.nf,1);
  elseif numel(parms.bvals)==parms.nf
    DTfit.bvals = reshape(parms.bvals,[parms.nf,1]);
  else
    error('number of bvals must match frames in vol');
  end;

  % check qmat is correct size
  if size(qmat,1)~=parms.nf
    error('number of directions in qmat must match frames in vol');
  end;
  if size(qmat,2)~=3
    error('qmat must have 3 columns');
  end;

  % ensure that qmat contains unit vectors
  qlength = max(sqrt(sum(qmat.^2,2)),eps);
  % set correct bvals if multishell
  if length(unique(DTfit.bvals))==1 && length(unique(qlength(qlength>0)))>1
    DTfit.bvals = DTfit.bvals.*(qlength.^2);
  end;
  DTfit.qmat_orig = qmat;
  DTfit.qmat = bsxfun(@rdivide,qmat,qlength);

  if parms.censor_niter>0
    if isempty(parms.censor_mat)
      DTfit.censor_mat = zeros(parms.nz,parms.nf);
    else
      DTfit.censor_mat = parms.censor_mat;
    end;
  else
    DTfit.censor_mat = [];
    DTfit.censor_err = [];
  end
  DTfit.censor_niter = parms.censor_niter;

  DTfit.M = parms.M;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  qlength = sqrt(sum(DTfit.qmat.^2,2));
  i_b0 = find(qlength<parms.smf | DTfit.bvals<=parms.b0_thresh);
  i_bx = setdiff([1:parms.nf],i_b0);
  nb0 = length(i_b0);
  volb0 = vol(:,:,:,i_b0);
  if numel(parms.bvals)==1
    DTfit.bvals(i_b0) = 0;
  end;
  % calculate scaling factors and apply to each frame
  if parms.scalefacts_flag && length(i_b0)>1
    DTfit.b0_scalefacts = dti_calc_b0_scale(volb0);
    fprintf('%s: scaling factors from %d b=0 images: %s\n',...
      mfilename,nb0,sprintf('%0.2f  ',DTfit.b0_scalefacts));
    vol = dti_apply_b0_scale(vol,i_b0,DTfit.b0_scalefacts);
  else
    DTfit.b0_scalefacts = 1;
  end;

  % average b=0 images
  if nb0>0
    DTfit.volb0 = mean(volb0,4);
  else
    DTfit.volb0 = [];
  end;
  clear volb0;

  if isempty(parms.volmask)
    if ~isempty(DTfit.volb0)
      % create approximate brain mask from average b=0 image
      DTfit.volmask = create_brain_mask(DTfit.volb0,parms,1);
    else
      % exclude voxels with zeros
      DTfit.volmask = ones(parms.nx,parms.ny,parms.nz);
      tmp = min(vol,[],4);
      DTfit.volmask(tmp==0) = 0;
    end;
  else
    DTfit.volmask = parms.volmask;
  end;

  DTfit.i_fit = [1:DTfit.volsz(4)];
  
  % exclude b=0 images
  if parms.nob0_flag
    i_keep = find(DTfit.bvals > parms.b0_thresh);
    parms.nf = length(i_keep);
    vol = vol(:,:,:,i_keep);
    DTfit.volsz = size(vol);
    DTfit.bvals = DTfit.bvals(i_keep);
    DTfit.qmat = DTfit.qmat(i_keep,:);
    DTfit.i_fit = DTfit.i_fit(i_keep);
  end;

  % exclude high b-value images
  if parms.max_bval < max(DTfit.bvals)
    i_keep = find(DTfit.bvals <= parms.max_bval);
    parms.nf = length(i_keep);
    vol = vol(:,:,:,i_keep);
    DTfit.volsz = size(vol);
    DTfit.bvals = DTfit.bvals(i_keep);
    DTfit.qmat = DTfit.qmat(i_keep,:);
    DTfit.i_fit = DTfit.i_fit(i_keep);
  end;

  fprintf('%s: including %d directions in tensor fit...\n',...
    mfilename,length(DTfit.i_fit));
  
  DTfit.i_b0 = find(DTfit.bvals <= parms.b0_thresh);
  DTfit.i_bx = find(DTfit.bvals > parms.b0_thresh);
  DTfit.nb0 = length(DTfit.i_b0);
  
  % create forward matrix for fitting tensor
  DTfit.Q = create_forward_matrix(DTfit);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTfit = fit_tensor(vol,DTfit,parms)
  DTfit.volDT = zeros(parms.nx,parms.ny,parms.nz,7);
  for z = 1:parms.nz
    data_mat = reshape(vol(:,:,z,:),[parms.nx*parms.ny,parms.nf])';
    T = DTfit.Q\log(max(1,data_mat));
    DTfit.volDT(:,:,z,:) = reshape(T',[parms.nx parms.ny 7]);
  end;
  % create approximate brain mask from synthesized b=0 image
  if isempty(parms.volmask)
    DTfit.volmask = create_brain_mask(DTfit.volDT(:,:,:,end),parms,0);
  end;
  % dilate brain mask
  DTfit.volmask_dilated = dilate_brain_mask(DTfit.volmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTfit = fit_tensor_censor(vol,DTfit,parms)
  for iter=1:parms.censor_niter
    DTfit.volDT = zeros(parms.nx,parms.ny,parms.nz,7);
    mean_err_mat = zeros(parms.nz,parms.nf);
    mean_sig_mat = zeros(parms.nz,parms.nf);
    for z = 1:parms.nz
      ind_keep = ~(DTfit.censor_mat(z,:)==1);
      nkeep = length(find(ind_keep));
      % do not censor if it means we have too few directions
      if nkeep < parms.censor_min_ndirs + DTfit.nb0
        nkeep = parms.nf;
        ind_keep = 1:parms.nf;
      end;
      data_mat = reshape(vol(:,:,z,ind_keep),[parms.nx*parms.ny,nkeep])';
      % clip values between 1 and b=0 value
      %  (to prevent log(0) and negative diffusion)
      if DTfit.nb0>0
        b0_mat = ones(nkeep,1)*reshape(DTfit.volb0(:,:,z),[parms.nx*parms.ny,1])';
        data_mat = max(1,min(b0_mat,data_mat));
      end;
      Qtmp = DTfit.Q(ind_keep,:);
      T = Qtmp\log(data_mat);
      DTfit.volDT(:,:,z,:) = reshape(T',[parms.nx parms.ny 7]);
      fit_mat = exp(DTfit.Q*T);
      data_mat = reshape(vol(:,:,z,:),[parms.nx*parms.ny,parms.nf])';
      err_mat = data_mat - fit_mat;
      mask_vec = mmil_rowvec(squeeze(DTfit.volmask(:,:,z)));
      for f=1:parms.nf
        sig_vec = fit_mat(f,:);
        err_vec = err_mat(f,:);
        sig_vec = sig_vec(mask_vec>0);
        err_vec = err_vec(mask_vec>0);
        if length(err_vec)>100
          mean_err_mat(z,f) = sqrt(mean(abs(err_vec).^2));
          mean_sig_mat(z,f) = sqrt(mean(abs(sig_vec).^2));
        end;
      end;
    end;

    % normalize by median error within slice across frames
    DTfit.raw_censor_err = mean_err_mat;
    med_err = median(DTfit.raw_censor_err (:,DTfit.i_bx),2);
    tmp = max(eps,med_err*ones(1,parms.nf));
    DTfit.censor_err = DTfit.raw_censor_err ./tmp;
    DTfit.censor_err(:,DTfit.i_b0) = 0;

    % calculate overall cost and update DTfit.censor_mat
    cost_vec = mmil_rowvec(DTfit.censor_err);
    DTfit.censor_mat = (DTfit.censor_err>parms.censor_thresh); 
    cost_cens = mean(cost_vec(find(DTfit.censor_mat==0 & DTfit.censor_err~=0)));
    cost = mean(cost_vec(find(DTfit.censor_err~=0)));
    fprintf('%s: iter=%d: cost=%f (cost_cens=%f)\n',mfilename,iter,cost,cost_cens);
    % create approximate brain mask from synthesized b=0 image
    if isempty(parms.volmask)
      DTfit.volmask = create_brain_mask(DTfit.volDT(:,:,:,end),parms,0);
    end;
  end;
  % dilate brain mask
  DTfit.volmask_dilated = dilate_brain_mask(DTfit.volmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTfit = nlfit_tensor(vol,DTfit,parms)
  % loop over all voxels inside mask
  ind_mask = find(DTfit.volmask_dilated);
  for m=1:length(ind_mask)
    [i,j,k]=ind2sub([parms.nx,parms.ny,parms.nz],ind_mask(m));
    % use censor information if available
    if parms.censor_niter > 0
      ind_keep = ~(DTfit.censor_mat(k,:)==1);
      nkeep = length(find(ind_keep));
      % do not censor if it means we have too few directions
      if nkeep < parms.censor_min_ndirs + DTfit.nb0
        nkeep = parms.nf;
        ind_keep = 1:parms.nf;
      end;
    else
      ind_keep = 1:parms.nf;
    end;
    % get data vector for this voxel
    data_vec = double(squeeze(vol(i,j,k,ind_keep)));
    % resize forward matrix for valid frames
    Qtmp = DTfit.Q(ind_keep,:);
    % linear fit as initial estimate
    beta_lin = squeeze(DTfit.volDT(i,j,k,:));
    % nonlinear fit
    beta_nonlin = lsqnonlin(@(beta) data_vec-exp(Qtmp*beta),...
                    beta_lin,[],[],parms.nlfit_options);
    DTfit.volDT(i,j,k,:) = beta_nonlin;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Q = create_forward_matrix(DTfit)
  ndirs = size(DTfit.qmat,1);
  Q = zeros(ndirs,7);
  for i=1:ndirs
    outerprod = -DTfit.bvals(i)*DTfit.qmat(i,:)'*DTfit.qmat(i,:);
    Q(i,1:3) = diag(outerprod)';
    Q(i,4) = 2*outerprod(1,2);
    Q(i,5) = 2*outerprod(1,3);
    Q(i,6) = 2*outerprod(2,3);
  end;
  Q(:,7) = 1;
return;

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

