function [vol_corr,DTfit] = dti_eddycurr_corr(vol,qmat,varargin)
%function [vol_corr,DTfit] = dti_eddycurr_corr(vol,qmat,[options])
%
% Required Parameters:
%   vol: 4d DTI volume
%     may be empty if fname is supplied
%   qmat: matrix of diffusion direction vectors
%
% Optional Parameters for input and output:
%   'fname_in': input file name
%     if supplied, vol is ignored and may be empty
%     {default = []}
%   'fname_out': output file name containing corrected volume
%     if not supplied, output file is not created
%     {default = []}
%   'fname_DT': output file name containing tensor fit structure
%     if not supplied, output file is not created
%     {default = []}
%   'forceflag': [0|1] overwrite existing fname_out
%     {default = 0}
%
% Optional Parameters for eddy current corrections:
%   'niter': number of iterations of eddy current correction
%     {default = 5}
%   'censor_niter': number of iterations of censoring on first ecc iteration
%     {default = 3}
%   'censor_thresh': error threshold for censoring bad frames
%     normalized to median error for each slice
%     higher values mean less censoring
%     {default = 3.2}
%   'censor_min_ndirs': minimum number of diffusion directions (not including
%     b=0 images) required for tensor fit after censoring
%     will not do censoring if it means reducing number of directions below min
%     {default = 6}
%   'nonlin_flag': [0|1] use nonlinear optimization for tensor fits
%     {default = 0}
%   'b0_thresh': threshold used for considering a b-value to be 0
%     {default = 10}
%   'bvals': vector of b values (one for each diffusion direction)
%     If single value supplied, will use same for all
%     {default = [1000]}
%   'nparams': number of parameters for fit
%     3=linear terms only, 9=linear plus cross-terms
%     {default = 3}
%   'phasedir': defines direction of expected distortions
%     1=x, 2=y, 3=z
%     {default = 2}
%   'driftcorr': [0|1] add drift correction term to qdiff
%     {default = 0}
%
% Created:  01/26/09 by Don Hagler
% Last Mod: 07/09/15 Don Hagler
%

% based on code from Anders Dale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_corr=[];
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms( varargin, { ...
  'fname_in',[],[],...
  'fname_out',[],[],...
  'fname_DT',[],[],...
  'forceflag',false,[false true],...
...
  'niter',5,[0,100],...
  'censor_niter',3,[0,100],...
  'censor_min_ndirs',6,[],...
  'censor_thresh',3.2,[],...
  'nonlin_flag',false,[false true],...
  'b0_thresh',10,[0,500],...
  'bvals',1000,[],...
  'nparams',3,[3,9],...
  'phasedir',2,[1,2,3],...
  'driftcorr',false,[false true],...
... % hidden parameters
  'smf',1e-5,[1e-100,1e-1],...
...
  'tensor_tags',{'bvals','censor_min_ndirs','censor_thresh',...
                 'nonlin_flag','b0_thresh'},[],...
});

if output_exists(parms), return; end;

if ~isempty(parms.fname_in)
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [vol,M] = fs_load_mgh(parms.fname_in);
else
  M = [];
end;

if isempty(vol)
  error('vol is empty');
end;

vol = single(vol);
% corrections below assume displacements are along 1st dimension
switch parms.phasedir
  case 2
    vol = permute(vol,[2,1,3,4]);
  case 3
    vol = permute(vol,[3,2,1,4]);
end;

% check input dimensions
[nx,ny,nz,nf] = size(vol);

if size(qmat,1)~=nf
  error('number of directions in qmat must match frames in vol');
end;
if size(qmat,2)~=3
  error('qmat must have 3 columns');
end;

qlength = sqrt(sum(qmat.^2,2));
i_b0 = find(qlength<parms.smf);
nb0 = length(i_b0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute unit deformation field for each deformation parameter
colnum = 0;
qdiff = zeros(nf,parms.nparams);
for i = 1:3 % cardinal axes
  colnum = colnum+1;
  qdiff(:,colnum) = qmat(:,i);
end
if parms.nparams==9
  for i = 1:3 % cross-terms (should check if these are needed)
    for j = i:3
      colnum = colnum+1;
      qdiff(:,colnum) = qmat(:,i).*qmat(:,j);
    end
  end
end
% add drift correction term 
if parms.driftcorr
  ramp = [0 1:(size(qdiff,1)-1)]';
  qdiff = cat(2,qdiff,ramp);
  parms.nparams = parms.nparams+1;
end

[rowvol,colvol,slicevol] = ndgrid(single(1:nx),single(1:ny),single(1:nz));
coordvol = cat(4,rowvol,colvol,slicevol,ones(size(rowvol)));

vol_corr = vol; % start with uncorrected data
paramvec = 0;
censor_niter = parms.censor_niter;
censor_mat = [];
for iter = 0:parms.niter
  % fit tensor, censoring bad frames for each slice
  args = mmil_parms2args(parms,parms.tensor_tags);
  DTfit = dti_fit_tensor(vol_corr,qmat,args{:},...
    'censor_niter',censor_niter,'censor_mat',censor_mat);
  % synthesize volume using tensor fit
  vol_err = dti_synth_vol(DTfit); % save memory by reusing variable later

  % replace censored slices in vol_corr
  if parms.censor_niter>0
    censor_mat = DTfit.censor_mat;
    censor_niter = 1; % next time use pre-determined censoring to init
    vol_corr = dti_censor_vol(vol_corr,vol_err,censor_mat,...
     'censor_min_ndirs',parms.censor_min_ndirs,...
      'volmask',DTfit.volmask,'nb0',nb0);
  end;
  % calculate error
  vol_err = vol_corr - vol_err;
  cost = sum(vol_err(:).^2)/sum(vol_corr(:).^2);
  fprintf('%s: iter=%d: cost=%0.6e\n',mfilename,iter,cost);

  dTdy = (cat(1,zeros(1,ny,nz,nf),diff(vol_corr,1,1))+...
          cat(1,diff(vol_corr,1,1),zeros(1,ny,nz,nf)))/2;

  nbetas = parms.nparams*size(coordvol,4); % total number of parameters to solve for
  g = 0;
  H = 0;
  for k = 1:nf
    colnum = 0;
    tmp0 = zeros(nx,ny,nz,nbetas);
    tmp1 = zeros(nx,ny,nz,nbetas);
    for i = 1:size(coordvol,4) % coordinates plus 1
      for j = 1:parms.nparams
        colnum = colnum+1;
        tmp1(:,:,:,colnum) = squeeze(dTdy(:,:,:,k).*coordvol(:,:,:,i)*qdiff(k,j));
        tmp0(:,:,:,colnum) = squeeze(vol_err(:,:,:,k)).*tmp1(:,:,:,colnum);
      end
    end
    tmp0 = reshape(tmp0,[nx*ny*nz nbetas]);
    tmp1 = reshape(tmp1,[nx*ny*nz nbetas]);
    g = g + sum(tmp0,1)';
    H = H + tmp1'*tmp1;
  end
  clear dTdy vol_err
  betahat = -H\g;
  paramvec = paramvec + betahat;
  dymat = zeros(nx,ny,nz,nf);
  colnum = 0;
  for i = 1:size(coordvol,4) % coordinates plus 1
    for j = 1:parms.nparams
      colnum = colnum+1;
      for k = 1:nf
        dymat(:,:,:,k) = dymat(:,:,:,k)+coordvol(:,:,:,i)*qdiff(k,j)*paramvec(colnum);
      end
    end
  end
  for i = 1:ny
    for j = 1:nz
      for k = 1:nf
        vol_corr(:,i,j,k) = mmil_fast_interp1(vol(:,i,j,k),rowvol(:,i,j)+dymat(:,i,j,k),'*linear',0);
      end
    end
  end
  clear dymat
end

switch parms.phasedir
  case 2
    vol_corr = permute(vol_corr,[2,1,3,4]);
  case 3
    vol_corr = permute(vol_corr,[3,2,1,4]);
end;

if ~isempty(parms.fname_out)
  fs_save_mgh(vol_corr,parms.fname_out,M);
end;

if ~isempty(parms.fname_DT)
  save(parms.fname_DT,'DTfit');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end;
  if ~isempty(parms.fname_out) && ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if ~isempty(parms.fname_DT) && ~exist(parms.fname_DT,'file')
    exist_flag = 0;
  end;
return;

