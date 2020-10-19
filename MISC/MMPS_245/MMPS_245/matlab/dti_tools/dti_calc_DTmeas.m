function DTmeas = dti_calc_DTmeas(DTfit,mask_flag)
%function DTmeas = dti_calc_DTmeas(DTfit,mask_flag)
%
% Required Input:
%   DTfit: output structure from dti_fit_tensor with fields:
%     volDT : tensor fit parameters    size = [nx,ny,nz,7]
%     volb0 : average b=0 image        size = [nx,ny,nz]
%     volmask : mask volume            size = [nx,ny,nz]
%     volmask_dilated : dilated mask   size = [nx,ny,nz]
%     volsz : size of input 4D vol = [nx,ny,nz,nframes]
%     Q: tensor fit forward matrix
%     Qinv: tensor fit inverse matrix
%     qmat: input q matrix
%     bvals: input b values
%
% Optional Parmaeters:
%   mask_flag: [0|1] apply dilated brain mask to output volumes
%     {default = 1}
%
% Output:
%   DTmeas: structure containing DTI measures with fields:
%     volb0:  average b=0 image        size = [nx,ny,nz]
%     volT : tensor matrix             size = [nx,ny,nz,3,3]
%     volV : eigen vectors             size = [nx,ny,nz,3,3]
%     volE : eigen values              size = [nx,ny,nz,3]
%     volFA : fractional anisotropy    size = [nx,ny,nz]
%     volmask : dilated brain mask     size = [nx,ny,nz] 
%     volsz : size of original data         = [nx,ny,nz,nf]
%
% Created:  02/15/10 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('mask_flag','var') || isempty(mask_flag)
  mask_flag = 1;
end;
DTmeas = [];


DTmeas.volsz = DTfit.volsz;
volsz = DTfit.volsz(1:3);

% initialize output
DTmeas.volT = zeros([volsz,3,3]);
DTmeas.volV = zeros([volsz,3,3]);
DTmeas.volE = zeros([volsz,3]);

DTmeas.volmask = DTfit.volmask_dilated;
if mask_flag
  % voxels in mask
  voxels = find(DTmeas.volmask);
  nvox = length(voxels);
else
  voxels = find(DTfit.volDT(:,:,:,1));
  nvox = length(voxels);
end;

for m=1:nvox
  % reshape tensor fit parameters into tensor
  [i,j,k]=ind2sub(volsz,voxels(m));
  T = dti_components2tensor(DTfit.volDT(i,j,k,1:6));
  DTmeas.volT(i,j,k,:,:) = T;

  % calculate sorted Eigenvalues and Eigenvectors
  [U,S] = eig(T);
  Evals = real(diag(S)');
  [Evals,ind]=sort(Evals,2,'descend'); % sort Eigenvalues
  Evals(Evals<0)=0; % set negative Eigenvalues to zero

  % Eigenvalues
  DTmeas.volE(i,j,k,:) = Evals;
 
  % Eigenvectors
  for vnum=1:3
    DTmeas.volV(i,j,k,vnum,:) = U(:,ind(vnum))';
  end;
end

% fractional anisotropy
volE_mean = mean(DTmeas.volE,4);
volE_diff = DTmeas.volE;
for i=1:3
  volE_diff(:,:,:,i) = volE_diff(:,:,:,i) - volE_mean;
end;
numer = sum(volE_diff.^2,4);
denom = sum(DTmeas.volE.^2,4) + eps;
DTmeas.volFA = sqrt((3/2)*(numer./denom));

if mask_flag
  DTmeas.volb0 = DTfit.volb0.*DTmeas.volmask;
else
  DTmeas.volb0 = DTfit.volb0;
end;


return;
