function DTmeas = dti_calc_eig(volT,volmask)
%function DTmeas = dti_calc_eig(volT,[volmask])
%
% Required Parameters:
%   volT:  5D volume containing tensor matrix
%       or 4D volume containing tensor vector
%
% Optional Parameters:
%   volmask: mask volume to restrict analysis to portion of image {default=[]}
%
% Output:
%   DTmeas: structure containing tensor calculations with fields:
%     volE : eigen values ; size = [nx,ny,nz,3]
%     volV : eigen vectors ; size = [nx,ny,nz,3,3]
%     volFA : fractional anisotropy ; size = [nx,ny,nz]
%
% Created:  08/27/07 Don Hagler
% Last Mod: 10/27/12 Don Hagler
%

%     volT : tensor matrix ; size = [nx,ny,nz,3,3]
%     volDR1 : dominance ratio (1st vs. 2nd eigen value) ; size = [nx,ny,nz]
%     volDR2 : dominance ratio (2nd vs. 3rd eigen value) ; size = [nx,ny,nz]
%     volADC : apparent difussion coeficient ; size = [nx,ny,nz]

if ~mmil_check_nargs(nargin,1), return; end;

smf = 10^-5;
nvox_msg = 1000;
if ~exist('volmask','var'), volmask = []; end;

DTmeas = [];

% check input dimensions
volsz = size(volT);
ndims = length(volsz);
switch ndims
  case {4,5}
    nx = volsz(1);
    ny = volsz(2);
    nz = volsz(3);
    if ndims == 4
      ncomp = volsz(4);
    else
      if volsz(4) == 3 & volsz(5) == 3
        ncomp = 9;
        volT = reshape(volT,[volsz(1:3),9]);
      else
        error('input volT is 5D but last two dims are not 3x3');
      end;
    end;
  otherwise
    error('input volT must be 4D or 5D (is %dD)',ndims);
end;
nvox = prod(volsz(1:3));

% initialize output
DTmeas.volE = zeros(nx,ny,nz,3); % eigen values
DTmeas.volV = zeros(nx,ny,nz,3,3); % eigen vectors
%DTmeas.volADC = zeros(nx,ny,nz); % apparent diffusion coefficient
DTmeas.volFA = zeros(nx,ny,nz); % fractional anisotropy
%DTmeas.volDR1 = zeros(nx,ny,nz); % (E1-E2)/E1
%DTmeas.volDR2 = zeros(nx,ny,nz); % (E2-E3)/E2


for k=1:nz
%  fprintf('%s: slice %d\n',mfilename,k);
  sliceT = squeeze(volT(:,:,k,:));
  % if original data are 0's, don't calculate anything for those voxels
  sliceT_std = sum(abs(sliceT),3);
  if ~isempty(volmask)
    voxels = find(sliceT_std>smf & squeeze(volmask(:,:,k)));
  else
    voxels = find(sliceT_std>smf);
  end;

  nvox_ROI = length(voxels);
%  fprintf('%s: calculating eigenvectors for %d selected voxels\n',mfilename,nvox_ROI);
%  tic
  for m=1:nvox_ROI
    if ~rem(m,nvox_msg)
%      dt = toc;
%      fprintf('%s:   voxel %d of %d   elapsed time = %0.1f seconds\n',...
%        mfilename,m,nvox_ROI,dt);
%      tic
    end;
    [i,j]=ind2sub([nx,ny],voxels(m));
    vec = squeeze(sliceT(i,j,:));
    T = dti_components2tensor(vec);

    if sum(abs(T))<smf, continue; end;

    % find sorted eigen vectors
    [U,S] = eig(T);
    evals = real(diag(S)');
    [evals,ind]=sort(evals,2,'descend');
    evals(evals<0)=0;
    emean = mean(evals);

    % store eigen values and vectors
    DTmeas.volE(i,j,k,:) = evals;
    for vnum=1:3
      DTmeas.volV(i,j,k,vnum,:) = U(:,ind(vnum))';
    end;

    % calculate FA
    evals_diff = evals - emean;
    numer = sum(evals_diff.^2);
    denom = sum(evals.^2);
    if denom~=0
      DTmeas.volFA(i,j,k) = sqrt((3/2)*(numer/denom));
    end;

    if 0
      % calculate ADC
      DTmeas.volADC(i,j,k) = sum(evals);

      % calculate dominance ratio 1
      numer = evals(1)-evals(2);
      denom = evals(1)+eps;
      if denom~=0
        DTmeas.volDR1(i,j,k) = numer/denom;
      end;

      % calculate dominance ratio 2
      numer = evals(2)-evals(3);
      denom = evals(2)+eps;
      if denom~=0
        DTmeas.volDR2(i,j,k) = numer/denom;
      end;
    end;
  end
end;

return;

