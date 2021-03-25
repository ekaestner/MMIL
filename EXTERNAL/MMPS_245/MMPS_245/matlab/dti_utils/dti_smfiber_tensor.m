function [vol_fiber_sm,vol_tensor_sm]=dti_smfiber_tensor(vol_fiber,vol_tensor,vol_weight,...
  sigmavec,win,init_niter,final_niter,first_only_flag)
%function [vol_fiber_sm,vol_tensor_sm]=dti_smfiber_tensor(vol_fiber,vol_tensor,[vol_weight],...
%  [sigmavec],win,init_niter,final_niter,[first_only_flag])
%
% Required Parameters:
%   vol_fiber: 3D volume specifying location of fiber bundle
%   vol_tensor: 5D volume containing diffusion tensors (size = [nx,ny,nz,3,3])
%            or 4D volume containing diffusion tensors (size = [nx,ny,nz,9])
%            or 4D volume containing diffusion tensors (size = [nx,ny,nz,6])
%
% Optional Parameters:
%   vol_weight: 3D volume with weightings specifying confidence in orientations
%     (e.g. FA, DR1)
%     if empty or omitted, all non-zero voxels in vol_fiber will have weight = 1
%     {default = []}
%   sigmavec: 3-member vector of smoothing sigmas for x, y, and z (voxels)
%     {default = [5,5,5]}
%   win: smoothing window for fiber orientations (voxels)
%     {default = 1}
%   init_niter: max number of iterations for inital orientation smoothing
%     (filling in values only for those not in vol_weight)
%     {default = 10}
%   final_niter: number of iterations for final orientation smoothing
%     {default = 1}
%   first_only_flag: output tensor with 2nd and 3rd eigen vectors removed
%     {default = 0}
%
%
% Created:  03/07/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%% todo: use varargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_fiber_sm = [];
vol_tensor_sm = [];

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('vol_weight','var'), vol_weight = []; end;
if ~exist('sigmavec','var') | isempty(sigmavec), sigmavec = [5,5,5]; end;
if ~exist('win','var') | isempty(win), win = 1; end;
if ~exist('init_niter','var') | isempty(init_niter), init_niter = 10; end;
if ~exist('final_niter','var') | isempty(final_niter), final_niter = 1; end;
if ~exist('first_only_flag','var') | isempty(first_only_flag), first_only_flag = 0; end;
if ~exist('verbose_flag','var') | isempty(verbose_flag), verbose_flag = 1; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volsz = size(vol_tensor);
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
        vol_tensor = reshape(vol_tensor,[volsz(1:3),9]);
      else
        error('input vol_tensor is 5D but last two dims are not 3x3');
      end;
    end;
  otherwise
    error('input vol_tensor must be 4D or 5D (is %dD)',ndims);
end;
nvox = prod(volsz(1:3));

if isempty(vol_weight)
  vol_weight = ones(size(vol_fiber));
end;
vol_weight(vol_fiber<=0) = 0;

% smooth fiber mask
sx = sigmavec(1);
sy = sigmavec(2);
sz = sigmavec(3);
vol_fiber_sm = mmil_smooth3d(vol_fiber,sx,sy,sz);
vol_fiber_sm = 1.0*(vol_fiber_sm>0.01);

% smooth orientations
i_trust = find(vol_weight>0);
i_unknown = find(vol_fiber_sm>0 & vol_weight==0);
if verbose_flag
  fprintf('%s: num unknown = %d\n',mfilename,length(i_unknown));
end;
vecT = reshape(vol_tensor,[nvox,ncomp]);
vecT_sm = zeros(size(vecT));
vecT_sm(i_trust,:) = vecT(i_trust,:);
vol_tensor_sm = reshape(vecT_sm,[nx,ny,nz,ncomp]);
still_unknown = zeros(size(i_unknown));
for i=1:init_niter
  if verbose_flag
    fprintf('%s: smooth step %d\n',mfilename,i);
  end;
  for k=1:length(i_unknown)
    [x,y,z]=ind2sub(size(vol_weight),i_unknown(k));
    x0 = max([1,x-win]);
    x1 = min([nx,x+win]);
    y0 = max([1,y-win]);
    y1 = min([ny,y+win]);
    z0 = max([1,z-win]);
    z1 = min([nz,z+win]);
    tmpT = squeeze(vol_tensor_sm(x0:x1,y0:y1,z0:z1,:));
    wt = squeeze(vol_weight(x0:x1,y0:y1,z0:z1));
    tmp_nvox = length(x0:x1)*length(y0:y1)*length(z0:z1);    
    tmpT = reshape(tmpT,[tmp_nvox,ncomp])';
    wt = reshape(wt,[tmp_nvox,1]);
    sum_wt = sum(wt);
    if sum(abs(tmpT(:)))<=0
      still_unknown(k) = 1;
      continue;
    else
      still_unknown(k) = 0;
    end;
    tmpT = (tmpT*wt)/sum(wt);
    vol_tensor_sm(x,y,z,:) = tmpT;
    vol_weight(x,y,z) = mean(wt);
  end;
  num_still_unknown = length(find(still_unknown));
  if verbose_flag
    fprintf('%s: num still unknown = %d\n',mfilename,num_still_unknown);
  end;
  if num_still_unknown == 0, break; end;
end;

i_unknown = find(vol_fiber>0);
if verbose_flag
  fprintf('%s: num voxels for final = %d\n',mfilename,length(i_unknown));
end;
tmp_vol_tensor_sm = vol_tensor_sm;
for i=1:final_niter
  if verbose_flag
    fprintf('%s: final smooth step %d\n',mfilename,i);
  end;
  for k=1:length(i_unknown)
    [x,y,z]=ind2sub(size(vol_weight),i_unknown(k));
    x0 = max([1,x-win]);
    x1 = min([nx,x+win]);
    y0 = max([1,y-win]);
    y1 = min([ny,y+win]);
    z0 = max([1,z-win]);
    z1 = min([nz,z+win]);
    tmpT = squeeze(vol_tensor_sm(x0:x1,y0:y1,z0:z1,:));
    wt = squeeze(vol_weight(x0:x1,y0:y1,z0:z1));
    tmp_nvox = length(x0:x1)*length(y0:y1)*length(z0:z1);    
    tmpT = reshape(tmpT,[tmp_nvox,ncomp])';
    wt = reshape(wt,[tmp_nvox,1]);
    sum_wt = sum(wt);
    if sum(abs(tmpT(:)))<=0
      continue;
    end;
    tmpT = (tmpT*wt)/sum(wt);
    if first_only_flag
      % reshape to tensor
      T = dti_components2tensor(tmpT);
      % calculate eigen vectors
      [U,S] = eig(T);
      lamdas = real(diag(S)');
      [sorted_lamdas,ind]=sort(lamdas,2,'descend');
      v = U(:,ind(1))';
      % use only first eigen vector
      T = v'*v;
      tmpT = dti_tensor2components(T,ncomp);
    end;
    tmp_vol_tensor_sm(x,y,z,:) = tmpT;
  end;
end;
vol_tensor_sm = tmp_vol_tensor_sm;

vol_fiber_sm(vol_fiber>0) = 2.0;
vol_fiber_sm(i_trust) = 3.0;


