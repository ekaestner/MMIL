function [vol_fiber_sm,vol_V_sm]=dti_smoothfiber(vol_fiber,vol_V,vol_mask,sigmavec,...
  win,init_niter,final_niter)
%function [vol_fiber_sm,vol_V_sm]=dti_smoothfiber(vol_fiber,vol_V,[vol_mask],[sigmavec],...
%  win,init_niter,final_niter)
%
% Required Parameters:
%   vol_fiber: 3D volume specifying location of fiber bundle
%   vol_V: 3D volume containing fiber orientations (size = [nx,ny,nz,3])
%
% Optional Parameters:
%   vol_mask: 3D volume specifying voxels with trusted orientations
%     if empty or omitted, will use all non-zero voxels in vol_fiber
%     {default = []}
%   sigmavec: 3-member vector of smoothing sigmas for x, y, and z (voxels)
%     {default = [5,5,5]}
%   win: smoothing window for fiber orientations (voxels)
%     {default = 1}
%   init_niter: max number of iterations for inital orientation smoothing
%     (filling in values only for those not in vol_mask)
%     {default = 100}
%   final_niter: number of iterations for final orientation smoothing
%     {default = 1}
%
%
% Created:  03/07/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_fiber_sm = [];
vol_V_sm = [];

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('vol_mask','var'), vol_mask = []; end;
if ~exist('sigmavec','var') | isempty(sigmavec), sigmavec = [5,5,5]; end;
if ~exist('win','var') | isempty(win), win = 1; end;
if ~exist('init_niter','var') | isempty(init_niter), init_niter = 10; end;
if ~exist('final_niter','var') | isempty(final_niter), final_niter = 1; end;
if ~exist('verbose_flag','var') | isempty(verbose_flag), verbose_flag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volsz = size(vol_fiber);
nx = volsz(1);
ny = volsz(2);
nz = volsz(3);
nvox = nx*ny*nz;

% find intersection of fiber and mask
if ~isempty(vol_mask)
  tmp = vol_mask;
else
  tmp = ones(size(vol_fiber));
end;
vol_mask = zeros(size(vol_fiber));
vol_mask(vol_fiber>0 & tmp>0) = 1;

% smooth fiber mask
sx = sigmavec(1);
sy = sigmavec(2);
sz = sigmavec(3);
vol_fiber_sm = mmil_smooth3d(vol_fiber,sx,sy,sz);
vol_fiber_sm = 1.0*(vol_fiber_sm>0.01);

% smooth orientations
i_trust = find(vol_mask>0);
i_unknown = find(vol_fiber_sm>0 & vol_mask<1);
if verbose_flag
  fprintf('%s: num unknown = %d\n',mfilename,length(i_unknown));
end;
vecV = reshape(vol_V,[nvox,3]);
vecV_sm = zeros(size(vecV));
vecV_sm(i_trust,:) = vecV(i_trust,:);
vol_V_sm = reshape(vecV_sm,[nx,ny,nz,3]);
still_unknown = zeros(size(i_unknown));
for i=1:init_niter
  if verbose_flag
    fprintf('%s: smooth step %d\n',mfilename,i);
  end;
  for k=1:length(i_unknown)
    [x,y,z]=ind2sub(size(vol_mask),i_unknown(k));
    x0 = max([1,x-win]);
    x1 = min([nx,x+win]);
    y0 = max([1,y-win]);
    y1 = min([ny,y+win]);
    z0 = max([1,z-win]);
    z1 = min([nz,z+win]);
    tmpV = vol_V_sm(x0:x1,y0:y1,z0:z1,:);
    if max(abs(tmpV(:)))==0
      still_unknown(k) = 1;
      continue;
    else
      still_unknown(k) = 0;
    end;
    tmpV = reshape(tmpV,[length(tmpV(:))/3,3]);
    C = tmpV'*tmpV;
    [U,S] = eig(C);
    lamdas = real(diag(S)');
    [sorted_lamdas,ind]=sort(lamdas,2,'descend');
    v = U(:,ind(1))';
    vol_V_sm(x,y,z,:) = v;
  end;
  num_still_unknown = length(find(still_unknown));
  if verbose_flag
    fprintf('%s: num still unknown = %d\n',mfilename,num_still_unknown);
  end;
  if num_still_unknown == 0, break; end;
end;

i_unknown = find(vol_fiber_sm>0);
tmp_vol_V_sm = vol_V_sm;
for i=1:final_niter
  if verbose_flag
    fprintf('%s: final smooth step %d\n',mfilename,i);
  end;
  for k=1:length(i_unknown)
    [x,y,z]=ind2sub(size(vol_mask),i_unknown(k));
    x0 = max([1,x-win]);
    x1 = min([nx,x+win]);
    y0 = max([1,y-win]);
    y1 = min([ny,y+win]);
    z0 = max([1,z-win]);
    z1 = min([nz,z+win]);
    tmpV = vol_V_sm(x0:x1,y0:y1,z0:z1,:);
    if max(abs(tmpV(:)))==0
      continue;
    end;
    tmpV = reshape(tmpV,[length(tmpV(:))/3,3]);
    C = tmpV'*tmpV;
    [U,S] = eig(C);
    lamdas = real(diag(S)');
    [sorted_lamdas,ind]=sort(lamdas,2,'descend');
    v = U(:,ind(1))';
    tmp_vol_V_sm(x,y,z,:) = v;
  end;
end;
vol_V_sm = tmp_vol_V_sm;

vol_fiber_sm(vol_fiber>0) = 2.0;
vol_fiber_sm(vol_mask>0) = 3.0;

