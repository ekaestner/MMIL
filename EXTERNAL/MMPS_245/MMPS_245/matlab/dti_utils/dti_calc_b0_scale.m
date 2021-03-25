function scalefacts = dti_calc_b0_scale(volb0,pthresh_min,pthresh_max)
%function scalefacts = dti_calc_b0_scale(volb0,[pthresh_min],[pthresh_max])
%
%  A histogram of the each frame is calculated and pthresh_min
%    and pthresh_max are used to select threshold values
%  Values between the min and max thresholds are averaged.
%  The average values for each frame are used to calculate scaling factors
%    to normalize each volume to the common mean
%
% Required Parameters:
%   volb0: 4D volume containing multiple b=0 volumes
%
% Optional Parameters:
%   pthresh_min: minimum cumulative probability threshold
%     {default=0.75}
%   pthresh_max: maximum cumulative probability threshold
%     {default=0.95}
%
% Output:
%    scalefacts : vector of scaling factors for each frame in volb0
%
% Created:  01/23/08 Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

plotflag = 0;
smf = 10^-5;
nbins = 1000;
if ~exist('pthresh_min','var') | isempty(pthresh_min), pthresh_min=0.75; end;
if ~exist('pthresh_max','var') | isempty(pthresh_max), pthresh_max=0.95; end;

scalefacts = [];

% check input dimensions
nx = size(volb0,1);
ny = size(volb0,2);
nz = size(volb0,3);
nf = size(volb0,4);

for f=1:nf
  vol = squeeze(volb0(:,:,:,f));
  [N,X] = hist(vol(:),nbins);
  P = cumsum(N);
  Pnorm = P/P(end);

  [tmp,ind_min] = min(abs(Pnorm-pthresh_min));
  [tmp,ind_max] = min(abs(Pnorm-pthresh_max));
  thresh_min = X(ind_min);
  thresh_max = X(ind_max);

  if plotflag
    figure(1);
    clf;
    subplot(2,1,1);
    plot(X,N/P(end));
    subplot(2,1,2);
    plot(X,Pnorm);

    vol_thresh = vol;
    vol_thresh(vol<thresh_min | vol>thresh_max)=0;
    tmp0 = squeeze(vol(:,:,round(nz/2)))';
    tmp1 = squeeze(vol_thresh(:,:,round(nz/2)))';
    figure(2); imagesc(tmp0,[thresh_min,thresh_max]); colorbar;
    figure(3); imagesc(tmp1,[thresh_min,thresh_max]); colorbar;
    
    fprintf('%s: press a key to continue...\n',mfilename);
    pause
  end;

  vals = vol(vol>thresh_min & vol<thresh_max);
  scalefacts(f) = mean(vals);
end;

scalefacts = mean(scalefacts)./scalefacts;

return;

