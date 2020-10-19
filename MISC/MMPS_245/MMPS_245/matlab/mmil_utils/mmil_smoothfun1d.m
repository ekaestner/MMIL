function [xvals,yvals,wtstruct,stats] = mmil_smoothfun1d(xvec,yvec,xvals,lam)
%function [xvals,yvals,wtstruct,stats] = mmil_smoothfun1d(xvec,yvec,xvals,lam)
%
% Created:  01/08/15 by Don Hagler
% Last Mod: 01/08/15 by Don Hagler
%

% based on code by Anders Dale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xvec = mmil_colvec(xvec);
yvec = mmil_colvec(yvec);
ivec = find(isfinite(xvec) & isfinite(yvec));
xvec = xvec(ivec);
yvec = yvec(ivec);
xvals = mmil_colvec(xvals);
xvali = interp1(xvals,[1:length(xvals)],xvec);
xvali(find(xvec<=xvals(1))) = 1+0.001;
xvali(find(xvec>=xvals(end))) = length(xvals)-0.001;
nobs = length(xvec);
nbins = length(xvals);
S = sparse([[1:nobs]'; [1:nobs]'],[floor(xvali); floor(xvali)+1],[1-(xvali-floor(xvali)); (xvali-floor(xvali))],nobs,nbins,2*nobs);
D = diff(speye(nbins,nbins));
L1 = D'*D;
L1(1,:) = 0;
L1(end,:) = 0;
L2 = L1'*L1;
W = inv(S'*S+lam*L2);
tmp = S'*yvec;
yvals = W*tmp;

if nargout>2
  wtstruct.S = S;
  wtstruct.W = W;
  wtstruct.xvals = xvals;
end

if nargout>3
  yvec_hat = S*W*(S'*yvec);
  m = nobs;
  A = S*W*S';
  I = eye(size(A));
  p_eff = trace(A); 
  stats.AIC = m*log(sum((yvec_hat-yvec).^2)) + 2*p_eff;
  stats.BIC = m*log(sum((yvec_hat-yvec).^2)) + p_eff*log(nobs);
end

