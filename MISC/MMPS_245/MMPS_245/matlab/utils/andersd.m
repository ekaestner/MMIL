function [dval,ci] = andersd(dvalsA, dvalsB, a, p, fraq)

if ~exist('fraq','var'), fraq = 0.25; end
if ~exist('p','var'),    p    = 0.9;  end
if ~exist('a','var'),    a    = 0.05; end

if isempty(dvalsB), dvalsB = 0; end
ivecA  = isfinite(dvalsA);
ivecB  = isfinite(dvalsB);
dvalsA = dvalsA(ivecA);
dvalsB = dvalsB(ivecB);

n1 = length(dvalsA);
n2 = length(dvalsB);
s1 = std(dvalsA);
s2 = std(dvalsB);
u = s2^2/s1^2;
mu   = mean(dvalsA)-mean(dvalsB);
tval = mu/sqrt(s1^2/n1+s2^2/n2);
dof   = (1/n1+u/n2)^2/(1/(n1^2*(n1-1))+u^2/(n2^2*(n2-1)));
dval = mu/sqrt((s1^2+s2^2)/2);
%dval = 2*tval/sqrt(dof); % This should work, according to uccs
sf = dval/tval;

if nargout>1 & n1>1 & n2>1
  dval0 = nctinv(0.025,dof,tval)*sf;
  dval1 = nctinv(0.975,dof,tval)*sf;
  ci = [dval0 dval1];
end
