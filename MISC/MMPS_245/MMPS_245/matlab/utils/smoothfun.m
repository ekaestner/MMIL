function [fvals,xbins] = smoothfun(xvec,fvec,nbins)

v0 = min(xvec);
v1 = max(xvec);

fvals = zeros(1,nbins);
for bin = 1:nbins
  vstep = (v1-v0)/nbins;
  v = v0+(bin-1)*vstep;
  ind = find((xvec>=v)&(xvec<=v+vstep));
  if length(ind)>0
    fvals(bin) = mean(fvec(ind));
  end
  xbins(bin) = v+vstep/2;
end
