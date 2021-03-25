function [sf] = compute_scalefact_hist1d(valvec,hc,binvals)

midvals = binvals+0.5*(binvals(2)-binvals(1));
meanval_target = sum(hc(:).*midvals(:))/sum(hc);
meanval = mean(valvec);
sf0 = meanval_target/meanval;

nbins = length(binvals);
min1 = min(binvals); 
max1 = max(binvals);
hc_logp = log((hc+1)/sum(hc(:)));
sfvec = sf0*(1+[-0.30:0.02:0.30]');
nsf = length(sfvec);
costvec = zeros(nsf,1);
for sfi = 1:nsf
  sf = sfvec(sfi);
  binvec = 1+floor((sf*valvec-min1)/(max1-min1)*nbins);
  binvec(find(binvec<1)) = 1;
  binvec(find(binvec>nbins)) = nbins;
  ind = binvec;
  cost = -sum(hc_logp(ind))/length(ind);
  costvec(sfi) = cost;
%  fprintf(1,'%f : %f\n',sf,cost);
end
%figure; plot(sfvec,costvec);
[minval,mini] = min(costvec);
mini = max(2,min(length(costvec)-1,mini));
ivec = mini+[-1:1];
s = sfvec(ivec);
M = [ones(size(s)) s s.^2];
beta = pinv(M)*costvec(ivec);
sf = -beta(2)/(2*beta(3));
