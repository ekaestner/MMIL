function [hcnt,binvals1,binvals2] = hist2d(valvec1,valvec2,nbins1,nbins2)

if max(size(nbins1))>1
  binvals1 = nbins1;
  nbins1 = length(binvals1);
end
if max(size(nbins2))>1
  binvals2 = nbins2;
  nbins2 = length(binvals2);
end
if ~exist('binvals1')
  min1 = min(valvec1); max1 = max(valvec1)*1.000001;
  binvals1 = min1+[0:nbins1-1]*(max1-min1)/nbins1;
end
if ~exist('binvals2')
  min2 = min(valvec2); max2 = max(valvec2)*1.000001;
  binvals2 = min2+[0:nbins2-1]*(max2-min2)/nbins2;
end
min1 = min(binvals1); 
max1 = max(binvals1);
min2 = min(binvals2); 
max2 = max(binvals2);
binvec1 = 1+floor((valvec1-min1)/(max1-min1)*nbins1);
binvec2 = 1+floor((valvec2-min2)/(max2-min2)*nbins2);
binvec1(find(binvec1<1)) = 1;
binvec1(find(binvec1>nbins1)) = nbins1;
binvec2(find(binvec2<1)) = 1;
binvec2(find(binvec2>nbins2)) = nbins2;
ind = sub2ind([nbins1,nbins2],binvec1,binvec2);
[hcntvec] = histc(ind,1:nbins1*nbins2);
hcnt = reshape(hcntvec,[nbins1,nbins2]);
%figure; imagesc(log10((hcnt+0.5)/length(valvec1))); colormap(hot); colorbar;
