function d = cohensd(vals1,vals2)

ivec1 = find(isfinite(sum(vals1,1)));
ivec2 = find(isfinite(sum(vals2,1)));
vals1 = vals1(:,ivec1);
vals2 = vals2(:,ivec2);
mu1 = mean(vals1,2);
mu2 = mean(vals2,2);
std1 = std(vals1,0,2);
std2 = std(vals2,0,2);
%stdall = sqrt((n1*std1.^2+n2*std2.^2)/(n1+n2));
stdall = sqrt((std1.^2+std2.^2)/2);
d = (mu1-mu2)./stdall;
