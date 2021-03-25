function vol_sm = SparseSmoothScalarVol(vol,volwt,niter)

lambda = 0.3;
lambda1 = 0.2;
lambda2 = 0.5;
momentum = 0.8;

dims = size(vol);
dims1 = dims(1);
dims2 = dims(2);
if(length(dims)==3)
  dims3 = dims(3);
else
  dims3 = 1; % single slice
end;

vol_sm = vol.*(volwt~=0);

options = optimset('GradObj','on');

delvol = zeros(size(vol),'single');

for iter = 1:niter

  tmp1 = lambda1*(vol_sm - vol).*volwt + lambda2*((cat(1,zeros(1,dims2,dims3,'single'),diff(vol_sm,1,1))+...
                                                   cat(2,zeros(dims1,1,dims3,'single'),diff(vol_sm,1,2))+...
                                                   cat(3,zeros(dims1,dims2,1,'single'),diff(vol_sm,1,3)))-...
                                                  (cat(1,diff(vol_sm,1,1),zeros(1,dims2,dims3,'single'))+...
                                                   cat(2,diff(vol_sm,1,2),zeros(dims1,1,dims3,'single'))+...
                                                   cat(3,diff(vol_sm,1,3),zeros(dims1,dims2,1,'single'))));

  diffvol = (vol - vol_sm).*volwt;
  lapvol = (cat(1,zeros(1,dims2,dims3,'single'),diff(vol_sm,1,1)).^2+...
            cat(2,zeros(dims1,1,dims3,'single'),diff(vol_sm,1,2)).^2+...
            cat(3,zeros(dims1,dims2,1,'single'),diff(vol_sm,1,3)).^2);

  c1 = floor(dims1/2);
  c2 = floor(dims2/2);
  c3 = floor(dims3/2);

  cost1 = (mean(diffvol(:).^2));
  cost2 = (mean(lapvol(:)));
  cost = (lambda1*cost1+lambda2*cost2)/(lambda1+lambda2);
%  fprintf('iter = %d: cost1 = %f (%f) cost2 = %f (%f) cost = %f\n',iter,cost1,lambda1*cost1,cost2,lambda2*cost2,cost);

  delvol = -lambda*tmp1 + momentum*delvol;
  delvol = max(-0.25,min(0.25,delvol));
  vol_sm = vol_sm + delvol;
end

