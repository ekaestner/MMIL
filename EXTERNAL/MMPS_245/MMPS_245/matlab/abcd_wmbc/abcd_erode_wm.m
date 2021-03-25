function vol_wm = abcd_erode_wm(vol_wm)
%function vol_wm = abcd_erode_wm(vol_wm)
%
% Created:  11/06/17 by Don Hagler
% Last Mod: 11/06/17 by Don Hagler
%

parms = mmil_args2parms({},{...
  'nvoxels',1,[1:100],...
  'sigma',1,[1e-2,10],...
  'thresh_init',1e-5,[0,Inf],...
  'thresh',0.99,[0.1,0.999],...
});

vol = 1.0 * (vol_wm.imgs>1-parms.thresh_init);
for j=1:parms.nvoxels
  vol = mmil_smooth3d(vol,parms.sigma,parms.sigma,parms.sigma);
  vol = 1.0*(vol>=parms.thresh);
end;
vol_wm.imgs = vol;

return;
