function [vol,volmask] = ComputeICV(vol_loFA,vol_hiFA)

if exist('vol_hiFA','var'), brainmesh = fastSkull(vol_loFA,vol_hiFA); else brainmesh = fastSkull(vol_loFA); end
[tmp, volmask]=getmaskvol(vol_loFA, brainmesh, eye(4));
vol = length(find(volmask.imgs)>0.5)*prod(svd(volmask.Mvxl2lph(1:3,1:3)))/1000;

