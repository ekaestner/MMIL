function vol_ctx = ctx_mgh2ctx(vol,M)
%function vol_ctx = ctx_mgh2ctx(vol,M)

[width,height,depth] = size(vol) ;
[xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s] = mat2mgh(M,width,height,depth);

vol_ctx.imgs = double(vol);
vol_ctx.Mvxl2lph = M_RAS_TO_LPH * M;
vol_ctx.dimc = size(vol,2);
vol_ctx.dimr = size(vol,1);
vol_ctx.dimd = size(vol,3);
vol_ctx.vx = xsize;
vol_ctx.vy = ysize;
vol_ctx.vz = zsize;
vol_ctx.lphcent=[-c_r;-c_a;c_s];
vol_ctx.DirCol=[-x_r;-x_a;x_s];
vol_ctx.DirRow=[-y_r;-y_a;y_s];
vol_ctx.DirDep=[-z_r;-z_a;z_s];
