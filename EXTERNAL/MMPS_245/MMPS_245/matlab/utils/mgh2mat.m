function M = mgh2mat(xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s,width,height,depth);

% function M = mgh2mat(xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s,width,height,depth);

M = eye(4);

M(1,1) = xsize * x_r ;
M(1,2) = ysize * y_r ;
M(1,3) = zsize * z_r ;

M(2,1) = xsize * x_a ;
M(2,2) = ysize * y_a ;
M(2,3) = zsize * z_a ;

M(3,1) = xsize * x_s ;
M(3,2) = ysize * y_s ;
M(3,3) = zsize * z_s ;

ci = (width)/2 ; cj = (height)/2 ; ck = (depth)/2 ;
M(1,4) = c_r - (M(1,1)*ci + M(1,2)*cj + M(1,3)*ck) ;
M(2,4) = c_a - (M(2,1)*ci + M(2,2)*cj + M(2,3)*ck);
M(3,4) = c_s - (M(3,1)*ci + M(3,2)*cj + M(3,3)*ck);
M(4,4) = 1 ;

M=c2matlab(M);

return
