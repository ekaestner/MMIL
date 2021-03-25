function r_new=mri2neuromag(r,r_LPA,r_RPA,r_NAS);
% function r_new=mri2neuromag(r,r_LPA,r_RPA,r_NAS);
% transfer the coordinate from mri 2 neuromag coordinates
% r: N by 3 locations in mri system
% r_LPA: 1 by 3 left PA location in mri system
% r_RPA: 1 by 3 right PA location in mri system
% r_NAS: 1 by 3 nasion location in mri system
% r_new: N by 3 locations in neuromag system
% 
% (c) Mingxiong Huang, Ph.D.

x1=r_LPA(1);
y1=r_LPA(2);
z1=r_LPA(3);

x2=r_RPA(1);
y2=r_RPA(2);
z2=r_RPA(3);

x3=r_NAS(1);
y3=r_NAS(2);
z3=r_NAS(3);

A=(z1-z3)*(z2-z1)*(x2-x1);
B=(y1-y3)*(y2-y1)*(x2-x1);
C=x1*((z2-z1)*(z2-z1)+(y2-y1)*(y2-y1))+x3*(x2-x1)*(x2-x1);
D=(z2-z1)*(z2-z1)+(y2-y1)*(y2-y1)+(x2-x1)*(x2-x1);

x0=(C-A-B)/D;
y0=y1+(x0-x1)*(y2-y1)/(x2-x1);
z0=z1+(x0-x1)*(z2-z1)/(x2-x1);


r_cen=[r(:,1)-x0, r(:,2)-y0, r(:,3)-z0];

newx=[x2-x0 y2-y0 z2-z0];
newy=[x3-x0 y3-y0 z3-z0];
newz=crossprod(newx,newy);

newx=newx/sqrt(newx(1)^2+newx(2)^2+newx(3)^2);
newy=newy/sqrt(newy(1)^2+newy(2)^2+newy(3)^2);
newz=newz/sqrt(newz(1)^2+newz(2)^2+newz(3)^2);

r_new=r_cen*[newx' newy' newz'];
