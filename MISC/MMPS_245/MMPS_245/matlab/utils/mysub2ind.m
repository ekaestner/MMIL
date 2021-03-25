function ind = mysub2ind(siz,svec1,svec2,svec3)
[I1,I2,I3] = meshgrid(svec1,svec2,svec3);
ind = sub2ind(siz,I1(:),I2(:),I3(:));