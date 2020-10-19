function vol_sm = smooth_sparse_vol(vol,wvol,p,ds1,ds2,ds3)

tic
vol = vol(1:ds1:end,1:ds2:end,1:ds3:end);
wvol = wvol(1:ds1:end,1:ds2:end,1:ds3:end);
nn = 6;
nvox = prod(size(vol));
dim1 = size(vol,1);
dim2 = size(vol,2);
dim3 = size(vol,3);
indvec = [1:nvox]';
[i1,i2,i3] = ind2sub(size(vol),indvec);
neigbormat = [i1-1 i2 i3 indvec; i1+1 i2 i3 indvec; i1 i2-1 i3 indvec; i1 i2+1 i3 indvec; i1 i2 i3-1 indvec; i1 i2 i3+1 indvec];
j1 = neigbormat(:,1);
j2 = neigbormat(:,2);
j3 = neigbormat(:,3);
j4 = neigbormat(:,4);
okind = find(~(j1<1|j1>dim1|j2<1|j2>dim2|j3<1|j3>dim3));
j1 = j1(okind);
j2 = j2(okind);
j3 = j3(okind);
j4 = j4(okind);
j5 = sub2ind(size(vol),j1,j2,j3);
[B,I,J] = unique(sort(j4));
nnvec = diff([0;I]);
M = sparse([indvec;j4],[indvec;j5],[p*double(wvol(indvec))+(1-p);-(1-p)*ones(size(j4))./nnvec(j4)],nvox,nvox);
b = double(p*wvol(indvec).*vol(indvec));
clear('B','I','J','i1','i2','i3','indvec','j1','j2','j3','j4','j5','neigbormat','nnvec','okind');
volvec_sm = M\b;
tmp = single(reshape(volvec_sm,size(vol)));
toc
vol_sm = upsample_volume(tmp,ds1,ds2,ds3);      

crange = [0.7 1.3];
figure(311);
subplot(2,2,1); imagesc(abs(squeeze(vol(floor(size(vol,1)/2)+1,:,:))),crange); colormap(hot); axis equal; axis image; colorbar;
subplot(2,2,2); imagesc(abs(squeeze(vol(:,floor(size(vol,2)/2)+1,:))),crange); colormap(hot); axis equal; axis image; colorbar;
subplot(2,2,3); imagesc(abs(squeeze(vol(end:-1:1,:,floor(size(vol,3)/2)+1))).',crange); colormap(hot); axis equal; axis image; colorbar;

figure(312);
subplot(2,2,1); imagesc(abs(squeeze(vol_sm(floor(size(vol_sm,1)/2)+1,:,:))),crange); colormap(hot); axis equal; axis image; colorbar;
subplot(2,2,2); imagesc(abs(squeeze(vol_sm(:,floor(size(vol_sm,2)/2)+1,:))),crange); colormap(hot); axis equal; axis image; colorbar;
subplot(2,2,3); imagesc(abs(squeeze(vol_sm(end:-1:1,:,floor(size(vol_sm,3)/2)+1))).',crange); colormap(hot); axis equal; axis image; colorbar;

% figure(312);
% subplot(2,2,1); imagesc(abs(squeeze(tmp(floor(size(tmp,1)/2)+1,:,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,2); imagesc(abs(squeeze(tmp(:,floor(size(tmp,2)/2)+1,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,3); imagesc(abs(squeeze(tmp(:,:,floor(size(tmp,3)/2)+1))).',crange); colormap(hot); axis equal; axis image; colorbar;

% figure(315);
% subplot(2,2,1); imagesc(abs(squeeze(wvol(floor(size(wvol,1)/2)+1,:,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,2); imagesc(abs(squeeze(wvol(:,floor(size(wvol,2)/2)+1,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,3); imagesc(abs(squeeze(wvol(:,:,floor(size(wvol,3)/2)+1))).',crange); colormap(hot); axis equal; axis image; colorbar;
