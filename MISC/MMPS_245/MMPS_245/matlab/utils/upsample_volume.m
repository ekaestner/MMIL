function vol_us = upsample_volume(vol,ds1,ds2,ds3)

tic
nvox = prod(size(vol));
dim1 = size(vol,1);
dim2 = size(vol,2);
dim3 = size(vol,3);
indvec = [1:nvox]';
[i1,i2,i3] = ind2sub(size(vol),indvec);
vol_us = zeros(dim1*ds1,dim2*ds2,dim3*ds3,'single');
vol_us(1+([1:dim1]-1)*ds1,1+([1:dim2]-1)*ds2,1+([1:dim3]-1)*ds3) = vol;
for i1 = 1:(dim1)
  for d1 = 1:(ds1-1)
    ii1 = 1+(i1-1)*ds1;
    if i1<dim1
      f = d1/ds1;
      vol_us(ii1+d1,:,:) = vol_us(ii1,:,:)*(1-f)+vol_us(ii1+ds1,:,:)*f;
    else
      vol_us(ii1+d1,:,:) = vol_us(ii1,:,:);
    end
  end
end
for i2 = 1:(dim2)
  for d2 = 1:(ds2-1)
    ii2 = 1+(i2-1)*ds2;
    if i2<dim2
      f = d2/ds2;
      vol_us(:,ii2+d2,:) = vol_us(:,ii2,:)*(1-f)+vol_us(:,ii2+ds2,:)*f;
    else
      vol_us(:,ii2+d2,:) = vol_us(:,ii2,:);
    end
  end
end
for i3 = 1:(dim3)
  for d3 = 1:(ds3-1)
    ii3 = 1+(i3-1)*ds3;
    if i3<dim3
      f = d3/ds3;
      vol_us(:,:,ii3+d3) = vol_us(:,:,ii3)*(1-f)+vol_us(:,:,ii3+ds3)*f;
    else
      vol_us(:,:,ii3+d3) = vol_us(:,:,ii3);
    end
  end
end
toc

% crange = [0.7 1.3];
% figure(511);
% subplot(2,2,1); imagesc(abs(squeeze(vol(floor(size(vol,1)/2)+1,:,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,2); imagesc(abs(squeeze(vol(:,floor(size(vol,2)/2)+1,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,3); imagesc(abs(squeeze(vol(:,:,floor(size(vol,3)/2)+1))).',crange); colormap(hot); axis equal; axis image; colorbar;
% 
% figure(512);
% subplot(2,2,1); imagesc(abs(squeeze(vol_us(floor(size(vol_us,1)/2)+1,:,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,2); imagesc(abs(squeeze(vol_us(:,floor(size(vol_us,2)/2)+1,:))).',crange); colormap(hot); axis equal; axis image; colorbar;
% subplot(2,2,3); imagesc(abs(squeeze(vol_us(:,:,floor(size(vol_us,3)/2)+1))).',crange); colormap(hot); axis equal; axis image; colorbar;
