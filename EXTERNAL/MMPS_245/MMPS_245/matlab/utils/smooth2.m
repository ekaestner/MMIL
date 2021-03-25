function im_out = smooth2(im_in,span1,span2)

im_out = im_in;

for i = 1:size(im_in,1)
  im_out(i,:) = smooth(squeeze(im_out(i,:)),span1,'loess');
end
for i = 1:size(im_in,2)
  im_out(:,i) = smooth(squeeze(im_out(:,i)),span2,'loess');
end
