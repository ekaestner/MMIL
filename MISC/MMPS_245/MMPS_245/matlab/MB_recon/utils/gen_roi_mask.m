function msk = gen_roi_mask(im, mag_thresh)
%
% function msk = gen_roi_mask(im, [mag_thresh=0.05])
%
% Generate binary masks for image ROIs, basing on the input magnitude threshold.
%
% Inputs
%   im     - Images to generate masks from. Each 2D image will be processed
%            independently and will correspond to one ROI mask.
%            Dim: [image_dimension1, image_dimension2, ...].
%   mag_thresh - Magnitude mag_threshold as percentage of the maximum magnitude in a 2D image.
%
% Output
%   msk    - Image masks, having the same matrix size as the input 'im'.
%
% (c) Kangrong Zhu      Stanford University     May 2014

if ~exist('mag_thresh', 'var') || isempty(mag_thresh)
    mag_thresh = 0.05;
end

sz = size(im);
nim = prod(sz(3 : end));
msk = zeros(sz);
erode_dilate_disk_sz = 3;
isolated_limit = 150;
for im_idx = 1 : nim
    im_2d = im(:, :, im_idx);
    m = (abs(im_2d) > mag_thresh*max(abs(im_2d(:))));
    
    % Erode and dilate: works well for phantom
% % %     se_erode = strel('disk', erode_dilate_disk_sz);
% % %     se_dilate = strel('disk', erode_dilate_disk_sz);
% % %     m = imerode(m, se_erode);
% % %     m = imdilate(m, se_dilate);
    
    % 2D filtering: works well for brain
    if size(m, 1) <= 64
        hsz = [3, 3];
    else
        hsz = [5, 5];
    end
    h = fspecial('gaussian', hsz, 0.2);
    m = filter2(h, m);
    
    m = (m > 0.5);
    m = ~elimiso(~m, isolated_limit);
    m = elimiso(m, isolated_limit);
    msk(:, :, im_idx) = m;
end

return