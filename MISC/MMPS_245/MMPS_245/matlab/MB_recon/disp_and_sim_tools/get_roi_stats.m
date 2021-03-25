function s = get_roi_stats(x, roi, prc)
%
% function s = get_roi_stats(x, [roi=ones(nx, ny, nec, nsl*nz)], [prc=[2, 5, 10, 20, 50, 80, 90, 95]])
%
% Get statistical performance of varialbe x in the ROI.
%
% Inputs
%   x   - Maps, e.g. G-factor maps or retained SNR maps. Dim: [X(=nx), Y(=ny), Echo(=nec), Slice->SimultaneousSlice Z(=nsl*nz), UndersamplePattern(=nmsk)]
%   roi - Binary ROI masks. Dim: [X(=nx), Y(=ny), Echo(=nec), Slice->Simultaneousslice Z(=nsl*nz)]
%   prc - A vector specifying the percentiles of x to calculate in the ROI.
%
% Output
%   s   - A structure for the statistical performance.
%         All fields have dimension: [Echo(=nec), Slice->Simultaneousslice Z(=nsl*nz), UndersamplePattern(=nmsk)]
%         Fields include:
%           s.min   - minimum
%           s.max   - maximum
%           s.mean  - mean
%           s.prcDD - DD percentile, e.g. s.prc10 is 10 percentile; s.prc50 is 50 percentile, i.e. median.
%
% (c) Kangrong Zhu  Stanford University     Jan 2015

if ~exist('prc', 'var') || isempty(prc)
    prc = [2, 5, 10, 20, 50, 80, 90, 95];
end
np = length(prc);

[nx, ny, nec, nsl_tot, nmsk] = size(x);
sz = [nec, nsl_tot, nmsk];

s.min  = zeros(sz);
s.max  = zeros(sz);
s.mean = zeros(sz);
for idx = 1 : np
    cmd = ['s.prc' num2str(prc(idx)) ' = zeros(sz);'];
    eval(cmd);
end

for ec = 1 : nec
    for sl = 1 : nsl_tot
        im_roi = reshape( roi(:, :, ec, sl), [nx*ny, 1]);
        
        for msk_idx = 1 : nmsk
            vals = reshape( x(:, :, ec, sl, msk_idx), [nx*ny, 1]);
            vals = vals(im_roi == 1);
            
            if ~isempty(vals)
                s.min(ec, sl, msk_idx)  = min(vals);
                s.max(ec, sl, msk_idx)  = max(vals);
                s.mean(ec, sl, msk_idx) = mean(vals);
                
                y = prctile(vals, prc);
                for idx = 1 : np
                    cmd = ['s.prc' num2str(prc(idx)) '(ec, sl, msk_idx) = y(' num2str(idx) ');'];
                    eval(cmd);
                end
            end
        end
    end
end

return
