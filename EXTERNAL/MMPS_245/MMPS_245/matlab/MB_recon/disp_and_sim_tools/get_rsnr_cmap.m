function [cmap, dispr] = get_rsnr_cmap(maxval)
%
% function [cmap, dispr] = get_rsnr_cmap([maxval=1.25])
%
% Generates gray and jet concatenated colormap similar to the one used in the blipped-CAIPI
% paper(Setsompop K. et al, MRM 2012;67(5):1210-24), for displaying retained SNR maps.
%
% Input
%   maxval - Maximum value in the display range.
%
% Outputs
%   cmap   - The colormap.
%   dispr  - The display range.
%
% (c) Kangrong Zhu  Stanford University     Nov 2014

if ~exist('maxval', 'var') || isempty(maxval)
    maxval = 1.25;
end
ncolors = 256;
dispr = [0, maxval];
cmap1 = jet(round(1/maxval * ncolors));
cmap2 = gray(ncolors - size(cmap1, 1));
cmap = cat(1, cmap1, cmap2);

return
