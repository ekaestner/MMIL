function area = us_msk_get_area(us_msk, ny, dispfig)
%
% Get the k-space area corresponding to each sample in the undersample masks.
%
% Inputs
%   us_msk  - Ky-omegaz undersample mask with fields 'ky', 'kz' and 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   ny      - Full prescribed matrix size in ky.
%   dispfig - True: display results in figure window (1000+dispfig).
%
% Output
%   area    - Output structure, with fields:
%             'area':    K-space area each sample takes. Dim: [nsamp, nmsk].
%             'var':     Variance of the k-space area for each mask. Dim: [1, nmsk].
%             'maxdiff': Maximum difference between the two k-space areas for each mask. Dim: [1, nmsk].
%
% (c) Kangrong Zhu      Stanford University     June 2014

if ~exist('dispfig', 'var') || isempty(dispfig)
    dispfig = 0;
end

nmsk = length(us_msk);
nsamp = length(us_msk(1).ky);
mtx = us_msk_get_mtx(us_msk, ny);

% Calculate undersample mask by undersample mask
area.area = zeros(nsamp, nmsk);
area.var = zeros(1, nmsk);
area.maxdiff = zeros(1, nmsk);
for midx = 1 : nmsk
    % ky, kz values, in [-pi, pi]
    ky = mtx(midx).omegay(us_msk(midx).ky);
    kz = mtx(midx).omegaz(us_msk(midx).kz);
    
    % To use voronoidens, repeat the sample pattern along ky and kz so that all sampled points have neighbours.
    ky = [ky(:); ky(:); ky(:); ky(:)-2*pi; ky(:)+2*pi; ky(:)-2*pi; ky(:)-2*pi; ky(:)+2*pi; ky(:)+2*pi];
    kz = [kz(:); kz(:)-2*pi; kz(:)+2*pi; kz(:); kz(:); kz(:)-2*pi; kz(:)+2*pi; kz(:)-2*pi; kz(:)+2*pi];
    area_this_msk = voronoidens(ky, kz);
    
    % Output
    area.area(:, midx) = reshape(area_this_msk(1:nsamp), [nsamp, 1]);
    area.var(midx) = var(area.area(:, midx));
    area.maxdiff(midx) = max(area.area(:, midx)) - min(area.area(:, midx));
end

% Display
if dispfig
    minarea = min(area.area(:));
    maxarea = max(area.area(:));
    [nrow, ncol] = get_montage_nrow_ncol(nmsk);
    figure(2000 + dispfig);
    for midx = 1 : nmsk
        subplot(nrow, ncol, midx);
        plot(1 : nsamp, area.area(:, midx), 'x--', 'MarkerSize', 5);
        ylim([minarea, maxarea]);
        xlabel('Sample');
        ylabel('K-space Area');
    end
end

return