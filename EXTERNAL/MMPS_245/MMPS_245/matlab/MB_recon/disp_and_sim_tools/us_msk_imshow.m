function us_msk_imshow(us_msk, ny, dispfig, axison, kydir, nomegaz)
%
% function us_msk_imshow(us_msk, ny, [dispfig=1], [axison=false], [kydir=BOTTOM_UP], [nomegaz=[]])
%
% Display undersample masks on the omegay-omegaz plane.
%
% Inputs
%   us_msk  - Ky-omegaz undersample mask with fields 'ky', 'kz' and 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   ny      - Full prescribed matrix size in ky.
%   dispfig - True: display results in figure window (1000+dispfig).
%   axison  - True: display axis on the figures.
%   kydir   - Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
%   nomegaz - Number of samples along omegaz axis.
%             If not empty, will define an omegaz grid with 'nomegaz' lines.
%             If empty, will use default number of omegaz lines (number of unique omegaz values for uniform sampling along omegaz; ny for nonuniform sampling along omegaz) to define the grid.
%
% (c) Kangrong Zhu      Stanford University     June 2014

BOTTOM_UP  = 2;

if ~exist('dispfig', 'var') || isempty(dispfig) || (dispfig <= 0)
    dispfig = 1;
end
if ~exist('axison', 'var') || isempty(axison)
    axison = false;
end
if ~exist('kydir', 'var') || isempty(kydir)
    kydir = BOTTOM_UP;
end
if ~exist('nomegaz', 'var') || isempty(nomegaz)
    nomegaz = [];
end

% Get matrix representation of the undersample mask
mtx = us_msk_get_mtx(us_msk, ny, kydir, nomegaz);

% Get number of columns and rows to display
nsubplots = length(mtx);
[nrow, ncol] = get_montage_nrow_ncol(nsubplots);

% Display
figure(1000+dispfig);
for msk_idx = 1 : nsubplots
    subplot(nrow, ncol, msk_idx);
    imagesc(mtx(msk_idx).omegaz, mtx(msk_idx).omegay, mtx(msk_idx).mtx);
    
    nxticks = 5;
    if length(mtx(msk_idx).omegaz) <= nxticks
        xticks = sort(mtx(msk_idx).omegaz);
    else
        xticks = linspace(min(mtx(msk_idx).omegaz), max(mtx(msk_idx).omegaz), nxticks);
    end
    set(gca,'XTick', xticks);
    nyticks = 10;
    yticks = linspace(min(mtx(msk_idx).omegay), max(mtx(msk_idx).omegay), nyticks);
    set(gca,'YTick', yticks);
    
    % Adjust colorbar for better contrast
    nclr = 256;
    lowval_cut_ratio = 0.15;  % This percent of the low values in the colorbar will be cut.
    highval_cut_ratio = 0.1;  % This percent of the high values in the colorbar will be cut.
    nclr_cut_lowval = round(lowval_cut_ratio * nclr);
    nclr_cut_highval = round(highval_cut_ratio * nclr);
    cmap = jet(nclr);
    cmap(1,:) = [0, 0, 0];   % Make lowest value black
    cmap = cat(1, cmap(1,:), cmap(nclr_cut_lowval+1 : nclr-nclr_cut_highval+1, :));
    colormap(cmap);
    if axison
        xlabel('\omega_z (\in [-\pi, \pi])'); ylabel('\omega_y (\in [-\pi, \pi])');
        h = colorbar;
        title(h, 'Echo Index');
    end
end

return