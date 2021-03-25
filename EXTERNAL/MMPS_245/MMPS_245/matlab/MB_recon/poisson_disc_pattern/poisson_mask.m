function [us_msk, us_msk_cand] = poisson_mask(nky, nkz, inplane_R, nz, ker_sup_dky, fov, kydir, non_repeat_dims, nmsk_tgt, debug)
%
% function [us_msk, us_msk_cand] = poisson_mask(nky, nkz, inplane_R, [nz=nkz], [ker_sup_dky=0.0259034741328851], [fov=200], [kydir=BOTTOM_UP], [non_repeat_dims=1], [nmsk_tgt=20], [debug=false])
%
% Generates a 2D Poisson-disc-like ky-kz undersampling mask.
%
% Inputs
%   nky              - Mask size in ky.
%   nkz              - Mask size in kz.
%   inplane_R        - Inplane acceleration factor.
%   nz               - Number of simultaneous slices. Used when calculating the point spread function of the undersampling mask.
%   ker_sup_dky      - Sensitivity map kernel support in ky, in cycles/mm. The default value is estimated for the NOVA 32 channel head coil.
%   fov              - Field of view, in mm.
%   kydir            - Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
%   non_repeat_dims  - The dimensions along which the sampled positions don't repeat. Value is 1 or [1,2].
%   nmsk_tgt         - Target number of masks to generate. Will generate 'nmsk_tgt' undersampling masks then choose the one with the best point spread function. 
%   debug            - True: print debugging messages.
%
% Outputs
%   us_msk           - Poisson-disc-like undersampling mask on the ky-kz plane, with fields 'ky', 'omegaz', 'kz'. See get_ky_omegaz_us_msk.m for details.
%   us_msk_cand      - 'nmsk_tgt' undersampling mask candidates generated in all the iterations. The final output us_msk is the one with the best point spread function in these 'nmsk_tgt' candidates.
%
% (c) Kangrong Zhu      Stanford University     2014

dispfig = 0; % >0: Display the mask generating process on a figure; =0: Don't display

if ~exist('nz', 'var') || isempty(nz)
    nz = nkz;
end
if ~exist('ker_sup_dky', 'var') || isempty(ker_sup_dky)
    ker_sup_dky = 0.0259034741328851;
end
if ~exist('non_repeat_dims', 'var') || isempty(non_repeat_dims)
    non_repeat_dims = 1;
end
if ~exist('nmsk_tgt', 'var') || isempty(nmsk_tgt)
    nmsk_tgt = 20;
end
if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end
if (inplane_R > 1) && (mod(nky, inplane_R) ~= 0)
    nky = round(nky/inplane_R) * inplane_R;
    if debug
        fprintf('   nky isn''t multiple of inplane_R. Set nky to %d.\n', nky);
    end
end

% -- For filling in gap in the ky dimension
delta_ky = 1 / fov;                                   % Ky sampling interval, in cycles/mm
ker_sup_nky = ker_sup_dky / delta_ky;                 % Kernel support along ky, in number of ky lines
ky_gap_tol = 1.8;                                     % The function will try to fill in gaps larger than ky_gap_tol*ker_sup_nky along the ky dimension.

% -- Probability density function
R = inplane_R * nkz;
pdf = (1/R) * ones(nky, nkz);

% -- Ky indices
p.TOP_DOWN   = 0;
p.CENTER_OUT = 1;
p.BOTTOM_UP  = 2;
if ~exist('kydir', 'var') || isempty(kydir)
    kydir = p.BOTTOM_UP;
end
p.kydir = kydir;
p.inplane_R = inplane_R;
ky_indices = 1 : 1 : nky;
if p.inplane_R > 1
    ky_indices = bypass_lower_ileaves(ky_indices, p); % Similar to the MICA case in function get_ky_omegaz_us_msk.
end
nsamp = length(ky_indices);                           % Total number of samples

% -- Normalize coordinates so that the distance between two points can be calculated more accurately.
% Longer dimension has range [0,1], index starts from 0.
normalize_factor = max(nky, nkz) - 1;
max_ky = (nky-1) / normalize_factor;
max_kz = (nkz-1) / normalize_factor;
ky_coord_all = (ky_indices - 1) / normalize_factor;   % Ky coordinates of all samples

% -- Distance tolerance between two sample points
area_tot = max_ky * max_kz;                           % Total area
area_per_samp = area_tot / nsamp;                     % Area each sample takes
asp_ratio = max_ky / max_kz;
kz_tol = sqrt(area_per_samp / asp_ratio / inplane_R); % Assume each sample takes a rectangular area. area_per_samp = (asp_ratio*kz_tol * inplane_R) * kz_tol 
ky_tol = asp_ratio * kz_tol * inplane_R;
dist_tol = [ky_tol, kz_tol];
dist_tol = dist_tol * 2/3;                            % Empirical value. Roughly set from the fact that if new points were generated around a sample point, they are generated in the range [dist_tol, 2*dist_tol], so 1.5*dist_tol roughly corresponds to the area each sample takes

% -- Generate 'nmsk_tgt' different undersampling masks as candidate masks
if dispfig
    figure(dispfig);
end
if debug
    fprintf('   Generating Poisson disc underample pattern...\n');
end

us_msk_cand(nmsk_tgt) = struct();
n_success_iters = 0;
while n_success_iters < nmsk_tgt
    if dispfig
        clf(dispfig);
    end
    
    sample_coord = zeros(nsamp, 2);                   % Coordinates for the sampled points, dim: [points, coordinates]. Preallocate for speed.
    nsamp_found = 0;                                  % Number of samples found
    sample_idx_looked = [];                           % Sample indices (in [1, nky*nkz]) that have been looked at

    while nsamp_found < nsamp                         % There are still ky candidates who don't have a sample yet
        
        if length(unique(sample_idx_looked)) == nky*nkz
            break;                                    % Already looked at all possible sample positions but couldn't find a sample for each of the ky candidates. So break.
        end
        
        % Randomly pick points according to the probability density function
        tmp = (rand(nky, nkz) < pdf);
        idx = find(tmp == 1);
        idx = idx(~ismember(idx, sample_idx_looked)); % Only look at indices that hasn't been looked at
        sample_idx_looked = [sample_idx_looked; idx];
        n_new_points = length(idx);
        idx = idx(randperm(n_new_points));            % Randomly permutes the indices, otherwise will look at indices on the first line and the first line will just become uniformly undersampled
        ky_coord_cand = (mod(idx-1, nky)) / normalize_factor;
        kz_coord_cand = (ceil(idx/nky)-1) / normalize_factor;
        if dispfig
            plot(kz_coord_cand, ky_coord_cand, '.');
            axis([0, max_kz, 0, max_ky]); hold on; drawnow;
        end
        
        % Check the new points and determine whether to include them or not
        for idx = 1 : n_new_points
            new_point = [ky_coord_cand(idx), kz_coord_cand(idx)];
            if ismember(new_point(1), ky_coord_all) % new_point is among the desired ky locations
                if ~inNeighbourhood(sample_coord(1:nsamp_found, :), new_point, dist_tol) % No other points exist in the new_point's neighbourhood
                    if ~acquiredBeforeAlongNonRepeatableDims(sample_coord(1:nsamp_found, :), new_point, non_repeat_dims, [nky, nkz], [max_ky, max_kz]) % The new_point's coordinates haven't been acquired before along the non repeatable dimensions
                    
                        % Update containers
                        nsamp_found = nsamp_found + 1;
                        sample_coord(nsamp_found, :) = new_point;
                        
                        if dispfig
                            plot(new_point(2), new_point(1), 'ro');
                            axis([0, max_kz, 0, max_ky]); hold on; drawnow;
                        end
                    end
                end
            end
        end
    end
    
    if debug
        fprintf('    dist_tol %.5f %.5f.\t', dist_tol(1), dist_tol(2));
    end
    
    % If successfully generated a mask, increase distance tolerance; otherwise, decrease distance tolerance.
    % In this way, the generated masks are always generated using a tolerance just small enough for generating a mask.
    if nsamp_found == nsamp % Successfully generated one undersample mask
        n_success_iters = n_success_iters + 1;
        
        % Calculate the indices of the sampled points on the mask
        sample_coord = sample_coord(1:nsamp, :);
        sample_grid_idx = coord2grid(sample_coord, [nky, nkz], [max_ky, max_kz]);
        
        us_msk_cand(n_success_iters).ky = sample_grid_idx(:, 1);
        us_msk_cand(n_success_iters).kz = sample_grid_idx(:, 2);
        us_msk_cand(n_success_iters).omegaz = 2*pi * ((0:1:nkz-1) - (nkz-1)/2) /nkz; % Similar to the MICA case in get_ky_omegaz_us_msk.m
        
        % Fill in big gaps along the ky dimension
        us_msk_cand(n_success_iters) = us_msk_fix_big_gap(us_msk_cand(n_success_iters), nky, ker_sup_nky, ky_gap_tol);
        
        % Increase distance tolerance
        dist_tol = dist_tol * 1.02;
        if debug
            fprintf('Generated 1 mask.\t Increased distance tolerance.\n');
        end
    else                    % Broke the loop because all possible positions have been examined but still couldn't generate a mask
        dist_tol = dist_tol * 0.98;
        if debug
            fprintf('Couldn''t generate mask.\t Decreased distance tolerance.\n');
        end
    end
end

% -- Choose a best undersample mask from all candidates using incoherency metric
[~, max_intr] = us_msk_psf(us_msk_cand, nky, nz, p.kydir, pdf);
[~, idx_sort_intr] = sort(max_intr(:), 1, 'ascend');

% Choose mask with minimum interference
us_msk_cand = us_msk_cand(idx_sort_intr);
us_msk = us_msk_cand(1);

% Adjust ordering so that ky is traversed monotonically
[us_msk.ky, idx_sort_ky] = sort(us_msk.ky(:), 1, 'ascend');
us_msk.kz = us_msk.kz(idx_sort_ky);

return

function grid_idx = coord2grid(point, grid_sz, max_coord)

grid_idx1 = 1 + round((grid_sz(1)-1) * point(:,1) / max_coord(1));
grid_idx2 = 1 + round((grid_sz(2)-1) * point(:,2) / max_coord(2));
grid_idx = [grid_idx1, grid_idx2];

return

function res = inNeighbourhood(sample_coord, new_point, dist_tol)

res = false;

if ~isempty(sample_coord)
    dist1 = abs(sample_coord(:,1) - new_point(1));
    dist2 = abs(sample_coord(:,2) - new_point(2));
    if any((dist1 < dist_tol(1)) & (dist2 < dist_tol(2))) % new_point inside at least one sample's adjacent rectangular region
        res = true;
    end
end

return

function res = acquiredBeforeAlongNonRepeatableDims(sample_coord, new_point, non_repeat_dims, grid_sz, max_coord)

res = false;

if ~isempty(non_repeat_dims) && ~isempty(sample_coord)
    sample_grid_idx = coord2grid(sample_coord, grid_sz, max_coord);
    newPos = coord2grid(new_point, grid_sz, max_coord);
    for dim_idx = 1 : length(non_repeat_dims)
        dim = non_repeat_dims(dim_idx);
        if any(sample_grid_idx(:, dim) == newPos(dim))
            res = true;
            break;
        end
    end
end

return