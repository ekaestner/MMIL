function ker = slice_grappa_kernel(dat_mux, dat_sing, us_msk, kersz, method, kspsz, use_split, odd_even_fit, kersz_min)
%
% function ker = slice_grappa_kernel(dat_mux, dat_sing, us_msk, kersz, [method='imSpace'], [kspsz], [use_split=true], [odd_even_fit=false], [kersz_min=[3,2]])
%
% Calculate reconstruction kernel for (Split-)Slice-GRAPPA.
% References: Standard Slice-GRAPPA: Kawin Setsompop, et al. MRM 2012;67(5):1210-24.
%             Split-Slice-GRAPPA:    Stephen F. Cauley, et al. MRM 2014;72(1):93-102.
%             Standard Slice-GRAPPA with odd even kernel fitting: Kawin Setsompop, et al. NeuroImage 2012;63:569-580.
%
% Inputs
%   dat_mux   - Slice-multiplexed k-space data. Only the acquired data, no zero
%               padding for in-plane acceleration. Dim: [Kx(=nx), Ky(=ny), Coil(=nc)].
%               Needed for standard slice-GRAPPA, not needed for split-slice-GRAPPA.
%   dat_sing  - Single-slice k-space data. Dim: [Kx(=nx), Ky(=ny), Coil(=nc),
%               SimultaneousSlice Z(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)].
%   us_msk    - Undersample mask with fields 'ky', 'kz', 'omegaz'.
%               See get_ky_omegaz_us_msk.m for details.
%   kersz     - KERnel SiZe. Dim: [Kx, Ky].
%   method    - How the reconstruction will be carried out.
%               'imSpace': Multiplication in image space.
%               'kSpace':  Convolution in k-space.
%   kspsz     - Full image matrix size, used only when 'method' is 'imSpace'.
%               Dim: [FE(=nx_full), PE(=ny_full)].
%   use_split - True: use split-slice-GRAPPA; False: use standard slice-GRAPPA.
%   odd_even_fit - True: fit different kernels to the odd and even ky lines.
%   kersz_min - If there aren't enough ACS points for the specified 'kersz',
%               'kersz' will be reduced automatically. 'kersz_min' is the
%               minimum kernel size allowed when 'kersz' is reduced. Dim: [Kx, Ky].
%
% Output
%   ker       - Structures for the interpolation kernels. If odd_even_fit is true, length(ker)=2; Otherwise length(ker)=1.
%               Fields include:
%               ker - The interpolation kernels. Dim: [FE(source), PE(source), Coil(=nc, source), Coil(=nc, target), SimultaneousSlice Z(=nz, target)].
%                     If 'method' is 'kSpace', size of 'ker' is [kersz(1), kersz(2), nc, nc, nz(z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%                     If 'method' is 'imSpace', size of 'ker' is [kspsz(1), kspsz(2), nc, nc, nz(z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%               ky_start - The starting target ky line for this kernel. Useful only when odd_even_fit is true (ky_start will be used to sort out the odd and even ky lines).
%
% (c) Kangrong Zhu,     Stanford University     Aug 2014

SPECIFIED_OMEGAZ_ENC = 0;

%% -- Check inputs
[nx, ny, nc, nz] = size(dat_sing);

if ~exist('use_split', 'var') || isempty(use_split)
    use_split = true;
end

if ~use_split
    if isempty(dat_mux)
        error('Must input slice-multiplexed data for standard slice-GRAPPA recon.');
    else
        if (size(dat_mux, 1) ~= nx) || (size(dat_mux, 2) ~= ny) || (size(dat_mux, 3) ~= nc)
            error('Sizes of slice-multiplexed data and single-slice data mismatch.');
        end
    end
end

if ~exist('method', 'var') || isempty(method)
    method = 'imSpace';
end

if ~exist('odd_even_fit', 'var') || isempty(odd_even_fit)
    odd_even_fit = false;
end

if ~exist('kersz_min', 'var') || isempty(kersz_min)
    kersz_min = [3, 2];
end

if ~all(kersz >= kersz_min)
    error('The specified kernel size is smaller than the minimum kernel size allowed.');
end

if nc < 2
    error('The number of coils is 1. Not a parallel reconstruction problem.');
end

if ~(strcmp(method, 'kSpace') || strcmp(method, 'imSpace'))
    error('''method'' must be either ''kSpace'' or ''imSpace''.');
end

if strcmp(method, 'imSpace') && (~exist('kspsz', 'var') || isempty(kspsz))
    error('Must specify the full image matrix size since the reconstruction will be conducted in the image space.');
end

acssz = [nx, ny];                                                                   % Size of ACS area. Dim: [Kx, Ky]
n_instances = calc_instances(kersz, acssz);

if any(acssz < kersz)
    fprintf('Not enough ACS lines in %d of the FE and PE dimensions.\n', sum(acssz < kersz) );
    error('Not enough ACS lines.');
end

%% -- Reduce the kernel size if there aren't enough ACS points for the specified kernel size
if n_instances < prod(kersz)*nc
    fprintf(['   The specified kernel size is [%d, %d] in [FE, PE]. ' ...
        'Not enough ACS points for the specified kernel size.\n'], kersz(1), kersz(2) );
    disp('   Reducing the kernel size...');
    
    for dim = 1 : 2                                                                 % First reduce kersz in FE, then in PE
        while (n_instances < prod(kersz)*nc) && (kersz(dim) > kersz_min(dim))
            kersz(dim) = kersz(dim) - 1;
            n_instances = calc_instances(kersz, acssz);
        end
    end
    
    if n_instances < prod(kersz)*nc
        fprintf('   The minimum kernel size allowed is [%d, %d] in [FE, PE].\n', kersz_min(1), kersz_min(2));
        fprintf('   Total number of instances available in the ACS area to fit the minimum kernel is %d.\n', n_instances);
        fprintf('   Total number of coefficients to fit for the minimum kernel is prod(kersz)*nc = %d.\n', prod(kersz)*nc);
        error('Unable to estimate the minimum kersz size allowed with the current ACS area size. Please increase the ACS area.');
    else
        fprintf('   The reduced kernel size is [%d, %d] in [FE, PE].\n', kersz(1), kersz(2));
    end
end

%% -- Add FTz encoding phase to each slice
us_msk = encode_ftz_pha(us_msk, nz, SPECIFIED_OMEGAZ_ENC);                          % us_msk.ftz_pha is the encoding phase added to each individual slice for each sample. Dim: [nsamp(=ny), nz(z indices -floor(nz/2):1:(ceil(nz/2)-1))])
for z = 1 : nz
    pha_mtx = repmat(us_msk.ftz_pha(:, z).', [nx, 1]);
    for coil = 1 : nc
        dat_sing(:, :, coil, z) = dat_sing(:, :, coil, z) .* pha_mtx;
    end
end

%% -- Solve kernel
if odd_even_fit
    nker = 2;
    n_instances_per_ky = acssz(1) - kersz(1) + 1;
    first_ky_tgt = floor(kersz(2)/2)+1; % First possible target point along ky
    last_ky_tgt = ny - first_ky_tgt + 1; % Last possible target point along ky
    ker{1}.ky_start = first_ky_tgt;
    ker{2}.ky_start = first_ky_tgt + 1;
    instances{nker} = [];
    for ker_idx = 1 : nker
        instances{ker_idx} = [];
        ky_indices = ker{ker_idx}.ky_start : 2 : last_ky_tgt;
        for ky_idx = 1 : length(ky_indices);
            instances{ker_idx} = [instances{ker_idx}, (ky_indices(ky_idx)-first_ky_tgt)*n_instances_per_ky + (1 : n_instances_per_ky)];
        end
    end
else
    nker = 1;
    instances{1} = 1 : n_instances;
end

for ker_idx = 1 : nker
    if ~use_split % Standard slice-GRAPPA
        % -- Instances to keep in the source and target matrices
        instances_to_keep = instances{ker_idx};
        
        % -- Source matrix
        src_mtx = setup_blk_mtx(dat_mux, kersz, n_instances); % Dim: [kersz(1), kersz(2), n_instances, nc]
        src_mtx = reshape(permute(src_mtx,[3,1,2,4]), [n_instances, kersz(1)*kersz(2)*nc]); % Dim: [n_instances, kersz(1)*kersz(2)*nc]
        src_mtx = src_mtx(instances_to_keep, :); % Dim: [length(instances_to_keep), kersz(1)*kersz(2)*nc]
        
        % -- Target matrix
        tgt_dat_mtx = setup_blk_mtx(dat_sing, kersz, n_instances); % Dim: [kersz(1), kersz(2), n_instances, nc, nz]
        
        tgt_fe_pos = floor(kersz(1)/2)+1;
        tgt_pe_pos = floor(kersz(2)/2)+1;
        tgt_mtx = tgt_dat_mtx(tgt_fe_pos, tgt_pe_pos, :, :, :); % Dim: [1, 1, n_instances, nc, nz]
        tgt_mtx = reshape(permute(tgt_mtx,[3,4,5,1,2]), [n_instances, nc*nz]); % Dim: [n_instances, nc*nz]
        tgt_mtx = tgt_mtx(instances_to_keep, :); % Dim: [length(instances_to_keep), nc*nz]
        
        % -- Raw kernel
        AtA = src_mtx' * src_mtx;
        lambda = norm(AtA, 'fro') / size(src_mtx,2) * 0.02; % Regularization parameter, from Miki Lustig's code
        %ker{ker_idx}.ker = inv(AtA + eye(size(AtA))*lambda) * src_mtx' * tgt_mtx; % Dim: [kersz(1)*kersz(2)*nc, nc*nz]
        ker{ker_idx}.ker = inv(AtA + eye(size(AtA))*lambda)\(src_mtx' * tgt_mtx); % Dim: [kersz(1)*kersz(2)*nc, nc*nz]
        
    else % Split-Slice-GRAPPA
        % -- Instances to keep in the source and target matrices
        instances_to_keep = zeros(1, length(instances{ker_idx})*nz);
        for z = 1 : nz
            instances_to_keep((z-1)*length(instances{ker_idx})+1 : z*length(instances{ker_idx})) = instances{ker_idx} + (z-1)*n_instances;
        end
        
        % -- Set up a matrix with all of the calibration data
        ker_blk_mtx = setup_blk_mtx(dat_sing, kersz, n_instances); % Dim: [kersz(1), kersz(2), n_instances, nc, nz]
        
        % -- Source matrix
        src_mtx = zeros(n_instances*nz, kersz(1)*kersz(2)*nc); % Dim: [n_instances*nz, kersz(1)*kersz(2)*nc]
        for z = 1 : nz
            src_mtx((z-1)*n_instances+1:z*n_instances, :) = reshape(permute(ker_blk_mtx(:, :, :, :, z), [3,1,2,4]), [n_instances, kersz(1)*kersz(2)*nc]);
        end
        src_mtx = src_mtx(instances_to_keep, :); % Dim: [length(instances_to_keep), kersz(1)*kersz(2)*nc]
        
        % -- Target matrix
        tgt_fe_pos = floor(kersz(1)/2)+1;
        tgt_pe_pos = floor(kersz(2)/2)+1;
        tgt_mtx_raw = zeros(n_instances*nz, nc, nz); % Dim: [n_instances*nz, nc, nz]
        for z = 1 : nz
            tgt_mtx_raw((z-1)*n_instances+1:z*n_instances, :, z) = permute(ker_blk_mtx(tgt_fe_pos, tgt_pe_pos, :, :, z), [3,4,1,2,5]);
        end
        tgt_mtx_raw = tgt_mtx_raw(instances_to_keep, :, :); % Dim: [length(instances_to_keep), nc, nz]
        
        % -- Raw kernel
        ker{ker_idx}.ker = zeros(kersz(1)*kersz(2)*nc, nc, nz);
        for z = 1 : nz
            tgt_mtx = tgt_mtx_raw(:, :, z); % Dim: [length(instances_to_keep), nc]
            AtA = src_mtx' * src_mtx;
            lambda = norm(AtA, 'fro') / size(src_mtx,2) * 0.02; % Regularization parameter, from Miki Lustig's code
            %ker{ker_idx}.ker(:, :, z) = inv(AtA + eye(size(AtA))*lambda) * src_mtx' * tgt_mtx; % Dim: [kersz(1)*kersz(2)*nc, nc, nz]
            ker{ker_idx}.ker(:, :, z) = inv(AtA + eye(size(AtA))*lambda)\(src_mtx' * tgt_mtx); % Dim: [kersz(1)*kersz(2)*nc, nc, nz]
        end
    end
    
    % -- Autocorrelation kernel
    autocorr_ker = reshape(ker{ker_idx}.ker, [kersz(1), kersz(2), nc, nc, nz]); % Dim: [kersz(1)(source), kersz(2)(source), nc(source), nc(target), nz(target)]
    
    % -- Convolution kernel
    conv_ker = FLIP(FLIP(autocorr_ker, 1), 2); % Dim: [Kx(src), Ky(src), nc(src), nc(tgt), nz(tgt)]
    
    if strcmp(method, 'kSpace')
        ker{ker_idx}.ker = conv_ker; % Dim: [Kx(=kersz(1), src), Ky(=kersz(2), src), nc(src), nc(tgt), nz(tgt)]
    end
    
    % -- Image space multiplication kernel
    if strcmp(method, 'imSpace')
        ker{ker_idx}.ker = zpad(sqrt(kspsz(1)*kspsz(2))*conv_ker, [kspsz(1), kspsz(2), nc, nc, nz]);
        ker{ker_idx}.ker = ifft2c(ker{ker_idx}.ker); % Dim: [X(=kspsz(1), src), Y(=kspsz(2), src), nc(src), nc(tgt), nz(tgt)]
    end
end

return

function n_instances = calc_instances(kersz, acssz)
%
% Calculates the total number of instances for kernel fitting in slice-GRAPPA.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2014

n_instances_indv = acssz - kersz + 1;  % # of instances for kernel fitting in [FE, PE]. The size of one block for kernel fitting is the same as the kernel size.
n_instances = prod(n_instances_indv);  % Total # of instances for kernel fitting

return;

function mtx = setup_blk_mtx(dat, blk_sz, n_instances)
%
% function mtx = setup_blk_mtx(dat, blk_sz, [n_instances=calc_instances(blk_sz, [nx, ny]);])
%
% Sets up a matrix for all instances for a block with size blk_sz.
%
% Inputs
%   dat         - K-space data. Dim: [Kx(=nx), Ky(=ny), Coil(=nc), SimultaneousSlice Z(=nz)].
%   blk_sz      - Block size. Dim: [Kx, Ky].
%   n_instances - Number of instances for the input block size.
%
% Output
%   mtx         - Data matrix for all the instances for this block size. Dim: [blk_sz(1), blk_sz(2), n_instances, nc, nz].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2014

[nx, ny, nc, nz] = size(dat);
if ~exist('n_instances', 'var') || isempty(n_instances)
    n_instances = calc_instances(blk_sz, [nx, ny]);
end

mtx = zeros(prod(blk_sz(1:2)), n_instances, nc, nz);
for z = 1 : nz
    for coil = 1 : nc
        mtx(:, :, coil, z) = im2col(dat(:,:,coil,z), blk_sz, 'sliding');  % Dim: [blk_sz(1)*blk_sz(2), n_instances, nc, nz]
    end
end
mtx = reshape(mtx, [blk_sz(1), blk_sz(2), n_instances, nc, nz]);          % Dim: [kersz(1), kersz(2), n_instances, nc, nz]

return