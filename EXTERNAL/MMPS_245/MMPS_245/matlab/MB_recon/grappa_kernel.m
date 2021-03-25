function ker = grappa_kernel(acs_dat, kersz, R, method, kspsz, kersz_min)
%
% function ker = grappa_kernel(acs_dat, kersz, R [, method, kspsz, kersz_min])
%
% Calculate reconstruction kernel for 1D GRAPPA.
%
% Inputs:
%   acs_dat   - AutoCalibrating Signal in k-space. Dim: [FE, PE, n_coil].
%   kersz     - KERnel SiZe, representing # of acquired points used in each
%               dimension to interpolate one missing point. Dim: [FE, PE].
%   R         - Reduction factor in the phase encoding direction.
%   method    - 'imSpace'(Default): the reconstruction will be carried out 
%                                   by multiplication in the image space.
%               'kSpace': the reconstruction will be carried out by 
%                         convolution in the k-space.
%   kspsz     - The full size of the k-space data, used only when 'method' 
%               is 'imSpace'. Dim: [n_fe(full), n_pe(full)].
%   kersz_min - If there aren't enough ACS points for the specified 'kersz', 
%               'kersz' will be reduced automatically. 'kersz_min' is the 
%               minimum kernel size allowed when 'kersz' is reduced. 
%               Default: [3, 2], dim: [FE, PE].
% Output:
%   ker       - The interpolation KERnel.
%               Dim: [FE(source), PE(source), n_coil(source), n_coil(target)].
%               If 'method' is 'kSpace', size of 'ker' is [kersz(1),
%               R*ceil(kersz(2)/2)*2-1, n_coil, n_coil].
%               If 'method' is 'imSpace', size of 'ker' is [kspsz(1),
%               kspsz(2), n_coil, n_coil].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

n_coil = size(acs_dat,3);

%% -- Check inputs.
if ~exist('method', 'var') || isempty(method)
    method = 'imSpace';
end

if ~exist('kersz_min', 'var') || isempty(kersz_min)
    kersz_min = [3, 2];
end

if ~all(kersz >= kersz_min)
    error('The specified kernel size is smaller than the minimum kernel size allowed.');
end

if n_coil < 2
    error('The number of coils is 1. Not a parallel reconstruction problem.');
end

if R <= 1
    error('Accelaration rate <= 1. Not a parallel reconstruction problem.');
end

if ~(strcmp(method, 'kSpace') || strcmp(method, 'imSpace'))
    error('''method'' must be either ''kSpace'' or ''imSpace''.');
end

if strcmp(method, 'imSpace') && (~exist('kspsz', 'var') || isempty(kspsz))
    error('Must specify the full k-space size if the reconstruction will be conducted in the image space.');
end

acs_sz = [size(acs_dat,1), size(acs_dat,2)]; % Dim: [FE, PE]
[n_instances, acs_blk_sz] = calc_instances(kersz, acs_sz, R);

if ( sum(acs_sz < acs_blk_sz) ) > 0
    fprintf('Not enough ACS lines in %d of the FE and PE dimensions.\n',...
        sum(acs_sz < acs_blk_sz) );
    error('Not enough ACS lines.');
end

%% -- Reduce the kernel size if there aren't enough ACS points for the specified kernel size.
if n_instances < prod(kersz)*n_coil
    fprintf(['   The specified kernel size is [%d, %d] in [FE, PE]. ' ...
        'Not enough ACS points for the specified kernel size.\n'], kersz(1), kersz(2) );
    disp('   Reducing the kernel size...');
    
    for dim = 1 : 2                 % First reduce kersz in FE, then in PE
        while (n_instances < prod(kersz)*n_coil) && (kersz(dim) > kersz_min(dim))
            kersz(dim) = kersz(dim) - 1;
            [n_instances, acs_blk_sz] = calc_instances(kersz, acs_sz, R);
        end
    end
    
    if n_instances < prod(kersz)*n_coil
        fprintf('   The minimum kernel size allowed is [%d, %d] in [FE, PE].\n', kersz_min(1), kersz_min(2));
        fprintf('   Total number of instances available in the ACS area to fit the minimum kernel is %d.\n', n_instances);
        fprintf('   Total number of coefficients to fit for the minimum kernel is prod(kersz)*n_coil = %d.\n', prod(kersz)*n_coil);
        error('Unable to estimate the minimum kersz size allowed with the current ACS area size. Please increase the ACS area.');
    else
        fprintf('   The reduced kernel size is [%d, %d] in [FE, PE].\n', kersz(1), kersz(2));
    end
end

%% -- Set up a matrix with all of the calibration data.
acs_mtx = zeros(prod(acs_blk_sz(1:2)), n_instances, n_coil);
for idx = 1 : n_coil
    acs_mtx(:, :, idx) = im2col(acs_dat(:,:,idx), acs_blk_sz, 'sliding');        % Dim: [acs_blk_sz(1)*acs_blk_sz(2), n_instances, n_coil]
end
acs_mtx = reshape(acs_mtx, [acs_blk_sz(1), acs_blk_sz(2), n_instances, n_coil]); % Dim: [acs_blk_sz(1)(i.e. kersz(1)), acs_blk_sz(2), n_instances, n_coil]

%% -- Source matrix.
src_mtx = acs_mtx(:, 1:R:end, :, :);                                             % Dim: [kersz(1), kersz(2), n_instances, n_coil]
src_mtx = reshape(permute(src_mtx,[3,1,2,4]), [n_instances, kersz(1)*kersz(2)*n_coil]); % Dim: [n_instances, kersz(1)*kersz(2)*n_coil]

%% -- Target matrix.
tgt_fe_pos = floor(kersz(1)/2)+1;
tgt_pe_pos = R*floor((kersz(2)-1)/2)+2 : R*floor((kersz(2)-1)/2)+R;
tgt_mtx = acs_mtx(tgt_fe_pos, tgt_pe_pos, :, :);                            % Dim: [1, R-1, n_instances, n_coil]
tgt_mtx = reshape(permute(tgt_mtx,[3,2,4,1]), [n_instances, (R-1)*n_coil]); % Dim: [n_instances, (R-1)*n_coil]

%% -- Solve the kernel.
% -- Raw kernel
AtA = src_mtx' * src_mtx; 
lambda = norm(AtA, 'fro') / size(src_mtx,2) * 0.02;            % Regularization parameter, from Miki Lustig's code                      
%ker = inv(AtA + eye(size(AtA))*lambda) * src_mtx' * tgt_mtx;   % Dim: [kersz(1)*kersz(2)*n_coil, (R-1)*n_coil]
ker = (AtA + eye(size(AtA))*lambda)\(src_mtx' * tgt_mtx);   % Dim: [kersz(1)*kersz(2)*n_coil, (R-1)*n_coil]
ker = reshape(ker, [kersz(1), kersz(2), n_coil, R-1, n_coil]); % Dim: [kersz(1)(source), kersz(2)(source), n_coil(source), (R-1)(target), n_coil(target)]

% -- Autocorrelation kernel
autocorr_ker = zeros(kersz(1), R*ceil(kersz(2)/2)*2-1, n_coil, n_coil); % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]
for tgt_coil = 1 : n_coil
    for tgt_pe = 1 : R-1
        autocorr_ker( :, R-tgt_pe:R:R*kersz(2)-1, :, tgt_coil ) = ker(:, :, :, tgt_pe, tgt_coil);
    end
    autocorr_ker( floor(kersz(1)/2)+1, R*floor((kersz(2)-1)/2)+R, ...
        tgt_coil, tgt_coil) = 1;
end

% -- Convolution kernel
conv_ker = flip(flip(autocorr_ker, 1), 2); % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]

if strcmp(method, 'kSpace')
    ker = conv_ker;                              % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]
end

% -- Image space multiplication kernel
if strcmp(method, 'imSpace')
    ker = zpad(sqrt(kspsz(1)*kspsz(2))*conv_ker, [kspsz(1), kspsz(2), n_coil, n_coil]);
    ker = ifft2c(ker);                           % Dim: [FE(src), PE(src), n_coil(src), n_coil(tgt)]
end


return;

function [n_instances, acs_blk_sz] = calc_instances(kersz, acs_sz, R)
%
% Calculates the total number of instances for kernel fitting in GRAPPA.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

acs_blk_sz = [kersz(1), R*(kersz(2)-1) + 1];   % Size of one block for kernel fitting, dim: [FE, PE]
n_instances_indv = acs_sz - acs_blk_sz + 1;    % # of instances for kernel fitting in [FE,PE]
n_instances = prod(n_instances_indv);          % Total # of instances for kernel fitting

return;