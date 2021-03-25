function dat = sense_map_espirit(dat, nmap, ncalib, ksize, eigThresh_1, eigThresh_2, partial_ky, ny_part, debug, tmp_dir)
%
% function dat = sense_map_espirit(dat, [nmap=1], [ncalib=[20,20]], [ksize=[6,6]], [eigThresh_1=round(1.6*prod(ksize))(if using Matlab code) or 0.02(if using compiled C code)], [eigThresh_2=0.95], [partial_ky=false], [ny_part=[]], [debug=false], [tmp_dir='tmp_data/'])
% 
% Calculates sensitivity maps using the approach in ESPIRiT. This is a wrapper
% function to the functions in the ESPIRiT package (Matlab code from Michael Lustig, C code from Martin Uecker).
% 
% Inputs
%   dat         - 2D k-space data with fully sampled calibration area. Zero-padded for partial ky acquisition. Dim: [UndersampledDim1, UndersampledDim2, Coil].
%   nmap        - Number of sensitivity map sets to return.
%   ncalib      - Size of the calibration area. Dim: [Dim1, Dim2].
%   ksize       - Kernel size, for forming the calibration matrix.
%   eigThresh_1 - Threshold to define the nullspace from the 1st SVD in the k-space. When using Matlab code, the largest eigThresh_1 singular values will remain; When using compiled C code, eigThresh_1 is in the range [0,1], and the cutoff threshold for the singular values will be maxSingularValue*eigThresh_1.
%   eigThresh_2 - Threshold of eigen vector decomposition in image space.
%   partial_ky  - True: partial ky acquisition, the compiled c code doesn't work in this case. False: full ky acquisition.
%   ny_part     - Number of acquired ky lines in a partial ky acquisition. Only needed if partial_ky is true.
%   debug       - True: project the coil images onto the maps to check whether the maps are good.
%   tmp_dir     - Temporary directory for saving intermediate data. Only used when Martin Uecker's compiled C code is used.
%
% Output
%   dat         - 2D sensitivity maps. Dim: [UndersampledDim1, UndersampledDim2, Coil, SetOfSensitivityMaps].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2013

%% Set defaults
if ~exist('nmap', 'var') || isempty(nmap)
    nmap = 1;
end

if ~exist('ncalib', 'var') || isempty(ncalib)
    ncalib = [20, 20];
end

if ~exist('ksize', 'var') || isempty(ksize)
    ksize = [6, 6];
end

if ~exist('eigThresh_2', 'var') || isempty(eigThresh_2)
    eigThresh_2 = 0.95;
end

if ~exist('partial_ky', 'var') || isempty(partial_ky)
    partial_ky = false;
end

if ~exist('ny_part', 'var')
    ny_part = [];
end

if isempty(ny_part) && partial_ky
    error('Must specify the number of acquired ky lines for a partial ky acquisition.');
end

if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

%% Size
[nx, ny, nc] = size(dat);

if ksize(1) > nx
    ksize(1) = nx;
end

if ksize(1) < 1
    ksize(1) = 1;
end

if ksize(2) > ny
    ksize(2) = ny;
end

if ksize(2) < 1
    ksize(2) = 1;
end

%% Calculate maps
if (ismac || isunix) && ~isempty(getenv('TOOLBOX_PATH')) && (ncalib(1) == ncalib(2)) && ~partial_ky % Use compiled C code
    
    if ~exist('eigThresh_1', 'var') || isempty(eigThresh_1)
        eigThresh_1 = 0.02;
    end
    
    if ~exist('tmp_dir', 'var') || isempty(tmp_dir)
        tmp_dir = 'tmp_data/';
    end
    if ~exist(tmp_dir, 'dir')
        mkdir(tmp_dir);
    end
    
    ksize = max(ksize);
    dat = reshape(dat, [1, nx, ny, nc]);
    dat_fname = [tmp_dir '/dat'];
    smap_fname = [tmp_dir '/smap'];
    
    % ESPIRiT calibration (using a maximum calibration region of size 'ncalib')
    % example command: ! ecalib -t 0.02 -c 0.95 -k 6 -r 64 -m 1 tmp_data/dat tmp_data/smap
    writecfl(dat_fname, dat);
    cmd = ['! ecalib -t ' num2str(eigThresh_1) ' -c ' num2str(eigThresh_2) ' -k ' num2str(ksize) ' -r ' num2str(ncalib) ' -m ' num2str(nmap) ' ' dat_fname ' ' smap_fname];
    eval(cmd);
    dat = readcfl(smap_fname);                                         % Dim: [1, UndersampledDim1(=nx), UndersampledDim2(=ny), nc, nmap]
    dat = reshape(dat, [nx, ny, nc, nmap]);
    
    rmdir(tmp_dir, 's');
    
else                                                                   % Use Matlab code
    
    if ~exist('eigThresh_1', 'var') || isempty(eigThresh_1)
        eigThresh_1 = round(1.6 * prod(ksize));
    end
    
    % K-space calibration area
    if debug
        dat_orig = dat;                                                % Keep the original data
    end
    
    if partial_ky && (ncalib(2) > (ny_part-ny/2)*2)
        dat = dat(:, ny_part-ncalib(2)+1 : ny_part, :);                % Keep ncalib(2) acquired ky lines
    end
    dat = crop(dat, [ncalib(1), ncalib(2), nc]);                       % dat: k-space calibration area

    % Compute calibration matrix, perform 1st SVD and convert singular vectors into k-space kernels
    [dat, S] = dat2Kernel(dat, ksize);                                 % dat: k-space kernels (Use the same variable name to reduce memory usage). S: singular values.
    
    % Crop kernels and compute eigen-value decomposition in image space to get maps
    [dat, W] = kernelEig(dat(:, :, :, 1:eigThresh_1), [nx, ny]);       % dat: sensitivity maps, Dim: [nx(UndersampledDim1), ny(UndersampledDim2), nc(Coils), nc(Set of maps)]. W: eigenvalues, Dim: [nx, ny, nc(Set of maps)]
    
    % Project coil images onto the maps. If the maps are good, sum_j(S_j^H * m_j)
    % should be equal to m, i.e. all the signal energy should live in the subspace
    % spanned by the eigenvectors with eigenvalue 1.
    if debug
        P = sum(repmat(ifft2c(dat_orig), [1,1,1,nc]) .* conj(dat), 3); % P: combined projections, Dim: [nx, ny, 1, nc(Set of maps)]
        raise2power = 1/3;
        figure; imshowALL((abs(P)).^raise2power);
        title(sprintf('Combined projection of the coil images onto the maps, raised to a power of %.2f, for all %d sets of maps.', raise2power, nc));
    end
    
    % Keep the desired number of map sets
    dat = dat(:, :, :, end-nmap+1:end);
    
    % Crop the sensitivity maps using soft-sense weights(i.e. Weight the eigenvectors with soft-senses eigen-values)
    weights = W(:, :, end-nmap+1:end);
    if eigThresh_2 == 0
        weights = (weights > eigThresh_2);
    else
        weights = (weights - eigThresh_2) ./ (1-eigThresh_2) .* (weights > eigThresh_2);
        weights = -cos(pi*weights)/2 + 1/2;
    end
    
    % Output sensitivity maps
    weights = repmat(reshape(weights, [nx,ny,1,nmap]), [1,1,nc,1]);
    if eigThresh_2 == 0
        dat = dat .* weights;
    else
        dat = dat .* sqrt(weights);
    end
end

return
