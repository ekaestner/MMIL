function im_noise = mux_calc_im_noise(us_msk, smap, psi_mtx, whiten_type, pha_coe, ramp_flt, ccmtx, im_noise_type, propagate_type, coil_noise_std, debug)
%
% function im_noise = mux_calc_im_noise(us_msk, smap, [psi_mtx=eye(nci)], [whiten_type='psi_mtx'], [pha_coe=zeros(2, ny, nsl, nci)], [ramp_flt], [ccmtx], [im_noise_type='actual'], [propagate_type='diagblk'], [coil_noise_std=sqrt(diag(psi_mtx))], [debug=false])
%
% Basing on the Hybrid-Space-SENSE reconstruction model, this function calculates the noise maps of the reconstructed mux images (Can also handle single-slice acquisitions).
% Starting with the coil noise covariance matrix, this function calculates the image noise maps by following the noise propagation in the reconstruction steps, all of which can be modeled as linear transformations.
% It is not necessary to consider the different eddy current effects in the simultaneous slices here, because these effects are not part of the encoding we want to apply but rather artifacts introduced in EPI-type acquisition.
%
% Inputs
%   us_msk         - Undersample mask on the ky-omegaz plane, must not be empty. A structure with
%                    fields 'ky', 'kz', 'omegaz'. length(us_msk) = nmsk.
%                    See get_ky_omegaz_us_msk.m for details.
%   smap           - Sensitivity maps, must not be empty. Dim: [X(=nx), Y(=ny), Echo(=nec=1), Slice(=nsl),
%                    SimultaneousSliceZ(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))),
%                    Coil(=nc, after coil compression), SetOfSensitivityMaps(=nmap)];
%                    For single-slice acquisition nz would be 1.
%   psi_mtx        - Receiver coil noise covariance matrix. Dim: [Coil(=nci, before coil compression),
%                    Coil(=nci, before coil compression)].
%   whiten_type    - Type of data used in the data whitening step. See get_whiten_mtx.m,
%                    epi_process_rawdata.m, mux_epi_params.m for details.
%                    'none':           No data whitening is applied.
%                    'coil_noise_std': Use coil noise standard deviation as an approximation.
%                    'psi_mtx':        Use the full noise covariance matrix.
%   pha_coe        - Coefficients for x-ky phase correction, loaded from ref.dat file or calculated from
%                    reference scan p-file or reference scan k-space data. See rawload_ref.m, epi_pha_correct.m for details.
%                    Dim: [num_coe_per_ky_line(=2, 0th order, 1st order), Ky(=ny), Slice(=nsl), Coil(=nci)].
%   ramp_flt       - Ramp sampling correction filter. Dim: [KxAfterRampSamplingCorrection(=nx),
%                    KxBeforeRampSamplingCorrection(=nxi)]. See rawload_vrgf.m, epi_vrgf_correct.m
%                    for details. When not exist or empty, assume no ramp sampling was used.
%   ccmtx          - Coil compression matrix. Dim: [CoilBeforeCompression(=nci), CoilAfterCompression(=nc),
%                    X(=1 if one set of compression is used for all x positions; =nx if GCC is used
%                    and different compression matrices are used for different x positions),
%                    Echo(=nec=1), Slice(=nsl)]. See coil_compress.m for details. When not exist or
%                    empty, assume no coil compression was used.
%   im_noise_type  - 'actual': Calculate noise propagation in the actual reconstruction which uses
%                              Tikhonov regularization.
%                    'snr_optimal': Calculate noise propagation in an SNR optimal reconstruction
%                              which assumes the sample noise covariance is known, as in the SENSE paper.
%   propagate_type - How to calculate the noise propagation.
%                    'cov': Keep track of the full sample noise covariance matrix.
%                           propagate_type must be 'cov' if im_noise_type is 'snr_optimal'.
%                    'diagblk': After the 1D iDFTx step, keep track of only the diagonal blocks of
%                               the sample noise covariance matrix. Each diagonal block corresponds
%                               to one x position. This is much faster than 'cov' when im_noise_type
%                               is 'actual', but this doesn't work when im_noise_type is 'snr_optimal'.
%   coil_noise_std - Coil noise standard deviation, i.e. sqrt of the diagonal elements in the coil
%                    noise covariance matrix. Needed only when whiten_type is 'coil_noise_std'.
%                    Although this can be calculated from the input 'psi_mtx', it is still listed as
%                    a separate input so that values used in an actual reconstruction can be passed
%                    to this function.
%   debug          - True: print debug messages.
%
% Output
%   im_noise       - Image noise (standard deviation) matrix. Dim: [X(=nx), Y(=ny), Echo(=nec=1),
%                    Slice*SimultaneousSliceZ(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)),
%                    MaskType(=nmsk)].
%
% (c) Kangrong Zhu,     Stanford University    March 2015

%% Parse inputs and set parameters

% us_msk
if ~exist('us_msk', 'var') || isempty(us_msk)
    error('No undersampling mask(s).');
end

if length(us_msk(1).ky) ~= length(us_msk(1).kz)
    error('Length of ''us_msk.ky'' does not match length of ''us_msk.kz''.');
end

nsamp = length(us_msk(1).ky); % Number of acquired points on the ky-omegaz plane
nmsk  = length(us_msk); % Number of undersample masks for the ky-omegaz plane
if nmsk == 0
    error('Undersampling mask empty.');
end

% smap
if ~exist('smap', 'var') || isempty(smap)
    error('No input sensitivity maps.');
end

[nx, ny, nec, nsl, nz, nc, nmap] = size(smap); % nc: Number of coils after coil compression

if nec ~= 1
    error('Number of echoes must be 1 in sensitivity map matrix.');
end

% ccmtx
coil_compress = (exist('ccmtx', 'var') && ~isempty(ccmtx)); % True: use coil compression

if coil_compress
    nci = size(ccmtx, 1); % nci: Initial number of coils before coil compression
    if nc ~= size(ccmtx, 2)
        error('Number of coils after coil compression mismatches in inputs ''smap'' and ''ccmtx''.');
    end
else
    nci = nc;
end

% psi_mtx
if ~exist('psi_mtx', 'var') || isempty(psi_mtx)
    psi_mtx = eye(nci);
end

% whiten_type
if ~exist('whiten_type', 'var') || isempty(whiten_type)
    whiten_type = 'psi_mtx';
end

if ~ismember(whiten_type, {'coil_noise_std', 'psi_mtx', 'none'})
    error('whiten_type must be either ''coil_noise_std'' or ''psi_mtx'' or ''none''.');
end

% pha_coe
if ~exist('pha_coe', 'var') || isempty(pha_coe)
    pha_coe = zeros(2, ny, nsl, nci);
end
pha_coe = pha_coe(:, us_msk(1).ky, :, :); % Dim: [num_coe_per_ky_line(=2, 0th order, 1st order), AcquiredSampleOnKyOmegazPlane(=nsamp), Slice(=nsl), Coil(=nci)]

% ramp_flt
vrgf = (exist('ramp_flt', 'var') && ~isempty(ramp_flt)); % True: ramp sampling was used

% im_noise_type
if ~exist('im_noise_type', 'var') || isempty(im_noise_type)
    im_noise_type = 'actual';
end

% propagate_type
if ~exist('propagate_type', 'var') || isempty(propagate_type)
    propagate_type = 'diagblk';
end

if strcmp(im_noise_type, 'snr_optimal') && ~strcmp(propagate_type, 'cov')
    error('propagate_type must be ''cov'' when im_noise_type is ''snr_optimal''.');
end

if ~ismember(propagate_type, {'cov', 'diagblk'})
    error('propagate_type must be either ''cov'' or ''diagblk''.');
end

% coil_noise_std
if strcmp(whiten_type, 'coil_noise_std') && (~exist('coil_noise_std', 'var') || isempty(coil_noise_std))
    coil_noise_std = sqrt(diag(psi_mtx));
end

% debug
if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

% Constants
p.KEEP_ORIG_SZ = 1;
p.ONLY_ACTUAL_COILS = 0;
p.ACTUAL_AND_VIRTUAL_COILS = 1;
p.ONLY_VIRTUAL_COILS = 2;

% Virtual coil concept
% Currently only works if p.add_vcc equals p.ONLY_ACTUAL_COILS or equals p.ONLY_VIRTUAL_COILs, but still keep this parameter so that the virtual coil concept might be added in the future (the noise in the actual and virtual coils are correlated).
p.add_vcc = p.ONLY_ACTUAL_COILS; % Currently just set to p.ONLY_ACTUAL_COILS
if p.add_vcc == p.ACTUAL_AND_VIRTUAL_COILS
    error('Function currently doesn''t work with reconstruction that uses both actual and virtual coils.');
end

% Others
p.debug = debug;
p.sense_lambda = []; % See mux_epi_params.m for details. The lambda for Tikhonov regularization in the SENSE recon. Pass an empty matrix to use the default of norm(encode_mtx' * encode_mtx, 'fro') / size(encode_mtx,2) * p.sense_lambda_default_ratio.
p.sense_lambda_default_ratio = 0.02; % See mux_epi_params.m for details. Default regularization coefficient for the SENSE recon is norm(encode_mtx' * encode_mtx, 'fro') / size(encode_mtx,2) * p.sense_lambda_default_ratio.

%% Calculate step by step the noise propagation from raw k-space data to image space
% -- Receiver coil noise covariance matrix (i.e. the input psi_mtx with dimension [nci, nci]) is the starting point of the noise propagation

% -- Data whitening
if ~strcmp(whiten_type, 'none')
    switch whiten_type
        case 'coil_noise_std'
            whiten_mtx = get_data_whiten_mtx(coil_noise_std, whiten_type);
        case 'psi_mtx'
            whiten_mtx = get_data_whiten_mtx(full(psi_mtx), whiten_type);
    end
    psi_mtx = cov_linear_transform(whiten_mtx, psi_mtx); % psi_mtx Dim: [Coil(=nci), Coil(=nci)]. Output from function cov_linear_transform is always a sparse matrix.
    clear whiten_mtx;
end

% -- Number of acquired kx points. If vrgf == true, this is before ramp sampling correction
if vrgf
    nxi = size(ramp_flt, 2);
else
    nxi = nx;
end

% -- Case 1: whiten_type is 'psi_mtx'. After data whitening step, psi_mtx becomes eye(nci).
if strcmp(whiten_type, 'psi_mtx')
    % - Default Nyquist ghosting correction: has no effects on sample noise covariance matrix
    % DO NOTHING
    
    % - Set up sample noise covariance matrix for all acquired kx samples
    psi_mtx = speye(nxi); % Dim: [Kx(=nxi), Kx(=nxi)]. The nxi kx samples are independent, and the nci coils are I.I.D.
    
    % - Ramp sampling correction
    if vrgf
        psi_mtx = cov_linear_transform(ramp_flt, psi_mtx); % Dim: [Kx(=nx), Kx(=nx)]
    end
    
    % - 1D iDFT along kx, corresponding to ifftc.m
    idftmtx = get_idftmtx(nx);
    psi_mtx = cov_linear_transform(idftmtx, psi_mtx); % Dim: [X(=nx), X(=nx)]
    clear idftmtx;
    
    % Set up sample noise covariance matrix for all coils and all x positions
    switch propagate_type
        case 'cov'
            psi_mtx = kron(psi_mtx, speye(nci)); % Dim: [Coil->X(nci*nx), Coil->X(nci*nx)]
        case 'diagblk' % Use a cell array to store the nx diagonal nci-by-nci blocks in psi_mtx. Will track the diagonal blocks instead of the full matrix from this point on.
            psi_mtx_diagblk = cell(nx); % psi_mtx_diagblk{x} (x = 1,2,...nx) Dim: [Coil(=nci), Coil(=nci)]. psi_mtx_diagblk{x} is the x-th diagonal block in the (nci*nx)-by-(nci*nx) sample noise covariance matrix
            for x = 1 : nx
                psi_mtx_diagblk{x} = psi_mtx(x, x) * speye(nci);
            end
            clear psi_mtx;
    end
end

% -- Case 2: whiten_type is 'coil_noise_std' or 'none'. After data whitening step, psi_mtx is not an identity matrix.
if ismember(whiten_type, {'coil_noise_std', 'none'})
    % - Default Nyquist ghosting correction, step 1: 1D iDFT along kx
    idftmtx = get_idftmtx(nxi);
    idftmtx = kron(speye(nci), idftmtx); % Dim: [X->Coil(=nxi*nci), Kx->Coil(=nxi*nci)]
    psi_mtx = kron(psi_mtx, speye(nxi)); % Dim: [Kx->Coil(=nxi*nci), Kx->Coil(=nxi*nci)]. Use the Kronecker tensor product to set up the sample noise covariance matrix since the nxi kx positions are independent. speye(nxi) is a sparse identity matrix of size [nxi, nxi]. Kronecker tensor product of a full matrix and a sparse matrix is a sparse matrix.
    psi_mtx = cov_linear_transform(idftmtx, psi_mtx); % Dim: [X->Coil(=nxi*nci), X->Coil(=nxi*nci)]
    clear idftmtx;
    
    % - Set up a few matrices which will be used when looping over samples
    
    % 1D DFT matrix for Nyquist ghosting correction
    fdftmtx = get_fdftmtx(nxi);
    fdftmtx = kron(speye(nci), fdftmtx); % Dim: [Kx->Coil(=nxi*nci), X->Coil(=nxi*nci)]
    
    % Ramp sampling correction matrix
    ramp_flt_mtx = kron(speye(nci), ramp_flt); % Dim: [KxAfterRampSampingCorrection->Coil(=nx*nci), KxBeforeRampSampingCorrection->Coil(=nxi*nci)]
    
    % 1D iDFT matrix
    idftmtx = get_idftmtx(nx);
    idftmtx = kron(speye(nci), idftmtx); % Dim: [X->Coil(=nx*nci), Kx->Coil(=nx*nci)]
    
    % Indices to change an ordering of nx*nci to nci*nx
    idx_permute_nx_nci = zeros(1, nx*nci);
    for coil = 1 : nci
        idx_permute_nx_nci(coil : nci : nci*nx) = (coil-1)*nx+1 : 1 : coil*nx;
    end
end

% -- Set up a matrix for combining images from multiple sensitivity maps at one x position. This matrix is the same for all x positions, so this is also a diagonal block in the combine_mtx which includes all x positions.
if nmap > 1
    combine_mtx_diagblk = repmat(speye(ny*nz), [p.KEEP_ORIG_SZ, nmap]); % Dim: [ny*nz, ny*nz*nmap].
end

% -- Calculate echo by echo, slice by slice
im_noise = zeros(nx, ny, nec, nsl, nz, nmsk);
for echo = 1 : nec % nec is always 1
    for slice = 1 : nsl
        if p.debug
            fprintf('    slice = %d/%d\n', slice, nsl);
        end
        
        if coil_compress
            % -- The coil compression matrix is block diagonal. Write its diagonal blocks into a cell array
            coil_compress_mtx_diagblk = cell(nx);
            for x = 1 : nx
                if size(ccmtx, 3) == nx % Different coil compression matrices for different x positions
                    coil_compress_mtx_diagblk{x} = ccmtx(:, :, x, echo, slice); % Dim: [nci, nc]
                else % Same coil compression matrix for all x positions
                    coil_compress_mtx_diagblk{x} = ccmtx(:, :, 1, echo, slice); % Dim: [nci, nc]
                end
                coil_compress_mtx_diagblk{x} = sparse(permute(coil_compress_mtx_diagblk{x}, [2, 1])); % coil_compress_mtx_diagblk{i} Dim: [nc, nci]. coil_compress_mtx_diagblk{i} (i = 1,2,...nx) is the i-th diagonal block in the matrix coil_compress_mtx.
            end
            if echo == nec && slice == nsl
                clear ccmtx;
            end
            
            % -- Set up the whole coil compression matrix if we'll keep track of the full sample noise covariance matrix
            if strcmp(propagate_type, 'cov')
                coil_compress_mtx = blkdiag(coil_compress_mtx_diagblk{:}); % coil_compress_mtx Dim: [CoilAfterCompression->X(=nc*nx), CoilBeforeCompression->X(=nci*nx)]
                clear coil_compress_mtx_diagblk;
            end
        end
        
        % -- Case 2 continued: whiten_type is 'coil_noise_std' or 'none'. After data whitening step, psi_mtx is not an identity matrix.
        if ismember(whiten_type, {'coil_noise_std', 'none'})
            switch propagate_type
                case 'cov' % Keep track of the full sample noise covariance matrix, psi_mtx_tilde
                    psi_mtx_tilde = spalloc(nsamp*nc*nx, nsamp*nc*nx, nc*nx*nc*nx*nsamp); % psi_mtx_tilde Dim: [Sample->Coil->X(=nsamp*nc*nx), Sample->Coil->X(=nsamp*nc*nx)]
                case 'diagblk'
                    psi_mtx_diagblk = cell(nx);
            end
            
            % - Loop over samples. Each sample(i.e. ky line) undergoes a different phase filtering, so need to keep track of the noise covariance for each sample
            for samp = 1 : nsamp
                % Default Nyquist ghosting correction, step 2: pha_flt correction
                pha_flt = get_pha_flt(pha_coe(:, samp, slice, :), nxi); % Dim: [X(=nxi), Ky/Sample(=1), Slice(=1), Coil(=nci)]
                pha_flt = sparse(1 : nxi*nci, 1 : nxi*nci, pha_flt(:), nxi*nci, nxi*nci, nxi*nci); % Dim: [X->Coil(=nxi*nci), X->Coil(=nxi*nci)]
                psi_mtx_tilde_samp = cov_linear_transform(pha_flt, psi_mtx); % psi_mtx_tilde_samp Dim: [X->Coil(=nxi*nci), X->Coil(=nxi*nci)]
                
                % Default Nyquist ghosting correction, step 3: 1D DFT along x
                psi_mtx_tilde_samp = cov_linear_transform(fdftmtx, psi_mtx_tilde_samp); % Dim: [Kx->Coil(=nxi*nci), Kx->Coil(=nxi*nci)]
                
                % Ramp sampling correction
                psi_mtx_tilde_samp = cov_linear_transform(ramp_flt_mtx, psi_mtx_tilde_samp); % Dim: [Kx->Coil(=nx*nci), Kx->Coil(=nx*nci)]
                
                % 1D iDFT along kx
                psi_mtx_tilde_samp = cov_linear_transform(idftmtx, psi_mtx_tilde_samp);% Dim: [X->Coil(=nx*nci), X->Coil(=nx*nci)]
                
                % Change index in psi_mtx_tilde_samp from an ordering of nx*nci to nci*nx
                psi_mtx_tilde_samp = psi_mtx_tilde_samp(idx_permute_nx_nci, idx_permute_nx_nci); % Dim: [Coil->X(=nci*nx), Coil->X(=nci*nx)]
                
                switch propagate_type
                    case 'cov'
                        % - Coil compression
                        if coil_compress
                            psi_mtx_tilde_samp = cov_linear_transform(coil_compress_mtx, psi_mtx_tilde_samp); % Dim: [CoilAfterCompression->X(=nc*nx), CoilAfterCompression->X(=nc*nx)]
                        end
                        
                        % - Write noise covariance of this sample into the full sample noise covariance matrix
                        idx = samp : nsamp : nsamp*nc*nx;
                        psi_mtx_tilde(idx, idx) = psi_mtx_tilde_samp; % psi_mtx_tilde Dim: [Sample->Coil->X(=nsamp*nc*nx), Sample->Coil->X(=nsamp*nc*nx)]
                        
                    case 'diagblk'
                        for x = 1 : nx
                            % - Write noise covariance of this sample into the cell array for different x positions
                            idx = (x-1)*nci+1 : x*nci;
                            psi_mtx_diagblk{x}{samp} = psi_mtx_tilde_samp(idx, idx); % psi_mtx_diagblk{x}{samp} Dim: [Coil(=nci), Coil(=nci)]. psi_mtx_diagblk{x}{samp} is the x-th nci-by-nci diagonal block in the matrix psi_mtx_tilde_samp.
                            
                            % - Coil compression
                            if coil_compress
                                psi_mtx_diagblk{x}{samp} = cov_linear_transform(coil_compress_mtx_diagblk{x}, psi_mtx_diagblk{x}{samp}); % psi_mtx_diagblk{x}{samp} Dim: [Coil(=nc), Coil(=nc)]
                            end
                        end
                end
            end

            switch propagate_type
                case 'cov'
                    clear coil_compress_mtx;
                case 'diagblk'
                    clear coil_comrress_mtx_diagblk;
            end
            if echo == nec && slice == nsl
                clear fdftmtx ramp_flt_mtx idftmtx idx_permute_nx_nci pha_coe;
            end
            
        else % Case 1 continued: whiten_type is 'psi_mtx', different samples on the ky-omegaz plane are I.I.D.. psi_mtx Dim: [Coil->X(nci*nx), Coil->X(nci*nx)]
            
            if strcmp(propagate_type, 'cov')
                % - Coil compression
                if coil_compress
                    psi_mtx_tilde = cov_linear_transform(coil_compress_mtx, psi_mtx); % Dim: [CoilAfterCompression->X(=nc*nx), CoilAfterCompression->X(=nc*nx)]
                end
                clear coil_compress_mtx;
                
                % - Set up sample noise covariance matrix for all samples on the ky-omeagaz plane
                psi_mtx_tilde = kron(psi_mtx_tilde, speye(nsamp)); % psi_mtx_tilde Dim: [Sample->Coil->X(=nsamp*nc*nx), Sample->Coil->X(=nsamp*nc*nx)]
            end
        end
        
        if strcmp(propagate_type, 'cov') % psi_mtx_tilde Dim: [Sample->Coil->X(=nsamp*nc*nx), Sample->Coil->X(=nsamp*nc*nx)]
            switch p.add_vcc
                case p.ONLY_ACTUAL_COILS
                    % DO NOTHING
                case p.ACTUAL_AND_VIRTUAL_COILS
                    % TODO: Set up sample noise covariance for both actual and virtual coils
                case p.ONLY_VIRTUAL_COILS
                    psi_mtx_tilde = conj(psi_mtx_tilde); % For 2 random variables X, Y, cov(conj(X), conj(Y)) = conj(cov(X, Y))
            end
        end
        
        % -- Solve ky-kz and sensitivity encoding
        % The encoding matrices are set up in the same way as in the sense reconstruction routine mux_recon_sense.m (In mux_recon_sense.m, include_sense_recon controls whether SENSE reconstruction is included, whereas here SENSE reconstruction is always included).
        for msk_idx = 1 : nmsk
            if p.debug
                fprintf('      undersample scheme %d/%d\n', msk_idx, nmsk);
            end
            
            % - ftyz_mtx: The DFTy & FTz(DFTz or randomly sampled DTFTz) encoding parts of the encoding matrix
            ftyz_mtx = encode_ftyz_mtx(us_msk(msk_idx), ny, nz, p.add_vcc); % Dim: [SamplesOnKyOmegazPlane(=nsamp), Y->Z(=ny*nz)]
            ftyz_mtx = repmat(ftyz_mtx, [nc, p.KEEP_ORIG_SZ]); % Dim: [SamplesOnKyOmegazPlane->Coil(=nsamp*nc), Y->Z(=ny*nz)].
            
            % - For propagate_type 'cov': Will find the solving/encoding matrices for each x position then put them into a big matrix which includes all x positions
            if strcmp(propagate_type, 'cov')
                ftyz_and_smap_mtx_diagblk = cell(nx); % ftyz_and_smap_mtx_diagblk{i}(i=1,2,...nx) is the solving/encoding matrix for x=i.
            end
            
            for x = 1 : nx
                
                % Encoding matrix including (1)DFTy & DTFTz
                encode_mtx = ftyz_mtx; % Dim: [nsamp*nc, ny*nz]
                
                % The sensitivity part in the encoding matrix
                smap_mtx = smap(x, :, echo, slice, :, :, :); % Dim: [X(=1), Y(=ny), Echo(=1), Slice(=1), SimultaneousSliceZ(=nz), Coil(=nc), SetOfSensitivityMaps(=nmap)
                smap_mtx = permute(smap_mtx, [6,2,5,7,1,3,4]); % Dim: [Coil(=nc), Y(=ny), SimultaneousSliceZ(=nz), SetOfSensitivityMaps(=nmap), X(=1), Echo(=1), Slice(=1)]
                smap_mtx = reshape(smap_mtx, [1, nc, ny*nz*nmap]); % Dim: [1, nc, ny*nz*nmap]
                smap_mtx = repmat(smap_mtx, [nsamp, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ]); % Dim: [nsamp, nc, ny*nz*nmap]
                switch p.add_vcc
                    case p.ONLY_ACTUAL_COILS
                        % DO NOTHING
                    case p.ACTUAL_AND_VIRTUAL_COILS
                        smap_mtx = cat(1, smap_mtx, conj(smap_mtx)); % Dim: [nsamp*(1+extra_dat_vcc), nc, ny*nz*nmap]
                    case p.ONLY_VIRTUAL_COILS
                        smap_mtx = conj(smap_mtx);
                end
                smap_mtx = reshape(smap_mtx, [nsamp*nc, ny*nz*nmap]); % Dim: [nsamp*nc, ny*nz*nmap]
                
                % Encoding matrix including (1)DFTy & DTFTz and (2)sensitivity
                encode_mtx = repmat(encode_mtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nmap]); % Dim: [nsamp*nc, ny*nz, nmap]
                encode_mtx = reshape(encode_mtx, [nsamp*nc, ny*nz*nmap]); % Dim: [nsamp*nc, ny*nz*nmap]
                encode_mtx = encode_mtx .* smap_mtx; % Dim: [nsamp*nc, ny*nz*nmap]
                
                if strcmp(im_noise_type, 'actual')
                    solve_mtx = get_solve_mtx(encode_mtx, p.sense_lambda, p.sense_lambda_default_ratio); % Dim: [ny*nz*nmap, nsamp*nc]
                    clear encode_mtx;
                end
                
                switch propagate_type
                    case 'diagblk' % Keep track of the nx diagonal nc-by-nc blocks in the full sample niose covariance matrix. In this case im_noise_type must be 'actual'.
                        
                        % Set up sample noise covariance matrix for all samples on the ky-omegaz plane that correspond to this x position
                        if ismember(whiten_type, {'coil_noise_std', 'none'}) % Case 2 continued
                            psi_mtx_tilde_diagblk = psi_mtx_diagblk{x}; % psi_mtx_tilde_diagblk{samp} Dim: [Coil(=nc), Coil(=nc)]. psi_mtx_diagblk{x}{samp}(x=1,2...nx; samp=1,2...nsamp) is the x-th nc-by-nc diagonal block in the sample noise covariance matrix for this sample on the ky-omegaz plane (i.e. psi_mtx_tilde_samp)
                            psi_mtx_diagblk{x} = []; % Set to empty to reduce memory use
                            
                            psi_mtx_tilde_diagblk = blkdiag(psi_mtx_tilde_diagblk{:}); % Dim: [Coil->Sample(=nc*nsamp), Coil->Sample(=nc*nsamp)]
                            
                            idx_permute_nc_nsamp = zeros(1, nc*nsamp);
                            for coil = 1 : nc
                                idx_permute_nc_nsamp((coil-1)*nsamp+1 : 1 : coil*nsamp) = coil : nc : nc*nsamp;
                            end
                            psi_mtx_tilde_diagblk = psi_mtx_tilde_diagblk(idx_permute_nc_nsamp, idx_permute_nc_nsamp); % Dim: [Sample->Coil(=nsamp*nc), Sample->Coil(=nsamp*nc)]
                        
                        else % Case 1 (whiten_type is 'psi_mtx') continued
                            psi_mtx_tilde_diagblk = psi_mtx_diagblk{x}; % psi_mtx_tilde_diagblk Dim: [Coil(=nci), Coil(=nci)]
                            if coil_compress
                                psi_mtx_tilde_diagblk = cov_linear_transform(coil_compress_mtx_diagblk{x}, psi_mtx_tilde_diagblk); % Dim: [Coil(=nc), Coil(=nc)]
                            end

                            psi_mtx_tilde_diagblk = kron(psi_mtx_tilde_diagblk, speye(nsamp)); % Dim: [Sample->Coil(nsamp*nc), Sample->Coil(=nsamp*nc)]
                        end
                        
                        switch p.add_vcc
                            case p.ONLY_ACTUAL_COILS
                                % DO NOTHING
                            case p.ACTUAL_AND_VIRTUAL_COILS
                                % TODO: Set up sample noise covariance for both actual and virtual coils
                            case p.ONLY_VIRTUAL_COILS
                                psi_mtx_tilde_diagblk = conj(psi_mtx_tilde_diagblk);
                        end
                        
                        % Solve sensitivity and ky-kz encoding
                        psi_mtx_tilde_diagblk = cov_linear_transform(solve_mtx, psi_mtx_tilde_diagblk); % Dim: [ny*nz*nmap, ny*nz*nmap]

                        % Combine images that correspond to multiple sensitivity maps
                        if nmap > 1
                            psi_mtx_tilde_diagblk = cov_linear_transform(combine_mtx_diagblk, psi_mtx_tilde_diagblk); % Dim: [ny*nz, ny*nz]
                        end

                        % Save results into image noise matrix
                        psi_mtx_tilde_diagblk = diag(psi_mtx_tilde_diagblk); % Image noise variance vector for this x position. Dim: [ny*nz, 1]
                        im_noise(x, :, echo, slice, :, msk_idx) = real(sqrt(reshape(full(psi_mtx_tilde_diagblk), [1, ny, 1, 1, nz, 1]))); % Image noise standard deviation for this x position

                    case 'cov'
                        switch im_noise_type
                            case 'actual' % Calculate the noise propagation in the actual reconstruction
                                ftyz_and_smap_mtx_diagblk{x} = sparse(solve_mtx); % ftyz_and_smap_mtx_diagblk{x} Dim: [ny*nz*nmap, nsamp*nc]

                            case 'snr_optimal' % Calculate the noise propagation of an SNR-optimal reconstruction
                                ftyz_and_smap_mtx_diagblk{x} = sparse(encode_mtx); % ftyz_and_smap_mtx_diagblk{x} Dim: [nsamp*nc, ny*nz*nmap]
                        end
                end
                
            end % x
            
            if strcmp(propagate_type, 'cov')
                % - Image noise covariance matrix
                switch im_noise_type
                    case 'actual'
                        solve_mtx = blkdiag(ftyz_and_smap_mtx_diagblk{:}); % Dim: [ny*nz*nmap*nx, nsamp*nc*nx]
                        psi_mtx_tilde = cov_linear_transform(solve_mtx, psi_mtx_tilde); % Image noise covariance matrix. Dim: [Y->Z->SetOfSensitivityMaps->X(ny*nz*nmap*nx), Y->Z->SetOfSensitivityMaps->X(ny*nz*nmap*nx)]. Didn't use Eq. [17] in Pruessmann's SENSE paper since we're not using the SNR-optimal SENSE solver here.
                        clear solve_mtx;
                        
                    case 'snr_optimal'
                        encode_mtx = blkdiag(ftyz_and_smap_mtx_diagblk{:}); % Dim: [nsamp*nc*nx, ny*nz*nmap*nx]
                        psi_mtx_tilde = inv(encode_mtx' * (psi_mtx_tilde \ encode_mtx)); % Image noise covariance matrix. Dim: [Y->Z->SetOfSensitivityMaps->X(ny*nz*nmap*nx), Y->Z->SetOfSensitivityMaps->X(ny*nz*nmap*nx)]. This corresponds to Eq.[17] in Pruessmann's SENSE paper(E: encode_mtx, psi_tilde: psi_mtx_tilde, X=inv(E'*inv(psi_tilde)*E)).
                        clear encode_mtx;
                end
                
                % - Combine images that correspond to multiple sensitivity maps
                if nmap > 1
                    combine_mtx = kron(speye(nx), combine_mtx_diagblk); % Dim: [Y->Z->X(=ny*nz*nx), Y->Z->SetOfSensitivityMaps->X(=ny*nz*nmap*nx)]
                    psi_mtx_tilde = cov_linear_transform(combine_mtx, psi_mtx_tilde); % Image noise covariance matrix. Dim: [Y->Z->X(=ny*nz*nx), Y->Z->X(=ny*nz*nx)].
                    clear combine_mtx;
                end
                
                % - Final image noise
                psi_mtx_tilde = diag(psi_mtx_tilde); % Image noise variance vector. Dim: [Y->Z->X(ny*nz*nx), 1]
                psi_mtx_tilde = permute(reshape(full(psi_mtx_tilde), [ny, nz, nx]), [3,1,4,5,2]); % Dim: [X(=nx), Y(=ny), Echo(=1), Slice(=1), Z(=nz), MaskType(=1)]
                im_noise(:, :, echo, slice, :, msk_idx) = real(sqrt(psi_mtx_tilde)); % Image noise standard deviation
            end % strcmp(propagate_type, 'cov')
            
        end % msk_idx
    end % slice
end % echo

im_noise = fftshift(im_noise, 5); % Dim: [X(=nx), Y(=ny), Echo(=nec=1), Slice(=nsl), SimultaneousSliceZ(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), MaskType(=nmsk)].
im_noise = reshape(im_noise, [nx, ny, nec, nsl*nz, nmsk]); % Dim: [X(=nx), Y(=ny), Echo(=nec=1), Slice*SimultaneousSliceZ(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), MaskType(=nmsk)].

return

function idftmtx = get_idftmtx(n)
%
% Generates a matrix for inverse DFT transform corresponding to function ifftc.m
%
% Input
%   n       - Length of the inverse DFT transform.
%
% Output
%   idftmtx - Matrix representation of the inverse DFT transform which corresponds to function ifftc.m. Dim: [n, n]
%

idftmtx = fftshift(fftshift(conj(dftmtx(n)), 2), 1) / sqrt(n);

return

function fdftmtx = get_fdftmtx(n)
%
% Generates a matrix for forward DFT transform corresponding to function fftc.m
%
% Input
%   n - Length of forward DFT transform.
%
% Output
%   fdftmtx - Matrix representation of the forward DFT transform which corresponds to function fftc.m. Dim: [n, n]
%

fdftmtx = fftshift(fftshift(dftmtx(n), 2), 1) / sqrt(n);

return

function sigma_AX = cov_linear_transform(A, sigma_X)
%
% Calculates covariance matrix of random vector AX, where X is a random vector
% with covariance matrix sigma_X. sigma_AX = A * sigma_X * A'.
%
% Inputs
%   A        - Matrix for the linear transform.
%   sigma_X  - Covariance matrix of random vector X.
%
% Output
%   sigma_AX - Covariance matrix of random vector AX. This is always a sparse matrix.
%
% Kangrong Zhu  Stanford University     March 2015

if ~issparse(A)
    A = sparse(A); % Make sure all matrices used in the multiplication are sparse so that the multiplication result will be sparse
end
if ~issparse(sigma_X)
    sigma_X = sparse(sigma_X);
end

sigma_AX = A * sigma_X * A';

return

function solve_mtx = get_solve_mtx(encode_mtx, lambda, lambda_default_ratio)
% Calculates the matrix for solving an encoding.

AtA = encode_mtx' * encode_mtx;
if ~exist('lambda', 'var') || isempty(lambda)
    lambda = norm(AtA, 'fro') / size(encode_mtx,2) * lambda_default_ratio; % Regularization parameter, from Miki Lustig's code
end

solve_mtx = pinv(AtA + eye(size(AtA))*lambda) * encode_mtx'; % Use pinv instead of inv so that no warning messages come up when cropped sensitivity maps are used

return