function [dat_out, ranks] = mux_recon_sense(dat_mux, us_msk, p, smap, eddy)
%
% function [dat_out, ranks] = mux_recon_sense(dat_mux, us_msk, p, [smap=[]], [eddy=[]])
% 
% Decode up to 3 types of encoding in simultaneous multislice imaging: 
%   (1) DFTy and FTz(DFTz for CAIPI, DTFTz for MICA) encoding.
%   (2)(optional) Sensitivity encoding(Can take multiple sets of sensitivity maps).
%   (3)(optional) eddy current effects induced phase. 
% For each group of multiplexed slices, the decoding is conducted x position by x position.
%
% Inputs
%   dat_mux - Slice-multiplexed k-space data. Dim: [Kx(=nx), Ky(=nsamp),
%             Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt)].
%   us_msk  - Undersample mask on the ky-omegaz plane. A structure with
%             fields 'ky', 'kz', 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   p       - Parameter structure. See mux_epi_params.m for details. The
%             following fields are used in this function: debug, mux, FE_DIM(=1),
%             PE_DIM(=2), EC_DIM(=3), SL_DIM(=4), C_DIM(=5), T_DIM(=6), KZ_DIM(=7),
%             KEEP_ORIG_SZ(=1), ny_pres, sense_lambda, sense_lambda_default_ratio,
%             add_vcc, ONLY_ACTUAL_COILS, ACTUAL_AND_VIRTUAL_COILS, ONLY_VIRTUAL_COILS.
%   smap    - If nonempty: Contains sensitivity maps. SENSE reconstruction
%                will be included in the decoding process. Dim: [X(=nx), Y(=ny), Echo(=1),
%                Slice(=nsl), SimultaneousSlice Z(=nz,z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil(=nc), SetOfSensitivityMaps(=nmap)];
%             If empty: SENSE reconstruction will not be included.
%   eddy    - If nonempty: Contains phase terms the eddy current effects cause.
%                Eddy current effect correction will be included in the decoding process.
%                xKyOmegaz-data-with-eddy-current = xKyOmegaz-data-without-eddy-current .* eddy.
%                If funtcion 'epi_pha_correct' is used to calculate the phase terms,
%                'eddy' should be reshaped from the output 'pha_flt'.
%                Dim: [X(=nx), PE(=nsamp), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))), Coil(=nc)].
%             If empty: eddy current effect correction will not be included.
%
% Outputs
%   dat_out - Reconstructed multislice images. Dim: [X(=nx), Y(=ny), Echo(=nec),
%             Slice*SimultaneousSlice(=nsl*nz), Coil(If SENSE recon included, =1;
%             If SENSE recon not included, =nc), Time(=nt),
%             SetOfSensitivityMap(=nmap if SENSE recon included, =1 if SENSE recon not included)].
%   ranks   - The rank of all encoding matrices which were (pseudo-)inverted
%             in the decoding process. Dim: [X(=nx), Echo(=nec), Slice(=nsl),
%             Coil(If SENSE recon included, =1; If SENSE recon not included, =nc), Undersampling mask].
%             Calculated only when p.debug==true. If p.debug==false, ranks = [].
%
% (c) Kangrong Zhu,     Stanford University    Aug 2013

if length(us_msk(1).ky) ~= length(us_msk(1).kz)
    error('Length of ''us_msk.ky'' does not match length of ''us_msk.kz''.');
end

if ~exist('smap', 'var')
    smap = [];
end

if ~exist('eddy', 'var')
    eddy = [];
end

ONE       = 1;                                               % For specifying matrix size of SENSE reconstructed data
NT_THRESH = 300;                                             % If nt > NT_THRESH, solve time point by time point; otherwise solve all time points together

calc_rank = p.debug;                                         % True: Calculate and record the rank of each encoding matrix
include_sense_recon = (~isempty(smap));                      % True: Include SENSE reconstruction in the decoding process
include_eddy_correct = (~isempty(eddy));                     % True: Include the eddy current effect correction in the decoding process

datsz = get_dat_sz(dat_mux, p);

nc    = datsz.c;                                             % Define some variables to make the later code easier to read
nt    = datsz.t;
ny    = p.ny_pres;
nz    = p.mux;
nsamp = length(us_msk(1).ky);                                % Number of acquired points on the ky-omegaz plane
nmsk  = length(us_msk);                                      % Number of undersample masks for the ky-omegaz plane. nmsk > 1: The undersample scheme changes with time

if ~include_sense_recon && (nsamp < ny*nz)                   % When SENSE recon is not included, must have sufficient number of sampled points on the ky-omegaz plane
    error('SENSE recon is not included, too few sampled points on the ky-omegaz plane.');
end
if nmsk == 0
    error('Undersample mask empty.');
end
if (nmsk > 1) && (nmsk ~= nt);
    error('Number of undersample masks doesn''t match number of time points');
end

dat_mux_xkyomegaz = ifftc(dat_mux, p.FE_DIM);                 % Slice multiplexed data in x-ky-omegaz space. Dim: [X(=nx), Ky(=nsamp), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt)]

if include_sense_recon                                        % If SENSE recon included in the decoding processs
    nc_rec = ONE;                                             % Number of coils in the reconstructed data is 1
    nmap = size(smap, 7);                                     % Number of sets of sensitivity maps
else
    nc_rec = nc;
    nmap = 1;
    p.add_vcc = p.ONLY_ACTUAL_COILS;                          % Virtual coil concept is of no use if SENSE recon is not included, so always turn it off
end

if p.debug
    fprintf('   Including data from ');
    switch p.add_vcc
        case p.ONLY_ACTUAL_COILS
            fprintf('only actual coils...\n');
        case p.ACTUAL_AND_VIRTUAL_COILS
            fprintf('both actual and virtual coils...\n');
        case p.ONLY_VIRTUAL_COILS
            fprintf('only virtual coils...\n');
    end
end

switch p.add_vcc
    case p.ONLY_ACTUAL_COILS
        extra_dat_vcc = 0;                                    % Number of additional sets of samples, for adding virtual coil concept
    case p.ACTUAL_AND_VIRTUAL_COILS
        extra_dat_vcc = 1;
    case p.ONLY_VIRTUAL_COILS
        extra_dat_vcc = 0;
end

dat_out = zeros(datsz.x, ny, datsz.ec, datsz.sl, nz, nc_rec, nt, nmap); % Dim: [X(=nx), Y(=ny), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz), Coil(If SENSE recon included, =1; If SENSE recon not included, =nc), Time(=nt), SetOfSensitivityMaps].
if calc_rank
    ranks = zeros(datsz.x, datsz.ec, datsz.sl, nc_rec, nmsk); % Dim: [X(=nx), Echo(=nec), Slice(=nsl), Coil(If SENSE recon included, =1; If SENSE recon not included, =nc), Undersampling masks]
else
    ranks = [];                                               % Avoid unnecessary computing when not in debug mode.
end

for msk_idx = 1 : nmsk
    if p.debug && mod(msk_idx-1, 1) == 0
        fprintf('   Calculating undersample scheme %d/%d\n', msk_idx, nmsk);
    end
    
    % ftyz_mtx: The DFTy & FTz(DFTz or randomly sampled DTFTz) encoding parts of the encoding matrix
    % Dim: If SENSE recon included: [nsamp*(1+extra_dat_vcc)*nc, ny*nz]. Row order: encoding index(acquisition order on the ky-omegaz plane) -> actual/virtual coils -> coil index; Column order: y -> z.
    %      If SENSE recon not included: [nsamp, ny*nz] (p.add_vcc and extra_dat_vcc are always 0). Row order: encoding index(acquisition order on the ky-omegaz plane); Column order: y -> z.
    ftyz_mtx = encode_ftyz_mtx(us_msk(msk_idx), ny, nz, p.add_vcc); % Dim: [nsamp*(1+extra_dat_vcc), ny*nz]
    if include_sense_recon
        ftyz_mtx = repmat(ftyz_mtx, [nc, p.KEEP_ORIG_SZ]); % Dim: [nsamp*(1+extra_dat_vcc)*nc, ny*nz]
    else
        ftyz_mtx = repmat(ftyz_mtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nc]); % When SENSE recon not included, Dim: [nsamp, ny*nz, nc]
    end
    
    for echo = 1 : datsz.ec
        for slice = 1 : datsz.sl
            for x = 1 : datsz.x
                
                if p.debug && mod(x-1, 40) == 0
                    fprintf('    x = %d/%d\n', x, datsz.x);
                end
                
                % Encoding matrix including (1)DFTy & DTFTz
                encode_mtx = ftyz_mtx; % Dim: If SENSE recon included: [nsamp*(1+extra_dat_vcc)*nc, ny*nz]; If SENSE recon not included: [nsamp, ny*nz, nc].
                
                % The eddy current effect part in the encoding matrix
                if include_eddy_correct
                    eddy_mtx = eddy(x, :, echo, slice, :, :); % Dim: [1, nsamp, 1, 1, nz, nc]
                    eddy_mtx = reshape(eddy_mtx, [nsamp, 1, nz, nc]); % Dim: [nsamp, 1, nz, nc]
                    eddy_mtx = repmat(eddy_mtx, [p.KEEP_ORIG_SZ, ny, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ]);
                    eddy_mtx = reshape(eddy_mtx, [nsamp, ny*nz, nc]); % Dim: [nsamp, ny*nz, nc]
                    if include_sense_recon
                        switch p.add_vcc
                            case p.ONLY_ACTUAL_COILS
                                % DO NOTHING
                            case p.ACTUAL_AND_VIRTUAL_COILS
                                eddy_mtx = cat(1, eddy_mtx, conj(eddy_mtx)); % Dim: [nsamp*(1+extra_dat_vcc), ny*nz, nc]
                            case p.ONLY_VIRTUAL_COILS
                                eddy_mtx = conj(eddy_mtx);
                        end
                    end
                    
                    % Encoding matrix including (1) DFTy & DTFTz and (2) eddy current effects
                    if include_sense_recon
                        eddy_mtx = permute(eddy_mtx, [1,3,2]); % Dim: [nsamp*(1+extra_dat_vcc), nc, ny*nz]
                        eddy_mtx = reshape(eddy_mtx, [nsamp*(1+extra_dat_vcc)*nc, ny*nz]); % Dim: [nsamp*(1+extra_dat_vcc)*nc, ny*nz]
                    end
                    encode_mtx = encode_mtx .* eddy_mtx; % Dim: If SENSE recon included: [nsamp*(1+extra_dat_vcc)*nc, ny*nz]; If SENSE recon not included: [nsamp, ny*nz, nc]
                end
                
                % The sensitivity part in the encoding matrix
                if include_sense_recon
                    smap_mtx = smap(x, :, echo, slice, :, :, :); % Dim: [X(=1), Y(=ny), Echo(=1), Slice(=1), SimultaneousSliceZ(=nz), Coil(=nc), SetOfSensitivityMaps(=nmap)
                    smap_mtx = permute(smap_mtx, [6,2,5,7,1,3,4]); % Dim: [Coil(=nc), Y(=ny), SimultaneousSliceZ(=nz), SetOfSensitivityMaps(=nmap), X(=1), Echo(=1), Slice(=1)]
                    smap_mtx = reshape(smap_mtx, [1, nc, ny*nz*nmap]); % Dim: [1, nc, ny*nz*nmap]
                    unknown_pixels = zeros(1, ny*nz*nmap); % Dim: [1, ny*nz*nmap]. 1 indicates a pixel is not in the background; 0 indicates a pixel is in the background and no sensitivity estimation exists for it.
                    for pixel_idx = 1 : ny*nz*nmap
                        if any(smap_mtx(1, :, pixel_idx) ~= 0)
                            unknown_pixels(pixel_idx) = 1;
                        end
                    end
                    unknown_pixels = find(unknown_pixels == 1); % Indices for the unknown pixels
                    num_unknown_pixels = length(unknown_pixels);
                    smap_mtx = smap_mtx(1, :, unknown_pixels); % Excludes background pixels from the reconstruction. Will directly set them to 0 later.
                    smap_mtx = repmat(smap_mtx, [nsamp, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ]); % Dim: [nsamp, nc, num_unknown_pixels]
                    switch p.add_vcc
                        case p.ONLY_ACTUAL_COILS
                            % DO NOTHING
                        case p.ACTUAL_AND_VIRTUAL_COILS
                            smap_mtx = cat(1, smap_mtx, conj(smap_mtx)); % Dim: [nsamp*(1+extra_dat_vcc), nc, num_unknown_pixels]
                        case p.ONLY_VIRTUAL_COILS
                            smap_mtx = conj(smap_mtx);
                    end
                    smap_mtx = reshape(smap_mtx, [nsamp*(1+extra_dat_vcc)*nc, num_unknown_pixels]); % Dim: [nsamp*(1+extra_dat_vcc)*nc, num_unknown_pixels]

                    % Encoding matrix including (1)DFTy & DTFTz (2) eddy current effects and (3)sensitivity
                    encode_mtx = repmat(encode_mtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nmap]); % Dim: [nsamp*(1+extra_dat_vcc)*nc, ny*nz, nmap]
                    encode_mtx = reshape(encode_mtx, [nsamp*(1+extra_dat_vcc)*nc, ny*nz*nmap]); % Dim: [nsamp*(1+extra_dat_vcc)*nc, ny*nz*nmap]
                    encode_mtx = encode_mtx(:, unknown_pixels); % Dim: [nsamp*(1+extra_dat_vcc)*nc, num_unknown_pixels]
                    encode_mtx = encode_mtx .* smap_mtx; % Dim: [nsamp*(1+extra_dat_vcc)*nc, num_unknown_pixels]
                end
                
                % Calculate the pseudo-inverse of the encoding matrices
                if include_sense_recon
                    solve_mtx = get_solve_mtx(encode_mtx, p.sense_lambda, p.sense_lambda_default_ratio); % Dim: [num_unknown_pixels, nsamp*(1+extra_dat_vcc)*nc]
                    if calc_rank
                        ranks(x, echo, slice, ONE, msk_idx) = rank(encode_mtx);
                    end
                else
                    solve_mtx = zeros(ny*nz, nsamp, nc); % Dim: [ny*nz, nsamp, nc]
                    for coil = 1 : nc
                        solve_mtx(:, :, coil) = get_solve_mtx(encode_mtx(:, :, coil), p.sense_lambda, p.sense_lambda_default_ratio);
                        if calc_rank
                            ranks(x, echo, slice, coil, msk_idx) = rank(encode_mtx(:, :, coil));
                        end
                    end
                end
                
                % Measured signal
                if nmsk == 1 % Undersample mask doesn't change with time
                    time_for_this_us_msk = 1 : nt; % Time points corresponding to the current undersample mask
                else % Undersample mask changes with time
                    time_for_this_us_msk = msk_idx;
                end
                nt_to_recon = length(time_for_this_us_msk); % Number of time points to reconstruct for this undersample mask. =nt if undersample scheme doesn't change with time; =1 if undersample scheme changes with time
                measured_sig = reshape(dat_mux_xkyomegaz(x, :, echo, slice, :, time_for_this_us_msk), [nsamp, nc, nt_to_recon]); % Dim: [nsamp, nc, nt_to_recon]
                
                % Solve images
                if (nmsk == 1) && (nt_to_recon <= NT_THRESH) % Undersample scheme doesn't change with time(nt_to_recon=nt) AND nt_to_recon no larger than NT_THRESH, reconstruct all time points together
                    if include_sense_recon % If SENSE recon included
                        switch p.add_vcc
                            case p.ONLY_ACTUAL_COILS
                                % DO NOTHING
                            case p.ACTUAL_AND_VIRTUAL_COILS
                                measured_sig = cat(1, measured_sig, conj(measured_sig)); % Dim: [nsamp*(1+extra_dat_vcc), nc, nt_to_recon]
                            case p.ONLY_VIRTUAL_COILS
                                measured_sig = conj(measured_sig);
                        end
                        measured_sig = reshape(measured_sig, [nsamp*(1+extra_dat_vcc)*nc, nt_to_recon]); % Dim: [nsamp*(1+extra_dat_vcc)*nc, nt_to_recon]
                        reconed_sig = solve_mtx * measured_sig; % Dim: [num_unknown_pixels, nt_to_recon]
                        if num_unknown_pixels < ny * nz * nmap
                            tmp = zeros(ny * nz * nmap, nt_to_recon);
                            tmp(unknown_pixels, :) = reconed_sig;
                            reconed_sig = tmp; % Dim: [ny*nz*nmap, nt_to_recon]
                        end
                        reconed_sig = permute(reshape(reconed_sig, [ny, nz, nmap, nt_to_recon]), [5,1,6,7,2,8,4,3]); % Dim: [X(=1), Y(=ny), Echo(=1), Slice(=1), SimultaneousSlice Z(=nz), Coil(=1), Time(=nt_to_recon), SetOfSensitivityMap(=nmap)]
                        dat_out(x, :, echo, slice, :, ONE, :, :) = reconed_sig; % Dim: [X(=nx), Y(=ny), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz), Coil(=1), Time(=nt), SetOfSensitivityMap(=nmap)]
                    else % If SENSE recon not included
                        for coil = 1 : nc
                            reconed_sig = solve_mtx(:, :, coil) * reshape(measured_sig(:, coil, :), [nsamp, nt_to_recon]); % Dim: [ny*nz, nt_to_recon]
                            reconed_sig = permute(reshape(reconed_sig, [ny, nz, nt_to_recon]), [4,1,5,6,2,7,3]); % Dim: [X(=1), Y(=ny), Echo(=1), Slice(=1), SimultaneousSlice Z(=nz), Coil(=1), Time(=nt_to_recon)]
                            dat_out(x, :, echo, slice, :, coil, :) = reconed_sig; % Dim: [X(=nx), Y(=ny), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz), Coil(=nc), Time(=nt)]
                        end
                    end
                else % Undersample scheme changes with time OR nt_to_recon larger than NT_THRESH
                    for time_idx = 1 : nt_to_recon
                        time = time_for_this_us_msk(time_idx);
                        sig = measured_sig(:, :, time_idx); % Dim: [nsamp, nc]
                        if include_sense_recon % If SENSE recon included
                            switch p.add_vcc
                                case p.ONLY_ACTUAL_COILS
                                    % DO NOTHING
                                case p.ACTUAL_AND_VIRTUAL_COILS
                                    sig = cat(1, sig, conj(sig)); % Dim: [nsamp*(1+extra_dat_vcc), nc]
                                case p.ONLY_VIRTUAL_COILS
                                    sig = conj(sig);
                            end
                            reconed_sig = solve_mtx * reshape(sig, [nsamp*(1+extra_dat_vcc)*nc, 1]); % Dim: [num_unknown_pixels, Time(=1)]
                            if num_unknown_pixels < ny * nz * nmap
                                tmp = zeros(ny * nz * nmap, 1);
                                tmp(unknown_pixels) = reconed_sig;
                                reconed_sig = tmp; % Dim: [ny*nz*nmap, Time(=1)]
                            end
                            reconed_sig = permute(reshape(reconed_sig, [ny, nz, nmap]), [4,1,5,6,2,7,8,3]); % Dim: [X(=1), Y(=ny), Echo(=1), Slice(=1), SimultaneousSlice Z(=nz), Coil(=1), Time(=1), SetOfSensitivityMaps(=nmap)]
                            dat_out(x, :, echo, slice, :, ONE, time, :) = reconed_sig; % Dim: [X(=nx), Y(=ny), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz), Coil(=1), Time(=nt), SetOfSensitivityMaps(=nmap)]
                        else % If SENSE recon not included
                            for coil = 1 : nc
                                reconed_sig = solve_mtx(:, :, coil) * sig(:, coil);% Dim: [ny*nz, Time(=1)]
                                reconed_sig = permute(reshape(reconed_sig, [ny, nz]), [3,1,4,5,2,6,7]); % Dim: [X(=1), Y(=ny), Echo(=1), Slice(=1), SimutlaneousSlice Z(=nz), Coil(=1), Time(=1)]
                                dat_out(x, :, echo, slice, :, coil, time) = reconed_sig; % Dim: [X(=nx), Y(=ny), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz), Coil(=nc), Time(=nt)]
                            end % coil
                        end % if include_sense_recon, else
                    end % time
                end % if (nmsk == 1) && (nt <= NT_THRESH), else
                
            end % x
        end % slice
    end % echo
end % msk_idx

% Output dat_out. Dim: [X(=nx), Y(=ny), Echo(=nec), Slice*SimultaneousSlice Z(=nsl*nz),
% Coil(If SENSE recon included, =1; If SENSE recon not included, =nc), Time(=nt)].
dat_out = fftshift(dat_out, 5); % Adjust slice ordering
dat_out = reshape(dat_out, [datsz.x, ny, datsz.ec, datsz.sl*nz, nc_rec, datsz.t, nmap]);

return

function solve_mtx = get_solve_mtx(encode_mtx, lambda, lambda_default_ratio)
% Calculates the matrix for solving an encoding.

AtA = encode_mtx' * encode_mtx;
if ~exist('lambda', 'var') || isempty(lambda)
    lambda = norm(AtA, 'fro') / size(encode_mtx,2) * lambda_default_ratio; % Regularization parameter, from Miki Lustig's code
end
% MJM - trying to overcome SVD error
%solve_mtx = pinv(AtA + eye(size(AtA))*lambda) * encode_mtx';               % Use pinv instead of inv so that no warning messages come up when cropped sensitivity maps are used
solve_mtx = (AtA + eye(size(AtA))*lambda)\encode_mtx';

return
