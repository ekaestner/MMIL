function [sense, oneD_grappa, slice_grappa, split_slice_grappa, pmr_cal, dat_mux, ccmtx, smap, noise, im_sing] = pmr_sim_mux_seq(p, us_msk, dat_cal, eddy, noise, im_sing)
%
% function [sense, oneD_grappa, slice_grappa, split_slice_grappa, pmr_cal, dat_mux, ccmtx, smap, noise, im_sing] = pmr_sim_mux_seq(p, us_msk, dat_cal, [eddy=[]], [noise], [im_sing])
%
% Pseudo Multiple Replica: Simulate slice-multiplexed acquisition.
% Reference: Philip M. Robson, et al. MRM 60:895-907 (2008).
%
% Inputs
%   p        - Parameter structure. See mux_epi_params.m for details.
%              The following fields are particularly needed for pseudo-multiple replica simulation:
%              p.psi_mtx                - Coil noise covariance matrix. Dim: [Coil(=nc), Coil(=nc)].
%              p.pmr_num_rep            - Number of replicas to simulate for calculating the noise in the reconstructed images. Default: 50.
%              In addition to fields defined in mux_epi_params.m, here p also contains fields:
%              p.sim_sense              - True: Simulate SENSE reconstruction.
%              p.sim_1d_grappa          - True: Simulate 1D GRAPPA reconstruction.
%              p.sim_slice_grappa       - True: Simulate Slice-GRAPPA reconstruction.
%              p.sim_split_slice_grappa - True: Simulate Split-Slice-GRAPPA reconstruction.
%              p.pmr_dat_ecc            - Raw ECC data. Input to function epi_process_rawdata. Needed only when p.cap_get_ecc is true. See epi_process_rawdata.m for details.
%              p.pmr_ref                - Path to 'ref.dat'. Input to function epi_process_rawdata. See epi_process_rawdata.m for details.
%              p.pmr_vrgf               - Path to 'vrgf.dat'. Input to function epi_process_rawdata. See epi_process_rawdata.m for details.
%              p.pmr_refp               - Path to reference scan pfile. Input to function epi_process_rawdata. See epi_process_rawdata.m for details.
%              p.pmr_slices             - Desired slices. Input to function epi_process_rawdata. See epi_process_rawdata.m for details.
%              p.pmr_slices_ref         - Slice indices in the reference scan. Input to function epi_process_rawdata. See epi_process_rawdata.m for details.
%              p.pmr_sim_cal            - True: Simulate multiple replica of the calibration acquisition.
%              p.pmr_recon_sing         - True: Reconstruct original single-slice images and compare them with reconstructed slice-multiplexed images.
%              p.ccmtx                  - Coil compression matrices. Used only when p.coil_compress is true, by function coil_compress. If not exist or is empty, will compute from calibration data.
%              p.smap                   - Sensitivity maps. Used only when p.sim_sense is true. If not exist or is empty, will be calculated from calibration data.
%   us_msk   - Undersampling mask on the ky-omegaz plane. When partial ky is used, length of us_msk.ky and us_msk.kz corresponds to data matrix size without zero padding for partial ky. See get_ky_omegaz_us_msk.m for details.
%   dat_cal  - Raw k-space calibration data. Dim: [Kx(=nxi, initial nx before ramp sampling correction), Ky(=nky, not zero-padded for partial ky), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nkz*p.num_mux_cycle)].
%   eddy     - If nonempty: Contains phase terms the eddy current effects cause, and eddy current effects will be included in the encoding process.
%                xKyOmegaz-data-with-eddy-current = xKyOmegaz-data-without-eddy-current .* eddy.
%                If funtcion 'epi_pha_correct' is used to calculate the phase terms, 'eddy' should be reshaped from the output 'pha_flt'.
%                Dim: [X(=p.nx_pres), PE(=nsamp), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))), Coil(=nc)].
%              If empty: eddy current effects will not be included.
%   noise, im_sing - These are included in the input list to enable use of precalculated values. Their dimensions are the same as their corresponding outputs.
%              noise    - Noise to add to the synthesized slice-multiplexed data.
%              im_sing  - SOS-coil-combined reconstructed original single-slice images. Used only when p.pmr_recon_sing is true. Default: Will be calculated from the input single-slice data.
%
% Outputs
%   sense/oneD_grappa/slice_grappa/split_slice_grappa - Structures holding pseudo-multiple replica simulation results for each undersampling mask when SENSE/1D-GRAPPA/Slice-GRAPPA/Split-Slice-GRAPPA reconstruction is used. Nonempty only when p.sim_sense/p.sim_1d_grappa/p.sim_slice_grappa/p.sim_split_slice_grappa is true.
%              Take 'sense' as an example. We have length(sense) = length(us_msk).
%              sense(i)         - Results corresponding to the i-th undersampling mask.
%              sense(i).snr     - SNR maps. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%              sense(i).noise   - Noise maps. Dim:  [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%              sense(i).signal  - The reconstructed magnitude image when no noise is added. Nonempty only when p.debug is true. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=1), Time(=1)].
%              sense(i).replica - The reconstructed magnitude images with added noise. Nonempty only when p.debug is true. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=1), Time(=p.pmr_num_rep)].
%              sense(i).rrms    - Relative root mean squared error of the 1st time point in sense(i).replica, with respect to the original single-slice images. Exists only when p.pmr_recon_sing is true. Dim: [1, 1, Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%              sense(i).diff    - Difference images between the 1st time point in sense(i).replica and the original single-slice images. Exists only when p.pmr_recon_sing is true. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%   pmr_cal  - A structure holding the pseudo-multiple replica simulation results for the calibration data. Nonempty only when p.pmr_sim_cal is true.
%              pmr_cal.snr, pmr_cal.noise, pmr_cal.signal, pmr_cal.replica are analogous to and having the same dimensions as the corresponding fields in the output 'sense'.
%   dat_mux  - Slice-multiplexed data synthesized from the input calibration data. Nonempty only when p.debug is true. Fields include:
%              dat_mux(i).dat   - Corresponding to the i-th undersampling mask. The 1st time point contains no added noise, the 2nd time point contains added noise that conforms to the input coil noise covariance matrix. Dim: [Kx(=p.nx_pres), PE(=nsamp), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=2)].
%   ccmtx, smap, noise, im_sing - Output these so that the same values may be used for future simulation.
%              ccmtx    - Coil compression matrices. Dim: [Coil(=p.num_coils, before compression), Coil(=p.num_vcoils, after compression), X(=1 if p.cc_method is 'single'; =p.nx_pres if p.cc_method is 'geometric'), Echo(=nec), Slice(=nsl)]
%              smap     - Sensitivity maps. Nonempty only when p.sim_sense is true. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil(=p.num_coils, after coil compression when applicable), SetOfSensitivityMaps(=p.espirit_nmap if p.smap_type is 'espirit'; =1 if p.smap_type is 'coil_over_sos')].
%              noise    - Noise added to the synthesized slice-multiplexed data. Dim: [Kx(=nxi), PE(=nsamp), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=1+p.pmr_num_rep)].
%              im_sing  - SOS-coil-combined reconstructed original single-slice images. Nonempty only when p.pmr_recon_sing is true. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%
% (c) Kangrong Zhu      Stanford University     Nov 2014

if ~isfield(p, 'pmr_num_rep') || isempty(p.pmr_num_rep)
    p.pmr_num_rep = 50;
end

if ~exist('dat_cal', 'var') || isempty(dat_cal)
    error('Must input calibration data.');
end

if ~exist('eddy', 'var')
    eddy = [];
end

OUTPUT_IDX_CCMTX = 7; % Output index for the coil compression matrices
OUTPUT_IDX_SMAP = 8; % Output index for the sensitivity maps
NT_TO_SIM_SIGNAL = 1; % Number of time points to use to simulate the signal (reconstruct with no added noise)
SPECIFIED_OMEGAZ_ENC = 0;

[nxi, nky, nec, nsl, nc, nt_cal] = size(dat_cal);
nz = p.mux;
nsamp = length(us_msk(1).kz); % Number of samples on the ky-omegaz plane
nt = NT_TO_SIM_SIGNAL + p.pmr_num_rep; % Use the time dimension to simulate the signal and all the replicas. nt is the total number of time points to synthesize.

%% -- Calibration data
nkz = p.mux;
num_mux_cycle_orig = nt_cal / nkz;
dat_cal = dat_cal(:, :, :, :, :, nkz*(num_mux_cycle_orig-1)+1 : nkz*num_mux_cycle_orig); % Keep only 1 set of calibration data to make it easier to synthesize multiple replica of calibration acquisition
p.num_mux_cycle = 1;

%% -- Original single-slice raw k-space data
dat_sing = mux_dftz(dat_cal, p.T_DIM, p.cap_fov_shift_cal, nz, 'decode');

%% -- Original single-slice images
if p.pmr_recon_sing
    if ~exist('im_sing', 'var') || isempty(im_sing)
        if p.debug
            fprintf(' Reconstructing original single-slice data used in pseudo-multiple replica simulation...\n');
        end
        
        % -- Match the EPI data processing steps
        im_sing = epi_process_rawdata(dat_sing, p); % Dim: [Kx(=p.nx_pres), Ky(=nky), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nz)]
        
        % -- Match partial k recon
        if p.partial_ky
            im_sing = cat(p.PE_DIM, im_sing, zeros(p.nx_pres, p.ny_pres-nky, nec, nsl, nc, nz));
            im_sing = homodyne(im_sing, nky, p.homodyne_ntran, p.homodyne_niter, 'kxky', 'xy');
        else
            im_sing = ifft2c(im_sing);
        end
        
        % -- Reshape and coil-combination
        im_sing = permute(im_sing, [1,2,3,4,6,5]); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil]
        im_sing = reshape(fftshift(im_sing, 5), [p.nx_pres, p.ny_pres, nec, nsl*nz, nc]); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=nc)]
        im_sing = sos(im_sing, p.C_DIM); % Coil-combined single-slice images. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))]
    end
else
    im_sing = [];
end

%% -- Simulate multiple replica of the calibration acquisition
if p.pmr_sim_cal
    if p.debug
        fprintf(' Simulating multiple replica of the calibration acquisition...\n');
    end

    % -- Form replicas of the calibration data
    dat_cal_reps = zeros(nxi, nky, nec, nsl, nc, nkz, nt);
    dat_cal_reps(:, :, :, :, :, :, NT_TO_SIM_SIGNAL) = dat_cal;
    for time = NT_TO_SIM_SIGNAL+1 : nt % Don't add noise to the 1st time point which is used to simulate the signal
        dat_cal_reps(:, :, :, :, :, :, time) = dat_cal + pmr_syn_noise(p.psi_mtx, [nxi, nky, nec, nsl, nc, nkz], p.C_DIM);
    end
    
    % -- Match the EPI data processing steps
    dat_cal_reps = reshape(dat_cal_reps, [nxi, nky, nec, nsl, nc, nkz*nt]); % Form a time series corresponding to p.num_mux_cycle=1
    dat_cal_reps = epi_process_rawdata(dat_cal_reps, p); % Dim: [Kx(=p.nx_pres), Ky(=nky), Echo(=nec), Slice(=nsl), Coil(=nc), Kz->Replica(=nkz*nt)]
    
    % -- Match coil compression
    if p.coil_compress
        [dat_cal_reps, p] = coil_compress(dat_cal_reps, p); % p.num_coils will be set to p.num_vcoils after coil compression
    end
    
    % -- Solve DFT encoding
    dat_cal_reps = reshape(dat_cal_reps, [p.nx_pres, nky, nec, nsl, p.num_coils, nkz, nt]);
    dat_cal_reps = mux_dftz(dat_cal_reps, p.T_DIM, p.cap_fov_shift_cal, nz, 'decode'); % Dim: [Kx, Ky, Echo, Slice, Coil, SimultaneousSlice Z, Time(i.e. Replica)]
    
    % -- Change to image space. Combine coil images using SENSE if needed. Refer to mux_epi_process_data_grappa.m and mux_epi_process_data_slice_grappa.m
    dat_cal_reps = permute(dat_cal_reps, [1,2,3,4,6,5,7]); % Dim: [Kx(=p.nx_pres), Ky(=nky), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))), Coil(=p.num_coils), Time(i.e. Replica, =nt)]
    dat_cal_reps = reshape(fftshift(dat_cal_reps, 5), [p.nx_pres, nky, nec, nsl*nz, p.num_coils, nt]); % Dim: [Kx(=p.nx_pres), Ky(=nky), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=p.num_coils), Time(i.e. Replica, =nt)]
    if ~p.return_sos_im && p.sense_coil_comb
        dat_cal_reps = sense_coil_comb(dat_cal_reps, p); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=1), Time(i.e. Replica, =nt)]
    else
        if p.partial_ky
            dat_cal_reps = cat(p.PE_DIM, dat_cal_reps, zeros(p.nx_pres, p.ny_pres-nky, nec, nsl*nz, p.num_coils, nt));
        end
        dat_cal_reps = ifft2c(dat_cal_reps);
    end
    
    % -- Match partial k recon
    if p.partial_ky
        dat_cal_reps = homodyne(dat_cal_reps, nky, p.homodyne_ntran, p.homodyne_niter, 'xy', 'xy');
    end

    % -- Zero pad reconstructed images
    if p.zpad_image
        dat_cal_reps = ifft2c(zpad(fft2c(dat_cal_reps), [p.zpad_size(1), p.zpad_size(2), nec, nsl*nz, size(dat_cal_reps, p.C_DIM), nt]));
    end
    
    % -- Final images
    if p.return_sos_im
        dat_cal_reps = sos(dat_cal_reps, p.C_DIM); % Coil-combined calibration images. Dim: [X(=p.nx_pres if p.zpad_image is false; =p.zpad_size(1) if p.zpad_image is true), Y(=p.ny_pres if p.zpad_image is false; =p.zpad_size(2) if p.zpad_image is true), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=1), Time(i.e. Replica, =nt)]
    end

    % -- Pseudo-multiple replica SNR
    [pmr_cal.signal, pmr_cal.replica] = separate_signal_and_replica(dat_cal_reps, NT_TO_SIM_SIGNAL);
    [pmr_cal.snr, pmr_cal.noise] = pmr_calc_snr(pmr_cal.signal, pmr_cal.replica, p.T_DIM);
    
    % -- For outputing results
    if ~p.debug
        pmr_cal.signal = [];
        pmr_cal.replica = [];
    end

    % -- Reduce memory use
    clear('dat_cal_reps');
else
    pmr_cal = [];
end

%% -- Simulate multiple replica of the slice-multiplexed acquisition
if p.debug
    fprintf(' Simulating multiple replica of the mux acquisition...\n');
end

% -- Noise to add to synthesized slice-multiplexed data. Calculate this up front so that the same noise will be added whening simulating different undersampling patterns.
if ~exist('noise', 'var') || isempty(noise)
    noise = zeros([nxi, nsamp, nec, nsl, nc, nt]);
    for time = NT_TO_SIM_SIGNAL+1 : nt % Don't add noise to the 1st time point which is used to simulate the signal
        noise(:, :, :, :, :, time) = pmr_syn_noise(p.psi_mtx, [nxi, nsamp, nec, nsl, nc], p.C_DIM);
    end
end

% -- Initialize result structures
nmsk = length(us_msk);
if p.debug
    dat_mux(nmsk).dat = [];
else
    dat_mux = [];
end
sense(nmsk).snr = [];
oneD_grappa(nmsk).snr = [];
slice_grappa(nmsk).snr = [];
split_slice_grappa(nmsk).snr = [];

% -- Simulate undersampling mask by undersampling mask
for msk_idx = 1 : nmsk
    if p.debug
        fprintf(' Simulating undersampling scheme %d/%d...\n', msk_idx, nmsk);
    end
    msk = us_msk(msk_idx);
    if p.sim_sense || p.sim_slice_grappa || p.sim_split_slice_grappa
        p.pmr_us_msk = msk;
    end
    
    % - Synthesize slice-multiplexed data
    if p.partial_ky && (msk_idx == 1)
        dat_sing = cat(p.PE_DIM, dat_sing, zeros(nxi, p.ny_pres-nky, nec, nsl, nc, nz));
    end
    mux_dat = mux_encode(dat_sing, msk, nt, noise, eddy); % Input dat_sing Dim: [Kx(=nxi), Ky(=p.ny_pres, zero-padded for partial ky), Echo(=nec), Slice(=nsl), Coil(=nc), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))]. Output mux_dat Dim: [Kx(=nxi), Ky(=nky, not zero-padded for partial ky), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt)]
    
    % - Put together calibration and slice-multiplexed data to get the synthesized time series
    dat_tseries = zeros(nxi, nky, nec, nsl, nc, nkz*p.num_mux_cycle+nt);
    dat_tseries(:, :, :, :, :, 1:nkz*p.num_mux_cycle) = dat_cal;
    dat_tseries(:, msk.ky, :, :, :, nkz*p.num_mux_cycle+1:end) = mux_dat;

    % - Reduce memory use
    clear('mux_dat');

    % - Match the EPI data processing steps
    p.num_coils = nc; % Set to original number of coils. p.num_coils might have been modified after coil compression when simulating multiple replica of the calibration acquisition. Also, when msk_idx > 1, it has been modified after coil compression in the previous loop.
    dat_tseries = epi_process_rawdata(dat_tseries, p); % Dim: [Kx(=p.nx_pres), Ky(=nky), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nkz*p.num_mux_cycle+nt)]

    % - Output slice-multiplexed data
    if p.debug
        dat_mux(msk_idx).dat = dat_tseries(:, :, :, :, :, nkz*p.num_mux_cycle+(1:NT_TO_SIM_SIGNAL + min(1, p.pmr_num_rep)));
    end

    % - Process data
    % Coil compression
    if p.coil_compress
        [dat_tseries, p] = coil_compress(dat_tseries, p);
        if nargout >= OUTPUT_IDX_CCMTX
            ccmtx = p.ccmtx;
            p = rmfield(p, 'ccmtx');
        end
    else
        if nargout >= OUTPUT_IDX_CCMTX
            ccmtx = [];
        end
    end
    
    % SENSE recon
    if p.sim_sense
        % Parameter
        p.mux_recon_method = 'sense'; % Directly set fields of p, to avoid copying the parameter structure
        p_orig.num_mux_cycle = p.num_mux_cycle; % Save a few fields that might be changed to p_orig so that they can be set back later
        
        % Sensitivity maps
        if isfield(p, 'smap') && ~isempty(p.smap)
            p.num_mux_cycle = 0; % Already have sensitivity maps, won't pass calibration data any more
        end

        % Recon the slice-multiplexed time series.
        [dat_rec, p] = mux_epi_process_data_sense(dat_tseries(:, :, :, :, :, (nkz*p_orig.num_mux_cycle*(p.num_mux_cycle==0)) + 1:end), p); % If input p.smap doesn''t exist or is empty, p.smap will be set to the sensitivity maps used in the recon.
        
        % Output sensitivity maps, if needed
        if nargout >= OUTPUT_IDX_SMAP
            smap = p.smap;
        end
        p = rmfield(p, 'smap'); % Reduce memmory use
        
        % Pseudo-multiple replica SNR
        [sense(msk_idx).signal, sense(msk_idx).replica] = separate_signal_and_replica(dat_rec(:, :, :, :, :, p.num_mux_cycle+1:end), NT_TO_SIM_SIGNAL);
        [sense(msk_idx).snr, sense(msk_idx).noise] = pmr_calc_snr(sense(msk_idx).signal, sense(msk_idx).replica, p.T_DIM);
        
        % Difference w.r.t original single-slice images
        if p.pmr_recon_sing
            [sense(msk_idx).rrms, sense(msk_idx).diff] = calc_rrms_diff(im_sing, abs(sense(msk_idx).replica(:, :, :, :, :, 1)));
        end
        
        % If param structure p will be used again, set back related fields
        if p.sim_1d_grappa || p.sim_slice_grappa || p.sim_split_slice_grappa || (nmsk > 1)
            p = copy_fields(p, p_orig);
        end
        
        % For output and reduce memory use
        if ~p.debug
            sense(msk_idx).signal = [];
            sense(msk_idx).replica = [];
        end
        
     else
        if nargout >= OUTPUT_IDX_SMAP
            smap = [];
        end
    end
    
    % 1D GRAPPA recon
    if p.sim_1d_grappa && (~p.use_gzblips || (p.use_gzblips && p.caipi)) % 1D GRAPPA only works with no-blip or CAIPI
        % Parameter. Directly set fields of p to avoid unnecessary copy of p. If p will be used again by (Split)-Slice-GRAPPA, its related fields will be reset.
        p.mux_recon_method = '1Dgrappa';
        p.cal_dat_tpoints = p.mux*(p.num_mux_cycle-1)+1 : p.mux*p.num_mux_cycle;
        p.pseq = 1;
        p.zpad_1Dgrappa = 'minZpad';
        if p.caipi && (abs(p.cap_fov_shift) ~= p.mux) && strcmp(p.zpad_1Dgrappa, 'minZpad');
            p.zpad_1Dgrappa = 'normalZpad';
        end
        
        % Recon
        dat_rec = mux_epi_process_data_grappa(dat_tseries, p);
        
        % Pseudo-multiple replica SNR
        [oneD_grappa(msk_idx).signal, oneD_grappa(msk_idx).replica] = separate_signal_and_replica(dat_rec(:, :, :, :, :, p.num_mux_cycle+1:end), NT_TO_SIM_SIGNAL);
        [oneD_grappa(msk_idx).snr, oneD_grappa(msk_idx).noise] = pmr_calc_snr(oneD_grappa(msk_idx).signal, oneD_grappa(msk_idx).replica, p.T_DIM);
        
        % Difference w.r.t original single-slice images
        if p.pmr_recon_sing
            [oneD_grappa(msk_idx).rrms, oneD_grappa(msk_idx).diff] = calc_rrms_diff(im_sing, abs(oneD_grappa(msk_idx).replica(:, :, :, :, :, 1)));
        end
        
        % For output and reduce memory use
        if ~p.debug
            oneD_grappa(msk_idx).signal = [];
            oneD_grappa(msk_idx).replica = [];
        end
    end
    
    % Slice-GRAPPA recon
    % Common parameters. Directly set fields of p to avoid unnecessary copy of p.
    if p.sim_slice_grappa || p.sim_split_slice_grappa
        p.mux_recon_method = 'slice-grappa';
        
        if p.partial_ky
            len = length(p.pmr_us_msk.ky);
            if len < p.ny_pres/p.inplane_R % pad undersampling mask size to p.ny_pres/p.inplane_R, with arbitrary ky,kz values. The undersampling mask size needs to be with partial ky when synthesizing mux data, but needs to be without partial ky for (Split-)Slice-GRAPPA recon.
                p.pmr_us_msk.ky = cat(2, p.pmr_us_msk.ky(:).', ones(1, p.ny_pres/p.inplane_R - len));
                p.pmr_us_msk.kz = cat(2, p.pmr_us_msk.kz(:).', ones(1, p.ny_pres/p.inplane_R - len));
                p.pmr_us_msk = encode_ftz_pha(p.pmr_us_msk, nz, SPECIFIED_OMEGAZ_ENC);
            end
        end
    end
    if p.sim_slice_grappa
        % Parameter
        p.use_split_slice_grappa = false;
        
        % Recon
        dat_rec = mux_epi_process_data_slice_grappa(dat_tseries, p);
        
        % Pseudo-multiple replica SNR
        [slice_grappa(msk_idx).signal, slice_grappa(msk_idx).replica] = separate_signal_and_replica(dat_rec(:, :, :, :, :, p.num_mux_cycle+1:end), NT_TO_SIM_SIGNAL);
        [slice_grappa(msk_idx).snr, slice_grappa(msk_idx).noise] = pmr_calc_snr(slice_grappa(msk_idx).signal, slice_grappa(msk_idx).replica, p.T_DIM);
        
        % Difference w.r.t original single-slice images
        if p.pmr_recon_sing
            [slice_grappa(msk_idx).rrms, slice_grappa(msk_idx).diff] = calc_rrms_diff(im_sing, abs(slice_grappa(msk_idx).replica(:, :, :, :, :, 1)));
        end
        
        % For output and reduce memory use
        if ~p.debug
            slice_grappa(msk_idx).signal = [];
            slice_grappa(msk_idx).replica = [];
        end
    end
    
    % Split-Slice-GRAPPA recon
    if p.sim_split_slice_grappa
        % Parameter
        p.use_split_slice_grappa = true;
        
        % Recon
        dat_rec = mux_epi_process_data_slice_grappa(dat_tseries, p);
        
        % Pseudo-multiple replica SNR
        [split_slice_grappa(msk_idx).signal, split_slice_grappa(msk_idx).replica] = separate_signal_and_replica(dat_rec(:, :, :, :, :, p.num_mux_cycle+1:end), NT_TO_SIM_SIGNAL);
        [split_slice_grappa(msk_idx).snr, split_slice_grappa(msk_idx).noise] = pmr_calc_snr(split_slice_grappa(msk_idx).signal, split_slice_grappa(msk_idx).replica, p.T_DIM);
        
        % Difference w.r.t original single-slice images
        if p.pmr_recon_sing
            [split_slice_grappa(msk_idx).rrms, split_slice_grappa(msk_idx).diff] = calc_rrms_diff(im_sing, abs(split_slice_grappa(msk_idx).replica(:, :, :, :, :, 1)));
        end
        
        % For output and reduce memory use
        if ~p.debug
            split_slice_grappa(msk_idx).signal = [];
            split_slice_grappa(msk_idx).replica = [];
        end
    end
end

return

function [sig, replica] = separate_signal_and_replica(dat, nt_sig)
%
% Seperate the signal and replicas from a reconstructed time series.
%

C_DIM = 5;
sig = sos(sos(dat(:, :, :, :, :, nt_sig, :), 7), C_DIM); % 7-th dimension: SetOfSensitivityMaps for SENSE; 1 for 1D-GRAPPA & Slice-GRAPPA & Split-Slice-GRAPPA
replica = sos(sos(dat(:, :, :, :, :, nt_sig+1:end, :), 7), C_DIM);

return

function [rrms_vals, diff] = calc_rrms_diff(Iref, Irecon)
%
% Calculates difference images and relative root mean squared error (RRMS) between reference images (Iref) and reconstructed images (Irecon).
%

rrms_vals = rrms(Iref, Irecon);
diff = abs(Iref - Irecon);

return

function p = copy_fields(p, copy_from)

if ~exist('copy_from', 'var') || isempty(copy_from)
    return;
end

fields = fieldnames(copy_from);
for idx = 1 : numel(fields)
    p.(fields{idx}) = copy_from.(fields{idx});
end

return