function [d,snr_d] = mb_mux_epi(pfile,outfile,ext_cal,slices,nt_to_recon,...
                                n_vcoils,debug,recon_method,apply_fermi,...
                                use_homodyne,notch_thresh)
%function [d,snr_d] = mb_mux_epi(pfile,outfile,ext_cal,slices,nt_to_recon,...
%                                n_vcoils,debug,recon_method,apply_fermi,...
%                                use_homodyne,notch_thresh)
%
% Purpose: slice-multiplexing EPI reconstruction
%
% Inputs
%   pfile       - Filename of the pfile, or the directory containing a pfile.
%   outfile     - The file for saving the results. If outfile is empty, then
%                 the results will not be saved to a file.
%   ext_cal     - External calibration. This can be
%                 (1) Filename of the pfile, or the directory containing
%                     a pfile. The corresponding external calibration scan
%                     can be either mux>1 or mux=1.
%                 (2) Calibration data matrix for all muxed slices.
%                     Dim: [FE, PE, Echo, MuxedSlice, Coil, Time(mux phase cycling)]
%                 (3) Empty.
%                 For external calibration, 'ext_cal' must not be empty.
%                 For internal calibration, if 'ext_cal'is empty internal
%                 calibration will be used, otherwise external calibration
%                 will be used.
%   slices      - Indices of the (muxed) slices to reconstruct. Default: All slices.
%   nt_to_recon - Number of time points to reconstruct, excluding the first
%                 few mux phase cycling time points. Not used if 'pfile' has
%                 sequential acquisition. Default: All time points.
%   n_vcoils    - Number of virtual coils for coil compression. No coil
%                 compression if n_vcoils is [] or 0 or the-number-of-physical-coils.
%   debug       - True: calculate intermediate results and print out messages.
%   recon_method- '1Dgrappa' or numbers other than 1: use 1D-GRAPPA;
%                 'sense' or number 1: use SENSE;
%                 'slice-grappa': use slice-GRAPPA;
%                 'split-slice-grappa': use split-slice-GRAPPA.
%                 Append '_sense1' to use sense coil combination. E.g., '1Dgrappa_sense1'.
%                 Default: 1Dgrappa.
%                 This option will be ignored for MICA-type acquisition (always uses SENSE).
%   apply_fermi - True: apply a circular Fermi filter using the rec.fermi parameters
%                 specified in the p-file header.
%   use_homodyne- True: use homodyne to fill partial-k acquisitions. False will use zero-filling.
%                 If use_homodyne is empty, then the p-file header will be used to determine whether
%                 to use homodyne or zero-filling.
%   notch_thresh- If non-zero, a denotching filter will be applied. Any
%                 point in k-space with a value lower than notch_thresh
%                 will be replaced by an adjacent time point.
%
% Outputs
%   d           - Reconstructed images. If p.return_sos_im == true, these are square-
%                 root-of-sum-of-squares (SOS) images; if p.return_sos_im == false,
%                 these are complex images. Dim: [1stImageDimension, 2ndImageDimension,
%                 Slice, Time, Coil(=nc_after_coil_compression, if (recon is GRAPPA-type &&
%                 p.return_sos_im == false && p.sense_coil_comb == false); =1, otherwise)].
%                 1stImageDimension - 2ndImageDimension:
%                             R/L - A/P (Axial and Axial-like Oblique)
%                             R/L - S/I (Coronal and Coronal-like Oblique)
%                             A/P - S/I (Sagittal and Sagittal-like Oblique)
%   snr_d       - SNR, noise, and retained SNR maps of the reconstructed images.
%                 Nonempty only when (p.calc_snr_amr == true || p.calc_snr_pmr ==true).
%                 snr_d.amr_snr   - SNR maps of the actually reconstructed SOS images,
%                                   corresponding to actual multiple replica method.
%                                   Dim: [1stImageDimension, 2ndImageDimension, Slice].
%                                   This is empty if p.calc_snr_amr == false, or if the acquisition is
%                                   a MICA-type whose undersampling mask changes with time,
%                                   or if the acquisition is a diffusion scan,
%                                   or if number of slice-multiplexed time points is smaller than p.amr_min_nt.
%                 snr_d.amr_noise - Noise maps of the actually reconstructed SOS images,
%                                   corresponding to actual multiple replica method.
%                                   Dimension same as snr_d.amr_snr. Nonempty only when snr_d.amr_snr is nonempty.
%                 snr_d.pmr_rsnr  - Retained SNR maps of the mux scan, from pseudo-multiple replica method.
%                                   Nonempty only when p.calc_snr_pmr == true.
%                                   Dim: [1stImageDimension,
%                                   2ndImageDimension, Slice,
%                                   UndersamplingMasks(=nmsk, number of unique undersampling masks used in the acquisition)].
%                 snr_d.pmr_noise - Noise maps of the mux scan, from pseudo-multiple replica method.
%                                   Nonempty only when p.calc_snr_pmr == true.
%                                   The noise maps may be scaled because the noise data in the noise.dat
%                                   file and the noise covariance matrix calculated from there are scaled.
%                                   Dim: [1stImageDimension, 2ndImageDimension, Slice, UndersamplingMasks(=nmsk)].
%                 snr_d.pmr_snr   - SNR maps of the SOS images, from pseudo-multiple replica method.
%                                   Nonempty only when (p.calc_snr_pmr == true && debug == true).
%                                   The SNR maps may be scaled.
%                                   Dim: [1stImageDimension, 2ndImageDimension, Slice, UndersamplingMasks(=nmsk)].
%                 snr_d.pmr_noise_ref   - Noise maps of the SOS reference images(either MUX1 images or slice-DFT-encoded calibration images), calculated by pseudo-multiple replica method.
%                                         Nonempty only when (p.calc_snr_pmr == true && debug == true).
%                                         Dim: [1stImageDimension, 2ndImageDimension, Slice].
%                 snr_d.sense_rsnr      - Retained SNR maps of the mux scan, from hybrid-space SENSE reconstruction model.
%                                         Nonempty only when (strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr).
%                                         Dim: [1stImageDimension, 2ndImageDimension, Slice, UndersamplingMasks(=nmsk)].
%                 snr_d.sense_noise     - Noise maps of the mux scan, from hybrid-space SENSE reconstruction model.
%                                         Nonempty only when (strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr).
%                                         The noise maps may be scaled.
%                                         Dim: [1stImageDimension, 2ndImageDimension, Slice, UndersamplingMasks(=nmsk)].
%                 snr_d.sense_noise_ref - SNR maps of the SENSE-reconstructed reference images(either MUX1 images or slice-DFT-encoded calibration images), calculated using the hybrid-space SENSE reconstruction model.
%                                         Nonempty only when (strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr && debug == true).
%                                         Dim: [1stImageDimension, 2ndImageDimension, Slice].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012
%
% Created:  08/01/12 by Kangrong Zhu (as mux_epi_main)
% Prev Mod: 02/03/16 by ?
% Prev Mod: 02/11/17 by Don Hagler
% Last Mod: 06/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters for the reconstruction.

% -- File directories and names.
if ~exist('pfile', 'var') || isempty(pfile)
    error('No pfile specified!');
end

mydir = fileparts(mfilename('fullpath'));
%addpath(fullfile(mydir, 'utils'), '-end');
%addpath(fullfile(mydir, 'pfile_reading'));
%addpath(fullfile(mydir, 'coil_compress'));
%addpath(fullfile(mydir, 'disp_and_sim_tools'));
%addpath(fullfile(mydir, 'correct_ref'));
%addpath(fullfile(mydir, 'correct_vrgf'));
%addpath(fullfile(mydir, 'espirit'));
%addpath(fullfile(mydir, 'pseudo_multiple_replica'));
%addpath(fullfile(mydir, 'poisson_disc_pattern'));
%addpath(fullfile(mydir, 'sense_rsnr'));

if(~exist('ext_cal', 'var'));      ext_cal     = [];     end
if(~exist('outfile', 'var'));      outfile     = [];     end
if(~exist('slices', 'var'));       slices      = [];     end
if(~exist('nt_to_recon', 'var'));  nt_to_recon = [];     end
if(~exist('n_vcoils', 'var'));     n_vcoils    = [];     end
if(~exist('debug', 'var'));        debug       = false;  end
if(~exist('recon_method', 'var')); recon_method= '1Dgrappa';end
if(~exist('apply_fermi', 'var'));  apply_fermi = false;  end
if(~exist('use_homodyne', 'var')); use_homodyne = true;  end
if(~exist('notch_thresh', 'var')); notch_thresh = 0;     end
if (~isnumeric(slices));           slices = str2double(slices);            end
if (~isnumeric(nt_to_recon));      nt_to_recon = str2double(nt_to_recon);  end
if (~isnumeric(n_vcoils));         n_vcoils = str2double(n_vcoils);        end
if (debug ~= 1)&&(debug ~= 0);     debug = str2double(debug);              end
if (~isnumeric(notch_thresh));     notch_thresh = str2double(notch_thresh);end
if (apply_fermi ~= 1)&&(apply_fermi ~= 0);     
    apply_fermi = str2double(apply_fermi);              
end
if (use_homodyne ~= 1)&&(use_homodyne ~= 0);     
    use_homodyne = str2double(use_homodyne);              
end

if debug && ~isempty(outfile)
    time_used = tic;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We do lots of ffts on arrays of the same size, so it's worth letting fftw
% measure the optimal algorithm the first time we run one.
% NOTE: this is broken on Octave versions < 3.6.4.
fftw('planner','measure');
if(exist('octave_config_info', 'builtin') && exist(fullfile(mydir,'fftw_wisdom.txt'), 'file'))
   tmp = load(fullfile(mydir,'fftw_wisdom.txt'));
   fftw('dwisdom',tmp.wisdom);
end

[pfile, ref, vrgf, refp, noise_fname, pfile_name] = get_pfile_name(pfile);

% -- Save all parameters of the mux epi scan and for the reconstruction in a structure.
if debug
    fprintf('Reading header from pfile %s...\n', pfile_name);
end
p = mux_epi_params(pfile, slices, nt_to_recon, n_vcoils, debug, recon_method, apply_fermi, use_homodyne, notch_thresh, noise_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.internal_cal && (~isempty(ext_cal))
    if p.debug
        fprintf('Will use external calibration data, although header indicates internal calibration data was acquired.\n');
    end
    p.internal_cal = false;
end

if (~p.internal_cal) && isempty(ext_cal)
    error('ext_cal must not be empty if external calibration is to be used.');
end

% -- External calibration
if ~p.internal_cal && ~isnumeric(ext_cal)                                % ~isnumeric(ext_cal): Input is pfile name or pfile directory, not calibration data matrix
    [ext_cal, ref_cal, vrgf_cal, refp_cal, noise_fname_cal] = get_pfile_name(ext_cal);

    % The ref.dat file might be missing or empty when external cal is used. If so, use the cal ref file.
    if ~exist(ref, 'file') || ~exist(vrgf, 'file')
        ref = ref_cal;
        vrgf = vrgf_cal;
    else                                                                 % it exists; make sure it's got data in it.
        d = dir(ref);
        if d.bytes < 64
            ref = ref_cal;
            vrgf = vrgf_cal;
        end
    end
end

% -- Check parameters
if ~exist(ref, 'file') && strcmp(p.ref_for_default_ecc, 'ref.dat')       % The ref.dat file is needed for all EPI scans
    error('Cannot find the ref.dat file!');
end

if ~exist(vrgf, 'file') && p.vrgf                                        % The vrgf.dat file is needed only when ramp sampling is on
    error('Cannot find the vrgf.dat file!');
end

if p.md_ecc
    if strcmp(p.ref_for_md_ecc, 'ref pfile') && ~exist(refp, 'file') && p.cap_get_ecc
        p.ref_for_md_ecc = 'ecc data';
        if p.debug
            fprintf('\nNo reference scan pfile. Will use ECC data collected by mux sequence in matrix-decoding ghosting correctrion.\n\n');
        end
    end
    if strcmp(p.ref_for_md_ecc, 'ecc data') && ~p.cap_get_ecc && exist(refp, 'file')
        p.ref_for_md_ecc = 'ref pfile';
        if p.debug
            fprintf('\nNo ECC data collected by mux sequence. Will use reference scan pfile in matrix-decoding ghosting correctrion.\n\n');
        end
    end
    if (strcmp(p.ref_for_md_ecc, 'ecc data') && ~p.cap_get_ecc) || (strcmp(p.ref_for_md_ecc, 'ref pfile') && ~exist(refp, 'file'))
        p.md_ecc = false;
        if p.debug
            fprintf('\nNo reference scan pfile and no ECC data collected by mux sequence, will not use matrix-decoding ghosting correction.\n\n');
        end
    end
end

if ~exist(refp, 'file') && (strcmp(p.ref_for_default_ecc, 'ref pfile') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')))
    error('Cannot find the reference scan pfile!');
end

if p.debug
    fprintf('Scan Parameters:\n');
    fprintf(' MUX excited %d; ARC %d; Number of mux cycling %d; Slice gap %gmm; Matrix [%d, %d];\n', p.mux_excited, p.inplane_R, p.num_mux_cycle, p.sldist, p.nx_pres, p.ny_pres);
    fprintf(' RF TBW=%g, pw=%gms, flip=%d;\n', p.multiband_TBW, p.multiband_pw, p.flip_angle);
    fprintf(' Acquired %d slices. Slice acquisition order (%s):', p.num_slices, p.acq_order); disp(p.sl_acq_order.');
    fprintf(' Acquired %d time points. ', p.num_passes);
    switch p.acq_order
        case 'interleaved'
            fprintf('Number of time points after the %d ecc and %d mux phase cycling to reconstruct is %d.\n', p.cap_get_ecc*p.mux_excited, p.mux_encoded*p.num_mux_cycle, p.nt_to_recon);
        case 'sequential'
            fprintf('Number of time points to reconstruct is 1, because each TR acquired one slice and equivalently the data contains only 1 time point.\n');
    end
    fprintf(' Internal calibration?\t %d;\n Gz blip on?\t %d;\n', p.internal_cal, p.use_gzblips);
    fprintf(' MUX encoded %d; MUX factor in recon will be %d;\n', p.mux_encoded, p.mux);
    if p.use_gzblips && (p.mux_excited > 1)
        if p.caipi
            fprintf(' CAIPI acquisition, CAIPI FOV shift for accelerated time points is %d;\n', p.cap_fov_shift);
        end
        if p.mica_br
            fprintf(' Bit-reversed MICA acquisition, temporal seed shift = %d;\n', p.cap_seed_shift);
        end
        if p.mica_rand
            fprintf(' Random MICA acquisition, temporal seed shift = %d;\n', p.cap_seed_shift);
        end
        if p.mica_perturbed_caipi
            fprintf(' Perburbed CAIPI MICA acquisition, temporal seed shift = %d, kz perturb range = %.4f;\n', p.cap_seed_shift, p.cap_kz_rand_pert_range);
        end
        if p.mica_poisson
            fprintf(' Poisson-disc MICA acquisition;\n');
        end
    end
    fprintf(' Ramp sampling on?\t %d;\n Partial ky?\t %d;\n descend_acq?\t %d;\n ychop?\t %d;\n\n', p.vrgf, p.partial_ky, p.descend_acq, p.ychop);

    if(p.rds_data_order)
        fprintf('RDS was used to acquire data; assuming rds data ordering for the raw data\n');
    end
end

% -- External calibration
if ~p.internal_cal && ~isnumeric(ext_cal)
    p_cal = mux_epi_params(ext_cal, slices, [], n_vcoils, debug, recon_method, apply_fermi, use_homodyne, notch_thresh, noise_fname_cal);

    % Sometimes the coil_noise_std in the target scan isn't set correctly.
    % In that case, use the values from the cal scan.
    if (~isfield(p, 'coil_noise_std') || isempty(p.coil_noise_std) || any(isnan(p.coil_noise_std))) && isfield(p_cal, 'coil_noise_std')
        p.coil_noise_std = p_cal.coil_noise_std;
    end

    if (p_cal.mux>1 && p_cal.num_mux_cycle <= 0); error('Ext Cal: num_mux_cycle <= 0.\n'); end
    if (p_cal.mux>1 && p_cal.sldist ~= p.sldist); error('Ext Cal: Band separation mismatch.\n'); end
    if (p_cal.nx_pres ~= p.nx_pres); error('Ext Cal: Prescribed size in FE mismatch.\n'); end
    if (p_cal.ny_pres ~= p.ny_pres); error('Ext Cal: Prescribed size in PE mismatch.\n'); end
    if (p_cal.inplane_R ~= p.inplane_R); disp('Ext Cal: Inplane acceleration factor mismatch.\n'); end
    if (p_cal.mux>1 && p_cal.multiband_TBW ~= p.multiband_TBW); disp('Ext Cal: TBW of RF pulse mismatch.\n'); end
    if (p_cal.mux>1 && p_cal.multiband_pw ~= p.multiband_pw); disp('Ext Cal: Pulse duration of RF pulse mismatch.\n'); end

    switch p_cal.mux_excited
        case 1                                                           % MUX=1 external calibration
            nt_to_recon_cal = [];                                        % Load all time points (including the first p_cal.num_mux_cycle time points)
        case p.mux_excited                                               % MUX=p.mux_excited external calibration
            nt_to_recon_cal = 0;                                         % Load all mux cycling data (don't load accelerated data)
            if p_cal.num_passes < (p_cal.cap_get_ecc*p_cal.mux_excited + p_cal.mux*p_cal.num_mux_cycle) % For external calibration whose actual number of time points is smaller than the number of time points needed for the specified num_mux_cycle
                p_cal.num_mux_cycle = floor((p_cal.num_passes-p_cal.cap_get_ecc*p_cal.mux_excited)/p_cal.mux); % Keep only fully sampled mux cyclings
            end
        otherwise
            error('Ext Cal: MUX factor must be either 1 or the same as the actual scan.\n');
    end
    p_cal = mux_epi_params_set_tpoints(p_cal, nt_to_recon_cal);

    % If external calibration has longer DFT encoding, need to update a few parameter fields
    if (p_cal.mux_excited == p.mux_excited) && (p_cal.mux_encoded > p.mux_encoded)
        p.recon_cal_use_mux_extra = p_cal.recon_cal_use_mux_extra;
        if p.mux ~= p_cal.mux
            p.mux = p_cal.mux;
            p.num_vcoils = p_cal.num_vcoils;
            p.num_unmuxed_slices = p_cal.num_unmuxed_slices;
            if strcmp(p.mux_recon_method, '1Dgrappa')
                p.mux_acssz_1d_grappa.y = p_cal.mux_acssz_1d_grappa.y;
                p.zpad_1Dgrappa = p_cal.zpad_1Dgrappa;
            end
        end
    end
end

% Slices(Normally ordered indices, not pfile indices) in reference scan corresponding to the slices(Normally ordered indices, not pfile indices) to reconstruct in actual mux scan
if strcmp(p.ref_for_default_ecc, 'ref pfile') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit && strcmp(p.ref_for_md_ecc, 'ref pfile'))
    if ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'})
        if p.slice_grappa_odd_even_fit
            md_ecc = true;
        else
            md_ecc = false;
        end
    else
        md_ecc = p.md_ecc;
    end
    ref_slices = get_ref_slices(refp, p.slices_to_recon, p.num_slices, p.mux, md_ecc, p.debug);
end

%% Call the reconstruction routines.
% We have to actually allocate this below, becuase we don't yet know the timeseries length.
d = [];

if strcmp(p.acq_order, 'sequential')                                     % Sequential acquisition is only used for inversion recovery measurement with mux_epi2(spin echo mux) sequence.
    if p.debug
        fprintf('Loading data from p-file...\n');
    end
    p_seq = p;
    p_seq = mux_epi_params_set_tpoints(p_seq, []);
    [dat_all, p] = epi_load_tseries(pfile, ref, vrgf, refp, [], [], [], p_seq); % Load all time points, i.e. load all slices (Each time point acquires one slice in this case)
end

for sl = p.slices_to_recon(:)'

    if p.debug
        fprintf('Processing slice %d/%d...\n', find(p.slices_to_recon == sl), length(p.slices_to_recon));
    end

    % -- Slice indices in reference scan
    if strcmp(p.ref_for_default_ecc, 'ref pfile') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit && strcmp(p.ref_for_md_ecc, 'ref pfile'))
        sl_idx = find(p.slices_to_recon == sl);
        sl_ref = ref_slices(sl_idx : length(p.slices_to_recon) : end);
    else
        sl_ref = [];
    end

    % -- Load the EPI time series data (Already corrected for EPI ghosting and ramp sampling).
    switch p.acq_order
        case 'interleaved'
            if p.debug
                fprintf(' Loading data from p-file...\n');
            end
            sl_in_pfile = p.sl_acq_order(sl);
            [dat, p_recon, dat_ecc] = epi_load_tseries(pfile, ref, vrgf, refp, sl, sl_in_pfile, sl_ref, p); % p_recon will be used as the parameter structure in the reconstruction. Several fields may have been changed in function epi_load_tseries. Parameters many change later if external calibration or coil compression is used. Can't modify values in p because it is reused by multiple slices.
        case 'sequential'
            dat = dat_all(:, :, :, sl, :, :);
            p_recon = p;
    end

    % -- For debugging: Check the effect of removing the encoding added to one of the simultaneous slices
    if p.debug && p.decode_each_slice
        if p.mica_br || p.mica_rand
            nt_decode = p.nt_to_recon;
        else
            nt_decode = 1;
        end
        dat_decoded_each_slice = mux_decode_each_slice(dat(:, 1:p.inplane_R:end, :, :, :, p.mux*p.num_mux_cycle+(1:nt_decode)), ...
            p.mux, get_ky_omegaz_us_msk(p, size(dat, p.PE_DIM), nt_decode, true));

        % Check synthesized mux data vs. actual mux data
        if p.use_gzblips
            us_msk = get_ky_omegaz_us_msk(p, size(dat, p.PE_DIM), nt_decode, true);
            dat_cal = dat(:, :, 1, 1, :, p.cal_dat_tpoints);
            dat_cal = 1/size(dat_cal, p.T_DIM) * fft(dat_cal, [], p.T_DIM);
            dat_cal_sz = get_dat_sz(dat_cal, p);
            nt_compare = size(dat, p.T_DIM) - p.mux*p.num_mux_cycle;
            ec = 1;
            sl_idx = find(p.slices_to_recon == sl);
            if ~exist('dat_acc', 'var')
                dat_acc = zeros(dat_cal_sz.x, dat_cal_sz.y, dat_cal_sz.ec, length(p.slices_to_recon), dat_cal_sz.c, nt_compare);
            end
            dat_acc(:, :, ec, sl_idx, :, :) = dat(:, :, 1, 1, :, p.mux*p.num_mux_cycle+1:end);
            if ~exist('dat_syn', 'var')
                dat_syn = zeros(dat_cal_sz.x, dat_cal_sz.y, dat_cal_sz.ec, length(p.slices_to_recon), dat_cal_sz.c, nt_compare);
            end
            for time = 1 : nt_compare
                if length(us_msk) == 1
                    us_msk_this_time = us_msk;
                else
                    us_msk_this_time = us_msk(time);
                end
                for coil = 1 : dat_cal_sz.c
                    for x = 1 : dat_cal_sz.x
                        tmp = dat_cal(x, 1:p.inplane_R:end, ec, 1, coil, :);
                        tmp = reshape(tmp, [dat_cal_sz.y/p.inplane_R, p.mux]);
                        tmp = tmp .* ifftshift(us_msk_this_time.ftz_pha, 2);
                        dat_syn(x, 1:p.inplane_R:end, ec, sl_idx, coil, time) = sum(tmp, 2);
                    end
                end
            end
        end
    end

    % -- External calibration
    if ~p.internal_cal
        if ~isnumeric(ext_cal)                                           % Input is pfile name or pfile directory
            if p.debug
                fprintf(' Loading data from external calibration p-file...\n');
                fprintf('  Number of MUX cycles in external calibration is %d.\n', p_cal.num_mux_cycle);
            end

            switch p_cal.mux_excited
                case 1                                                   % MUX=1 external calibration
                    sl_cal = sl : p.num_slices : sl+(p.mux-1)*p.num_slices;
                    sl_in_pfile_cal = p_cal.sl_acq_order(sl_cal);
                case p.mux_excited                                       % MUX=p.mux_excited external calibration
                    sl_cal = sl;
                    sl_in_pfile_cal = p_cal.sl_acq_order(sl);
            end

            [dat_cal, p_cal_recon] = epi_load_tseries(ext_cal, ref_cal, vrgf_cal, refp_cal, sl_cal, sl_in_pfile_cal, sl_ref, p_cal);

            %if (size(dat_cal, p.FE_DIM) ~= size(dat, p.FE_DIM)) || (size(dat_cal, p.PE_DIM) ~= size(dat, p.PE_DIM))
            %    error('Ext Cal: matrix size mismatch.\n');
            %end

            if (p_cal_recon.mux_excited == p_recon.mux_excited) && (p_cal_recon.mux_encoded ~= p_recon.mux_encoded)
                p_recon.cap_fov_shift_cal = p_cal_recon.cap_fov_shift_cal; % p_cal_recon.cap_fov_shift_cal and p_cal_recon.mux_encoded may have been changed in function epi_load_tseries
                p_recon.mux_encoded = p_cal_recon.mux_encoded;
            end

            switch p_cal.mux_excited
                case 1                                                   % MUX=1 external calibration
                    dat_cal_num_mux_cycle = 1;                           % Number of mux cycles in the synthesized multiplexed calibration data
                    if p_cal.inplane_R > 1
                        dat_cal = dat_cal(:, :, :, :, :, p_cal.num_mux_cycle); % Use the last fully sampled time point since all later time points are undersampled.
                    else
                        dat_cal = dat_cal(:, :, :, :, :, p_cal.num_mux_cycle+1:end); % Exclude the first few non steady state time points
                        dat_cal = mean(dat_cal, p.T_DIM);                % Average all time points for higher SNR
                    end
                    dat_cal = ifftshift(dat_cal, p.SL_DIM);              % Make the middle slice the 0th indexed slice. Dim: [Kx, Ky, Echo, SimultaneousSlice Z(=p.mux), Coil, Time(=1)]
                    dat_cal = mux_dftz(dat_cal, p.SL_DIM, p.cap_fov_shift_cal, p.mux, 'encode');
                    dat_cal = permute(dat_cal, [1,2,3,6,5,4]);           % Dim: [Kx, Ky, Echo, Slice(=1), Coil, Time(=p.mux)]
                case p.mux_excited                                       % MUX=p.mux_excited external calibration
                    dat_cal_num_mux_cycle = p_cal.num_mux_cycle;
            end
        else                                                             % Input is calibration data matrix
            dat_cal_num_mux_cycle = size(ext_cal, p.T_DIM) / p.mux;
            dat_cal = ext_cal(:, :, :, sl, :, :);
        end

        dat = dat(:, :, :, :, :, p.mux_encoded*p.num_mux_cycle+1:end);   % Discard internal calibration data, if there were any
        dat = cat(p.T_DIM, dat_cal, dat);
        p_recon.num_mux_cycle = dat_cal_num_mux_cycle;                   % num_mux_cycle must be changed after loading the data from the actual scan and before reconstruction
        p_recon.cal_dat_tpoints = p_recon.mux*(p_recon.num_mux_cycle-1)+1 : p_recon.mux*p_recon.num_mux_cycle; % Time points corresponding to the calibration data
    end

    % -- Coil Compression.
    if p.coil_compress
        if p.debug
            fprintf(' Conducting coil compression...\n');
        end

        [dat, p_recon] = coil_compress(dat, p_recon);
    end

    % -- Recon the data.
    if p.debug
        fprintf(' Reconstructing...\n');
    end
    switch p.mux_recon_method
        case 'sense'
            [im, p_recon] = mux_epi_process_data_sense(dat, p_recon);
        case '1Dgrappa'
            im = mux_epi_process_data_grappa(dat, p_recon);
        case 'slice-grappa'
            im = mux_epi_process_data_slice_grappa(dat, p_recon);
        case 'split-slice-grappa'
            im = mux_epi_process_data_slice_grappa(dat, p_recon);
        otherwise
            error(['Unknown recon method "' p.mux_recon_method '". Aborting.']);
    end

    % -- Add to results the ECC data collected by the mux sequence
    if p.output_ecc_data && ~p.zpad_image
        if size(im, p.C_DIM) == 1
            dat_ecc = sos(dat_ecc, p.C_DIM);
        end
        if p.partial_ky
            dat_ecc_sz = get_dat_sz(dat_ecc, p);
            dat_ecc = cat(p.PE_DIM, dat_ecc, zeros(dat_ecc_sz.x, p.ny_pres-dat_ecc_sz.y, dat_ecc_sz.ec, dat_ecc_sz.sl, dat_ecc_sz.c, dat_ecc_sz.t));
        end
        im = cat(p.T_DIM, dat_ecc, im);
    end

    % -- Common parameters for pseudo-multiple replica and hybrid-space SENSE model based calculations
    if p.calc_snr_pmr || (strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr)
        nt = size(dat, p.T_DIM) - p.mux_excited*p.cap_get_ecc - p.mux_encoded*p.num_mux_cycle; % Number of slice-accelerated time points
    end

    % -- Pseudo-multiple replica: Calculate (scaled) SNR maps and retained SNR maps, then add them to reconed time series so that they will undergo the same dimension adjustment.
    if p.calc_snr_pmr
        if p.debug
            fprintf(' Pseudo-multiple replica simulation...\n');
        end
        if p.coil_compress
            p_recon.num_coils = p.num_coils; % Change number of coils to before coil compression for pseudo-multiple replica simulation
        end
        if p.debug
            tpmr = tic;
        end
        [pmr_rsnr, pmr_noise, pmr_snr, pmr_noise_ref, nmsk] = simulate_pmr_snr(p_recon, nt, p_recon.pmr_raw_data); % pmr_rsnr Dim: [X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=nec), Slice(=nz), Coil(=1), UndersamplingMask(=nmsk)].
        if p.debug
            tpmr = toc(tpmr);
            fprintf(' Pseudo-multiple replica simulation used %f s.\n', tpmr);
        end
        if p.coil_compress
            p_recon.num_coils = p_recon.num_vcoils; % Change number of coils back to after coil compression
        end

        datsz = get_dat_sz(im, p);
        if p.debug
            im = cat(p.T_DIM, zeros(datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, 3*nmsk+1), im);
        else
            im = cat(p.T_DIM, zeros(datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, 2*nmsk), im);
        end
        im(:, :, :, :, 1, 1:nmsk) = pmr_rsnr; % Add retained SNR maps to the first coil dimension, time points 1~nmsk
        im(:, :, :, :, 1, nmsk+1:2*nmsk) = pmr_noise; % Add (scaled) noise maps to the first coil dimension, time points nmsk+1~2*nmsk
        if p.debug
            im(:, :, :, :, 1, 2*nmsk+1:3*nmsk) = pmr_snr; % Add (scaled) SNR maps to the first coil dimension, time points 2*nmsk+1:3*nmsk
            im(:, :, :, :, 1, 3*nmsk+1) = pmr_noise_ref; % Add (scaled) reference noise maps to the first coil dimension, time point 3*nmsk+1
        end
    end

    % Hybrid-space SENSE model based retained SNR maps
    if strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr
        if p.debug
            fprintf(' Calculating hybrid-space SENSE model based retained SNR...\n');
        end
        if p.debug
            t_sense = tic;
        end
        [sense_rsnr, sense_noise, sense_noise_ref, nmsk] = calc_sense_rsnr(p_recon, nt, p_recon.pmr_raw_data); % sense_rsnr Dim: [X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=1), Slice(=nz), Coil(=1), UndersamplingMask(=nmsk)].
        if p.debug
            t_sense = toc(t_sense);
            fprintf(' Hybrid-space SENSE model based retained SNR used %f s.\n', t_sense);
        end
        p_recon = rmfield(p_recon, 'pmr_raw_data'); % pmr_raw_data Won't be used any more

        datsz = get_dat_sz(im, p);
        if p.debug
            im = cat(p.T_DIM, zeros(datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, 2*nmsk+1), im);
        else
            im = cat(p.T_DIM, zeros(datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, 2*nmsk), im);
        end
        im(:, :, :, :, 1, 1:nmsk) = sense_rsnr; % Add retained SNR maps to the first coil dimension, time points 1~nmsk
        im(:, :, :, :, 1, nmsk+1:2*nmsk) = sense_noise; % Add (scaled) noise maps to the first coil dimension, time points nmsk+1~2*nmsk
        if p.debug
            im(:, :, :, :, 1, 2*nmsk+1) = sense_noise_ref; % Add (scaled) reference noise maps to the first coil dimension, time point 2*nmsk+1
        end
    end

    if p.debug
        fprintf(' Preparing outputs...\n');
    end

    % -- Fix phase-encode direction for pepolar scans
    if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
        im = flip(im, p.PE_DIM);
        pepolar = 1;
    else
        pepolar = 0;
    end

    % -- Fix slice ordering.
    if p.descend_acq                                                     % If the multiplexed slices are acquired in descending order. Need this because the unmuxed slices are always in ascending order (i.e., from bottom to top) for the current DFT encoding scheme in the PSD.
        sl_loc = (p.num_slices*(p.mux-1)+sl) : -p.num_slices : 1;
    else
        sl_loc = sl : p.num_slices : p.num_unmuxed_slices;
    end

    % -- Save the recon results.
    % Note that size(im,6) might be less than the prescribed timeseries length
    % (i.e., p.num_mux_cycle+p.nt_to_recon) if the scan was stopped prematurely.
    if p.zpad_image
        nx_image = p.zpad_size(1);
        ny_image = p.zpad_size(2);
    else
        nx_image = p.nx_pres;
        ny_image = p.ny_pres;
    end
    if p.return_sos_im                                                                      % Return SOS coil-combined magnitude image
        d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(im, p.T_DIM)];             % Dim: [X, Y, Slice, Time]
    else                                                                                    % Return complex images
        if strcmp(p.mux_recon_method, 'sense')                                              % SENSE recon, return coil-combined complex image
            d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(im, p.T_DIM)];         % Dim: [X, Y, Slice, Time]
        else                                                                                % GRAPPA recon, return complex single-coil images (when p.sense_coil_comb is false) or coil-combined complex image (when p.sense_coil_comb is true)
            if ~p.sense_coil_comb
                d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(im, p.T_DIM), size(im, p.C_DIM)]; % Dim: [X, Y, Slice, Time, Coil]
            else
                d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(im, p.T_DIM)];     % Dim: [X, Y, Slice, Time]
            end
        end
    end
    if p.isdifscan && (p.next2 > 1)
        im = mean(im, 3);
    end
    if p.return_sos_im
        im = reshape(im, [nx_image, ny_image, p.mux, size(im, p.T_DIM)]);
    else
        if strcmp(p.mux_recon_method, 'sense')
            im = reshape(im, [nx_image, ny_image, p.mux, size(im, p.T_DIM)]);
        else
            if ~p.sense_coil_comb
                im = permute(reshape(im, [nx_image, ny_image, p.mux, size(im, p.C_DIM), size(im, p.T_DIM)]), [1,2,3,5,4]); % Dim: [X, Y, Slice, Time, Coil]
            else
                im = reshape(im, [nx_image, ny_image, p.mux, size(im, p.T_DIM)]);
            end
        end
    end
    if numel(p.slices_to_recon) > 1
        if isempty(d)
            d = zeros(d_size);
        end
        d(:, :, sl_loc, :, :) = im;
    else
        d = im;
    end
end

% -- Adjust orientation of reconed images
% Axial and Axial-like Oblique:         R/L - A/P (1stImageDimension - 2ndImageDimension)
% Coronal and Coronal-like Oblique:     R/L - S/I
% Sagittal and Sagittal-like Oblique:   A/P - S/I
if p.swappf
    % RFD: when phase/freq are swapped, axial & axial-like obliques are left-right flipped.
    % TODO: check that this flipdim is nexessary for other orientations.
    if p.scan_orient(1)=='a'
        % Only axials are flipped
        d = FLIP(d, p.PE_DIM);
    elseif p.scan_orient(1)=='c'
        % for coronals, the freq dim seems to be flipp!?!?!?
        d = FLIP(d, p.FE_DIM);
    end
    d = permute(d, [2,1,3,4,5]); % Dim: [1stImageDimension, 2ndImageDimension, Slice, Time, Coil(=nc_after_coil_compression if d contains complex images from GRAPPA recon; =1 otherwise)]
    d_size([1 2])=d_size([2 1]);
end

% -- SNR maps
snr_d = [];
% Hybrid-space SENSE model based retained SNR and noise: Separate maps from reconstructed images.
if strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr
    snr_d.sense_rsnr = d(:, :, :, 1:nmsk, 1); % retained SNR map
    snr_d.sense_noise = d(:, :, :, nmsk+1:2*nmsk, 1); % (scaled) noise map
    if p.debug
        snr_d.sense_noise_ref = d(:, :, :, 2*nmsk+1, 1); % (scaled) noise map for reference scan(either a MUX1 scan or a slice-DFT-encoded calibration scan) that has the same inplane acceleration as the mux scan
        d = d(:, :, :, 2*nmsk+2:end, :);
    else
        d = d(:, :, :, 2*nmsk+1:end, :);
    end
end
% Pseudo-multiple replica: Separate SNR maps and reconstructed images.
if p.calc_snr_pmr
    snr_d.pmr_rsnr = d(:, :, :, 1:nmsk, 1); % retained SNR map
    snr_d.pmr_noise = d(:, :, :, nmsk+1:2*nmsk, 1); % (scaled) noise map
    if p.debug
        snr_d.pmr_snr = d(:, :, :, 2*nmsk+1:3*nmsk, 1); % (scaled) SNR map
        snr_d.pmr_noise_ref = d(:, :, :, 3*nmsk+1, 1); % (scaled) noise map for reference scan(either a MUX1 scan or a slice-DFT-encoded calibration scan) that has the same inplane acceleration as the mux scan
        d = d(:, :, :, 3*nmsk+2:end, :);
    else
        d = d(:, :, :, 2*nmsk+1:end, :);
    end
end

% -- Scale the images according to the receive gains
%scaling = 10 * 2^((p.r2 - 15) + 0.5 * (p.r1 - 14));
if p.edr_flag
    scaling = 1/2^16;
else
    scaling = 100;
end
d = d * scaling;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- SNR maps
% Scale pseudo-multiple replica noise maps
if p.calc_snr_pmr
    snr_d.pmr_noise = snr_d.pmr_noise * scaling;
    if p.debug
        snr_d.pmr_noise_ref = snr_d.pmr_noise_ref * scaling;
    end
end
% Actual multiple replica
if p.calc_snr_amr
    nt_dat_cal = (p.output_ecc_data && ~p.zpad_image)*1 + p.num_mux_cycle; % Number of time points for the ECC data and calibration data
    if size(d, 4) >= nt_dat_cal + p.amr_min_nt
        [snr_d.amr_snr, ~, snr_d.amr_noise] = imsnr(sos(d(:, :, :, nt_dat_cal + 1:end, :), p.C_DIM), 4); % snr_d.amr_snr: SNR map; snr_d.amr_noise: noise map
    end
end
% Hybrid-space SENSE model based noise maps
if strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr
    snr_d.sense_noise = snr_d.sense_noise * scaling;
    if p.debug
        snr_d.sense_noise_ref = snr_d.sense_noise_ref * scaling;
    end
end

% -- Timing
if debug && ~isempty(outfile)
    time_used = toc(time_used);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Save results to file
if exist('d', 'var') && ~isempty(d) && ~isempty(outfile)
    if exist('octave_config_info', 'builtin')
        save(outfile, 'd', 'sl_loc', 'd_size', '-v7');
    else
        if ~(p.debug && p.decode_each_slice && p.use_gzblips)
            dat_acc = [];
            dat_syn = [];
            us_msk = [];
        end
        if p.debug && p.decode_each_slice
            im_decoded_each_slice = sos(ifft2c(dat_decoded_each_slice), p.C_DIM);
        else
            im_decoded_each_slice = [];
        end
        if ~p.debug
            time_used = 'not recorded';
        end
        save(outfile, 'd', 'sl_loc', 'd_size', 'snr_d', 'p', 'dat_acc', 'dat_syn', 'im_decoded_each_slice', 'us_msk', 'time_used', 'ext_cal', 'pfile', '-v7');
    end
end

return
