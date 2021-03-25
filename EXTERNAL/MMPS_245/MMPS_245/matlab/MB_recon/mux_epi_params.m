function p = mux_epi_params(pfile, slices, nt_to_recon, n_vcoils, debug, recon_method, apply_fermi, use_homodyne, notch_thresh, noise_fname)
%
% p = mux_epi_params(pfile, [slices=[]], [nt_to_recon=[]], [n_vcoils=[]], [debug=false], [recon_method='1Dgrappa'], [apply_fermi=false], [use_homodyne=true], [notch_thresh=0], [noise_fname])
%
% Generates all parameters needed for reconstructing slice-multiplexed EPI data.
%
% Inputs:
%   pfile       - The p-file (need information in p-file header).
%   slices      - Indices of the (muxed) slices to reconstruct. Default: all slices.
%   nt_to_recon - Number of time points to reconstruct, excluding the first
%                 few mux phase cycling time points. Default: all time points.
%   n_vcoils    - Number of virtual coils for coil compression. Set to [] for no coil compression.
%   debug       - True: calculate intermediate results and print out messages.
%   recon_method- '1Dgrappa' or numbers other than 1: use 1D-GRAPPA;
%                 'sense' or number 1: use SENSE;
%                 'slice-grappa': use slice-GRAPPA;
%                 'split-slice-grappa': use split-slice-GRAPPA.
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
%   noise_fname - Filename of the noise.dat (noise data) file.
%
% Output:
%   p           - A structure with all parameters.
%
% (c) Bob Dougherty <bobd@stanford.edu> Stanford University     September 2012
% Modified by Kangrong Zhu              Stanford University     September 2012

% -- Inputs and constants
if ~exist('pfile', 'var') || ~exist(pfile,'file')
    error('Need a p-file!');
end

if(~exist('slices', 'var'));                               slices      = [];    end
if(~exist('nt_to_recon', 'var'));                          nt_to_recon = [];    end
if(~exist('n_vcoils', 'var'));                             n_vcoils    = [];    end
if(~exist('debug', 'var')) || isempty(debug);              debug       = false; end
if(~exist('recon_method', 'var') || isempty(recon_method));recon_method= '1Dgrappa'; end
if(~exist('apply_fermi', 'var'));                          apply_fermi = false; end
if(~exist('use_homodyne', 'var'));                         use_homodyne = true; end
if(~exist('notch_thresh', 'var'));                         notch_thresh = 0;    end

% Dimension of the k-/image-space data
p.FE_DIM = 1;        % Dimension for frequency encoding
p.PE_DIM = 2;        % Dimension for phase encoding
p.EC_DIM = 3;        % Dimension for echoes
p.SL_DIM = 4;        % Dimension for slices
p.C_DIM  = 5;        % Dimension for coils
p.T_DIM  = 6;        % Dimension for time points
p.KZ_DIM = 7;        % Dimension for Kz, after reformulating the slice-multiplexed data into a 3D data set

% For the 'repmat' function
p.KEEP_ORIG_SZ = 1;  % Keep the original size in a dimension when repeating a matrix.

% For ky traversal direction
p.TOP_DOWN   = 0;
p.CENTER_OUT = 1;
p.BOTTOM_UP  = 2;

% Number of coefficients per ky line for EPI eddy current correction (i.e. phase correction in x-ky space)
p.NUM_COE_PER_KY_LINE = 2;

% For virtual coil concept
p.ONLY_ACTUAL_COILS        = 0;
p.ACTUAL_AND_VIRTUAL_COILS = 1;
p.ONLY_VIRTUAL_COILS       = 2;

% For Homodyne partial-k recon
RHTYP1HOMODYNE = 1;  % The bit indicating homodyne recon in rhtype1(int16). GE software manual for rhtype1: bit0 (1).

% -- Load full p-file header.
hdr = read_MR_headers(pfile, 'all', 'raw');

% -- Set parameters using p-file header.
% Save some header fields into a structure to avoid loading them again when loading the k-space data using load_raw_tseries.m
p.pfile_header.version    = hdr.rdb_hdr.rdbm_rev;
p.pfile_header.npasses    = hdr.rdb_hdr.npasses;
p.pfile_header.nslices    = hdr.rdb_hdr.nslices;
p.pfile_header.nframes    = hdr.rdb_hdr.nframes;
p.pfile_header.nechoes    = hdr.rdb_hdr.nechoes;
p.pfile_header.hnover     = hdr.rdb_hdr.hnover;
p.pfile_header.frsize     = hdr.rdb_hdr.frame_size;
p.pfile_header.ncoils     = hdr.rdb_hdr.dab(2)-hdr.rdb_hdr.dab(1)+1;
p.pfile_header.ptsize     = hdr.rdb_hdr.point_size;
p.pfile_header.rawhdrsize = hdr.rdb_hdr.off_data;
p.pfile_header.rawsize    = hdr.rdb_hdr.raw_pass_size;

% Scan params
p.num_echoes      = hdr.rdb_hdr.nechoes;                     % Total number of echoes
p.num_slices      = hdr.image.slquant;                       % Total number of multiplexed slices
p.num_coils       = hdr.rdb_hdr.dab(2)-hdr.rdb_hdr.dab(1)+1; % Total number of receiving coils
p.num_passes      = hdr.rdb_hdr.npasses;                     % Total number of time points in the mux epi scan, including the first few mux phase cycling time points
p.inplane_R       = hdr.rdb_hdr.ileaves;                     % In-plane reduction factor (e.g. ARC 2 corresponds to inplane_R = 2)
p.partial_ky      = (hdr.rdb_hdr.hnover ~= 0);               % True: used partial ky acquisition
p.dacq_ctrl       = hdr.rdb_hdr.dacq_ctrl;                   % Bits contianing information about the DACQ read-out
p.kissoff_views   = hdr.rdb_hdr.kissoff_views;               % Not sure what this is, but seems important for RDS PE-line interleaving
p.vrgf            = hdr.rdb_hdr.vrgf;                        % True: ramp sampling on
p.nx_pres         = hdr.rdb_hdr.vrgfxres;                    % Prescribed size in FE. Pass an empty matrix to use the default value. Default: the same as the size in PE of the raw k-space data
p.kydir           = hdr.rdb_hdr.kydir;                       % Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
p.ny_pres         = hdr.image.dim_Y;                         % Prescribed size in PE
p.fov             = hdr.image.dfov;                          % Image field of view, in mm.
p.slthick         = hdr.image.slthick;                       % Slice thickness, in mm
p.slspacing       = hdr.image.scanspacing;                   % Slice spacing (gap between two adjacent prescribed slices), in mm.
p.start_loc       = hdr.series.start_loc;                    % The z location of the 1st prescribed slices, in mm

% MJM
p.end_loc       = hdr.series.end_loc;                        % The z location of the last prescribed slices, in mm

p.mux_excited     = hdr.rdb_hdr.user6;                       % Number of simultaneously excited slices
p.num_mux_cycle   = hdr.rdb_hdr.user7;                       % Number of complete slice phase cycling at the beginning of the slice-multiplexed scan
p.sldist          = hdr.rdb_hdr.user8;                       % Distance between two adjacent simultaneously excited slices, in mm
p.use_gzblips     = hdr.rdb_hdr.user13;                      % True: used Gz blips for slice encoding; False: no Gz blips, direct aliasing of all simultaneous slices
p.cap_random_flag = hdr.rdb_hdr.user17;                      % Type of Gz blip phase encoding. 0: CAIPI; 1: Bit-reversed MICA; 2: Random MICA.
p.swappf          = ~hdr.image.swappf;                       % Whether to swap the phase and frequency encoding directions in the reconed images. True: swap; False: don't swap.
p.flip_angle      = hdr.image.mr_flip;                       % Flip angle, in degrees
p.slice_shuffling = hdr.rdb_hdr.user30;                      % True: slice shuffling during acquisition; false: normal slice order
p.extra_tr        = hdr.rdb_hdr.user31;                      % Wait time after the last slice in each pass
p.pseq            = hdr.rdb_hdr.user32;                      % Type of pulse sequence
p.r1              = hdr.psc.aps_r1;                          % Receive gain R1
p.r2              = hdr.psc.aps_r2;                          % Receive gain R2

% Fermi filter params
p.fermi_radius   = hdr.rdb_hdr.fermi_radius;
p.fermi_width    = hdr.rdb_hdr.fermi_width;
p.fermi_ecc      = hdr.rdb_hdr.fermi_ecc;
p.fermi_rad_scale = 1.0;                                     % Used to apply a more aggressive fermi filter
p.fermi_width_scale = 1.0;

% Slice acquisition order params
switch p.num_slices
    case hdr.rdb_hdr.nslices/hdr.rdb_hdr.npasses
        p.acq_order = 'interleaved';                         % Interleaved acquisition, hdr.data_acq_tab.pass_number is always 0.
        p.sl_acq_order = hdr.data_acq_tab.slice_in_pass(1 : p.num_slices); % The acquisition order of the prescribed slices
    case hdr.rdb_hdr.nslices/hdr.rdb_hdr.reps                % For sequential acquisition, reps = opfphases is the number of time points.
        p.acq_order = 'sequential';                          % Sequential acquisition, hdr.data_acq_tab.slice_in_pass is always 1.
        p.sl_acq_order = hdr.data_acq_tab.slice_in_pass(1) : p.num_slices; % Hardcoded for now for sequential IR scans
        p.pfile_header.npasses = hdr.rdb_hdr.reps;
        p.num_passes = hdr.rdb_hdr.reps;
    otherwise
        error('Inconsistent slice information in header.');
end

% Contains internal calibration or not
p.internal_cal = hdr.rdb_hdr.user5 && (p.num_mux_cycle > 0); % True: acquired mux phase cycling data as internal calibration; False: need external calibration from a separate calibration scan.

% Orientation params
if hdr.series.start_loc > hdr.series.end_loc                 % Header indicates the 1st slice has a higher location in the slice select direction than the last slice.
    p.descend_acq = true;                                    % True: acquired slices in descending order, i.e. from higher frequency to lower frequency in the slice select direction; False: acquired in ascending order.
else
    p.descend_acq = false;
end

if hdr.series.start_ras=='S' || hdr.series.start_ras=='I'
    p.scan_orient = 'axial';
elseif hdr.series.start_ras=='A' || hdr.series.start_ras=='P'
    p.scan_orient = 'coronal';
else
    p.scan_orient = 'sagittal';
end

% RF attributes
p.multiband_TBW = hdr.rdb_hdr.user9;                         % Time-bandwidth product of the multiband RF pulse
p.multiband_pw  = hdr.rdb_hdr.user10/1000.0;                 % Duration of the multiband RF pulse

% Flag indicating the raw data order. When we use the RDS to save a pfile, the data are ordered differently than a normal p-file.
% (Note: for rdsdata saved before 2013.10.23, rhuser16=1 indicates rds data order.)
p.rds_data_order = ( (hdr.rdb_hdr.recon>=9000) || (hdr.rdb_hdr.user16==1) );

% Is diffusion scan or not
if (isfield(hdr.rdb_hdr, 'numdifdirs'))
    p.isdifscan = (hdr.rdb_hdr.numdifdirs > 1);              % True: The scan was a diffusion scan acquired with the muxarcepi2 sequence; False: The scan was acquired with the muxarcepi sequence.
    p.next2     = hdr.rdb_hdr.difnext2;                      % NEX
else
    p.isdifscan = 0;
end

% Params for calibration time points
p.dftz_type = 'inverse';                                     % Type of DFTz encoding conducted on the simultaneous slices. 'inverse': inverse DFT; 'forward': forward DFT.
p.cap_blip_start_cal = hdr.rdb_hdr.user22;                   % Starting index of the kz blips for the calibration time points. 0~(abs(cap_fov_shift_cal)-1) corrrespond to -kzmax~kzmax.
if hdr.rdb_hdr.user25 == 0                                   % Deal with older version sequence
    if p.isdifscan
        % This is set to the # of T2 images for diffusion scans, so default to the mux factor if this is a dwi scan.
        p.cap_fov_shift_cal = p.mux_excited;
    else
        p.cap_fov_shift_cal = hdr.rdb_hdr.user23;                % Older version sequence used rhuser23 to save cap_fov_shift_cal
    end
else
    p.cap_fov_shift_cal = hdr.rdb_hdr.user25;                % CAIPI FOV shift for the calibration time points. A positive integer.
end
if p.cap_fov_shift_cal == 0
    p.cap_fov_shift_cal = p.mux_excited;
end
if (p.mux_excited == 1) && (p.cap_fov_shift_cal ~= 1)
    p.cap_fov_shift_cal = 1;
end
p.mux_encoded = p.cap_fov_shift_cal;
if p.mux_encoded < p.mux_excited
    error('mux_encoded cannot be smaller than mux_excited.');
end
p.recon_cal_use_mux_extra = 0;                               % True: when mux_encoded>mux_excited, use extra encoded slices in recon calibration; False: only use the desired mux_excited slices in recon calibration and ignore the extra slices even if there were any.
if ~(p.use_gzblips && (p.mux_encoded > p.mux_excited))       % Use Gz blips && Encode more slices than excited in calibration time points
    p.recon_cal_use_mux_extra = 0;
end
if p.recon_cal_use_mux_extra
    p.mux = p.mux_encoded;                                   % p.mux is the mux factor in recon
else
    p.mux = p.mux_excited;
end
p.cap_fov_shift_cal = (-1)^strcmp(p.dftz_type, 'inverse') * p.cap_fov_shift_cal; % The type of DFTz encoding for calibration data will be passed as the sign of p.cap_fov_shift_cal.
if ~p.isdifscan
    p.cal_gzblips = hdr.rdb_hdr.user33;                      % True: use gzblip for calibration; false: use RF phase cycling for calibration
end

% Type of acquisition
p.caipi = (p.use_gzblips && p.cap_random_flag == 0);         % CAIPI
p.mica_br = (p.use_gzblips && p.cap_random_flag == 1);       % Bit-reversed MICA
p.mica_rand = (p.use_gzblips && p.cap_random_flag == 2);     % Random MICA
p.mica_perturbed_caipi = (p.use_gzblips && p.cap_random_flag == 3); % MICA, whose kz encoding is CAIPI kz encoding wih random kz perturbation added
p.mica_poisson = (p.use_gzblips && p.cap_random_flag == 4);  % MICA with Poisson-disc-like undersampling pattern on the ky-kz plane

% Params needed only for CAIPI, MICA with randomly perturbed CAIPI scheme, and MICA with Poisson-disc ky-kz pattern
if p.caipi || p.mica_perturbed_caipi || p.mica_poisson
    p.cap_blip_start = hdr.rdb_hdr.user14;                   % Starting index of the kz blips for the accelerated time points. 0~(abs(cap_fov_shift)-1) correspond to -kzmax~kzmax.
    p.cap_blip_inc   = hdr.rdb_hdr.user15;                   % Increment of the kz blip index for adjacent acquired ky lines
    p.cap_fov_shift  = hdr.rdb_hdr.user21;                   % CAIPI FOV shift for the accelerated time points. A positive integer.
    if p.cap_fov_shift == 0
        p.cap_fov_shift = p.mux_excited;
    end
    p.cap_fov_shift = (-1)^strcmp(p.dftz_type, 'inverse') * p.cap_fov_shift; % The shifting in the PE direction of successive slices is FOVy/p.cap_fov_shift. The type of DFTz encoding for accelerated data will be passed as the sign of p.cap_fov_shift.
end

% Params needed only for MICA
if p.mica_br || p.mica_rand || p.mica_perturbed_caipi || p.mica_poisson
    p.cap_seed_shift = hdr.rdb_hdr.user18;                   % Seed shift for subsequent images in subsequent time points
end
if p.mica_perturbed_caipi
    p.cap_kz_rand_pert_range = hdr.rdb_hdr.user26;
end

% Params for eddy current correction
p.md_ecc = 0;                                                % True: use matrix-decoding eddy current correction; False: use the default correction (either a single-slice or a slice-averaged correction).
p.ref_for_default_ecc = 'ref.dat';                           % Which data to use for the default eddy current correction. 'ref.dat': use ref.dat file; 'ecc data': use ECC data collected by the mux sequence; 'ref pfile': use reference scan pfile.
if p.md_ecc
    p.ref_for_md_ecc = 'ref pfile';                          % Which data to use for the matrix-decoding eddy current correction. 'ecc data': use ECC data collected by the mux sequence; 'ref pfile': use reference scan pfile.
end
p.cap_get_ecc = hdr.rdb_hdr.user19;
if p.cap_get_ecc
    p.cap_get_ecc_z = hdr.rdb_hdr.user20;
    p.pcslice       = hdr.rdb_hdr.pcspacial;                 % Slice averaging when calculating x-ky phase correction coefficients. 0: no averaging across slice; >=1 && <=nslics: use one of the slices' coefficients for all slices; -1: average across slices.
    if p.recon_cal_use_mux_extra
        p.ref_for_default_ecc = 'ref.dat';                   % TODO: Enable use of ecc data
        if p.md_ecc
            p.md_ecc = false;                                % TODO: Enable matrix-decoding ghosting correction
        end
    else
        if p.md_ecc
            p.ref_for_default_ecc = 'ecc data';
            p.ref_for_md_ecc = 'ecc data';
        else
            p.ref_for_default_ecc = 'ref.dat';
        end
    end
end

% Receiver coil noise standard deviation from pfile header
if isfield(hdr, 'psc')
    p.coil_noise_std = hdr.psc.rec_std(1:p.num_coils);       % Receiver coil noise standard deviation. Dim: [1, p.num_coils]
    p.coil_noise_std = 1./ ((1./p.coil_noise_std) ./ sqrt(sum((1./p.coil_noise_std).^2))); % Normalize so that sum_k(1/sigma_k.^2)=1, where sigma_k is the noise standard deviation in the k-th coil. Reference: Matt A. Bernstein et al, MRM 1994, 32:330-334, Eq[1]
end

% flag for Extended Dynamic Range
if (p.pfile_header.ptsize == 2)
    p.edr_flag = false;
elseif (p.pfile_header.ptsize == 4)
    p.edr_flag = true;
end

% Params for Homodyne partial-k recon
if p.partial_ky
    if isempty(use_homodyne)
        rhtype1 = hdr.rdb_hdr.data_collect_type1;
        use_homodyne = ( bitand(uint32(RHTYP1HOMODYNE), rhtype1) == RHTYP1HOMODYNE ); % True: use Homodyne; False: use zero-filling
    end
    if use_homodyne
        p.homodyne_niter = 4;                                % Number of iterations for homodyne partial-k recon
    else
        p.homodyne_niter = 0;
    end
    p.homodyne_ntran = hdr.rdb_hdr.ntran;                    % Transition width for homodyne partial-k recon
end

% -- Specify Y chopping
p.ychop = true;                                              % True: used y-chopping in the acquisition; False: no y-chopping. TODO: Get ychop from the pfile header (Whether to use ychop was not set by bit 0 of rhtype)

% -- Specify whether to remove the encoding added to each individual slice for debugging
p.decode_each_slice = false;
p.decode_each_slice = (p.decode_each_slice && debug);        % Set to false when not in debug mode to avoid unnecessary computing.

% -- Specify params for coil compression
if isempty(n_vcoils) || (n_vcoils == 0) || (n_vcoils >= p.num_coils)
    p.coil_compress = false;                                 % True: conduct coil compression before parallel imaging recon
else
    p.coil_compress = true;
    p.cc_method = 'geometric';                               % 'single' or 'geometric'(default). 'single': Use a single compression matrix for all data. 'geometric': Use different compression matrices for different x positions (i.e. geometric decomposition coil compression with alignment).
    p.cc_dat = 'muxcal';                                     % 'muxcal'(default) or 'smap'. 'muxcal': Use calibration k-space data to calculate coil compression matrices; 'smap': Use sensitivity maps to calculate coil compression matrices.
    if strcmp(p.cc_method, 'geometric')
        p.gcc_slwin     = 5;                                 % Odd number of window size in space for computing the coil compression matrices in Geometric decomposition Coil Compression (GCC)
    end
    if n_vcoils >= p.mux * p.inplane_R
        p.num_vcoils = n_vcoils;                             % Number of virtual coils for coil compression
    else
        p.num_vcoils = p.mux * p.inplane_R;                  % To conduct parallel imaging, must have num_vcoils >= Total_acceleration_factor.
    end
end

% -- Specify params for reconstruction
p.frames             = [];                                   % Always pass an empty matrix to reconstruct all frames.
p.coils              = [];                                   % Always pass an empty matrix to reconstruct all coils.
p.echoes             = [];                                   % Always pass an empty matrix to reconstruct all echoes.
p.num_unmuxed_slices = p.num_slices * p.mux;                 % Total # of unmuxed slices.
p.cal_dat_tpoints    = p.mux*(p.num_mux_cycle-1)+1 : p.mux*p.num_mux_cycle; % Time points corresponding to the calibration data, i.e. the last group of mux phase cycling
p.zpad_image         = 0;                                    % True: Zero pad reconstructed images to size p.zpad_size
if p.zpad_image
    p.zpad_size      = [128, 128];                           % Matrix size to pad the reconstructed images to
end

% -- Specify reconstruction algorithm to use
if ~ischar(recon_method)
    % Support the old calling convention: 0=1Dgrappa, 1=sense
    if recon_method==1
        recon_method = 'sense';
    else
        recon_method = '1Dgrappa';
    end
end

% Check for sense1 flag
sense1 = 0;
tmp = strsplit(recon_method,'_');
if numel(tmp)>1
    if strcmp(tmp{2}, 'sense1')
        sense1 = 1;
    else
        error(['Unknown recon flag "' tmp{2} '". Aborting.']);
    end
end
recon_method = tmp{1};

if sum(strcmp(recon_method, {'1Dgrappa', 'sense', 'slice-grappa', 'split-slice-grappa'})) > 0
    p.mux_recon_method = recon_method; % 'sense', '1Dgrappa', 'slice-grappa' or 'split-slice-grappa'
else
    error(['Unknown recon method "' recon_method '". Aborting.']);
end
if (p.mica_br || p.mica_rand || p.mica_perturbed_caipi || p.mica_poisson) && strcmp(p.mux_recon_method, '1Dgrappa')
    if debug
        fprintf('1D-GRAPPA doesn''t work with MICA. Set reconstruction method to Slice-GRAPPA.\n');
    end
    p.mux_recon_method = 'slice-grappa';
end

% -- Noise related parameters
% Specify SNR map params
p.calc_snr_pmr = 0;                                          % True: calculate SNR maps using pseudo-multiple replica method.
if p.mux_excited == 1
    p.calc_snr_pmr = 0;                                      % Donot need to calculate pseudo-multiple replica SNR or retained SNR maps for single-slice acquisition.
end
if p.md_ecc
    p.calc_snr_pmr = 0;
end
if strcmp(p.mux_recon_method, 'sense')
    p.calc_sense_rsnr = 1;                                   % True: Calculate retained SNR and noise maps basing on the hybrid-space SENSE reconstruction model
    if p.mux_excited == 1
        p.calc_sense_rsnr = 0;
    end
    if p.partial_ky && p.homodyne_niter > 0
        p.calc_sense_rsnr = 0;                               % Model-based rSNR calculation only works with zero-filling (doesn't work with Homodyne) for partial ky acquisition
    end
    if p.md_ecc
        p.calc_sense_rsnr = 0;
    end
else
    p.calc_sense_rsnr = 0;
end

% Receiver coil noise covariance matrix from noise.dat file
if (p.calc_snr_pmr || (strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr)) && exist('noise_fname', 'var') && ~isempty(noise_fname) && exist(noise_fname, 'file')
    p.psi_mtx = noise_file_to_cov(noise_fname, p.num_coils);

    % Empirical scaling
    % TODO: Fix this by understanding the scaling in the data acquisition.
    if p.isdifscan % muxarcepi2 sequence
        p.psi_mtx = p.psi_mtx .* 14;
    else % muxarcepi sequence
        p.psi_mtx = p.psi_mtx ./ 10^8;
    end
end

% Specify what type of data to use for data whitening
p.whiten_type = 'coil_noise_std'; % 'psi_mtx': Use (p.psi_mtx)^(-1/2); 'coil_noise_std'(This option is needed because many scans don't measure psi_mtx): Use diag(1./p.coil_noise_std); 'none': don't conduct data whitening
if strcmp(p.whiten_type, 'coil_noise_std')
    if ~isfield(p, 'coil_noise_std') || isempty(p.coil_noise_std) || any(isnan(p.coil_noise_std))
        p.whiten_type = 'none';
    end
end
if strcmp(p.whiten_type, 'psi_mtx');
    if ~isfield(p, 'psi_mtx') || isempty(p.psi_mtx)
        p.whiten_type = 'none';
    end
end

% -- Specify what type of images to return
if p.isdifscan
    p.return_sos_im = 0;                                     % True: return the square-root-of-sum-of-squares (SOS) images; False: if using SENSE, return the coil-combined complex images, if using GRAPPA, the type of images to return is determined by p.sense_coil_comb.
    % Also, apply a more aggressive fermi filter for diffusion scans
    p.fermi_rad_scale = 0.9;
    p.fermi_width_scale = 1.5;
else
    if sense1
        p.return_sos_im = 0;
    else
        p.return_sos_im = 1;
    end
end
p.sense_coil_comb = 0;
if ~p.return_sos_im
    if sum(strcmp(p.mux_recon_method, {'1Dgrappa', 'slice-grappa', 'split-slice-grappa'})) > 0 % If using GRAPPA
        p.sense_coil_comb = 1;                               % True: return the coil-combined (combined using SENSE) complex image; False: return the complex single-coil images.
    end
end

if p.sense_coil_comb
    if debug
        fprintf('Doing SENSE1 coil combination.\n');
    end
end

% -- Specify params for SENSE reconstruction
if strcmp(p.mux_recon_method, 'sense') || (~p.return_sos_im && p.sense_coil_comb)
    p.crop_smap          = 0;                                % True: crop the sensitivity maps to make SENSE recon work better.
    p.sense_lambda       = [];                               % The lambda for Tikhonov regularization in the SENSE recon. Pass an empty matrix to use the default of norm(encode_mtx' * encode_mtx, 'fro') / size(encode_mtx,2) * p.sense_lambda_default_ratio.
    p.sense_lambda_default_ratio = 0.02;                     % Default regularization coefficient for the SENSE recon is norm(encode_mtx' * encode_mtx, 'fro') / size(encode_mtx,2) * p.sense_lambda_default_ratio.

    % Algorithm to calculate sensitivity maps
    p.smap_type = 'espirit';                                 % Algorithm to calculate the sensitivity maps. 'espirit' or 'coil_over_sos'.

    % Params for ESPIRiT sensitivity map calculations
    if strcmp(p.smap_type, 'espirit')
        p.espirit_nmap = 1;                                  % Number of sensitivity map sets to use in the recon
        p.espirit_catim = 0;                                 % True: concatenate simultaneous slices into one image when calculating sensitivity maps
        if p.mux_excited == 1
            p.espirit_catim = 0;
        end
        if (~p.return_sos_im && p.sense_coil_comb)           % GRAPPA-type recon with SENSE coil combination
            p.espirit_catim = 0;
        end
        if p.espirit_catim
            p.smap_acssz = struct('x', {32}, 'y', {p.ny_pres * (p.mux-1)}); % Enlarge the ACS area
            p.espirit_ksize = [7, 4*p.mux];                 % Enlarge the kernel size
        else
            p.smap_acssz = struct('x', {64}, 'y', {64});    % K-space size to use to calculate the sensitivity maps. Maximum size possible will be used if empty matrices are passed or if the specified size is larger than the maximum size available.
            p.espirit_ksize = [6, 6];                       % Kernel size to calculate the ESPIRiT sensitivity maps. Dim: [Kx, Ky] (The sensitivity maps will be calculated for each individual slice)
        end

        p.espirit_use_c_code = 0;                            % True: Use Martin Uecker's compiled C code for ESPIRiT
        if ispc || p.partial_ky || (~isempty(p.smap_acssz.x) && ~isempty(p.smap_acssz.y) && (p.smap_acssz.x ~= p.smap_acssz.y))
            p.espirit_use_c_code = 0;                        % No compiled C code for Windows. Compiled C code doesn't work if acquisition uses partial ky or if calibration area is not square (if p.smap_acssz.x or p.smap_acssz.y is empty and maximum possible calibration size for x and y doesn't equal, c code won't be used either.).
        end
        if p.espirit_use_c_code
            p.espirit_c_code_dir = '~/Dropbox/ESPIRiT_matlab_and_c_code/recon_v0.1.10';% Directory to the compiled C code for ESPIRiT
            if ~exist(p.espirit_c_code_dir, 'dir')
                p.espirit_use_c_code = 0;
            end
        end
        if ~p.espirit_use_c_code                             % Use Michael Lustig's matlab code for ESPIRiT
            p.espirit_eigThresh_1 = round(1.6 * prod(p.espirit_ksize)); % Threshold to define the nullspace from the 1st SVD in the k-space
            setenv('TOOLBOX_PATH', '');                      % Set environmental variable used by compiled C code to be empty
            p.espirit_tmp_dir = '';                          % Set variable used by compiled C code to be empty
        else                                                 % Use Martin Uecker's compiled C code for ESPIRiT
            p.espirit_eigThresh_1 = 0.001;
            addpath([p.espirit_c_code_dir '/matlab']);       % Add directory to the matlab interface functions
            setenv('TOOLBOX_PATH', p.espirit_c_code_dir);    % Set environmental variable to work around a Matlab bug
            setenv('PATH', strcat(getenv('TOOLBOX_PATH'), ':', getenv('PATH'))); % Set environmental variable to work around a Matlab bug
            if ismac                                         % On Mac
                setenv('DYLD_LIBRARY_PATH', '');
            else                                             % On unix
                setenv('LD_LIBRARY_PATH', '');
            end
            p.espirit_tmp_dir = 'tmp_data/';                 % Temporary directory to store ESPIRiT intermediate results
        end

        if (p.espirit_nmap > 1) || p.espirit_use_c_code
            p.crop_smap = 1;
        end
        if p.crop_smap
            p.espirit_eigThresh_2 = 0.95;                    % Threshold on the eigenvalues for cropping the sensitivity maps
        else
            p.espirit_eigThresh_2 = 0;
        end
    end

    % Params for coil-over-sos sensitivity map calculations
    if strcmp(p.smap_type, 'coil_over_sos')
        p.smap_acssz = struct('x', {[]}, 'y', {[]});         % Always use all data available, better performance than using only the central part of k-space
    end
end

% -- Specify params for GRAPPA reconstruction
% The kernel size represents the number of acquired points used in each
% dimension to interpolate one missing point. The ACS size represents the
% number of points used in each dimension for the ACS area.
if strcmp(p.mux_recon_method, '1Dgrappa') || strcmp(p.mux_recon_method, 'slice-grappa') || strcmp(p.mux_recon_method, 'split-slice-grappa')
    p.grappa_domain = 'imSpace';                             % 'imSpace' or 'kSpace', specifies the domain in which the GRAPPA recon will be carried out.

    p.inplane_kersz = struct( 'x',{7},  'y',{4}  );          % Kernel size for solving inplane acceleration. Pass an empty matrix for a field to use the default value. Default: [7, 4] for [x, y].
    p.inplane_acssz = struct( 'x',{[]}, 'y',{[]} );          % ACS size for solving inplane acceleration. Pass an empty matrix for a field to use the default value. Default: use all available calibration k-space data.

    if strcmp(p.mux_recon_method, '1Dgrappa')
        p.mux_kersz_1d_grappa = struct( 'x',{7},  'y',{4}  ); % Kernel size for solving slice multiplexing. Pass an empty matrix for a field to use the default value. Default: [7, 4] for [x, y].
        p.mux_acssz_1d_grappa = struct( 'x',{32}, 'y',{p.ny_pres * (p.mux-1)} ); % ACS size for solving slice multiplexing. Pass an empty matrix for a field to use the default value. Default: [32, p.ny_pres * (p.mux-1)] for [x, y].

        p.zpad_1Dgrappa = 'minZpad';                         % How to zero pad in 1D GRAPPA recon. 'minZpad': pad minimum size zero matrices between adjacent slices; 'normalZpad': pad zero matrices of the same size as the original image matrix.
        if p.caipi && (abs(p.cap_fov_shift) ~= p.mux) && strcmp(p.zpad_1Dgrappa, 'minZpad');
            p.zpad_1Dgrappa = 'normalZpad';                  % Minimum zero padding does not work for CAIPI if abs(p.cap_fov_shift) ~= p.mux
        end
    else
        p.mux_kersz_slice_grappa = struct( 'x',{7}, 'y',{7} ); % Kernel size for (split-)slice-GRAPPA. Pass an empty matrix for a field to use the default value. Default: [7, 7] for [x, y]
        p.mux_acssz_slice_grappa = struct( 'x',{[]}, 'y',{[]} ); % ACS size for (split-)slice-GRAPPA. Pass an empty matrix for a field to use the default value. Default: Use all available k-space data.
        p.slice_grappa_odd_even_fit = 0; % True: Fit and apply different kernels to the odd and even ky lines in slice-GRAPPA.
        if p.slice_grappa_odd_even_fit
            p.ref_for_md_ecc = 'ref pfile'; % Will use the md_ecc function to calculate residual ghost
            if p.cap_get_ecc
                p.ref_for_md_ecc = 'ecc data';
            end
        end
        switch p.mux_recon_method
            case 'slice-grappa'
                p.use_split_slice_grappa = false;            % True: Use split-slice-GRAPPA; False: Use standard slice-GRAPPA.
            case 'split-slice-grappa'
                p.use_split_slice_grappa = true;
        end
    end
end

% -- More for Nyquist ghost correction
if p.md_ecc && ~strcmp(p.mux_recon_method, 'sense')
    if ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit
        p.md_ecc = false;
        if debug
            fprintf('Set p.md_ecc to false since recon is (split)-slice-grappa with odd even kernel fitting.\n');
        end
    else
        if debug
            fprintf('Reconstruction is 1D-GRAPPA or (split)-slice-grappa without odd even kernel fitting, will only apply matrix decoding eddy current correction to calibration data.\n');
        end
    end
end

if ~strcmp(p.ref_for_default_ecc, 'ref.dat') || p.md_ecc || (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit)
    if p.use_gzblips && p.cap_get_ecc && p.cap_get_ecc_z     % Use Gz blips && Measured eddy current effects with Gz blips on
        p.smooth_pha_coe_acc = 0;                            % True: smooth the x-ky phase correction coefficients along the ky direction for the accelerated data; False: no smoothing along ky.
    else
        p.smooth_pha_coe_acc = 1;
    end
end

p.output_ecc_data = (p.cap_get_ecc && (strcmp(p.ref_for_default_ecc, 'ecc data') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ecc data')) || (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit && strcmp(p.ref_for_md_ecc, 'ecc data')))); % True: output ECC data collected by the mux sequence along with the reconstructed images; False: don't output the ECC data.

% -- Specify whether to add virtual coil concept in recon
% The virtual concept(Reference: Martin Blaimer et al. MRM 2009; 61:93-102)
% is different from coil compression: It synthesizes an additional set of
% data to better utilize the encoding power of the phase(including background phase) in the sensitivity profiles.
% Encoding scheme for virtual coils: Take conjugate of measured data <=>
% Take conjugate of sensitivity, use negative frequencies, and take conjugate of eddy current effect phase term.
p.add_vcc = p.ONLY_ACTUAL_COILS;                             % p.ONLY_ACTUAL_COILS: Use data in only actual coils; p.ACTUAL_AND_VIRTUAL_COILS: Use data in both actual and virtual coils; p.ONLY_VIRTUAL_COILS: Use data in only virtual coils.
if p.add_vcc
    if ~p.internal_cal
        if debug
            fprintf('Not using virtual coils because they currently don''t work with external calibration.');
        end
        p.add_vcc = p.ONLY_ACTUAL_COILS;                     % TODO: external calibration may work with VCC if phase mismatch is corrected either in PSD or in recon
    end
    switch p.mux_recon_method
        case 'sense'
            switch p.smap_type
                case 'coil_over_sos'
                    p.smap_acssz = struct('x', {[]}, 'y', {[]}); % Must use all data available, must capture rapid phase change for virtual coil concept to work well
                case 'espirit'
                    if debug
                        fprintf('Not using virtual coils because they currently don''t work with ESPIRiT sensitivity maps.');
                    end
                    p.add_vcc = p.ONLY_ACTUAL_COILS;         % TODO: add virtual coil concept for ESPIRiT maps
            end
        case '1Dgrappa'
            if debug
                fprintf('Not using virtual coils because they currently don''t work with GRAPPA.');
            end
            p.add_vcc = p.ONLY_ACTUAL_COILS;                 % TODO: add virtual coil concept for GRAPPA
    end
end

% -- Reconstruction params input from the user
p.debug = debug;                                             % True: calculate intermediate results and print out messages.
p.apply_fermi = apply_fermi;                                 % True: apply a circular Fermi filter using the rec.fermi parameters specified in the p-file header.
p.notch_thresh = notch_thresh;                               % If non-zero, a denotching filter will be applied. Any point in k-space with a value lower than p.notch_thresh will be replaced by an adjacent time point.

if isempty(slices)
    slices = 1 : p.num_slices;
end
[~, nslices_inrange] = checkrange(slices, 1, p.num_slices);
if nslices_inrange < length(slices)
    error('Acquired %d slices. Slices to recon contain out of range slice indices.', p.num_slices);
end
p.slices_to_recon = slices;                                  % Indices of the (muxed) slices to reconstruct

p = mux_epi_params_set_tpoints(p, nt_to_recon);

% -- For SNR calculation
% Actual multiple replica method
p.calc_snr_amr = 0;                                          % True: calculate SNR maps using actual multiple replica method.
if (p.mica_br || p.mica_rand || p.mica_perturbed_caipi || p.mica_poisson) && p.cap_seed_shift ~= 0 % If the MICA undersampling pattern is changing with time
    p.calc_snr_amr = 0;
end
if p.isdifscan                                               % If the scan is a diffusion scan where the diffusion encoding direction changes with time
    p.calc_snr_amr = 0;
end
if p.calc_snr_amr
    p.amr_min_nt = 10;                                       % Minimum number of time points needed for calculate SNR maps of the reconstructed images.
end

% Pseudo-multiple replica method
if ~isfield(p, 'psi_mtx') || isempty(p.psi_mtx)
    p.calc_snr_pmr = 0;                                      % Cannot conduct pseudo-multiple replica simulation if no coil noise covariance matrix exists.
end
if p.calc_snr_pmr
    p.pmr_num_rep = 50;                                      % Number of replicas to simulate using the pseudo-multiple replica method.
    if ~strcmp(p.mux_recon_method, 'sense') && (p.inplane_R == 1)
        p.pmr_ref_type = 'muxcal';                           % Type of reference image in the pseudo-multiple replica simulation, 'muxcal'(faster) or 'mux1'(slower). 'muxcal': Use fully DFT encoded slice-multiplexed calibration acquisition; 'mux1': Use single-slice acquisition. For p.inplane_R=1, The simulated retained SNR maps are the same using either type.
    else
        p.pmr_ref_type = 'mux1';                             % When p.inplane_R > 1, the simulated retained SNR are different using 'muxcal' or 'mux1' as reference, because 'muxcal' data have no inplane acceleration and 'mux1' data have inplane acceleration. Since we want the retained SNR w.r.t. a reference scan which has the same inplane acceleration as the actual mux scan, we set the reference type to 'mux1' here.
    end
end

% Hybrid-Space-SENSE model based retained SNR
if strcmp(p.mux_recon_method, 'sense') && (~isfield(p, 'psi_mtx') || isempty(p.psi_mtx))
    p.calc_sense_rsnr = 0;                                   % Cannot calculate hybrid-space-SENSE model based retained SNR if no coil noise covariance matrix exists.
end
if strcmp(p.mux_recon_method, 'sense') && p.calc_sense_rsnr
    p.whiten_type = 'psi_mtx';                               % Set to 'psi_mtx' since 'coil_noise_std' or 'none' would be too slow.
    p.sense_rsnr_im_noise_type = 'actual';                   % 'actual': Calculate noise propagation in the actual reconstruction which uses Tikhonov regularization. 'snr_optimal'(much slower): Calculate noise propagation in an SNR optimal reconstruction which assumes the sample noise covariance is known, as in the SENSE paper.
    p.sense_rsnr_propagate_type = 'diagblk';                 % How to calculate the noise propagation. 'cov': Keep track of the full sample noise covariance matrix. 'diagblk'(much faster): After the 1D iDFTx step, keep track of only the diagonal blocks of the sample noise covariance matrix.
    if strcmp(p.sense_rsnr_im_noise_type, 'snr_optimal')
        p.sense_rsnr_propagate_type = 'cov';                 % propagate_type must be 'cov' if im_noise_type is 'snr_optimal'.
    end
    p.sense_rsnr_ref_type = 'mux1';                          % Type of reference image, 'muxcal'(slower) or 'mux1'(faster). 'muxcal': Use fully DFT encoded slice-multiplexed calibration acquisition with matching inplane acceleration; 'mux1': Use single-slice acquisition with matching inplane acceleration.
end

return
