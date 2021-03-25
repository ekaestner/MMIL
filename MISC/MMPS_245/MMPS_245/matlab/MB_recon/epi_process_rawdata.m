function [dat, p, dat_ecc] = epi_process_rawdata(dat, p, dat_ecc, ref, vrgf, refp, slices, slices_ref)
%
% function [dat, p, dat_ecc] = epi_process_rawdata(dat, p, [dat_ecc=[]], [ref=[]], [vrgf=[]], [refp=[]], [slices=[]], [slices_ref=[]])
%
% Process EPI raw k-space data.
% Steps include: Data whitening, Y chopping correction, Nyquist ghost correction, ramp sampling correction, Fermi filtering and notch filtering.
% Write this into a function so that it can be used in both data reconstruction and pseudo-multiple replica simulation.
%
% Inputs
%   dat        - EPI raw k-space data. Dim: [Kx(=nx_in_raw_data_matrix), Ky(=ny_in_raw_data_matrix), Echo, Slice, Coil, Time].
%   p          - Parameter structure. See mux_epi_params and epi_load_tseries.m for details.
%   dat_ecc    - Reference data for Nyquist ghost correction. This is nonempty only when p.cap_get_ecc == true.
%   ref        - Path to 'ref.dat'.
%   vrgf       - Path to 'vrgf.dat'.
%   refp       - Path to reference scan pfile.
%   slices     - Desired slices(Normally ordered indices) to load. Default: all slices.
%   slices_ref - Slice indices(Normally ordered indices, not pfile indices)
%                in the reference scan corresponding to 'slices' in the actual scan.
%                Only needed if the eddy current correction algorithm will use reference scan pfile.
%
% Outputs
%   dat        - Processed k-space data. Dim: [Kx(=p.nx_pres), Ky(=ny_in_raw_data_matrix), Echo, Slice, Coil, Time].
%   p          - Same as the input 'p', with the following field added if p.md_ecc == true:
%                pha_coe(EPI x-ky phase correction coefficients, will be passed to
%                mux_epi_process_data_sense.m. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo,
%                Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])),
%                Coil(=1, if p.coil_compress == true; =nc, if p.coil_compress == false)].
%   dat_ecc    - The input dat_ecc which has gone through ramp sampling correction. This is nonempty only when p.cap_get_ecc == true.
%
% (c) Kangrong Zhu      Stanford University     Dec 2014

if ~exist('dat_ecc', 'var');    dat_ecc = [];    end
if ~exist('ref', 'var');        ref = [];        end
if ~exist('vrgf', 'var');       vrgf = [];       end
if ~exist('refp', 'var');       refp = [];       end
if ~exist('slices', 'var');     slices = [];     end
if ~exist('slices_ref', 'var'); slices_ref = []; end

% -- For pseudo-multiple replica simulation: Save to or load from the paramter structure
if p.calc_snr_pmr || p.calc_sense_rsnr
    if p.cap_get_ecc
        [p, dat_ecc] = pmr_save_load_params(p, dat_ecc);
    end
    [p, ref] = pmr_save_load_params(p, ref);
    if p.vrgf
        [p, vrgf] = pmr_save_load_params(p, vrgf);
    end
    if (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || strcmp(p.ref_for_default_ecc, 'ref pfile')
        [p, refp] = pmr_save_load_params(p, refp);
    end
    [p, slices] = pmr_save_load_params(p, slices);
    if p.md_ecc || (~p.md_ecc && strcmp(p.ref_for_default_ecc, 'ref pfile'))
        [p, slices_ref] = pmr_save_load_params(p, slices_ref);
    end
end

% Check params
if p.cap_get_ecc && isempty(dat_ecc)
    error('Header indicates ECC data was acquired with the sequence, but the inputs don''t contain it.');
end
if isempty(ref)
    error('No ref file');
end
if p.vrgf && isempty(vrgf)
    error('No vrgf file.');
end
if ((p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || strcmp(p.ref_for_default_ecc, 'ref pfile')) && isempty(refp)
    error('No reference scan pfile.');
end
if isempty(slices)
    error('No slices specified.');
end
if ( (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || (~p.md_ecc && strcmp(p.ref_for_default_ecc, 'ref pfile')) ) && isempty(slices_ref)
    error('No slices for reference scan specified.');
end

datsz = get_dat_sz(dat, p);

if p.vrgf && (datsz.x == p.nx_pres)
    error('Header indicates ramp sampling on, but data are not.');
end

% -- Data whitening
% After the whitening transformation, the coil noise covariance matrix will become identity and can be ignored in further data processing
% (SOS will directly be SNR optimal coil combination and the sample covariance matrix in SENSE recon will become identity)
if ~strcmp(p.whiten_type, 'none')
    switch p.whiten_type
        case 'coil_noise_std' % Approximation using coil noise standard deviation
            whiten_mtx = get_data_whiten_mtx(p.coil_noise_std, p.whiten_type);
        case 'psi_mtx' % Use full coil noise covariance matrix
            whiten_mtx = get_data_whiten_mtx(p.psi_mtx, p.whiten_type);
    end
    
    dat = reshape(permute(dat, [5,1,2,3,4,6]), [datsz.c, datsz.x*datsz.y*datsz.ec*datsz.sl*datsz.t]); % Dim: [Coil, Kx-Ky-Echo-Slice-Time]
    dat = whiten_mtx * dat; % If covariance matrix M of the random vector X is positive definite, then M^(-1/2)X has covariance matrix I. Reference: Matt A. Bernstein et al, MRM 1994, 32:330-334, Eq[1]
    dat = permute(reshape(dat, [datsz.c, datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.t]), [2,3,4,5,1,6]); % Dim: [Kx, Ky, Echo, Slice, Coil, Time]
end

% -- Y chopping correction
if p.ychop
    dat(:, 2:2:end, :, :, :, :) = -dat(:, 2:2:end, :, :, :, :);
    if p.cap_get_ecc
        dat_ecc(:, 2:2:end, :, :, :, :) = -dat_ecc(:, 2:2:end, :, :, :, :);
    end
end

% -- Nyquist ghost correction
if p.md_ecc || (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit) % Matrix-decoding ghost correction or (split-)slice-GRAPPA with odd even kernel fitting
    if strcmp(p.ref_for_md_ecc, 'ref pfile') && ~exist(refp, 'file')
        error('No reference scan p-file.');
    end
    if strcmp(p.ref_for_md_ecc, 'ecc data') && ~p.cap_get_ecc
        error('No ECC data collected by the mux sequence.');
    end
    if p.debug
        if p.md_ecc
            fprintf('  Matrix-decoding ghost correction...\n');
        else
            fprintf('  (Split)-slice-GRAPPA with odd even kernel fitting ghost correction...\n');
        end
    end
    if strcmp(p.ref_for_md_ecc, 'ref pfile')
        md_ref = refp;
    else
        md_ref = dat_ecc;
    end
    if p.vrgf
        % Apply default (single-slice or slice-averaged) correction as a first-pass correction
        switch p.ref_for_default_ecc
            case 'ecc data'
                [dat, pha_coe_default] = default_ecc(dat, p, dat_ecc);  % Use reference k-space data collected by the mux sequence in default correction
            case 'ref pfile'
                [dat, pha_coe_default] = default_ecc(dat, p, refp, slices+(floor((p.mux-1)/2))*p.num_slices); % Use reference scan pfile in default correction
            case 'ref.dat'
                [dat, pha_coe_default] = default_ecc(dat, p, ref, slices); % Use ref.dat file in default correction
        end
    else
        if p.md_ecc
            [dat, p.pha_coe] = md_ecc(dat, p, md_ref, [], [], slices_ref);
        else % (split-)slice-GRAPPA with odd even kernel fitting
            [~, p.pha_coe] = md_ecc([], p, md_ref, [], [], slices_ref);
        end
    end
else                                                                    % Default(single-slice or slice-averaged) ghost correction
    if p.debug
        fprintf('  Default ghost correction...\n');
    end
    switch p.ref_for_default_ecc
        case 'ecc data'
            [dat, pha_coe_default] = default_ecc(dat, p, dat_ecc); % Use reference k-space data collected by the mux sequence in default correction
        case 'ref pfile'
            [dat, pha_coe_default] = default_ecc(dat, p, refp, slices_ref); % Use reference scan pfile in default correction
        case 'ref.dat'
            [dat, pha_coe_default] = default_ecc(dat, p, ref, slices); % Use ref.dat file in default correction
    end
    if p.calc_sense_rsnr
        p.pha_coe_default = pha_coe_default;
    end
end

% -- Ramp sampling correction
if p.vrgf
    
    % Ramp sampling correction
    if isfield(p, 'pmr_ramp_flt') && ~isempty(p.pmr_ramp_flt)
        ramp_flt = p.pmr_ramp_flt;
    else
        ramp_flt = rawload_vrgf(datsz.x, p.nx_pres, vrgf);              % Load the filter contained in the 'vrgf.dat' file.
    end
    
    if size(ramp_flt, 1) ~= p.nx_pres;
        error('Wrong size of ramp filter.');
    end
    dat = epi_vrgf_correct(dat, ramp_flt);
    datsz.x = p.nx_pres;
    
    if p.cap_get_ecc
        dat_ecc = epi_vrgf_correct(dat_ecc, ramp_flt);
    end
    
    % Nyquist ghost correction: matrix decoding correction for ramp on data
    if p.md_ecc
        [dat, p.pha_coe] = md_ecc(dat, p, md_ref, pha_coe_default, ramp_flt, slices_ref);
    end
    
    % Nyquist ghost correction: (split-)slice-GRAPPA with odd even kernel fitting
    if ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit
        [~, p.pha_coe] = md_ecc([], p, md_ref, pha_coe_default, ramp_flt, slices_ref);
    end

    % For pseudo-multiple replica simulation: Save ramp_flt in the paramter structure
    if p.calc_snr_pmr && ~isfield(p, 'pmr_ramp_flt')
        p.pmr_ramp_flt = ramp_flt;
    end
else
    % For pseudo-multiple replica simulation: Save empty ramp_flt in the paramter structure
    if p.calc_snr_pmr && ~isfield(p, 'pmr_ramp_flt')
        p.pmr_ramp_flt = [];
    end
end

% -- Fermi filtering
if p.apply_fermi
    fermi_filt = gen_fermi_filter([datsz.x, p.ny_pres], p.fermi_radius*p.fermi_rad_scale, p.fermi_width*p.fermi_width_scale, 'circ');
    fermi_filt = fermi_filt(:,1:datsz.y);
    nim = datsz.ec * datsz.sl * datsz.c * datsz.t;
    for im_idx = 1 : nim
        dat(:, :, im_idx) = dat(:, :, im_idx) .* fermi_filt;
    end
end

% -- Notch filtering
if p.notch_thresh ~= 0
    nskip = max(p.cal_dat_tpoints);
    if ~p.internal_cal
        % if internal cal was not used, skip the first couple of time points
        nskip = nskip + 2;
    end
    dat = notch_filter(dat, nskip, p.notch_thresh);
end

return

function [dat, pha_coe_acc] = md_ecc(dat, p, ref, pha_coe_1stpass, ramp_flt, slices_ref)
%
% Matrix-decoding eddy current correction. Corrects fully sampled calibration data
% and returns the x-ky phase correction coefficients for accelerated data.
%
% Inputs
%   dat             - K-space data, including both mux phase cycling and accelerated data.
%                     Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=mux*num_mux_cycle+num_accelerated_tpoints)].
%   p               - Parameter structure. See mux_epi_params.m for details.
%   ref             - Reference scan k-space data matrix(Dim: [Kx(=nx), Ky(=ny), Echo(=1), Slice(=nsl*mux), Coil(=nc)]),
%                     or ref.dat file name('*.dat'), or reference scan pfile name('*refscan.7').
%   pha_coe_1stpass - First-pass phase correction coefficients. See rawload_ref_data.m for details.
%   ramp_flt        - First-pass ramp sampling correction filter. See rawload_ref_data.m for details.
%   slices_ref      - Slice indices(Normally order indices, not pfile indices) in reference scan which correspond to the input 'dat'.
%                     Only needed if the type of 'ref' is reference scan pfile name.
%
% Outputs
%   dat             - K-space data with the calibration data corrected by the matrix decoding correction.
%                     Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=mux*num_mux_cycle+num_accelerated_tpoints)].
%   pha_coe_acc     - EPI x-ky phase correction coefficients for the accelerated data. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky(=ny), Echo(=nec),
%                     Slice(=nsl), SimultaneousSlice Z(=mux, indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])),
%                     Coil(=1 if p.coil_compress==true; =nc Otherwise)].
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

MUX_DIM = 5;

datsz = get_dat_sz(dat, p);
ref_type = get_ref_type(ref);

if (~exist('slices_ref', 'var') || isempty(slices_ref)) && strcmp(ref_type, 'ref_pfile')
    error('No slice indices for reference scan pfile');
end

if strcmp(ref_type, 'ksp_data')
    if p.inplane_R > 1                            % Multi-shot
        ny_pershot = size(ref, p.PE_DIM) / p.inplane_R;
        for shot = 2 : p.inplane_R
            ref(:, shot:ny_pershot:end, :, :, :) = ref(:, 1:ny_pershot:end, :, :, :); % Repeat because process_refdat.m averages over shots
        end
    end
end

% Load eddy current correction coefficients
pccoil_cal = 0;
do_quad_final_cal = true;
if p.coil_compress
    pccoil_acc = -1;                              % Average coefficients across coils to get one set of coefficients for all coils if coil compression will be used
else
    pccoil_acc = 0;                               % No averaging across coil
end
do_quad_final_acc = p.smooth_pha_coe_acc;
debug = false;
if p.debug
    if p.md_ecc
        method = 'matrix-decoding';
    else % ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit
        method = '(split)-slice-GRAPPA with odd even kernel fitting';
    end
end
switch ref_type
    case 'ref_pfile'
        if p.debug
            fprintf(['   Using reference scan pfile in ' method ' ghost correction.\n']);
        end
        [pha_coe_cal, descend_acq_ref] = rawload_ref_pfile(p.frames, slices_ref, p.coils, ref, pha_coe_1stpass, ramp_flt, pccoil_cal, do_quad_final_cal, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=nc)]
        pha_coe_acc = rawload_ref_pfile(p.frames, slices_ref, p.coils, ref, pha_coe_1stpass, ramp_flt, pccoil_acc, do_quad_final_acc, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=1 if p.coil_compress==true; =nc Otherwise)]
    case 'ksp_data'
        if p.debug
            fprintf(['   Using ecc data collected by mux sequence in ' method ' ghost correction.\n']);
        end
        pha_coe_cal = rawload_ref_data(ref, p, pha_coe_1stpass, ramp_flt, pccoil_cal, do_quad_final_cal, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=nc)]
        pha_coe_acc = rawload_ref_data(ref, p, pha_coe_1stpass, ramp_flt, pccoil_acc, do_quad_final_acc, debug); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice(num_muxed_slice*mux), Coil(=1 if p.coil_compress==true; =nc Otherwise]
end

% Eddy current correction coefficients for accelerated data
pha_coe_acc = reshape(pha_coe_acc, [p.NUM_COE_PER_KY_LINE, size(pha_coe_acc, 2), 1, size(pha_coe_acc,3)/p.mux, p.mux, size(pha_coe_acc,4)]); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [-floor(mux/2):1:ceil(mux/2)-1] for ascending acquisition; [ceil(mux/2)-1:-1:-floor(mux/2)] for descending acquisition), Coil(=1 if p.coil_compress==true; =nc Otherwise)]
if strcmp(ref_type, 'refp_file') && descend_acq_ref
    pha_coe_acc = FLIP(pha_coe_acc, MUX_DIM);  % Make the slice indices from negative to positive, for mux=3, the slice ordering is [-1, 0, 1].
end
pha_coe_acc = ifftshift(pha_coe_acc, MUX_DIM);    % Make the 0-th indexed slice the first slice, for mux=3, the slice ordering is [0, 1, -1]. pha_coe_acc Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])), Coil(=1 if p.coil_compress==true; =nc Otherwise)]

% Eddy current correction coefficients for calibration data
pha_coe_cal = reshape(pha_coe_cal, [p.NUM_COE_PER_KY_LINE, size(pha_coe_cal, 2), 1, size(pha_coe_cal,3)/p.mux, p.mux, size(pha_coe_cal,4)]); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [-floor(mux/2):1:ceil(mux/2)-1]), Coil(=nc)]
if strcmp(ref_type, 'refp_file') && descend_acq_ref
    pha_coe_cal = FLIP(pha_coe_cal, MUX_DIM);
end
pha_coe_cal = ifftshift(pha_coe_cal, MUX_DIM);    % For mux=3, the slice ordering is [0, 1, -1]. pha_coe_cal Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])), Coil(=nc)]

% Correct fully sampled calibration data
if exist('dat', 'var') && ~isempty(dat)
    eddy = get_pha_flt( - pha_coe_cal, datsz.x); % Dim: [X, Ky, Echo, Slice, SimultaneousSlice Z, Coil]
    cal_samp_ky = 1 : 1 : datsz.y; % All ky lines are sampled

    cal_dat = dat(:, :, :, :, :, p.cal_dat_tpoints);
    cal_dat_nt = length(p.cal_dat_tpoints)/p.mux;
    cal_dat = reshape(cal_dat, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux, cal_dat_nt]); % Calibration data for recon. Dim: [Kx, Ky, Echo, Slice, Coil, Kz(=p.mux), Time]
    cal_dat = permute(cal_dat, [1,2,3,4,5,7,6]); % Dim: [Kx, Ky, Echo, Slice, Coil, Time, Kz(=p.mux)]

    [cal_dat, ~] = mux_epi_md_ecc_full_kz(cal_dat, cal_samp_ky, eddy, p);

    cal_dat = reshape(permute(cal_dat, [1,2,3,4,5,7,6]), [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux*cal_dat_nt]);
    dat(:, :, :, :, :, p.cal_dat_tpoints) = cal_dat;
end

return

function [dat, pha_coe] = default_ecc(dat, p, ref, slices)
%
% Default(single-slice or slice-averaged) eddy current correction on both calibration and accelerated data.
%
% Inputs
%   dat     - K-space data to correct. Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time].
%   p       - Parameter structure. See mux_epi_params.m for details.
%   ref     - Reference scan k-space data matrix(Dim: [Kx(=nx), Ky(=ny), Echo(=1), Slice(=nsl*mux), Coil(=nc)]),
%             or ref.dat file name('*.dat'), or reference scan pfile name('*refscan.7').
%   slices  - Slice indices(Normally order indices, not pfile indices) the
%             input 'dat' corresponds to. Only needed if the type of 'ref'
%             is ref.dat file name or reference scan pfile name.
%
% Outputs
%   dat     - Corrected k-space data. Size same as the input 'dat'.
%   pha_coe - EPI x-ky phase correction coefficients. Dim: [2(0th order, 1st order), ny, nsl, nc].
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

ref_type = get_ref_type(ref);

if strcmp(ref_type, 'ref.dat_file') || strcmp(ref_type, 'ksp_data')
    datsz = get_dat_sz(dat, p);
end
if (~exist('slices', 'var') || isempty(slices)) && (strcmp(ref_type, 'ref.dat_file') || strcmp(ref_type, 'ref_pfile'))
    error('No slice indices');
end

if strcmp(ref_type, 'ref.dat_file') % Load coefficients from ref.dat file
    if p.debug
        fprintf('   Using ref.dat file in default ghosting correction.\n');
    end
    pha_coe = rawload_ref(datsz.y, p.num_slices, p.num_coils, p.frames, slices, p.coils, ref);
else                                % Calculate coefficients from reference scan k-space data or from reference scan pfile
    pha_coe_1stpass = [];           % []: No first-pass phase correction
    ramp_flt = [];                  % []: No first-pass ramp sampling correction
    pccoil = 0;                     % 0: No averaging across coil
    do_quad_final = true;           % True: Conduct smoothing along ky on the least squares fitted coefficients
    debug = false;                  % True: Display correction coefficient maps
    
    switch ref_type
        case 'ref_pfile'
            if p.debug
                fprintf('   Using reference scan pfile in default ghosting correction.\n');
            end
            pha_coe = rawload_ref_pfile(p.frames, slices, p.coils, ref, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug);
        case 'ksp_data'
            if p.debug
                fprintf('   Using ecc data collected by mux sequence in default ghosting correction.\n');
            end
            ref = ref(:, :, :, datsz.sl*floor(p.mux/2)+1 : datsz.sl*(floor(p.mux/2)+1), :); % Reference scan data for the middle band. Dim: [Kx, Ky, Echo(=1), Slice(=nsl), Coil(=nc)].
            pha_coe = rawload_ref_data(ref, p, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug);
    end
end

dat = epi_pha_correct(dat, pha_coe, p);

return

function ref_type = get_ref_type(ref)
%
% Determines the the data type for eddy current correction.
%
% Input
%   ref      - Reference scan k-space data, or '*.dat' file name, or '*refscan.7' file name.
%
% Output
%   ref_type - What the input 'ref' contains.
%              'ksp_data'     - Reference scan k-space data.
%              'ref.dat_file' - Name of ref.dat file.
%              'ref_pfile'    - Name of reference scan pfile.
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

if isnumeric(ref)
    ref_type = 'ksp_data';
else
    if strcmp(ref(end-3:end), '.dat')
        ref_type = 'ref.dat_file';
    else
        if strcmp(ref(end-1:end), '.7')
            ref_type = 'ref_pfile';
        else
            error('ref must be kspace data or *.dat file or *.7 file.');
        end
    end
end

return

function [p, v] = pmr_save_load_params(p, v)
%
% function [p, v] = pmr_save_load_params(p, v)
%
% Use this function to
%   (1) save value of a variable into the parameter structure when loading
%       data or conducting reconstruction for the actual scan.
%   (2) load previously saved variable value when conducting pseudo-multiple
%       replica simulation.
%
% Inputs
%   p - Input parameter structure. See mux_epi_params.m for details.
%   v - Input value of a variable.
%
% Outputs
%   p - Output parameter structure.
%   v - Output value of the variable.
%
% Usage example
%   [p, ref] = pmr_save_load_params(p, ref);
%   If p doesn't contain field 'pmr_ref', will set p.pmr_ref = ref;
%   If p contains field 'pmr_ref', will set ref = p.pmr_ref.
%
% (c) Kangrong Zhu      Stanford University     Jan 2015

field = inputname(2); % This returns string 'ref' for function call [p, ref] = pmr_save_load_params(p, ref).
field = ['pmr_' field];
if ~isfield(p, field)
    cmd = ['p.' field ' = v;'];
    eval(cmd);
else
    cmd = ['v = p.' field ';'];
    eval(cmd);
end

return