function [reconed, sz_out] = mux_recon_1Dgrappa(ksp, p)
%
% function [reconed, sz_out] = mux_recon_1Dgrappa(ksp, p)
% 
% Reconstruct slice-multiplexed data using 1D GRAPPA.
% The last group of mux phase cycling is used as the reference images.
%
% Referenes
%   1. Blaimer M. et al. JMRI 2006; 24(2): 444-450.
%   2. Blaimer M. et al. MRM 2013; 69(4): 974-980.
%
% Inputs
%   ksp     - Slice-multiplexed k-space data. The first 'p.mux_encoded
%             *p.num_mux_cycle' time points are the mux phase cycling time points.
%             Dim: [Kx, Ky, Echo, Slice(with slice multiplexing), Coil, Time].
%   p       - Parameter structure. See mux_epi_params.m for details.
%             The following fields are used in this function: mux_kersz_1d_grappa,
%             mux_acssz_1d_grappa, mux, num_mux_cycle, PE_DIM, T_DIM, KEEP_ORIG_SZ,
%             caipi, cap_fov_shift, grappa_domain, zpad_1Dgrappa, debug.
%
% Outputs
%   reconed - Reconctructed data. This is k-space data if p.grappa_domain is
%             'kSpace', and is image space data if p.grappa_domain is 'imSpace'.
%             Dim: [FE, PE, Echo, Slice(solved for slice multiplexing), Coil, Time].
%   sz_out  - Size of 'reconed'. A structure array with fields: x, y, ec, sl, c, t, kz.
%
% (c) Kangrong Zhu,     Stanford University     July 2012

% -- Variables deducted from inputs.
if p.caipi
    if (abs(p.cap_fov_shift) ~= p.mux) && strcmp(p.zpad_1Dgrappa, 'minZpad');
        p.zpad_1Dgrappa = 'normalZpad';
        if p.debug
            fprintf('    mux_recon_1Dgrappa: Minimum zero padding does not work when abs(p.cap_fov_shift)~=p.mux. Changed to normal zero padding');
        end
    end
    if strcmp(p.zpad_1Dgrappa, 'minZpad')   % Use minimum zero padding, only zero pad a total of 1*datsz.y pixels in the PE direction.
        mux_R = p.mux + 1;                  % Reduction factor for the 1D GRAPPA recon.
    else                                    % Use an easy way for zero padding, zero pad (p.mux-1)*datsz.y pixels in the PE direction.
        mux_R = 1 + (p.mux-1)*2;
    end
else
    mux_R = p.mux;
end
datsz = get_dat_sz(ksp, p);
nt_solve = datsz.t - p.mux*p.num_mux_cycle; % # of temporal phases in 'ksp' to solve for.

if p.caipi
    if strcmp(p.zpad_1Dgrappa, 'minZpad')   % Minimum zero padding
        ny_zpad = floor(datsz.y/p.mux);                       % # of y lines in the zero matrices for zero padding. Zero padding using matrices can only account for shifts with interger number of pixels.
        ny_zpad_residual_total = datsz.y - p.mux*ny_zpad;     % # of y lines lacking after the zero matrix padding.
        ny_zpad_residual_each = ny_zpad_residual_total/p.mux; % The difference between the fractional number and the floored integer. This will be used to correct shifts with fractional number of pixels.
        slice_shift_amt = sign(p.cap_fov_shift) * (0:p.mux-1) * ny_zpad_residual_each/(datsz.y+ny_zpad); % The amount of shift to apply to each slice, as a fraction of the FOV in y. This will be used to correct shifts with fractional number of pixels.
    else                                    % Easy way for zero padding
        ny_zpad = datsz.y;                                    % # of y lines in the zero matrices for zero padding.
        slice_indices = - floor(p.mux/2) : 1 : (ceil(p.mux/2)-1);       % Slice indices, must center around index 0. [-1,0,1] for mux=3, [-2,-1,0,1] for mux=4. ([-1,0,1] and [2,0,1] are the same when abs(cap_fov_shift)==p.mux, but different when abs(cap_fov_shift)~=p.mux).
        slice_indices = ifftshift(slice_indices);                       % Make the first slice the 0-th indexed slice.
        slice_shift_amt_raw = mod(slice_indices, abs(p.cap_fov_shift)); % [0,1,2] for p.mux=3 and abs(p.cap_fov_shift)=3; [0,1,1] for p.mux=3 and abs(p.cap_fov_shift)=2; [0,1,2,3] for p.mux=4 and abs(p.cap_fov_shift)=4; [0,1,1,2] for p.mux=4 and abs(p.cap_fov_shift)=3.
        slice_shift_amt = slice_shift_amt_raw /(p.cap_fov_shift*2);     % The amount of shift to apply to each slice, as a fraction of the FOV in y. This will be used to synthesize the shifting for the reference image.
    end
end

% -- Check inputs.
if p.mux <= 1
    error('No slice multiplexing to solve.');
end

if nt_solve < 0
    error('The number of temporal phases to solve after mux phase cycling can''t be negative.');
end

% -- The temporal phases for mux phase cycling.
ref_dat = ksp(:, :, :, :, :, 1 : p.mux_encoded*p.num_mux_cycle);
ref_dat = reshape(ref_dat, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux_encoded, p.num_mux_cycle]);
ref_dat = mux_dftz(ref_dat, p.T_DIM, p.cap_fov_shift_cal, p.mux, 'decode'); % Dim: [Kx, Ky, Echo, Slice, Coil, SimultaneousSlice, MultiplexingPhaseCycle]
if (p.pseq == 2) && mod(p.mux, 2)==0 && (~p.caipi || (p.caipi && ~p.cal_gzblips))         % Adjust the slice ordering for the non-caipi case, if the number of simultaneous bands is even. Somehow this is needed for the muxarcepi GRE sequence, but not for the muxarcepi2 sequence.
    ref_dat = circshift(ref_dat, [0, 0, 0, 0, 0, -1, 0]);
end

% -- GRAPPA reconstruction.
reconed = zeros(datsz.x, datsz.y, datsz.ec, datsz.sl, p.mux, datsz.c, nt_solve);

if nt_solve > 0
    % Raw reference images.
    refic_raw = ifft2c(ref_dat(:, :, :, :, :, :, p.num_mux_cycle)); % Dim: [FE, PE, Echo, Slice, Coil, SimultaneousSlice, MultiplexingPhaseCycle(=1)]
    if p.caipi
        if p.cap_fov_shift > 0
            refic_raw = cat(p.PE_DIM, zeros(datsz.x, ny_zpad, datsz.ec, datsz.sl, datsz.c, p.mux), refic_raw);
        else % p.cap_fov_shift <0
            refic_raw = cat(p.PE_DIM, refic_raw, zeros(datsz.x, ny_zpad, datsz.ec, datsz.sl, datsz.c, p.mux));
        end
        refic_raw = mux_shift(refic_raw, p.T_DIM, slice_shift_amt, 'im'); % Note the negative y direction in matlab is the positive y direction in the physical image space. If (CAIPI && Minimum zero padding), this step is the 1st step to correct the shifts with fractional number of pixels.
    end
    
    % The position for the ACS data.
    acspos = grappa_acspos( struct('x',{datsz.x},'y',{datsz.y*mux_R}), p.mux_acssz_1d_grappa,...
        struct('x',{[1,datsz.x]}, 'y',{[1,datsz.y*mux_R]}), p.debug);
    
    % Recon.
    for slice = 1 : datsz.sl
        for echo = 1 : datsz.ec
            refic = refic_raw( :, :, echo, slice, :, :); % Dim: [FE, PE, Echo, Slice, Coil, SimultaneousSlice]
            refic = permute(refic, [1,2,6,5,3,4]);       % Dim: [FE, PE, SimultaneousSlice, Coil]
            if p.caipi && strcmp(p.zpad_1Dgrappa, 'minZpad') && p.cap_fov_shift > 0
                refic = FLIP(refic, 3);
            end
            refic = reshape( refic, [datsz.x, size(refic,2)*size(refic,3), datsz.c]); % Dim: CAIPI && (Minimum zero padding): [nx, ny*mux+(ny-ny_zpad_residual_total), nc]; CAIPI && (easy way for zero padding): [nx, ny*mux*2, nc]; Non-CAIPI: [nx, ny*mux, nc].
            if p.caipi
                if strcmp(p.zpad_1Dgrappa, 'minZpad')    % Minimum zero padding
                    if p.cap_fov_shift > 0
                        refic = cat(p.PE_DIM, zeros(datsz.x, ny_zpad_residual_total, datsz.c), refic); % 2nd step to correct the shifts with fractional number of pixels.
                    else % p.cap_fov_shift < 0
                        refic = cat(p.PE_DIM, refic, zeros(datsz.x, ny_zpad_residual_total, datsz.c));
                    end
                else                                     % Use an easy way for zero padding
                    if p.cap_fov_shift > 0
                        refic = refic(:, datsz.y+1:end, :);
                    else % p.cap_fov_shift < 0
                        refic = cat(p.PE_DIM, refic(:, 1:datsz.y, :), refic(:, datsz.y*2+1:end, :));
                    end
                end
            end
            refksp = fft2c(refic); 
            
            acs_dat = refksp(acspos.x, acspos.y, :);
            
            % Calculate the kernel.
            ker = grappa_kernel(acs_dat, [p.mux_kersz_1d_grappa.x, p.mux_kersz_1d_grappa.y], ...
                mux_R, p.grappa_domain, [datsz.x, datsz.y*mux_R]);
            
            % Solve the multiplexed slices for each temporal frame.
            for tphase = 1 : nt_solve
                us_dat = zeros(datsz.x, datsz.y*mux_R, datsz.c);
                us_dat(:, 1:mux_R:end, :) = reshape( ksp(:, :, echo, slice, :, p.mux*p.num_mux_cycle+tphase), ...
                    [datsz.x, datsz.y, datsz.c]);        % Dim: [FE, PE, Coil]
                
                if mod(mux_R, 2) == 0                    % For even 'mux_R'(must be a non-caipi case or (CAIPI && Minimum zero padding) case), need to adjust the aliasing pattern to match the concatenated reference image.
                    us_dat(:, 1 : mux_R*2 : end, :) = -us_dat(:, 1 : mux_R*2 : end, :); % Dim: [FE, PE, Coil]
                end
                
                rec_dat = grappa_recon(us_dat, ker, p.grappa_domain); % Reconed data, Dim: [FE, PE, Coil]
                
                if p.caipi
                    if strcmp(p.grappa_domain, 'kSpace')
                        rec_dat = ifft2c(rec_dat);                    % Reconstructed coil images
                    end
                    
                    if strcmp(p.zpad_1Dgrappa, 'minZpad')% Minimum zero padding
                        if p.cap_fov_shift >0
                            rec_dat = rec_dat(:, ny_zpad_residual_total+1 : end, :); % Undo the effects of the 2nd step to correct the shifts with fractional number of pixels.
                        else % p.cap_fov_shift <0
                            rec_dat = rec_dat(:, 1 : end-ny_zpad_residual_total, :);
                        end
                        rec_dat = reshape(rec_dat, [datsz.x, datsz.y + ny_zpad, p.mux, datsz.c]);
                        if p.cap_fov_shift > 0
                            rec_dat = FLIP(rec_dat, 3);
                        end
                    else                                 % Use an easy way for zero padding
                        if p.cap_fov_shift > 0
                            rec_dat = cat(p.PE_DIM, zeros(datsz.x, datsz.y, datsz.c), rec_dat);
                        else % p.cap_fov_shift < 0
                            rec_dat = cat(p.PE_DIM, cat(p.PE_DIM, rec_dat(:, 1:datsz.y, :), zeros(datsz.x, datsz.y, datsz.c)), rec_dat(:, datsz.y+1:end, :));
                        end
                        rec_dat = reshape(rec_dat, [datsz.x, 2*datsz.y, p.mux, datsz.c]);
                    end
                    
                    rec_dat = fftshift(mux_shift(rec_dat, 3, -slice_shift_amt, 'im'), 3);           % If minimum zero padding, this undos the effects of the 1st step to correct the shifts with fractional number of pixels.
                    
                    if p.cap_fov_shift > 0
                        reconed(:, :, echo, slice, :, :, tphase) = rec_dat(:, ny_zpad+1:end, :, :); % reconed Dim: [FE, PE, Echo, Slice, Mux, Coil, TemporalPhase]
                    else % p.cap_fov_shift < 0
                        reconed(:, :, echo, slice, :, :, tphase) = rec_dat(:, 1:datsz.y, :, :);
                    end
                    if strcmp(p.grappa_domain, 'kSpace')
                        reconed(:, :, echo, slice, :, :, tphase) = fft2c(reconed(:, :, echo, slice, :, :, tphase));
                    end
                else     % Non-CAIPI
                    if strcmp(p.grappa_domain, 'imSpace')
                        reconed(:, :, echo, slice, :, :, tphase) = ...                              % reconed Dim: [FE, PE, Echo, Slice, Mux, Coil, TemporalPhase]
                            fftshift(reshape(rec_dat, [datsz.x, datsz.y, p.mux, datsz.c]), 3);
                    else % p.grappa_domain is 'kSpace'
                        reconed(:, :, echo, slice, :, :, tphase) = ...
                            fft2c( fftshift(reshape(ifft2c(rec_dat), [datsz.x, datsz.y, p.mux, datsz.c]), 3) );
                    end
                end

            end
        end
    end
end

reconed = reconed .* sqrt(p.mux/mux_R) ./ sqrt(p.mux);  % sqrt(p.mux/mux_R): Adjust scaling in ifft2c to account for the upsampling. The scaling in ifft2c should be consistent with the scaling in fft2c; sqrt(p.mux): Account for the fact that the slice direction is not using orthonormal DFT encoding.
ref_dat = permute(fftshift(ref_dat,6), [1, 2, 3, 4, 6, 5, 7]);
if strcmp(p.grappa_domain, 'imSpace')
    ref_dat = ifft2c(ref_dat);
end
reconed = cat(7, ref_dat, reconed);

reconed = reshape( reconed, [datsz.x, datsz.y, datsz.ec, datsz.sl*p.mux, ...
    datsz.c, p.num_mux_cycle+nt_solve]); % Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase]

sz_out = get_dat_sz(reconed, p);

return