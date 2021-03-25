function [dat_out, ranks] = mux_epi_md_ecc_full_kz(dat_mux, samp_ky, eddy, p)
%
% function [dat_out, ranks] = mux_epi_md_ecc_full_kz(dat_mux, samp_ky, eddy, p)
% 
% Mux EPI matrix-decoding eddy current correction for data with fully sampled kz.
% The input data can be fully sampled 3D slice-muxed data or undersampled slice-muxed
% data with fully sampled kz and undersampled ky. Currently only works for DFT encoding along z.
%
% Inputs
%   dat_mux - Slice-multiplexed k-space data. Dim: [Kx(=nx), Ky(=ny), 
%             Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt), Kz(=nz)].
%   samp_ky - A vector of length 'nsamp_ky'(nsamp_ky <= ny), specifying the
%             acquired ky indices(range 1,2...ny).
%   eddy    - The phase terms induced by the eddy current effects,
%             xKyData-with-eddy-current = xKyData-without-eddy-current .* eddy.
%             Use rawload_ref_pfile.m(or rawload_ref.m or rawload_ref_data.m) and
%             get_pha_flt.m to generate this(this should be reshaped from the output 
%             'pha_flt' of get_pha_flt.m). If size(eddy, 7)==1, the input matrix 'eddy' 
%             will be used for all kz lines. If 'eddy' is empty: The eddy
%             current effect correction will not be included. Dim: [X(=nx), Ky(=nsamp_ky),
%             Echo(=nec), Slice(=nsl), SilmultaneousSlice Z(=nz), Coil(=nc), Kz(=nz or =1)].
%   p       - Parameter structure. See mux_epi_params.m for details.
%             The following fields are used in this function: FE_DIM(=1),
%             PE_DIM(=2), EC_DIM(=3), SL_DIM(=4), C_DIM(=5), T_DIM(=6),
%             KZ_DIM(=7), KEEP_ORIG_SZ, debug, cap_fov_shift_cal.
%
% Outputs
%   dat_out - Slice-multiplexed data with eddy current effects corrected.
%             Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc),
%             Time(=nt), Kz(=nz)].
%   ranks   - The rank of all encoding matrices which were (pseudo)inverted
%             in the decoding process. Dim: [X(=nx), Ky(=nsamp_ky), 
%             Echo(=nec), Slice(=nsl), Coil(=nc)].
%
% (c) Kangrong Zhu,     Stanford University    April 2013

if size(dat_mux, p.KZ_DIM) ~= p.mux_encoded
    error('data size in kz and p.mux_encoded mismatch)');
end

NT_THRESH = 300;                                              % If there are more time points than this, solve time point by time point; otherwise solve all time points together.
calc_rank = p.debug;                                          % True: Calculate and record the rank of each encoding matrix
datsz     = get_dat_sz(dat_mux, p);
nkz       = size(dat_mux, p.KZ_DIM);                          % Define nkz, nz, nsamp_ky to make the later code easier to read
nz        = nkz;
nsamp_ky  = length(samp_ky);                                  % Number of acquired points along ky

include_eddy_correct = (~isempty(eddy));                      % True: Include eddy current correction in the decoding process
if include_eddy_correct
    eddy_kz_sz = size(eddy, p.KZ_DIM);
    if (eddy_kz_sz ~= nkz) && (eddy_kz_sz ~= 1)
        error('Wrong size in the kz dimension of the input ''eddy'' matrix.');
    end
end

dat_mux = ifftc(dat_mux, p.FE_DIM);                           % Slice multiplexed data in x-ky-kz space. Dim: [X(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt), Kz(=nz)]

dat_out = zeros( size(dat_mux) );                             % Dim: [X(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt), Z(=nz)]
if calc_rank
    ranks = zeros(datsz.x, nsamp_ky, datsz.ec, datsz.sl, datsz.c);
else
    ranks = [];                                               % Avoid unnecessary computing when not in debug mode
end

dftz_mtx = encode_dftz_mtx(p.cap_fov_shift_cal, nz);             % The DFTz encoding part in the encoding matrix. Dim: [nkz, nz]

for echo = 1 : datsz.ec
    for slice = 1 : datsz.sl
        for coil = 1 : datsz.c
            if p.debug && (mod(coil-1, 8) == 0)
                fprintf('   mux_epi_md_ecc_full_kz: Echo %d/%d, Slice %d/%d, Coil %d/%d...\n', ...
                        echo, datsz.ec, slice, datsz.sl, coil, datsz.c);
            end
            
            for x = 1 : datsz.x
                for ky = 1 : nsamp_ky
                    % The eddy current part in the encoding matrix
                    if include_eddy_correct
                        eddy_mtx = reshape(eddy(x, ky, echo, slice, :, coil, :), [nz, eddy_kz_sz]); % Dim: [Z(=nz), Kz(=nkz or =1)]
                        if eddy_kz_sz == 1                    % The same eddy current effects for all kz
                            eddy_mtx = repmat(eddy_mtx, [p.KEEP_ORIG_SZ, nkz]); % Dim: [nz, nkz]
                        end
                        eddy_mtx = permute(eddy_mtx, [2, 1]); % Dim: [nkz, nz]
                    else
                        eddy_mtx = ones(size(dftz_mtx));
                    end
                    
                    % Encoding matrix including (1)DFTz and (2) eddy current effects
                    encode_mtx = dftz_mtx .* eddy_mtx;        % Dim: [nkz, nz]
                    
                    % Decoding matrix
                    solve_mtx = inv(encode_mtx);              % Dim: [nz, nkz]
                    if calc_rank
                        ranks(x, ky, echo, slice, coil) = rank(encode_mtx);
                    end
                    
                    % Decode
                    if datsz.t <= NT_THRESH
                        % Measured signal
                        measured_sig = dat_mux(x, samp_ky(ky), echo, slice, coil, :, :); % Dim: [1, 1, 1, 1, 1, nt, nkz]
                        measured_sig = permute(reshape(measured_sig, [datsz.t, nkz]), [2,1]);  % Dim: [nkz, nt]
                        
                        dat_out(x, samp_ky(ky), echo, slice, coil, :, :) = permute(solve_mtx * measured_sig, [2,1]); % solve_mtx * measured_sig Dim: [nz, nt]
                    else
                        for time = 1 : datsz.t
                            % Measured signal
                            measured_sig = dat_mux(x, samp_ky(ky), echo, slice, coil, time, :);
                            measured_sig = reshape(measured_sig, [nkz, 1]);                    % Dim: [nkz, 1]
                            
                            dat_out(x, samp_ky(ky), echo, slice, coil, time, :) = solve_mtx * measured_sig; % solve_mtx * measured_sig Dim: [nz, 1]; dat_out Dim: [X(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt), Z(=nz)]
                        end % time
                    end % if datsz.t <= NT_THRESH, else
                    
                end % ky
            end % x
        end % coil
    end % slice
end % echo

% Change x-z dimensions in dat_out into kx-kz dimensions
dat_out = mux_dftz(dat_out, p.KZ_DIM, p.cap_fov_shift_cal, nz, 'encode'); % z => Kz
dat_out = fftc(dat_out, p.FE_DIM);                            % x => Kx. dat_out Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt), Kz(=nz)]

return