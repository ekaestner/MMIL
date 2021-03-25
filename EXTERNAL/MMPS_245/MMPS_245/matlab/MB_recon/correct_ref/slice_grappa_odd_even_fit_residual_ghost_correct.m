function dat = slice_grappa_odd_even_fit_residual_ghost_correct(dat, pha_coe, p)
%
% Correct residual Nyquist ghost for (split)-slice-GRAPPA with odd even fitting.
%
% Inputs
%   dat     - k-space data with residual Nyquist ghost. Dim: [Kx, Ky(=ny_full), Echo, Slice(=nsl*nz), Coil(=nc)]
%   pha_coe - Coefficients describing the residual Nyquist ghost. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky(=ny_part if p.partial_ky; =ny_full if ~p.partial_ky), Echo, Slice(=nsl), SimultaneousSlice Z(indices ifftshift(-floor(mux/2):1:ceil(mux/2)-1)), Coil(=1 if p.coil_compress; =nc Otherwise)]
%   p       - Parameter structure. See mux_epi_params.m for details.
%
% Output
%   dat     - k-space data with residual Nyquist ghost corrected. Dimension same as the input 'dat'.
%
% (c) Kangrong Zhu  June 2015

[nx, ny_full, nec, nsl_tot, nc, nt] = size(dat);
coesz = size(pha_coe);
if length(coesz) < 6
    coesz(end+1 : 6) = 1;
end
if (coesz(3) ~= nec) || (coesz(4)*coesz(5) ~= nsl_tot)
    error('pha_coe and data size mismatch.');
end
pha_coe = fftshift(pha_coe, 5); % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo, Slice, SimultaneousSlice Z(indices (-floor(mux/2):1:ceil(mux/2)-1)), Coil(=1 if p.coil_compress; =nc Otherwise)]
pha_coe = reshape(pha_coe, [p.NUM_COE_PER_KY_LINE, coesz(2), nec, nsl_tot, coesz(6)]);
if p.partial_ky
    pha_coe = cat(2, pha_coe, zeros(p.NUM_COE_PER_KY_LINE, ny_full-coesz(2), nec, nsl_tot, coesz(6)));
end
if size(pha_coe, p.C_DIM) == 1
    pha_coe = repmat(pha_coe, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nc]);
end
dat = epi_pha_correct(dat, pha_coe, p);

return