function [ksp, pha_flt] = epi_pha_correct(ksp, pha_coe, p)
%
% function [ksp, pha_flt] = epi_pha_correct(ksp, pha_coe, p)
%
% Correct the 0-th and 1st order phases in the x-ky space for EPI data.
% Correspondingly, the constant phases of and the misalignment between the
% odd and even echoes in the k-space are corrected.
%
% Inputs
%   ksp     - Uncorrected k-space data. Dim: [Kx, Ky, Echo, Slice, Coil, Time]
%   pha_coe - Coefficients for x-ky phase correction. See get_pha_flt.m for details.
%             Dim: [2(0th order, 1st order), other dimensions(e.g. ny, nsl, nc].
%   p       - Parameter structure. See mux_epi_params.m for details. The following
%             fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM.
%
% Outputs
%   ksp     - Corrected k-space data (note that the correction is done in-place).
%   pha_flt - X-ky phase filter. xKy_corrected = xKy_uncorrected .* pha_flt.
%             Dim: [X, Ky, Echo(=1), Slice, Coil].
%
% (c) Kangrong Zhu,     Stanford University     July 2012

datsz = get_dat_sz(ksp, p);
pha_flt = get_pha_flt(pha_coe, datsz.x);                              % Dim: [X, Ky, Slice, Coil]
pha_flt = reshape(pha_flt, [datsz.x, datsz.y, 1, datsz.sl, datsz.c]); % Dim: [X, Ky, Echo(=1), Slice, Coil]
for echo = 1 : datsz.ec
    for t = 1 : datsz.t
        ksp(:, :, echo, :, :, t) = fftc( ifftc(ksp(:, :, echo, :, :, t), p.FE_DIM) .* pha_flt, p.FE_DIM); % The same coefficients are used for every echo and every time point.
    end
end

return