function pha_flt = get_pha_flt(pha_coe, nx)
% function pha_flt = get_pha_flt(pha_coe, nx)
%
% Calculates the phase filter in the x-ky space for EPI ghosting correction.
%
% Inputs
%   pha_coe - Coefficients for x-ky phase correction, loaded from ref.dat 
%             file or calculated from reference scan p-file or reference
%             scan k-space data. Dim: [num_coe_per_ky_line(0th order, 1st order, 2nd order...),
%             other dimensions(e.g. ny, nsl, nc)].
%   nx      - Matrix size in x.
%
% Output
%   pha_flt - X-ky phase filter. xKy_corrected = xKy_uncorrected .* pha_flt.
%             Dim: [nx, other dimensions(same as other dimensions in the input pha_coe)].
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

sz = size(pha_coe);
num_coe_per_ky_line = sz(1);                                        % Number of coefficients for each ky line.

xidx = (-nx/2 : 1 : (nx/2-1)).';
pix_idx = zeros(nx, num_coe_per_ky_line);                           % Dim: [X, PixelIndex(0thOrder-1stOrder-2ndOrder...)]
for order = 0 : 1 : num_coe_per_ky_line-1
    pix_idx(:, order+1) = xidx.^order;
end
pha_coe = reshape(pha_coe, [num_coe_per_ky_line, prod(sz(2:end))]); % Dim: [PhaseCoefficients(0thOrder-1stOrder-2ndOrder...), OtherDimensions(e.g. Ky-Slice-Coil)]
pha_flt = pix_idx * pha_coe;                                        % Dim: [X, OtherDimensions(e.g. Ky-Slice-Coil)]
pha_flt = exp(1i * pha_flt);
pha_flt = reshape(pha_flt, [nx, sz(2:end)]);                        % Dim: [X, OtherDimensions(e.g. Ky, Slice, Coil)]
