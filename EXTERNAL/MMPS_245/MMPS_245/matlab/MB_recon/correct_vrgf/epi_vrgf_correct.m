function dat_out = epi_vrgf_correct(dat, ramp_flt)
% function dat_out = epi_vrgf_correct(dat, ramp_flt)
%
% Correct ramp sampling in EPI data. Correction will be conducted image by image.
%
% Inputs
%   dat      - K-space data with ramp sampling. Dim: [Kx(=nxi), Ky, ...]
%   ramp_flt - Correction filter. When not exist or empty, no correction will be conducted. Dim: [nxo, nxi]
%
% Output
%   dat_out  - Corrected k-space data. Dim: [Kx(=nxo), Ky, ...]
%
% (c) Kangrong Zhu,     Stanford University     Sep 2013

if ~exist('ramp_flt', 'var') || isempty(ramp_flt)
    dat_out = dat;
    return;
end

[nxo, nxi] = size(ramp_flt); % [Prescribed size in X, Ramp sampled size in X]

sz = size(dat);
if sz(1) ~= nxi
    error('Data size in x doesn''t match the vrgf filter.');
end

szn = sz;
szn(1) = nxo;

dat_out = zeros(szn);
for idx = 1 : prod(szn(3:end))
    dat_out(:, :, idx) = ramp_flt * dat(:, :, idx);
end