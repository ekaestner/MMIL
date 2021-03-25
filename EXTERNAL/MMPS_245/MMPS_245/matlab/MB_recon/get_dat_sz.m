function sz = get_dat_sz(dat, p)
%
% function sz = get_dat_sz(dat, p)
%
% Returns the size of the input matrix 'dat' in the structure array 'sz'.
%
% Inputs
%   dat - Input data. Dim: [FE, PE, Echo, Slice, Coil, Time, Kz].
%   p   - Parameter structure. See mux_epi_params.m for details. The following
%         fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM, KZ_DIM.
%
% Output
%   sz  - A structure array with the following fields: x(FE), y(PE),
%         ec(Echo), sl(Slice), c(Coil), t(Time), kz(Kz).
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

sz = struct('x', {size(dat, p.FE_DIM)}, 'y', {size(dat,p.PE_DIM)}, ...
    'ec', {size(dat, p.EC_DIM)}, 'sl', {size(dat, p.SL_DIM)}, ...
    'c', {size(dat, p.C_DIM)}, 't', {size(dat, p.T_DIM)}, ...
    'kz', {size(dat, p.KZ_DIM)});
