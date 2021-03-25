function pha_coe = rawload_ref_data(dat, p, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug)
%
% function pha_coe = rawload_ref_data(dat, p, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug)
%
% Calls process_refdata.m to calculate EPI x-ky phase correction coefficients from reference scan k-space data.
%
% Inputs
%   dat             - Reference scan k-space data. Dim: [Kx(=nx), Ky(=ny), Echo(=nec), Slice(=nsl), Coil(=nc)].
%   p               - Parameter structure. See mux_epi_params.m for details. Fields
%                     used in this function: KEEP_ORIG_SZ, nshot(if not exist,
%                     use inplane_R), frames, coils. Fields used in process_refdata.m:
%                     debug, pccoil, pcslice, nshot. Fields used in epi_pha_correct:
%                     FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM, KZ_DIM.
%   pha_coe_1stpass - First-pass phase correction coefficients. If not empty,
%                     correct the x-ky space data with these coefficients before
%                     fitting the output coefficients 'pha_coe'. Dim: [2(0thOrderPhase, 1stOrderPhase),
%                     Ky, Slice(=nsl OR =1(the same coefficients will be used for all slices)),
%                     Coil(=nc OR =1(the same coefficients will be used for all coils))].
%                     If empty(Default), skip this step.
%   ramp_flt        - First-pass ramp sampling correction filter.
%                     If not empty, it is used to correct ramp sampling in the k-space data after the
%                     first-pass phase correction and before fitting the output coefficients 'pha_coe'.
%                     If empty(Default), skip this step.
%   pccoil          - 0(Default):      No averaging across coil; 
%                     >=1 && <=ncoils: Use one of the coils' coefficients for all coils; 
%                     -1:              Average across coils.
%   do_quad_final   - True(Default): Conduct quadsmooth along ky after the x-ky phase correction coefficients are calculated by least squares fitting.
%                     False:         No smoothing along ky.
%   debug           - If exist and not empty, this will overwrite p.debug.
%                     p.debug - True: Display per slice per coil static phase correction coefficients. Default: False.
%
% Output
%   pha_coe         - EPI x-ky phase correction coefficients.
%                     Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice, Coil].
%
% (c) Kangrong Zhu,     Stanford University     Sep 2013

%% Parameters
if ~exist('pha_coe_1stpass', 'var'); pha_coe_1stpass = []; end;
if ~exist('ramp_flt', 'var');        ramp_flt = [];        end;
if ~exist('pccoil', 'var');          pccoil = 0;           end;
if ~exist('do_quad_final', 'var') || isempty(do_quad_final)
    do_quad_final = true;
end

if exist('debug', 'var') && ~isempty(debug)
    p.debug = debug;
end
p.pccoil = pccoil;
p.do_quad_final = do_quad_final;
if ~isfield(p, 'nshot')          % Number of shots
    if isfield(p, 'inplane_R')
        p.nshot = p.inplane_R;
    else
        error('No info for number of shots.');
    end
end
[nx, ny, nec, nsl, nc] = size(dat);
if nec ~= 1
    error('There are %d echos in the reference scan k-space data. There should be only one echo.', nec);
end

%% Prepare data
% Conduct a 1st-pass phase correction in x-ky space if the input filter is nonempty
if ~isempty(pha_coe_1stpass)
    switch size(pha_coe_1stpass, 3)
        case 1
            pha_coe_1stpass = repmat(pha_coe_1stpass, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nsl, p.KEEP_ORIG_SZ]);
        case nsl
            
        otherwise
            error('Number of slices in ''pha_coe_1stpass'' must be either 1 or be equal to number of slices in the data.');
    end
    switch size(pha_coe_1stpass, 4)
        case 1
            pha_coe_1stpass = repmat(pha_coe_1stpass, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nc]);
        case nc
            
        otherwise
            error('Number of coils in ''pha_coe_1stpass'' must be either 1 or be equal to number of coils in the data.');
    end
    dat = epi_pha_correct(dat, pha_coe_1stpass, p);
end

% Correct ramp sampling if the input filter is nonempty
if ~isempty(ramp_flt)
    dat = epi_vrgf_correct(dat, ramp_flt);
end

% Permute data
dat = permute(dat, [2,1,4,5,3]); % Dim: [Ky, Kx, Slice, Coil, Echo(=1)]

%% Calculate phase correction coefficients
if p.kydir == p.BOTTOM_UP
    dat = flipdim(dat, 1);
end

pha_coe= process_refdata(dat, p);

if p.kydir == p.BOTTOM_UP
    pha_coe = flipdim(pha_coe, 2);
end

%% Output desired frames and coils
if isempty(p.frames)
    p.frames = 1 : ny;
end
if isempty(p.coils)
    p.coils = 1 : size(pha_coe, 4);
end
pha_coe = pha_coe(:, p.frames, :, p.coils);

return