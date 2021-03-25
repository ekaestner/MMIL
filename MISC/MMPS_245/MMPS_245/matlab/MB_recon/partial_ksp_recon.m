function ksp_part = partial_ksp_recon(ksp_part, ny_part, method, niter)
%
% function ksp = partial_ksp_recon(ksp_part, ny_part [, method='homodyne', niter=4] )
%
% Reconstruct partially acquired k-space.
%
% References:
% [1] John Pauly, handout on partial k-space reconstruction, Stanford EE369C course, Aut07.
% [2] Matt A. Bernstein et al., Handbook of MRI Pulse Squences, section 13.4, 2004.
% [3] G. McGibney et al., Quantitative Evaluation of Several Partial Fourier
%     Reconstruction Algorithms, MRM 30:51-59 (1993).
%
% Inputs:
%   ksp_part - Partially acquired k-space, with full matrix size. The
%              acquired k-space data is FOLLOWED by zeros in the PE direction.
%              Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase].
%   ny_part  - # of acquired ky lines in the partial k-space acquisition.
%   method   - A string specifying the algorithm to use.
%              'pocs'                     : An iterative POCS algorithm.
%              'h' or 'homo' or 'homodyne': The Homodyne algorithm.
%              'zero': zero-fill
%   niter    - Number of iterations.
%
% Output:
%   ksp      - Reconstructed k-space data, having the same size as 'ksp_part'.
%
% (c) Kangrong Zhu,     Stanford University     Sep 2012

% -- Parse inputs.
KY_DIM = 2;

if ~exist('method', 'var')
    method = 'homodyne';
end

if ~exist('niter', 'var')
    niter = 4;
end

possible_methods = {'h', 'homo', 'homodyne', 'pocs', 'zero'};
if ~ismember(method, possible_methods)
    error('The input ''method'' must be ''h'', ''homo'', ''homodyne'',  ''pocs'', or ''zero''.');
end

if strcmp(method, 'zero')
    % Nothing to do...
    return
end

kspsz   = size(ksp_part);
if length(kspsz) < 7
    kspsz(end+1 : 7) = 1;                    % To make sure kspsz(7) exists when it is called below
end
ny_full = kspsz(KY_DIM);                     % Full matrix size in PE.

if ny_full < ny_part
    error('Total # of ky lines is fewer than the # of acquired ky lines.');
end

if ny_full == ny_part
    disp('No partial k-space reconstruction is needed. The input k-space data is returned.');
    return;
end

if ny_part <= ny_full/2
    error('No symmetrically acquired ky lines to estimate the phase.');
end

% -- The merging filter and the lowpass filter in Homodyne.
DC_OFF    = 1/2;                             % For calculating the transition band, suppose the (N/2)-th / (N/2+1)-th ky lines are at k = (-/+)DC_OFF*delta_ky.
N_TRANS   = 8;                               % (N_TRANS+1) is the length of the transition band.
N_TRANS   = N_TRANS + mod(N_TRANS, 2);       % N_TRANS must be even.
W_ASYMM   = 2;                               % Weighting for the asymmetrically acquried ky lines in the Homodyne merging filter.
W_SYMM    = 1;                               % Weighting for the symmetrically acquired ky lines.

ny_asymm  = ny_full - ny_part;               % # of asymmetrically acquired ky lines.
ny_symm   = (ny_part - ny_full/2) * 2;       % # of symmetrically acquired ky lines.
nf_asymm  = ny_asymm - (N_TRANS/2+1);        % length of the flat portion of the Homodyne merging filter for the asymmetrically acquired ky lines.
nf_symm   = ny_symm - N_TRANS;               % length of the flat portion for the symmetrically acquired ky lines.
nf_not    = nf_asymm;                        % length of the flat portion for the not acquired ky lines.

k0        = ny_symm/2 - DC_OFF;              % The acquired ky lines are in the range [-kmax, k0]*delta_ky.
k_trans   = k0-N_TRANS/2 : 1 : k0+N_TRANS/2; % Ky indices for the transition band.
trans     = ( cos(pi.*(abs(k_trans)-(k0-N_TRANS/2))./2./N_TRANS) ).^2; % Amplitude for the transition band, corresponding to Eqn 13.94 in section 13.4 in the Handbook of MRI Pulse Sequences.

low_flt   = [zeros(1, nf_asymm), fliplr(trans), ones(1, nf_symm), trans, zeros(1, nf_not) ];      % Low-pass filter.
merge_flt = [W_ASYMM*ones(1, nf_asymm), W_SYMM+trans, ones(1, nf_symm), trans, zeros(1, nf_not)]; % Merging filter for Homodyne.

% RFD: loop over time points to reduce RAM usage. Very important for parallelization.
for i7 = 1 : kspsz(7)
for tp = 1:kspsz(6)
    % -- MoFIR reconstruction (MRM 30:51-59 (1993), by G. McGibney).
    %    Homodyne: merging filter   => phase correction;
    %    MoFIR   : phase correction => merging filter (performs a bit better than Homodyne).
    im_low  = zpad( ksp_part(:, ny_asymm+1 : ny_asymm+ny_symm, :, :, :, tp, i7), kspsz(1:5) );
    for fe = 1:kspsz(1), for ec = 1:kspsz(3), for sl = 1:kspsz(4), for cl = 1:kspsz(5)
        im_low(fe,:,ec,sl,cl) = im_low(fe,:,ec,sl,cl) .* low_flt;
    end; end; end; end
    im_low    = ifft2c( im_low );    % The low-passed images.

    ksp        = ifft2c(ksp_part(:,:,:,:,:,tp,i7)) .* conj(im_low) ./ abs(im_low); % Correct the image phase in the image domain.
    ksp = fft2c(ksp);
    for fe = 1:kspsz(1), for ec = 1:kspsz(3), for sl = 1:kspsz(4), for cl = 1:kspsz(5)
        ksp(fe,:,ec,sl,cl) = ksp(fe,:,ec,sl,cl) .* merge_flt;
    end; end; end; end

    ksp_conj  = conj( FLIP(FLIP(ksp, 1),2) );      % F*(-u,-v) <=> f*(x,y)
    ksp_conj  = circshift(ksp_conj, [1, 1, 0, 0, 0]);    % A k-space matrix is centered at (N/2+1, N/2+1). The reflected k-space needs to be shifted to align its center with the original one.
    ksp_conj(1, :, :, :, :, :, :) = 0;                   % Since the k-space is centered at (N/2+1, N/2+1), the samples conjugate to the current 1st kx and 1st ky lines were not acquired.
    ksp_conj(:, 1, :, :, :, :, :) = 0;

    ksp = (ksp + ksp_conj)/2;                            % (F(u,v)+F*(-u,-v))/2 <=> { (f(x,y)+f*(x,y))/2 = real(f(x,y) }. Taking the real part of the image is performed in k-space for superior speed.

    % -- The iterative POCS algorithm.
    % RFD: modified to always do four iterations regardless of mse improvement for increased speed.
    if strcmp(method, 'pocs')
        merge_flt  = [ones(1, ny_asymm + ny_symm - N_TRANS/2), trans, zeros(1, nf_not)]; % Merging filter for POCS.
        ksp         = ifft2c(ksp);
        for iter = 1 : niter
            ksp  = ksp .* im_low ./ abs(im_low);         % im .* exp(1i * phase)
            ksp  = fft2c(ksp);
            for fe = 1:kspsz(1), for ec = 1:kspsz(3), for sl = 1:kspsz(4), for cl = 1:kspsz(5)
                ksp(fe,:,ec,sl,cl) = ksp_part(fe,:,ec,sl,cl,tp) .* merge_flt + ksp(fe,:,ec,sl,cl) .* (1-merge_flt);
            end; end; end; end
            ksp   = ifft2c(ksp);
            ksp   = ksp .* conj(im_low) ./ abs(im_low);   % im .* exp(-1i * phase)
            ksp   = real(ksp);
        end
        ksp = fft2c(ksp);
    end
    ksp_part(:,:,:,:,:,tp,i7) = ksp;
end
end

return;


