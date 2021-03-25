function [ref_slices, ref_nslices] = get_ref_slices(refp, mux_slices, mux_nslices, mux, md_ecc, debug)
% function [ref_slices, ref_nslices] = get_ref_slices(refp, mux_slices, mux_nslices, mux, md_ecc, debug)
%
% Get slice indices in reference scan which correspond to the specified 'mux_slices' in the actual mux scan.
%
% Inputs
%   refp        - File name of reference scan pfile.
%   mux_slices  - Slice indices in actual mux scan (Normally ordered indices, not pfile indices).
%   mux_nslices - Total number of muxed slices in actual mux scan.
%   mux         - Number of simultaneous slices.
%   md_ecc      - True: get slice indices for matrix-decoding Eddy Current Correction(ECC); False: get slice indices for default ECC.
%   debug       - True: print debugging messages.
%
% Outputs
%   ref_slices  - Slice indices in reference scan corresponding to the input 'mux_slices' in actual mux scan (Normally ordered indices, not pfile indices).
%   ref_nslices - Total number of slices in reference scan pfile.
%
% (c) Kangrong Zhu,     Stanford University     Nov 2013

header = rawheadX(refp);            % Reference scan pfile header
ref_nslices = header.image_slquant; % slquant in image header

mux_nslices_to_exam = length(mux_slices);

switch ref_nslices
    case mux_nslices * mux          % Reference scan covered the entire volume
        if debug
            fprintf(' Reference scan: Covered the entire volume.\n');
        end
        if md_ecc                   % Matrix-decoding ghosting correction
            band_base = (0 : mux-1) .* mux_nslices;
            ref_slices = zeros(1, mux_nslices_to_exam * mux);
            for idx = 1 : length(band_base)
                ref_slices((idx-1)*mux_nslices_to_exam+1 : idx*mux_nslices_to_exam) = band_base(idx) + mux_slices; % For mux=3, ref_slices = [mux_slices, mux_slices+mux_nslices, mux_slices+2*mux_nslices]
            end
        else                        % Default ghosting correction
            ref_slices = floor((mux-1)/2)*mux_nslices + mux_slices;
        end
        
    case mux_nslices                % Reference scan covered one of the mux bands
        if debug
            fprintf(' Reference scan: Covered one of the mux bands.\n');
        end
        if md_ecc                   % Matrix-decoding ghosting correction
            error('Reference scan covered only one of the mux bands. Cannot use such pfile for matrix-decoding ECC or (split-)slice-GRAPPA with odd even kernel fitting.');
        else                        % Default ghosting correction
            ref_slices  = mux_slices;
        end
        
    otherwise
        error('Total number of slices in reference scan must be either mux_nslices or mux_nslices*mux.');
end

return