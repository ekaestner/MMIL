function mtx = us_msk_get_mtx(us_msk, ny, kydir, nomegaz)
%
% function mtx = us_msk_get_mtx(us_msk, ny, [kydir=BOTTOM_UP], [nomegaz])
%
% Get matrix representation of undersample masks.
%
% Inputs
%   us_msk  - Ky-omegaz undersample mask with fields 'ky', 'kz' and 'omegaz'. nmsk=length(us_msk). See get_ky_omegaz_us_msk.m for details.
%   ny      - Full prescribed matrix size in ky.
%   kydir   - Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
%   nomegaz - Number of samples along omegaz axis.
%             If not empty, will define an omegaz grid with 'nomegaz' lines.
%             If empty, will use default number of omegaz lines (number of unique omegaz values for uniform sampling along omegaz; ny for nonuniform sampling along omegaz) to define the grid.
%
% Output
%   mtx     - Output structure with length nmsk, with fields:
%             'mtx': Matrix representation of undersample masks. Dim: [ny, nomegaz]
%             'omegay': DFTy encoding frequencies corresponding to the 1st dimension of 'mtx'. Dim: a vector with length ny.
%             'omegaz': DFTz encoding frequencies corresponding to the 2nd dimension of 'mtx'. Dim: a vector with length nomegaz.
%
% (c) Kangrong Zhu      Stanford University     June 2014

TOP_DOWN   = 0;
CENTER_OUT = 1;
BOTTOM_UP  = 2;

if ~exist('kydir', 'var')
    kydir = BOTTOM_UP;
end
if ~exist('nomegaz', 'var')
    nomegaz = [];
end

nmsk = length(us_msk);
nsamp = length(us_msk(1).ky);
mtx(nmsk) = struct('mtx', {[]}, 'omegay', {[]}, 'omegaz', {[]});
for msk_idx = 1 : nmsk
    % Omegay grid
    mtx(msk_idx).omegay =  - 2*pi * (-ny/2 : 1 : ny/2-1) / ny;
    
    % Omegaz grid
    if isempty(nomegaz) || (nomegaz <= 0) % Use default nomegaz
        omegaz_all = sort(unique(us_msk(msk_idx).omegaz));
        
        if length(omegaz_all) == 1
            nomegaz = 1;
            mtx(msk_idx).omegaz = omegaz_all;
        else
            % Check whether the sampling along the omegaz axis is uniform or not
            unique_delta_omegaz = unique(diff(omegaz_all));
            
            % Deal with roundoff errors
            tmp = unique_delta_omegaz;
            for i1 = 1 : length(tmp)
                for i2 = i1+1 : length(tmp)
                    if abs(tmp(i2)-tmp(i1)) < 10^-4
                        unique_delta_omegaz = unique_delta_omegaz(unique_delta_omegaz ~= tmp(i2));
                    end
                end
            end
            
            % Set up omegaz grid
            if length(unique_delta_omegaz) == 1 % Sampling along omegaz axis is uniform
                nomegaz = length(omegaz_all);
                mtx(msk_idx).omegaz = omegaz_all;
            else % Sampling along omegaz axis is nonuniform
                nomegaz = ny;
                mtx(msk_idx).omegaz = get_omegaz_from_nomegaz(nomegaz);
            end
        end
    else % Use specified nomegaz
        mtx(msk_idx).omegaz = get_omegaz_from_nomegaz(nomegaz);
    end
    
    % Set up the undersampling mask matrix
    mtx(msk_idx).mtx = zeros(ny, nomegaz);
    samp_order = 1 : 1 : nsamp;
    if kydir == BOTTOM_UP
        samp_order = flipdim(samp_order, 2);
    end
    for samp = 1 : nsamp
        ky_idx = us_msk(msk_idx).ky(samp);
        omegaz_this_samp = us_msk(msk_idx).omegaz(us_msk(msk_idx).kz(samp));
        omegaz_diff = abs(mtx(msk_idx).omegaz - omegaz_this_samp);
        [~, kz_idx] = min(omegaz_diff);
        mtx(msk_idx).mtx(ky_idx, kz_idx) = samp_order(samp);
    end
end

return

function omegaz = get_omegaz_from_nomegaz(nomegaz)

min_omegaz = - (2*pi/nomegaz) * (nomegaz-1) * (1/2);
omegaz = min_omegaz : (2*pi/nomegaz) : pi;

return