function fname = poisson_mask_save(us_msk, inplane_R, kydir)
%
% function fname = poisson_mask_save(us_msk, inplane_R, [kydir=BOTTOM_UP])
%
% Saves a Poisson-disc-like undersampling mask into a file.
%
% Inputs
%   us_msk    - Poisson disk undersample mask with fields 'ky', 'kz', 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   inplane_R - Inplane reduction factor.
%   kydir     - Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
%
% Output
%   fname     - File name for the saved file.
%
% (c) Kangrong Zhu      Stanford University     July 2014

TOP_DOWN   = 0;
CENTER_OUT = 1;
BOTTOM_UP  = 2;

if ~exist('kydir', 'var') || isempty(kydir)
    kydir = BOTTOM_UP;
end

dftz_type = 'inverse';
blip_polarity = (-1)^strcmp(dftz_type, 'inverse');

% -- Number of kz lines
% Unique omegaz values
unique_omegaz = blip_polarity * sort(unique(us_msk.omegaz), 'descend');
% Deal with roundoff errors
tmp = unique_omegaz;
for i1 = 1 : length(tmp)
    for i2 = i1+1 : length(tmp)
        if abs(tmp(i2)-tmp(i1)) < 10^-4
            unique_omegaz = unique_omegaz(unique_omegaz ~= tmp(i2));
        end
    end
end
% number of kz lines
cap_fov_shift = length(unique_omegaz);

% -- kz indices
if kydir == BOTTOM_UP
    us_msk.kz = flipdim(us_msk.kz(:).', 2);
end
nsamp = length(us_msk.kz);
kz_indices = zeros(nsamp, 1);
for samp = 1 : nsamp
    kz_indices(samp) = find(unique_omegaz == us_msk.omegaz(us_msk.kz(samp))); % kz_indices 1,2,3...cap_fov_shift corresponds to unique_omegaz(1), unique_omegaz(2), unique_omegaz(3)...unique_omegaz(cap_fov_shift)
end
kz_indices = kz_indices - 1; % kz_indices 0,1,2...cap_fov_shift-1 corresponds to unique_omegaz(1), unique_omegaz(2), unique_omegaz(3)...unique_omegaz(cap_fov_shift)(i.e., -kzmax ... kzmax)

% -- Save file
maxetl = nsamp;              % Maximum echo train length this undersample mask allows
fname = ['poisson_mask_arc' num2str(inplane_R) '_capfovshift' num2str(cap_fov_shift) '_maxetl' num2str(maxetl) '.raw'];
fid = fopen(fname, 'w', 'l');
fwrite(fid, inplane_R, 'int32');
fwrite(fid, cap_fov_shift, 'int32');
fwrite(fid, maxetl, 'int32');
fwrite(fid, kz_indices, 'int32');
fclose(fid);

return