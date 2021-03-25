function [us_msk, inplane_R, cap_fov_shift, maxetl] = poisson_mask_read(fname, nky, p)
%
% function [us_msk, inplane_R, cap_fov_shift, maxetl] = poisson_mask_read(fname, nky, [p])
%
% Reads in a Poisson-disc-like undersampling mask from a file.
%
% Inputs
%   fname         - File name.
%   nky           - Size of the acquired k-space data (raw k-space data loaded
%                   from p-file. Inplane acceleration interleaved with zeros,
%                   partial ky not zero-padded) in the PE dimension.
%   p             - Parameter structure. See mux_epi_params.m for details.
%                   The following fields are used in this function:
%                   cap_blip_start, inplane_R, kydir, BOTTOM_UP, cap_fov_shift.
%
% Outputs
%   us_msk        - Poisson disk undersample mask with fields 'ky', 'kz', 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   inplane_R     - Inplane reduction factor.
%   cap_fov_shift - CAIPI FOV shift, i.e. number of omegaz lines.
%   maxetl        - Maximum echo train length the file allows.
%
% (c) Kangrong Zhu      Stanford University     July 2014

TOP_DOWN   = 0;
CENTER_OUT = 1;
BOTTOM_UP  = 2;

if ~exist('p', 'var')
    p = struct();
end
if ~isfield(p, 'TOP_DOWN') || isempty(p.TOP_DOWN)
    p.TOP_DOWN = TOP_DOWN;
end
if ~isfield(p, 'CENTER_OUT') || isempty(p.CENTER_OUT)
    p.CENTER_OUT = CENTER_OUT;
end
if ~isfield(p, 'BOTTOM_UP') || isempty(p.BOTTOM_UP)
    p.BOTTOM_UP = BOTTOM_UP;
end

% -- Read file
fid = fopen(fname, 'r', 'l');
inplane_R = fread(fid, 1, 'int32');
cap_fov_shift = fread(fid, 1, 'int32');
maxetl = fread(fid, 1, 'int32');
if maxetl*inplane_R < nky
    error('Number of ky lines exceeds maximum allowed.');
end
if mod(nky, inplane_R) ~= 0
    error('nky must be integer multiples of inplane_R.');
end
etl = nky / inplane_R;
kz_indices = fread(fid, etl, 'int32');
fclose(fid);

% -- Set up parameter structure
% p.cap_blip_start
if ~isfield(p, 'cap_blip_start') || isempty(p.cap_blip_start)
    p.cap_blip_start = floor(cap_fov_shift/2);
end
% p.inplane_R
if isfield(p, 'inplane_R') && (p.inplane_R ~= inplane_R)
    error('Inplane acceleration factor in input parameter structure not the same as in %s.', fname);
end
p.inplane_R = inplane_R;
% p.kydir
if ~isfield(p, 'kydir') || isempty(p.kydir)
    p.kydir = p.BOTTOM_UP;
end
% p.cap_fov_shift
if isfield(p, 'cap_fov_shift') && (abs(p.cap_fov_shift) ~= cap_fov_shift)
    error('cap-fov-shift in input parameter structure not the same as in %s.', fname);
end
if ~isfield(p, 'cap_fov_shift') || isempty(p.cap_fov_shift)
    if ~isfield(p, 'dftz_type') || isempty(p.dftz_type)
        p.dftz_type = 'inverse';
    end
    p.cap_fov_shift = (-1)^strcmp(p.dftz_type, 'inverse') * cap_fov_shift;
end
blip_polarity = sign(p.cap_fov_shift);

% -- Set up undersample mask
% Get a CAIPI undersample mask first
inc_inp = true;
us_msk_caipi = get_caipi_us_msk(p, nky, inc_inp);
% Set up the Poisson disk mask
us_msk.ky = us_msk_caipi.ky;
us_msk.omegaz = blip_polarity * sort(unique(us_msk_caipi.omegaz), 'descend');
us_msk.kz = kz_indices(:).' + 1; % us_msk.kz 1,2,3...cap_fov_shift corresponds to -kzmax ... kzmax
if p.kydir == p.BOTTOM_UP
    us_msk.kz = flipdim(us_msk.kz, 2);
end

return