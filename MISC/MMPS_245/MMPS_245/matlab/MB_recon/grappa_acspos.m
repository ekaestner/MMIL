function [acs_pos, sz_acs] = grappa_acspos(sz_dat, sz_acs, sz_lim, debug)
%
% function [acs_pos, sz_acs] = grappa_acspos(sz_dat, sz_acs[, sz_lim, debug])
%
% Calculates the indices of the Autocalibration area for GRAPPA reconstruction. 
%
% Inputs:
%   sz_dat  - Size of the whole k-space matrix. A structure array with fields 'x' and 'y'.
%   sz_acs  - Size of the ACS area. A structure array with fields 'x' and 'y'.
%   sz_lim  - Limits on acs_pos. A structure array with fields 'x' and/or 'y'.
%             Will check the x/y indices if field x/y exists and is non-empty.
%   debug   - True: print debug messages. Only needed when sz_lim has a
%             non-empty 'x' and/or 'y' field.
%
% Output:
%   acs_pos - Indices for the ACS area. A structure array with fields 'x' and 'y'.
%   sz_acs  - Size of the ACS area, having the same fields as the input 'sz_acs'.
%             The values might have been changed in this function if partial k is used.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

if ~exist('sz_lim', 'var')
    sz_lim = [];
end

if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

acs_pos = struct('x', {[]}, 'y', {[]});

acs_pos.x = get_acs_pos(sz_dat.x, sz_acs.x);
acs_pos.y = get_acs_pos(sz_dat.y, sz_acs.y);

% Remove the out-of-range k-space indices. Useful for partial k acquisitions.
if ~isempty(sz_lim)
    if isfield(sz_lim, 'x') && ~isempty(sz_lim.x)
        [acs_pos.x, sz_acs.x] = checkrange(acs_pos.x, sz_lim.x(1), sz_lim.x(2), 'ACS kx positions', debug); % acs_pos.x, sz_acs.x change if partial kx 
    end
    if isfield(sz_lim, 'y') && ~isempty(sz_lim.y)
        [acs_pos.y, sz_acs.y] = checkrange(acs_pos.y, sz_lim.y(1), sz_lim.y(2), 'ACS ky positions', debug); % acs_pos.y, sz_acs.y change if partial ky
    end
end

function pos = get_acs_pos(n_total, n_acs)

pos = floor(n_total/2) + ( -floor(n_acs/2)+1 : ceil(n_acs/2) );
