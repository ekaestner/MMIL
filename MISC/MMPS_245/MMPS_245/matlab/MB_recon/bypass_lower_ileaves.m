function indices = bypass_lower_ileaves(indices, p)
%
% The PSD conducts in-plane acceleration by by-passing the lower interleaves.
% Do the same thing here for a 1D vector of indices.
%
% Inputs
%   indices - 1D vector of indices.
%   p       - Parameter structure. See mux_epi_params.m for details.
%             The following fields are used in this function: inplane_R,
%             kydir, BOTTOM_UP, TOP_DOWN, CENTER_OUT.
%
% Output
%   indices - Vector of indices with the by-passed indices taken out.
%
% (c) Kangrong Zhu      Stanford University     Dec 2013

idx = 1 : 1 : length(indices);
switch p.kydir
    case p.BOTTOM_UP
        for ileave = 0 : p.inplane_R-2
            idx(end-ileave : -p.inplane_R : 1) = 0;
        end
    case p.TOP_DOWN
        for ileave = 0 : p.inplane_R-2
            idx(ileave+1 : p.inplane_R : end) = 0;
        end
    case p.CENTER_OUT
        error('Center-Out k-space acquisition is not supported in recon yet.'); % TODO: support center-out k-space acquisition
end
indices = indices(idx ~= 0);

return