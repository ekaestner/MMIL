function ff = gen_fermi_filter(np, fr, fw, ty)
%
% function ff = gen_fermi_filter(np, fr, fw [, ty])
%
% Generates Fermi filter.
%
% Inputs
%   np - Number of pixels, [np(1), np(2)(optional)].
%   fr - Fermi radius, in pixels.
%   fw - Fermi width, in pixels.
%   ty - Filter type to generate when ( numel(np)==2 && np(2)>1 ).
%        'circ'(default): Circularly symmetric.
%        'sepa':          Separable.
%
% Output
%   ff - Fermi filter.
%        If (numel(np)==1 || (numel(np)==2 && np(2)==1)), ff is a 1D Fermi filter of size [np(1), 1].
%        If (numel(np)==2 && np(2)>1 && strcmp(ty, 'circ')), ff is a 2D circularly symmetric Fermi filter of size [np(1), np(2)].
%        If (numel(np)==2 && np(2)>1 && strcmp(ty, 'sepa')), ff is a 2D separable Fermi filter of size [np(1), np(2)].
%
% (c) Kangrong Zhu,     Stanford University     Sep 2013

if ~exist('ty', 'var') || isempty(ty)
    ty = 'circ';
end

if (numel(np)==1) || ( (numel(np)==2) && (np(2)==1) ) % 1D Fermi filter
    ff = fermf_1d(np(1), fr, fw);
end

if (numel(np) == 2) && (np(2) > 1)                    % 2D Fermi filter
    switch ty
        case 'circ'                                   % Circularly symmetric filter
            ff = zeros(np(1), np(2));
            or = np/2 + 1;                            % Origin
            for x = 1 : np(1)
                for y = 1 : np(2)
                    d = sqrt((x-or(1))^2 + (y-or(2))^2);
                    ff(x, y) = ( 1 + exp((d-fr)/fw) ) .^ (-1);
                end
            end
        case 'sepa'                                   % Separable filter
            ff = fermf_1d(np(1), fr, fw);
            ff2 = fermf_1d(np(2), fr, fw);
            ff = repmat(ff, [1, np(2)]) .* repmat(ff2.', [np(1), 1]);
    end
end

return

function f = fermf_1d(n, r, w)
%
%   n - Number of pixels.
%   r - Fermi radius, in pixels.
%   w - Fermi width, in pixels.
%
%   f - 1D Fermi filter.
%

d = abs( -floor(n/2) : 1 : ceil(n/2-1) );
d = d(:);
f = ( 1 + exp((d-r)/w) ) .^ (-1);

return