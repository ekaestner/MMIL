function us_msk = us_msk_fix_big_gap(us_msk, nky, ker_sup_nky, ky_gap_tol)
%
% function us_msk = us_msk_fix_big_gap(us_msk, nky, [ker_sup_nky=5], [ky_gap_tol=2])
%
% Fix big gaps along ky for an undersampling mask.
%
% Inputs
%   us_msk      - Undersample masks with possible big holes. Having fields 'ky', 'kz', 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   nky         - Total number of ky lines in the undersampling mask.
%   ker_sup_nky - Size of kernel support along the ky dimension, in number of ky lines.
%   ky_gap_tol  - Tolerance for the gap size along ky, with respect to the size of the kernel support. The function will try to fill in holes larger than ky_gap_tol*ker_sup_nky along the ky dimension.
%
% Output
%   us_msk      - Undersample masks with big holes (possibly) filled.
%
% (c) Kangrong Zhu  Stanford University     July 2014

NITER = 3; % Conduct NITER iterations

if ~exist('ker_sup_nky', 'var') || isempty(ker_sup_nky)
    ker_sup_nky = 5;
end

if ~exist('ky_gap_tol', 'var') || isempty(ky_gap_tol)
    ky_gap_tol = 2;
end

nmsk = length(us_msk);
nsamp = length(us_msk(1).ky);
for iter = 1 : NITER
    for msk_idx = 1 : nmsk
        % Put ky indices into monotonically increasing order
        [us_msk(msk_idx).ky, idx] = sort(us_msk(msk_idx).ky);
        us_msk(msk_idx).kz = us_msk(msk_idx).kz(idx);
        
        % Get gap before and after each sample
        [nky_gap_before, nky_gap_after] = get_nky_gap(us_msk(msk_idx), nky);
        
        % Fill in big gaps
        for samp = 1 : nsamp
            % Fill in big gap before a sample
            if nky_gap_before(samp) >= ky_gap_tol*ker_sup_nky % Gap is big enough(larger than or equal to ky_gap_tol*ker_sup_nky) to be filled
                [us_msk(msk_idx), nky_gap_before, nky_gap_after] = fill_one_gap(us_msk(msk_idx), nky_gap_before, nky_gap_after, nky, samp, 'before_samp');
            end
            
            % Fill in big gap after a sample
            later_ky = us_msk.ky(samp) + nky_gap_after(samp) + 1; % ky for the later sample on the same kz line
            if later_ky > nky % Only consider the case when this sample is at the edge, all other cases have already been covered before
                if nky_gap_after(samp) >= ky_gap_tol*ker_sup_nky % Gap is big enough to be filled
                    [us_msk(msk_idx), nky_gap_before, nky_gap_after] = fill_one_gap(us_msk(msk_idx), nky_gap_before, nky_gap_after, nky, samp, 'after_samp');
                end
            end

        end % for samp = 1 : nsamp
    end % for msk_idx = 1 : nmsk
end % for iter = 1 : NITER

return


function [us_msk, nky_gap_before, nky_gap_after] = fill_one_gap(us_msk, nky_gap_before, nky_gap_after, nky, samp, before_or_after_samp)
%
% Fill in one gap in the undersampling mask.
%
% Inputs
%   us_msk         - Undersampling mask. length(us_msk.ky) = nsamp. length(us_msk.omegaz) = nomegaz.
%   nky_gap_before - Gap size before each sample, in number of ky lines. Dim: [nsamp, 1].
%   nky_gap_after  - Gap size after each sample, in number of ky lines. Dim: [nsamp, 1].
%   nky            - Total number of ky lines in the undersampling mask.
%   samp           - Current sample index.
%   before_or_after_samp - 'before_samp': Fill in gap before the current sample.
%                          'after_samp':  Fill in gap after the current sample.
%
% Outputs
%   us_msk         - Undersampling mask, with the current gap (possibly) filled.
%   nky_gap_before - Gap size before each sample, in number of ky lines, for the gap-filled undersampling mask. Dim: [nsamp, 1].
%   nky_gap_after  - Gap size after each sample, in number of ky lines, for the gap-filled undersampling mask. Dim: [nsamp, 1].
%
% (c) Kangrong Zhu  Stanford University     Feb 2015

nsamp = length(us_msk.ky);
nomegaz = length(us_msk.omegaz);

current_ky = us_msk.ky(samp);
previous_ky = current_ky - nky_gap_before(samp) - 1; % ky for the previous sample on the same kz line
later_ky = current_ky + nky_gap_after(samp) + 1; % ky for the later sample on the same kz line

switch before_or_after_samp
    case 'before_samp'
        tgt_ky_pos = (previous_ky + current_ky) / 2;
        current_gap = nky_gap_before(samp);
    case 'after_samp'
        tgt_ky_pos = (later_ky + current_ky) / 2;
        current_gap = nky_gap_after(samp);
end

% Find ky candidates
if nomegaz > 1
    switch before_or_after_samp
        case 'before_samp'

            if previous_ky == 0
                previous_ky_samp_idx = 0;
            else
                previous_ky_samp_idx = find(us_msk.ky == previous_ky);
            end
            ky_cand_samp_idx = previous_ky_samp_idx+1 : samp-1;

        case 'after_samp'
            ky_cand_samp_idx = samp+1 : nsamp;
    end
else % Only one kz line at omegaz=0
    ky_cand_samp_idx = [1 : samp-1, samp+1 : nsamp];
end

ky_cand = us_msk.ky(ky_cand_samp_idx);
gap_to_open = nky_gap_before(ky_cand_samp_idx) + nky_gap_after(ky_cand_samp_idx) + 1; % Gaps that will be opened if these ky lines are moved
[min_gap_to_open, min_gap_to_open_idx] = min(gap_to_open);
ky_cand = ky_cand(min_gap_to_open_idx); % Only consider moving samples that will leave minimum gap

% Move best ky candidate to fill in one gap
if ~isempty(ky_cand) && (min_gap_to_open < current_gap) % Only consider moving a sample if the gap to be opend is smaller than the current one
    if nomegaz > 1
        [us_msk, nky_gap_before, nky_gap_after] = fill_one_gap_for_nomegaz_greater_than_one(us_msk, ky_cand, tgt_ky_pos, nky, samp, current_ky, previous_ky, later_ky, before_or_after_samp);
    else
        [us_msk, nky_gap_before, nky_gap_after] = fill_one_gap_for_nomegaz_equals_one(us_msk, ky_cand, tgt_ky_pos, nky);
    end
end

return


function [us_msk, nky_gap_before, nky_gap_after] = fill_one_gap_for_nomegaz_greater_than_one(us_msk, ky_cand, tgt_ky_pos, nky, samp, current_ky, previous_ky, later_ky, before_or_after_samp)
%
% Fills in one gap for nomegaz > 1 (more than 1 omegaz lines).
%
% Inputs
%   us_msk         - Undersampling mask, before gap filling.
%   ky_cand        - ky candidates that will leave minimum gap if being moved.
%   tgt_ky_pos     - Target ky position, i.e. middle point of the current gap.
%   nky            - Total number of ky lines.
%   samp           - Current sample index.
%   current_ky     - Current ky position.
%   previous_ky    - ky for the previous sample on the same kz line
%   later_ky       - ky for the later sample on the same kz line
%   before_or_after_samp - 'before_samp': Fill in gap before the current sample.
%                          'after_samp':  Fill in gap after the current sample.
%
% Outputs
%   us_msk         - Undersampling mask, with the target gap filled.
%   nky_gap_before - Gap size before each sample, in number of ky lines, for the gap-filled undersampling mask. Dim: [nsamp, 1].
%   nky_gap_after  - Gap size after each sample, in number of ky lines, for the gap-filled undersampling mask. Dim: [nsamp, 1].
%
% (c) Kangrong Zhu  Stanford University     Feb 2015

switch before_or_after_samp
    case 'before_samp'
        new_gap_before = ky_cand - previous_ky - 1; % New gap before the candidate sample if the candidate sample is moved from another kz line to the current kz line
        new_gap_after = current_ky - ky_cand - 1; % New gap after the candidate sample if the candidate sample is moved from another kz line to the current kz line
    case 'after_samp'
        new_gap_before = ky_cand - current_ky - 1;
        new_gap_after = later_ky - ky_cand - 1;
end
new_gap_max = max([new_gap_before(:), new_gap_after(:)], [], 2); % Maximum new gap if the candidate sample is moved from another kz line to the current kz line
ky_cand = ky_cand(new_gap_max == min(new_gap_max)); % All ky indices that will leave a minimum gap and will create minimum maximum new gap
if length(ky_cand) > 1 % If multiple candidates
    dist_to_mid_gap = abs(ky_cand - tgt_ky_pos);
    [~, idx] = sort(dist_to_mid_gap);
    ky_to_move = ky_cand(idx(1)); % Choose the one that's closest to the target ky position (i.e. the middle of the current gap)
else
    ky_to_move = ky_cand;
end

% Move the sample that will leave a minimum gap and will create minimum maximum new gap and is closet to the middle of the current gap
us_msk.kz(us_msk.ky == ky_to_move) = us_msk.kz(samp); % Only change us_msk.kz, us_msk.ky is never changed

[nky_gap_before, nky_gap_after] = get_nky_gap(us_msk, nky);

return


function [us_msk, nky_gap_before, nky_gap_after] = fill_one_gap_for_nomegaz_equals_one(us_msk, ky_cand, tgt_ky_pos, nky)
%
% Fills in one gap for nomegaz = 1 (only 1 omegaz lines).
%
% Inputs
%   us_msk         - Undersampling mask, before gap filling.
%   ky_cand        - ky candidates that will leave minimum gap if being moved.
%   tgt_ky_pos     - Target ky position, i.e. middle point of the current gap.
%   nky            - Total number of ky lines.
%
% Outputs
%   us_msk         - Undersampling mask, with the target gap filled.
%   nky_gap_before - Gap size before each sample, in number of ky lines, for the gap-filled undersampling mask. Dim: [nsamp, 1].
%   nky_gap_after  - Gap size after each sample, in number of ky lines, for the gap-filled undersampling mask. Dim: [nsamp, 1].
%
% (c) Kangrong Zhu  Stanford University     Feb 2015

if isempty(ky_cand)
    error('Empty ky candidates.');
end

ky_to_move = ky_cand( randi(length(ky_cand), [1,1]) );
us_msk.ky = us_msk.ky(us_msk.ky ~= ky_to_move);
us_msk.ky = [us_msk.ky(:); round(tgt_ky_pos)];
us_msk.ky = sort(us_msk.ky);

[nky_gap_before, nky_gap_after] = get_nky_gap(us_msk, nky);

return


function [nky_gap_before, nky_gap_after] = get_nky_gap(us_msk, nky)
%
% Calcualtes the gap along the ky dimension before and after each sample.
%
% Inputs
%   us_msk         - The undersampling mask with fields 'ky', 'kz', 'omegaz'.
%   nky            - Total number of ky lines in the undersampling mask.
%
% Outputs
%   nky_gap_before - Gap size before each sample, in number of ky lines. Dim: [nsamp, 1].
%   nky_gap_after  - Gap size after each sample, in number of ky lines. Dim: [nsamp, 1].
%
% (c) Kangrong Zhu  Stanford University     Feb 2015

BOTTOM_UP  = 2;
nomegaz = length(us_msk(1).omegaz);
us_msk_mtx = us_msk_get_mtx(us_msk, nky, BOTTOM_UP, nomegaz);
nsamp = length(us_msk(1).ky);
nky_gap_before = zeros(nsamp, 1);
nky_gap_after = zeros(nsamp, 1);
for kz_idx = 1 : nomegaz
    kz_line = us_msk_mtx.mtx(:, kz_idx);
    samp_ky_idx = find(kz_line >= 1);
    diff_ky_idx = diff(samp_ky_idx);
    nky_gap_before_on_this_kz_line = [samp_ky_idx(1); diff_ky_idx] - 1;
    nky_gap_after_on_this_kz_line = [diff_ky_idx; nky-samp_ky_idx(end)+1] - 1;
    for this_kz_line_samp_idx = 1 : length(samp_ky_idx)
        samp_idx = find(us_msk.ky == samp_ky_idx(this_kz_line_samp_idx)); % To find the sample index corresponding to the input undersampling mask. Use the fact that each ky line is acquired only once
        nky_gap_before(samp_idx) = nky_gap_before_on_this_kz_line(this_kz_line_samp_idx);
        nky_gap_after(samp_idx) = nky_gap_after_on_this_kz_line(this_kz_line_samp_idx);
    end
end

return