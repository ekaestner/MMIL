function [t_str, t_num] = get_current_time(rs, fd)
%
% Get the current time, return both strings and numbers.
%
% Inputs
%   rs    - True: Round the seconds to an integer value. Default: true.
%   fd    - True: Use fixed number of digits in the time strings. Default: true.
%
% Outputs
%   t_str - Strings for the current time.
%           Element  - field, number of digits if 'fd' is true.
%           t_str{1} - year,    4
%           t_str{2} - month,   2
%           t_str{3} - day,     2
%           t_str{4} - hour,    2
%           t_str{5} - minute,  2
%           t_str{6} - seconds, 2
%   t_num - Numbers for the current time. A 1-by-6 vector with elements
%           [year month day hour minute seconds].
%
% (c) Kangrong Zhu,     Stanford University,    June 2013

if ~exist('rs', 'var') || isempty(rs)
    rs = true;
end

if ~exist('fd', 'var') || isempty(fd)
    fd = true;
end

ne = 6;                      % Number of elements in the time record
t_num = clock;               % current time

if rs                        % If round the seconds to an integer value
    t_num = fix(t_num);
end

if fd                        % If use fixed number of digits in the time strings
    nd = [4, 2, 2, 2, 2, 2]; % Number of expected digits in each string
end

t_str{ne} = '';
for idx = 1 : ne
    t_str{idx} = num2str(t_num(idx));
    
    if fd
        nds = nd(idx) - length(t_str{idx}); % How many digits short
        if nds > 0
            t_str{idx} = [repmat('0', [1, nds]), t_str{idx}];
        end
    end
end

