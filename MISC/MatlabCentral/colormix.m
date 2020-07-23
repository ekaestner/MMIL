function CLR = colormix(CLRSTR,MODE)
% COLORMIX - Combine basic colors with string syntax to create an RGB color 
% 
% CLR = colormix(CLRSTR)
% Convert string CLRSTR with standard MATLAB RGB shorthands to
% an RGB vector. Combines multiple colors, defaults to linear combination.
% I.e returns the mathematical 'average' color.
% 
% CLR = colormix(CLRSTR,MODE)
% Combines colors based on mixing mode MODE. MODE can be 'lincomb' for 
% linear combination, 'add' for additive light mixing (normalized), or 
% 'sub' for subtractive "paint" colormixing.
% See note below.
%
% Available colors: r,g,b,c,m,y,k,w  (Caps insensitive)
%
%--- Usage Example:
% CLR = colormix('ckkr')
% Creates a 1x3 vector CLR describing a color that is one part cyan, 
% two parts black and one part red.
% Alternatively put an integer before the color to specify the amount.
% CLR = colormix('c2kr') is equivalent to the above example.
%
%--- Example 2:
% CLR = colormix('by','sub')
% Produces black: [0 0 0] in the subtractive model, because blue and yellow 
% are compliments of eachother. 
% Similarly, CLR = colormix('by','add') produces white [1 1 1].
%
%------
% Note: Linear is distinct from additve or subtractive color mixing.
% It is similar to subtractive in that adding blacks and low intensity 
% colors lowers total luminosity. It is similar to additive in how 
% combining colors can never make blacks, only desaturating grays.
%
% E.g. 'rc' would be [1 0 0] + [0 1 1] = [1 1 1] for additive
% and [0 0 0] for subtractive. Linear combination gives [0.5 0.5 0.5].
%
% 'rrc' would be [2 1 1] (or [1 0.5 0.5] when normalized) for additive.
% Subtractive gives [0.33 0 0], linear is [0.66 0.33 0.33].
%
% 'rrgb' = [2 1 1] or [0.25 0 0], linear is [0.5 0.25 0.25].


assert(ischar(CLRSTR),...
    'Incorrect input, require a character array as first argument')

CLRSTR = lower(CLRSTR);

%--- Check if CLRSTR contains other characters
if ~all(ismember(CLRSTR,'1234567890rgbcmykw'))
    error('Input string CLRSTR contains illegal character(s)');
end

if nargin >1
    assert(ischar(MODE),...
         'Incorrect input, MODE must be lincomb, add or sub. See help.');
    MODE = lower(MODE);
else
    MODE = 'lincomb'; %default mode
end

%--- Check if using numbers:
if any(ismember(CLRSTR,'1234567890'))
    numbers = true;
else
    numbers = false;
end

%--- Count rgbcmykw occurences
R = 0+sum(ismember(CLRSTR,'r'));
G = 0+sum(ismember(CLRSTR,'g'));
B = 0+sum(ismember(CLRSTR,'b'));
C = 0+sum(ismember(CLRSTR,'c'));
M = 0+sum(ismember(CLRSTR,'m'));
Y = 0+sum(ismember(CLRSTR,'y'));
K = 0+sum(ismember(CLRSTR,'k'));
W = 0+sum(ismember(CLRSTR,'w'));

%--- Parse numbers
if numbers
    count   = zeros(8,1);
    values  = regexp(CLRSTR,'\d+.?','match'); %returns {'44k','2b','9r'} form
    for ii = 1:numel(values)
        input = values{ii};
        switch input(end)
            case 'r'
                count(1) = count(1)-1+ str2double(input(1:end-1));
            case 'g'
                count(2) = count(2)-1+str2double(input(1:end-1));
            case 'b'
                count(3) = count(3)-1+str2double(input(1:end-1));
            case 'c'
                count(4) = count(4)-1+str2double(input(1:end-1));
            case 'm'
                count(5) = count(5)-1+str2double(input(1:end-1));
            case 'y'
                count(6) = count(6)-1+str2double(input(1:end-1));
            case 'k'
                count(7) = count(7)-1+str2double(input(1:end-1));
            case 'w'
                count(8) = count(8)-1+str2double(input(1:end-1));
            otherwise
            error(['Number on the last position of CLRSTR?'...
                ' Or there is a bug in parsing numbers. Sorry!']);
        end
    end        
    R = R + count(1);
    G = G + count(2);
    B = B + count(3);
    C = C + count(4);
    M = M + count(5);
    Y = Y + count(6);
    K = K + count(7);
    W = W + count(8); 
end

%--- Calculate CLR
N   = sum([R G B C M Y K W]); %no. of mix components
CLR = cat(1,R*[1 0 0],...
    G*[0 1 0],...
    B*[0 0 1],...
    C*[0 1 1],...
    Y*[1 1 0],...
    M*[1 0 1],...
    K*[0 0 0],...
    W*[1 1 1]);

if strcmp(MODE,'lincomb')
    CLR = sum(CLR,1)/N;
elseif strcmp(MODE,'add')    
    CLR = sum(CLR,1);    
    CLR = CLR/max(CLR);
elseif strcmp(MODE,'sub')
    CLR = sum(CLR,1);
    CLR = (CLR - min(CLR))/N;      %basic color
    CLR = CLR + ([1 1 1] * W / N); %add whites
else
    error(['Unrecognized input for MODE,'...
        ' expect ''lincomb'', ''add'' or ''sub''.']);
end


%--- uint16 codes
%Numbers 0-9 are: 48-57 
%characters rgbcmykw are: 114,103,98,99,109,121,107,119

%--- Calculation examples
% subtractive mixing: I separate black first.
% 'rrg' = [0 -1 -1][0 -1 -1][-1 0 -1]
%       = [1 0 0][1 1 0][-1 -1 -1]
%       = [0.66 .33.0]
% 'yr'  = [0 0 -1][0 -1 -1] = [1 0.5 0]
%  Notice how mixing green and red produces a darker yellow, compared to 
%  using a pure yellow pigment.

% Additive mixing uses normalization to screen intensity (0-1).
% 'rrg' = [1 0 0][1 0 0][0 1 0] = [2 1 0] = [1 0.5 0]
% 'yr'  = [1 1 0][1 0 0] = [2 1 0] = [1 0.5 0]

% Linear workout:
% 'rrg' = [2 1 0]/3 = [0.66 0.33 0]
% 'yr'  = [2 1 0]/2 = [1 0.5 0]
% Distinct from subtractive in that lower luminosity is not due to 
% inclusion of blacks, but because of averaging over more colors.