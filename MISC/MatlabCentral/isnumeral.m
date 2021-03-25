function flg = isnumeral(s)
 %ISNUMERAL True for numeral 0-9.
 % For a string S, ISNUMERAL(S) returns an array that
 % contains 1 for numbers and 0 where not.
 %
 % See also ISLETTER, ISSPACE.

 if ~ischar(s)
 error('isnumeral argument not character')
 end

 flg = ismember(s,'0123456789');